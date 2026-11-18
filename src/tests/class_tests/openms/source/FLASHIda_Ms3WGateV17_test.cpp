// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// ADR-0079: the MS3 gate W of package 17 against the contract of its reference values (CONTRACT.md of the ms3_esm package):
// level 1 (the eight features, the 38 inputs and W per recorded MS3) and level 2 (the decision at 0.3, 0.5 and 0.8 per
// recorded MS3), on all 12,440 MS3 of the 14 recorded LC runs, with the engine's own classes (Ms1Buffer, planPeakFeatures,
// w17ModelInputs, WitnessModel). The data is PRIVATE and lives in the parent repo's test-data-private/ms3_w_gate_v17/
// (README.md there); a missing file FAILS, never skips. level12_full_decisions writes the decisions for the level-3 step
// of the CI (contract/level3.py), and level12_full_buffer_features prints the cost per plan.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/LogCell.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms1Buffer.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3WitnessModel.h>
#include <OpenMS/FORMAT/GzipIfstream.h>

#include "FLASHIda_TestHelpers.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace OpenMS;

namespace
{
  const std::string kData = "../../test-data-private/ms3_w_gate_v17/";
  const std::string kModelSha256 = "5d3f44a99ed580c6557605ee6433ea44f4a8b6cfeac857072ca09f8224f755d3";
  const double kThresholds[3] = {0.3, 0.5, 0.8};
  const char* const kThresholdNames[3] = {"0.3", "0.5", "0.8"};

  std::string readAll(const std::string& path)
  {
    if (path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0)
    {
      GzipIfstream gz(path.c_str());
      if (!gz.isOpen()) return std::string();
      std::string out;
      std::vector<char> buf(1 << 20);
      while (!gz.streamEnd())
      {
        const size_t n = gz.read(buf.data(), buf.size());
        out.append(buf.data(), n);
        if (n == 0) break;
      }
      return out;
    }
    std::ifstream in(path, std::ios::binary);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
  }

  /// A tab-separated table with a header; CSV-quoted cells (tracking ids that begin with a quote) unquoted.
  struct Table
  {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    std::unordered_map<std::string, size_t> col;

    static Table read(const std::string& path)
    {
      Table t;
      const std::string text = readAll(path);
      size_t pos = 0;
      bool first = true;
      while (pos < text.size())
      {
        size_t eol = text.find('\n', pos);
        if (eol == std::string::npos) eol = text.size();
        std::string line = text.substr(pos, eol - pos);
        pos = eol + 1;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> cells;
        size_t s = 0;
        while (true)
        {
          const size_t tab = line.find('\t', s);
          cells.push_back(unquoteLogCell(line.substr(s, tab == std::string::npos ? std::string::npos : tab - s)));
          if (tab == std::string::npos) break;
          s = tab + 1;
        }
        if (first)
        {
          t.header = cells;
          for (size_t i = 0; i < cells.size(); ++i) t.col[cells[i]] = i;
          first = false;
        }
        else
          t.rows.push_back(std::move(cells));
      }
      return t;
    }
    const std::string& at(const std::vector<std::string>& r, const std::string& name) const
    {
      static const std::string empty;
      auto it = col.find(name);
      if (it == col.end() || it->second >= r.size()) return empty;
      return r[it->second];
    }
    double num(const std::vector<std::string>& r, const std::string& name) const
    {
      const std::string& c = at(r, name);
      if (c.empty()) return std::numeric_limits<double>::quiet_NaN();
      try
      {
        return std::stod(c);
      }
      catch (const std::exception&)
      {
        return std::numeric_limits<double>::quiet_NaN();
      }
    }
  };

  /// a and b equal as 32-bit floats (both NaN counts as equal), or within @p ulps units in the last place.
  bool f32Equal(float a, float b, int ulps = 0)
  {
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
    if (a == b) return true;
    if (ulps == 0) return false;
    int32_t ia, ib;
    std::memcpy(&ia, &a, 4);
    std::memcpy(&ib, &b, 4);
    if ((ia < 0) != (ib < 0)) return false;
    return std::abs(ia - ib) <= ulps;
  }

  bool near(double a, double b, double tol)
  {
    if (std::isnan(a) || std::isnan(b)) return std::isnan(a) && std::isnan(b);
    return std::fabs(a - b) <= tol;
  }

  const WitnessModel& model()
  {
    static const WitnessModel m = WitnessModel::load(kData + "w_v13.json", kW17InputCount);
    return m;
  }

  /// The 29 raw values of a level12 / reference row (the package's RAW names; frag_type N/C -> frag_N).
  std::array<double, kW17RawCount> rawOf(const Table& t, const std::vector<std::string>& r)
  {
    std::array<double, kW17RawCount> raw{};
    const std::string& ft = t.at(r, "frag_type");
    raw[kW17FragN] = ft == "N" ? 1.0 : (ft == "C" ? 0.0 : std::numeric_limits<double>::quiet_NaN());
    for (int i = 1; i < kW17RawCount; ++i) raw[static_cast<size_t>(i)] = t.num(r, w17RawName(i));
    return raw;
  }

  const char* const kPk3[8] = {"pk3_since_start", "pk3_since_apex", "pk3_since_last", "pk3_missed",
                               "pk3_rel_max",     "pk3_log_max",    "pk3_trend",      "pk3_n_seen"};

  bool isLog10Input(int i)
  {
    const std::string n = w17InputName(i);
    return n == "id_score" || n == "x_prec_intensity" || n == "x_ms1_plan" || n == "x_ms2_mass_count";
  }

  /// The per-MS3 gate of the reference (w_gate_v17.decide), for a computed W.
  bool sends(double w, int n_uncovered, double threshold) { return n_uncovered >= 1 && !std::isnan(w) && w >= threshold; }

  /// The pruned surveys of all 14 runs, by run, in resolve_ts order.
  std::map<std::string, std::vector<BufferedSurvey>> loadSurveys()
  {
    std::map<std::string, std::vector<BufferedSurvey>> by_run;
    const Table s = Table::read(kData + "pruned_surveys/surveys.tsv.gz");
    std::map<std::string, std::vector<size_t>> index_by_run;   // survey_index -> position in by_run[run]
    for (const auto& r : s.rows)
    {
      BufferedSurvey b;
      b.tracking_id = s.at(r, "tracking_id");
      b.resolve_ts = static_cast<uint64_t>(std::stoull(s.at(r, "resolve_ts")));
      b.kind = SurveyKind{s.at(r, "scan_type"), s.at(r, "faims_cv"), s.at(r, "faims_enabled"), s.at(r, "first_mass"), s.at(r, "last_mass")};
      b.ion_current = std::stod(s.at(r, "ion_current"));
      b.centroid_count = static_cast<uint64_t>(std::stoull(s.at(r, "centroids_full")));
      const std::string& run = s.at(r, "run");
      const size_t idx = static_cast<size_t>(std::stoul(s.at(r, "survey_index")));
      auto& v = by_run[run];
      auto& ix = index_by_run[run];
      if (ix.size() <= idx) ix.resize(idx + 1, static_cast<size_t>(-1));
      ix[idx] = v.size();
      v.push_back(std::move(b));
    }
    for (auto& [run, v] : by_run)
    {
      const std::string text = readAll(kData + "pruned_surveys/centroids_" + run + ".tsv.gz");
      std::istringstream in(text);
      std::string line;
      std::getline(in, line);   // header
      while (std::getline(in, line))
      {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const size_t t1 = line.find('\t'), t2 = line.find('\t', t1 + 1);
        const size_t idx = static_cast<size_t>(std::stoul(line.substr(0, t1)));
        BufferedSurvey& b = v[index_by_run[run][idx]];
        b.mz.push_back(std::stod(line.substr(t1 + 1, t2 - t1 - 1)));
        b.intensity.push_back(std::stod(line.substr(t2 + 1)));
      }
      std::stable_sort(v.begin(), v.end(), [](const BufferedSurvey& a, const BufferedSurvey& b) { return a.resolve_ts < b.resolve_ts; });
    }
    return by_run;
  }

  double percentile(std::vector<double> v, double q)
  {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t k = static_cast<size_t>(std::min<double>(static_cast<double>(v.size() - 1), std::floor(q * static_cast<double>(v.size() - 1) + 0.5)));
    return v[k];
  }
} // namespace

START_TEST(FLASHIda_Ms3WGateV17, "$Id$")

START_SECTION(model_sha_and_structure)
{
  TEST_EQUAL(model().treeCount(), 658)
  TEST_EQUAL(model().numFeatures(), 38)
  TEST_STRING_EQUAL(model().sha256(), kModelSha256)
  TEST_STRING_EQUAL(model().sha8(), "5d3f44a9")
  // The file declares 38 inputs: read as a 36-input model it is refused.
  TEST_EXCEPTION(std::invalid_argument, WitnessModel::load(kData + "w_v13.json", 36))
}
END_SECTION

START_SECTION(reference_rows_give_w)
{
  const Table t = Table::read(kData + "tests/reference_rows.tsv");
  TEST_EQUAL(t.rows.size(), 329)
  double worst = 0.0;
  int decisions = 0, agree = 0;
  for (const auto& r : t.rows)
  {
    std::array<float, kW17InputCount> x{};
    for (int i = 0; i < kW17InputCount; ++i) x[static_cast<size_t>(i)] = toFloat32(t.num(r, w17InputName(i)));
    const double w = model().w(x);
    worst = std::max(worst, std::fabs(w - t.num(r, "expected_w")));
    const int nu = static_cast<int>(t.num(r, "n_uncovered"));
    for (int k = 0; k < 3; ++k)
    {
      ++decisions;
      if (static_cast<int>(sends(w, nu, kThresholds[k])) == static_cast<int>(t.num(r, std::string("send_") + kThresholdNames[k]))) ++agree;
    }
  }
  std::cout << "[W17] reference_rows: largest |W - expected| = " << worst << ", decisions " << agree << " / " << decisions << std::endl;
  TEST_EQUAL(worst <= 1e-6, true)
  TEST_EQUAL(agree, decisions)
}
END_SECTION

START_SECTION(level12_full_inputs_to_w_and_transforms)
{
  const Table t = Table::read(kData + "contract/level12_full.tsv.gz");
  TEST_EQUAL(t.rows.size(), 12440)
  double worst_w = 0.0;
  int input_mismatch = 0, ulp_differences = 0;
  for (const auto& r : t.rows)
  {
    // Level 1, the trees: the 38 recorded inputs -> W.
    std::array<float, kW17InputCount> x{};
    for (int i = 0; i < kW17InputCount; ++i) x[static_cast<size_t>(i)] = toFloat32(t.num(r, std::string("in_") + w17InputName(i)));
    worst_w = std::max(worst_w, std::fabs(model().w(x) - t.num(r, "w")));
    // Level 1, the transforms (SPEC 7): the 29 raw values and the 64-bit features -> the 38 inputs, equal as 32-bit
    // floats; 1 ulp allowed for the four log10 inputs (tolerances.json).
    std::array<double, 8> pk3{};
    for (size_t k = 0; k < 8; ++k) pk3[k] = t.num(r, std::string(kPk3[k]) + "_f64");
    const std::array<float, kW17InputCount> mine = w17ModelInputs(rawOf(t, r), pk3);
    for (int i = 0; i < kW17InputCount; ++i)
    {
      if (f32Equal(mine[static_cast<size_t>(i)], x[static_cast<size_t>(i)])) continue;
      if (isLog10Input(i) && f32Equal(mine[static_cast<size_t>(i)], x[static_cast<size_t>(i)], 1))
      {
        ++ulp_differences;
        continue;
      }
      if (input_mismatch < 10)
        std::cout << "[W17] input mismatch " << t.at(r, "run") << " " << t.at(r, "tracking_id") << " " << w17InputName(i) << ": "
                  << mine[static_cast<size_t>(i)] << " vs " << x[static_cast<size_t>(i)] << std::endl;
      ++input_mismatch;
    }
  }
  std::cout << "[W17] level12_full: largest |W - reference| = " << worst_w << "; inputs off: " << input_mismatch
            << "; 1-ulp log10 differences: " << ulp_differences << " (the package: 43)" << std::endl;
  TEST_EQUAL(worst_w <= 1e-6, true)
  TEST_EQUAL(input_mismatch, 0)
}
END_SECTION

START_SECTION(testvectors_buffer_to_decision)
{
  const nlohmann::json j = nlohmann::json::parse(readAll(kData + "tests/testvectors.json.gz"));
  const auto& vectors = j.at("vectors");
  TEST_EQUAL(vectors.size(), 81)
  int feature_off = 0, input_off = 0, w_off = 0, decision_off = 0;
  for (const auto& v : vectors)
  {
    // The vector's buffer IS the selection (already the kind and horizon of the plan).
    std::vector<BufferedSurvey> owned;
    owned.reserve(v.at("buffer").size());
    Ms1Buffer::Selection sel;
    for (const auto& b : v.at("buffer"))
    {
      BufferedSurvey s;
      s.resolve_ts = static_cast<uint64_t>(b.at("resolve_ts_ms").get<double>());
      s.ion_current = b.at("ion_current").get<double>();
      s.mz = b.at("mz").get<std::vector<double>>();
      s.intensity = b.at("intensity").get<std::vector<double>>();
      sel.lo = b.at("scan_range")[0].get<double>();
      sel.hi = b.at("scan_range")[1].get<double>();
      owned.push_back(std::move(s));
    }
    for (const auto& s : owned) sel.surveys.push_back(&s);
    const uint64_t t = static_cast<uint64_t>(v.at("enqueue_ts_ms").get<double>());
    const PeakFeatures f = planPeakFeatures(sel, v.at("precursor_mono_mass").get<double>(), t);
    const auto exp_f = v.at("expected_features_f64").get<std::vector<double>>();
    for (size_t k = 0; k < 8; ++k)
      if (!near(f[k], exp_f[k], 1e-9) || !f32Equal(toFloat32(f[k]), toFloat32(exp_f[k]))) ++feature_off;

    std::array<double, kW17RawCount> raw{};
    const auto& ri = v.at("raw_inputs");
    const std::string ft = ri.at("frag_type").get<std::string>();
    raw[kW17FragN] = ft == "N" ? 1.0 : 0.0;
    for (int i = 1; i < kW17RawCount; ++i)
    {
      const auto& c = ri.at(w17RawName(i));
      raw[static_cast<size_t>(i)] = c.is_number() ? c.get<double>() : std::numeric_limits<double>::quiet_NaN();
    }
    const std::array<float, kW17InputCount> x = w17ModelInputs(raw, f);
    const auto exp_x = v.at("expected_inputs_f32").get<std::vector<double>>();
    for (int i = 0; i < kW17InputCount; ++i)
      if (!f32Equal(x[static_cast<size_t>(i)], toFloat32(exp_x[static_cast<size_t>(i)]), isLog10Input(i) ? 1 : 0)) ++input_off;
    const double w = model().w(x);
    if (!near(w, v.at("expected_w").get<double>(), 1e-6)) ++w_off;
    const int nu = v.at("n_uncovered").get<int>();
    for (int k = 0; k < 3; ++k)
      if (static_cast<int>(sends(w, nu, kThresholds[k])) != v.at("expected_send").at(kThresholdNames[k]).get<int>()) ++decision_off;
  }
  TEST_EQUAL(feature_off, 0)
  TEST_EQUAL(input_off, 0)
  TEST_EQUAL(w_off, 0)
  TEST_EQUAL(decision_off, 0)
}
END_SECTION

// Level 1 of the features and level 2, on all 12,440 recorded MS3: the pruned surveys of each run fed in resolve order
// into the engine's Ms1Buffer at their recorded resolve_ts, each MS3 selected at its recorded plan time (enqueue_ts) from
// the survey its parent MS2 was picked from, the features computed by planPeakFeatures, the inputs by w17ModelInputs.
START_SECTION(level12_full_buffer_features_and_decisions)
{
  const Table t = Table::read(kData + "contract/level12_full.tsv.gz");
  const auto surveys = loadSurveys();
  TEST_EQUAL(surveys.size(), 14)

  // Rows by run, in plan-time order (stable: a plan's MS3 keep their file order).
  std::map<std::string, std::vector<size_t>> by_run;
  for (size_t i = 0; i < t.rows.size(); ++i) by_run[t.at(t.rows[i], "run")].push_back(i);
  std::vector<double> w_of(t.rows.size(), std::numeric_limits<double>::quiet_NaN());
  int descriptor_off = 0, feature_off = 0, input_off = 0, w_off = 0, rows_seen = 0;
  double worst_feature = 0.0;
  std::vector<double> plan_ms;   // features once + W of every MS3, per plan
  for (auto& [run, idx] : by_run)
  {
    std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return t.num(t.rows[a], "enqueue_ts") < t.num(t.rows[b], "enqueue_ts"); });
    const auto sit = surveys.find(run);
    if (sit == surveys.end())
    {
      std::cout << "[W17] no surveys for " << run << std::endl;
      ++descriptor_off;
      continue;
    }
    const std::vector<BufferedSurvey>& sv = sit->second;
    Ms1Buffer buffer;
    size_t next = 0;
    size_t i = 0;
    while (i < idx.size())
    {
      // One plan: the MS3 that share the plan time, the parent and the precursor mass.
      const auto& r0 = t.rows[idx[i]];
      const uint64_t plan_t = static_cast<uint64_t>(t.num(r0, "enqueue_ts"));
      size_t end = i + 1;
      while (end < idx.size() && static_cast<uint64_t>(t.num(t.rows[idx[end]], "enqueue_ts")) == plan_t
             && t.at(t.rows[idx[end]], "parent_tracking_id") == t.at(r0, "parent_tracking_id"))
        ++end;
      while (next < sv.size() && sv[next].resolve_ts <= plan_t) buffer.add(sv[next++]);

      const auto t0 = std::chrono::steady_clock::now();
      const Ms1Buffer::Selection sel = buffer.select(t.at(r0, "survey_tracking_id"), plan_t);
      const PeakFeatures f = planPeakFeatures(sel, t.num(r0, "mass_ms2"), plan_t);
      for (size_t q = i; q < end; ++q)
      {
        const auto& r = t.rows[idx[q]];
        ++rows_seen;
        w_of[idx[q]] = model().w(w17ModelInputs(rawOf(t, r), f));
      }
      plan_ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());

      for (size_t q = i; q < end; ++q)
      {
        const auto& r = t.rows[idx[q]];
        // The buffer descriptor first, so a selection error shows before the features.
        const bool desc_ok = sel.kind_from == t.at(r, "buffer_kind_from") && sel.surveys.size() == static_cast<size_t>(t.num(r, "buffer_surveys"))
                             && sel.lo == t.num(r, "scan_range_lo") && sel.hi == t.num(r, "scan_range_hi")
                             && sel.centroids == static_cast<uint64_t>(t.num(r, "buffer_centroids"))
                             && (sel.surveys.empty()
                                 || (static_cast<double>(sel.surveys.front()->resolve_ts) == t.num(r, "buffer_first_resolve_ts")
                                     && static_cast<double>(sel.surveys.back()->resolve_ts) == t.num(r, "buffer_last_resolve_ts")
                                     && sel.surveys.front()->tracking_id == t.at(r, "buffer_first_survey_id")
                                     && sel.surveys.back()->tracking_id == t.at(r, "buffer_last_survey_id")));
        if (!desc_ok)
        {
          if (descriptor_off < 10)
            std::cout << "[W17] buffer descriptor off " << run << " " << t.at(r, "tracking_id") << ": surveys " << sel.surveys.size()
                      << " vs " << t.at(r, "buffer_surveys") << ", kind_from " << sel.kind_from << std::endl;
          ++descriptor_off;
        }
        for (size_t k = 0; k < 8; ++k)
        {
          const double e = t.num(r, std::string(kPk3[k]) + "_f64");
          worst_feature = std::max(worst_feature, std::fabs(f[k] - e));
          if (!near(f[k], e, 1e-9) || !f32Equal(toFloat32(f[k]), toFloat32(e)))
          {
            if (feature_off < 10)
              std::cout << "[W17] feature off " << run << " " << t.at(r, "tracking_id") << " " << kPk3[k] << ": " << f[k] << " vs " << e << std::endl;
            ++feature_off;
          }
        }
        std::array<float, kW17InputCount> x_ref{};
        for (int k = 0; k < kW17InputCount; ++k) x_ref[static_cast<size_t>(k)] = toFloat32(t.num(r, std::string("in_") + w17InputName(k)));
        const std::array<float, kW17InputCount> x = w17ModelInputs(rawOf(t, r), f);
        for (int k = 0; k < kW17InputCount; ++k)
          if (!f32Equal(x[static_cast<size_t>(k)], x_ref[static_cast<size_t>(k)], isLog10Input(k) ? 1 : 0)) ++input_off;
        if (!near(w_of[idx[q]], t.num(r, "w"), 1e-6)) ++w_off;
      }
      i = end;
    }
  }
  std::cout << "[W17] level 1 on " << rows_seen << " MS3 in " << plan_ms.size() << " plans: descriptors off " << descriptor_off
            << ", features off " << feature_off << " (largest |difference| " << worst_feature << "), inputs off " << input_off
            << ", W off " << w_off << std::endl;
  std::cout << "[W17] runtime per plan (features once + W of every MS3), ms: median " << percentile(plan_ms, 0.5) << ", p99 "
            << percentile(plan_ms, 0.99) << ", max " << percentile(plan_ms, 1.0) << std::endl;
  TEST_EQUAL(rows_seen, 12440)
  TEST_EQUAL(descriptor_off, 0)
  TEST_EQUAL(feature_off, 0)
  TEST_EQUAL(input_off, 0)
  TEST_EQUAL(w_off, 0)

  // Level 2. Per recorded MS3 with the reference n_uncovered: the decision of the engine's W at each threshold equals
  // the reference's, 100%. Then the engine's rule for a CE sweep (ADR-0079 decision 5): one verdict per target on the
  // largest W of its variants; it may differ from the per-MS3 reference only on a target whose variants the reference
  // decides apart, and every such target is listed.
  const std::string out_dir = freshLogDir("w17_contract");
  for (int k = 0; k < 3; ++k)
  {
    const double th = kThresholds[k];
    const std::string tn = kThresholdNames[k];
    std::map<std::string, double> target_w;   // run|parent|ion_type|ion_index -> max W
    auto key = [&](const std::vector<std::string>& r) {
      return t.at(r, "run") + "|" + t.at(r, "parent_tracking_id") + "|" + t.at(r, "ion_type") + "|" + t.at(r, "ion_index");
    };
    for (size_t q = 0; q < t.rows.size(); ++q)
    {
      double& m = target_w.emplace(key(t.rows[q]), -1.0).first->second;
      m = std::max(m, w_of[q]);
    }
    std::map<std::string, std::set<int>> ref_by_target;
    for (const auto& r : t.rows) ref_by_target[key(r)].insert(static_cast<int>(t.num(r, "send_" + tn)));

    int per_ms3_agree = 0, edge = 0, target_differs = 0, target_differs_unexplained = 0;
    std::ofstream wdec(out_dir + "/w_decisions_" + tn + ".tsv"), edec(out_dir + "/engine_decisions_" + tn + ".tsv");
    wdec << "run\ttracking_id\tsend\n";
    edec << "run\ttracking_id\tsend\n";
    for (size_t q = 0; q < t.rows.size(); ++q)
    {
      const auto& r = t.rows[q];
      const int nu = static_cast<int>(t.num(r, "n_uncovered"));
      const int ref = static_cast<int>(t.num(r, "send_" + tn));
      const int per_ms3 = sends(w_of[q], nu, th) ? 1 : 0;
      if (per_ms3 == ref) ++per_ms3_agree;
      if (std::fabs(w_of[q] - th) <= 1e-6) ++edge;
      const int engine = sends(target_w[key(r)], nu, th) ? 1 : 0;
      if (engine != ref)
      {
        ++target_differs;
        if (ref_by_target[key(r)].size() < 2) ++target_differs_unexplained;
        else std::cout << "[W17] level 2 at " << tn << ": " << t.at(r, "run") << " " << t.at(r, "tracking_id") << " (" << t.at(r, "ion_type")
                       << t.at(r, "ion_index") << ") sent with its sweep's largest W " << target_w[key(r)] << "; its own W " << w_of[q] << std::endl;
      }
      wdec << quoteLogCell(t.at(r, "run")) << "\t" << quoteLogCell(t.at(r, "tracking_id")) << "\t" << per_ms3 << "\n";
      edec << quoteLogCell(t.at(r, "run")) << "\t" << quoteLogCell(t.at(r, "tracking_id")) << "\t" << engine << "\n";
    }
    std::cout << "[W17] level 2 at " << tn << ": per MS3 " << per_ms3_agree << " / " << t.rows.size() << " (edge band " << edge
              << "); one verdict per target: " << target_differs << " MS3 differ, all on targets the reference decides apart: "
              << (target_differs_unexplained == 0 ? "yes" : "NO") << std::endl;
    TEST_EQUAL(per_ms3_agree, static_cast<int>(t.rows.size()))
    TEST_EQUAL(target_differs_unexplained, 0)
  }
}
END_SECTION

START_SECTION(buffer_horizon_kind_fallback)
{
  const SurveyKind a{"survey", "0", "0", "500.0000", "1800.0000"};
  const SurveyKind b{"survey", "-45", "1", "400.0000", "2000.0000"};
  auto survey = [](const std::string& id, uint64_t ts, const SurveyKind& k) {
    BufferedSurvey s;
    s.tracking_id = id;
    s.resolve_ts = ts;
    s.kind = k;
    s.mz = {600.0};
    s.intensity = {1.0};
    s.ion_current = 1.0;
    s.centroid_count = 7;
    return s;
  };
  Ms1Buffer buf;
  TEST_EQUAL(buf.select("!!#", 1000).kind_from, "none")   // nothing resolved: an empty selection, (500, 1800)
  TEST_EQUAL(buf.select("!!#", 1000).surveys.size(), 0)
  TEST_REAL_SIMILAR(buf.select("!!#", 1000).lo, 500.0)
  buf.add(survey("!!#", 100000, a));
  buf.add(survey("!!$", 150000, b));
  buf.add(survey("!!%", 280000, a));
  // Both ends of [t - 180000, t] are included.
  Ms1Buffer::Selection s = buf.select("!!#", 280000);
  TEST_EQUAL(s.kind_from, "parent survey")
  TEST_EQUAL(s.surveys.size(), 2)                         // !!# at t - 180000 and !!% at t; !!$ is of kind b
  TEST_EQUAL(s.centroids, 14)
  s = buf.select("!!#", 280001);
  TEST_EQUAL(s.surveys.size(), 1)                         // !!# is 180001 ms old
  // The parent's survey picks the kind and the range.
  s = buf.select("!!$", 280000);
  TEST_EQUAL(s.surveys.size(), 1)
  TEST_REAL_SIMILAR(s.lo, 400.0)
  TEST_REAL_SIMILAR(s.hi, 2000.0)
  // An unknown parent survey: the kind of the last survey resolved at or before t.
  s = buf.select("zzz", 279999);
  TEST_EQUAL(s.kind_from, "last survey")
  TEST_EQUAL(s.surveys.size(), 1)                         // the last at or before t is !!$ (kind b): itself
  TEST_REAL_SIMILAR(s.lo, 400.0)
  // Eviction: a survey more than 180 s older than the newest is dropped; its kind stays known.
  buf.add(survey("!!&", 330001, b));
  TEST_EQUAL(buf.size(), 2)                               // !!# and !!$ are gone
  s = buf.select("!!#", 330001);
  TEST_EQUAL(s.kind_from, "parent survey")
  TEST_EQUAL(s.surveys.size(), 1)                         // !!% only
  // No sighting: 180, 180, 180, the number of surveys, 0, 0, 0, 0; an empty buffer gives 0 surveys.
  const PeakFeatures none = peakFeatures({1.0, 2.0, 3.0}, {0.0, 0.0, 0.0}, 10.0);
  TEST_REAL_SIMILAR(none[0], 180.0)
  TEST_REAL_SIMILAR(none[3], 3.0)
  TEST_REAL_SIMILAR(none[7], 0.0)
  const PeakFeatures empty = peakFeatures({}, {}, 10.0);
  TEST_REAL_SIMILAR(empty[2], 180.0)
  TEST_REAL_SIMILAR(empty[3], 0.0)
  // A new elution after a gap of more than 30 s with a survey between: the earlier sightings do not count.
  const PeakFeatures e = peakFeatures({0.0, 10000.0, 20000.0, 60000.0, 61000.0}, {5.0, 0.0, 0.0, 2.0, 4.0}, 62000.0);
  TEST_REAL_SIMILAR(e[0], 2.0)    // since the first sighting of the current elution (60 s)
  TEST_REAL_SIMILAR(e[1], 1.0)    // since its apex (61 s)
  TEST_REAL_SIMILAR(e[7], 2.0)    // two sightings
}
END_SECTION

END_TEST
