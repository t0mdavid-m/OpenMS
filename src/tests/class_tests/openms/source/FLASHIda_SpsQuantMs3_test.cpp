// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

// ADR-0080: the SPS quantification MS3 on hardware data. The fixture (test-data-private/sps_quant_ms3/, in the PRIVATE
// parent repo; made by test-scripts/make_sps_quant_fixture.py) holds every identification scan of the 2026-11-10
// E. coli + Caco-2 iodoTMT mix (fi-ms3-mix_261110) that bought a quantification MS3, 371 cases, each with its decision
// survey and the protein the live search named. Every case is replayed through runInterleaved twice, at sps_windows
// 1 and 10: the survey, strict inclusion on the commanded mass, the identification scan as the MS2's return, the
// protein inline as the reference, the synthetic reporter spectrum as the quantification MS3's return.
//
//   replay_fi_ms3_mix         (T9) the anchor at 10 is the anchor at 1 (sps_windows 1 is ADR-0064, on hardware data),
//                             and every window is a matched labelled fragment of the identification with window SNR > 1,
//                             logged once in quant_ms3.tsv; the windows are disjoint, at most 10, no peak twice.
//   one_real_scan_pool_table  (T8) the case with the most SNR-passing pool peaks: its whole quant_ms3.tsv table at 10,
//                             printed for the owner's review (ADR-0080 d.7), with plausibility asserts only.
//
// The live anchor (the run's own 'Q' command) is printed beside the replay's, never asserted: the live build predates
// ADR-0074..0079, so a difference there is not the SPS fill's (owner, 2026-11-19).

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanRole.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include "FLASHIda_TestHelpers.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace OpenMS;

namespace
{
  const std::string kDir = "../../test-data-private/sps_quant_ms3/";
  const std::string kReporterMs3 = "../../FlashIDA/test-data/spectra/ms2_quant_tmt.txt";   // synthetic, six TMT peaks

  struct Case
  {
    int id = 0;
    size_t survey_index = 0, r_index = 0;
    std::string r_id, mono_mass, protein, sequence, live_anchor;
  };

  std::vector<Case> loadCases()
  {
    std::vector<Case> out;
    const TSVFile t = TSVFile::parse(kDir + "cases.tsv");
    for (const auto& r : t.rows)
    {
      Case c;
      c.id = std::stoi(r[t.colIndex("case")]);
      c.survey_index = std::stoul(r[t.colIndex("survey_index")]);
      c.r_index = std::stoul(r[t.colIndex("r_index")]);
      c.r_id = r[t.colIndex("r_id")];
      c.mono_mass = r[t.colIndex("mono_mass")];
      c.protein = r[t.colIndex("protein")];
      c.sequence = r[t.colIndex("sequence")];
      // As a 'Q' description ends: "@<charge><ion type><index>".
      c.live_anchor = "@" + r[t.colIndex("live_charge")] + r[t.colIndex("live_ion_type")] + r[t.colIndex("live_ion_index")];
      out.push_back(c);
    }
    return out;
  }

  ScanData toScanData(const MSSpectrum& s)
  {
    ScanData d;
    for (const auto& p : s)
    {
      d.mzs.push_back(p.getMZ());
      d.ints.push_back(p.getIntensity());
    }
    d.rt = s.getRT();   // seconds, as every C++ drive feeds (ADR-0051)
    const std::string nid = s.getNativeID();
    const size_t at = nid.find("scan=");
    d.scan_id = at == std::string::npos ? std::string() : nid.substr(at + 5);
    return d;
  }

  /// The run's own method (deconvolution, the identification scan, the search) with scan ms3, strict inclusion on the
  /// case's commanded mass and its protein inline.
  std::string caseJson(const std::string& dir, const Case& c, int sps)
  {
    const std::string inclusion = dir + "/inclusion.txt";
    {
      std::ofstream f(inclusion);
      f << "mass\tcharge\trt_start\trt_end\tpriority\n" << c.mono_mass << "\t-1\t0\t1000000\t1\n";
    }
    nlohmann::json j = nlohmann::json::parse(buildJsonWithLogDir(dir));
    j["deconvolution"]["min_mass"] = 2000;
    j["precursor_selection"]["max_precursors"] = 1;
    j["precursor_selection"]["targeting"] = "inclusion";
    j["precursor_selection"]["strict_inclusion"] = true;
    j["files"]["inclusion_list"] = inclusion;
    j["ms_settings"]["ms2"] = {{"analyzer", "Orbitrap"}, {"activation", "CID"}, {"collision_energy", 35},
                               {"first_mass", 200}, {"last_mass", 2000}, {"resolution", 120000}};
    j["ms_settings"]["ms3_quant"] = {{"analyzer", "Orbitrap"}, {"activation", "HCD"}, {"collision_energy", 80},
                                     {"first_mass", 120}, {"last_mass", 140}, {"resolution", 30000}};
    j["identification"] = {
      {"reference", {{"sequence", c.sequence}}},
      {"flashtnt", {{"min_length", 3}, {"max_length", 8}, {"max_blind_mod_count", 5}, {"max_mod_mass", 500}}},
      {"threshold", 50}};
    j["quantification"] = {{"enabled", true}, {"labelling", "tmt6plex"}, {"scan", "ms3"}, {"ms3_targets", 1},
                           {"sps_windows", sps},
                           {"conditions", {{{"name", "glucose"}, {"channels", {"126", "127", "128"}}},
                                           {{"name", "acetate"}, {"channels", {"129", "130", "131"}}}}}};
    return j.dump();
  }

  struct Replay
  {
    std::vector<ScanCommand> quant;               ///< the quantification MS3s, in dequeue order
    TSVFile pool;                                 ///< quant_ms3.tsv
    std::set<std::pair<std::string, int>> matched;   ///< (ion type, index) of every MS2 identification row's fragments
  };

  Replay replay(const Case& c, const ScanData& survey, const ScanData& ms2, const ScanData& reporter, int sps)
  {
    Replay out;
    const std::string dir = freshLogDir("sps_replay/c" + std::to_string(c.id) + "_n" + std::to_string(sps));
    const std::string cfg = caseJson(dir, c, sps);
    const std::map<std::string, std::vector<ScanData>> no_char_fixtures;
    AcqResult acq;
    {
      FLASHIda ida(const_cast<char*>(cfg.c_str()));
      acq = runInterleaved(&ida, {survey}, {ms2}, &no_char_fixtures, 600, false, 0.0, nullptr, nullptr, &reporter);
    }
    for (const auto& cmd : acq.ms3_cmds)
      if (roleOf(cmd) == ScanRole::Quantification) out.quant.push_back(cmd);
    out.pool = TSVFile::parse(dir + "/quant_ms3.tsv");
    const TSVFile idf = TSVFile::parse(dir + "/identification.tsv");
    for (const auto& row : idf.rows)
    {
      if (row[idf.colIndex("ms_level")] != "2") continue;
      std::stringstream ss(row[idf.colIndex("ms2_fragments")]);
      std::string ion;
      while (std::getline(ss, ion, ';'))
      {
        size_t k = 0;
        while (k < ion.size() && std::isalpha(static_cast<unsigned char>(ion[k]))) ++k;
        if (k > 0 && k < ion.size()) out.matched.emplace(ion.substr(0, k), std::stoi(ion.substr(k)));
      }
    }
    return out;
  }

  std::string anchorOf(const ScanCommand& c) { return std::string(c.scan_description).substr(3); }

  /// The pool rows that are windows of the MS3 @p id, in window order.
  std::vector<const std::vector<std::string>*> windowRows(const TSVFile& pool, const std::string& id)
  {
    std::vector<const std::vector<std::string>*> rows;
    for (const auto& r : pool.rows)
      if (r[pool.colIndex("tracking_id")] == id && std::stoi(r[pool.colIndex("window")]) >= 1) rows.push_back(&r);
    std::sort(rows.begin(), rows.end(), [&pool](const auto* a, const auto* b) {
      return std::stoi((*a)[pool.colIndex("window")]) < std::stoi((*b)[pool.colIndex("window")]);
    });
    return rows;
  }

  int chosen_case = -1;   // T9 picks the case T8 prints: the most SNR-passing pool peaks, the lowest case on a tie
}

START_TEST(FLASHIda_SpsQuantMs3, "$Id$")

const std::vector<Case> cases = loadCases();
MSExperiment slice;
MzMLFile().load(kDir + "slice.mzML", slice);
const std::vector<ScanData> reporter = loadTsvScans(kReporterMs3);

START_SECTION(replay_fi_ms3_mix)
{
  TEST_EQUAL(cases.size(), 371)
  TEST_EQUAL(slice.size(), 725)
  ABORT_IF(cases.size() != 371 || slice.size() != 725 || reporter.empty())
  const auto t0 = std::chrono::steady_clock::now();
  int not_bought = 0, live_same = 0, best_passing = -1;
  std::map<int, int> window_hist;
  for (const Case& c : cases)
  {
    const ScanData survey = toScanData(slice[c.survey_index]);
    const ScanData ms2 = toScanData(slice[c.r_index]);
    const Replay r1 = replay(c, survey, ms2, reporter[0], 1);
    const Replay r10 = replay(c, survey, ms2, reporter[0], 10);
    TEST_EQUAL(r10.quant.size(), r1.quant.size())
    if (r1.quant.empty() || r10.quant.empty())
    {
      ++not_bought;
      continue;
    }
    const ScanCommand& a = r1.quant[0];
    const ScanCommand& s = r10.quant[0];
    // sps_windows 1 is ADR-0064's MS3: one window. At 10 the anchor is the same.
    TEST_EQUAL(notchesForStage(a, 1).second, 0)
    TEST_STRING_EQUAL(anchorOf(s), anchorOf(a))
    TEST_REAL_SIMILAR(s.stages[1].precursor_mz, a.stages[1].precursor_mz)
    const std::string aa = anchorOf(a);
    if (aa.size() >= c.live_anchor.size() && aa.compare(aa.size() - c.live_anchor.size(), std::string::npos, c.live_anchor) == 0) ++live_same;

    const auto notches = notchesForStage(s, 1);
    const int windows = 1 + notches.second;
    ++window_hist[windows];
    TEST_EQUAL(windows <= 10, true)
    // Every window is one quant_ms3.tsv row of this MS3: a matched fragment of the identification, with a cysteine,
    // window SNR > 1; no (ion, charge) twice; the windows disjoint.
    const auto rows = windowRows(r10.pool, ScanCommandQueue::encode(s.scan_id));
    TEST_EQUAL(static_cast<int>(rows.size()), windows)
    std::set<std::tuple<std::string, int, int>> seen;
    std::vector<std::pair<double, double>> w;
    for (const auto* row : rows)
    {
      const std::string ion = (*row)[r10.pool.colIndex("ion_type")];
      const int index = std::stoi((*row)[r10.pool.colIndex("ion_index")]);
      const int z = std::stoi((*row)[r10.pool.colIndex("charge")]);
      TEST_EQUAL(r10.matched.count({ion, index}), 1)
      TEST_EQUAL((*row)[r10.pool.colIndex("cys_positions")].empty(), false)
      TEST_EQUAL(std::stod((*row)[r10.pool.colIndex("window_snr")]) > 1.0, true)
      TEST_EQUAL(seen.insert({ion, index, z}).second, true)
      const double mz = std::stod((*row)[r10.pool.colIndex("mz")]);
      const double width = std::stod((*row)[r10.pool.colIndex("width")]);
      w.emplace_back(mz - width / 2.0, mz + width / 2.0);
    }
    bool disjoint = true;
    for (size_t i = 0; i < w.size(); ++i)
      for (size_t k = i + 1; k < w.size(); ++k) disjoint = disjoint && !(w[i].first < w[k].second && w[k].first < w[i].second);
    TEST_EQUAL(disjoint, true)

    int passing = 0;
    for (const auto& row : r10.pool.rows)
      if (row[r10.pool.colIndex("status")] != "snr" && row[r10.pool.colIndex("status")] != "yield") ++passing;
    if (passing > best_passing)
    {
      best_passing = passing;
      chosen_case = c.id;
    }
  }
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::cout << "[SPS-REPLAY] cases=" << cases.size() << " not_bought=" << not_bought
            << " anchor_equals_live=" << live_same << " seconds=" << secs << "\n[SPS-REPLAY] windows per MS3:";
  for (const auto& [n, k] : window_hist) std::cout << " " << n << ":" << k;
  std::cout << "\n[SPS-REPLAY] chosen_case=" << chosen_case << " passing_pool_peaks=" << best_passing << std::endl;
  TEST_EQUAL(not_bought < 371, true)   // the replay measures something
}
END_SECTION

START_SECTION(one_real_scan_pool_table)
{
  ABORT_IF(chosen_case < 0)
  const Case& c = cases[static_cast<size_t>(chosen_case)];
  const ScanData survey = toScanData(slice[c.survey_index]);
  const ScanData ms2 = toScanData(slice[c.r_index]);
  const Replay r1 = replay(c, survey, ms2, reporter[0], 1);
  const Replay r10 = replay(c, survey, ms2, reporter[0], 10);
  ABORT_IF(r1.quant.empty() || r10.quant.empty())
  std::cout << "[SPS-TABLE] case=" << c.id << " r_id=" << c.r_id << " protein=" << c.protein << " mono_mass=" << c.mono_mass
            << " live_anchor=" << c.live_anchor << "\n[SPS-TABLE]";
  for (const auto& h : r10.pool.headers) std::cout << "\t" << h;
  std::cout << "\n";
  for (const auto& row : r10.pool.rows)
  {
    std::cout << "[SPS-TABLE]";
    for (const auto& v : row) std::cout << "\t" << v;
    std::cout << "\n";
  }
  std::cout.flush();
  // Plausibility only until the owner has reviewed this table (ADR-0080 d.7); the reviewed list is pinned after.
  TEST_STRING_EQUAL(anchorOf(r10.quant[0]), anchorOf(r1.quant[0]))
  TEST_EQUAL(windowRows(r10.pool, ScanCommandQueue::encode(r10.quant[0].scan_id)).size() >= 2, true)
}
END_SECTION

END_TEST
