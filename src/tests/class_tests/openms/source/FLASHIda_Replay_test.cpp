// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MassOffset.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FixedModifications.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include "FLASHIda_EventCapture.h"

#include <nlohmann/json.hpp>

#include <chrono>
#include <filesystem>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// Engine build 1: the replay of a recorded run through the engine's own MS3 path -- the level-3 deconvolution, the
// target fragment's ladder, the match at the run's MS3 mass offset, the offset state and the fallback vote (ADR-0081).
// Each MS3 of the run's mzML is joined to its scan_commands.tsv row by tracking id and taken in file order (the order the
// MS3 returned); its ladder is rebuilt from the row's ms3_proteoform (FLASHIda's own fragment ProForma, the fixed
// modifications in its prefix) with the row's ion type.
//
// In ctest the sections below run on synthetic inputs. With FLASHIDA_REPLAY_SPEC naming a spec file, the last section also
// replays every run it lists and writes, per run, <out_dir>/<name>_ms3.tsv (one row per MS3) and <name>_events.txt (the
// [MASS-OFFSET*] lines). The data stay on this machine; the CI step that sets the variable runs only where they exist.
//
// Spec: {"out_dir": "...", "runs": [{"name": "...", "mzml": "...", "scan_commands": "...", "method": "...",
//        "ms3_start_ppm": -25, "ms3_vote": true, "fixed_offset_ppm": null}]}

using namespace OpenMS;
using json = nlohmann::json;

namespace
{
  /// The sites of FLASHIda's fragment ProForma (FragmentAnalysis::toProForma): "K[+79.9663]" on residue p (1-based),
  /// "(TKI)[+15.9949]" over a range, a leading "<...>" prefix of fixed modifications (read by FixedModifications).
  bool parseFragmentProForma(const std::string& proforma, std::string& seq, std::vector<FragmentAnalysis::PTMSite>& sites)
  {
    seq.clear();
    sites.clear();
    size_t i = 0;
    while (i < proforma.size() && proforma[i] == '<')
    {
      const size_t close = proforma.find('>', i);
      if (close == std::string::npos) return false;
      i = close + 1;
    }
    int open_start = -1;
    for (; i < proforma.size(); ++i)
    {
      const char c = proforma[i];
      if (c == '(')
      {
        open_start = static_cast<int>(seq.size()) + 1;
      }
      else if (c == ')')
      {
        // the range's mass follows
        if (open_start < 0 || i + 1 >= proforma.size() || proforma[i + 1] != '[') return false;
        const size_t close = proforma.find(']', i + 1);
        if (close == std::string::npos) return false;
        FragmentAnalysis::PTMSite s;
        s.start_position = open_start;
        s.end_position = static_cast<int>(seq.size());
        s.position = (s.start_position + s.end_position) / 2;
        s.mass_shift = std::stod(proforma.substr(i + 2, close - i - 2));
        sites.push_back(s);
        open_start = -1;
        i = close;
      }
      else if (c == '[')
      {
        const size_t close = proforma.find(']', i);
        if (close == std::string::npos) return false;
        FragmentAnalysis::PTMSite s;
        s.start_position = s.end_position = s.position = static_cast<int>(seq.size());
        s.mass_shift = std::stod(proforma.substr(i + 1, close - i - 1));
        sites.push_back(s);
        i = close;
      }
      else if (c >= 'A' && c <= 'Z')
      {
        seq.push_back(c);
      }
      else
      {
        return false;
      }
    }
    return !seq.empty() && open_start < 0;
  }

  /// One TSV cell, RFC 4180 unquoted (quoteLogCell quotes an id that starts with a quote).
  std::string unquote(const std::string& cell)
  {
    if (cell.size() >= 2 && cell.front() == '"' && cell.back() == '"')
    {
      std::string out;
      for (size_t i = 1; i + 1 < cell.size(); ++i)
      {
        out.push_back(cell[i]);
        if (cell[i] == '"' && i + 1 < cell.size() - 1 && cell[i + 1] == '"') ++i;
      }
      return out;
    }
    return cell;
  }

  std::vector<std::string> splitTab(const std::string& line)
  {
    std::vector<std::string> out;
    std::string cur;
    for (const char c : line)
    {
      if (c == '\t')
      {
        out.push_back(cur);
        cur.clear();
      }
      else if (c != '\r')
        cur.push_back(c);
    }
    out.push_back(cur);
    return out;
  }

  /// The stage-1 value of a "s0;s1" cell, its first notch ("a,b").
  std::string stage1(const std::string& cell)
  {
    const size_t semi = cell.find(';');
    std::string s = semi == std::string::npos ? cell : cell.substr(semi + 1);
    const size_t comma = s.find(',');
    return comma == std::string::npos ? s : s.substr(0, comma);
  }

  struct Ms3Row
  {
    std::string ion_type;
    int ion_index = 0;
    std::string proforma;
    double mass = 0.0;
    int charge = 0;
  };

  std::map<std::string, Ms3Row> readMs3Rows(const std::string& path)
  {
    std::map<std::string, Ms3Row> rows;
    std::ifstream in(path);
    std::string line;
    if (!std::getline(in, line)) return rows;
    const std::vector<std::string> head = splitTab(line);
    std::map<std::string, size_t> col;
    for (size_t i = 0; i < head.size(); ++i) col[head[i]] = i;
    for (const char* need : {"tracking_id", "ms_level", "mono_mass", "charge", "ion_type", "ion_index", "ms3_proteoform"})
      if (col.find(need) == col.end()) return rows;
    while (std::getline(in, line))
    {
      const std::vector<std::string> c = splitTab(line);
      if (c.size() != head.size() || c[col["ms_level"]] != "3") continue;
      Ms3Row r;
      r.ion_type = c[col["ion_type"]];
      r.ion_index = std::atoi(c[col["ion_index"]].c_str());
      r.proforma = c[col["ms3_proteoform"]];
      r.mass = std::atof(stage1(c[col["mono_mass"]]).c_str());
      r.charge = std::abs(std::atoi(stage1(c[col["charge"]]).c_str()));
      rows[unquote(c[col["tracking_id"]])] = r;
    }
    return rows;
  }

  /// The ladder of a target fragment from its row: the subsequence and sites of its ProForma, the ion types of its class.
  bool ladderOf(const Ms3Row& r, std::vector<MS3FragmentMatcher::TheoreticalMass>& theo, int& length)
  {
    if (r.ion_type.empty() || !MS3FragmentMatcher::isAnnotatableMs3Precursor(r.ion_type[0])) return false;
    std::string seq;
    std::vector<FragmentAnalysis::PTMSite> sites;
    if (!parseFragmentProForma(r.proforma, seq, sites)) return false;
    const FixedModifications fixed = FixedModifications::fromProFormaPrefix(String(r.proforma));
    theo = MS3FragmentMatcher::computeTheoreticalMasses(seq, MS3FragmentMatcher::getMS3IonTypes(r.ion_type[0]), sites, fixed);
    length = static_cast<int>(seq.size());
    return !theo.empty();
  }

  /// One recorded run replayed: returns the number of MS3 replayed.
  int replayRun(const json& run, const std::string& out_dir)
  {
    const std::string name = run.at("name").get<std::string>();
    json method = json::object();
    {
      std::ifstream in(run.at("method").get<std::string>());
      if (in) in >> method;
    }
    // The run's deconvolution settings, as its method authored them; the level-3 cosine is this build's default
    // unless the spec gives one.
    json dc = json::object();
    const json md = method.value("deconvolution", json::object());
    for (const char* k : {"min_charge", "max_charge", "min_mass", "max_mass", "tol"})
      if (md.contains(k)) dc[k] = md[k];
    if (run.contains("ms3_min_cos")) dc["ms3_min_cos"] = run["ms3_min_cos"];
    // No ms_settings: rank_by "none" is the selection that needs none (the replay selects no precursor).
    Config cfg(json {{"deconvolution", dc}, {"precursor_selection", {{"rank_by", "none"}}}}.dump());
    const DoubleList tol = cfg.toleranceList();
    Deconvolution deconv(cfg, tol);
    const double tol3 = tol.size() >= 3 ? tol[2] : 10.0;

    const bool fixed = run.contains("fixed_offset_ppm") && !run["fixed_offset_ppm"].is_null();
    const double fixed_ppm = fixed ? run["fixed_offset_ppm"].get<double>() : 0.0;
    MassOffsets offsets(0.0, run.value("ms3_start_ppm", 0.0), run.value("ms3_vote", false));

    const std::map<std::string, Ms3Row> rows = readMs3Rows(run.at("scan_commands").get<std::string>());
    PeakMap exp;
    MzMLFile().load(run.at("mzml").get<std::string>(), exp);

    std::filesystem::create_directories(out_dir);
    std::ofstream out(out_dir + "/" + name + "_ms3.tsv");
    out << "tracking_id\tscan\tion\tn_masses\tmatched\tbonds\tbond_indices\tconfident\tused_ppm\testimate_after\tdeconv_us\n";
    int replayed = 0;
    EventCapture cap;
    offsets.logStart();
    for (Size si = 0; si < exp.size(); ++si)
    {
      const MSSpectrum& sp = exp[si];
      if (sp.getMSLevel() != 3) continue;
      const std::string desc = sp.metaValueExists("scan description") ? sp.getMetaValue("scan description").toString() : "";
      if (desc.size() < 3) continue;
      const std::string tid = desc.substr(0, 3);
      const auto row = rows.find(tid);
      if (row == rows.end()) continue;
      std::vector<MS3FragmentMatcher::TheoreticalMass> theo;
      int length = 0;
      if (!ladderOf(row->second, theo, length)) continue;

      std::vector<double> mzs, ints;
      for (const auto& p : sp)
      {
        mzs.push_back(p.getMZ());
        ints.push_back(p.getIntensity());
      }
      const auto t0 = std::chrono::steady_clock::now();
      deconv.deconvolveMSn(mzs.data(), ints.data(), static_cast<int>(mzs.size()), sp.getRT(), row->second.mass,
                           row->second.charge, 3);
      const auto t1 = std::chrono::steady_clock::now();
      const DeconvolvedSpectrum spec = deconv.hasStoredMS2() ? deconv.storedMS2() : DeconvolvedSpectrum(0);

      const double used = fixed ? fixed_ppm : offsets.ms3.estimate();
      std::vector<MS3FragmentMatcher::MatchDetail> details;
      const MS3FragmentMatcher::OffsetEvidence ev = MS3FragmentMatcher::offsetEvidence(spec, theo, used, tol3, length, &details);
      if (!fixed) offsets.onMs3(tid, used, ev.n_masses, ev.matched, ev.bonds, ev.raw_errors, ev.wide_errors);

      std::set<int> bond_set;
      for (const auto& d : details)
      {
        if (d.ion_type.empty()) continue;
        const char t = d.ion_type[0];
        const int b = (t == 'a' || t == 'b' || t == 'c') ? d.position : (length - d.position);
        if (b >= 1 && b <= length - 1) bond_set.insert(b);
      }
      std::string bonds;
      for (const int b : bond_set) bonds += (bonds.empty() ? "" : ";") + std::to_string(b);
      out << tid << "\t" << sp.getNativeID() << "\t" << row->second.ion_type << row->second.ion_index << "\t" << ev.n_masses
          << "\t" << ev.matched << "\t" << ev.bonds << "\t" << bonds << "\t" << (ev.bonds >= kMassOffsetConfidentMs3Bonds ? 1 : 0)
          << "\t" << std::fixed << std::setprecision(4) << used << "\t" << (fixed ? fixed_ppm : offsets.ms3.estimate()) << "\t"
          << std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() << "\n";
      ++replayed;
    }
    std::ofstream(out_dir + "/" + name + "_events.txt") << cap.str();
    return replayed;
  }
} // namespace

START_TEST(FLASHIda_Replay, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION(parses_flashidas_fragment_proforma)
{
  std::string seq;
  std::vector<FragmentAnalysis::PTMSite> sites;
  TEST_EQUAL(parseFragmentProForma("<[+329.2266]@C>ACK[+42.0106]DE(FGH)[-17.0265]IK", seq, sites), true)
  TEST_EQUAL(seq, "ACKDEFGHIK")
  ABORT_IF(sites.size() != 2)
  TEST_EQUAL(sites[0].start_position, 3)   // on K, residue 3
  TEST_EQUAL(sites[0].end_position, 3)
  TEST_REAL_SIMILAR(sites[0].mass_shift, 42.0106)
  TEST_EQUAL(sites[1].start_position, 6)   // F..H, residues 6-8
  TEST_EQUAL(sites[1].end_position, 8)
  TEST_REAL_SIMILAR(sites[1].mass_shift, -17.0265)
  TEST_EQUAL(parseFragmentProForma("PEPTIDE", seq, sites), true)
  TEST_EQUAL(sites.empty(), true)
  TEST_EQUAL(parseFragmentProForma("PEP(TIDE", seq, sites), false)   // an open range
  TEST_EQUAL(parseFragmentProForma("", seq, sites), false)
}
END_SECTION

START_SECTION(reads_the_stage_1_cells_and_quoted_ids)
{
  TEST_EQUAL(stage1("12345.6;5183.55"), "5183.55")
  TEST_EQUAL(stage1("12;7,6"), "7")
  TEST_EQUAL(unquote("\"!\"\"a\""), "!\"a")
  TEST_EQUAL(unquote("!!a"), "!!a")
}
END_SECTION

START_SECTION(a_synthetic_ms3_replays_at_the_run_offset)
{
  // The ladder of b10 of ACDEFGHIKLMNPQRSTVWY, every entry observed 25 ppm low: at the run offset -25 they match.
  Ms3Row r;
  r.ion_type = "b";
  r.ion_index = 10;
  r.proforma = "ACDEFGHIKL";
  std::vector<MS3FragmentMatcher::TheoreticalMass> theo;
  int length = 0;
  TEST_EQUAL(ladderOf(r, theo, length), true)
  TEST_EQUAL(length, 10)
  DeconvolvedSpectrum spec(0);
  for (size_t i = 0; i < theo.size() && i < 6; ++i)
  {
    PeakGroup pg(1, 1, true);
    pg.setMonoisotopicMass(theo[i].mass * (1.0 - 25e-6));
    spec.push_back(pg);
  }
  const auto at_start = MS3FragmentMatcher::offsetEvidence(spec, theo, -25.0, 10.0, length);
  TEST_EQUAL(at_start.matched >= 3, true)
  const auto at_zero = MS3FragmentMatcher::offsetEvidence(spec, theo, 0.0, 10.0, length);
  TEST_EQUAL(at_zero.matched, 0)

  // The replay's engine is built from a config holding only the run's deconvolution section, as replayRun builds it.
  json dc = {{"min_charge", 1}, {"max_charge", 50}, {"tol", {10, 10, 10}}};
  Config cfg(json {{"deconvolution", dc}, {"precursor_selection", {{"rank_by", "none"}}}}.dump());
  TEST_EQUAL(cfg.toleranceList().size(), 3)
  TEST_REAL_SIMILAR(cfg.deconvolution().ms3_min_cos, 0.5)
  Deconvolution deconv(cfg, cfg.toleranceList());
  TEST_EQUAL(deconv.hasStoredMS2(), false)
}
END_SECTION

START_SECTION(replays_the_recorded_runs_of_the_spec_when_asked)
{
  // ctest leaves FLASHIDA_REPLAY_SPEC unset: nothing to replay, the section passes. The CI replay step sets it.
  const char* spec_path = std::getenv("FLASHIDA_REPLAY_SPEC");
  if (spec_path == nullptr || *spec_path == '\0')
  {
    TEST_EQUAL(true, true)
  }
  else
  {
    std::ifstream in(spec_path);
    ABORT_IF(!in)
    json spec;
    in >> spec;
    const std::string out_dir = spec.value("out_dir", std::string("testlogs/replay"));
    for (const json& run : spec.at("runs"))
    {
      const int n = replayRun(run, out_dir);
      std::cout << "[REPLAY] " << run.at("name").get<std::string>() << " ms3=" << n << std::endl;
      TEST_EQUAL(n > 0, true)
    }
  }
}
END_SECTION

END_TEST
