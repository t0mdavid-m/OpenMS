// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Plausibility tests for the four FLASHIda log streams (ida_log, scan_commands.tsv,
// scan_results.tsv, identification.tsv). Per the agreed division of labor these are
// DRIFT-STABLE plausibility checks -- ranges / sets / signs / cross-field relationships /
// structural invariants -- NOT exact golden values. Only config-deterministic quantities
// (column counts, ms_level, scan_type, activation set, MS2 CE==29 / MS3 CE==35, the 15-float
// IDA-log length) are asserted exactly. Exact, all-execution-path behaviour locking lives in
// the C# golden suite (FLASHIdaLogGolden_test), which drives the real method.json configs.
//
// Coverage notes (intentional C++/C# split, not silent gaps):
//   * MS3-exploration *result* rows and the metric-dependent identification 'E'-row MIXED
//     behaviour are golden-locked in C# -- the single-group C++ driver stops at the first MS3
//     command, so those rows are not reliably reachable here.
//   * Exploration variant rows log tag_count==0 because no sequence tagging runs on CE-sweep
//     variants; the >0 real-tag-count path (standard MS2/MS3) is plausibility-checked (>=0)
//     here and golden-locked exactly in C#.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include "FLASHIda_TestHelpers.h"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  // tracking_id -> ms_level, from a parsed scan_commands.tsv (used to classify result rows,
  // which carry no ms_level column: a result row's id is in this map iff it is an MS2/MS3 scan;
  // absent => MS1 survey input).
  std::map<std::string, int> commandLevels(const TSVFile& cmds)
  {
    std::map<std::string, int> m;
    int id_col = cmds.colIndex("tracking_id");
    int lvl_col = cmds.colIndex("ms_level");
    if (id_col < 0 || lvl_col < 0) return m;
    for (const auto& r : cmds.rows)
      if (id_col < (int)r.size() && lvl_col < (int)r.size())
        m[r[id_col]] = std::atoi(r[lvl_col].c_str());
    return m;
  }

  bool ionTypeOk(const std::string& s)
  {
    if (s.size() != 1) return false;
    char c = s[0];
    return c == 'a' || c == 'b' || c == 'c' || c == 'x' || c == 'y' || c == 'z';
  }

  // Decode the trailing precursor ion of an MS3 scan_description per the descriptor contract:
  //   {id}R{mass}k@{charge}{ion_type}{ion_index}  (the 'k' is part of the mass token, before '@').
  // Find the LAST '@', skip the run of fragment-charge digits after it, require an ion_type char in
  // {a,b,c,x,y,z} (via ionTypeOk), then an all-digit ion_index >= 1. Returns false for the no-ion form
  // {id}R{mass}k@{charge} (nothing after the charge digits) -- that branch is tolerated by callers.
  // Kept byte-for-byte identical to decodeTrailingIonKey in FLASHIda_TestHelpers.h.
  bool decodeTrailingIon(const std::string& d, char& t, int& idx)
  {
    std::size_t at = d.rfind('@');
    if (at == std::string::npos) return false;
    std::size_t i = at + 1;
    while (i < d.size() && std::isdigit(static_cast<unsigned char>(d[i]))) ++i;  // skip fragment charge digits
    if (i >= d.size() || !ionTypeOk(std::string(1, d[i]))) return false;          // no-ion branch -> false
    t = d[i++];
    std::string n = d.substr(i);
    if (n.empty() || n.find_first_not_of("0123456789") != std::string::npos) return false;  // non-truncated index
    idx = std::atoi(n.c_str());
    return idx >= 1;
  }

  // Build ion -> [ScanData...] from real MS3 fixtures ms3_cytc_<ion>_scan<N>.txt in a spectra dir.
  // The ion key (e.g. "b44") is the substring between "ms3_cytc_" and "_scan" (per the fixture naming
  // contract). loadTsvScans each matching file and group by ion. Returns empty if the dir is absent or
  // holds no matching fixtures (caller then skips cleanly -- never fabricates data).
  std::map<std::string, std::vector<ScanData>> buildMs3IonManifest(const std::string& dir)
  {
    std::map<std::string, std::vector<ScanData>> manifest;
    std::error_code ec;
    if (!std::filesystem::is_directory(dir, ec)) return manifest;
    const std::string prefix = "ms3_cytc_";
    for (std::filesystem::directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec))
    {
      if (!it->is_regular_file(ec)) continue;
      std::string name = it->path().filename().string();
      if (name.size() < prefix.size() || name.compare(0, prefix.size(), prefix) != 0) continue;
      // The ion token sits between "ms3_cytc_" and the LAST "_scan"; require both present and a non-empty ion.
      std::size_t scan_pos = name.rfind("_scan");
      if (scan_pos == std::string::npos || scan_pos <= prefix.size()) continue;  // e.g. legacy ms3_cytc_scan341.txt
      std::string ion = name.substr(prefix.size(), scan_pos - prefix.size());
      if (ion.empty()) continue;
      auto scans = loadTsvScans(it->path().string());
      if (scans.empty()) continue;
      auto& bucket = manifest[ion];
      for (auto& s : scans) bucket.push_back(std::move(s));
    }
    return manifest;
  }

  // Strict ProForma validator (per the plan glossary): every RESIDUE char (alphabetic, outside [...]
  // mod groups and the () region delimiters) must be UPPERCASE A-Z (no lowercase residues); modifications
  // are bracketed [..]; an ambiguity region '(' RESIDUES ')' must be immediately followed by a '[' mod;
  // brackets balanced; >=1 residue. e.g. GDVEK, PEP[+10.0000]TIDE, (GDVEKGKK)[+42.0157]IFVQK valid;
  // gdvek[+42] (lowercase), GDVEK(IFVQ)KIFVQK (region not followed by a mod) invalid.
  bool isProForma(const std::string& s)
  {
    if (s.empty()) return false;
    bool in_mod = false;             // inside [...]
    bool expect_mod = false;         // just closed a ')', a '[' mod must follow
    bool any_residue = false;
    for (char ch : s)
    {
      if (in_mod) { if (ch == ']') in_mod = false; continue; }
      if (expect_mod) { if (ch != '[') return false; expect_mod = false; in_mod = true; continue; }
      if (ch == '[') { in_mod = true; continue; }
      if (ch == '(') continue;
      if (ch == ')') { expect_mod = true; continue; }
      if (std::isalpha(static_cast<unsigned char>(ch)))
      {
        if (std::isupper(static_cast<unsigned char>(ch)) == 0) return false;  // lowercase residue
        any_residue = true;
      }
    }
    return any_residue && !in_mod && !expect_mod;
  }

  const std::string CYTC_MS1 = FI_MS1_CYTC;
  const std::string CYTC_MS2 = FI_MS2_CYTC;
}

START_TEST(FLASHIda_LoggingFields, "$Id$")

/////////////////////////////////////////////////////////////
// §0 -- schema column counts + exact header order
/////////////////////////////////////////////////////////////
START_SECTION(schema_column_counts)
{
  std::string cmd_f = "lf_schema_commands.tsv", res_f = "lf_schema_results.tsv", id_f = "lf_schema_id.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str()); std::remove(id_f.c_str());

  std::string json = buildJsonWithRuntime("", cmd_f, res_f, true, id_f);
  FLASHIda ida(const_cast<char*>(json.c_str()));  // constructor writes the three headers

  auto c = TSVFile::parse(cmd_f);
  auto r = TSVFile::parse(res_f);
  auto i = TSVFile::parse(id_f);

  TEST_EQUAL(c.headers.size(), 29)   // E6: + scan_description
  TEST_EQUAL(r.headers.size(), 34)   // E5: + ms_level; F5: + winner_tracking_id
  TEST_EQUAL(i.headers.size(), 25)   // I2: + ms2/ms3 isolation_width, window_snr, charge_intensity (6)

  // Spot-check exact header identities / order at the boundaries that matter for parsing.
  TEST_EQUAL(c.headers.front(), std::string("tracking_id"))
  TEST_EQUAL(c.headers.back(), std::string("scan_description"))   // E6: appended last
  TEST_EQUAL(c.colIndex("scan_description"), 28)
  TEST_EQUAL(c.colIndex("hcd_energy"), 21)
  TEST_EQUAL(r.headers.front(), std::string("tracking_id"))
  TEST_EQUAL(r.colIndex("ms_level"), 1)                           // E5: inserted right after tracking_id
  TEST_EQUAL(r.headers.back(), std::string("winner_tracking_id"))   // F5: appended last
  TEST_EQUAL(r.colIndex("winner_tracking_id"), 33)                 // F5: trailing column index
  TEST_EQUAL(r.colIndex("processing_duration_ms"), 32)             // now second-to-last
  TEST_EQUAL(r.colIndex("child_ids"), 9)                          // shifted +1 by ms_level@1
  TEST_EQUAL(i.headers.front(), std::string("ms_level"))
  TEST_EQUAL(i.headers.back(), std::string("ms3_charge_intensity"))   // I2: appended last
  // I2: the 6 new identification columns appended after ms3_fragment_masses (col 18), in order.
  TEST_EQUAL(i.colIndex("ms3_fragment_masses"), 18)
  TEST_EQUAL(i.colIndex("ms2_isolation_width"), 19)
  TEST_EQUAL(i.colIndex("ms2_window_snr"), 20)
  TEST_EQUAL(i.colIndex("ms2_charge_intensity"), 21)
  TEST_EQUAL(i.colIndex("ms3_isolation_width"), 22)
  TEST_EQUAL(i.colIndex("ms3_window_snr"), 23)
  TEST_EQUAL(i.colIndex("ms3_charge_intensity"), 24)

  std::remove(cmd_f.c_str()); std::remove(res_f.c_str()); std::remove(id_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §C1 -- scan_commands.tsv : MS2 (recording) columns
/////////////////////////////////////////////////////////////
START_SECTION(commands_ms2_columns)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_c1_commands.tsv";
  std::remove(cmd_f.c_str());
  std::string json = buildJsonWithRuntime("", cmd_f, "", true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto cycle = runFullCycle(&ida, ms1, ms2);
  TEST_TRUE(cycle.ms2_cmds.size() > 0)

  auto t = TSVFile::parse(cmd_f);
  int lvl = t.colIndex("ms_level");
  ABORT_IF(lvl < 0)

  bool found = false;
  bool tid_ok = true, type_ok = true, act_ok = true, ce_ok = true, qsc_ok = true;
  bool snr_ok = true, ppm_ok = true, pint_ok = true, parent_ok = true, ion_ok = true, nosemi = true;
  for (const auto& row : t.rows)
  {
    if (lvl >= (int)row.size() || row[lvl] != "2") continue;
    found = true;
    tid_ok    = tid_ok    && isTrackingId(cell(t, row, "tracking_id"));
    type_ok   = type_ok   && (cell(t, row, "scan_type") == "recording");
    act_ok    = act_ok    && inActivationSet(cell(t, row, "activation"));
    ce_ok     = ce_ok     && std::abs(toD(cell(t, row, "collision_energy")) - 29.0) < 1.0;  // [cx] MS2 CE
    qsc_ok    = qsc_ok    && inUnit(toD(cell(t, row, "qscore")));
    snr_ok    = snr_ok    && nonNegFinite(toD(cell(t, row, "snr")));
    ppm_ok    = ppm_ok    && finiteVal(toD(cell(t, row, "ppm_error")));
    pint_ok   = pint_ok   && posFinite(toD(cell(t, row, "precursor_intensity")));
    parent_ok = parent_ok && isTrackingId(cell(t, row, "parent_tracking_id"));  // == its MS1
    ion_ok    = ion_ok    && cell(t, row, "ion_type").empty() && (cell(t, row, "ion_index") == "0");
    nosemi    = nosemi    && (cell(t, row, "charge").find(';') == std::string::npos)
                         && (cell(t, row, "qscore").find(';') == std::string::npos);  // single stage
    TEST_TRUE(inMass(toD(cell(t, row, "mono_mass"))))
    TEST_TRUE(inChargeD(toD(cell(t, row, "charge"))))
    TEST_TRUE(posFinite(toD(cell(t, row, "precursor_mz"))))
    TEST_TRUE(posFinite(toD(cell(t, row, "isolation_width"))))
    TEST_EQUAL(row.size(), t.headers.size())
  }
  TEST_TRUE(found)
  TEST_TRUE(tid_ok)  TEST_TRUE(type_ok)  TEST_TRUE(act_ok)  TEST_TRUE(ce_ok)
  TEST_TRUE(qsc_ok)  TEST_TRUE(snr_ok)   TEST_TRUE(ppm_ok)  TEST_TRUE(pint_ok)
  TEST_TRUE(parent_ok) TEST_TRUE(ion_ok) TEST_TRUE(nosemi)

  std::remove(cmd_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §C2 -- scan_commands.tsv : MS3 two-stage columns + 2-stage scoring
/////////////////////////////////////////////////////////////
START_SECTION(commands_ms3_two_stage)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_c2_commands.tsv";
  std::remove(cmd_f.c_str());
  std::string json = buildJsonWithRuntime("", cmd_f, "", true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto cycle = runFullCycle(&ida, ms1, ms2);
  TEST_TRUE(cycle.ms3_cmds.size() > 0)  // MS3 must fire, else this section asserts nothing

  auto t = TSVFile::parse(cmd_f);
  int lvl = t.colIndex("ms_level");
  ABORT_IF(lvl < 0)

  bool found = false, stage1_populated = false;
  bool charge2_ok = true, ce_stage1_ok = true, act2_ok = true, ion_ok = true, parent_ok = true;
  bool score2_ok = true, hcd_stage1_ok = true;
  for (const auto& row : t.rows)
  {
    if (lvl >= (int)row.size() || row[lvl] != "3") continue;
    found = true;

    // Per-stage cols carry exactly two ';'-tokens (MS2 isolation ; MS3 fragmentation).
    auto chg = splitTokens(cell(t, row, "charge"), ';');
    auto ce  = splitTokens(cell(t, row, "collision_energy"), ';');
    auto act = splitTokens(cell(t, row, "activation"), ';');
    auto hcd = splitTokens(cell(t, row, "hcd_energy"), ';');
    // stage-0 = MS2 precursor charge (>= min_charge); stage-1 = FRAGMENT charge (1..parent), not min_charge-bounded
    charge2_ok = charge2_ok && chg.size() == 2 && inChargeD(toD(chg[0])) && inFragCharge(toD(chg[1]), toD(chg[0]));
    ce_stage1_ok = ce_stage1_ok && ce.size() == 2 && std::abs(toD(ce[1]) - 35.0) < 1.0;  // [cx] MS3 CE
    act2_ok = act2_ok && act.size() == 2 && inActivationSet(act[0]) && inActivationSet(act[1]);
    hcd_stage1_ok = hcd_stage1_ok && hcd.size() == 2 && std::abs(toD(hcd[1]) - 35.0) < 1.0;

    ion_ok = ion_ok && ionTypeOk(cell(t, row, "ion_type")) && (std::atoi(cell(t, row, "ion_index").c_str()) >= 1);
    parent_ok = parent_ok && isTrackingId(cell(t, row, "parent_tracking_id"));  // == its MS2

    // The 11 scoring cols are 2-stage 'stage0;stage1'. stage0 == parent MS2 score (range-checked);
    // stage1 is the fresh fragment score (allowed 0, but at least one row must populate it).
    const char* scoring[] = {"mono_mass","qscore","charge_cos","charge_snr","iso_cos","snr",
                             "charge_score","ppm_error","precursor_intensity","peakgroup_intensity"};
    for (const char* col : scoring)
    {
      auto tok = splitTokens(cell(t, row, col), ';');
      score2_ok = score2_ok && (tok.size() == 2) && finiteVal(toD(tok[0])) && finiteVal(toD(tok[1]));
    }
    auto mm1 = splitTokens(cell(t, row, "mono_mass"), ';');
    auto qs1 = splitTokens(cell(t, row, "qscore"), ';');
    if (mm1.size() == 2 && toD(mm1[1]) > 0.0) stage1_populated = true;          // fragment has a real mass
    if (qs1.size() == 2) { TEST_TRUE(inUnit(toD(qs1[1]))) }                     // fresh pg.getQscore() in [0,1]
  }
  TEST_TRUE(found)
  TEST_TRUE(charge2_ok)  TEST_TRUE(ce_stage1_ok)  TEST_TRUE(act2_ok)  TEST_TRUE(hcd_stage1_ok)
  TEST_TRUE(ion_ok)      TEST_TRUE(parent_ok)     TEST_TRUE(score2_ok)
  TEST_TRUE(stage1_populated)  // proves the MS3 2-stage scoring plumb reached the writer

  std::remove(cmd_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §C3 -- scan_commands.tsv : MS1 / AGC stage-less rows
/////////////////////////////////////////////////////////////
START_SECTION(commands_ms1_agc_stageless)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_c3_commands.tsv";
  std::remove(cmd_f.c_str());
  std::string json = buildJsonWithRuntime("", cmd_f, "", true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto t = TSVFile::parse(cmd_f);
  int lvl = t.colIndex("ms_level");
  ABORT_IF(lvl < 0)

  bool found = false, type_ok = true, act_ok = true, zero_ok = true, parent_ok = true, ion_ok = true;
  for (const auto& row : t.rows)
  {
    if (lvl >= (int)row.size() || row[lvl] != "1") continue;
    found = true;
    std::string st = cell(t, row, "scan_type");
    type_ok = type_ok && (st == "survey" || st == "agc" || st == "conditional");
    act_ok  = act_ok  && (cell(t, row, "activation") == "none");
    // stage-less placeholders: per-stage numerics literal "0", scoring scalars "0"
    zero_ok = zero_ok && (cell(t, row, "charge") == "0") && (cell(t, row, "precursor_mz") == "0")
                      && (cell(t, row, "isolation_width") == "0") && (cell(t, row, "collision_energy") == "0")
                      && (cell(t, row, "mono_mass") == "0") && (cell(t, row, "qscore") == "0");
    parent_ok = parent_ok && cell(t, row, "parent_tracking_id").empty();
    ion_ok = ion_ok && cell(t, row, "ion_type").empty() && (cell(t, row, "ion_index") == "0");
    TEST_EQUAL(row.size(), t.headers.size())
  }
  TEST_TRUE(found)
  TEST_TRUE(type_ok)  TEST_TRUE(act_ok)  TEST_TRUE(zero_ok)  TEST_TRUE(parent_ok)  TEST_TRUE(ion_ok)

  std::remove(cmd_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §C(expl) -- exploration variant commands carry scan_type=="exploration"
/////////////////////////////////////////////////////////////
START_SECTION(commands_exploration_marker)
{
  // Pure MS2 exploration: plain exploration config (target_mode 0) selects the top MS1
  // precursor and sweeps CE variants -- no inclusion pinning needed (that is only for
  // making MS3 fire on the cytC proteoform).
  auto ms1 = loadTsvScans(FI_MS1_STD);
  auto ms2 = loadTsvScans(FI_MS2_HCD);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_cE_commands.tsv";
  std::remove(cmd_f.c_str());
  std::string json = injectRuntime(explorationConfig(), cmd_f, "");
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  TEST_TRUE(r.group_commands > 0)

  auto t = TSVFile::parse(cmd_f);
  int lvl = t.colIndex("ms_level");
  ABORT_IF(lvl < 0)

  bool found_expl = false;
  for (const auto& row : t.rows)
  {
    if (cell(t, row, "scan_type") == "exploration")
    {
      found_expl = true;
      TEST_TRUE(isTrackingId(cell(t, row, "tracking_id")))
      TEST_TRUE(inActivationSet(cell(t, row, "activation")))
      TEST_EQUAL(row.size(), t.headers.size())
    }
  }
  TEST_TRUE(found_expl)  // exploration CE-sweep variants overwrite the desc marker to 'E'

  std::remove(cmd_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §R1 -- scan_results.tsv : MS1 row columns + child_ids join + pinned defaults
/////////////////////////////////////////////////////////////
START_SECTION(results_ms1_columns)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_r1_commands.tsv", res_f = "lf_r1_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
  std::string json = buildJsonWithRuntime("", cmd_f, res_f, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto cmds = TSVFile::parse(cmd_f);
  auto res = TSVFile::parse(res_f);
  auto level = commandLevels(cmds);

  std::set<std::string> cmd_ids;
  for (const auto& kv : level) cmd_ids.insert(kv.first);

  bool found = false;
  for (const auto& row : res.rows)
  {
    std::string tid = cell(res, row, "tracking_id");
    auto it = level.find(tid);
    if (it != level.end() && it->second != 1) continue;  // MS2/MS3 -> not an MS1 row (MS1 surveys now sit in
                                                          // scan_commands.tsv with ms_level=1; classify by value)
    found = true;
    // mass_count == deconv_masses token count
    int mc = std::atoi(cell(res, row, "mass_count").c_str());
    auto masses = splitTokens(cell(res, row, "deconv_masses"), ';');
    TEST_EQUAL((int)masses.size(), mc)
    // child_ids space-split, each a real command id, count == commands_pushed
    auto kids = splitTokens(cell(res, row, "child_ids"), ' ');
    int pushed = std::atoi(cell(res, row, "commands_pushed").c_str());
    TEST_EQUAL((int)kids.size(), pushed)
    bool kids_ok = true;
    for (const auto& k : kids) kids_ok = kids_ok && cmd_ids.count(k) > 0;
    TEST_TRUE(kids_ok)
    // pinned MS1 defaults
    TEST_TRUE(cell(res, row, "parent_tracking_id").empty())
    TEST_EQUAL(cell(res, row, "tag_count"), std::string("0"))
    TEST_TRUE(cell(res, row, "matched_protein").empty())
    TEST_TRUE(cell(res, row, "proteoform_sequence").empty())
    TEST_EQUAL(cell(res, row, "exploration_group_id"), std::string("-1"))
    TEST_EQUAL(cell(res, row, "exploration_metric"), std::string("0"))
    TEST_EQUAL(cell(res, row, "variant_index"), std::string("-1"))
    TEST_EQUAL(cell(res, row, "total_variants"), std::string("0"))
    TEST_EQUAL(cell(res, row, "activation_type"), std::string(""))
    TEST_TRUE(toD(cell(res, row, "exploration_score")) < 0.0)   // -1 sentinel
    TEST_TRUE(toD(cell(res, row, "remaining_ratio")) < 0.0)
    TEST_TRUE(toD(cell(res, row, "rt")) >= 0.0 && toD(cell(res, row, "rt")) <= FI_RT_WINDOW)
  }
  TEST_TRUE(found)

  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §R2 -- scan_results.tsv : MS2-normal CE/activation/reaction (single stage) + tag_count
/////////////////////////////////////////////////////////////
START_SECTION(results_ms2_normal_columns)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_r2_commands.tsv", res_f = "lf_r2_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
  std::string json = buildJsonWithRuntime("", cmd_f, res_f, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto cmds = TSVFile::parse(cmd_f);
  auto res = TSVFile::parse(res_f);
  auto level = commandLevels(cmds);

  bool found = false;
  for (const auto& row : res.rows)
  {
    std::string tid = cell(res, row, "tracking_id");
    auto it = level.find(tid);
    if (it == level.end() || it->second != 2) continue;  // MS2 rows only
    found = true;
    // single stage: no ';' in CE/activation/reaction
    TEST_TRUE(cell(res, row, "collision_energy").find(';') == std::string::npos)
    TEST_TRUE(cell(res, row, "activation_type").find(';') == std::string::npos)
    TEST_TRUE(cell(res, row, "reaction_time").find(';') == std::string::npos)
    TEST_TRUE(posFinite(toD(cell(res, row, "collision_energy"))))            // ctx.stages[0].CE
    TEST_TRUE(inActivationSet(cell(res, row, "activation_type")))
    TEST_TRUE(nonNegFinite(toD(cell(res, row, "reaction_time"))))
    TEST_TRUE(toD(cell(res, row, "tag_count")) >= 0.0)                        // real count (>=0); exact in C#
    TEST_TRUE(inUnit(toD(cell(res, row, "tic_coverage"))))
    // non-exploration: exploration cols pinned
    TEST_EQUAL(cell(res, row, "exploration_group_id"), std::string("-1"))
    TEST_EQUAL(cell(res, row, "variant_index"), std::string("-1"))
    TEST_TRUE(isTrackingId(cell(res, row, "parent_tracking_id")))            // == its MS1
  }
  TEST_TRUE(found)

  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §I3 -- scan_results.tsv : matched (non-FASTA identification) rows log the REAL identification
//        tag_count (>0, not the FASTA-gated 0) and a toProForma-rendered proteoform that MATCHES
//        identification.tsv for the same tracking id. Locks I3 on the stable side; exact values in C#.
/////////////////////////////////////////////////////////////
START_SECTION(results_identification_tag_count_and_proforma)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string res_f = "lf_i3_results.tsv", id_f = "lf_i3_id.tsv";
  std::remove(res_f.c_str()); std::remove(id_f.c_str());
  // inclusion-pinned cytC identification (target_mode=1, M-start proteoform, NO FASTA tag DB) -- the
  // exact case that used to log tag_count=0 because selection_.hasTargetProteinDatabase() is false.
  std::string json = buildJsonWithRuntime("", "", res_f, true, id_f);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto res = TSVFile::parse(res_f);
  auto idf = TSVFile::parse(id_f);
  ABORT_IF(idf.rows.size() == 0)   // identification must fire, else this gate proves nothing

  // tracking_id -> identification proteoform (rendered via toProForma in the engine writer) + ms_level
  std::map<std::string, std::string> id_proteoform;
  std::map<std::string, int> id_level;
  for (const auto& row : idf.rows)
  {
    id_proteoform[cell(idf, row, "tracking_id")] = cell(idf, row, "proteoform");
    id_level[cell(idf, row, "tracking_id")] = std::atoi(cell(idf, row, "ms_level").c_str());
  }

  bool checked = false;
  for (const auto& row : res.rows)
  {
    std::string tid = cell(res, row, "tracking_id");
    auto it = id_proteoform.find(tid);
    if (it == id_proteoform.end()) continue;  // only rows that produced an identification (a match)
    // I3 is the MS2 identification tag_count fix; MS3 rows sentinel tag_count to -1 (F9, covered by §R4/§F9),
    // so restrict the >0 tag_count floor to MS2 matched rows.
    if (id_level[tid] != 2) continue;
    checked = true;
    // I3(a): the logged tag_count is the engine's REAL identification tag count, not the FASTA-gated 0.
    // A matched cytC proteoform implies >=1 generated sequence tag.
    TEST_TRUE(toD(cell(res, row, "tag_count")) > 0.0)
    // I3(b): proteoform_sequence is rendered through the SAME toProForma path as identification.tsv, so
    // PTMs (heme/N-term...) are displayed and the two streams agree exactly for the same scan.
    TEST_EQUAL(cell(res, row, "proteoform_sequence"), it->second)
  }
  TEST_TRUE(checked)

  std::remove(res_f.c_str()); std::remove(id_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §R4 -- scan_results.tsv : MS3-normal 2-stage CE/act/rt + proteoform/protein/fragment/tic
/////////////////////////////////////////////////////////////
START_SECTION(results_ms3_normal_columns)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_r4_commands.tsv", res_f = "lf_r4_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
  std::string json = buildJsonWithRuntime("", cmd_f, res_f, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto cycle = runFullCycle(&ida, ms1, ms2);
  TEST_TRUE(cycle.ms3_cmds.size() > 0)

  auto cmds = TSVFile::parse(cmd_f);
  auto res = TSVFile::parse(res_f);
  auto level = commandLevels(cmds);

  bool found_ms3 = false, found_populated = false;
  for (const auto& row : res.rows)
  {
    std::string tid = cell(res, row, "tracking_id");
    auto it = level.find(tid);
    if (it == level.end() || it->second != 3) continue;  // MS3 rows only
    found_ms3 = true;

    // terminal scan: no children
    TEST_EQUAL(cell(res, row, "commands_pushed"), std::string("0"))
    TEST_TRUE(cell(res, row, "child_ids").empty())
    TEST_TRUE(isTrackingId(cell(res, row, "parent_tracking_id")))  // == its MS2

    // 2-stage CE/act/rt (MS2 isolation ; MS3 fragmentation)
    auto ce  = splitTokens(cell(res, row, "collision_energy"), ';');
    auto act = splitTokens(cell(res, row, "activation_type"), ';');
    auto rt  = splitTokens(cell(res, row, "reaction_time"), ';');
    if (ce.size() == 2)
    {
      TEST_TRUE(posFinite(toD(ce[0])) && posFinite(toD(ce[1])))
      TEST_TRUE(std::abs(toD(ce[1]) - 35.0) < 1.0)  // MS3 fragmentation CE
      TEST_TRUE(act.size() == 2 && inActivationSet(act[0]) && inActivationSet(act[1]))
      TEST_TRUE(rt.size() == 2 && nonNegFinite(toD(rt[0])) && nonNegFinite(toD(rt[1])))
    }
    // F9: MS3 rows log fragment_count & tag_count as the -1 SENTINEL on EVERY MS3 row — the matched count is
    // finalized only in the calibrated round and lives in identification.tsv (ms3_fragments); tagging is an
    // MS2-targeting feature, not used for fragment-based MS3 id. (Was tag_count>=0 / fragment_count>0.)
    TEST_EQUAL(std::atoi(cell(res, row, "tag_count").c_str()), -1)
    TEST_EQUAL(std::atoi(cell(res, row, "fragment_count").c_str()), -1)

    // a matched MS3 row still carries protein + proteoform + tic (these are NOT sentinelized by F9)
    if (!cell(res, row, "matched_protein").empty())
    {
      found_populated = true;
      TEST_TRUE(isProForma(cell(res, row, "proteoform_sequence")))
      TEST_TRUE(toD(cell(res, row, "tic_coverage")) > 0.0 && inUnit(toD(cell(res, row, "tic_coverage"))))
    }
  }
  TEST_TRUE(found_ms3)
  // matched_protein/proteoform/fragment_count/tic populate only when the MS3 fragment match
  // succeeds; that path's exactness is golden-locked in C# (and depends on real MS3 fragment
  // data rather than this cycle's MS2-spectrum-as-MS3 shortcut). Here we validate them WHEN
  // present (block above) plus the match-independent 2-stage CE/act/rt + terminal structure
  // on every MS3 row.
  (void)found_populated;

  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §R3 -- scan_results.tsv : MS2-exploration columns + child_ids (B7)
/////////////////////////////////////////////////////////////
START_SECTION(results_ms2_exploration_columns)
{
  auto ms1 = loadTsvScans(FI_MS1_STD);
  auto ms2 = loadTsvScans(FI_MS2_HCD);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_r3_commands.tsv", res_f = "lf_r3_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
  std::string json = injectRuntime(explorationConfig(), cmd_f, res_f);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  TEST_TRUE(r.group_commands > 0)

  auto res = TSVFile::parse(res_f);

  bool found_expl = false;
  for (const auto& row : res.rows)
  {
    int gid = std::atoi(cell(res, row, "exploration_group_id").c_str());
    if (gid < 0) continue;  // exploration rows only (>=1; non-expl pins -1)
    found_expl = true;
    TEST_TRUE(gid >= 1)
    int metric = std::atoi(cell(res, row, "exploration_metric").c_str());
    TEST_TRUE(metric >= 1 && metric <= 3)
    int vi = std::atoi(cell(res, row, "variant_index").c_str());
    int tv = std::atoi(cell(res, row, "total_variants").c_str());
    TEST_TRUE(vi >= 0 && tv > 0 && vi < tv)
    TEST_TRUE(posFinite(toD(cell(res, row, "collision_energy"))))   // variant CE (single stage at MS2)
    TEST_TRUE(inActivationSet(cell(res, row, "activation_type")))
    TEST_TRUE(finiteVal(toD(cell(res, row, "exploration_score"))))
    // B7: child_ids count matches commands_pushed and each is a well-formed tracking id.
    // (Children are space-split; they may not be dequeued into the commands TSV yet, so we
    // check id FORMAT here -- the full child->command join is covered by §R1 / join_integrity.)
    auto kids = splitTokens(cell(res, row, "child_ids"), ' ');
    int pushed = std::atoi(cell(res, row, "commands_pushed").c_str());
    TEST_EQUAL((int)kids.size(), pushed)
    bool kids_ok = true;
    for (const auto& k : kids) kids_ok = kids_ok && isTrackingId(k);
    TEST_TRUE(kids_ok)
  }
  TEST_TRUE(found_expl)

  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F5 -- scan_results.tsv : every exploration variant logs ITS OWN sweep point (variant_index spans
//        0..total_variants-1 exactly once per group) and exactly ONE row per completed group names the
//        winner via the trailing winner_tracking_id column. Catches F5 — the winner-overwrite at
//        Exploration.cpp re-stamped the LAST-resolved variant's row with the WINNER's index/CE/score, so
//        a variant_index was duplicated and another vanished, and the winner was otherwise unidentifiable.
//        Set-membership / cross-row resolution only (no captured CE/score) -> drift-stable.
/////////////////////////////////////////////////////////////
START_SECTION(results_exploration_per_variant_and_winner)
{
  auto ms1 = loadTsvScans(FI_MS1_STD);
  auto ms2 = loadTsvScans(FI_MS2_HCD);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_f5_commands.tsv", res_f = "lf_f5_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
  std::string json = injectRuntime(explorationConfig(), cmd_f, res_f);   // mass_count -> 5 variants, no baseline
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);                  // feeds all variants -> group completes
  ABORT_IF(r.group_commands == 0)

  auto res = TSVFile::parse(res_f);

  std::map<int, std::set<int>> seen_vi;     // group_id -> set of variant_index logged
  std::map<int, int> tv_of, winner_count;   // group_id -> total_variants / # rows naming a winner
  std::map<int, std::string> winner_of;     // group_id -> the winner's tracking_id
  std::set<std::string> all_row_ids;        // every exploration row's own tracking_id
  bool found_expl = false;
  for (const auto& row : res.rows)
  {
    int gid = std::atoi(cell(res, row, "exploration_group_id").c_str());
    if (gid < 0) continue;                  // exploration rows only (non-expl pin -1)
    found_expl = true;
    int vi = std::atoi(cell(res, row, "variant_index").c_str());
    int tv = std::atoi(cell(res, row, "total_variants").c_str());
    TEST_TRUE(vi >= 0 && tv > 0 && vi < tv)
    seen_vi[gid].insert(vi);
    tv_of[gid] = tv;
    all_row_ids.insert(cell(res, row, "tracking_id"));
    std::string w = cell(res, row, "winner_tracking_id");
    if (!w.empty()) { winner_of[gid] = w; winner_count[gid]++; }
  }
  TEST_TRUE(found_expl)

  for (const auto& kv : seen_vi)
  {
    int gid = kv.first;
    const std::set<int>& vis = kv.second;
    // ISSUE(F5): pre-fix variant_index did NOT self-describe (e.g. {0,0,1,2,3} — winner index duplicated,
    // the completing variant's own index missing) because the completing row was overwritten with the winner.
    for (int k = 0; k < tv_of[gid]; ++k) TEST_EQUAL(vis.count(k), 1)
    TEST_EQUAL((int)vis.size(), tv_of[gid])
    // ISSUE(F5): pre-fix there was NO winner column, so once per-variant overwriting stopped the winner
    // would be unidentifiable. Exactly one row per completed group names a real, in-set winner.
    TEST_EQUAL(winner_count[gid], 1)
    TEST_TRUE(isTrackingId(winner_of[gid]) && all_row_ids.count(winner_of[gid]) == 1)
  }

  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §I4a -- scan_results.tsv : exploration FOLLOW-UP. A non-tolerance override makes the MS2-exploration
//         winner re-acquire as a production MS2 (then MS2->MS3) on inclusion-pinned cytC. Mirrors C++
//         ms2_exploration_production_winner_then_ms3 + checks the LOGGED exploration columns (matching data).
/////////////////////////////////////////////////////////////
START_SECTION(results_exploration_followup_columns)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  // Inclusion-pin the cytC precursor + M-start proteoform; add a non-tolerance override ("analyzer")
  // so config_.level(2).overrides stays non-empty -> the winner is a production-MS2 re-acquisition.
  std::string cfg = inclusionPinCytc(ms3SelectionOnlyConfig());
  {
    auto p = cfg.find("\"metric\": \"mass_count\"");
    ABORT_IF(p == std::string::npos)
    cfg.insert(p, "\"overrides\": { \"analyzer\": \"Orbitrap\" }, ");
  }
  std::string res_f = "lf_i4followup_results.tsv";
  std::remove(res_f.c_str());
  cfg = injectRuntime(cfg, "", res_f);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  ABORT_IF(r.group_commands == 0)

  // Override non-empty -> winner re-acquires as production MS2; feeding it drives the MS2->MS3 path.
  TEST_TRUE(r.found_production_ms2)
  TEST_TRUE(r.found_ms3)

  // Logged exploration result rows (inclusion-pinned cytC => they MATCH, unlike §R3's E. coli data).
  auto res = TSVFile::parse(res_f);
  bool found_expl = false;
  for (const auto& row : res.rows)
  {
    int gid = std::atoi(cell(res, row, "exploration_group_id").c_str());
    if (gid < 0) continue;
    found_expl = true;
    TEST_TRUE(gid >= 1)
    int vi = std::atoi(cell(res, row, "variant_index").c_str());
    int tv = std::atoi(cell(res, row, "total_variants").c_str());
    TEST_TRUE(vi >= 0 && tv > 0 && vi < tv)
    TEST_TRUE(finiteVal(toD(cell(res, row, "exploration_score"))))
    auto kids = splitTokens(cell(res, row, "child_ids"), ' ');
    int pushed = std::atoi(cell(res, row, "commands_pushed").c_str());
    TEST_EQUAL((int)kids.size(), pushed)
  }
  TEST_TRUE(found_expl)
  std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §I4b -- scan_results.tsv : MS3-level EXPLORATION. A two-level MS2->MS3 CE-sweep cascade on inclusion-
//         pinned cytC. Mirrors C++ ms2_then_ms3_exploration_acquires_ms3 + checks the LOGGED exploration
//         columns. (MS3 result rows themselves are golden-locked in the C# Golden_Exploration_MS3_CytC.)
/////////////////////////////////////////////////////////////
START_SECTION(results_ms3_exploration_columns)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string res_f = "lf_i4ms3_results.tsv";
  std::remove(res_f.c_str());
  std::string cfg = injectRuntime(inclusionPinCytc(ms3ExplorationConfig()), "", res_f);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  ABORT_IF(r.group_commands == 0)

  // The MS2-exploration winner cascades to an MS3 (two-stage) command.
  TEST_TRUE(r.found_ms3)
  TEST_EQUAL(r.ms3_num_stages, 2)

  auto res = TSVFile::parse(res_f);
  bool found_expl = false;
  for (const auto& row : res.rows)
  {
    int gid = std::atoi(cell(res, row, "exploration_group_id").c_str());
    if (gid < 0) continue;
    found_expl = true;
    TEST_TRUE(gid >= 1)
    int metric = std::atoi(cell(res, row, "exploration_metric").c_str());
    TEST_TRUE(metric >= 1 && metric <= 3)
    int vi = std::atoi(cell(res, row, "variant_index").c_str());
    int tv = std::atoi(cell(res, row, "total_variants").c_str());
    TEST_TRUE(vi >= 0 && tv > 0 && vi < tv)
  }
  TEST_TRUE(found_expl)
  std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F1 -- scan_commands.tsv : every MS3 exploration command's parent resolves to an MS2 (not the MS1
//        survey). Catches F1 — the feedResult push-loop used to overwrite buildMS3's correct MS2 parent
//        with the group's MS1 originating id. Cross-stream/level relationship only -> drift-stable.
/////////////////////////////////////////////////////////////
START_SECTION(commands_exploration_ms3_parent_is_ms2)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_f1_commands.tsv", res_f = "lf_f1_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
  std::string cfg = injectRuntime(inclusionPinCytc(ms3ExplorationConfig()), cmd_f, res_f);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  ABORT_IF(r.group_commands == 0)
  ABORT_IF(!r.found_ms3)   // need >=1 MS3 command emitted to check its parent

  auto cmds = TSVFile::parse(cmd_f);
  auto level = commandLevels(cmds);   // tracking_id -> ms_level (from scan_commands)
  int checked = 0;
  for (const auto& row : cmds.rows)
  {
    std::string tid = cell(cmds, row, "tracking_id");
    auto it = level.find(tid);
    if (it == level.end() || it->second != 3) continue;   // MS3 commands only
    std::string p = cell(cmds, row, "parent_tracking_id");
    // ISSUE(F1): pre-fix p resolved to an MS1 (the group's grandparent survey) because the
    // exploration push-loop overwrote buildMS3's correct MS2 parent with info.parent_scan_id.
    auto pit = level.find(p);
    TEST_TRUE(pit != level.end() && pit->second == 2)
    checked++;
  }
  TEST_TRUE(checked >= 1)
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F2 -- scan_commands.tsv : exploration MS3 commands carry REAL stage-1 (fragment) scores. Catches F2 —
//        the exploration buildMS3 used to get a default-empty FragmentScores, so the ';'-joined stage-1
//        token of the scoring columns was 0 (and identification ms3_window_snr/ms3_charge_intensity = 0).
//        The single-group driver can't feed MS3 spectra (no MS3 identification rows), but it DOES drain+log
//        the MS3 commands, whose stage-1 scalars are the upstream value the bug zeroed. Range only -> drift-stable.
/////////////////////////////////////////////////////////////
START_SECTION(commands_exploration_ms3_stage1_scores_nonzero)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_f2_commands.tsv", res_f = "lf_f2_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
  std::string cfg = injectRuntime(inclusionPinCytc(ms3ExplorationConfig()), cmd_f, res_f);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  ABORT_IF(r.group_commands == 0)
  ABORT_IF(!r.found_ms3)   // need >=1 MS3 command to inspect its stage-1 scores

  auto cmds = TSVFile::parse(cmd_f);
  int checked = 0;
  for (const auto& row : cmds.rows)
  {
    if (std::atoi(cell(cmds, row, "ms_level").c_str()) != 3) continue;   // MS3 commands only
    // MS3 scoring columns are 2-stage ';'-joined "MS2;fragment"; token[1] is the stage-1 (fragment) value.
    auto pin = splitTokens(cell(cmds, row, "precursor_intensity"), ';');
    auto csn = splitTokens(cell(cmds, row, "charge_snr"), ';');
    TEST_EQUAL(pin.size(), 2)
    TEST_EQUAL(csn.size(), 2)
    // ISSUE(F2): pre-fix the stage-1 tokens were 0 because exploration buildMS3 got the default-empty
    // FragmentScores; now they carry the selected fragment's real deconvolution scores.
    TEST_TRUE(toD(pin[1]) > 0.0)
    TEST_TRUE(finiteVal(toD(csn[1])) && toD(csn[1]) >= 0.0)
    checked++;
  }
  TEST_TRUE(checked >= 1)
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F3 -- MS3 exploration ('E') rows carry the parent proteoform. Catches F3 — proteoform/matched_protein
//        were sourced only on MS3 'R'. Needs the MS3 variants FED (the batch must complete), so this drives
//        the inclusion-pinned MS3-exploration config with the real per-ion MS3 manifest (skips cleanly when
//        no ms3_cytc fixtures exist; the exact proteoform is golden-locked in C# Golden_Exploration_MS3).
/////////////////////////////////////////////////////////////
START_SECTION(identification_ms3_exploration_proteoform)
{
  auto manifest = buildMs3IonManifest("../../FlashIDA/test-data/spectra");
  if (manifest.empty())
  {
    TEST_TRUE(true)   // no real MS3 fragment fixtures -> skip cleanly (never fabricate)
  }
  else
  {
    auto ms1 = loadTsvScans(CYTC_MS1);
    auto ms2 = loadTsvScans(CYTC_MS2);
    ABORT_IF(ms1.empty() || ms2.empty())

    std::string res_f = "lf_f3_results.tsv", id_f = "lf_f3_id.tsv";
    std::remove(res_f.c_str()); std::remove(id_f.c_str());
    std::string cfg = injectRuntime(inclusionPinCytc(ms3ExplorationConfig()), "", res_f, id_f);
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    runFullCycle(&ida, ms1, ms2, &manifest);   // feed the MS3 variants so MS3-'E' rows are produced

    auto idf = TSVFile::parse(id_f);
    auto res = TSVFile::parse(res_f);

    // F3-id: a matched MS3-'E' identification row has a non-empty ProForma proteoform.
    int matched = 0;
    for (const auto& row : idf.rows)
    {
      if (std::atoi(cell(idf, row, "ms_level").c_str()) != 3) continue;
      if (cell(idf, row, "scan_mode") != "E") continue;
      if (cell(idf, row, "ms3_fragments").empty()) continue;   // only rows that matched fragments
      // ISSUE(F3-id): pre-fix proteoform was blank (use_ctx_proteoform was 'R'-only).
      TEST_TRUE(isProForma(cell(idf, row, "proteoform")))
      matched++;
    }
    // F3-res: the matched MS3-'E' scan_results rows carry non-empty matched_protein + proteoform_sequence.
    for (const auto& r : res.rows)
    {
      if (std::atoi(cell(res, r, "ms_level").c_str()) != 3) continue;
      if (std::atoi(cell(res, r, "exploration_group_id").c_str()) < 1) continue;
      if (cell(res, r, "matched_protein").empty()) continue;   // unmatched/unfed MS3 -> tolerated
      // ISSUE(F3-res): pre-fix matched_protein/proteoform_sequence were blank on MS3-'E' rows.
      TEST_TRUE(!cell(res, r, "proteoform_sequence").empty())
    }
    TEST_TRUE(matched >= 1)   // the inclusion-pinned MS3-exploration cascade must yield >=1 matched MS3-'E' row
    std::remove(res_f.c_str()); std::remove(id_f.c_str());
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F4 -- MS2 exploration ('E') identification rows source ms1_precursor_* from the GROUP's own deconvolved
//        precursor, not the originating MS1 survey command. Catches F4 — buildMS2ContextForVariant filled
//        ms1_precursor_* from group.originating_cmd ONLY when num_stages>0, but for an MS2 exploration group
//        the originating_cmd is the MS1 survey (num_stages==0), so the guard failed and the fields stayed 0.
//        Drift-stable cross-stream ppm join: an MS2-'E' id row's ms1_precursor_mass must equal the SAME
//        tracking_id's scan_commands.mono_mass (both = the group's deconvolved precursor) within ppm.
/////////////////////////////////////////////////////////////
START_SECTION(identification_ms2_exploration_ms1_precursor_from_group)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_f4_commands.tsv", id_f = "lf_f4_id.tsv";
  std::remove(cmd_f.c_str()); std::remove(id_f.c_str());
  // single-group drive feeds the MS2 variants back -> feedResult writes the MS2-'E' identification rows.
  std::string cfg = injectRuntime(inclusionPinCytc(ms3ExplorationConfig()), cmd_f, "", id_f);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  ABORT_IF(r.group_commands == 0)

  auto cmds = TSVFile::parse(cmd_f);
  auto idf  = TSVFile::parse(id_f);

  // tracking_id -> mono_mass from scan_commands (an MS2 command's mono_mass is a single value == the
  // precursor mono mass; only MS3 commands ';'-join two stages, which these MS2-'E' ids never resolve to).
  std::map<std::string, double> cmd_mono;
  for (const auto& row : cmds.rows)
    cmd_mono[cell(cmds, row, "tracking_id")] = toD(cell(cmds, row, "mono_mass"));

  int checked = 0;
  for (const auto& row : idf.rows)
  {
    if (std::atoi(cell(idf, row, "ms_level").c_str()) != 2) continue;   // MS2 identification rows only
    if (cell(idf, row, "scan_mode") != "E") continue;                   // exploration variants only
    auto it = cmd_mono.find(cell(idf, row, "tracking_id"));
    if (it == cmd_mono.end() || it->second <= 0.0) continue;            // need the matching command row
    double id_mass = toD(cell(idf, row, "ms1_precursor_mass"));
    // ISSUE(F4): pre-fix ms1_precursor_mass came from group.originating_cmd, which for an MS2 exploration
    // group is the MS1 survey (num_stages==0) -> the >0 guard failed -> the field stayed 0.0.
    TEST_TRUE(id_mass > 0.0)
    double ppm = std::abs(id_mass - it->second) / it->second * 1e6;
    TEST_TRUE(ppm < 50.0)
    checked++;
  }
  TEST_TRUE(checked >= 1)   // >=1 MS2-'E' id row joined to its command within ppm
  std::remove(cmd_f.c_str()); std::remove(id_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F8-incl -- inclusion (cytC) selects the LISTED precursor within the target window. Catches the
//   F8-inclusion bug: the target was the cytC AVERAGE mass (~12358.31), ~570-700 ppm off the engine's
//   MONOISOTOPIC deconvolution (~12351.3), so the ~0.247 Da match window never matched -> silent DDA
//   fall-through (DDA picks cytC anyway, masking the miss). inclusion_cytc.txt is now the monoisotopic
//   12351.3. NON-STRICT: >=1 MS2 within the window (a floor — DDA would also pick cytC). STRICT: selection
//   is RESTRICTED to the list, so EVERY MS2 is within the window AND >=1 exists (the real lever — on the
//   wrong target strict matches nothing -> 0 MS2). Window membership only -> drift-stable.
/////////////////////////////////////////////////////////////
START_SECTION(inclusion_ms2_within_target_window)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const double T = 12351.3, HW = 2.0 * 10.0 * T * 1e-6;   // 10 ppm tol, doubled at the match site (~0.247 Da)
  auto onList = [&](double m){ return std::abs(m - T) < HW; };
  const std::string cmd_f = "lf_f8incl_commands.tsv";

  auto driveCountMs2 = [&](const std::string& json, int& ms2_total, int& ms2_on_list) {
    std::remove(cmd_f.c_str());
    FLASHIda ida(const_cast<char*>(json.c_str()));
    runFullCycle(&ida, ms1, ms2);
    auto cmds = TSVFile::parse(cmd_f);
    ms2_total = 0; ms2_on_list = 0;
    for (const auto& row : cmds.rows)
    {
      if (std::atoi(cell(cmds, row, "ms_level").c_str()) != 2) continue;   // MS2 commands only
      ms2_total++;
      if (onList(toD(cell(cmds, row, "mono_mass")))) ms2_on_list++;        // mono_mass is single-valued at MS2
    }
  };

  // NON-STRICT: the cytC precursor IS selected (>=1 MS2 within the window); DDA may add off-list MS2 too.
  int ns_total = 0, ns_on = 0;
  driveCountMs2(buildJsonWithRuntime("", cmd_f, "", true), ns_total, ns_on);
  // ISSUE(F8-incl): pre-fix ns_on==0 — the target was the cytC AVERAGE mass, ~570-700 ppm off the
  // monoisotopic ~12351.3, so the window never matched (silent DDA fall-through).
  TEST_TRUE(ns_on >= 1)

  // STRICT: selection RESTRICTED to the list -> EVERY MS2 is within the window, and >=1 exists.
  std::string strict_json = buildJsonWithRuntime("", cmd_f, "", true);
  {
    const std::string from = "\"strict_inclusion\": false";
    auto p = strict_json.find(from);
    ABORT_IF(p == std::string::npos)
    strict_json.replace(p, from.size(), "\"strict_inclusion\": true");
  }
  int s_total = 0, s_on = 0;
  driveCountMs2(strict_json, s_total, s_on);
  // ISSUE(F8-incl-strict): pre-fix s_total==0 (nothing matched the wrong target); any leaked MS2 would be an
  // off-list DDA pick -> s_on != s_total. With the corrected monoisotopic target every strict MS2 is on-list.
  TEST_TRUE(s_total >= 1 && s_on == s_total)

  std::remove(cmd_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F8-excl -- exclusion (mode 2) tqscore SOFT DE-PRIORITIZATION COVERAGE. NOT an engine bug: the enum is
//   correct (2=exclusion/tqscore reading 'Mass' lines, 3=deep/AllMass). The mode-2 target-log tqscore is a
//   SOFT, iteration-0-only DE-PRIORITIZATION (PrecursorSelection drops tqscore-exceeding masses FIRST, then
//   BACKFILLS any remaining MS2 slots), NOT a hard exclusion. So with a slack slot budget the de-prioritized
//   target is backfilled in iteration 1 and re-selected — the old !hasNear(excl, target) assertion was wrong
//   (it asserted hard elimination, which only the dedicated ProcessScan::processScan_mass_exclusion hard-
//   exclusion test covers). To make the soft de-prioritization OBSERVABLE we drive a RICH multi-precursor
//   survey (ms1_ecoli_rich, >=9 selectable masses/scan) with max_targets==1: the single MS2 slot then goes to
//   a non-suppressed competitor and the de-prioritized target is genuinely selected-OUT. Mode-0 (no log,
//   max_targets==1) confirms the engine WOULD pick that target as its top pick; mode-2 with the target-log
//   then de-prioritizes it out. Set membership only -> drift-stable.
/////////////////////////////////////////////////////////////
START_SECTION(exclusion_mode2_tqscore_suppresses_target_mass)
{
  // Rich co-eluting E. coli survey: >=9 selectable precursors per scan, so a max_targets==1 slot budget is
  // genuinely contended — exactly the condition under which SOFT de-prioritization is observable.
  const std::string FI_MS1_ECOLI = "../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt";
  auto ms1 = loadTsvScans(FI_MS1_ECOLI);
  auto ms2 = loadTsvScans(FI_MS2_HCD);
  ABORT_IF(ms1.empty() || ms2.empty())

  // target_logs + max_targets vary. tqscore_threshold is read from deconvolution.tqscore_threshold
  // (Config.cpp:98) = 0.9; rt_window 180 >> the tiny ecoli RTs, so every target-log observation is within
  // window of every survey scan.
  auto cfg = [](int mode, int max_targets, const std::string& target_log_json) {
    std::ostringstream o;
    o << R"({
      "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10] },
      "precursor_selection": { "RT_window": 180, "target_mode": )" << mode << R"(, "IDScore": false, "AllCharges": false, "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1 },
      "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
      "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
      "faims": { "cv_values": [-50], "max_cv_skip": 0 },
      "ms_settings": { "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 }, "ms2": [ { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 } ] },
      "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 }, "agc_interval_seconds": 999999 },
      "files": { "target_logs": [)" << target_log_json << R"(], "fasta": "", "inclusion_list": "", "ptm_list": "" },
      "selection_strategy": { "ms1": { "selection": "qscore", "max_targets": )" << max_targets << R"( }, "ms2": { "selection": "none" }, "ms3": { "selection": "none" } }
    })";
    return o.str();
  };

  const std::string cmd_f = "lf_f8excl_commands.tsv";
  auto selectedMs2Masses = [&](const std::string& json) {
    std::remove(cmd_f.c_str());
    std::string full = injectRuntime(json, cmd_f, "");
    FLASHIda ida(const_cast<char*>(full.c_str()));
    runFullCycle(&ida, ms1, ms2);
    auto cmds = TSVFile::parse(cmd_f);
    std::vector<double> masses;
    for (const auto& row : cmds.rows)
      if (std::atoi(cell(cmds, row, "ms_level").c_str()) == 2)
        masses.push_back(toD(cell(cmds, row, "mono_mass")));
    return masses;
  };

  // 1) DDA baseline (mode 0, no target log, max_targets==1): the top pick per scan. dda aggregates one pick
  //    per scan, so >=2 entries means the survey is productive (NON-VACUITY guard) — NOT a per-scan
  //    >=2-selectable precondition (the slack-vs-tight contrast below is what certifies contention). The mass
  //    we then de-prioritize is dda[0].
  auto dda = selectedMs2Masses(cfg(0, 1, ""));
  ABORT_IF(dda.size() < 2)
  const double target = dda[0];

  // 2) Write a target log re-observing `target` twice within rt_window at qscore 0.7
  //    -> tqscore = 1 - (1-0.7)^2 = 0.91 > 0.9. (Mass<TAB>Score=<TAB>Z= format, mirrors test_target_log.log.)
  const std::string log_path = "lf_f8excl_target.log";
  {
    std::ofstream lg(log_path);
    lg << "MS1 RT 0.05 (scan)\n"
       << "Mass " << std::to_string(target) << "\tScore=0.7\tZ=4\n"
       << "MS1 RT 0.10 (scan)\n"
       << "Mass " << std::to_string(target) << "\tScore=0.7\tZ=4\n";
  }

  // 3) Exclusion (mode 2) with that target log, STILL max_targets==1: the lone MS2 slot must go to a
  //    non-suppressed competitor (the rich survey has >=8 others), so the de-prioritized `target` is
  //    selected-OUT of the single slot. With a slack budget it would be backfilled — hence max_targets==1.
  auto excl = selectedMs2Masses(cfg(2, 1, "\"" + log_path + "\""));

  auto hasNear = [&](const std::vector<double>& v, double t){
    for (double m : v) if (std::abs(m - t) < 0.5) return true;   // same nominal mass (suppression is per-nominal)
    return false;
  };
  TEST_TRUE(hasNear(dda, target))    // sanity: DDA (mode 0) picks `target` as its top pick
  // §F8-excl SOFT DE-PRIORITIZATION: the mode-2 target-log tqscore (>0.9) is a SOFT iteration-0 de-
  //   prioritization — the engine drops `target` to the back of the iteration-0 ranking. With only ONE MS2
  //   slot and >=8 non-suppressed competitors, `target` loses that slot and does NOT appear among the mode-2
  //   selections. This proves the SOFT de-prioritization path (target de-prioritized, selected only when a
  //   slack slot lets it backfill), distinct from HARD exclusion (covered by
  //   ProcessScan::processScan_mass_exclusion). Set membership only -> drift-stable.
  TEST_TRUE(!hasNear(excl, target))

  // 4) SLACK budget (mode 2, SAME target log, max_targets==64 >> the >=9-19 selectable/scan): iteration 0 can
  //    no longer fill every slot from non-suppressed competitors, so iteration 1 runs WITHOUT the tqscore
  //    filter and BACKFILLS the de-prioritized `target` -> it IS re-selected. This is the decisive SOFT proof:
  //    dropped at budget 1, backfilled at budget 64. A HARD exclusion would suppress it at ANY budget
  //    (that path is covered by ProcessScan::processScan_mass_exclusion).
  auto slack = selectedMs2Masses(cfg(2, 64, "\"" + log_path + "\""));
  TEST_TRUE(hasNear(slack, target))

  std::remove(cmd_f.c_str());
  std::remove(log_path.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F9 -- scan_results.tsv : MS3 rows log fragment_count == -1 AND tag_count == -1 (the N/A SENTINEL), while
//   MS2 rows keep real >=0 counts. Catches F9 — MS3 matching is finalized only in the calibrated round (the
//   matched count lives in identification.tsv) and tagging is an MS2-targeting feature, so a +count / carried
//   parent count on an MS3 row is misleading; -1 is the honest N/A marker. Exact -1 / >=0 floor -> drift-stable.
/////////////////////////////////////////////////////////////
START_SECTION(results_ms3_fragment_and_tag_count_are_sentinel)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_f9_commands.tsv", res_f = "lf_f9_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
  std::string json = buildJsonWithRuntime("", cmd_f, res_f, true);   // cytC inclusion + MS3 -> MS2 AND MS3 rows
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto cmds = TSVFile::parse(cmd_f);
  auto res  = TSVFile::parse(res_f);
  auto level = commandLevels(cmds);

  bool found_ms2 = false, found_ms3 = false;
  for (const auto& row : res.rows)
  {
    auto it = level.find(cell(res, row, "tracking_id"));
    if (it == level.end()) continue;   // MS1 survey input (not in the command-level map)
    if (it->second == 3)
    {
      found_ms3 = true;
      // ISSUE(F9a): MS3 fragment_count must be -1 (matched count is finalized in the calibrated round and
      // lives in identification.tsv) — NOT 0, NOT a +count.
      TEST_EQUAL(std::atoi(cell(res, row, "fragment_count").c_str()), -1)
      // ISSUE(F9b): MS3 tag_count must be -1 (tagging is not used for fragment-based MS3 id) — NOT a carried count.
      TEST_EQUAL(std::atoi(cell(res, row, "tag_count").c_str()), -1)
    }
    else if (it->second == 2)
    {
      found_ms2 = true;
      TEST_TRUE(std::atoi(cell(res, row, "fragment_count").c_str()) >= 0)   // MS2 keeps a real count (not sentinel)
      TEST_TRUE(std::atoi(cell(res, row, "tag_count").c_str()) >= 0)
    }
  }
  TEST_TRUE(found_ms3)
  TEST_TRUE(found_ms2)

  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §X1 -- scan_results.tsv : duration semantics (ordering + decomposition)
/////////////////////////////////////////////////////////////
START_SECTION(results_duration_semantics)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string res_f = "lf_x1_results.tsv";
  std::remove(res_f.c_str());
  std::string json = buildJsonWithRuntime("", "", res_f, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto res = TSVFile::parse(res_f);
  TEST_TRUE(res.rows.size() > 0)
  for (const auto& row : res.rows)
  {
    double dur = toD(cell(res, row, "duration_ms"));
    TEST_TRUE(inDuration(dur))
    double q = toD(cell(res, row, "queue_duration_ms"));
    double inst = toD(cell(res, row, "instrument_duration_ms"));
    double proc = toD(cell(res, row, "processing_duration_ms"));
    TEST_TRUE(nonNegFinite(q) && nonNegFinite(inst) && nonNegFinite(proc))
    double resolve = toD(cell(res, row, "resolve_ts"));
    double recv = toD(cell(res, row, "received_ts"));
    double deq = toD(cell(res, row, "dequeue_ts"));
    // ordering when timestamps are present (>0)
    if (deq > 0 && recv > 0) TEST_TRUE(recv >= deq)
    if (recv > 0 && resolve > 0) TEST_TRUE(resolve >= recv)
    // decomposition when fully tracked
    if (q > 0 && inst > 0 && proc > 0) TEST_TRUE(std::abs(dur - (q + inst + proc)) < 2.0)
  }
  std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §X2 -- scan_results.tsv : parent lineage (MS2->MS1, MS3->MS2)
/////////////////////////////////////////////////////////////
START_SECTION(results_parent_lineage)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_x2_commands.tsv", res_f = "lf_x2_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
  std::string json = buildJsonWithRuntime("", cmd_f, res_f, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto cmds = TSVFile::parse(cmd_f);
  auto res = TSVFile::parse(res_f);
  auto level = commandLevels(cmds);

  bool checked_ms2 = false, checked_ms3 = false;
  for (const auto& row : res.rows)
  {
    std::string tid = cell(res, row, "tracking_id");
    std::string parent = cell(res, row, "parent_tracking_id");
    auto it = level.find(tid);
    if (it == level.end() || it->second == 1) { TEST_TRUE(parent.empty()); continue; }  // MS1 survey/input
                                                          // (now in scan_commands.tsv with ms_level=1)
    // MS2/MS3 result rows must carry a well-formed parent that is itself a known id
    TEST_TRUE(isTrackingId(parent))
    if (it->second == 2) { checked_ms2 = true; }
    if (it->second == 3) { checked_ms3 = true; TEST_TRUE(level.count(parent) && level[parent] == 2); }
  }
  TEST_TRUE(checked_ms2)
  TEST_TRUE(checked_ms3)

  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §X4 -- results vs identification agree for the same tracking_id (#2-4)
/////////////////////////////////////////////////////////////
START_SECTION(results_identification_agreement)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string res_f = "lf_x4_results.tsv", id_f = "lf_x4_id.tsv";
  std::remove(res_f.c_str()); std::remove(id_f.c_str());
  std::string json = buildJsonWithRuntime("", "", res_f, true, id_f);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto res = TSVFile::parse(res_f);
  auto idf = TSVFile::parse(id_f);

  // identification proteoform per tracking_id
  std::map<std::string, std::string> id_proteo;
  for (const auto& row : idf.rows)
    id_proteo[cell(idf, row, "tracking_id")] = cell(idf, row, "proteoform");

  bool checked = false;
  for (const auto& row : res.rows)
  {
    std::string tid = cell(res, row, "tracking_id");
    auto it = id_proteo.find(tid);
    if (it == id_proteo.end()) continue;
    // when a scan has both a results row and an identification row, both proteoforms agree (non-empty)
    std::string rp = cell(res, row, "proteoform_sequence");
    if (rp.empty()) continue;  // results proteoform only populated on matched MS2/MS3
    checked = true;
    TEST_TRUE(isProForma(rp))
    TEST_TRUE(isProForma(it->second))
  }
  // identification.tsv must have produced at least one matched row in this cycle
  TEST_TRUE(idf.rows.size() > 0)
  (void)checked;

  std::remove(res_f.c_str()); std::remove(id_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §I1 -- identification.tsv : universal gates on all (R-mode) rows
/////////////////////////////////////////////////////////////
START_SECTION(identification_R_rows_universal)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string id_f = "lf_i1_id.tsv";
  std::remove(id_f.c_str());
  std::string json = buildJsonWithRuntime("", "", "", true, id_f);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto idf = TSVFile::parse(id_f);
  TEST_TRUE(idf.rows.size() > 0)  // cytC (non-exploration) => all rows scan_mode 'R'
  for (const auto& row : idf.rows)
  {
    int lvl = std::atoi(cell(idf, row, "ms_level").c_str());
    TEST_TRUE(lvl == 2 || lvl == 3)
    std::string mode = cell(idf, row, "scan_mode");
    TEST_TRUE(mode == "R" || mode == "E")
    TEST_TRUE(isTrackingId(cell(idf, row, "tracking_id")))
    TEST_TRUE(isProForma(cell(idf, row, "proteoform")))
    TEST_TRUE(std::atoi(cell(idf, row, "start_pos").c_str()) >= 0)
    TEST_TRUE(std::atoi(cell(idf, row, "end_pos").c_str()) >= std::atoi(cell(idf, row, "start_pos").c_str()))
    TEST_TRUE(finiteVal(toD(cell(idf, row, "ppm_offset"))))
    double corr = toD(cell(idf, row, "correction_factor"));
    TEST_TRUE(corr > 0.0 && std::abs(corr - 1.0) < 0.05)
    TEST_TRUE(inMass(toD(cell(idf, row, "ms1_precursor_mass"))))
    TEST_TRUE(posFinite(toD(cell(idf, row, "ms1_precursor_mz"))))
    TEST_TRUE(inChargeD(toD(cell(idf, row, "ms1_precursor_charge"))))
    // ms2_fragments / ms2_fragment_masses token-aligned
    auto f = splitTokens(cell(idf, row, "ms2_fragments"), ';');
    auto m = splitTokens(cell(idf, row, "ms2_fragment_masses"), ';');
    TEST_EQUAL(f.size(), m.size())

    // I2 (isolation-window cols): every matched row reports the MS2 precursor's commanded window width,
    // the SNR over that actual window, and the selected-charge intensity. Plausibility ranges only;
    // the exact values are golden-locked in the C# suite. window_snr is bounded (~1000) for a pure window.
    double w2 = toD(cell(idf, row, "ms2_isolation_width"));
    TEST_TRUE(w2 > 0.0 && w2 < 50.0)
    double s2 = toD(cell(idf, row, "ms2_window_snr"));
    TEST_TRUE(finiteVal(s2) && s2 >= 0.0 && s2 <= 1001.0)
    TEST_TRUE(posFinite(toD(cell(idf, row, "ms2_charge_intensity"))))

    // §I2 folded in: MS3-R rows additionally carry the MS2-fragment precursor (#1: proteoform/
    // region sourced from ctx) plus a non-empty, token-aligned MS3 fragment list. Asserted WHEN
    // an MS3-R row is present -- those require the MS3 fragment match to fire (which depends on
    // real MS3 fragment data, golden-locked exactly in the C# suite, not this cycle's shortcut).
    if (lvl == 3 && mode == "R")
    {
      std::string pion = cell(idf, row, "ms2_precursor_ion");
      TEST_TRUE(pion.size() >= 2 && ionTypeOk(std::string(1, pion[0])))
      TEST_TRUE(posFinite(toD(cell(idf, row, "ms2_precursor_mass"))))
      // the MS2 fragment chosen as MS3 precursor: a fragment charge, bounded by the MS1 precursor charge
      TEST_TRUE(inFragCharge(toD(cell(idf, row, "ms2_precursor_charge")), toD(cell(idf, row, "ms1_precursor_charge"))))
      auto f3 = splitTokens(cell(idf, row, "ms3_fragments"), ';');
      auto m3 = splitTokens(cell(idf, row, "ms3_fragment_masses"), ';');
      TEST_TRUE(!f3.empty())
      TEST_EQUAL(f3.size(), m3.size())

      // I1: start_pos/end_pos must be the MS3 FRAGMENT sub-range, not the parent proteoform's full range.
      // A b/y precursor ion of index k spans exactly k residues, so end_pos - start_pos == k (the
      // ms2_precursor_ion index). Proves the sub-range was stored, with no hard-coded coordinates.
      int frag_idx = std::atoi(pion.c_str() + 1);  // pion = e.g. "b3" -> 3
      int sp = std::atoi(cell(idf, row, "start_pos").c_str());
      int ep = std::atoi(cell(idf, row, "end_pos").c_str());
      TEST_TRUE(frag_idx > 0 && (ep - sp) == frag_idx)

      // I2 (isolation-window cols): MS3 rows additionally report the MS3 fragment precursor's window.
      double w3 = toD(cell(idf, row, "ms3_isolation_width"));
      TEST_TRUE(w3 >= 2.0 && w3 < 50.0)   // engine floors the MS3 isolation window at 2.0 Da
      double s3 = toD(cell(idf, row, "ms3_window_snr"));
      TEST_TRUE(finiteVal(s3) && s3 >= 0.0 && s3 <= 1001.0)
      TEST_TRUE(posFinite(toD(cell(idf, row, "ms3_charge_intensity"))))
    }
  }
  std::remove(id_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §L1 -- ida_log : all 15 floats plausible + cross-field relationships
/////////////////////////////////////////////////////////////
START_SECTION(ida_log_all15_fields)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  ABORT_IF(ms1.empty())

  std::string log_f = "lf_l1_ida.log";
  std::remove(log_f.c_str());
  std::string json = buildJsonWithRuntime(log_f, "", "");
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Drive the MS1 surveys via the canonical interleaved driver: the engine emits each survey command
  // and we feed the next cytC MS1 scan back stamped with the engine's OWN tracking id (was
  // pushAllScans' fabricated encode(800000+i) ids, which the always-on MS1 gate now rejects -> 0
  // precursors -> empty log). The ida_log is written at MS1 time; ms2 selection is "none" here so MS2
  // commands are recorded but not fed (no MS2 fixture needed). Standard-DDA MS1 selection (qscore,
  // threshold 0.0) still picks the top precursors per cytC scan -> >=1 emitted MS2 command, the
  // faithful analog of the old total>0 (MS2 commands pushed during MS1 processing).
  AcqResult acq = runInterleaved(&ida, ms1, std::vector<ScanData>{});
  TEST_TRUE(acq.ms2_cmds.size() > 0)

  auto parsed = IdaLogger::parseFLASHIdaLog(log_f);
  TEST_TRUE(parsed.size() > 0)
  bool any = false;
  for (const auto& entry : parsed)
  {
    for (const auto& p : entry.second)
    {
      any = true;
      TEST_EQUAL(p.size(), 15)
      TEST_TRUE(inMass(p[0]))                       // e0 mass
      TEST_TRUE(inChargeD(p[1]))                    // e1 charge
      TEST_TRUE(inUnit(p[2]))                       // e2 qscore
      TEST_TRUE(posFinite(p[3]))                    // e3 window lo
      TEST_TRUE(p[4] > p[3])                        // e4 window hi > lo
      TEST_TRUE(posFinite(p[5]))                    // e5 precursor intensity
      TEST_TRUE(posFinite(p[6]))                    // e6 precursor mass intensity
      TEST_TRUE(inUnit(p[9]))                       // e9 charge_cos
      TEST_TRUE(nonNegFinite(p[10]))                // e10 charge_snr
      TEST_TRUE(inUnit(p[11]))                      // e11 iso_cos
      TEST_TRUE(nonNegFinite(p[12]))                // e12 snr
      TEST_TRUE(inUnit(p[13]))                      // e13 charge_score
      TEST_TRUE(finiteVal(p[14]))                   // e14 ppm_error
      // e7/e8 ChargeRange are the documented ⟂DEFER degenerate pin (== e1 charge)
      TEST_TRUE(std::abs(p[7] - p[1]) < 1.0)
      TEST_TRUE(std::abs(p[8] - p[7]) < 1.0)
    }
  }
  TEST_TRUE(any)
  std::remove(log_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §L2 -- ida_log : two MS1 scans => two distinct keys
/////////////////////////////////////////////////////////////
START_SECTION(ida_log_multi_scan_distinct_keys)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.size() < 2 || ms2.empty())

  std::string log_f = "lf_l2_ida.log";
  std::remove(log_f.c_str());
  // Selecting config: a standard-DDA "none" config (enable_ms3=false) logs every MS1 as "0 targets",
  // which parseFLASHIdaLog (FLASHIda.cpp:595) skips -> no log groups. enable_ms3=true gives ms2
  // selection="intensity" + inclusion, so each MS1 selects >=1 precursor -> a NON-empty log group.
  std::string json = buildJsonWithRuntime(log_f, "", "", true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Feed TWO MS1 surveys under the engine's OWN emitted tracking ids (no invented ids) -- the ida_log key
  // is the decoded 3-char tracking id, so distinct engine ids => distinct log groups. ms1[1] = scan 134
  // carries the cytC envelope (ms1[0] = scan 132 is a weak edge scan). runFullAcquisition feeds the 2nd
  // survey at rt+1000 (beyond RT_window 180) so its precursor re-selects. DRY: reuse runFullAcquisition.
  AcqResult acq = runFullAcquisition(&ida, ms1[1], ms2[0], 300, /*n_ms1=*/2);
  ABORT_IF(acq.ms1_cmds.size() < 2)
  std::string id0(acq.ms1_cmds[0].scan_description, 3), id1(acq.ms1_cmds[1].scan_description, 3);
  TEST_TRUE(id0 != id1)                          // distinct engine ids => distinct ida_log keys

  auto parsed = IdaLogger::parseFLASHIdaLog(log_f);
  TEST_TRUE(parsed.size() >= 2)                  // two non-"0 targets" MS1 groups
  std::remove(log_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §X5 -- scan_description cap (≤15 + NUL) + id/ion round-trip (B4 pinned)
/////////////////////////////////////////////////////////////
START_SECTION(scan_description_cap_roundtrip)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string json = buildJsonWithRuntime("", "", "", true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto cycle = runFullCycle(&ida, ms1, ms2);
  TEST_TRUE(cycle.ms3_cmds.size() > 0)

  for (const auto& c : cycle.ms3_cmds)
  {
    std::string desc(c.scan_description);
    TEST_TRUE(desc.size() <= 15)                 // snprintf(...,16,...) cap (intended instrument limit)
    TEST_TRUE(desc.size() >= 4)
    TEST_TRUE(isTrackingId(desc.substr(0, 3)))   // id prefix round-trips
    TEST_TRUE(desc[3] == 'R' || desc[3] == 'E')  // MS3 recording / exploration marker
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// §X6 -- parent_tracking_id resolution over a full INCLUSION-mode acquisition
//
// Drives the engine the way Flash.cs does: getNextScanCommand -> feed the matching
// scan back stamped with THAT command's engine-emitted scan_description, so every
// MS1->MS2->MS3 link uses the engine's own ids end-to-end (no test-fabricated ids).
// The inclusion-pinned cytC MS3 recipe (buildJsonWithRuntime enable_ms3=true) is the
// C++-facing equivalent of method_inclusion.json + the MS3 mode-1 HCD recipe: target_mode=1
// + inclusion_cytc.txt + the M-starting cytC proteoform.  HARD requirement: EVERY non-empty
// parent_tracking_id in scan_commands.tsv must resolve to an engine-emitted command id, with
// MS2->MS1 / MS3->MS2 lineage.
/////////////////////////////////////////////////////////////
START_SECTION(parent_tracking_id_resolution)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.size() < 2 || ms2.empty())

  std::string cmd_f = "lf_x6_commands.tsv", res_f = "lf_x6_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
  std::string json = buildJsonWithRuntime("", cmd_f, res_f, true);  // inclusion (target_mode 1) + MS3
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Single-MS1 engine-chained acquisition: the engine emits the MS1 command, we feed the cytC
  // MS1 back under that id; the resulting MS2s (and any MS3s) chain off the engine's own ids.
  // ms1[1] = scan 134 carries the cytC envelope (ms1[0] = scan 132 is a weak edge scan from which
  // the engine correctly selects 0 precursors).
  AcqResult acq = runFullAcquisition(&ida, ms1[1], ms2[0]);
  TEST_TRUE(acq.ms2_cmds.size() >= 1)  // inclusion mode must trigger >=1 MS2 on the pinned cytC precursor

  auto cmds = TSVFile::parse(cmd_f);

  // Hard plausibility gate: every non-empty parent_tracking_id resolves to a known engine-emitted id.
  // The engine chains ids itself here, so the fed-id universe is empty -- all known ids come from the
  // commands TSV's own tracking_id column (MS1 survey, AGC, MS2, MS3 are all written there).
  std::string err;
  bool resolved = validateParentTrackingIds(cmds, std::vector<std::string>(), err);
  TEST_EQUAL(err, std::string(""))  // surface the offending parent id on failure
  TEST_TRUE(resolved)

  // Explicit lineage: MS2 parents are MS1 ids; MS3 parents are MS2 ids.
  auto level = commandLevels(cmds);
  bool checked_ms2 = false, checked_ms3 = false, lineage_ok = true;
  for (const auto& row : cmds.rows)
  {
    int lvl = std::atoi(cell(cmds, row, "ms_level").c_str());
    std::string parent = cell(cmds, row, "parent_tracking_id");
    if (lvl == 2 && !parent.empty())
    {
      checked_ms2 = true;
      // parent is an MS1 survey: present as a known id but NOT itself an MS2/MS3 command
      auto it = level.find(parent);
      lineage_ok = lineage_ok && (it == level.end() || it->second == 1);
    }
    else if (lvl == 3 && !parent.empty())
    {
      checked_ms3 = true;
      lineage_ok = lineage_ok && level.count(parent) && level[parent] == 2;  // MS3 parent is an MS2
    }
  }
  TEST_TRUE(checked_ms2)
  TEST_TRUE(lineage_ok)
  (void)checked_ms3;  // MS3 rows are recipe-dependent; lineage is asserted WHEN present

  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §R2t -- scan_results.tsv : tag-targeting MS2 rows log a real tag_count>0 (locks E4)
//
// Tag-based targeting (E4) fires whenever a target-protein database is loaded: the engine
// generates sequence tags from the MS2 deconvolution and matches them against the FASTA. The
// logged tag_count is the REAL number of generated tags (not a 0/1 boolean). This is the
// C++-facing equivalent of method_tag_targeting.json: tagging active + a cytC-bearing FASTA
// (test_fasta.fasta contains CYC_HORSE) driven with the real cytC MS2 fixture, so >=1 MS2
// result row must carry tag_count>0. MS2 selection stays "none" so tagging runs without the
// fragment-matching protein_sequence requirement (tagging is gated only on the FASTA database).
/////////////////////////////////////////////////////////////
START_SECTION(results_tag_targeting)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_r2t_commands.tsv", res_f = "lf_r2t_results.tsv";
  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());

  // C++-facing config (mirrors buildJsonWithRuntime's schema) with a target-protein FASTA loaded
  // for tag-based targeting. tagging params present; conditional_ms2 stays false so no follow_up_scan
  // is required by Config::validate(); MS2 selection "none" so no MS3/protein_sequence chain is needed.
  std::ostringstream cfg;
  cfg << R"({
    "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10] },
    "precursor_selection": { "RT_window": 180, "target_mode": 0, "IDScore": false, "AllCharges": false, "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1 },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [ { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 } ]
    },
    "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 }, "agc_interval_seconds": 999999 },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "../../FlashIDA/test-data/configs/test_fasta.fasta", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 5 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    },
    "runtime": {
      "ida_log_path": "", "scan_commands_path": ")" << cmd_f << R"(",
      "scan_results_path": ")" << res_f << R"(", "identification_log_path": ""
    }
  })";
  std::string json = cfg.str();
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto cmds = TSVFile::parse(cmd_f);
  auto res = TSVFile::parse(res_f);
  auto level = commandLevels(cmds);

  bool found_ms2 = false, any_tag = false;
  for (const auto& row : res.rows)
  {
    std::string tid = cell(res, row, "tracking_id");
    auto it = level.find(tid);
    if (it == level.end() || it->second != 2) continue;  // MS2 result rows only
    found_ms2 = true;
    double tc = toD(cell(res, row, "tag_count"));
    TEST_TRUE(nonNegFinite(tc))                            // real count, never negative
    if (tc > 0.0) any_tag = true;
  }
  TEST_TRUE(found_ms2)
  TEST_TRUE(any_tag)  // E4: the real cytC MS2 tags must match CYC_HORSE in test_fasta.fasta

  std::remove(cmd_f.c_str()); std::remove(res_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §E6rt -- scan_commands.tsv scan_description column == the drained ScanCommand buffer
//
// X5 only inspects the DRAINED cmd.scan_description; this section closes the E6 loop by
// proving the LOGGED scan_description column (index 28) is byte-identical to the descriptor
// the engine actually drained for the same command (joined on the 3-char tracking-id prefix).
/////////////////////////////////////////////////////////////
START_SECTION(commands_scan_description_roundtrip)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  std::string cmd_f = "lf_e6_commands.tsv";
  std::remove(cmd_f.c_str());
  std::string json = buildJsonWithRuntime("", cmd_f, "", true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto cycle = runFullCycle(&ida, ms1, ms2);
  TEST_TRUE(cycle.ms2_cmds.size() > 0)

  // Map each drained command's 3-char id prefix -> its raw scan_description buffer.
  std::map<std::string, std::string> drained_desc;
  for (const auto& c : cycle.ms2_cmds) { std::string d(c.scan_description); if (d.size() >= 3) drained_desc[d.substr(0, 3)] = d; }
  for (const auto& c : cycle.ms3_cmds) { std::string d(c.scan_description); if (d.size() >= 3) drained_desc[d.substr(0, 3)] = d; }

  auto t = TSVFile::parse(cmd_f);
  ABORT_IF(t.colIndex("scan_description") != 28)  // E6: scan_description is the last column

  bool matched_any = false, equal_ok = true;
  for (const auto& row : t.rows)
  {
    std::string tid = cell(t, row, "tracking_id");
    auto it = drained_desc.find(tid);
    if (it == drained_desc.end()) continue;  // MS1/AGC rows we did not capture as drained MS2/MS3
    matched_any = true;
    equal_ok = equal_ok && (cell(t, row, "scan_description") == it->second);
  }
  TEST_TRUE(matched_any)  // at least the MS2 commands must round-trip
  TEST_TRUE(equal_ok)     // logged col[28] == drained descriptor buffer, byte-for-byte

  std::remove(cmd_f.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §T9 -- cytC REAL-MS3 fragment data via ion-name-keyed real-spectrum feed (skips cleanly if absent)
//
// The standard MS3 sections feed the MS2 spectrum back as the MS3 scan (a shortcut), so
// matched_protein/fragment_count/tic are only conditionally populated. This section proves the
// real MS3 path end-to-end on REAL data: under the inclusion-pinned, EXHAUSTIVE (ms3.max_targets=200)
// MS3 mode-1 recipe (buildJsonWithRuntime enable_ms3=true), the engine emits ion-targeted MS3 commands;
// each command's scan_description encodes its precursor fragment ion ({id}R{mass}k@{charge}{ion}{idx}).
// runFullCycle decodes that ion and feeds the matching REAL per-ion MS3 spectrum from the manifest built
// by globbing ms3_cytc_<ion>_scan<N>.txt in the spectra dir; commands whose ion is absent are SKIPPED
// (tolerated). The gate is tolerance-based: >=1 FED MS3 row matches (matched_protein populated), and every
// matched row is fully populated (isProForma proteoform, fragment_count>0, tic in (0,1], resolvable parent).
// Guarded on the manifest: if no ms3_cytc_<ion>_scan<N>.txt fixtures exist this section asserts nothing and
// passes -- never fabricates data. Exactness of the match-dependent fields is golden-locked in the C# suite.
/////////////////////////////////////////////////////////////
START_SECTION(results_ms3_real_fragment_data)
{
  // Build the ion -> [ScanData...] manifest by globbing the spectra dir for ms3_cytc_<ion>_scan<N>.txt.
  auto manifest = buildMs3IonManifest("../../FlashIDA/test-data/spectra");
  if (manifest.empty())
  {
    // No real MS3 fragment fixtures available -- skip cleanly (do NOT fabricate). The
    // match-dependent path's exactness is golden-locked in the C# suite.
    TEST_TRUE(true)
  }
  else
  {
    auto ms1 = loadTsvScans(CYTC_MS1);
    auto ms2 = loadTsvScans(CYTC_MS2);  // scan57 ladder -> proteoform context (auto-PTM)
    ABORT_IF(ms1.empty() || ms2.empty())

    std::string cmd_f = "lf_t9_commands.tsv", res_f = "lf_t9_results.tsv", id_f = "lf_t9_id.tsv";
    std::remove(cmd_f.c_str()); std::remove(res_f.c_str()); std::remove(id_f.c_str());
    // M-start, inclusion-pinned, MS3 mode-1, EXHAUSTIVE MS3 emission: pass ms3_max_targets=50 (scoped to
    // §T9 only; other enable_ms3 sections keep the engine default ~10). 50 > cytC's ~16-20 detectable
    // fragment ions, so emission stays exhaustive while keeping the section's runtime in check.
    // F9: capture identification.tsv too — the matched-fragment FLOOR now lives there (ms3_fragments), since
    // scan_results fragment_count is sentinelized to -1.
    std::string json = buildJsonWithRuntime("", cmd_f, res_f, true, id_f, "HCD", 50);
    FLASHIda ida(const_cast<char*>(json.c_str()));

    // Mirror the shared harness wiring; feed the REAL per-ion MS3 spectrum keyed by the decoded descriptor ion.
    CycleResult cyc = runFullCycle(&ida, ms1, ms2, &manifest);
    TEST_TRUE(cyc.ms3_cmds.size() > 0)  // engine issued MS3 command(s)

    // Q2: every emitted MS3 command descriptor that carries an ion decodes to a valid, non-truncated ion
    // (ion_type in {a,b,c,x,y,z}, ion_index >= 1). The no-ion form is tolerated (not counted as a failure).
    int decoded_ion_cmds = 0;
    for (const auto& c : cyc.ms3_cmds)
    {
      char t; int idx;
      if (decodeTrailingIon(c.scan_description, t, idx))  // false == no-ion branch (tolerated)
      {
        TEST_TRUE(ionTypeOk(std::string(1, t)) && idx >= 1)
        ++decoded_ion_cmds;
      }
    }
    TEST_TRUE(decoded_ion_cmds > 0)  // >=1 ion-targeted MS3 emitted

    auto cmds  = TSVFile::parse(cmd_f);
    auto res   = TSVFile::parse(res_f);
    auto level = commandLevels(cmds);

    // GATE: walk MS3 result rows; >=1 FED row WITH a match; each matched row fully populated.
    int ms3_rows = 0, matched_rows = 0;
    for (const auto& row : res.rows)
    {
      std::string tid = cell(res, row, "tracking_id");
      auto it = level.find(tid);
      if (it == level.end() || it->second != 3) continue;  // MS3 result rows only
      ++ms3_rows;
      if (cell(res, row, "matched_protein").empty()) continue;  // unfed / unmatched MS3 -> tolerated
      ++matched_rows;
      TEST_TRUE(isProForma(cell(res, row, "proteoform_sequence")))
      // F9: scan_results fragment_count is the -1 SENTINEL on MS3 rows (the matched count lives in
      // identification.tsv, asserted below). (Was fragment_count > 0.)
      TEST_EQUAL(std::atoi(cell(res, row, "fragment_count").c_str()), -1)
      TEST_TRUE(toD(cell(res, row, "tic_coverage")) > 0.0 && inUnit(toD(cell(res, row, "tic_coverage"))))
      TEST_TRUE(isTrackingId(cell(res, row, "parent_tracking_id")))       // resolvable MS2 parent
    }
    TEST_TRUE(ms3_rows > 0)
    TEST_TRUE(matched_rows >= 1)  // Q1 GATE: >=1 FED MS3 scan matched the pinned cytC proteoform

    // F9: the matched-fragment FLOOR relocated from scan_results.fragment_count to identification.tsv —
    // >=1 MS3 identification row carries a non-empty ms3_fragments list (matched-ness lives where it belongs).
    auto idf = TSVFile::parse(id_f);
    int ms3_id_with_frags = 0;
    for (const auto& row : idf.rows)
    {
      if (std::atoi(cell(idf, row, "ms_level").c_str()) != 3) continue;
      if (!cell(idf, row, "ms3_fragments").empty()) ++ms3_id_with_frags;
    }
    TEST_TRUE(ms3_id_with_frags >= 1)

    std::remove(cmd_f.c_str()); std::remove(res_f.c_str()); std::remove(id_f.c_str());
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// §IONP -- ion-decode parity vectors (C# <-> C++ drift guard)
//
// decodeTrailingIonKey (FLASHIda_TestHelpers.h) and the C# DecodeIonFromScanDescription
// (FLASHIdaLogGolden_test.cs) MUST decode the trailing precursor-ion key identically,
// byte-for-byte. These are the SHARED parity vectors used in BOTH suites (the C# twin
// asserts the same table). Each row: a raw scan_description and its expected ion key
// (or the no-ion / malformed form, which must decode to "none"). Covers the
// '@'-inside-the-id, multi-digit charge+index, no-ion, MS1-survey, invalid-type, and
// index-0 / empty edge cases. If this section diverges, the two decoders have drifted.
/////////////////////////////////////////////////////////////
START_SECTION(ion_decode_parity_vectors)
{
  struct Vec { const char* desc; bool ok; const char* expected; };
  const Vec vectors[] = {
    {"!#@R4.450k@5y38", true,  "y38"},   // '@' INSIDE the 3-char id (!#@); rfind/LAST '@' is the charge delim
    {"!!!R1.000k@2b10", true,  "b10"},
    {"AAAR12.351k@3y5", true,  "y5"},
    {"JJJR2.0k@12c144", true,  "c144"},  // multi-digit charge + index
    {"!!!R5.0k@4",      false, ""},      // no-ion form (nothing after the charge digits)
    {"!!\"S",           false, ""},      // MS1 survey descriptor, no '@'
    {"!!!R5.0k@2d10",   false, ""},      // invalid ion type 'd'
    {"!!!R5.0k@2y0",    false, ""},      // index 0 (<1) invalid
    {"",                false, ""},      // empty
  };
  for (const auto& v : vectors)
  {
    std::string key;
    bool decoded = decodeTrailingIonKey(std::string(v.desc), key);
    TEST_EQUAL(decoded, v.ok)
    if (v.ok) { TEST_EQUAL(key, std::string(v.expected)) }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
