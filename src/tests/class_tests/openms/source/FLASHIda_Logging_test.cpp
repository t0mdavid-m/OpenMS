// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Phase 10 unit tests: IDA logging and scan tracking TSV files.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandJoin.h>

#include <cmath>     // std::abs on the parsed float masses
#include <iomanip>   // fixed4 in scan_commands_mono_mass_is_written_at_four_decimals
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cstdio>
#include <iostream>  // the two lock-split guards below report their counts
#include <thread>    // concurrent_drain_writes_one_wellformed_row_per_call
#include <map>       // precursor_id inheritance check
#include <set>       // tracking-id uniqueness check
#include <utility>   // std::pair

#include "FLASHIda_TestHelpers.h"
#include "FLASHIda_TestAccess.h"  // the target-log loader's maps (PrecursorSelection is a friend)

using namespace OpenMS;

namespace
{
  // Spectrum paths for these logging tests. The shared loaders, cycle drivers, and
  // JSON config builder now live in FLASHIda_TestHelpers.h (included above).
  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_standard.txt";
  const std::string ms2_tsv_path = "../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt";

  // ---- Fixture for the ADR-0046 sections (scan_commands_parse_* / scan_commands_join_*) ----------
  //
  // REAL data, not invented: the header and three rows of a run folder's own scan_commands.tsv
  // (Eclipse, 2026-09-15 -- a PRE-ADR-0046 file, mono_mass at six significant digits), with the
  // spectra facts from that run's mzML: survey !!= is scan 32, survey !!> is scan 33, and MS2 !!B is
  // scan 37 isolating 531.790649. A newer survey really does sit between the command's survey and
  // its scan -- which is the whole reason the join exists.
  using Cells = std::vector<std::string>;

  const Cells kCmdHeader = {"tracking_id", "scan_type", "ms_level", "parent_tracking_id", "precursor_id", "priority",
                            "mono_mass", "charge", "precursor_mz", "isolation_width", "qscore", "charge_cos",
                            "charge_snr", "iso_cos", "snr", "charge_score", "activation", "collision_energy",
                            "hcd_energy", "reaction_time", "reagent_max_it", "reagent_agc_target", "ppm_error",
                            "precursor_intensity", "peakgroup_intensity", "ion_type", "ion_index", "ms3_proteoform",
                            "scan_description", "faims_cv", "faims_enabled", "first_mass", "last_mass", "enqueue_ts"};

  Cells surveyRow(const std::string& id, const std::string& enqueue_ts)
  {
    return {id, "survey", "1", "", "0", "3",
            "0", "0", "0", "0", "0", "0", "0", "0", "0", "0",   // mono_mass .. charge_score: stage-less placeholders
            "none",
            "0", "0", "0", "0", "0", "0", "0", "0",             // collision_energy .. peakgroup_intensity
            "", "0", "", id + "S", "0", "0", "500.0000", "2000.0000", enqueue_ts};
  }

  const Cells kS32 = surveyRow("!!=", "267351474");
  const Cells kS33 = surveyRow("!!>", "267351740");
  const Cells kM37 = {"!!B", "recording", "2", "!!=", "1", "2", "5305.33", "10", "531.791", "0.901196", "0.443474",
                      "0.738412", "5.70422", "0.859759", "3.90517", "1", "HCD", "29", "29", "0", "0", "0", "0.608297",
                      "2134.16", "2947.43", "", "0", "", "!!BR5.30533k@10", "0", "0", "200.0000", "2000.0000",
                      "267352272"};

  // @p row with one named cell replaced.
  Cells withCell(Cells row, const std::string& column, const std::string& value)
  {
    for (size_t i = 0; i < kCmdHeader.size(); ++i)
      if (kCmdHeader[i] == column) row[i] = value;
    return row;
  }

  // A scan_commands.tsv written by hand into a fresh dir; returns its path.
  std::string writeCommandsFile(const std::string& tag, const std::vector<Cells>& lines)
  {
    const std::string path = freshLogDir(tag) + "/scan_commands.tsv";
    std::ofstream f(path);
    for (const auto& cells : lines)
    {
      for (size_t i = 0; i < cells.size(); ++i) f << (i ? "\t" : "") << cells[i];
      f << "\n";
    }
    return path;
  }

  ScanCommandJoin::Scan scanOf(int number, int level, const std::string& description, const std::vector<double>& targets = {})
  {
    ScanCommandJoin::Scan s;
    s.scan_number = number;
    s.ms_level = level;
    s.scan_description = description;
    s.isolation_targets = targets;
    return s;
  }

  // The three spectra the three fixture rows belong to.
  std::vector<ScanCommandJoin::Scan> realScans()
  {
    return {scanOf(32, 1, "!!=S"), scanOf(33, 1, "!!>S"), scanOf(37, 2, "!!BR5.30533k@10", {531.790649})};
  }
}

START_TEST(FLASHIda_Logging, "$Id$")

/////////////////////////////////////////////////////////////

// deconvolveMS1 writes nothing to stdout.
//
// It used to print every mass of every survey, then std::endl -- one flush per survey under
// analysis_mutex_, the console pathology of ADR-0025 / ADR-0033, for masses ida.log (AllMass=) and
// scan_results (deconv_masses) already record. Not vacuous on any input: the old code wrote at least
// that newline on every call, whatever the spectrum held. No fixture: one synthetic 10 kDa species.
START_SECTION(deconvolveMS1_writes_nothing_to_stdout)
{
  Config cfg{buildJsonWithLogDir("")};
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});

  // Charges 12 down to 8, so the m/z come out ascending; a Poisson envelope (lambda 5.4) per charge.
  std::vector<double> mzs, ints;
  const double mass = 10000.0, lambda = 5.4;
  for (int z = 12; z >= 8; --z)
  {
    double p = std::exp(-lambda);
    for (int k = 0; k < 14; ++k)
    {
      if (k > 0) { p *= lambda / k; }
      mzs.push_back((mass + k * Constants::ISOTOPE_MASSDIFF_55K_U + z * Constants::PROTON_MASS_U) / z);
      ints.push_back(1e6 * p);
    }
  }

  std::ostringstream captured;
  std::streambuf* const real_cout = std::cout.rdbuf(captured.rdbuf());
  try
  {
    deconv.deconvolveMS1(mzs.data(), ints.data(), static_cast<int>(mzs.size()), 60.0, 0.0);
  }
  catch (...)
  {
    std::cout.rdbuf(real_cout);
    throw;
  }
  std::cout.rdbuf(real_cout);
  TEST_EQUAL(captured.str(), std::string())
}
END_SECTION

// Test 0: a >= 10 target count must not be mistaken for "0 targets"
//
// The skip test was `line.find("0 targets")`, an unanchored SUBSTRING, so it also matched
// "10 targets" and "20 targets". Those headers were dropped and their Mass= rows then inherited the
// PREVIOUS header's scan key -- or, at the head of a file, an uninitialised int. This is not
// hypothetical: the committed separate_charges golden carries eight "- 10 targets" headers and one
// "- 20 targets", because command fan-out (per MS2 config x CE-sweep variant x charge under
// `separate`) makes double-digit counts routine.
//
// Driven off a hand-written log rather than an acquisition: the defect is in the READER, and a
// synthetic file states the exact input that triggers it without needing a config that happens to
// fan out past ten.
START_SECTION(ida_log_double_digit_target_count_is_not_skipped)
{
  const std::string dir = freshLogDir("logging_ida_double_digit_targets");
  const std::string path = dir + "/ida.log";
  {
    std::ofstream f(path);
    // Entry A: 10 targets. Must be kept, under its own key.
    f << "MS1 Scan# 7 RT 1.0000 (Access ID !!\") - 10 targets\n";
    f << "Mass=12351.3933\tZ=15\tScore=0.90000\tWindow=[824.0356-825.8985]"
         "\tPrecursorIntensity=1.00000\tPrecursorMassIntensity=2.00000"
         "\tFeatures=[0.9,1.0,0.9,1.0,0.9,1.0]\tChargeRange=[12-19]\tHCD=0\n";
    f << "AllMass=12351.3933\n";
    // Entry B: genuinely 0 targets. Must still be skipped -- the anchor must not over-correct.
    f << "MS1 Scan# 8 RT 2.0000 (Access ID !!#) - 0 targets\n";
    f << "AllMass=\n";
    // Entry C: 20 targets, the other count the substring swallowed.
    f << "MS1 Scan# 9 RT 3.0000 (Access ID !!$) - 20 targets\n";
    f << "Mass=8604.9191\tZ=8\tScore=0.80000\tWindow=[1076.6249-1077.6249]"
         "\tPrecursorIntensity=3.00000\tPrecursorMassIntensity=4.00000"
         "\tFeatures=[0.8,1.0,0.8,1.0,0.8,1.0]\tChargeRange=[6-10]\tHCD=0\n";
    f << "AllMass=8604.9191\n";
  }

  auto parsed = IdaLogger::parseFLASHIdaLog(path);

  TEST_EQUAL(parsed.count(7), 1)          // the "10 targets" header survived
  TEST_EQUAL(parsed.count(9), 1)          // and the "20 targets" one
  TEST_EQUAL(parsed.count(8), 0)          // while a real "0 targets" entry is still skipped
  TEST_EQUAL(parsed.size(), 2)

  // Each row landed under ITS OWN header, not the previous one. Before the fix scan 7's row went to
  // an uninitialised key and scan 9's was appended to whatever came before it.
  TEST_EQUAL(parsed[7].size(), 1)
  TEST_EQUAL(parsed[9].size(), 1)
  TEST_TRUE(std::abs(parsed[7][0][0] - 12351.3933) < 0.001)
  TEST_TRUE(std::abs(parsed[9][0][0] - 8604.9191) < 0.001)
}
END_SECTION

// The OTHER reader of this grammar: PrecursorSelection's target-log loader (targeting.target_log_files),
// which files every AllMass= line under the retention time of the MS1 header above it. It used to skip
// "- 0 targets" headers before reading their RT, copying parseFLASHIdaLog's skip -- but a survey with no
// target still writes its AllMass= line, and exclusion_masses (mode 3) reads exactly those lines. So the
// masses of every target-less survey were filed under the PREVIOUS header's RT, and those ahead of the
// first header with targets under 0.0, which filterAndRank applies at EVERY RT. MS1-only replays of the
// Eclip lysates (2026-10-20) filed 17-41 % of the non-empty AllMass= lines of their logs under another
// survey's RT, up to 92 s away.
//
// Driven off a hand-written log for the same reason as the section above: the defect is in the reader.
START_SECTION(target_log_loader_files_allmass_under_its_own_survey)
{
  const std::string dir = freshLogDir("logging_target_log_loader_rt");
  const std::string path = dir + "/ida.log";
  {
    std::ofstream f(path);
    // Entry A, 0.5 min: no target, ahead of the first header that has one. Was filed under rt 0.
    f << "MS1 Scan# 1 RT 0.5000 (Access ID !!\") - 0 targets\n";
    f << "AllMass=5000.0000 6000.0000\n";
    // Entry B, 1 min: one target.
    f << "MS1 Scan# 2 RT 1.0000 (Access ID !!#) - 1 targets\n";
    f << "Mass=12351.3933\tZ=15\tScore=0.90000\tWindow=[824.0356-825.8985]"
         "\tPrecursorIntensity=1.00000\tPrecursorMassIntensity=2.00000"
         "\tFeatures=[0.9,1.0,0.9,1.0,0.9,1.0]\tChargeRange=[12-19]\tHCD=0\n";
    f << "AllMass=12351.3933 7000.0000\n";
    // Entry C, 9 min: no target. Was filed under B's 1 min, eight minutes from where it was seen.
    f << "MS1 Scan# 3 RT 9.0000 (Access ID !!$) - 0 targets\n";
    f << "AllMass=8604.9191\n";
  }

  // The loader runs in PrecursorSelection's constructor; the targeting mode decides which lines it keeps.
  auto loaderConfig = [&](const std::string& targeting) {
    std::string cfg = buildJsonWithLogDir(dir);
    auto rep = [&cfg](const std::string& from, const std::string& to) {
      auto p = cfg.find(from);
      TEST_EQUAL(p != std::string::npos, true)   // a missing literal would leave the config silently unchanged
      if (p != std::string::npos) cfg.replace(p, from.size(), to);
    };
    rep("\"targeting\": \"none\"", "\"targeting\": \"" + targeting + "\"");
    rep("\"target_logs\": []", "\"target_logs\": [\"" + path + "\"]");
    return cfg;
  };

  {
    Config config(loaderConfig("exclusion_masses"));
    Deconvolution deconv(config, config.toleranceList());
    PrecursorSelection selection(config, deconv);
    // survey rt in SECONDS (the log writes minutes, the loader multiplies by 60) -> that survey's masses
    const auto excl = FLASHIdaTestAccess::exclusionRtMasses(selection);

    TEST_EQUAL(excl.size(), 3)
    TEST_EQUAL(excl.count(0.0), 0)       // nothing is excluded at every retention time
    TEST_EQUAL(excl.count(30.0), 1)      // A under its own 0.5 min
    TEST_EQUAL(excl.count(60.0), 1)
    TEST_EQUAL(excl.count(540.0), 1)     // C under its own 9 min
    if (excl.count(30.0) == 1) { TEST_EQUAL(excl.at(30.0).size(), 2) }
    if (excl.count(60.0) == 1) { TEST_EQUAL(excl.at(60.0).size(), 2) }  // B's own two, C's no longer appended
    if (excl.count(540.0) == 1)
    {
      TEST_EQUAL(excl.at(540.0).size(), 1)
      TEST_REAL_SIMILAR(excl.at(540.0)[0], 8604.9191)
    }
  }
  {
    // in_depth (mode 2) reads Mass= lines only, and a 0-target survey writes none: its load is unchanged.
    Config config(loaderConfig("in_depth"));
    Deconvolution deconv(config, config.toleranceList());
    PrecursorSelection selection(config, deconv);
    const auto targets = FLASHIdaTestAccess::targetMassRts(selection);

    TEST_EQUAL(targets.size(), 1)
    if (targets.size() == 1)
    {
      TEST_REAL_SIMILAR(targets.begin()->first, 12351.3933)
      TEST_EQUAL(targets.begin()->second.size(), 1)
      if (targets.begin()->second.size() == 1) { TEST_REAL_SIMILAR(targets.begin()->second[0], 60.0) }
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////

// Test 1: IDA Log contract -- write + parseFLASHIdaLog roundtrip
START_SECTION(ida_log_contract_roundtrip)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  // Own (wiped + created) log folder for this engine; the IDA log is dir + "/ida.log"
  const std::string dir = freshLogDir("logging_ida_log_contract");

  std::string json = buildJsonWithLogDir(dir);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Drive the MS1 surveys via the canonical interleaved driver: the engine emits each survey
  // command and we feed the next ms1_standard scan back stamped with the engine's OWN tracking id
  // (was pushAllScans' fabricated encode(800000+i) ids, which the always-on MS1 gate now rejects ->
  // 0 precursors). No MS2 fixture needed: the ida_log is written at MS1 time and ms2 selection is
  // "none" here, so MS2 commands are recorded but not fed. The engine selecting >=1 precursor per
  // MS1 (ms1_standard) yields >=1 emitted MS2 command -- the faithful analog of the old
  // total_commands>0 (which counted MS2 commands pushed during MS1 processing).
  AcqResult acq = runInterleaved(&ida, ms1_scans, std::vector<ScanData>{});
  TEST_TRUE(acq.ms2_cmds.size() > 0);

  // Parse the IDA log back using parseFLASHIdaLog
  auto parsed = IdaLogger::parseFLASHIdaLog(dir + "/ida.log");

  // Verify: at least one scan group with precursors
  TEST_TRUE(parsed.size() > 0);

  // Verify each precursor has exactly 15 floats
  for (const auto& entry : parsed)
  {
    for (const auto& precursor : entry.second)
    {
      TEST_EQUAL(precursor.size(), 15);
      // mass (index 0) should be > 0
      TEST_TRUE(precursor[0] > 0);
      // charge (index 1) should be >= 4 (min_charge in config)
      TEST_TRUE(precursor[1] >= 4);
      // qscore (index 2) should be >= 0
      TEST_TRUE(precursor[2] >= 0);
      // window (indices 3,4) should be > 0
      TEST_TRUE(precursor[3] > 0);
      TEST_TRUE(precursor[4] > precursor[3]);
    }
  }
}
END_SECTION

// Test 2: ScanCommands TSV -- full MS1->MS2->MS3 cycle
START_SECTION(scan_commands_tsv_format)
{
  // Real CytC MS1+MS2 so the MS3 path actually fires (MS3 enabled in the config below);
  // generic data deconvolves to no proteoform-matching fragments and emits zero MS3.
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  const std::string dir = freshLogDir("logging_scan_commands");

  // Enable MS3 so we get MS3 commands in the TSV
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Full MS1->MS2->MS3 cycle
  auto cycle = runFullCycle(&ida, ms1_scans, ms2_scans);
  TEST_TRUE(cycle.ms2_cmds.size() > 0);

  // Parse and verify TSV
  auto tsv = TSVFile::parse(dir + "/scan_commands.tsv");

  // Header check
  TEST_TRUE(tsv.colIndex("tracking_id") >= 0);
  TEST_TRUE(tsv.colIndex("ms_level") >= 0);
  TEST_TRUE(tsv.colIndex("scan_type") >= 0);
  TEST_TRUE(tsv.colIndex("enqueue_ts") >= 0);
  TEST_TRUE(tsv.colIndex("qscore") >= 0);
  TEST_TRUE(tsv.colIndex("ion_type") >= 0);
  TEST_TRUE(tsv.colIndex("ion_index") >= 0);

  int ms_level_col = tsv.colIndex("ms_level");
  int charge_col = tsv.colIndex("charge");
  int activation_col = tsv.colIndex("activation");
  int precursor_mz_col = tsv.colIndex("precursor_mz");
  int iso_width_col = tsv.colIndex("isolation_width");
  int col_energy_col = tsv.colIndex("collision_energy");

  // Fail closed on a dropped/renamed format column: these must exist, otherwise the per-row
  // semicolon-format checks below would silently no-op (colIndex returns -1) and a schema
  // regression would pass unnoticed. (The per-row `< row.size()` guards remain for bounds.)
  TEST_TRUE(ms_level_col >= 0);
  TEST_TRUE(charge_col >= 0);
  TEST_TRUE(activation_col >= 0);
  TEST_TRUE(precursor_mz_col >= 0);
  TEST_TRUE(iso_width_col >= 0);
  TEST_TRUE(col_energy_col >= 0);

  bool found_ms2 = false;
  bool found_ms3 = false;
  for (const auto& row : tsv.rows)
  {
    if (ms_level_col < 0 || ms_level_col >= (int)row.size())
      continue;

    if (row[ms_level_col] == "2")
    {
      found_ms2 = true;
      // MS2 rows: single stage, no semicolons
      if (charge_col >= 0 && charge_col < (int)row.size())
        TEST_TRUE(row[charge_col].find(';') == std::string::npos);
      if (activation_col >= 0 && activation_col < (int)row.size())
        TEST_TRUE(row[activation_col].find(';') == std::string::npos);
    }

    if (row[ms_level_col] == "3")
    {
      found_ms3 = true;
      // MS3 rows: two stages, semicolons present
      if (charge_col >= 0 && charge_col < (int)row.size())
        TEST_TRUE(row[charge_col].find(';') != std::string::npos);
      if (activation_col >= 0 && activation_col < (int)row.size())
        TEST_TRUE(row[activation_col].find(';') != std::string::npos);
      if (precursor_mz_col >= 0 && precursor_mz_col < (int)row.size())
        TEST_TRUE(row[precursor_mz_col].find(';') != std::string::npos);
      if (iso_width_col >= 0 && iso_width_col < (int)row.size())
        TEST_TRUE(row[iso_width_col].find(';') != std::string::npos);
      if (col_energy_col >= 0 && col_energy_col < (int)row.size())
        TEST_TRUE(row[col_energy_col].find(';') != std::string::npos);
    }
  }
  TEST_TRUE(found_ms2);
  // MS3 must actually be produced and appear in the TSV with two-stage semicolon fields.
  TEST_TRUE(cycle.ms3_cmds.size() > 0);
  TEST_TRUE(found_ms3);

  // Every row should have the same number of columns as the header
  for (const auto& row : tsv.rows)
  {
    TEST_EQUAL(row.size(), tsv.headers.size());
  }
}
END_SECTION

// Test 3: ScanResults TSV -- full MS1->MS2->MS3 cycle with duration tracking
START_SECTION(scan_results_tsv_format)
{
  // Real CytC MS1+MS2 so MS3 result rows + child_ids are actually produced (MS3 enabled
  // in the config below); generic data yields zero MS3 and would leave the MS3 checks dead.
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  const std::string dir = freshLogDir("logging_scan_results");

  // Enable MS3 so we get MS3 result rows
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Full MS1->MS2->MS3 cycle (including feeding MS3 back)
  auto cycle = runFullCycle(&ida, ms1_scans, ms2_scans);

  // Parse and verify
  auto tsv = TSVFile::parse(dir + "/scan_results.tsv");
  TEST_TRUE(tsv.colIndex("tracking_id") >= 0);
  TEST_TRUE(tsv.colIndex("resolve_ts") >= 0);
  TEST_TRUE(tsv.colIndex("duration_ms") >= 0);
  TEST_TRUE(tsv.colIndex("mass_count") >= 0);
  TEST_TRUE(tsv.colIndex("commands_pushed") >= 0);
  TEST_TRUE(tsv.colIndex("child_ids") >= 0);

  // Should have MS1, MS2, and (if MS3 commands were created) MS3 result rows.
  // MS1 results come from engine-emitted survey commands (echoed back via runFullCycle/runInterleaved);
  // MS2 + MS3 from feeding their commands back. runFullCycle's iteration budget guarantees EVERY input MS1
  // is fed, so a short-feed (fewer surveys driven than scans) fails LOUDLY here rather than under-counting.
  ABORT_IF((int)cycle.ms1_cmds.size() != (int)ms1_scans.size())
  int expected_min_rows = (int)ms1_scans.size() + (int)cycle.ms2_cmds.size();
  if (cycle.ms3_cmds.size() > 0)
    expected_min_rows += (int)cycle.ms3_cmds.size();
  TEST_TRUE((int)tsv.rows.size() >= expected_min_rows);

  // Every row should have correct column count
  for (const auto& row : tsv.rows)
  {
    TEST_EQUAL(row.size(), tsv.headers.size());
  }

  // duration_ms should be non-negative
  int dur_col = tsv.colIndex("duration_ms");
  for (const auto& row : tsv.rows)
  {
    if (dur_col >= 0 && dur_col < (int)row.size())
    {
      uint64_t dur = std::stoull(row[dur_col]);
      TEST_TRUE(dur < 3600000ULL);  // < 1 hour: a real upper bound, not a tautology
    }
  }

  // MS3 must be produced, and MS2 result rows must carry child_ids linking to MS3 commands.
  TEST_TRUE(cycle.ms3_cmds.size() > 0);
  {
    int child_col = tsv.colIndex("child_ids");
    bool found_ms2_with_children = false;
    for (const auto& row : tsv.rows)
    {
      if (child_col >= 0 && child_col < (int)row.size() && ! row[child_col].empty())
        found_ms2_with_children = true;
    }
    TEST_TRUE(found_ms2_with_children);
  }
}
END_SECTION

// Test 4: Join integrity -- every child_id in results exists in commands, full MS3 cycle
START_SECTION(join_integrity)
{
  // Real CytC MS1+MS2 with MS3 enabled so the parent-child join graph is actually
  // populated; with generic data + MS3 off every child_ids cell is empty and the join
  // loop below never runs, letting the section pass having asserted nothing.
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  const std::string dir = freshLogDir("logging_join");

  // Enable MS3 for full parent-child graph testing
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Full MS1->MS2->MS3 cycle (with MS3 fed back). runFullCycle already drains the queue to idle,
  // so there is no leftover command to drain here -- the previous trailing
  // `while (getNextScanCommand) { if (is_agc) break; }` loop was redundant (and is exactly the
  // unbounded-drain hazard the harness exists to remove). The TSVs are written incrementally
  // during the cycle, so the parses below see every command/result the cycle produced.
  auto cycle = runFullCycle(&ida, ms1_scans, ms2_scans);

  // Parse both files
  auto cmd_tsv = TSVFile::parse(dir + "/scan_commands.tsv");
  auto res_tsv = TSVFile::parse(dir + "/scan_results.tsv");

  // MS3 must actually have fired, else the join graph is empty and the loop below
  // would validate nothing.
  TEST_TRUE(cycle.ms3_cmds.size() > 0);

  // Build set of all command tracking_ids
  std::set<std::string> cmd_ids;
  int cmd_id_col = cmd_tsv.colIndex("tracking_id");
  for (const auto& row : cmd_tsv.rows)
  {
    if (cmd_id_col >= 0 && cmd_id_col < (int)row.size())
      cmd_ids.insert(row[cmd_id_col]);
  }
  TEST_TRUE(! cmd_ids.empty());

  // Strict join: every child_id in a results row must resolve to a commands-TSV tracking_id,
  // and commands_pushed must equal the child count. child_ids are space-separated (the
  // separator is outside the 0x21-0x7E tracking-id alphabet so ids can never collide with it;
  // ';' would, e.g. the id "!!;"). e.g. MS2 parent !!! -> child_ids "!"B .. !"K" == its 10 MS3 ids.
  int child_col = res_tsv.colIndex("child_ids");
  int pushed_col = res_tsv.colIndex("commands_pushed");
  bool checked_any_child = false;
  for (const auto& row : res_tsv.rows)
  {
    if (child_col >= 0 && child_col < (int)row.size() && ! row[child_col].empty())
    {
      std::istringstream child_ss(row[child_col]);
      std::string child_id;
      int child_count = 0;
      while (std::getline(child_ss, child_id, ' '))
      {
        TEST_TRUE(cmd_ids.count(child_id) > 0);
        child_count++;
      }
      if (pushed_col >= 0 && pushed_col < (int)row.size())
      {
        TEST_EQUAL(std::stoi(row[pushed_col]), child_count);
      }
      checked_any_child = true;
    }
  }
  TEST_TRUE(checked_any_child);

  // ----------------------------------------------------------------------------
  // Backward edge (ADDITIVE -- the forward edge above checks every results child_id
  // resolves to a commands tracking_id; here we walk the lineage the OTHER way).
  // Build tracking_id -> ms_level from the commands TSV so result rows (which classify
  // by their id appearing as an emitted command) can be resolved to a level.
  // ----------------------------------------------------------------------------
  std::map<std::string, int> cmd_level;
  {
    int id_col = cmd_tsv.colIndex("tracking_id");
    int lvl_col = cmd_tsv.colIndex("ms_level");
    TEST_TRUE(id_col >= 0 && lvl_col >= 0);
    for (const auto& row : cmd_tsv.rows)
      if (id_col < (int)row.size() && lvl_col < (int)row.size())
        cmd_level[row[id_col]] = std::atoi(row[lvl_col].c_str());
  }

  int res_id_col = res_tsv.colIndex("tracking_id");
  int res_parent_col = res_tsv.colIndex("parent_tracking_id");
  TEST_TRUE(res_id_col >= 0 && res_parent_col >= 0);

  // (a) No MS1 results-row carries the "~~~" sentinel tracking_id. An MS1 result row is one
  //     whose id is NOT an emitted MS2/MS3 command (i.e. absent from cmd_level OR level==1).
  //     "~~~" is the survey-MS1 sentinel the engine must never echo into a real results row.
  bool checked_ms1_sentinel = false;
  for (const auto& row : res_tsv.rows)
  {
    if (res_id_col >= (int)row.size()) continue;
    const std::string& tid = row[res_id_col];
    auto it = cmd_level.find(tid);
    bool is_ms1 = (it == cmd_level.end()) || (it->second == 1);  // MS1 survey input row
    if (!is_ms1) continue;
    checked_ms1_sentinel = true;
    TEST_TRUE(tid != "~~~");
  }
  TEST_TRUE(checked_ms1_sentinel);

  // (b) Every MS2 results-row parent_tracking_id resolves to an emitted MS1-LEVEL command id.
  //     An MS2 result row is one whose id is an emitted command at level 2; its parent must be
  //     present in the commands level map at level 1 (the survey MS1 that spawned it).
  bool checked_ms2_parent = false;
  bool ms2_parent_ok = true;
  for (const auto& row : res_tsv.rows)
  {
    if (res_id_col >= (int)row.size() || res_parent_col >= (int)row.size()) continue;
    auto it = cmd_level.find(row[res_id_col]);
    if (it == cmd_level.end() || it->second != 2) continue;  // MS2 result rows only
    const std::string& parent = row[res_parent_col];
    checked_ms2_parent = true;
    auto pit = cmd_level.find(parent);
    ms2_parent_ok = ms2_parent_ok && (pit != cmd_level.end()) && (pit->second == 1);
  }
  TEST_TRUE(checked_ms2_parent);
  TEST_TRUE(ms2_parent_ok);
}
END_SECTION

// Test 5: Crash safety -- the command/results TSV files stay valid across a full MS3 cycle
START_SECTION(crash_safety_valid_tsv)
{
  // Real CytC MS1+MS2 with MS3 enabled so the MS2 and MS3 crash-safety paths are actually
  // reached; with MS3 off the MS3 commands never fire and the section would validate only the
  // MS1 path. (Driven via runInterleaved -- see the per-stage rationale at the cycle call below.)
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  const std::string dir = freshLogDir("logging_crash");

  // Enable MS3 for full cycle
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // After constructor: headers should exist (pre-cycle: no scans driven yet)
  {
    auto cmd_tsv = TSVFile::parse(dir + "/scan_commands.tsv");
    TEST_TRUE(cmd_tsv.headers.size() > 0);
    auto res_tsv = TSVFile::parse(dir + "/scan_results.tsv");
    TEST_TRUE(res_tsv.headers.size() > 0);
  }

  // Drive the full MS1->MS2->MS3 cycle via the canonical interleaved driver. The old staged feed
  // (pushAllScans + manual MS2 drain + manual MS2-feed + manual MS3 drain + manual MS3-feed) fed MS1
  // under fabricated encode(800000+i) ids, which the always-on MS1 gate now rejects -> 0 precursors ->
  // the MS2/MS3 stages were never reached. runInterleaved pulls one command at a time and feeds the
  // matching cytC scan back stamped with the ENGINE's own descriptor (MS3 via the MS2-as-MS3 shortcut,
  // null manifest), so MS1->MS2->MS3 chains off the engine's ids end-to-end. The files are written
  // incrementally during the call, so the post-cycle validity parse below confirms no partial/corrupt
  // row was ever emitted -- the same crash-safety invariant the staged checks asserted per operation.
  auto cycle = runFullCycle(&ida, ms1_scans, ms2_scans);
  TEST_TRUE(cycle.ms2_cmds.size() > 0);  // MS2 commands produced (was: ms2_cmds.size() > 0 after drain)
  TEST_TRUE(cycle.ms3_cmds.size() > 0);  // MS3 must have fired (was: ms3_cmds.size() > 0 after drain)

  // results file is valid: >=1 row and every row has the header column count (no torn writes)
  {
    auto res_tsv = TSVFile::parse(dir + "/scan_results.tsv");
    TEST_TRUE(res_tsv.rows.size() >= 1);
    for (const auto& row : res_tsv.rows)
      TEST_EQUAL(row.size(), res_tsv.headers.size());
  }

  // commands file is valid: >=1 row and every row has the header column count (no torn writes)
  {
    auto cmd_tsv = TSVFile::parse(dir + "/scan_commands.tsv");
    TEST_TRUE(cmd_tsv.rows.size() >= 1);
    for (const auto& row : cmd_tsv.rows)
      TEST_EQUAL(row.size(), cmd_tsv.headers.size());
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// Lock-split guards.
//
// Neither of these is a golden. They assert RELATIONSHIPS between logged values and structural
// invariants of the writer -- never an absolute number -- so they cannot drift into being a second,
// unmanaged copy of the log goldens.
/////////////////////////////////////////////////////////////

// Every dequeued MSn command carries a real precursor_id, and children inherit their parent's.
START_SECTION(every_dequeued_command_logs_a_row_with_the_right_precursor_id)
{
  // WHAT THIS PINS, and why it is worth its lines even though it passes today.
  //
  // The precursor_id written on a scan_commands row is read at dequeue time from a map that
  // processScan populates. The lock split moves that read behind its own mutex and routes the six
  // writes through a helper. Two mis-applications of that change are silent:
  //
  //   1. Deleting the whole braced block at FLASHIda.cpp:707-716 instead of just the lock_guard
  //      LINE inside it. That removes writeScanCommandRow entirely and every non-AGC row with it.
  //      Nothing else in this suite would notice -- the log goldens would, but they are deliberately
  //      out of scope here.
  //   2. Converting a write site to the helper but dropping or mis-keying its argument, so a whole
  //      family of children silently logs precursor_id 0.
  //
  // Both are caught by asserting the ROW EXISTS and its value RELATES correctly to its siblings.
  // The cytC + MS3 recipe, same as scan_commands_tsv_format above. MS3 is REQUIRED here, not
  // incidental: inheritance is only observable when a child's parent itself carries a non-zero
  // precursor_id, and an MS2's parent is an MS1, which logs 0 by definition. With MS3 off the
  // inheritance loop below would iterate, skip every pair, and assert nothing -- vacuous and green.
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  const std::string dir = freshLogDir("logging_precursor_id_pin");
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  auto cycle = runFullCycle(&ida, ms1_scans, ms2_scans);
  TEST_TRUE(cycle.ms2_cmds.size() > 0);
  TEST_TRUE(cycle.ms3_cmds.size() > 0);  // no MS3 -> the inheritance check below is vacuous

  auto tsv = TSVFile::parse(dir + "/scan_commands.tsv");
  TEST_TRUE(tsv.colIndex("precursor_id") >= 0);

  // FAIL-CLOSED. Zero MSn rows is a failure, not a skip -- that is precisely mis-application (1),
  // and a test that quietly passes on an empty file is the thing this whole exercise exists to avoid.
  int msn_rows = 0;
  int msn_rows_with_precursor = 0;
  std::map<std::string, std::string> precursor_by_tracking;  // tracking_id -> precursor_id
  std::vector<std::pair<std::string, std::string>> child_parent;  // (tracking_id, parent_tracking_id)

  for (const auto& row : tsv.rows)
  {
    const std::string lvl = cell(tsv, row, "ms_level");
    const std::string tid = cell(tsv, row, "tracking_id");
    const std::string pid = cell(tsv, row, "precursor_id");
    const std::string par = cell(tsv, row, "parent_tracking_id");

    if (!tid.empty()) precursor_by_tracking[tid] = pid;

    if (lvl == "2" || lvl == "3")
    {
      msn_rows++;
      if (!pid.empty() && pid != "0") msn_rows_with_precursor++;
      if (!par.empty() && par != "0") child_parent.emplace_back(tid, par);
    }
  }

  TEST_TRUE(msn_rows > 0);                              // mis-application (1) fails here
  TEST_EQUAL(msn_rows_with_precursor, msn_rows);        // mis-application (2) fails here

  // A child command inherits its parent's precursor_id verbatim -- that inheritance is the whole
  // reason five of the six write sites exist, and it is invisible in any single row.
  int checked_inheritance = 0;
  for (const auto& cp : child_parent)
  {
    auto par_it = precursor_by_tracking.find(cp.second);
    if (par_it == precursor_by_tracking.end()) continue;  // parent not in this file (AGC/idle parents)
    if (par_it->second.empty() || par_it->second == "0") continue;
    TEST_EQUAL(precursor_by_tracking[cp.first], par_it->second);
    checked_inheritance++;
  }
  std::cout << "[PID-PIN] msn_rows=" << msn_rows << " inheritance_pairs_checked=" << checked_inheritance << std::endl;

  // FAIL-CLOSED on the check itself. Without this, a config or fixture change that stopped producing
  // parent-carrying children would turn the loop above into a no-op and this test would keep passing
  // while asserting nothing -- the exact failure mode it was written to catch elsewhere.
  TEST_TRUE(checked_inheritance > 0);
}
END_SECTION

// Concurrent drains must each write exactly one intact row.
START_SECTION(concurrent_drain_writes_one_wellformed_row_per_call)
{
  // WHAT THIS PINS: scan_commands.tsv stays intact when getNextScanCommand runs on several threads.
  //
  // OBSERVED, not hypothesised. This test was landed green, the coarse lock was then removed while
  // the per-stream logger mutexes were deliberately held back one commit, and CI reported:
  //
  //     [DRAIN-CONCURRENCY] calls=1000 rows=1000 wellformed=17 unique_ids=156
  //
  // 17 intact rows out of 1000. That run is why IdaLogger owns a mutex per stream rather than
  // relying on the streams being thread-disjoint in practice.
  //
  // No processScan, no spectra, no id echo. The MS1 admission gate makes fabricated ids useless, but
  // the DRAIN needs none of that -- it manufactures its own work. That is what makes this cheap and
  // immune to fixture drift.
  //
  // THREE ASSERTIONS, and the one that actually caught it was NOT the one expected to. Every path
  // through getNextScanCommand writes exactly one row and returns 1 (Step 1, Step 4 and Step 5 each
  // write once; Step 2 pushes and falls through rather than returning; Step 3 writes nothing), so
  // rows == calls. The original reasoning was that a torn row MERGES two lines and drops that count.
  // It does not: the count held at exactly 1000 and the damage was entirely intra-line -- interleaved
  // field writes within correctly-terminated rows. The row count alone would have passed. Keep all
  // three; the field-count and unique-id checks are what have actually earned their keep.
  const std::string dir = freshLogDir("logging_concurrent_drain");
  std::string json = buildJsonWithLogDir(dir);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  const int kThreads = 4;
  const int kCallsPerThread = 250;
  const int kTotalCalls = kThreads * kCallsPerThread;

  std::vector<std::thread> workers;
  for (int t = 0; t < kThreads; t++)
  {
    workers.emplace_back([&ida, kCallsPerThread]() {
      for (int i = 0; i < kCallsPerThread; i++)
      {
        ScanCommand cmd {};
        ida.getNextScanCommand(cmd);
      }
    });
  }
  for (auto& w : workers) w.join();

  auto tsv = TSVFile::parse(dir + "/scan_commands.tsv");

  // Exactly one row per call. Held even under full tearing, so on its own this proves little -- but
  // it is the assertion that would catch a path that stopped writing, or started writing twice.
  TEST_EQUAL((int)tsv.rows.size(), kTotalCalls);

  // Every row has the full column count. THIS is the one that catches interleaved writes: fields
  // from two threads land in one line and the row tokenizes to the wrong width.
  int wellformed = 0;
  for (const auto& row : tsv.rows)
    if (row.size() == tsv.headers.size()) wellformed++;
  TEST_EQUAL(wellformed, (int)tsv.rows.size());

  // Tracking ids are allocated under queue_mutex_, so they must all be distinct even across threads.
  // Independent of the width check: a spliced row can tokenize to the right width and still carry a
  // duplicated or garbled id, which is what the observed run showed (156 distinct ids across 1000
  // rows).
  std::set<std::string> ids;
  for (const auto& row : tsv.rows) ids.insert(cell(tsv, row, "tracking_id"));
  TEST_EQUAL((int)ids.size(), (int)tsv.rows.size());

  std::cout << "[DRAIN-CONCURRENCY] calls=" << kTotalCalls << " rows=" << tsv.rows.size()
            << " wellformed=" << wellformed << " unique_ids=" << ids.size() << std::endl;
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0046 -- ScanCommandJoin: reading a scan_commands.tsv back, and joining it to a data file
//
// FLASHDeconv locates the precursor of a commanded MS2 through this join instead of searching for
// it. The join is pure (file -> rows, data -> data), which is why it is tested HERE: no CI job runs
// FLASHDeconv. Each section states the bug under which it fails.
/////////////////////////////////////////////////////////////

// T1 -- every file the engine has ever written must parse (ADR-0046 decision 9).
// Fails if a required column is misnamed, or the reader demands a column old files do not have.
START_SECTION(scan_commands_parse_reads_a_pre_adr0046_run_folder_file)
{
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t01", {kCmdHeader, kS32, kS33, kM37}));
  TEST_EQUAL(rows.size(), 3)
  ABORT_IF(rows.count("!!B") == 0 || rows.count("!!=") == 0)

  const auto& m = rows["!!B"];
  TEST_EQUAL(m.ms_level, 2)
  TEST_EQUAL(m.parent_tracking_id, std::string("!!="))
  TEST_EQUAL(m.charge, 10)
  TEST_TRUE(std::abs(m.mono_mass - 5305.33) < 1e-9)
  TEST_TRUE(std::abs(m.precursor_mz - 531.791) < 1e-9)
  TEST_TRUE(std::abs(m.qscore - 0.443474) < 1e-9)
  TEST_TRUE(std::abs(m.charge_snr - 5.70422) < 1e-9)
  TEST_TRUE(std::abs(m.ppm_error - 0.608297) < 1e-9)
  TEST_TRUE(std::abs(m.precursor_intensity - 2134.16) < 1e-9)
  TEST_TRUE(std::abs(m.peakgroup_intensity - 2947.43) < 1e-9)

  // A non-MS2 row keeps id / parent / level and nothing else.
  const auto& s = rows["!!="];
  TEST_EQUAL(s.ms_level, 1)
  TEST_TRUE(s.parent_tracking_id.empty())
  TEST_EQUAL(s.mono_mass, 0.0)
}
END_SECTION

// T2 -- the drift guard: whatever the WRITER emits, the reader must take.
// Fails if IdaLogger renames or drops a column the reader needs -- so the break lands in a test, not
// in somebody's analysis.
START_SECTION(scan_commands_parse_roundtrips_an_engine_written_file)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  const std::string dir = freshLogDir("scj_t02");
  std::string json = buildJsonWithLogDir(dir);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  AcqResult acq = runInterleaved(&ida, ms1_scans, std::vector<ScanData>{});
  ABORT_IF(acq.ms2_cmds.empty())

  auto t = TSVFile::parse(dir + "/scan_commands.tsv");
  auto rows = ScanCommandJoin::parse(dir + "/scan_commands.tsv");
  TEST_EQUAL(rows.size(), t.rows.size())   // every row the engine wrote is a Row

  int ms2_rows = 0;
  bool fields_ok = true, parent_ok = true;
  for (const auto& row : t.rows)
  {
    if (cell(t, row, "ms_level") != "2") continue;
    ms2_rows++;
    auto it = rows.find(cell(t, row, "tracking_id"));
    if (it == rows.end()) { fields_ok = false; continue; }
    const auto& r = it->second;
    fields_ok = fields_ok && std::abs(r.mono_mass - toD(cell(t, row, "mono_mass"))) < 1e-6
                          && r.charge == (int)toD(cell(t, row, "charge"))
                          && std::abs(r.precursor_mz - toD(cell(t, row, "precursor_mz"))) < 1e-6
                          && std::abs(r.qscore - toD(cell(t, row, "qscore"))) < 1e-9
                          && r.mono_mass > 0 && r.charge > 0 && r.precursor_mz > 0;
    // a production MS2's parent is its survey
    auto p = rows.find(r.parent_tracking_id);
    parent_ok = parent_ok && p != rows.end() && p->second.ms_level == 1;
  }
  TEST_TRUE(ms2_rows > 0)   // vacuity guard
  TEST_TRUE(fields_ok)
  TEST_TRUE(parent_ok)
}
END_SECTION

// T3 -- a multiplexed row carries "anchor,notch,notch"; the join wants the anchor (ADR-0016 / 0046 d7).
// Fails if the whole cell is parsed ("10,9" is not a number) or the LAST notch is taken.
START_SECTION(scan_commands_parse_takes_the_anchor_of_a_multiplexed_row)
{
  Cells msx = withCell(withCell(withCell(kM37, "charge", "10,9"), "precursor_mz", "531.791,590.878"), "isolation_width", "0.9,1.0");
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t03", {kCmdHeader, kS32, msx}));
  ABORT_IF(rows.count("!!B") == 0)
  TEST_EQUAL(rows["!!B"].charge, 10)
  TEST_TRUE(std::abs(rows["!!B"].precursor_mz - 531.791) < 1e-9)
}
END_SECTION

// T4 -- columns are resolved by NAME, so a reorder (and a column from the future) is free.
// Fails if any index is positional.
START_SECTION(scan_commands_parse_resolves_columns_by_name)
{
  // mono_mass moved to the front, plus an unknown trailing column.
  auto permute = [](const Cells& c, const std::string& extra) {
    Cells out;
    out.push_back(c[6]);
    for (size_t i = 0; i < c.size(); ++i) if (i != 6) out.push_back(c[i]);
    out.push_back(extra);
    return out;
  };
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t04", {permute(kCmdHeader, "future_col"), permute(kS32, "x"), permute(kM37, "y")}));
  ABORT_IF(rows.count("!!B") == 0)
  const auto& m = rows["!!B"];
  TEST_EQUAL(m.ms_level, 2)
  TEST_EQUAL(m.parent_tracking_id, std::string("!!="))
  TEST_EQUAL(m.charge, 10)
  TEST_TRUE(std::abs(m.mono_mass - 5305.33) < 1e-9)
  TEST_TRUE(std::abs(m.precursor_mz - 531.791) < 1e-9)
  TEST_TRUE(std::abs(m.peakgroup_intensity - 2947.43) < 1e-9)
}
END_SECTION

// T5-T9 -- a bad file FAILS CLOSED (ADR-0046 decision 5). Each is one way a bad file could otherwise
// slip through as "a few rows fewer".
START_SECTION(scan_commands_parse_rejects_a_missing_file)
{
  const std::string dir = freshLogDir("scj_t05");
  TEST_EXCEPTION(Exception::FileNotFound, ScanCommandJoin::parse(dir + "/absent.tsv"))
}
END_SECTION

START_SECTION(scan_commands_parse_rejects_a_missing_required_column)
{
  auto drop = [](const Cells& c) { Cells out = c; out.erase(out.begin() + 3); return out; };   // parent_tracking_id
  const std::string path = writeCommandsFile("scj_t06", {drop(kCmdHeader), drop(kS32), drop(kM37)});
  TEST_EXCEPTION(Exception::ParseError, ScanCommandJoin::parse(path))
}
END_SECTION

START_SECTION(scan_commands_parse_rejects_a_short_row)
{
  Cells cut(kM37.begin(), kM37.begin() + 20);
  const std::string path = writeCommandsFile("scj_t07", {kCmdHeader, kS32, cut});
  TEST_EXCEPTION(Exception::ParseError, ScanCommandJoin::parse(path))
}
END_SECTION

START_SECTION(scan_commands_parse_rejects_a_non_numeric_mass)
{
  const std::string path = writeCommandsFile("scj_t08", {kCmdHeader, kS32, withCell(kM37, "mono_mass", "n/a")});
  TEST_EXCEPTION(Exception::ParseError, ScanCommandJoin::parse(path))
}
END_SECTION

START_SECTION(scan_commands_parse_rejects_a_duplicate_tracking_id)
{
  const std::string path = writeCommandsFile("scj_t09", {kCmdHeader, kS32, kM37, kM37});
  TEST_EXCEPTION(Exception::InvalidValue, ScanCommandJoin::parse(path))
}
END_SECTION

// T10 -- THE assertion of ADR-0046: the survey is the one the COMMAND names.
// Fails if the join takes "the MS1 before me" -- which is what FLASHDeconv did, and is wrong for 88 %
// of commanded MS2 (scan 33 is a newer survey sitting between the command's survey 32 and its scan 37).
START_SECTION(scan_commands_join_names_the_survey_the_command_came_from)
{
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t10", {kCmdHeader, kS32, kS33, kM37}));
  auto located = ScanCommandJoin::join(realScans(), rows);
  TEST_EQUAL(located.size(), 1)
  ABORT_IF(located.count(37) == 0)
  TEST_EQUAL(located[37].parent_scan_number, 32)
  TEST_EQUAL(located[37].row.tracking_id, std::string("!!B"))
  TEST_TRUE(std::abs(located[37].row.mono_mass - 5305.33) < 1e-9)
}
END_SECTION

// T11 -- a follow-up MS2's parent is the MS2 that TRIGGERED it (ScanCommandQueue::buildFollowUp), not
// the survey; the join walks the parent chain up to the MS1 row.
// Fails if the parent is assumed to be the survey.
START_SECTION(scan_commands_join_walks_a_follow_up_up_to_its_survey)
{
  Cells c40 = withCell(withCell(withCell(withCell(kM37, "tracking_id", "!!C"), "scan_type", "conditional"),
                                "parent_tracking_id", "!!B"), "scan_description", "!!CC5.30533k@10");
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t11", {kCmdHeader, kS32, kS33, kM37, c40}));
  auto scans = realScans();
  scans.push_back(scanOf(40, 2, "!!CC5.30533k@10", {531.790649}));
  auto located = ScanCommandJoin::join(scans, rows);
  TEST_EQUAL(located.size(), 2)
  ABORT_IF(located.count(40) == 0)
  TEST_EQUAL(located[40].parent_scan_number, 32)   // through !!B, up to !!=
}
END_SECTION

// T12 -- an uncommanded scan is not an error (ADR-0046 decision 4).
// Fails if a blank or unknown id aborts the join, or if the three-blank non-id the instrument's own
// scans have been seen to carry trips the duplicate check when it appears twice.
START_SECTION(scan_commands_join_passes_uncommanded_scans_through)
{
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t12", {kCmdHeader, kS32, kS33, kM37}));
  auto scans = realScans();
  scans.push_back(scanOf(34, 1, ""));
  scans.push_back(scanOf(35, 1, "   "));
  scans.push_back(scanOf(36, 1, "   "));
  scans.push_back(scanOf(38, 2, "~~~R9.99k@9", {700.0}));   // another acquisition's id
  auto located = ScanCommandJoin::join(scans, rows);
  TEST_EQUAL(located.size(), 1)
  TEST_EQUAL(located.count(37), 1)
  TEST_EQUAL(located.count(38), 0)
}
END_SECTION

// T13 -- the survey can be legitimately absent from the data file (an RT crop), and a root MS2 names
// no parent at all. Neither is an error: parent_scan_number is -1 and FLASHDeconv keeps its search.
// Fails if a cropped input aborts, or if -1 is ever taken for a scan number.
START_SECTION(scan_commands_join_tolerates_a_survey_missing_from_the_data_file)
{
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t13a", {kCmdHeader, kS32, kS33, kM37}));
  std::vector<ScanCommandJoin::Scan> cropped = {scanOf(33, 1, "!!>S"), scanOf(37, 2, "!!BR5.30533k@10", {531.790649})};
  auto located = ScanCommandJoin::join(cropped, rows);
  ABORT_IF(located.count(37) == 0)
  TEST_EQUAL(located[37].parent_scan_number, -1)

  auto root_rows = ScanCommandJoin::parse(writeCommandsFile("scj_t13b", {kCmdHeader, kS32, withCell(kM37, "parent_tracking_id", "")}));
  auto root_located = ScanCommandJoin::join(realScans(), root_rows);
  ABORT_IF(root_located.count(37) == 0)
  TEST_EQUAL(root_located[37].parent_scan_number, -1)
}
END_SECTION

// T14 / T15 -- tracking ids restart in every run, so a FOREIGN scan_commands.tsv joins every scan by
// id (measured: replicate R2's file joins 13,873 of 13,873 of R1's MS2). The row-vs-spectrum check is
// what tells them apart. Fail if a foreign file joins silently.
START_SECTION(scan_commands_join_rejects_a_foreign_file_by_mz)
{
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t14", {kCmdHeader, kS32, kS33, withCell(kM37, "precursor_mz", "600.123")}));
  TEST_EXCEPTION(Exception::InvalidValue, ScanCommandJoin::join(realScans(), rows))
}
END_SECTION

START_SECTION(scan_commands_join_rejects_a_foreign_file_by_ms_level)
{
  // In the other run, id !!= was an MS2; here the spectrum carrying it is a survey.
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t15", {kCmdHeader, withCell(kS32, "ms_level", "2"), kS33, kM37}));
  TEST_EXCEPTION(Exception::InvalidValue, ScanCommandJoin::join(realScans(), rows))
}
END_SECTION

// T16 -- a file that joins NOTHING is the extreme foreign file, and also what a converter that drops
// the scan description produces. -FD:scan_commands was given, so an uncoupled run is never silent.
START_SECTION(scan_commands_join_rejects_a_file_that_joins_nothing)
{
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t16", {kCmdHeader, kS32, kS33, kM37}));
  std::vector<ScanCommandJoin::Scan> bare = {scanOf(32, 1, ""), scanOf(33, 1, ""), scanOf(37, 2, "", {531.790649})};
  TEST_EXCEPTION(Exception::InvalidValue, ScanCommandJoin::join(bare, rows))
}
END_SECTION

// T17 -- one commanded id on two spectra makes every parent lookup ambiguous.
START_SECTION(scan_commands_join_rejects_one_id_on_two_spectra)
{
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t17", {kCmdHeader, kS32, kS33, kM37}));
  auto scans = realScans();
  scans.push_back(scanOf(41, 2, "!!BR5.30533k@10", {531.790649}));
  TEST_EXCEPTION(Exception::InvalidValue, ScanCommandJoin::join(scans, rows))
}
END_SECTION

// T18 -- an MSX scan has one precursor element per notch; the anchor may be any of them.
// Fails if only the first isolation target is compared.
START_SECTION(scan_commands_join_accepts_any_isolation_target_of_an_msx_scan)
{
  auto rows = ScanCommandJoin::parse(writeCommandsFile("scj_t18", {kCmdHeader, kS32, kS33, kM37}));
  std::vector<ScanCommandJoin::Scan> scans = {scanOf(32, 1, "!!=S"), scanOf(33, 1, "!!>S"),
                                              scanOf(37, 2, "!!BR5.30533k@10", {590.878, 531.790649})};
  auto located = ScanCommandJoin::join(scans, rows);
  TEST_EQUAL(located.count(37), 1)
}
END_SECTION

// T19 -- the locator tolerance and the isotope ladder, in one place (ADR-0046 decision 3).
START_SECTION(scan_commands_massRank)
{
  double res = -1;
  // A six-significant-digit logged mass (every pre-ADR-0046 file) still locates the species.
  TEST_EQUAL(ScanCommandJoin::massRank(12351.3933, 12351.4, res), 0)
  TEST_TRUE(res >= 0 && res < 0.01)   // the tie-break input is the distance to the ACCEPTED isotope

  TEST_EQUAL(ScanCommandJoin::massRank(12352.3957, 12351.3933, res), 1)    // +1 isotope
  TEST_TRUE(res < 0.001)
  TEST_EQUAL(ScanCommandJoin::massRank(12349.3886, 12351.3933, res), 2)    // -2 isotopes
  TEST_TRUE(res < 0.001)
  TEST_EQUAL(ScanCommandJoin::massRank(12354.4004, 12351.3933, res), -1)   // three isotopes: another call
  TEST_EQUAL(ScanCommandJoin::massRank(12351.8933, 12351.3933, res), -1)   // half a dalton: another mass
}
END_SECTION

// T20 -- scan_commands.tsv's mono_mass is the DECISION value, at four decimals (ADR-0046 decision 10).
//
// It used to go through sc(), i.e. the stream default of six SIGNIFICANT digits: "12351.4" for a
// commanded 12351.3933 -- the defect ADR-0035 decision 5 fixed in ida.log. This section is the ONLY
// gate on the fix: the C# golden comparer accepts 12351.4 against 12351.3933 with ~1800x headroom
// (RelTol 1e-3), so no golden can see the precision either way.
// Fails if mono_mass goes back through sc(), or if the stage-less "0" becomes "0.0000" -- which would
// revalue every MS1 / AGC row and break FLASHIda_LoggingFields_test::commands_ms1_agc_stageless.
START_SECTION(scan_commands_mono_mass_is_written_at_four_decimals)
{
  auto ms1 = loadTsvScans(FI_MS1_CYTC);
  auto ms2 = loadTsvScans(FI_MS2_CYTC);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("scj_t20");
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  const int budget = 256 + 64 * static_cast<int>(ms1.size() + ms2.size());
  AcqResult acq = runInterleaved(&ida, ms1, ms2, nullptr, budget);
  ABORT_IF(acq.ms2_cmds.empty() || acq.ms3_cmds.empty())   // vacuity guard: all three row shapes must occur

  auto t = TSVFile::parse(dir + "/scan_commands.tsv");
  std::map<std::string, std::string> text_of;   // tracking_id -> the mono_mass cell AS WRITTEN
  for (const auto& row : t.rows) text_of[cell(t, row, "tracking_id")] = cell(t, row, "mono_mass");

  auto fixed4 = [](double v) { std::ostringstream os; os << std::fixed << std::setprecision(4) << v; return os.str(); };

  int stageless = 0, ms2_rows = 0, ms3_rows = 0, above_10k = 0, missing = 0, wrong = 0;
  for (const auto& c : acq.all_cmds)
  {
    auto it = text_of.find(ScanCommandQueue::encode(c.scan_id));
    if (it == text_of.end()) { missing++; continue; }
    std::string want;
    if (c.num_stages == 0)      { want = "0"; stageless++; }
    else if (c.msn_level == 3)  { want = fixed4(c.mono_mass) + ";" + fixed4(c.mono_mass_s1); ms3_rows++; }
    else                        { want = fixed4(c.mono_mass); ms2_rows++; if (c.mono_mass >= 10000.0) above_10k++; }
    if (it->second != want)     // byte for byte
    {
      wrong++;
      std::cout << "[MONO-MASS] id=" << it->first << " wrote '" << it->second << "' want '" << want << "'" << std::endl;
    }
  }
  TEST_EQUAL(missing, 0)
  TEST_EQUAL(wrong, 0)
  TEST_TRUE(stageless > 0)
  TEST_TRUE(ms2_rows > 0)
  TEST_TRUE(ms3_rows > 0)
  TEST_TRUE(above_10k > 0)   // where six significant digits lose a decimal
}
END_SECTION

// ---- Tracking ids that begin with '"' (LogCell.h) ----------------------------------------------------
//
// Tracking ids 8,836 to 17,671 begin with '"' -- 5,052 of the 13,887 scan_commands rows of an Eclip_1480
// MS1 replay -- and a CSV-style reader takes it as an opening quote: Python's csv read that file as 11,400
// rows, pandas stopped with a tokenizing error. IdaLogger now quotes exactly the cells that begin with '"'
// (RFC 4180: enclose, double the quotes); every other cell, and so every log of a shorter run, is unchanged.

// The rule itself, and its inverse, over every id that begins with '"' and every id that does not.
// Fails if a cell without a leading '"' is touched (that would move every golden: id 1, "!!\"", opens every
// run), if a quoted id does not read back, or if a raw id of an old file is taken for a quoted one.
START_SECTION(log_cells_are_quoted_exactly_when_they_begin_with_a_quote)
{
  TEST_EQUAL(quoteLogCell(""), std::string())
  TEST_EQUAL(quoteLogCell("!!\""), std::string("!!\""))                 // id 1: a '"', not a leading one
  TEST_EQUAL(quoteLogCell("!!\"R5.30533k@10"), std::string("!!\"R5.30533k@10"))
  TEST_EQUAL(quoteLogCell("\"!!"), std::string("\"\"\"!!\""))          // id 8,836
  TEST_EQUAL(quoteLogCell("\"!\""), std::string("\"\"\"!\"\"\""))      // id 8,837: every quote doubled
  TEST_EQUAL(quoteLogCell("\"!# !!$"), std::string("\"\"\"!# !!$\""))  // an id list led by such an id

  int touched = 0, lost = 0;
  for (int v = 1; v < 3 * 8836; ++v)
  {
    const std::string id = ScanCommandQueue::encode(v);
    const std::string desc = id + "R5.30533k@10";
    if (id.front() != '"')
    {
      if (quoteLogCell(id) != id || quoteLogCell(desc) != desc || unquoteLogId(id) != id) { ++touched; }
      continue;
    }
    const std::string q = quoteLogCell(id);
    if (q.size() < 6 || unquoteLogId(q) != id || unquoteLogCell(q) != id || unquoteLogCell(quoteLogCell(desc)) != desc) { ++lost; }
  }
  TEST_EQUAL(touched, 0)
  TEST_EQUAL(lost, 0)

  // A file re-saved by pandas or Excel quotes every cell that CONTAINS a '"'.
  TEST_EQUAL(unquoteLogCell("\"!!\"\"\""), std::string("!!\""))
  TEST_EQUAL(unquoteLogId("\"!!\"\"\""), std::string("!!\""))
  // Not quoted cells: a lone inner quote (the raw id 8,931), an unclosed one.
  TEST_EQUAL(unquoteLogCell("\"\"\""), std::string("\"\"\""))
  TEST_EQUAL(unquoteLogCell("\"!#"), std::string("\"!#"))
  // A raw id of a file written before the rule: unquoteLogCell would take "!" for a quoted !, unquoteLogId never.
  TEST_EQUAL(unquoteLogCell("\"!\""), std::string("!"))
  TEST_EQUAL(unquoteLogId("\"!\""), std::string("\"!\""))
}
END_SECTION

// All four TSV streams, written by IdaLogger with ids that begin with '"', read three ways: by a CSV-style
// reader (the defect), by the class tests' TSVFile, and by ScanCommandJoin, which is what FLASHDeconv's
// -FD:scan_commands runs. No fixture and no acquisition: the rows are written directly.
// Fails if a stream writes such an id raw (the CSV reader merges rows), or if a reader in this repository
// does not unquote it (the ids or the join no longer match the engine's).
START_SECTION(ids_that_begin_with_a_quote_survive_every_stream_and_reader)
{
  const std::string s_id = ScanCommandQueue::encode(8836);   // "!!  survey
  const std::string m_id = ScanCommandQueue::encode(8837);   // "!"  MS2
  const std::string c_id = ScanCommandQueue::encode(8838);   // "!#  a child
  const std::string w_id = ScanCommandQueue::encode(8839);   // "!$  a winner
  const std::string one = ScanCommandQueue::encode(1);       // !!"  a '"' that does not lead
  ABORT_IF(s_id != "\"!!" || m_id != "\"!\"" || one != "!!\"")

  const std::string dir = freshLogDir("logcell_t02");
  Config cfg{buildJsonWithLogDir(dir)};
  {
    IdaLogger logger(cfg);

    ScanCommand survey{};
    survey.scan_id = 8836;
    survey.msn_level = 1;
    ScanCommandQueue::writeMs1Description(survey, ScanRole::Survey);
    logger.writeScanCommandRow(survey);

    ScanCommand ms2{};
    ms2.scan_id = 8837;
    ms2.msn_level = 2;
    ms2.priority = 2;
    ms2.num_stages = 1;
    ms2.mono_mass = 5305.33;
    ms2.qscore = 0.443474;
    ms2.stages[0].charge_state = 10;
    ms2.stages[0].precursor_mz = 531.791;
    ms2.stages[0].isolation_width = 0.9;
    ms2.stages[0].collision_energy = 29;
    std::strncpy(ms2.stages[0].activation_type, "HCD", sizeof(ms2.stages[0].activation_type) - 1);
    std::strncpy(ms2.parent_scan_id, s_id.c_str(), sizeof(ms2.parent_scan_id) - 1);
    std::snprintf(ms2.scan_description, sizeof(ms2.scan_description), "%sR5.30533k@10", m_id.c_str());
    logger.writeScanCommandRow(ms2, 1);

    ScanCommand early = ms2;   // the same command as id 1: written exactly as before
    early.scan_id = 1;
    std::snprintf(early.scan_description, sizeof(early.scan_description), "%sR5.30533k@10", one.c_str());
    logger.writeScanCommandRow(early, 1);

    IdaLogger::ScanRowDescriptor res;
    res.tracking_id = m_id;
    res.ms_level = 2;
    res.parent_tracking_id = s_id;
    res.child_ids = {c_id, one};
    res.winner_tracking_id = w_id;
    logger.writeScanResultRow(res);

    IdaLogger::IdRowDescriptor idr{};
    idr.tracking_id = m_id;
    idr.ms_level = 2;
    idr.scan_mode = 'R';
    FragmentAnalysis::ProteoformMatch::FragmentMatch frag;
    frag.ion_type = "b";
    frag.ion_index = 5;
    frag.observed_mass = 500.25;
    idr.match.fragments.push_back(frag);
    logger.writeIdentificationRow(idr);

    IdaLogger::PooledModelDescriptor pooled;
    pooled.trigger_scan_id = m_id;
    pooled.trigger = "MS2";
    pooled.contributing_scan_ids = {8837, 1};
    logger.writePooledModelRow(pooled);
  }

  // (1) The bytes: a leading '"' is quoted, a '"' elsewhere is not.
  auto lines = [](const std::string& path) {
    std::vector<std::string> out;
    std::ifstream in(path);
    std::string l;
    while (std::getline(in, l)) out.push_back(l);
    return out;
  };
  const auto cmd_lines = lines(dir + "/scan_commands.tsv");
  ABORT_IF(cmd_lines.size() != 4)
  TEST_EQUAL(cmd_lines[1].substr(0, cmd_lines[1].find('\t')), std::string("\"\"\"!!\""))
  TEST_EQUAL(cmd_lines[2].substr(0, cmd_lines[2].find('\t')), std::string("\"\"\"!\"\"\""))
  TEST_EQUAL(cmd_lines[3].substr(0, cmd_lines[3].find('\t')), std::string("!!\""))
  TEST_TRUE(cmd_lines[3].find("\t!!\"R5.30533k@10\t") != std::string::npos)   // id 1's description, raw

  // (2) A CSV-style reader -- Python's csv / pandas at their defaults: a '"' opens a quoted cell only at
  // the start of a cell, a doubled '"' inside one is a quote, and a quoted cell runs across tabs and line
  // ends. Every file must come back with one row per line, each as wide as the header.
  auto csvRows = [](const std::string& path) {
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    const std::string text = ss.str();
    std::vector<std::vector<std::string>> rows(1);
    std::string c;
    bool quoted = false, at_start = true;
    for (size_t i = 0; i < text.size(); ++i)
    {
      const char ch = text[i];
      if (quoted)
      {
        if (ch != '"') { c += ch; }
        else if (i + 1 < text.size() && text[i + 1] == '"') { c += '"'; ++i; }
        else { quoted = false; }
        continue;
      }
      if (ch == '"' && at_start) { quoted = true; at_start = false; continue; }
      if (ch == '\t' || ch == '\n')
      {
        rows.back().push_back(c);
        c.clear();
        at_start = true;
        if (ch == '\n') { rows.emplace_back(); }
        continue;
      }
      c += ch;
      at_start = false;
    }
    if (rows.back().empty() && c.empty()) { rows.pop_back(); }
    else { rows.back().push_back(c); }
    return rows;
  };
  for (const std::string name : {"scan_commands.tsv", "scan_results.tsv", "identification.tsv", "pooled_identification.tsv"})
  {
    const auto rows = csvRows(dir + "/" + name);
    const size_t n_lines = lines(dir + "/" + name).size();
    TEST_EQUAL(rows.size(), n_lines)
    bool widths_ok = !rows.empty();
    for (const auto& r : rows) { widths_ok = widths_ok && r.size() == rows.front().size(); }
    TEST_TRUE(widths_ok)
  }
  const auto csv_cmds = csvRows(dir + "/scan_commands.tsv");
  ABORT_IF(csv_cmds.size() != 4)
  TEST_EQUAL(csv_cmds[1][0], s_id)
  TEST_EQUAL(csv_cmds[2][0], m_id)
  TEST_EQUAL(csv_cmds[2][3], s_id)
  TEST_EQUAL(csv_cmds[3][0], one)

  // (3) TSVFile, the class tests' reader: every id-bearing cell equals the id the engine minted.
  auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
  ABORT_IF(cmds.rows.size() != 3)
  TEST_EQUAL(cell(cmds, cmds.rows[0], "tracking_id"), s_id)
  TEST_EQUAL(cell(cmds, cmds.rows[0], "scan_description"), s_id + "S")
  TEST_EQUAL(cell(cmds, cmds.rows[1], "tracking_id"), m_id)
  TEST_EQUAL(cell(cmds, cmds.rows[1], "parent_tracking_id"), s_id)
  TEST_EQUAL(cell(cmds, cmds.rows[1], "scan_description"), m_id + "R5.30533k@10")
  TEST_EQUAL(cell(cmds, cmds.rows[2], "tracking_id"), one)
  TEST_EQUAL(cell(cmds, cmds.rows[2], "scan_description"), one + "R5.30533k@10")
  auto results = TSVFile::parse(dir + "/scan_results.tsv");
  ABORT_IF(results.rows.size() != 1)
  TEST_EQUAL(cell(results, results.rows[0], "tracking_id"), m_id)
  TEST_EQUAL(cell(results, results.rows[0], "parent_tracking_id"), s_id)
  TEST_EQUAL(cell(results, results.rows[0], "child_ids"), c_id + " " + one)
  TEST_EQUAL(cell(results, results.rows[0], "winner_tracking_id"), w_id)
  auto ident = TSVFile::parse(dir + "/identification.tsv");
  ABORT_IF(ident.rows.size() != 1)
  TEST_EQUAL(cell(ident, ident.rows[0], "tracking_id"), m_id)
  auto pool = TSVFile::parse(dir + "/pooled_identification.tsv");
  ABORT_IF(pool.rows.size() != 1)
  TEST_EQUAL(cell(pool, pool.rows[0], "trigger_scan_id"), m_id)
  TEST_EQUAL(cell(pool, pool.rows[0], "contributing_scan_ids"), m_id + " " + one)

  // (4) ScanCommandJoin, i.e. FLASHDeconv: the ids parse and the MS2 joins to its survey.
  auto rows = ScanCommandJoin::parse(dir + "/scan_commands.tsv");
  TEST_EQUAL(rows.size(), 3)
  ABORT_IF(rows.count(m_id) == 0 || rows.count(s_id) == 0 || rows.count(one) == 0)
  TEST_EQUAL(rows[m_id].parent_tracking_id, s_id)
  TEST_EQUAL(rows[m_id].charge, 10)
  auto located = ScanCommandJoin::join({scanOf(40, 1, s_id + "S"), scanOf(41, 2, m_id + "R5.30533k@10", {531.790649})}, rows);
  ABORT_IF(located.count(41) == 0)
  TEST_EQUAL(located[41].parent_scan_number, 40)
}
END_SECTION

// A file written BEFORE the rule holds those ids raw, and every such file must keep parsing (ADR-0046 d9).
// "!" (id 8,837) begins and ends with '"' and would read as a quoted ! to a general unquote.
// Fails if ScanCommandJoin unquotes a raw 3-character id.
START_SECTION(scan_commands_parse_reads_raw_quote_ids_of_an_older_file)
{
  const Cells survey = surveyRow("\"!!", "267351474");
  const Cells ms2 = withCell(withCell(withCell(kM37, "tracking_id", "\"!\""), "parent_tracking_id", "\"!!"),
                             "scan_description", "\"!\"R5.30533k@10");
  auto rows = ScanCommandJoin::parse(writeCommandsFile("logcell_t03", {kCmdHeader, survey, ms2}));
  TEST_EQUAL(rows.size(), 2)
  ABORT_IF(rows.count("\"!\"") == 0 || rows.count("\"!!") == 0)
  TEST_EQUAL(rows["\"!\""].parent_tracking_id, std::string("\"!!"))
  auto located = ScanCommandJoin::join({scanOf(32, 1, "\"!!S"), scanOf(37, 2, "\"!\"R5.30533k@10", {531.790649})}, rows);
  ABORT_IF(located.count(37) == 0)
  TEST_EQUAL(located[37].parent_scan_number, 32)
}
END_SECTION

END_TEST
