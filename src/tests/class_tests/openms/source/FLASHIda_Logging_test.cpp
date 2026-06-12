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

#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <cstdio>

#include "FLASHIda_TestHelpers.h"

using namespace OpenMS;

namespace
{
  // Spectrum paths for these logging tests. The shared loaders, cycle drivers, and
  // JSON config builder now live in FLASHIda_TestHelpers.h (included above).
  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_standard.txt";
  const std::string ms2_tsv_path = "../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt";
}

START_TEST(FLASHIda_Logging, "$Id$")

/////////////////////////////////////////////////////////////

// Test 1: IDA Log contract -- write + parseFLASHIdaLog roundtrip
START_SECTION(ida_log_contract_roundtrip)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  // Use temp file for IDA log
  std::string ida_log_file = "test_ida_log_contract.log";
  std::remove(ida_log_file.c_str());

  std::string json = buildJsonWithRuntime(ida_log_file, "", "");
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Push all MS1 scans
  int total_commands = pushAllScans(&ida, ms1_scans);
  TEST_TRUE(total_commands > 0);

  // Parse the IDA log back using parseFLASHIdaLog
  auto parsed = FLASHIda::parseFLASHIdaLog(ida_log_file);

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

  // Cleanup
  std::remove(ida_log_file.c_str());
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

  std::string commands_file = "test_scan_commands.tsv";
  std::remove(commands_file.c_str());

  // Enable MS3 so we get MS3 commands in the TSV
  std::string json = buildJsonWithRuntime("", commands_file, "", true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Full MS1->MS2->MS3 cycle
  auto cycle = runFullCycle(&ida, ms1_scans, ms2_scans);
  TEST_TRUE(cycle.ms2_cmds.size() > 0);

  // Parse and verify TSV
  auto tsv = TSVFile::parse(commands_file);

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

  std::remove(commands_file.c_str());
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

  std::string results_file = "test_scan_results.tsv";
  std::remove(results_file.c_str());

  // Enable MS3 so we get MS3 result rows
  std::string json = buildJsonWithRuntime("", "", results_file, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Full MS1->MS2->MS3 cycle (including feeding MS3 back)
  auto cycle = runFullCycle(&ida, ms1_scans, ms2_scans);

  // Parse and verify
  auto tsv = TSVFile::parse(results_file);
  TEST_TRUE(tsv.colIndex("tracking_id") >= 0);
  TEST_TRUE(tsv.colIndex("resolve_ts") >= 0);
  TEST_TRUE(tsv.colIndex("duration_ms") >= 0);
  TEST_TRUE(tsv.colIndex("mass_count") >= 0);
  TEST_TRUE(tsv.colIndex("commands_pushed") >= 0);
  TEST_TRUE(tsv.colIndex("child_ids") >= 0);

  // Should have MS1, MS2, and (if MS3 commands were created) MS3 result rows
  // MS1 results from pushAllScans, MS2 from feeding back, MS3 from feeding back
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

  std::remove(results_file.c_str());
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

  std::string commands_file = "test_join_commands.tsv";
  std::string results_file = "test_join_results.tsv";
  std::remove(commands_file.c_str());
  std::remove(results_file.c_str());

  // Enable MS3 for full parent-child graph testing
  std::string json = buildJsonWithRuntime("", commands_file, results_file, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // Full MS1->MS2->MS3 cycle (with MS3 fed back)
  auto cycle = runFullCycle(&ida, ms1_scans, ms2_scans);

  // Drain any remaining commands (break on the idle-cycle AGC; the engine emits idle
  // AGC/MS1 commands indefinitely once the real queue is empty, so an unbounded drain
  // would loop forever and OOM).
  ScanCommand cmd;
  while (ida.getNextScanCommand(cmd) > 0) { if (cmd.is_agc) break; }

  // Parse both files
  auto cmd_tsv = TSVFile::parse(commands_file);
  auto res_tsv = TSVFile::parse(results_file);

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

  std::remove(commands_file.c_str());
  std::remove(results_file.c_str());
}
END_SECTION

// Test 5: Crash safety -- files are valid TSV after each operation (including MS3)
START_SECTION(crash_safety_valid_tsv)
{
  // Real CytC MS1+MS2 with MS3 enabled so the MS2 and MS3 crash-safety paths below are
  // actually reached; with MS3 off the MS3 block was unreachable and the MS2 block was a
  // no-failing-else conditional, so the section validated only the MS1 path.
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  std::string commands_file = "test_crash_commands.tsv";
  std::string results_file = "test_crash_results.tsv";
  std::remove(commands_file.c_str());
  std::remove(results_file.c_str());

  // Enable MS3 for full cycle
  std::string json = buildJsonWithRuntime("", commands_file, results_file, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // After constructor: headers should exist
  {
    auto cmd_tsv = TSVFile::parse(commands_file);
    TEST_TRUE(cmd_tsv.headers.size() > 0);
    auto res_tsv = TSVFile::parse(results_file);
    TEST_TRUE(res_tsv.headers.size() > 0);
  }

  // Push all MS1 scans, check results file is valid
  pushAllScans(&ida, ms1_scans);
  {
    auto res_tsv = TSVFile::parse(results_file);
    TEST_TRUE(res_tsv.rows.size() >= 1);
    for (const auto& row : res_tsv.rows)
      TEST_EQUAL(row.size(), res_tsv.headers.size());
  }

  // Dequeue all MS2 commands (break on the idle-cycle AGC); commands file must stay valid
  // and at least one MS2 command must have been produced (positive expectation, not a
  // no-failing-else conditional).
  std::vector<ScanCommand> ms2_cmds;
  ScanCommand cmd{};
  while (ida.getNextScanCommand(cmd) > 0)
  {
    if (cmd.msn_level == 2) ms2_cmds.push_back(cmd);
    if (cmd.is_agc) break;
  }
  TEST_TRUE(ms2_cmds.size() > 0);
  {
    auto cmd_tsv = TSVFile::parse(commands_file);
    TEST_TRUE(cmd_tsv.rows.size() >= 1);
    for (const auto& row : cmd_tsv.rows)
      TEST_EQUAL(row.size(), cmd_tsv.headers.size());
  }

  // Feed every MS2 result back, check results file is still valid
  for (const auto& m : ms2_cmds)
    ida.processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                    (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 2,
                    m.scan_description);
  {
    auto res_tsv = TSVFile::parse(results_file);
    for (const auto& row : res_tsv.rows)
      TEST_EQUAL(row.size(), res_tsv.headers.size());
  }

  // Dequeue MS3 commands (break on idle AGC); commands file valid and MS3 must have fired.
  std::vector<ScanCommand> ms3_cmds;
  while (ida.getNextScanCommand(cmd) > 0)
  {
    if (cmd.msn_level == 3) ms3_cmds.push_back(cmd);
    if (cmd.is_agc) break;
  }
  TEST_TRUE(ms3_cmds.size() > 0);
  {
    auto cmd_tsv = TSVFile::parse(commands_file);
    for (const auto& row : cmd_tsv.rows)
      TEST_EQUAL(row.size(), cmd_tsv.headers.size());
  }

  // Feed every MS3 result back, check results file is still valid
  for (const auto& m : ms3_cmds)
    ida.processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                    (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 3,
                    m.scan_description);
  {
    auto res_tsv = TSVFile::parse(results_file);
    for (const auto& row : res_tsv.rows)
      TEST_EQUAL(row.size(), res_tsv.headers.size());
  }

  std::remove(commands_file.c_str());
  std::remove(results_file.c_str());
}
END_SECTION

END_TEST
