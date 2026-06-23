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
  auto parsed = IdaLogger::parseFLASHIdaLog(ida_log_file);

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

  // Full MS1->MS2->MS3 cycle (with MS3 fed back). runFullCycle already drains the queue to idle,
  // so there is no leftover command to drain here -- the previous trailing
  // `while (getNextScanCommand) { if (is_agc) break; }` loop was redundant (and is exactly the
  // unbounded-drain hazard the harness exists to remove). The TSVs are written incrementally
  // during the cycle, so the parses below see every command/result the cycle produced.
  auto cycle = runFullCycle(&ida, ms1_scans, ms2_scans);

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

  std::remove(commands_file.c_str());
  std::remove(results_file.c_str());
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

  std::string commands_file = "test_crash_commands.tsv";
  std::string results_file = "test_crash_results.tsv";
  std::remove(commands_file.c_str());
  std::remove(results_file.c_str());

  // Enable MS3 for full cycle
  std::string json = buildJsonWithRuntime("", commands_file, results_file, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // After constructor: headers should exist (pre-cycle: no scans driven yet)
  {
    auto cmd_tsv = TSVFile::parse(commands_file);
    TEST_TRUE(cmd_tsv.headers.size() > 0);
    auto res_tsv = TSVFile::parse(results_file);
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
    auto res_tsv = TSVFile::parse(results_file);
    TEST_TRUE(res_tsv.rows.size() >= 1);
    for (const auto& row : res_tsv.rows)
      TEST_EQUAL(row.size(), res_tsv.headers.size());
  }

  // commands file is valid: >=1 row and every row has the header column count (no torn writes)
  {
    auto cmd_tsv = TSVFile::parse(commands_file);
    TEST_TRUE(cmd_tsv.rows.size() >= 1);
    for (const auto& row : cmd_tsv.rows)
      TEST_EQUAL(row.size(), cmd_tsv.headers.size());
  }

  std::remove(commands_file.c_str());
  std::remove(results_file.c_str());
}
END_SECTION

END_TEST
