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
#include <iostream>  // the two lock-split guards below report their counts
#include <thread>    // concurrent_drain_writes_one_wellformed_row_per_call
#include <map>       // precursor_id inheritance check
#include <set>       // tracking-id uniqueness check
#include <utility>   // std::pair

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
  // Today analysis_mutex_ serialises the three write sites, so this passes. It is expected to go RED
  // the moment those guards are deleted and BEFORE the per-stream logger mutexes land -- that
  // intermediate state is deliberate, and this test is the evidence those mutexes are load-bearing
  // rather than decorative.
  //
  // No processScan, no spectra, no id echo. The MS1 admission gate makes fabricated ids useless, but
  // the DRAIN needs none of that -- it manufactures its own work. That is what makes this cheap and
  // immune to fixture drift.
  //
  // THE INVARIANT IS EXACT, not "roughly right": every path through getNextScanCommand writes
  // exactly one row and returns 1. Step 1 (:668), Step 4 (:715) and Step 5 (:759) each write once;
  // Step 2 pushes and falls through rather than returning; Step 3 writes nothing. So rows == calls,
  // and a torn row merges two lines and DROPS the count -- deterministic detection of the common
  // tearing mode rather than a probabilistic one.
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

  // Exactly one row per call. Fewer means rows were merged (torn writes); more would mean a path
  // started writing twice.
  TEST_EQUAL((int)tsv.rows.size(), kTotalCalls);

  // Every row has the full column count. A row torn mid-write tokenizes short.
  int wellformed = 0;
  for (const auto& row : tsv.rows)
    if (row.size() == tsv.headers.size()) wellformed++;
  TEST_EQUAL(wellformed, (int)tsv.rows.size());

  // Tracking ids are allocated under queue_mutex_, so they must all be distinct even across threads.
  // This catches a torn row that happens to tokenize to the right width by splicing two half-rows.
  std::set<std::string> ids;
  for (const auto& row : tsv.rows) ids.insert(cell(tsv, row, "tracking_id"));
  TEST_EQUAL((int)ids.size(), (int)tsv.rows.size());

  std::cout << "[DRAIN-CONCURRENCY] calls=" << kTotalCalls << " rows=" << tsv.rows.size()
            << " wellformed=" << wellformed << " unique_ids=" << ids.size() << std::endl;
}
END_SECTION

END_TEST
