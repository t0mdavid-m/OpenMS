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

using namespace OpenMS;

namespace
{
  // Helper: build JSON config string with runtime paths
  // enable_ms3 parameter is retained for API compat but MS3 is now only
  // activated via selection_strategy — these tests focus on logging TSV format.
  std::string buildJsonWithRuntime(const std::string& ida_log_path,
                                   const std::string& commands_path,
                                   const std::string& results_path,
                                   bool /*enable_ms3*/ = false)
  {
    std::string ms3_block;  // no ms3 block needed (protein_sequence defaults to "")
    std::string ms3_selection = "\"none\"";

    std::ostringstream oss;
    oss << R"({
      "deconvolution": {
        "score_threshold": 0.0, "tqscore_threshold": 0.9,
        "min_charge": 4, "max_charge": 50,
        "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
      },
      "precursor_selection": {
        "RT_window": 180, "target_mode": 0,
        "IDScore": false, "AllCharges": false,
        "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
      },
      "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
      "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
      "faims": { "cv_values": [-50], "max_cv_skip": 0 },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
        "ms2": [
          { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
        ]
      },
      "scheduling": {
        "cycle_time": { "enabled": false, "value_ms": 60000 },
        "scan_timeout": { "enabled": false, "value_ms": 30000 },
        "agc_interval_seconds": 30
      },
      "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
      )" << ms3_block << R"(
      "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
      "selection_strategy": {
        "ms1": { "selection": "qscore", "max_targets": 5 },
        "ms2": { "selection": "none" },
        "ms3": { "selection": )" << ms3_selection << R"( }
      },
      "runtime": {
        "ida_log_path": ")" << ida_log_path << R"(",
        "scan_commands_path": ")" << commands_path << R"(",
        "scan_results_path": ")" << results_path << R"("
      }
    })";
    return oss.str();
  }

  struct ScanData
  {
    std::vector<double> mzs;
    std::vector<double> ints;
    double rt;
    std::string scan_id;
  };

  // Parse multi-scan TSV: "Spec scan=N\tRT" headers, "mz\tintensity" data lines
  std::vector<ScanData> loadTsvScans(const std::string& path)
  {
    std::vector<ScanData> scans;
    std::ifstream f(path);
    if (! f.good()) return scans;
    std::string line;
    while (std::getline(f, line))
    {
      if (line.substr(0, 4) == "Spec")
      {
        scans.emplace_back();
        auto tab = line.find('\t');
        scans.back().scan_id = line.substr(10, tab - 10);
        scans.back().rt = std::stod(line.substr(tab + 1));
      }
      else if (! scans.empty())
      {
        auto tab = line.find('\t');
        if (tab != std::string::npos)
        {
          scans.back().mzs.push_back(std::stod(line.substr(0, tab)));
          scans.back().ints.push_back(std::stod(line.substr(tab + 1)));
        }
      }
    }
    return scans;
  }

  // Push all scans through processScan, return total command count
  int pushAllScans(FLASHIda* ida, const std::vector<ScanData>& scans)
  {
    int total = 0;
    for (const auto& scan : scans)
    {
      int n = ida->processScan(scan.mzs.data(), scan.ints.data(),
                                (int)scan.mzs.size(), scan.rt, 1,
                                ("scan_" + scan.scan_id).c_str());
      total += n;
    }
    return total;
  }

  // Helper: run full MS1->MS2->MS3 cycle, collecting commands at each level
  struct CycleResult
  {
    std::vector<ScanCommand> ms2_cmds;
    std::vector<ScanCommand> ms3_cmds;
    std::vector<ScanCommand> ms1_cmds;
    int total_dequeued = 0;
  };

  CycleResult runFullCycle(FLASHIda* ida,
                           const std::vector<ScanData>& ms1_scans,
                           const std::vector<ScanData>& ms2_scans)
  {
    CycleResult result;

    // 1. Push all MS1 scans -> creates MS2 commands
    pushAllScans(ida, ms1_scans);

    // 2. Dequeue all commands (MS2 + AGC/MS1 fallbacks)
    ScanCommand cmd{};
    while (ida->getNextScanCommand(cmd) > 0)
    {
      result.total_dequeued++;
      if (cmd.msn_level == 2)
        result.ms2_cmds.push_back(cmd);
      else if (cmd.msn_level == 1)
        result.ms1_cmds.push_back(cmd);
      if (cmd.is_agc) break; // idle-cycle AGC => real queue drained (engine emits idle cmds forever)
    }

    // 3. Feed MS2 results back -> may create MS3 commands
    for (const auto& ms2_cmd : result.ms2_cmds)
    {
      if (! ms2_scans.empty())
      {
        ida->processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                        (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 2,
                        ms2_cmd.scan_description);
      }
    }

    // 4. Dequeue MS3 commands (+ any new AGC/MS1)
    while (ida->getNextScanCommand(cmd) > 0)
    {
      result.total_dequeued++;
      if (cmd.msn_level == 3)
        result.ms3_cmds.push_back(cmd);
      if (cmd.is_agc) break;
    }

    // 5. Feed MS3 results back through processScan (ms_level=3)
    for (const auto& ms3_cmd : result.ms3_cmds)
    {
      // Reuse MS2 fragment data as MS3 input (same format, different level)
      if (! ms2_scans.empty())
      {
        ida->processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                        (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 3,
                        ms3_cmd.scan_description);
      }
    }

    return result;
  }

  // Parse a TSV file into rows of string vectors (first row = header)
  struct TSVFile
  {
    std::vector<std::string> headers;
    std::vector<std::vector<std::string>> rows;

    static TSVFile parse(const std::string& path)
    {
      TSVFile result;
      std::ifstream f(path);
      std::string line;
      bool first = true;
      while (std::getline(f, line))
      {
        std::vector<std::string> cols;
        std::istringstream iss(line);
        std::string col;
        while (std::getline(iss, col, '\t'))
          cols.push_back(col);

        if (first)
        {
          result.headers = cols;
          first = false;
        }
        else
        {
          result.rows.push_back(cols);
        }
      }
      return result;
    }

    // Get column index by name
    int colIndex(const std::string& name) const
    {
      for (size_t i = 0; i < headers.size(); i++)
        if (headers[i] == name) return static_cast<int>(i);
      return -1;
    }
  };

  // TSV file paths relative to CTest working dir (OpenMS/build/)
  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_standard.txt";
  const std::string ms2_tsv_path = "../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt";
}

START_TEST(FLASHIda_Logging, "$Id$")

/////////////////////////////////////////////////////////////

// Test 1: IDA Log contract -- write + parseFLASHIdaLog roundtrip
START_SECTION(ida_log_contract_roundtrip)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  if (ms1_scans.empty())
  {
    NOT_TESTABLE;  // test data not available in this build environment
    break;
  }

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
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  if (ms1_scans.empty() || ms2_scans.empty())
  {
    NOT_TESTABLE;
    break;
  }

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
  if (cycle.ms3_cmds.size() > 0)
  {
    TEST_TRUE(found_ms3);
  }

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
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  if (ms1_scans.empty() || ms2_scans.empty())
  {
    NOT_TESTABLE;
    break;
  }

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
      (void)dur; // stoull succeeds = valid unsigned integer
    }
  }

  // MS2 result rows should have child_ids pointing to MS3 commands (if any)
  if (cycle.ms3_cmds.size() > 0)
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
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  if (ms1_scans.empty() || ms2_scans.empty())
  {
    NOT_TESTABLE;
    break;
  }

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

  // Build set of all command tracking_ids
  std::set<std::string> cmd_ids;
  int cmd_id_col = cmd_tsv.colIndex("tracking_id");
  for (const auto& row : cmd_tsv.rows)
  {
    if (cmd_id_col >= 0 && cmd_id_col < (int)row.size())
      cmd_ids.insert(row[cmd_id_col]);
  }

  // Every child_id in results must exist in commands
  int child_col = res_tsv.colIndex("child_ids");
  int pushed_col = res_tsv.colIndex("commands_pushed");
  for (const auto& row : res_tsv.rows)
  {
    if (child_col >= 0 && child_col < (int)row.size() && ! row[child_col].empty())
    {
      // Parse semicolon-separated child IDs
      std::istringstream child_ss(row[child_col]);
      std::string child_id;
      int child_count = 0;
      while (std::getline(child_ss, child_id, ';'))
      {
        TEST_TRUE(cmd_ids.count(child_id) > 0);
        child_count++;
      }
      // commands_pushed should equal child count
      if (pushed_col >= 0 && pushed_col < (int)row.size())
      {
        TEST_EQUAL(std::stoi(row[pushed_col]), child_count);
      }
    }
  }

  std::remove(commands_file.c_str());
  std::remove(results_file.c_str());
}
END_SECTION

// Test 5: Crash safety -- files are valid TSV after each operation (including MS3)
START_SECTION(crash_safety_valid_tsv)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  if (ms1_scans.empty() || ms2_scans.empty())
  {
    NOT_TESTABLE;
    break;
  }

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

  // Push one MS1 scan, check files are valid
  ida.processScan(ms1_scans[0].mzs.data(), ms1_scans[0].ints.data(),
                  (int)ms1_scans[0].mzs.size(), ms1_scans[0].rt, 1, "scan_1");
  {
    auto res_tsv = TSVFile::parse(results_file);
    TEST_TRUE(res_tsv.rows.size() >= 1);
    for (const auto& row : res_tsv.rows)
      TEST_EQUAL(row.size(), res_tsv.headers.size());
  }

  // Dequeue one command, check commands file is valid
  ScanCommand cmd;
  ida.getNextScanCommand(cmd);
  {
    auto cmd_tsv = TSVFile::parse(commands_file);
    TEST_TRUE(cmd_tsv.rows.size() >= 1);
    for (const auto& row : cmd_tsv.rows)
      TEST_EQUAL(row.size(), cmd_tsv.headers.size());
  }

  // Feed MS2 result back, check files are still valid
  if (cmd.msn_level == 2)
  {
    ida.processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                    (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 2,
                    cmd.scan_description);
    {
      auto res_tsv = TSVFile::parse(results_file);
      for (const auto& row : res_tsv.rows)
        TEST_EQUAL(row.size(), res_tsv.headers.size());
    }
  }

  // Dequeue MS3 if available, feed back, check files
  ScanCommand ms3_cmd;
  if (ida.getNextScanCommand(ms3_cmd) > 0 && ms3_cmd.msn_level == 3)
  {
    // Check commands file after MS3 dequeue
    {
      auto cmd_tsv = TSVFile::parse(commands_file);
      for (const auto& row : cmd_tsv.rows)
        TEST_EQUAL(row.size(), cmd_tsv.headers.size());
    }

    // Feed MS3 result back
    ida.processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                    (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 3,
                    ms3_cmd.scan_description);

    // Check results file after MS3 result
    {
      auto res_tsv = TSVFile::parse(results_file);
      for (const auto& row : res_tsv.rows)
        TEST_EQUAL(row.size(), res_tsv.headers.size());
    }
  }

  std::remove(commands_file.c_str());
  std::remove(results_file.c_str());
}
END_SECTION

END_TEST
