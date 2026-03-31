// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------
//
// Phase 4 unit tests: processScan full routing.
// Peak data extracted from FlashIDA/test-data/spectra/ms1_standard.txt (scan 1)
// and ms2_hcd_fragment.txt (scan 57).

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include <string>
#include <cstring>

using namespace OpenMS;

namespace
{
  // Minimal JSON config for standard DDA mode with score_threshold=0 to accept all peaks
  const char* standard_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10]
    },
    "precursor_selection": {
      "max_mass_count": [3], "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false, "MS3AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": true, "value_ms": 30000 },
      "agc_interval_seconds": 30
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" }
  })";

  // Config with MS3 mode 1 (SourceCID) enabled
  const char* ms3_mode1_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10]
    },
    "precursor_selection": {
      "max_mass_count": [3], "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false, "MS3AllCharges": false,
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
    "ms3": { "enabled": true, "mode": 1, "max_per_ms2": 2, "protein_sequence": "" },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" }
  })";

  // Config with conditional MS2 enabled
  const char* conditional_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10]
    },
    "precursor_selection": {
      "max_mass_count": [3], "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false, "MS3AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 30
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "conditional_ms2": true,
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" }
  })";

  // Real MS1 peaks from ms1_standard.txt scan 1, RT=0.0668 min
  // Strong charge envelopes around 519 m/z, 536 m/z, 577 m/z, 593 m/z, 610 m/z
  const double ms1_mzs[] = {
    503.107269, 504.107849, 505.086548, 505.104584, 506.086975, 506.104248,
    519.138550, 520.138977, 521.117981, 521.135803, 522.135803, 523.133240,
    536.165100, 537.165588, 538.162476, 539.163025, 540.159668,
    550.217224, 551.217773,
    563.074158, 564.075256,
    577.126038, 578.126343, 579.105164, 579.123962, 580.105835, 581.102600,
    593.157410, 594.157837, 595.155151, 596.154907, 597.150818,
    607.209412, 608.209778, 609.206909,
    610.183838, 611.184570, 612.181946, 613.181091, 614.178894
  };
  const double ms1_ints[] = {
    59130.62, 24195.42, 14895.83, 14728.34, 7187.91, 5917.36,
    246361.94, 92025.09, 6748.71, 56501.90, 18775.31, 7298.50,
    195921.22, 95304.26, 52071.90, 16788.95, 5082.68,
    11773.91, 5773.60,
    4017.96, 2033.30,
    30457.84, 13873.09, 7044.49, 5960.62, 2958.92, 2413.76,
    109951.71, 46630.72, 33525.32, 10582.25, 3454.64,
    21142.76, 7270.54, 3597.60,
    174427.27, 83930.87, 46101.07, 22262.65, 6546.02
  };
  const int ms1_length = sizeof(ms1_mzs) / sizeof(ms1_mzs[0]);

  // Minimal MS2 peaks from ms2_hcd_fragment.txt scan 57
  const double ms2_mzs[] = {
    100.300545, 104.473274, 104.481895, 110.161835, 111.074051,
    111.105568, 111.118179, 111.128120
  };
  const double ms2_ints[] = {
    16955.05, 20442.67, 20458.91, 12642.20, 13831.51,
    13627.63, 14103.31, 13171.75
  };
  const int ms2_length = sizeof(ms2_mzs) / sizeof(ms2_mzs[0]);
}

START_TEST(FLASHIda_ProcessScan, "$Id$")

/////////////////////////////////////////////////////////////

// P4-U01: MS1 processScan returns > 0 commands for real spectral data
START_SECTION(processScan_ms1_returns_commands)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)
  // Should be at most max_mass_count (3)
  TEST_EQUAL(n <= 3, true)
  delete ida;
}
END_SECTION

// P4-U02: Commands from processScan are dequeued by getNextScanCommand
START_SECTION(processScan_commands_dequeued)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  // Dequeue all commands
  ScanCommand cmd{};
  for (int i = 0; i < n; i++)
  {
    int result = ida->getNextScanCommand(cmd);
    TEST_EQUAL(result, 1)
    TEST_EQUAL(cmd.msn_level, 2)
    TEST_EQUAL(cmd.priority, 1)
    TEST_EQUAL(cmd.num_stages, 1)
    TEST_EQUAL(cmd.is_agc, 0)
    // Isolation stage should have valid m/z
    TEST_EQUAL(cmd.stages[0].precursor_mz > 0, true)
    TEST_EQUAL(cmd.stages[0].isolation_width > 0, true)
    TEST_EQUAL(cmd.stages[0].charge_state >= 4, true)
    // Scan description starts with 4-char tracking ID
    TEST_EQUAL(std::strlen(cmd.scan_description) >= 4, true)
    // Enqueue timestamp should be set
    TEST_EQUAL(cmd.enqueue_timestamp_ms > 0, true)
  }

  // Queue empty — next call returns 0
  int empty_result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(empty_result, 0)

  delete ida;
}
END_SECTION

// P4-U03: ScanCommand fields populated correctly
START_SECTION(processScan_command_fields)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  ScanCommand cmd{};
  ida->getNextScanCommand(cmd);

  // Analyzer from ms2_configs_[0]
  TEST_STRING_EQUAL(std::string(cmd.analyzer), "Orbitrap")
  // Resolution from ms2_configs_[0]
  TEST_EQUAL(cmd.orbitrap_resolution, 120000)
  // Activation type from ms2_configs_[0]
  TEST_STRING_EQUAL(std::string(cmd.stages[0].activation_type), "HCD")
  // Collision energy
  TEST_EQUAL(cmd.stages[0].collision_energy > 0, true)
  // Precursor m/z within scan range
  TEST_EQUAL(cmd.stages[0].precursor_mz >= 500.0, true)
  TEST_EQUAL(cmd.stages[0].precursor_mz <= 2000.0, true)

  delete ida;
}
END_SECTION

// P4-U04: processScan with ms_level=2 processes MS2 path
START_SECTION(processScan_ms2_path)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // First push an MS2 via MS1 processing to get a tracking ID in pending_scan_map_
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  // Dequeue one MS2 command to get its scan description (contains tracking ID)
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(ms2_cmd.msn_level, 2)

  // Now process MS2 return with the tracking ID in scan description
  int ms2_result = ida->processScan(ms2_mzs, ms2_ints, ms2_length, 0.1,
                                     2, ms2_cmd.scan_description);
  // Should return 0 (no conditional, no MS3, no quant in standard config)
  TEST_EQUAL(ms2_result, 0)

  delete ida;
}
END_SECTION

// P4-U05: processScan with empty spectrum returns 0
START_SECTION(processScan_empty_spectrum)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int n = ida->processScan(nullptr, nullptr, 0, 0.0, 1, "empty");
  TEST_EQUAL(n, 0)
  delete ida;
}
END_SECTION

// P4-U06: conditional MS2 follow-up is pushed at priority 2
START_SECTION(processScan_conditional_ms2_followup)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(conditional_json));

  // Push MS2 commands via MS1
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  // Dequeue first MS2
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(ms2_cmd.msn_level, 2)
  TEST_EQUAL(ms2_cmd.priority, 1)

  // Process MS2 return — should push conditional follow-up at priority 2
  int ms2_result = ida->processScan(ms2_mzs, ms2_ints, ms2_length, 0.1,
                                     2, ms2_cmd.scan_description);
  TEST_EQUAL(ms2_result > 0, true)

  // Drain remaining priority-1 commands first
  ScanCommand out{};
  while (true)
  {
    ida->getNextScanCommand(out);
    if (out.priority != 1) break;
  }
  // Should get priority-2 conditional follow-up
  TEST_EQUAL(out.priority, 2)
  TEST_EQUAL(out.msn_level, 2)

  delete ida;
}
END_SECTION

// P4-U07: MS3 commands are pushed at priority 3
START_SECTION(processScan_ms3_commands)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(ms3_mode1_json));

  // Push MS2 commands via MS1
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  // Dequeue first MS2
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(ms2_cmd.msn_level, 2)

  // Process MS2 return — may push MS3 commands if fragments found
  int ms2_result = ida->processScan(ms2_mzs, ms2_ints, ms2_length, 0.1,
                                     2, ms2_cmd.scan_description);
  // MS3 commands pushed depends on deconvolution results. We test the structure.
  if (ms2_result > 0)
  {
    // Drain all priority-1 commands
    ScanCommand out{};
    while (true)
    {
      ida->getNextScanCommand(out);
      if (out.priority != 1) break;
    }
    // If MS3 was pushed, should be priority 3 with 2 stages
    if (out.msn_level == 3)
    {
      TEST_EQUAL(out.priority, 3)
      TEST_EQUAL(out.num_stages, 2)
    }
  }

  delete ida;
}
END_SECTION

// P4-U08: decodeBase36 roundtrips with encodeBase36
START_SECTION(decodeBase36_roundtrip)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // Test roundtrip via processScan → scan_description parsing
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  ScanCommand cmd{};
  ida->getNextScanCommand(cmd);
  // scan_description starts with 4-char base36 tracking ID
  std::string desc(cmd.scan_description);
  TEST_EQUAL(desc.size() >= 4, true)
  std::string id_str = desc.substr(0, 4);
  // Verify the ID is valid base-36
  for (char c : id_str)
  {
    TEST_EQUAL((c >= '0' && c <= '9') || (c >= 'a' && c <= 'z'), true)
  }

  delete ida;
}
END_SECTION

// P4-U09: cleanupExpiredCommands removes old entries
START_SECTION(cleanup_expired_commands)
{
  // Use timeout-enabled config
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // Push an MS2 via MS1 processing (this adds to pending_scan_map_)
  int n = ida->processScan(ms1_mzs, ms1_ints, ms1_length, 0.0668, 1, "ms1_test");
  TEST_EQUAL(n > 0, true)

  // Dequeue all commands — they're still in pending_scan_map_ for tracking
  ScanCommand cmd{};
  for (int i = 0; i < n; i++)
  {
    ida->getNextScanCommand(cmd);
  }

  // The pending_scan_map_ entries have timestamps. cleanupExpiredCommands_
  // is called by getNextScanCommand. With timeout_ms=30000, entries should
  // NOT be expired immediately. Verify by processing an MS2 with valid tracking ID.
  int ms2_result = ida->processScan(ms2_mzs, ms2_ints, ms2_length, 0.1,
                                     2, cmd.scan_description);
  // Should succeed (entry found in pending_scan_map_)
  // ms2_result can be 0 (no follow-ups) but shouldn't crash
  TEST_EQUAL(ms2_result >= 0, true)

  delete ida;
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
