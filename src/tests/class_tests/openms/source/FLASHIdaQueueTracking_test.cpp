// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include <set>
#include <string>

using namespace OpenMS;

namespace
{
  // Minimal JSON config matching Phase 1 schema — just enough to construct FLASHIda
  const char* minimal_json_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0
    },
    "ms_settings": {
      "ms1": {
        "analyzer": "Orbitrap",
        "first_mass": 500,
        "last_mass": 2000,
        "resolution": 120000,
        "agc_target": 800000,
        "max_it": 246
      },
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "ETD",
          "collision_energy": 0,
          "resolution": 120000
        }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "exploration": {
      "enabled": false,
      "max_depth": 1,
      "max_variants": 5
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 1 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  FLASHIda* createTestInstance()
  {
    // const_cast needed because constructor takes char* (legacy strtok compatibility)
    return new FLASHIda(const_cast<char*>(minimal_json_config));
  }
}

START_TEST(FLASHIdaQueueTracking, "$Id$")

/////////////////////////////////////////////////////////////

// P3-U05: encodeTracking_ correctness (base-94, 3-char, alphabet '!' to '~')
START_SECTION(encodeTracking_correctness)
{
  // encode(0) = "!!!" (all zeros in base-94)
  TEST_STRING_EQUAL(ScanCommandQueue::encode(0), "!!!")
  // encode(1) = "!!"" (last char advances to index 1 = '"')
  TEST_STRING_EQUAL(ScanCommandQueue::encode(1), "!!\"")
  // encode(94) = "!"!" (middle char advances to index 1 = '"')
  TEST_STRING_EQUAL(ScanCommandQueue::encode(94), "!\"!")
  // encode(830583) = "~~~" (max value: 94^3 - 1)
  TEST_STRING_EQUAL(ScanCommandQueue::encode(830583), "~~~")
  // encode(8836) = '"!!' (94*94 = 8836, first char advances to index 1 = '"')
  TEST_STRING_EQUAL(ScanCommandQueue::encode(8836), "\"!!")
}
END_SECTION

// P3-U06: tracking IDs are sequential and unique
START_SECTION(tracking_ids_sequential_unique)
{
  FLASHIda* ida = createTestInstance();
  std::set<int> seen;
  int prev = -1;
  for (int i = 0; i < 10000; ++i)
  {
    int id = ida->getNextTrackingId();
    TEST_EQUAL(seen.count(id), 0)
    TEST_EQUAL(id > prev, true)
    // Verify ID stays within base-94 3-char range (94^3 = 830584)
    TEST_EQUAL(id < 830584, true)
    seen.insert(id);
    prev = id;
  }
  delete ida;
}
END_SECTION

// P3-U07: empty queue returns idle AGC (never returns 0)
START_SECTION(empty_queue_returns_idle_agc)
{
  FLASHIda* ida = createTestInstance();
  ScanCommand cmd{};
  int result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)  // Queue empty, idle cycle returns AGC
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.is_agc, 1)
  delete ida;
}
END_SECTION

// P3-U08: queue priority dequeue order
START_SECTION(queue_priority_dequeue_order)
{
  FLASHIda* ida = createTestInstance();

  // Push commands at different priorities (out of order)
  ScanCommand cmd3{};
  cmd3.msn_level = 3;
  cmd3.priority = 3;
  cmd3.scan_id = 100;
  ida->pushCommandForTest(cmd3);

  ScanCommand cmd1{};
  cmd1.msn_level = 2;
  cmd1.priority = 1;
  cmd1.scan_id = 101;
  ida->pushCommandForTest(cmd1);

  ScanCommand cmd2{};
  cmd2.msn_level = 2;
  cmd2.priority = 2;
  cmd2.scan_id = 102;
  ida->pushCommandForTest(cmd2);

  // Dequeue should return priority 1, then 2, then 3
  ScanCommand out{};
  ida->getNextScanCommand(out);
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.scan_id, 101)  // priority 1
  TEST_EQUAL(out.priority, 1)

  ida->getNextScanCommand(out);
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.scan_id, 102)  // priority 2
  TEST_EQUAL(out.priority, 2)

  ida->getNextScanCommand(out);
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.scan_id, 100)  // priority 3
  TEST_EQUAL(out.priority, 3)

  // Queue empty — idle cycle returns AGC (never returns 0)
  int idle_result = ida->getNextScanCommand(out);
  TEST_EQUAL(idle_result, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.is_agc, 1)

  delete ida;
}
END_SECTION

// P3-U09: AGC scan is dequeued first
START_SECTION(agc_scan_is_dequeued_first)
{
  // Create FLASHIda with agc_interval_seconds = 0 (always needs AGC)
  const char* agc_config = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10]
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
      "ms2": [{ "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 0
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 1 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  FLASHIda* ida = new FLASHIda(const_cast<char*>(agc_config));

  // Push an MS2 command into the queue
  ScanCommand ms2_cmd{};
  ms2_cmd.msn_level = 2;
  ms2_cmd.priority = 2;
  ms2_cmd.scan_id = 200;
  ida->pushCommandForTest(ms2_cmd);

  // First dequeue: with agc_interval_seconds=0, AGC is always needed.
  // AGC should be returned before any queued command.
  ScanCommand out{};
  ida->getNextScanCommand(out);
  TEST_EQUAL(out.is_agc, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_STRING_EQUAL(std::string(out.analyzer), "IonTrap")
  TEST_EQUAL(out.msn_level, 1)

  delete ida;
}
END_SECTION

// P3-U10: timeout cleanup does not crash when disabled
START_SECTION(timeout_cleanup_no_crash)
{
  FLASHIda* ida = createTestInstance();
  // timeout_enabled_ defaults to false from config above.
  // Calling getNextScanCommand (which calls cleanupExpiredCommands_) should not crash.
  ScanCommand cmd{};
  int result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)  // Queue empty, idle cycle returns AGC
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.is_agc, 1)
  delete ida;
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
