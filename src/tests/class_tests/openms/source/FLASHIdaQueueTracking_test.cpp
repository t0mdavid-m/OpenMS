// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

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
      "max_mass_count": [1],
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "MS3AllCharges": true,
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

// P3-U05: encodeBase36_ correctness
START_SECTION(encodeBase36_correctness)
{
  TEST_STRING_EQUAL(FLASHIda::encodeBase36ForTest(0), "0000")
  TEST_STRING_EQUAL(FLASHIda::encodeBase36ForTest(35), "000z")
  TEST_STRING_EQUAL(FLASHIda::encodeBase36ForTest(36), "0010")
  TEST_STRING_EQUAL(FLASHIda::encodeBase36ForTest(1679615), "zzzz")
  // Mid-range value
  TEST_STRING_EQUAL(FLASHIda::encodeBase36ForTest(1296), "0100")
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
    seen.insert(id);
    prev = id;
  }
  delete ida;
}
END_SECTION

// P3-U07: empty queue returns MS1 scan
START_SECTION(empty_queue_returns_ms1)
{
  FLASHIda* ida = createTestInstance();
  ScanCommand cmd{};
  int result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_EQUAL(cmd.is_agc, 0)
  // Analyzer should be "Orbitrap" from config
  TEST_STRING_EQUAL(std::string(cmd.analyzer), "Orbitrap")
  // Resolution from config
  TEST_EQUAL(cmd.orbitrap_resolution, 120000)
  delete ida;
}
END_SECTION

// P3-U08: queue priority dequeue order (deferred to Phase 4)
START_SECTION(queue_priority_dequeue_order_stub)
{
  // Phase 4 will test enqueueing commands at different priorities
  // and verifying dequeue order matches priority.
  NOT_TESTABLE
}
END_SECTION

// P3-U09: AGC scan is dequeued first (deferred to Phase 4)
START_SECTION(agc_scan_is_dequeued_first_stub)
{
  // Phase 4 will test that AGC scans have highest priority.
  NOT_TESTABLE
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
  TEST_EQUAL(result, 1)
  TEST_EQUAL(cmd.msn_level, 1)
  delete ida;
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
