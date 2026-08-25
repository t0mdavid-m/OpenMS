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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include "FLASHIda_TestAccess.h"   // FLASHIdaTestAccess::push (private-state access)

#include <chrono>
#include <set>
#include <string>
#include <thread>  // std::this_thread::sleep_for -- reaching the interval-AGC path deterministically

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
    "tol": [
      10,
      10,
      10
    ]
  },
  "flashtnt": {
    "min_length": 3,
    "max_length": 8
  },
  "faims": {
    "cv_values": [],
    "max_cv_skip": 0
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": false,
      "value_ms": 30000
    },
    "agc_interval_seconds": 9999999
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 1
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
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
    "ms2": {
      "analyzer": "Orbitrap",
      "activation": "ETD",
      "collision_energy": 0,
      "reaction_time": 10.0,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)";

  FLASHIda* createTestInstance()
  {
    // const_cast needed because constructor takes char* (legacy strtok compatibility)
    return new FLASHIda(const_cast<char*>(minimal_json_config));
  }

  // Used by the two builder sections at the end of this file. It lives here, at file scope, and
  // not next to them: START_TEST expands to `int main(...) {` (ClassTest.h:453), so everything
  // after it is inside a function body and a namespace there is ill-formed (C2870).
  //
  // ms1 carries a full source region and microscans 4; ms2 has TWO configs whose agc_target and
  // max_it deliberately differ, which is what makes the buildMS2 section non-vacuous.
  const char* scan_sourcing_config = R"({
  "deconvolution": {
    "min_charge": 4,
    "max_charge": 50,
    "min_mass": 500,
    "max_mass": 50000,
    "tol": [
      10,
      10,
      10
    ]
  },
  "flashtnt": {
    "min_length": 3,
    "max_length": 8
  },
  "faims": {
    "cv_values": [],
    "max_cv_skip": 0
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": false,
      "value_ms": 30000
    },
    "agc_interval_seconds": 9999999
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "rank_by": "qscore",
    "max_precursors": 1,
    "additional_scans": [
      "secondary"
    ]
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
  },
  "ms_settings": {
    "ms1": {
      "analyzer": "Orbitrap",
      "first_mass": 500,
      "last_mass": 2000,
      "resolution": 120000,
      "agc_target": 800000,
      "max_it": 246,
      "microscans": 4,
      "data_type": "Centroid",
      "rf_lens": 60,
      "source_cid": 15,
      "source_cid_scaling": 0
    },
    "ms2": {
      "analyzer": "Orbitrap",
      "activation": "HCD",
      "collision_energy": 29,
      "resolution": 120000,
      "agc_target": 500000,
      "max_it": 150
    },
    "additional_ms2": {
      "secondary": {
        "analyzer": "Orbitrap",
        "activation": "HCD",
        "collision_energy": 35,
        "resolution": 60000,
        "agc_target": 300000,
        "max_it": 100
      }
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false
  }
}
)";
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

// P3-U07: empty queue returns an idle survey MS1 (never returns 0)
//
// The drained queue emits a survey, NOT an AGC prescan (ADR-0031). Prescans are scheduled by
// scheduling.agc_interval_seconds alone, and this fixture pins that at 9999999 so the timer cannot
// fire here. priority == 3 is asserted rather than incidental: three C# drain loops break on
// exactly `msn_level == 1 && priority == 3`, so it is the contract, not a detail.
START_SECTION(empty_queue_returns_idle_survey)
{
  FLASHIda* ida = createTestInstance();
  ScanCommand cmd{};
  int result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)  // Queue empty, idle cycle returns a survey MS1
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.is_agc, 0)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_EQUAL(cmd.priority, 3)
  std::string desc(cmd.scan_description);
  TEST_EQUAL(desc.size() >= 4, true)
  TEST_EQUAL(desc[3], 'S')
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
  FLASHIdaTestAccess::push(*ida,cmd3);

  ScanCommand cmd1{};
  cmd1.msn_level = 2;
  cmd1.priority = 1;
  cmd1.scan_id = 101;
  FLASHIdaTestAccess::push(*ida,cmd1);

  ScanCommand cmd2{};
  cmd2.msn_level = 2;
  cmd2.priority = 2;
  cmd2.scan_id = 102;
  FLASHIdaTestAccess::push(*ida,cmd2);

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

  // Queue empty — idle cycle returns a fresh survey MS1 (never returns 0).
  // Note the p3 command dequeued just above was TEST-INJECTED via FLASHIdaTestAccess::push, which
  // bypasses the builders; scan_id 100 distinguishes it from the survey the engine mints here.
  // Nothing in production emits a priority-3 command except Step 5 —
  // only_the_idle_survey_is_emitted_at_priority_3 pins that on engine-emitted commands.
  int idle_result = ida->getNextScanCommand(out);
  TEST_EQUAL(idle_result, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.is_agc, 0)
  TEST_EQUAL(out.msn_level, 1)
  TEST_EQUAL(out.priority, 3)
  TEST_NOT_EQUAL(out.scan_id, 100)

  delete ida;
}
END_SECTION

// P3-U09: AGC scan is dequeued first
START_SECTION(agc_scan_is_dequeued_first)
{
  // Create FLASHIda with agc_interval_seconds = 0 (always needs AGC)
  const char* agc_config = R"({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.9,
    "min_charge": 4,
    "max_charge": 50,
    "min_mass": 500,
    "max_mass": 50000,
    "tol": [
      10,
      10,
      10
    ]
  },
  "flashtnt": {
    "min_length": 3,
    "max_length": 8
  },
  "faims": {
    "cv_values": [],
    "max_cv_skip": 0
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": false,
      "value_ms": 30000
    },
    "agc_interval_seconds": 0
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 1
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
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
    "ms2": {
      "analyzer": "Orbitrap",
      "activation": "ETD",
      "collision_energy": 0,
      "reaction_time": 10.0,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)";

  FLASHIda* ida = new FLASHIda(const_cast<char*>(agc_config));

  // Push an MS2 command into the queue
  ScanCommand ms2_cmd{};
  ms2_cmd.msn_level = 2;
  ms2_cmd.priority = 2;
  ms2_cmd.scan_id = 200;
  FLASHIdaTestAccess::push(*ida,ms2_cmd);

  // First dequeue: with agc_interval_seconds=0, an AGC prescan is always due.
  // It should be returned before any queued command — Step 1 sits ahead of the dequeue.
  //
  // needsAGC() is `elapsed > agc_interval_ms`, so with the interval at 0 any non-zero elapsed
  // qualifies. Sleep past the millisecond boundary the steady_clock difference is measured in:
  // this used to rely on FLASHIda construction happening to take >=1 ms, and the drained-queue
  // path fabricating a prescan anyway if it did not. Both crutches are gone (ADR-0031).
  std::this_thread::sleep_for(std::chrono::milliseconds(5));

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
  TEST_EQUAL(result, 1)  // Queue empty, idle cycle returns a survey MS1
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.is_agc, 0)
  TEST_EQUAL(cmd.priority, 3)
  delete ida;
}
END_SECTION

/////////////////////////////////////////////////////////////
// What the builders put into a command: config sourcing, not queueing.
/////////////////////////////////////////////////////////////

// The AGC prescan is split by physics (ADR-0011): what decides WHICH ions arrive comes from the
// survey's config, what decides HOW they are measured is fixed to the fast-prescan identity.
//
// makeAGC used to read only first_mass/last_mass, leaving the source region at 0. That was
// invisible while ScanFactory guarded those fields with `> 0` -- the keys were simply omitted. Now
// that the group is emitted unconditionally, a zero here would actively command RF lens 0 on the
// scan whose ion-flux estimate gains every scan that follows it.
START_SECTION(makeAGC_takes_source_region_from_config_but_not_speed)
{
  Config cfg{std::string(scan_sourcing_config)};
  ScanCommandQueue queue(cfg);
  ScanCommand cmd = queue.makeAGC();

  // From the survey's config: the prescan must sample the same ion population it is calibrating.
  TEST_REAL_SIMILAR(cmd.rf_lens, 60.0)
  TEST_REAL_SIMILAR(cmd.source_cid, 15.0)
  TEST_REAL_SIMILAR(cmd.source_cid_scaling, 0.0)

  // Fixed to the fast-prescan identity, NOT copied from ms1. microscans is the one that bites:
  // ms1 asks for 4, and an AGC scan at 4 microscans is four times as long -- on a priority-0 scan
  // with a 1 ms max IT. Asserting 1 against a config that says 4 is what makes this non-vacuous.
  TEST_EQUAL(cmd.microscans, 1)
  TEST_STRING_EQUAL(std::string(cmd.data_type), "Profile")
  TEST_STRING_EQUAL(std::string(cmd.analyzer), "IonTrap")
  TEST_STRING_EQUAL(std::string(cmd.scan_rate), "Turbo")
  TEST_EQUAL(cmd.agc_target, 30000)
  TEST_REAL_SIMILAR(cmd.max_it, 1.0)
  TEST_EQUAL(cmd.is_agc, 1)
}
END_SECTION

// buildMS2 must use the ScanConfig it was handed (ADR-0009).
//
// It took a ScanConfig& and then ignored it for exactly two fields, reading level(2).scans[0]
// directly. FLASHIda loops over every level-2 scan config, so ms_settings.ms2[1..N] acquired at
// ms2[0]'s AGC target and max IT, and an exploration override of either key was inert at level 2.
//
// The section is written against scans[1] with values that differ from scans[0]; under the old
// code it reports 500000/150 instead of 300000/100. Run against scans[0] it would pass either way.
START_SECTION(buildMS2_uses_the_scan_config_it_was_given)
{
  Config cfg{std::string(scan_sourcing_config)};
  ScanCommandQueue queue(cfg);

  TEST_EQUAL(cfg.level(2).scans.size(), 2)   // guards the premise: two DIFFERENT configs exist

  PeakGroup pg(10, 10, true);
  ScanCommand cmd = queue.buildMS2(pg, 10, cfg.level(2).scans[1], 2, 0);

  TEST_EQUAL(cmd.agc_target, 300000)          // ISSUE(pre-fix): 500000, i.e. scans[0]'s
  TEST_REAL_SIMILAR(cmd.max_it, 100.0)        // ISSUE(pre-fix): 150.0, i.e. scans[0]'s
  TEST_EQUAL(cmd.orbitrap_resolution, 60000)  // already correct; pins that the fix did not regress it
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
