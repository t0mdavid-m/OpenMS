// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include <vector>
#include <set>
#include <cmath>

using namespace OpenMS;

namespace
{
  // JSON config with 3 CVs, no skipping
  const char* faims_3cv_config = R"({
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
      "cv_values": [-40, -50, -60],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
          "activation": "HCD",
          "collision_energy": 29,
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
      "ms1": { "selection": "qscore", "max_precursors": 1 },
      "ms2": { "selection": "intensity" },
      "ms3": { "selection": "none" }
    }
  })";

  // JSON config with 3 CVs, adaptive skipping enabled
  const char* faims_skip_config = R"({
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
      "cv_values": [-40, -50, -60],
      "max_cv_skip": 2,
      "cv_precursor_threshold": 15
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
          "activation": "HCD",
          "collision_energy": 29,
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
      "ms1": { "selection": "qscore", "max_precursors": 1 },
      "ms2": { "selection": "intensity" },
      "ms3": { "selection": "none" }
    }
  })";

  // JSON config with no FAIMS (single CV = non-FAIMS mode)
  const char* non_faims_config = R"({
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
          "activation": "HCD",
          "collision_energy": 29,
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
      "ms1": { "selection": "qscore", "max_precursors": 1 },
      "ms2": { "selection": "intensity" },
      "ms3": { "selection": "none" }
    }
  })";

  FLASHIda* createFaims3CV()
  {
    return new FLASHIda(const_cast<char*>(faims_3cv_config));
  }

  FLASHIda* createFaimsSkip()
  {
    return new FLASHIda(const_cast<char*>(faims_skip_config));
  }

  FLASHIda* createNonFaims()
  {
    return new FLASHIda(const_cast<char*>(non_faims_config));
  }
}

START_TEST(FLASHIdaFAIMS, "$Id$")

/////////////////////////////////////////////////////////////

// P6-U01: CV cycling order matches config [-40, -50, -60] with wrap-around
START_SECTION(cv_cycling_order_matches_config)
{
  FLASHIda* ida = createFaims3CV();

  // Each processScan at MS1 level pushes a CV-transition MS1 into the queue.
  // current_cv_index_ starts at 0; advanceToNextCV_ increments first:
  //   1st processScan -> advance to index 1 -> CV-transition MS1 with CV -50
  //   2nd processScan -> advance to index 2 -> CV-transition MS1 with CV -60
  //   3rd processScan -> advance to index 0 (wrap) -> CV-transition MS1 with CV -40
  //   4th processScan -> advance to index 1 -> CV-transition MS1 with CV -50

  // Input: what the instrument scanned at (starts at initial CV -40, then follows transitions)
  std::vector<double> input_cvs = {-40.0, -50.0, -60.0, -40.0};
  // Output: what advanceToNextCV_ pushes as CV-transition MS1
  std::vector<double> expected_cvs = {-50.0, -60.0, -40.0, -50.0};
  for (size_t i = 0; i < expected_cvs.size(); ++i)
  {
    ida->processScan(nullptr, nullptr, 0, (double)(i + 1), 1, "ms1", input_cvs[i]);
    ScanCommand cmd{};
    int result = ida->getNextScanCommand(cmd);
    TEST_EQUAL(result, 1)
    TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
    TEST_EQUAL(cmd.msn_level, 1)
    TEST_REAL_SIMILAR(cmd.faims_cv, expected_cvs[i])
    // Drain queued commands; stop at idle AGC + consume following idle MS1
    while (ida->getNextScanCommand(cmd) == 1) {
      TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
      if (cmd.is_agc) { ida->getNextScanCommand(cmd); TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true) break; }
    }
  }

  delete ida;
}
END_SECTION

// P6-U02: Adaptive CV skipping — low precursor count activates skip
START_SECTION(adaptive_cv_skip_low_precursor)
{
  FLASHIda* ida = createFaimsSkip();

  // processScan with empty spectrum -> 0 precursors < threshold 15
  // -> updateCVSkip_ doubles CVSkipAmount[0] from 0 to 1
  // -> advanceToNextCV_: index 1 -> CV -50 (CVSkipAmount[1]=0, no skip)
  ida->processScan(nullptr, nullptr, 0, 1.0, 1, "ms1", -40.0);

  ScanCommand cmd{};
  int result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_REAL_SIMILAR(cmd.faims_cv, -50.0)  // advanced to next CV

  // Drain remaining; consume idle AGC + following idle MS1
  while (ida->getNextScanCommand(cmd) == 1) {
    TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
    if (cmd.is_agc) { ida->getNextScanCommand(cmd); TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true) break; }
  }

  // Second empty scan at CV=-50 -> CVSkipAmount[1] doubles 0->1
  // advanceToNextCV_ from index 1: index 2 -> CV -60 (amount=0, use it)
  ida->processScan(nullptr, nullptr, 0, 2.0, 1, "ms1", -50.0);

  result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_REAL_SIMILAR(cmd.faims_cv, -60.0)  // advanced to CV -60

  delete ida;
}
END_SECTION

// P6-U02b: Threshold boundary — precursor_count at threshold does NOT trigger skip
START_SECTION(adaptive_cv_skip_threshold_boundary)
{
  FLASHIda* ida = createFaimsSkip();  // threshold=15, max_cv_skip=2

  // Precursor count = 14 (below threshold) -> should double skip amount 0->1
  ida->updateCVSkipForTest(-40.0, 14);
  TEST_EQUAL(ida->getCVSkipAmountForTest(0), 1)

  // Precursor count = 15 (at threshold, NOT strictly less) -> should RESET to 0
  ida->updateCVSkipForTest(-40.0, 15);
  TEST_EQUAL(ida->getCVSkipAmountForTest(0), 0)

  // Precursor count = 14 -> doubles 0->1 again
  ida->updateCVSkipForTest(-40.0, 14);
  TEST_EQUAL(ida->getCVSkipAmountForTest(0), 1)

  // Precursor count = 14 -> doubles 1->2
  ida->updateCVSkipForTest(-40.0, 14);
  TEST_EQUAL(ida->getCVSkipAmountForTest(0), 2)

  // Precursor count = 14 -> would double to 4 but capped at max_cv_skip=2
  ida->updateCVSkipForTest(-40.0, 14);
  TEST_EQUAL(ida->getCVSkipAmountForTest(0), 2)

  delete ida;
}
END_SECTION

// P6-U03: CV skip limit enforced — after max_cv_skip, CV is still visited
START_SECTION(cv_skip_limit_enforced)
{
  FLASHIda* ida = createFaimsSkip();  // max_cv_skip=2, 3 CVs

  // Drive multiple cycles with empty spectra to build up skip amounts for all CVs.
  // After enough cycles, all CVSkipAmounts should be capped at 2.
  // But the skip counters exhaust, so all CVs are still reachable.
  std::set<double> seen_cvs;
  for (int cycle = 0; cycle < 15; ++cycle)
  {
    // Rotate through CVs
    double cv = (cycle % 3 == 0) ? -40.0 : (cycle % 3 == 1) ? -50.0 : -60.0;
    ida->processScan(nullptr, nullptr, 0, (double)(cycle + 1), 1, "ms1", cv);

    // Drain and record CV values
    ScanCommand drain{};
    while (ida->getNextScanCommand(drain) == 1)
    {
      TEST_EQUAL(std::strlen(drain.scan_description) <= 15, true)
      if (drain.is_agc) { ida->getNextScanCommand(drain); TEST_EQUAL(std::strlen(drain.scan_description) <= 15, true) break; } // consume idle AGC + MS1
      seen_cvs.insert(drain.faims_cv);
    }
  }

  // All 3 CVs should appear — none permanently blocked
  TEST_EQUAL(seen_cvs.count(-40.0), 1)
  TEST_EQUAL(seen_cvs.count(-50.0), 1)
  TEST_EQUAL(seen_cvs.count(-60.0), 1)

  delete ida;
}
END_SECTION

// P6-U04: MS2 commands carry parent MS1's CV
START_SECTION(ms2_carries_parent_cv)
{
  FLASHIda* ida = createFaims3CV();

  // Push an MS2 with parent CV=-40 set at build time
  ScanCommand ms2{};
  ms2.msn_level = 2;
  ms2.priority = 1;
  ms2.scan_id = 999;
  ms2.faims_cv = -40.0;
  ida->pushCommandForTest(ms2);

  // Dequeue the MS2 — it should retain its parent CV
  ScanCommand out{};
  int result = ida->getNextScanCommand(out);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.msn_level, 2)
  TEST_REAL_SIMILAR(out.faims_cv, -40.0)  // parent CV preserved

  // Queue is now empty — idle cycle returns AGC
  result = ida->getNextScanCommand(out);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.is_agc, 1)

  delete ida;
}
END_SECTION

// P6-U05: CV transition injects MS1 with new CV; MS2s carry parent CV
START_SECTION(cv_transition_ms1_before_ms2s)
{
  FLASHIda* ida = createFaims3CV();

  // Push two MS2s at priority 2 with parent CV=-40
  ScanCommand ms2_a{};
  ms2_a.msn_level = 2;
  ms2_a.priority = 2;
  ms2_a.scan_id = 500;
  ms2_a.faims_cv = -40.0;
  ida->pushCommandForTest(ms2_a);

  ScanCommand ms2_b{};
  ms2_b.msn_level = 2;
  ms2_b.priority = 2;
  ms2_b.scan_id = 501;
  ms2_b.faims_cv = -40.0;
  ida->pushCommandForTest(ms2_b);

  // processScan at CV=-40 pushes a CV-transition MS1 at priority 0 for next CV
  ida->processScan(nullptr, nullptr, 0, 1.0, 1, "ms1", -40.0);

  // Dequeue order: CV-transition MS1 (prio 0), then MS2s (prio 2)
  ScanCommand out{};

  ida->getNextScanCommand(out);
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.msn_level, 1)
  TEST_EQUAL(out.priority, 0)
  TEST_REAL_SIMILAR(out.faims_cv, -50.0)  // CV transition to next CV

  ida->getNextScanCommand(out);
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.msn_level, 2)
  TEST_REAL_SIMILAR(out.faims_cv, -40.0)  // parent CV preserved

  ida->getNextScanCommand(out);
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.msn_level, 2)
  TEST_REAL_SIMILAR(out.faims_cv, -40.0)  // parent CV preserved

  delete ida;
}
END_SECTION

// P6-U06: Non-FAIMS mode — processScan does not push CV-transition MS1
START_SECTION(non_faims_no_cv_transition)
{
  FLASHIda* ida = createNonFaims();

  // Call processScan with non-FAIMS config (single CV = faims_enabled_=false).
  // With null spectrum data, produces 0 MS2 commands.
  // Key check: faims_enabled_=false means NO CV-transition MS1 is pushed.
  ida->processScan(nullptr, nullptr, 0, 1.0, 1, "ms1", 0.0);

  // Queue should be empty — no CV-transition MS1, no MS2 commands
  ScanCommand out{};
  int result = ida->getNextScanCommand(out);
  TEST_EQUAL(result, 1)  // idle cycle returns AGC = non-FAIMS behavior confirmed
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.is_agc, 1)

  delete ida;
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
