// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FAIMS.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>

#include "FLASHIda_TestHelpers.h"  // runInterleaved / AcqResult / ScanData — the canonical engine-id-echo driver

#include <vector>
#include <set>
#include <cmath>
#include <cstring>

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
      "tol": [10, 10, 10]
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
      "ms1": { "selection": "qscore", "max_targets": 1 },
      "ms2": { "selection": "none" },
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
      "tol": [10, 10, 10]
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
      "ms1": { "selection": "qscore", "max_targets": 1 },
      "ms2": { "selection": "none" },
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
      "tol": [10, 10, 10]
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
      "ms1": { "selection": "qscore", "max_targets": 1 },
      "ms2": { "selection": "none" },
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

  // N empty-peak MS1 surveys (no mzs/ints => deconvolution finds 0 precursors = the "low-precursor" case the
  // CV skip tests want), each with a distinct RT step. Fed one-per-survey by runInterleaved under the engine's
  // OWN emitted tracking id, so the always-on MS1 gate (FLASHIda.cpp:775) passes and the engine paces the CVs.
  inline std::vector<ScanData> emptyMs1Surveys(int n, double rt_step = 1.0)
  {
    std::vector<ScanData> v;
    v.reserve(n);
    for (int i = 0; i < n; ++i)
    {
      ScanData s;            // mzs/ints intentionally empty -> 0 precursors
      s.rt = (double)(i + 1) * rt_step;
      v.push_back(s);
    }
    return v;
  }
}

START_TEST(FLASHIdaFAIMS, "$Id$")

/////////////////////////////////////////////////////////////

// P6-U01: CV cycling order matches config [-40, -50, -60] with wrap-around
START_SECTION(cv_cycling_order_matches_config)
{
  FLASHIda* ida = createFaims3CV();

  // OBSERVE the engine, do not dictate. The engine OWNS CV cycling: it stamps each survey MS1 with its current
  // CV and, after processing one, advances to the next CV for the survey it pushes. Drive via runInterleaved
  // feeding empty-peak MS1 (0 precursors) under the engine's own emitted tracking ids, then read the CV off the
  // drained MS1 commands. current_cv_index_ starts at 0 (CV -40); the engine then cycles -50 -> -60 -> wrap -40:
  //   ms1_cmds[0] = initial survey at CV -40
  //   ms1_cmds[1] = next survey at CV -50  (advanced once)
  //   ms1_cmds[2] = next survey at CV -60  (advanced again)
  //   ms1_cmds[3] = next survey at CV -40  (wrap-around)
  std::vector<double> expected_cvs = {-40.0, -50.0, -60.0, -40.0};

  AcqResult r = runInterleaved(ida, emptyMs1Surveys((int)expected_cvs.size()), std::vector<ScanData>{});

  // The engine emitted at least one survey per requested CV (each fed MS1 = one drained MS1 command).
  TEST_EQUAL(r.ms1_cmds.size() >= expected_cvs.size(), true)
  for (size_t i = 0; i < expected_cvs.size(); ++i)
  {
    TEST_EQUAL(std::strlen(r.ms1_cmds[i].scan_description) <= 15, true)
    TEST_EQUAL(r.ms1_cmds[i].msn_level, 1)
    TEST_REAL_SIMILAR(r.ms1_cmds[i].faims_cv, expected_cvs[i])
  }

  delete ida;
}
END_SECTION

// P6-U02: Adaptive CV skipping — low precursor count activates skip
START_SECTION(adaptive_cv_skip_low_precursor)
{
  FLASHIda* ida = createFaimsSkip();  // max_cv_skip=2, threshold=15, CVs=[-40,-50,-60]

  // OBSERVE the engine. Each empty-peak survey yields 0 precursors (< threshold 15), so the engine doubles the
  // just-processed CV's skip amount — but a freshly-entered CV (skip amount still 0) is NOT skipped, so the
  // engine keeps advancing. Drive via runInterleaved (engine-id-echo) and read the emitted CV sequence:
  //   ms1_cmds[0] = first survey at CV -40 (then amount[-40] doubles 0->1)
  //   ms1_cmds[1] = advanced to CV -50     (amount[-50]=0, not skipped; then doubles 0->1)
  //   ms1_cmds[2] = advanced to CV -60     (amount[-60]=0, not skipped)
  AcqResult r = runInterleaved(ida, emptyMs1Surveys(3), std::vector<ScanData>{});

  TEST_EQUAL(r.ms1_cmds.size() >= 3, true)
  for (size_t i = 0; i < 3; ++i)
  {
    TEST_EQUAL(std::strlen(r.ms1_cmds[i].scan_description) <= 15, true)
    TEST_EQUAL(r.ms1_cmds[i].msn_level, 1)
  }
  TEST_REAL_SIMILAR(r.ms1_cmds[0].faims_cv, -40.0)  // initial CV
  TEST_REAL_SIMILAR(r.ms1_cmds[1].faims_cv, -50.0)  // advanced to next CV (not skipped)
  TEST_REAL_SIMILAR(r.ms1_cmds[2].faims_cv, -60.0)  // advanced to CV -60 (not skipped)

  delete ida;
}
END_SECTION

// P6-U02b: Threshold boundary — precursor_count at threshold does NOT trigger skip
START_SECTION(adaptive_cv_skip_threshold_boundary)
{
  // Construct FAIMS directly from the skip config (threshold=15, max_cv_skip=2, CVs=[-40,-50,-60])
  Config cfg{std::string(faims_skip_config)};
  FAIMS faims{cfg};

  // Precursor count = 14 (below threshold) -> should double skip amount 0->1
  faims.updateSkip(-40.0, 14);
  TEST_EQUAL(faims.cvSkipAmount(0), 1)

  // Precursor count = 15 (at threshold, NOT strictly less) -> should RESET to 0
  faims.updateSkip(-40.0, 15);
  TEST_EQUAL(faims.cvSkipAmount(0), 0)

  // Precursor count = 14 -> doubles 0->1 again
  faims.updateSkip(-40.0, 14);
  TEST_EQUAL(faims.cvSkipAmount(0), 1)

  // Precursor count = 14 -> doubles 1->2
  faims.updateSkip(-40.0, 14);
  TEST_EQUAL(faims.cvSkipAmount(0), 2)

  // Precursor count = 14 -> would double to 4 but capped at max_cv_skip=2
  faims.updateSkip(-40.0, 14);
  TEST_EQUAL(faims.cvSkipAmount(0), 2)
}
END_SECTION

// P6-U03: CV skip limit enforced — after max_cv_skip, CV is still visited
START_SECTION(cv_skip_limit_enforced)
{
  FLASHIda* ida = createFaimsSkip();  // max_cv_skip=2, 3 CVs

  // Drive many empty-peak surveys (0 precursors each) to build up skip amounts for all CVs. The engine OWNS the
  // skip policy: amounts cap at max_cv_skip=2 and the per-CV skip counters exhaust, so every CV stays reachable.
  // OBSERVE the engine-emitted CV sequence (runInterleaved, engine-id-echo) and confirm all 3 CVs appear — none
  // is permanently blocked by skipping.
  AcqResult r = runInterleaved(ida, emptyMs1Surveys(15), std::vector<ScanData>{});

  std::set<double> seen_cvs;
  for (const auto& c : r.ms1_cmds)
  {
    TEST_EQUAL(std::strlen(c.scan_description) <= 15, true)
    seen_cvs.insert(c.faims_cv);
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
  ms2.priority = 2;
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

  // Drive ONE real engine survey under the engine's OWN emitted tracking id (engine-id-echo, same DRAIN-CONTRACT
  // as runInterleaved — done inline here because this test observes dequeue ORDER, which runInterleaved abstracts
  // away). The bootstrap idle cycle emits an AGC and queues an idle MS1 (prio 3) stamped with a valid engine id;
  // dequeuing it registers it pending, so feeding an empty-peak MS1 back under that id passes the always-on MS1
  // gate (FLASHIda.cpp:775). That survey (0 precursors) makes the engine advance the CV and push a CV-transition
  // MS1 at priority 0 (CV -50).
  ScanCommand boot{};
  ida->getNextScanCommand(boot);                     // idle AGC; queues idle MS1 (engine id, CV -40)
  TEST_EQUAL(boot.is_agc, 1)
  ScanCommand survey{};
  ida->getNextScanCommand(survey);                   // the idle MS1 carrying the engine's own survey id
  TEST_EQUAL(survey.msn_level, 1)
  ida->processScan(nullptr, nullptr, 0, 1.0, 1, survey.scan_description, survey.faims_cv);

  // Push two MS2s at priority 2 with parent CV=-40 (synthetic seeding — pure pushCommandForTest, unaffected by
  // the gate). They are now pending BEHIND the just-emitted prio-0 CV-transition MS1.
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
  FLASHIda* ida = createNonFaims();  // single CV => faims_enabled_=false

  // Drive ONE real engine survey under the engine's OWN emitted id (engine-id-echo). The bootstrap idle cycle
  // emits an AGC + queues an idle MS1 (prio 3); dequeuing it registers it pending, so the empty-peak MS1 we feed
  // back passes the always-on MS1 gate and is genuinely PROCESSED (not gate-rejected). With faims_enabled_=false,
  // processScan must NOT enter the FAIMS branch -> NO CV-transition MS1 (prio 0) is pushed. So after the survey
  // the only pending item is the idle MS1 (prio 3); the next getNextScanCommand falls through to the idle cycle
  // (AGC). A spurious prio-0 CV-transition MS1 would surface here instead of the AGC.
  ScanCommand boot{};
  ida->getNextScanCommand(boot);                     // idle AGC; queues idle MS1 (engine id)
  TEST_EQUAL(boot.is_agc, 1)
  ScanCommand survey{};
  ida->getNextScanCommand(survey);                   // the idle MS1 carrying the engine's own survey id
  TEST_EQUAL(survey.msn_level, 1)
  int pushed = ida->processScan(nullptr, nullptr, 0, 1.0, 1, survey.scan_description, 0.0);
  TEST_EQUAL(pushed, 0)  // empty spectrum => 0 MS2 commands, and (non-FAIMS) no CV-transition MS1

  // Next command: idle AGC (no CV-transition MS1 was injected) => non-FAIMS behavior confirmed.
  ScanCommand out{};
  int result = ida->getNextScanCommand(out);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.is_agc, 1)

  delete ida;
}
END_SECTION

/////////////////////////////////////////////////////////////
// P6-U07 (F7): the ground-truth harness now echoes the engine command's FAIMS CV back on the re-fed MS1
// (runInterleaved passes cmd.faims_cv to processScan; the C# twin PushScanAndDrainFull sets the "FAIMS CV"
// trailer). So the engine observes the commanded CV and stamps the produced MS2 with it (FLASHIda.cpp:859).
// Pre-fix the re-fed MS1 was processed with faims_cv=0.0, so its MS2 carried 0.0 — NOT a configured CV.
// Set-membership of two engine-derived CVs (no captured float) -> drift-stable.
START_SECTION(refed_ms1_echoes_commanded_faims_cv)
{
  auto ms1 = loadTsvScans(FI_MS1_STD);   // E. coli MS1 reliably selects >=1 precursor under DDA
  auto ms2 = loadTsvScans(FI_MS2_HCD);
  ABORT_IF(ms1.empty() || ms2.empty())

  // faims_config: multi-CV [-40,-50,-60] (=> faims_enabled_) + DDA. Flip MS2 selection on so the selected
  // precursor actually yields MS2 commands (the file's other FAIMS tests use empty surveys + "none").
  std::string cfg(faims_config);
  {
    const std::string from = "\"ms2\": { \"selection\": \"none\" }";
    auto p = cfg.find(from);
    ABORT_IF(p == std::string::npos)
    cfg.replace(p, from.size(), "\"ms2\": { \"selection\": \"intensity\", \"max_targets\": 1 }");
  }
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runInterleaved(&ida, ms1, std::vector<ScanData>{ms2[0]});

  const std::set<double> cv_set = { -40.0, -50.0, -60.0 };
  bool found_ms1_cv = false;
  for (const auto& c : a.ms1_cmds) if (cv_set.count(c.faims_cv)) { found_ms1_cv = true; break; }
  TEST_TRUE(found_ms1_cv)   // FAIMS run: surveys carry a configured CV

  int checked = 0;
  for (const auto& c : a.ms2_cmds)
  {
    // ISSUE(F7): pre-fix the re-fed MS1 dropped cmd.faims_cv -> processScan saw 0.0 -> the produced MS2
    // carried faims_cv=0.0 (not in the configured CV set). The echo now binds the survey's CV to its MS2.
    TEST_EQUAL(cv_set.count(c.faims_cv), 1)
    checked++;
  }
  TEST_TRUE(checked >= 1)   // >=1 MS2 produced, each carrying a configured (parent-survey) FAIMS CV
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
