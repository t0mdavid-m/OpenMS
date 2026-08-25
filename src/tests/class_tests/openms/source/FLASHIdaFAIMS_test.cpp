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
#include "FLASHIda_TestAccess.h"   // FLASHIdaTestAccess::push (private-state access)

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
    "cv_values": [
      -40,
      -50,
      -60
    ],
    "max_cv_skip": 0,
    "cv_precursor_threshold": 15
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
      "activation": "HCD",
      "collision_energy": 29,
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

  // JSON config with 3 CVs, adaptive skipping enabled
  const char* faims_skip_config = R"({
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
    "cv_values": [
      -40,
      -50,
      -60
    ],
    "max_cv_skip": 2,
    "cv_precursor_threshold": 15
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
      "activation": "HCD",
      "collision_energy": 29,
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

  // JSON config with no FAIMS. An EMPTY cv_values is the only way to say that (ADR-0012).
  // This fixture used to carry [-50] and rely on faims_.enabled being `cv_values.size() > 1`,
  // i.e. it expressed "no FAIMS" as "one CV, which is not enough to cycle". That conflation was
  // the defect: a real fixed-CV FAIMS method was silently treated as no FAIMS at all.
  const char* non_faims_config = R"({
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
      "activation": "HCD",
      "collision_energy": 29,
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

  // non_faims_config with a different cv_values, so enablement is the ONLY variable. Building the
  // string this way rather than hand-writing a config keeps every other key identical by
  // construction -- an enablement test that also perturbed selection or tolerances would prove
  // nothing about enablement.
  std::string withCvValues(const std::string& cv_list)
  {
    std::string s(non_faims_config);
    const std::string needle = "\"cv_values\": []";
    const size_t p = s.find(needle);
    if (p == std::string::npos) return s;   // fixture drifted; the section's TEST_EQUALs will catch it
    return s.replace(p, needle.size(), "\"cv_values\": " + cv_list);
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
  FLASHIdaTestAccess::push(*ida,ms2);

  // Dequeue the MS2 — it should retain its parent CV
  ScanCommand out{};
  int result = ida->getNextScanCommand(out);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.msn_level, 2)
  TEST_REAL_SIMILAR(out.faims_cv, -40.0)  // parent CV preserved

  // Queue is now empty — idle cycle returns a survey MS1, not an AGC prescan (ADR-0031).
  result = ida->getNextScanCommand(out);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.is_agc, 0)
  TEST_EQUAL(out.msn_level, 1)
  TEST_EQUAL(out.priority, 3)

  delete ida;
}
END_SECTION

// P6-U05: CV transition injects an MS1 with the next CV at priority 0 (so it is sent before pending MS2s)
START_SECTION(cv_transition_ms1_before_ms2s)
{
  // Drive a REAL, guaranteed-selecting survey through the canonical interleaved harness (engine-id-echo) and
  // observe the CROSS-LEVEL dequeue order over AcqResult.all_cmds (the additive raw-order capability — the
  // per-level ms1/ms2/ms3 buckets lose cross-level interleave). ms1_ecoli_rich[0] deterministically selects
  // >=1 precursor at the -40 STARTUP CV (cf. ProcessScan A1[0]==1 under a value-identical selection config),
  // so processScan pushes prio-2 MS2 commands at CV -40 (FLASHIda.cpp:858) and THEN a prio-0 CV-transition MS1
  // carrying the NEXT CV -50 (FLASHIda.cpp:898-915). The queue dequeues by priority 0->3, so the prio-0
  // CV-transition MS1 is drained BEFORE the prio-2 MS2s — deterministically. Feeding exactly ONE survey
  // (n_ms1==1) keeps a single selection event: the transition MS1 then arrives as an idle tick (recorded in
  // all_cmds, not re-fed), so no second selection perturbs the trace -> exact -50/-40 are pinnable.
  const std::string FI_MS1_ECOLI = "../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt";
  auto ecoli = loadTsvScans(FI_MS1_ECOLI);
  ABORT_IF(ecoli.empty())

  FLASHIda ida(const_cast<char*>(faims_3cv_config));
  AcqResult a = runInterleaved(&ida, std::vector<ScanData>{ecoli[0]}, std::vector<ScanData>{});

  // First prio-0 NON-AGC MS1 (the CV transition) and first MS2, in raw dequeue order. Exclude the startup idle
  // The CV-transition MS1 is the ONLY is_agc==0, prio-0, level-1 command the engine emits; the idle
  // survey MS1 is priority 3. A scheduled AGC prescan would also be msn_level 1 / priority 0 but
  // carries is_agc==1 — and this fixture pins agc_interval_seconds at 9999999, so none can fire.
  int cv_idx = -1, ms2_idx = -1;
  for (int i = 0; i < (int)a.all_cmds.size(); ++i)
  {
    if (cv_idx < 0 && a.all_cmds[i].msn_level == 1 && a.all_cmds[i].priority == 0 && a.all_cmds[i].is_agc == 0) cv_idx = i;
    if (ms2_idx < 0 && a.all_cmds[i].msn_level == 2) ms2_idx = i;
  }

  TEST_TRUE(cv_idx >= 0)                 // a prio-0 CV-transition MS1 was emitted
  TEST_TRUE(ms2_idx >= 0)                // the survey selected >=1 precursor -> a prio-2 MS2 was emitted
  TEST_TRUE(cv_idx < ms2_idx)            // prio-0 CV-transition MS1 drained BEFORE the prio-2 MS2s
  TEST_REAL_SIMILAR(a.all_cmds[cv_idx].faims_cv, -50.0)   // CV-transition MS1 carries the NEXT CV (-40 -> -50)
  TEST_REAL_SIMILAR(a.all_cmds[ms2_idx].faims_cv, -40.0)  // MS2 carries the -40 startup (parent survey) CV
}
END_SECTION

// Enablement and cycling are two different questions (ADR-0012).
//
// faims_.enabled used to be `cv_values.size() > 1`, which answered "is there anything to cycle
// between?" and then used that as the answer to "is FAIMS in use?". A fixed-CV method -- one CV,
// perfectly ordinary -- therefore reported no FAIMS: faims_cv stayed 0, ScanFactory's magnitude
// test failed, and the run silently acquired at whatever FAIMS state the instrument method held.
START_SECTION(faims_enablement_and_cycling_are_separate)
{
  {
    Config c(withCvValues("[]"));
    TEST_EQUAL(c.faims().enabled, false)          // empty is the ONLY way to say "no FAIMS"
    FAIMS f(c);
    TEST_EQUAL(f.isEnabled(), false)
    TEST_EQUAL(f.isCycling(), false)
  }
  {
    Config c(withCvValues("[-45]"));
    FAIMS f(c);
    TEST_EQUAL(f.isEnabled(), true)               // ISSUE(pre-fix): was false -- the whole defect
    TEST_EQUAL(f.isCycling(), false)              // one CV: nothing to rotate between
    TEST_REAL_SIMILAR(f.currentCV(), -45.0)
  }
  {
    // CV 0 is a legitimate compensation voltage, and is now distinguishable from "no FAIMS".
    Config c(withCvValues("[0]"));
    FAIMS f(c);
    TEST_EQUAL(f.isEnabled(), true)
    TEST_EQUAL(f.isCycling(), false)
    TEST_REAL_SIMILAR(f.currentCV(), 0.0)
  }
  {
    Config c(withCvValues("[-40, -50]"));
    FAIMS f(c);
    TEST_EQUAL(f.isEnabled(), true)
    TEST_EQUAL(f.isCycling(), true)
  }
}
END_SECTION

// A fixed-CV run carries its CV but must not transition.
//
// The CV-transition MS1 push is guarded on isCycling(), not isEnabled(). Guarding it on
// isEnabled() would make a single-CV run push a priority-0 MS1 after EVERY MS1 -- advanceToNextCV
// would keep returning the one CV it has -- silently doubling the survey rate. This is the
// regression the split exists to prevent, so it is asserted rather than assumed.
START_SECTION(single_cv_is_enabled_but_never_transitions)
{
  std::string cfg = withCvValues("[-45]");
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult r = runInterleaved(&ida, emptyMs1Surveys(3), std::vector<ScanData>{});

  TEST_EQUAL(r.ms1_cmds.empty(), false)
  for (const auto& c : r.ms1_cmds)
  {
    TEST_EQUAL(c.msn_level, 1)
    TEST_NOT_EQUAL(c.priority, 0)              // no CV-transition MS1 was injected
    TEST_REAL_SIMILAR(c.faims_cv, -45.0)       // ...but the configured CV still travels
    TEST_EQUAL(c.faims_enabled, 1)             // ...and the instrument is told FAIMS is on
  }
}
END_SECTION

// P6-U06: Non-FAIMS mode — processScan does not push a CV-transition MS1
START_SECTION(non_faims_no_cv_transition)
{
  FLASHIda* ida = createNonFaims();  // empty cv_values => faims_.enabled false

  // Drive via the canonical interleaved harness (engine-id-echo). Feed empty-peak surveys (0 precursors) the same
  // way the FAIMS sections do; the difference under test is purely the engine response. With faims_enabled_=false,
  // processScan must NOT enter the FAIMS branch (FLASHIda.cpp:899) -> NO priority-0 CV-transition MS1 is ever
  // pushed. "No CV-transition MS1" is a per-level property of ms1_cmds (not a cross-level interleave), so it is
  // faithfully re-expressible over AcqResult.
  AcqResult r = runInterleaved(ida, emptyMs1Surveys(3), std::vector<ScanData>{});

  TEST_EQUAL(r.ms1_cmds.empty(), false)  // the engine still emitted idle survey MS1s

  // Every emitted MS1 is an ordinary (idle, prio 3) survey — none is a prio-0 CV-transition MS1. With
  // faims_enabled_=false the engine reports current_faims_cv_ as 0.0 (FLASHIda.cpp:919,1285), so the idle
  // survey MS1s carry faims_cv 0.0 (the FAIMS branch never stamps a configured CV). A spurious FAIMS branch
  // would surface as a prio-0 MS1 carrying a non-zero (configured) CV.
  for (const auto& c : r.ms1_cmds)
  {
    TEST_EQUAL(std::strlen(c.scan_description) <= 15, true)
    TEST_EQUAL(c.msn_level, 1)
    TEST_NOT_EQUAL(c.priority, 0)            // non-FAIMS: no priority-0 CV-transition MS1 injected
    TEST_REAL_SIMILAR(c.faims_cv, 0.0)       // non-FAIMS reports CV 0.0; no CV transition occurs
  }

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

  // faims_3cv_config AS-IS: multi-CV [-40,-50,-60] (=> faims_enabled_) + DDA. MS1->MS2 selection is already on
  // (selection_strategy.ms1 = qscore); feeding a REAL survey yields MS2 commands carrying the parent survey's
  // FAIMS CV. Do NOT enable ms2 (MS2->MS3) selection: it defaults to fragment matching, which Config::validate()
  // requires a protein_sequence for — and this DDA config has none (would throw at construction).
  FLASHIda ida(const_cast<char*>(faims_3cv_config));
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
