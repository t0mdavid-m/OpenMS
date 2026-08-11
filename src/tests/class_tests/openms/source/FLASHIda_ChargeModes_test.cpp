// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Behavioural tests for precursor_selection.precursor_charges -- the MS2 acquisition GEOMETRY axis
// (ADR-0016): single | separate | multiplexed.
//
// These assert what each value PROMISES, under a config that sets nothing else unusual. That
// distinction is the whole reason the file exists. "separate" was previously covered by exactly one
// section, and that section configured charge_based_exclusion: true -- an unrelated exclusion-KEYING
// developer flag which, as a side effect, was the only thing that made the mode's charge list
// multi-valued. With the flag at its default the mode emitted one scan, i.e. behaved exactly like
// "single", in every shipped config, and the suite was green throughout.
//
// So: no config here sets charge_based_exclusion. Geometry must come from precursor_charges alone.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include "FLASHIda_TestHelpers.h"  // ground-truth harness: ScanData/loadTsvScans/runInterleaved

#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;

namespace
{
  // cytC (~12356 Da neutral) resolves at 9-10 charge states per species in this fixture, so the
  // difference between "one charge" and "the envelope" is large and unambiguous. ms1_standard would
  // not do: it yields <=1 selectable mass per scan.
  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_cytc.txt";

  // One config, one varying value. max_precursors is 10 so the budget never masks the geometry:
  // the species guard bounds selection by SPECIES, and a whole envelope costs ONE slot (ADR-0016),
  // but a cramped budget would still make "did it fan out?" ambiguous.
  std::string cfgFor(const std::string& precursor_charges)
  {
    return std::string(R"JSON({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.3,
    "min_charge": 4,
    "max_charge": 50,
    "min_mass": 500,
    "max_mass": 50000,
    "tol": [10, 10, 10]
  },
  "flashtnt": { "min_length": 3, "max_length": 8 },
  "faims": { "cv_values": [], "max_cv_skip": 0 },
  "scheduling": {
    "cycle_time": { "enabled": false, "value_ms": 60000 },
    "scan_timeout": { "enabled": true, "value_ms": 30000 },
    "agc_interval_seconds": 30
  },
  "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 10,
    "precursor_charges": ")JSON")
           + precursor_charges + R"JSON("
  },
  "characterization": { "mode": "off", "max_targets": 10 },
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
  "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 }
}
)JSON";
  }

  // Every charge the emitted MS2s isolate, keyed by (survey, mass).
  //
  // Keying per SURVEY, not per run, is the point: charge-keyed fallback legitimately reaches other
  // charges of a mass on LATER surveys, so a run-wide grouping would conflate "acquired the envelope
  // at once" with "drifted onto another charge eventually" -- the exact conflation that let the
  // defect survive. parent_scan_id identifies the spawning survey; the first survey's MS2s carry an
  // EMPTY parent_scan_id (buildMS2 stamps it only for a non-zero parent tracking id), which is still
  // one coherent group.
  //
  // @p with_notches folds a multiplexed command's co-isolation notches in alongside its anchor, so
  // the two geometries become comparable as SETS.
  std::map<std::pair<std::string, long long>, std::set<int>>
  chargesPerSurveyMass(const AcqResult& a, bool with_notches)
  {
    std::map<std::pair<std::string, long long>, std::set<int>> out;
    for (const auto& c : a.ms2_cmds)
    {
      if (c.num_stages < 1) continue;
      auto& charges = out[{std::string(c.parent_scan_id), std::llround(c.mono_mass)}];
      charges.insert(c.stages[0].charge_state);
      if (with_notches)
      {
        auto [first, count] = notchesForStage(c, 0);
        for (int i = 0; i < count; ++i) charges.insert(first[i].charge_state);
      }
    }
    return out;
  }

  AcqResult runMode(FLASHIda& ida, const std::vector<ScanData>& scans)
  {
    return runInterleaved(&ida, scans, {}, nullptr, /*max_iters=*/4000);
  }
}

START_TEST(FLASHIda_ChargeModes, "$Id$")

// CM-01: the assertion that was missing. "separate" means one MS2 per SNR-positive charge state, so
// within a SINGLE survey a mass must appear at two or more charges -- under a config that does
// nothing but ask for it.
//
// Under the defect every group has exactly one charge and this fails, which is the point.
START_SECTION(separate_fans_out_under_a_default_config)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("separate");
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)  // non-vacuity: with no MS2 the assertion is trivial

  const auto groups = chargesPerSurveyMass(a, /*with_notches=*/false);
  int masses_with_multiple_charges = 0;
  for (const auto& kv : groups)
    if (kv.second.size() >= 2) ++masses_with_multiple_charges;

  if (masses_with_multiple_charges == 0)
    std::cout << "[CM-01] no mass acquired at >=2 charges within one survey -- "
                 "precursor_charges: \"separate\" is not fanning out" << std::endl;
  TEST_EQUAL(masses_with_multiple_charges >= 1, true)
}
END_SECTION

// CM-02: ADR-0016's "the two modes differ only in SCAN COUNT, never in which charges" as an
// executable invariant. Both must resolve the same acquisition charge set from the same
// peakGroupNotchCandidates + selectNotches pair; sourcing them differently is precisely how
// "separate" ended up reading an exclusion flag while "multiplexed" read the PeakGroup.
//
// Scoped to the FIRST survey (empty parent_scan_id), which is identical in both runs. Later surveys
// are not comparable across runs: "separate" emits more commands, so it consumes more tracking ids
// and the two runs' survey ids diverge. Comparing them would test id allocation, not geometry.
START_SECTION(separate_and_multiplexed_acquire_the_same_charge_set)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string sep_cfg = cfgFor("separate");
  const std::string mux_cfg = cfgFor("multiplexed");

  FLASHIda sep_ida(const_cast<char*>(sep_cfg.c_str()));
  FLASHIda mux_ida(const_cast<char*>(mux_cfg.c_str()));
  AcqResult sep = runMode(sep_ida, scans);
  AcqResult mux = runMode(mux_ida, scans);

  TEST_EQUAL(sep.ms2_cmds.size() > 0, true)
  TEST_EQUAL(mux.ms2_cmds.size() > 0, true)

  // separate: N commands, anchor only. multiplexed: 1 command, anchor + notches.
  const auto sep_groups = chargesPerSurveyMass(sep, /*with_notches=*/false);
  const auto mux_groups = chargesPerSurveyMass(mux, /*with_notches=*/true);

  int compared = 0, mismatches = 0;
  for (const auto& kv : mux_groups)
  {
    if (!kv.first.first.empty()) continue;  // first survey only
    auto it = sep_groups.find(kv.first);
    if (it == sep_groups.end()) continue;
    ++compared;
    if (it->second == kv.second) continue;
    ++mismatches;
    std::cout << "[CM-02] mass " << kv.first.second << " separate={";
    for (int z : it->second) std::cout << " " << z;
    std::cout << " } multiplexed={";
    for (int z : kv.second) std::cout << " " << z;
    std::cout << " }" << std::endl;
  }

  TEST_EQUAL(compared > 0, true)  // non-vacuity: an empty intersection would pass trivially
  TEST_EQUAL(mismatches, 0)

  // The modes must differ in scan count, or "separate" collapsed into "multiplexed" instead of
  // fanning out and the set comparison above would agree for the wrong reason.
  TEST_EQUAL(sep.ms2_cmds.size() > mux.ms2_cmds.size(), true)
}
END_SECTION

// CM-03: the default mode's within-survey invariant -- a mass is acquired at exactly ONE charge.
// Guards the other direction: a fan-out that fires regardless of mode would spend the whole budget
// on one proteoform's charge states and never fragment P2 or P3.
START_SECTION(single_acquires_one_charge_per_mass_per_survey)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("single");
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)

  const auto groups = chargesPerSurveyMass(a, /*with_notches=*/false);
  int violations = 0;
  for (const auto& kv : groups)
  {
    if (kv.second.size() <= 1) continue;
    ++violations;
    std::cout << "[CM-03] survey '" << kv.first.first << "' mass " << kv.first.second
              << " acquired at " << kv.second.size() << " charges:";
    for (int z : kv.second) std::cout << " " << z;
    std::cout << std::endl;
  }
  TEST_EQUAL(violations, 0)

  // "single" must also command no notches -- one charge means one isolation window.
  int notched = 0;
  for (const auto& c : a.ms2_cmds)
    if (notchesForStage(c, 0).second > 0) ++notched;
  TEST_EQUAL(notched, 0)
}
END_SECTION

END_TEST
