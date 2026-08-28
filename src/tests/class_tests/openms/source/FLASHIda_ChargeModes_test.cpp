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
// That flag is now gone (ADR-0021), so geometry has exactly one source and these sections exercise
// it directly. The lesson generalises past this key: a mode wired at parse time and asserted only
// for round-trip is indistinguishable from a mode that does nothing.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include "FLASHIda_TestHelpers.h"  // ground-truth harness: ScanData/loadTsvScans/runInterleaved

#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
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
  // The AUTHORED CHARGE SET fixtures (ADR-0028). cytC resolves at EVERY charge z=8..20 in this
  // spectrum, so naming three of them is a large, visible restriction rather than a no-op, and the
  // envelope's most intense charge (z=15) is deliberately NOT among them -- an anchor that ignored
  // the set would land outside it.
  const std::string charge_set_list      = "../../FlashIDA/test-data/configs/inclusion_cytc_charges.txt";
  const std::string charge_set_rows_list = "../../FlashIDA/test-data/configs/inclusion_cytc_charges_rows.txt";
  const std::string unrestricted_list    = "../../FlashIDA/test-data/configs/inclusion_cytc.txt";
  const long long   authored_mass        = 12351;   // llround of the authored monoisotopic mass
  const std::set<int> authored_charges{10, 13, 16};

  // The WIDE authored row (ADR-0036), naming a contiguous span the winning PeakGroup often cannot
  // cover on its own. `10;13;16` above cannot exercise split-envelope completion at all: its winner
  // carries all three in 18 of the 25 productive surveys, and in the 7 where it does not, the
  // missing charges were already spent by an earlier survey -- so no sibling ever has anything to
  // complete. Against `12..18` the same spectrum behaves differently: the winner carries all seven
  // in 12 surveys, carries a strict non-empty subset in several more (survey 6 is the first), and
  // the species' PeakGroups collectively carry all seven in 23 of 25.
  const std::string wide_set_list      = "../../FlashIDA/test-data/configs/inclusion_cytc_charges_wide.txt";
  const std::string wide_set_rows_list = "../../FlashIDA/test-data/configs/inclusion_cytc_charges_wide_rows.txt";
  const std::set<int> wide_charges{12, 13, 14, 15, 16, 17, 18};

  bool isSubsetOfAuthored(const std::set<int>& charges)
  {
    for (int z : charges)
      if (authored_charges.count(z) == 0) return false;
    return true;
  }

  // @p max_precursors and @p strict_inclusion default to what every CM-01..CM-08 case already used,
  // so adding them changes no existing config text. CM-09 needs both: the split-envelope completion
  // it exercises is reachable only while the candidate loop is still running (ADR-0036 leaves the
  // cap hard), and under strict inclusion a non-target never spends a slot, which is what makes a
  // budget of 2 enough for one target and its siblings.
  std::string cfgFor(const std::string& precursor_charges,
                     const std::string& targeting = "none",
                     const std::string& inclusion_list = "",
                     int max_precursors = 10,
                     bool strict_inclusion = false)
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
  "files": { "target_logs": [], "fasta": "", "inclusion_list": ")JSON")
           + inclusion_list + R"JSON(", "ptm_list": "" },
  "precursor_selection": {
    "rt_window": 180,
    "targeting": ")JSON" + targeting + R"JSON(",
    "consider_all_charges": false,
    "strict_inclusion": )JSON" + (strict_inclusion ? "true" : "false") + R"JSON(,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": )JSON" + std::to_string(max_precursors) + R"JSON(,
    "precursor_charges": ")JSON" + precursor_charges + R"JSON("
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
  // defect survive. parent_scan_id identifies the spawning survey.
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

  // The FIRST survey's acquisitions only, as mass -> charges.
  //
  // "First survey" is read from the run itself -- the parent of the first MS2 command emitted -- and
  // NOT assumed to be the empty string. buildMS2 stamps parent_scan_id only when the parent tracking
  // id is non-zero (ScanCommandQueue.cpp:345), and by the time the first productive survey runs, the
  // idle AGC/MS1 ticks have already consumed ids, so every MS2 here carries a real parent. Assuming
  // "" matched nothing at all and compared zero groups.
  //
  // Reading it per-run is not incidental: the two runs mint different ids once one emits several
  // times as many commands, so the id is only meaningful within its own run. The first survey itself
  // is comparable across runs because both start from the same engine state and the same MS1.
  std::map<long long, std::set<int>> firstSurveyCharges(const AcqResult& a, bool with_notches)
  {
    std::map<long long, std::set<int>> out;
    if (a.ms2_cmds.empty()) return out;

    const std::string first_survey(a.ms2_cmds.front().parent_scan_id);
    for (const auto& c : a.ms2_cmds)
    {
      if (c.num_stages < 1) continue;
      if (std::string(c.parent_scan_id) != first_survey) continue;
      auto& charges = out[std::llround(c.mono_mass)];
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
// Scoped to the FIRST survey, which both runs reach from the same engine state and the same MS1.
// Later surveys are not comparable across runs: "separate" emits several times as many commands, so
// it consumes more tracking ids and the two runs' survey ids -- and then their exclusion state --
// diverge. Comparing those would test id allocation, not geometry.
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
  const auto sep_groups = firstSurveyCharges(sep, /*with_notches=*/false);
  const auto mux_groups = firstSurveyCharges(mux, /*with_notches=*/true);

  int compared = 0, mismatches = 0;
  for (const auto& kv : mux_groups)
  {
    auto it = sep_groups.find(kv.first);
    if (it == sep_groups.end()) continue;
    ++compared;
    if (it->second == kv.second) continue;
    ++mismatches;
    std::cout << "[CM-02] mass " << kv.first << " separate={";
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

// CM-04: an AUTHORED CHARGE SET is a RESTRICTION (ADR-0028 decision 1). Under multiplexed the scan
// co-isolates the envelope, so if the set were ignored the command would carry ~10 charges, most of
// them unnamed. It must carry only named ones -- and still more than one, or "restrict" quietly
// became "acquire a single charge" and the co-isolation was thrown away with the rest.
START_SECTION(authored_charge_set_restricts_the_acquisition)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("multiplexed", "inclusion", charge_set_list);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)  // non-vacuity

  const auto groups = chargesPerSurveyMass(a, /*with_notches=*/true);
  int checked = 0, violations = 0, co_isolated = 0;
  for (const auto& kv : groups)
  {
    if (kv.first.second != authored_mass) continue;
    ++checked;
    if (kv.second.size() >= 2) ++co_isolated;
    if (isSubsetOfAuthored(kv.second)) continue;
    ++violations;
    std::cout << "[CM-04] survey '" << kv.first.first << "' isolated unauthored charges:";
    for (int z : kv.second)
      if (authored_charges.count(z) == 0) std::cout << " " << z;
    std::cout << " (authored 10;13;16)" << std::endl;
  }

  TEST_EQUAL(checked > 0, true)      // the authored mass was selected at all
  TEST_EQUAL(violations, 0)          // FAILS if the notch filter is skipped
  TEST_EQUAL(co_isolated > 0, true)  // FAILS if the filter emptied the notch set
}
END_SECTION

// CM-05: the anchor is drawn FROM the set (ADR-0028 decision 2). cytC resolves at 13 charges here
// and three are named, so an anchor picked from the envelope rather than the set lands outside it
// almost surely -- which is exactly the regression this guards: the old code reassigned the charge
// inside the inclusion matcher, and that reassignment is now gone.
START_SECTION(authored_charge_set_moves_the_anchor)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("single", "inclusion", charge_set_list);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)

  const auto first = firstSurveyCharges(a, /*with_notches=*/false);
  auto it = first.find(authored_mass);
  TEST_EQUAL(it != first.end(), true)  // non-vacuity: the authored mass must have been acquired
  ABORT_IF(it == first.end())

  // "single" acquires exactly one charge per species per survey, and it must be a named one.
  TEST_EQUAL(it->second.size(), 1)
  if (!isSubsetOfAuthored(it->second))
    std::cout << "[CM-05] anchor " << *it->second.begin()
              << " is not one of the authored charges 10;13;16 -- the set did not reach the anchor"
              << std::endl;
  TEST_EQUAL(isSubsetOfAuthored(it->second), true)
}
END_SECTION

// CM-06: per-charge exclusion walks the set across surveys (ADR-0028 decision 3b) -- the load-
// bearing one. Under mass-keyed exclusion the species is retired by its own first acquisition and
// the second and third named charges are unreachable, so this fails with exactly one anchor.
START_SECTION(authored_charge_set_is_walked_across_surveys)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("single", "inclusion", charge_set_list);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)

  // Anchors for the authored mass, one entry per survey that acquired it.
  std::vector<int> anchors;
  std::set<std::string> surveys_seen;
  for (const auto& c : a.ms2_cmds)
  {
    if (c.num_stages < 1) continue;
    if (std::llround(c.mono_mass) != authored_mass) continue;
    const std::string survey(c.parent_scan_id);
    if (!surveys_seen.insert(survey).second) continue;
    anchors.push_back(c.stages[0].charge_state);
  }

  std::cout << "[CM-06] authored mass acquired at charges:";
  for (int z : anchors) std::cout << " " << z;
  std::cout << " (authored 10;13;16)" << std::endl;

  const std::set<int> distinct(anchors.begin(), anchors.end());
  TEST_EQUAL(anchors.size() >= 2, true)           // FAILS under mass-keyed exclusion
  TEST_EQUAL(distinct.size(), anchors.size())     // FAILS if a spent charge is re-acquired
  TEST_EQUAL(isSubsetOfAuthored(distinct), true)  // FAILS if the walk escapes the set
  TEST_EQUAL(anchors.size() <= 3, true)           // FAILS if the set never terminates
}
END_SECTION

// CM-07: several rows naming one mass union their sets (ADR-0028 decision 5), so the one-row and
// three-row spellings are the same method. It is also the regression guard for the matcher that
// took the first ACTIVE row regardless of its mass: under that bug the three-row file hands over
// the charge of whichever row std::sort happened to leave first, and the two runs diverge.
START_SECTION(rows_naming_one_mass_union_their_charge_sets)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string one_cfg  = cfgFor("multiplexed", "inclusion", charge_set_list);
  const std::string rows_cfg = cfgFor("multiplexed", "inclusion", charge_set_rows_list);

  FLASHIda one_ida(const_cast<char*>(one_cfg.c_str()));
  FLASHIda rows_ida(const_cast<char*>(rows_cfg.c_str()));
  AcqResult one = runMode(one_ida, scans);
  AcqResult rows = runMode(rows_ida, scans);

  TEST_EQUAL(one.ms2_cmds.size() > 0, true)
  TEST_EQUAL(rows.ms2_cmds.size() > 0, true)

  // Scoped to the first survey for the same reason CM-02 is: later surveys are only comparable
  // while both runs have minted the same tracking ids.
  const auto one_first  = firstSurveyCharges(one, /*with_notches=*/true);
  const auto rows_first = firstSurveyCharges(rows, /*with_notches=*/true);

  auto a_it = one_first.find(authored_mass);
  auto b_it = rows_first.find(authored_mass);
  TEST_EQUAL(a_it != one_first.end(), true)  // non-vacuity on both sides
  TEST_EQUAL(b_it != rows_first.end(), true)
  ABORT_IF(a_it == one_first.end() || b_it == rows_first.end())

  if (a_it->second != b_it->second)
  {
    std::cout << "[CM-07] one-row={";
    for (int z : a_it->second) std::cout << " " << z;
    std::cout << " } three-row={";
    for (int z : b_it->second) std::cout << " " << z;
    std::cout << " }" << std::endl;
  }
  TEST_EQUAL(a_it->second == b_it->second, true)
}
END_SECTION

// CM-08: a row naming NO charge stays unrestricted. This is the guard behind the claim that every
// committed config -- all of which write -1 -- is byte-identical. If authoredChargesFor_ ever
// returned a non-empty set for such a row, the whole envelope would collapse onto three charges
// here and every existing inclusion golden would move with it.
START_SECTION(a_row_naming_no_charge_is_unrestricted)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("multiplexed", "inclusion", unrestricted_list);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)

  const auto first = firstSurveyCharges(a, /*with_notches=*/true);
  auto it = first.find(authored_mass);
  TEST_EQUAL(it != first.end(), true)
  ABORT_IF(it == first.end())

  // cytC co-isolates ~10 of its 13 resolved charges here, so an unrestricted row must reach well
  // outside {10,13,16}.
  TEST_EQUAL(isSubsetOfAuthored(it->second), false)
}
END_SECTION

// CM-09: a SPLIT ENVELOPE is completed within one survey (ADR-0036), and the budget cap is what
// bounds it.
//
// The discriminator has to be "two charges of one mass in ONE survey", not "the named set completed
// eventually". Across surveys this fixture's winner varies -- z8..18, then z8..22, then z9,10,11 --
// so the set is walked to completion with or without sibling admission, and an eventual-completion
// assertion would pass against the defect. Within a single survey under `single`, a mass reaching
// two charges is only possible if a second PeakGroup of it was admitted.
//
// The max_precursors 1 arm pins the other half of the decision: the cap stays HARD, so a sibling is
// not examined once the loop has ended, and the run says so rather than silently acquiring less.
START_SECTION(a_split_authored_envelope_is_completed_within_one_survey)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  std::ostringstream captured;
  std::streambuf* const real_cout = std::cout.rdbuf(captured.rdbuf());
  AcqResult roomy, cramped;
  try
  {
    const std::string roomy_cfg = cfgFor("single", "inclusion", wide_set_list, /*max_precursors=*/2,
                                         /*strict_inclusion=*/true);
    FLASHIda roomy_ida(const_cast<char*>(roomy_cfg.c_str()));
    roomy = runMode(roomy_ida, scans);
  }
  catch (...)
  {
    std::cout.rdbuf(real_cout);
    throw;
  }
  const std::string roomy_log = captured.str();
  captured.str("");

  try
  {
    const std::string cramped_cfg = cfgFor("single", "inclusion", wide_set_list, /*max_precursors=*/1,
                                           /*strict_inclusion=*/true);
    FLASHIda cramped_ida(const_cast<char*>(cramped_cfg.c_str()));
    cramped = runMode(cramped_ida, scans);
  }
  catch (...)
  {
    std::cout.rdbuf(real_cout);
    throw;
  }
  const std::string cramped_log = captured.str();
  std::cout.rdbuf(real_cout);

  TEST_EQUAL(roomy.ms2_cmds.size() > 0, true)
  TEST_EQUAL(cramped.ms2_cmds.size() > 0, true)

  // Non-vacuity first: the authored mass must actually have been acquired in both runs, or the
  // counts below would both be zero and the test would pass having measured nothing.
  int roomy_surveys_with_the_mass = 0, cramped_surveys_with_the_mass = 0;
  int roomy_completing_surveys = 0, cramped_completing_surveys = 0;
  std::set<int> roomy_union;

  for (const auto& kv : chargesPerSurveyMass(roomy, /*with_notches=*/false))
  {
    if (kv.first.second != authored_mass) continue;
    ++roomy_surveys_with_the_mass;
    roomy_union.insert(kv.second.begin(), kv.second.end());
    if (kv.second.size() >= 2) ++roomy_completing_surveys;
  }
  for (const auto& kv : chargesPerSurveyMass(cramped, /*with_notches=*/false))
  {
    if (kv.first.second != authored_mass) continue;
    ++cramped_surveys_with_the_mass;
    if (kv.second.size() >= 2) ++cramped_completing_surveys;
  }

  TEST_EQUAL(roomy_surveys_with_the_mass > 0, true)
  TEST_EQUAL(cramped_surveys_with_the_mass > 0, true)

  if (roomy_completing_surveys == 0)
    std::cout << "[CM-09] no survey acquired mass " << authored_mass << " at >=2 charges -- a "
              << "sibling PeakGroup was never admitted, so the split envelope was not completed"
              << std::endl;
  TEST_EQUAL(roomy_completing_surveys > 0, true)

  // Everything acquired is a NAMED charge: completion must not become "acquire more".
  for (int z : roomy_union)
  {
    if (wide_charges.count(z) == 0)
      std::cout << "[CM-09] charge " << z << " is outside the authored set 12..18" << std::endl;
    TEST_EQUAL(wide_charges.count(z) > 0, true)
  }

  // The hard cap: at max_precursors 1 the loop ends on the first species and no sibling is reached.
  TEST_EQUAL(cramped_completing_surveys, 0)

  // ...and it is not silent about it. This is the whole reason the cap was left hard rather than
  // special-cased: the configuration remains legal and unchanged, and the run reports what it cost.
  TEST_EQUAL(cramped_log.find("reason=budget_exhausted") != std::string::npos, true)
  TEST_EQUAL(roomy_log.find("reason=budget_exhausted") == std::string::npos, true)
}
END_SECTION

// CM-10: the union rule (ADR-0028) survives a SPLIT. CM-07 already pins that one row naming
// `10;13;16` and three rows naming one charge each acquire identically -- but its winner resolves
// the whole named set every survey, so no sibling is ever admitted and the equivalence never
// crosses the completion path. Against the wide row it does.
START_SECTION(rows_naming_one_mass_union_their_charge_sets_under_a_split)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string one_cfg  = cfgFor("single", "inclusion", wide_set_list, 2, true);
  const std::string rows_cfg = cfgFor("single", "inclusion", wide_set_rows_list, 2, true);

  FLASHIda one_ida(const_cast<char*>(one_cfg.c_str()));
  FLASHIda rows_ida(const_cast<char*>(rows_cfg.c_str()));
  AcqResult one  = runMode(one_ida, scans);
  AcqResult rows = runMode(rows_ida, scans);

  TEST_EQUAL(one.ms2_cmds.size() > 0, true)
  TEST_EQUAL(one.ms2_cmds.size(), rows.ms2_cmds.size())
  ABORT_IF(one.ms2_cmds.size() != rows.ms2_cmds.size())

  const auto one_groups  = chargesPerSurveyMass(one, /*with_notches=*/true);
  const auto rows_groups = chargesPerSurveyMass(rows, /*with_notches=*/true);
  TEST_EQUAL(one_groups.size(), rows_groups.size())
  TEST_EQUAL(one_groups == rows_groups, true)
}
END_SECTION

// CM-11: an UNRESTRICTED row is untouched by any of this (ADR-0036 as amended; ADR-0037 withdrawn).
//
// This is the regression guard for the scope decision, and it is the test whose absence let a
// change aimed at charge completeness move thirteen log goldens -- every one of them a `-1` row,
// while the two authored-charge modes stayed still. Two things must hold for a `-1` row: no sibling
// is ever admitted (so a mass reaches at most one anchor per survey under `single`), and matching
// an inclusion row does not lift dynamic exclusion (so the mass is acquired once per rt window
// rather than on every survey that beats the stored qscore by any margin).
START_SECTION(an_unrestricted_row_is_unaffected_by_split_envelope_completion)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  // Under `single`: one anchor per mass per survey, exactly as CM-03 requires of a non-targeted run.
  {
    const std::string cfg = cfgFor("single", "inclusion", unrestricted_list);
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    AcqResult a = runMode(ida, scans);
    TEST_EQUAL(a.ms2_cmds.size() > 0, true)

    int violations = 0, surveys_with_the_mass = 0;
    for (const auto& kv : chargesPerSurveyMass(a, /*with_notches=*/false))
    {
      if (kv.first.second != authored_mass) continue;
      ++surveys_with_the_mass;
      if (kv.second.size() <= 1) continue;
      ++violations;
      std::cout << "[CM-11] survey '" << kv.first.first << "' acquired the unrestricted mass at "
                << kv.second.size() << " charges -- a sibling was admitted for a `-1` row"
                << std::endl;
    }
    TEST_EQUAL(surveys_with_the_mass > 0, true)   // non-vacuity
    TEST_EQUAL(violations, 0)

    // Dynamic exclusion still governs: the target is acquired ONCE inside the rt window, not on
    // every survey that resolves it a fraction better. rt_window is 180 s and this fixture is a
    // 63 s run, so nothing expires -- a count above one here means the bars were exempted.
    TEST_EQUAL(surveys_with_the_mass, 1)
  }

  // Under `multiplexed`: still one scan for the mass per survey, however many PeakGroups of it the
  // survey produced. The notch set is the winner's own envelope and is NOT extended by a sibling.
  {
    const std::string cfg = cfgFor("multiplexed", "inclusion", unrestricted_list);
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    AcqResult a = runMode(ida, scans);
    TEST_EQUAL(a.ms2_cmds.size() > 0, true)

    int scans_for_the_mass = 0;
    for (const auto& c : a.ms2_cmds)
      if (c.num_stages >= 1 && std::llround(c.mono_mass) == authored_mass) ++scans_for_the_mass;

    TEST_EQUAL(scans_for_the_mass, 1)
  }
}
END_SECTION

END_TEST
