// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// characterization.mode: exhaustive (ADR-0023) -- the mode whose MS3 target pool is EVERY deconvolved
// mass of the winner MS2 scan, not only the masses that matched the winning proteoform.
//
// Everything here is driven through the tracker's public sequence -- feedScan -> finalizeMS2 ->
// planNextScans -- exactly as the neighbouring ProteoformTracker_* tests do. There are no *ForTest
// hooks and none may be added: test scaffolding lives at the test location.
//
// CONSTRUCTION NOTE, and it constrains what can honestly be asserted here. A hand-built PeakGroup has
// never been through a scoring pass, so PeakGroup::getChargeIntensity() returns 0 for every charge --
// which is what feedScan copies into PeakRecord::intensity. EVERY pool entry in this file therefore
// carries intensity 0, and the intensity ORDER of the plan is not observable from C++. It is pinned
// by the C# golden suite instead, per the standing division of labour (C++ asserts plausibility and
// structure, C# asserts exact numbers). What IS observable, and is what section 3 pins, is that the
// pool is ranked as ONE list rather than mapped-first-then-unassigned: with the intensities all tied
// the planner's stable_sort leaves the winner scan's own peak order, so a mapped target must still be
// found sandwiched between two unassigned ones.

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Exploration.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms2Params.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include "FLASHIda_TestAccess.h"   // ExplorationTestAccess::group -- Exploration's group map is private

#include <cmath>
#include <cstring>
#include <set>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  // 8-residue winner proteoform, the same one the other ProteoformTracker tests use.
  const char* WINNER_SEQ = "PEPTIDEK";

  // One ambiguous modification both mapped fragments CONTAIN and neither BRACKETS, so it survives
  // narrowModifications_ and the Ambiguity objective still has work to do in section 7:
  //   b6 covers region [1,6] -> contains [4,5]; brackets iff 4 <= 6 < 5  -> false
  //   y6 covers region [3,8] -> contains [4,5]; brackets iff 4 <  3 <= 5 -> false
  const int AMB_MOD_START = 4;
  const int AMB_MOD_END = 5;

  // The two masses the winner match claims as real fragments...
  const double MASS_B6 = 700.0;
  const double MASS_Y6 = 900.0;
  // ...and five masses nothing claims. Deliberately far apart, so no two share a NOMINAL mass
  // (round(m * 0.999497), the key the dispatch memory uses) and none lands within the 10 ppm level-2
  // tolerance of a mapped one.
  const double MASS_U1 = 1234.5;
  const double MASS_U2 = 2345.6;
  const double MASS_U3 = 3456.7;
  const double MASS_U4 = 4567.8;
  const double MASS_U5 = 5678.9;

  // 1.5 Da off b6: ~2100 ppm away, i.e. 200x the level-2 tolerance, yet still the CLOSEST mapped
  // fragment to it. The probe for ADR-0023 D-c -- a "closest overall" label lookup (the rule
  // mapScanOntoModel_ deliberately uses for its own, opposite question) would stamp this peak `b6`.
  const double MASS_NEAR_MISS = MASS_B6 + 1.5;

  // A config that differs from the neighbouring tracker fixtures only in the characterization block.
  // ms_settings.ms3 is present with a CID activation because Config::validate() requires both as soon
  // as mode != off.
  //
  // The raw strings carry an explicit JSON delimiter. An undelimited R"(...)" ends at the first `)"`
  // in the text, and a single parenthesised value dropped into this JSON later would silently truncate
  // it mid-document -- that has already cost one CI round trip in this repo.
  const char* kConfigHead = R"JSON({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.9,
    "min_charge": 1,
    "max_charge": 50,
    "min_mass": 100,
    "max_mass": 50000,
    "tol": [10, 10, 10]
  },
  "flashtnt": {"min_length": 3, "max_length": 8},
  "faims": {"cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15},
  "scheduling": {
    "cycle_time": {"enabled": false, "value_ms": 60000},
    "scan_timeout": {"enabled": false, "value_ms": 30000}
  },
  "files": {"target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": ""},
  "conditional_ms2": false,
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3
  },
  "characterization": {
    "mode": ")JSON";

  const char* kConfigTail = R"JSON(
  },
  "ms_settings": {
    "ms1": {"analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000,
            "resolution": 120000, "agc_target": 800000, "max_it": 246},
    "ms2": {"analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000},
    "ms3": {"analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25, "resolution": 120000}
  },
  "tagging": {},
  "quantification": {"enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4}
}
)JSON";

  std::string trackerConfig(const std::string& mode, int max_targets,
                            double min_target_mass, int min_fragment_charge)
  {
    return std::string(kConfigHead) + mode
         + "\",\n    \"protein_sequence\": \"" + WINNER_SEQ
         + "\",\n    \"max_targets\": " + std::to_string(max_targets)
         + ",\n    \"min_target_mass\": " + std::to_string(min_target_mass)
         + ",\n    \"min_fragment_charge\": " + std::to_string(min_fragment_charge)
         + kConfigTail;
  }

  // One synthetic deconvolved mass. Single peak, so getMzRange(charge) is a zero-width window and
  // PeakRecord::by_charge comes out empty -- which is exactly the fragment_charges == Single state
  // every section here runs in.
  PeakGroup makeSyntheticPeakGroup(double mono_mass, int charge)
  {
    PeakGroup pg(charge, charge, true);
    pg.setMonoisotopicMass(mono_mass);
    FLASHHelperClasses::LogMzPeak lp;
    lp.mz = mono_mass / charge + 1.00728;
    lp.abs_charge = charge;
    lp.intensity = 1000.0f;   // recorded on the peak; getChargeIntensity stays 0 (see the file header)
    pg.push_back(lp);
    return pg;
  }

  struct FragSpec
  {
    std::string ion_type;
    int ion_index;
    double observed_mass;
  };

  // A winner ProteoformMatch: score >= 0 plus a non-empty sequence makes finalizeMS2 accept it, and
  // its FragmentMatches are pooled VERBATIM (the winner scan is not re-matched), so the masses need
  // not sit on a theoretical ladder for this fixture to behave like a real identification.
  FragmentAnalysis::ProteoformMatch makeWinnerMatch(const std::vector<FragSpec>& frags, bool with_ambiguous_ptm)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.score = 1.0;
    m.region_start = -1;   // full-sequence region -> identity frame mapping
    m.region_end = -1;
    m.matched_protein = "synthetic";
    m.proteoform_sequence = WINNER_SEQ;
    for (const FragSpec& f : frags)
    {
      FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
      fm.ion_type = f.ion_type;
      fm.ion_index = f.ion_index;
      fm.observed_mass = f.observed_mass;
      m.fragments.push_back(fm);
    }
    if (with_ambiguous_ptm)
    {
      FragmentAnalysis::PTMSite site;
      site.start_position = AMB_MOD_START;
      site.end_position = AMB_MOD_END;
      site.position = (AMB_MOD_START + AMB_MOD_END) / 2;
      site.mass_shift = 42.0106;
      m.ptm_sites.push_back(site);
    }
    return m;
  }

  ScanCommand makeMs2Ctx()
  {
    ScanCommand ctx{};
    ctx.scan_id = 1;
    ctx.msn_level = 2;
    ctx.num_stages = 1;
    ctx.mono_mass = 12000.0;
    ctx.stages[0].precursor_mz = 1001.0;
    ctx.stages[0].isolation_width = 2.0;
    ctx.stages[0].charge_state = 12;
    ctx.stages[0].collision_energy = 29.0;
    std::strncpy(ctx.stages[0].activation_type, "HCD", sizeof(ctx.stages[0].activation_type) - 1);
    return ctx;
  }

  struct MassCharge
  {
    double mass;
    int charge;
  };

  // Feed ONE MS2 scan carrying @p peaks and identifying @p frags, then finalize. One scan means it is
  // the winner by construction, so mapScanOntoModel_ takes the winner-pooling path (a non-winner scan
  // would be re-matched against the winner ladder instead, which these synthetic masses would fail).
  void feedWinnerScan(ProteoformTracker& tracker, int precursor_id, int scan_id,
                      const std::vector<MassCharge>& peaks, const std::vector<FragSpec>& frags,
                      const std::string& activation, bool with_ambiguous_ptm = false)
  {
    Ms2Params p;
    p.collision_energy = 29.0;
    p.activation_type = activation;
    p.reaction_time = 0.0;

    DeconvolvedSpectrum d(scan_id);
    for (const MassCharge& mc : peaks) d.push_back(makeSyntheticPeakGroup(mc.mass, mc.charge));

    ScanCommand ms2_ctx = makeMs2Ctx();
    tracker.feedScan(precursor_id, 2, p, scan_id, d, makeWinnerMatch(frags, with_ambiguous_ptm), 1.0, ms2_ctx);
    tracker.finalizeMS2(precursor_id);
  }

  const Ms3Target* targetAtMass(const std::vector<Ms3Target>& targets, double mass)
  {
    for (const Ms3Target& t : targets)
      if (std::abs(t.frag_mass - mass) < 1e-6) return &t;
    return nullptr;
  }

  int countUnassigned(const std::vector<Ms3Target>& targets)
  {
    int n = 0;
    for (const Ms3Target& t : targets)
      if (t.ion_type == "u" && t.ion_index == 0) ++n;
    return n;
  }

  bool dispatched(const ProteoformModel& m, double mass)
  {
    return m.dispatched_nominal_masses.count(SpectralDeconvolution::getNominalMass(mass)) > 0;
  }

  // The seven-mass winner spectrum used by most sections: two masses the match claims, five it does
  // not. Interleaved so no section can accidentally depend on mapped masses coming first.
  const std::vector<MassCharge> SEVEN_MASSES{
    {MASS_U1, 2}, {MASS_B6, 2}, {MASS_U2, 2}, {MASS_Y6, 2}, {MASS_U3, 2}, {MASS_U4, 2}, {MASS_U5, 2}};
  const std::vector<FragSpec> TWO_FRAGS{{"b", 6, MASS_B6}, {"y", 6, MASS_Y6}};
}

START_TEST(ProteoformTracker_Exhaustive, "$Id$")

/////////////////////////////////////////////////////////////
// 1. The core claim: the pool is the winner scan's deconvolved masses, not its mapped fragments.
//
// FAILS UNDER: any implementation that still walks m.fragments -- it would return the 2 mapped
// fragments instead of all 7 masses. Also fails if the exhaustive fork is never taken (a mode that
// assigns CharacterizationMode but not CharacterizationObjective ships byte-identical to ambiguity,
// ADR-0023 D-a: that would return b6 + y6 here, i.e. 2).
/////////////////////////////////////////////////////////////
START_SECTION(pool_is_every_deconvolved_mass)
{
  const int precursor_id = 21;
  Config cfg{trackerConfig("exhaustive", 7, 0.0, 0)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);

  feedWinnerScan(tracker, precursor_id, 501, SEVEN_MASSES, TWO_FRAGS, "HCD");

  const ProteoformModel* mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL(mdl->proteoform_sequence, std::string(WINNER_SEQ))
  // Exactly 2 of the 7 masses mapped -- the very asymmetry the mode exists to exploit.
  TEST_EQUAL(static_cast<int>(mdl->fragments.size()), 2)

  std::vector<Ms3Target> plan = tracker.planNextScans(precursor_id);
  TEST_EQUAL(static_cast<int>(plan.size()), 7)   // N, not M
  ABORT_IF(plan.size() != 7)

  // Every fed mass is targeted exactly once, mapped or not.
  for (const MassCharge& mc : SEVEN_MASSES) TEST_TRUE(targetAtMass(plan, mc.mass) != nullptr)
  std::set<int> distinct;
  for (const Ms3Target& t : plan) distinct.insert(SpectralDeconvolution::getNominalMass(t.frag_mass));
  TEST_EQUAL(static_cast<int>(distinct.size()), 7)

  // Every target is isolatable and carries the winner scan's own acquisition as its MS3 stage[0]
  // (ADR-0003): for an unassigned mass there is no per-ion best observation to source it from.
  for (const Ms3Target& t : plan)
  {
    TEST_TRUE(t.frag_mz > 0.0)
    TEST_EQUAL(t.frag_charge, 2)
    TEST_EQUAL(t.stage0_params.activation_type, std::string("HCD"))
    TEST_TRUE(std::abs(t.stage0_params.collision_energy - 29.0) < 1e-6)
    TEST_EQUAL(static_cast<int>(t.notches.size()), 0)   // fragment_charges defaults to single
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// 2. Labelling: a mass that mapped keeps its real ion identity, a mass that did not gets 'u'/0.
//
// FAILS UNDER: (a) labelling everything 'u' -- which is what keying the lookup on
// MappedFragment::theoretical_mass produces, since that member is declared and never assigned;
// (b) a "closest overall" label lookup copied from mapScanOntoModel_ (ADR-0023 D-c) -- then the
// near-miss peak 1.5 Da off b6, and in fact every unassigned mass, would inherit a real ion label and
// the unassigned count would collapse; (c) pairing the sentinel with a non-zero index, which would
// put a 'u' on the wire instead of down buildMS3's no-ion branch (D-f).
/////////////////////////////////////////////////////////////
START_SECTION(unassigned_targets_are_labelled_u0)
{
  const int precursor_id = 22;
  Config cfg{trackerConfig("exhaustive", 10, 0.0, 0)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);

  std::vector<MassCharge> peaks = SEVEN_MASSES;
  peaks.push_back({MASS_NEAR_MISS, 2});   // the D-c probe
  feedWinnerScan(tracker, precursor_id, 502, peaks, TWO_FRAGS, "HCD");

  std::vector<Ms3Target> plan = tracker.planNextScans(precursor_id);
  TEST_EQUAL(static_cast<int>(plan.size()), 8)
  ABORT_IF(plan.size() != 8)

  const Ms3Target* mapped_b = targetAtMass(plan, MASS_B6);
  const Ms3Target* mapped_y = targetAtMass(plan, MASS_Y6);
  TEST_TRUE(mapped_b != nullptr)
  TEST_TRUE(mapped_y != nullptr)
  ABORT_IF(mapped_b == nullptr || mapped_y == nullptr)
  TEST_EQUAL(mapped_b->ion_type, std::string("b"))
  TEST_EQUAL(mapped_b->ion_index, 6)
  TEST_EQUAL(mapped_y->ion_type, std::string("y"))
  TEST_EQUAL(mapped_y->ion_index, 6)

  // The other six -- the five unclaimed masses plus the near-miss -- carry the sentinel and index 0.
  TEST_EQUAL(countUnassigned(plan), 6)
  const Ms3Target* near_miss = targetAtMass(plan, MASS_NEAR_MISS);
  TEST_TRUE(near_miss != nullptr)
  ABORT_IF(near_miss == nullptr)
  TEST_EQUAL(near_miss->ion_type, std::string("u"))
  TEST_EQUAL(near_miss->ion_index, 0)

  // The whole-plan invariant: a target is either a KNOWN ion class with a positive index, or exactly
  // the 'u'/0 sentinel. Nothing in between, because index and class travel independently and only the
  // pair (sentinel class, index 0) is safe at both the matcher and the wire (ADR-0023 5a/5b).
  for (const Ms3Target& t : plan)
  {
    if (t.ion_type == "u") { TEST_EQUAL(t.ion_index, 0) }
    else
    {
      TEST_TRUE(!t.ion_type.empty())
      TEST_TRUE(t.ion_index > 0)
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// 3. One flat ranking: mapped and unassigned masses compete in the SAME list.
//
// The intensities are all tied here (see the file header), so what this pins is that the planner ranks
// one pool with a stable order rather than concatenating two class-segregated lists.
//
// FAILS UNDER: the natural "extend the existing path" implementation -- emit the mapped fragments
// first, then append the leftovers. That gives [b6, y6, U1] for a budget of 3; this asserts
// [U1, b6, U2], i.e. a mapped target sandwiched between two unassigned ones. It also fails under any
// re-ordering of the pool by mass or by class.
/////////////////////////////////////////////////////////////
START_SECTION(ranking_is_one_pool_not_segregated_by_class)
{
  const int precursor_id = 23;
  Config cfg{trackerConfig("exhaustive", 3, 0.0, 0)};   // budget below the pool -> order decides
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);

  // Winner-scan peak order: U1, b6, U2, y6, U3.
  const std::vector<MassCharge> interleaved{
    {MASS_U1, 2}, {MASS_B6, 2}, {MASS_U2, 2}, {MASS_Y6, 2}, {MASS_U3, 2}};
  feedWinnerScan(tracker, precursor_id, 503, interleaved, TWO_FRAGS, "HCD");

  std::vector<Ms3Target> plan = tracker.planNextScans(precursor_id);
  TEST_EQUAL(static_cast<int>(plan.size()), 3)
  ABORT_IF(plan.size() != 3)

  TEST_TRUE(std::abs(plan[0].frag_mass - MASS_U1) < 1e-6)
  TEST_EQUAL(plan[0].ion_type, std::string("u"))
  TEST_TRUE(std::abs(plan[1].frag_mass - MASS_B6) < 1e-6)
  TEST_EQUAL(plan[1].ion_type, std::string("b"))     // a MAPPED target, in the middle of the plan
  TEST_EQUAL(plan[1].ion_index, 6)
  TEST_TRUE(std::abs(plan[2].frag_mass - MASS_U2) < 1e-6)
  TEST_EQUAL(plan[2].ion_type, std::string("u"))

  // y6 is mapped and still LOST its slot to two unassigned masses ranked ahead of it -- decision 2,
  // "exhaustive is not a superset of ambiguity".
  TEST_TRUE(targetAtMass(plan, MASS_Y6) == nullptr)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 4. Filters run BEFORE the dispatch memory is stamped (ADR-0023 D-d).
//
// The stake is not just "the filtered mass is absent from this plan" -- it is that the mass is still
// REACHABLE afterwards. The memory is monotone and dispatched-but-never-returned counts as done, so a
// stamp made before the filters burns that species for the rest of the Precursor's life.
//
// FAILS UNDER: stamping the nominal mass while building the pool (or before applying the floors),
// which is the shape the ADR's pseudocode invites -- then the filtered masses appear in
// dispatched_nominal_masses even though nothing was ever acquired for them.
/////////////////////////////////////////////////////////////
START_SECTION(filters_apply_before_the_dispatch_memory_is_stamped)
{
  // --- min_target_mass ---------------------------------------------------------------------------
  {
    const int precursor_id = 24;
    Config cfg{trackerConfig("exhaustive", 10, 3000.0, 0)};   // floor above b6, y6, U1, U2
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedWinnerScan(tracker, precursor_id, 504, SEVEN_MASSES, TWO_FRAGS, "HCD");

    std::vector<Ms3Target> plan = tracker.planNextScans(precursor_id);
    TEST_EQUAL(static_cast<int>(plan.size()), 3)   // only U3, U4, U5 clear 3000 Da
    TEST_TRUE(targetAtMass(plan, MASS_B6) == nullptr)
    TEST_TRUE(targetAtMass(plan, MASS_Y6) == nullptr)
    TEST_TRUE(targetAtMass(plan, MASS_U1) == nullptr)
    TEST_TRUE(targetAtMass(plan, MASS_U2) == nullptr)
    TEST_TRUE(targetAtMass(plan, MASS_U3) != nullptr)

    const ProteoformModel* mdl = tracker.getModel(precursor_id);
    TEST_TRUE(mdl != nullptr)
    ABORT_IF(mdl == nullptr)
    // THE ASSERTION THAT MATTERS: the four masses the floor rejected were never stamped.
    TEST_EQUAL(static_cast<int>(mdl->dispatched_nominal_masses.size()), 3)
    TEST_EQUAL(dispatched(*mdl, MASS_B6), false)
    TEST_EQUAL(dispatched(*mdl, MASS_Y6), false)
    TEST_EQUAL(dispatched(*mdl, MASS_U1), false)
    TEST_EQUAL(dispatched(*mdl, MASS_U2), false)
    TEST_EQUAL(dispatched(*mdl, MASS_U3), true)
  }

  // --- min_fragment_charge, enforced downstream at Exploration.cpp:950 and therefore the harder case
  {
    const int precursor_id = 25;
    Config cfg{trackerConfig("exhaustive", 10, 0.0, 4)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);

    // Charges 2 and 3 are below the floor; 4, 5 and 6 clear it.
    const std::vector<MassCharge> charged{
      {MASS_U1, 2}, {MASS_B6, 2}, {MASS_U2, 3}, {MASS_Y6, 3},
      {MASS_U3, 4}, {MASS_U4, 5}, {MASS_U5, 6}};
    feedWinnerScan(tracker, precursor_id, 505, charged, TWO_FRAGS, "HCD");

    std::vector<Ms3Target> plan = tracker.planNextScans(precursor_id);
    TEST_EQUAL(static_cast<int>(plan.size()), 3)
    TEST_TRUE(targetAtMass(plan, MASS_B6) == nullptr)   // mapped, but charge 2
    TEST_TRUE(targetAtMass(plan, MASS_U2) == nullptr)   // unassigned, charge 3
    TEST_TRUE(targetAtMass(plan, MASS_U3) != nullptr)

    const ProteoformModel* mdl = tracker.getModel(precursor_id);
    TEST_TRUE(mdl != nullptr)
    ABORT_IF(mdl == nullptr)
    TEST_EQUAL(static_cast<int>(mdl->dispatched_nominal_masses.size()), 3)
    TEST_EQUAL(dispatched(*mdl, MASS_B6), false)
    TEST_EQUAL(dispatched(*mdl, MASS_U2), false)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// 5. The dispatch memory is what terminates the mode (decision 7).
//
// Exhaustive has no feedback that shrinks its pool -- a localized mod leaves the ambiguous list and a
// witnessed bond leaves the uncovered set, but a deconvolved mass stays deconvolved. This set IS the
// termination condition.
//
// FAILS UNDER: no memory at all (call 2 returns the same three masses), a memory keyed on something
// other than the nominal mass, or a memory that is reset per plan.
/////////////////////////////////////////////////////////////
START_SECTION(dispatch_memory_suppresses_duplicates_and_terminates)
{
  const int precursor_id = 26;
  Config cfg{trackerConfig("exhaustive", 3, 0.0, 0)};   // budget 3, pool 7
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  // The ambiguous modification is seeded ON PURPOSE here, and it is what makes the objectiveUnmet
  // assertion at the end of this section non-vacuous: it is never localized, so the AMBIGUITY rule
  // would answer "unmet" forever. Only the Exhaustive arm can answer "met" once the pool is drained.
  feedWinnerScan(tracker, precursor_id, 506, SEVEN_MASSES, TWO_FRAGS, "HCD", /*with_ambiguous_ptm=*/true);

  auto nominalsOf = [](const std::vector<Ms3Target>& ts) {
    std::set<int> out;
    for (const Ms3Target& t : ts) out.insert(SpectralDeconvolution::getNominalMass(t.frag_mass));
    return out;
  };

  std::vector<Ms3Target> first = tracker.planNextScans(precursor_id);
  TEST_EQUAL(static_cast<int>(first.size()), 3)
  std::vector<Ms3Target> second = tracker.planNextScans(precursor_id);
  TEST_EQUAL(static_cast<int>(second.size()), 3)

  // Disjoint: the second plan re-picked nothing the first one took.
  std::set<int> a = nominalsOf(first);
  std::set<int> b = nominalsOf(second);
  TEST_EQUAL(static_cast<int>(a.size()), 3)
  TEST_EQUAL(static_cast<int>(b.size()), 3)
  for (int n : b) TEST_EQUAL(a.count(n) == 0, true)

  // 7 masses, 3 + 3 taken -> one left, then the pool is exhausted and stays exhausted.
  std::vector<Ms3Target> third = tracker.planNextScans(precursor_id);
  TEST_EQUAL(static_cast<int>(third.size()), 1)
  std::vector<Ms3Target> fourth = tracker.planNextScans(precursor_id);
  TEST_EQUAL(static_cast<int>(fourth.size()), 0)   // pool_exhausted
  std::vector<Ms3Target> fifth = tracker.planNextScans(precursor_id);
  TEST_EQUAL(static_cast<int>(fifth.size()), 0)    // and it stays exhausted

  const ProteoformModel* mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL(static_cast<int>(mdl->dispatched_nominal_masses.size()), 7)

  // NOTE: the matching objectiveUnmet assertion -- that a drained pool reports "met" so the escalation
  // ladder stops -- is deliberately NOT here. objectiveUnmet belongs to the ladder (ADR-0022), which is
  // not in this branch; the assertion lands with the ladder's own push. `pool_exhausted` above is the
  // same termination fact stated through the API that does exist here.
}
END_SECTION

/////////////////////////////////////////////////////////////
// 6. The winner scan's MS3-capability gate, and its deliberate asymmetry (decision 3 + D-b).
//
// An MS3's stage[0] REPLAYS the MS2 that made the fragment (ADR-0003), so an ETD winner makes every
// mass in its scan MS3-incapable. But an EMPTY activation is "not recorded", not "incapable".
//
// FAILS UNDER: no gate at all (the ETD case plans 7), or the gate written as the ADR's wording reads
// -- fail-closed on "" -- which returns zero targets for every hand-built fixture and for any scan
// config that omits the key. Both halves are asserted, because each is the other's regression.
/////////////////////////////////////////////////////////////
START_SECTION(etd_winner_plans_nothing_but_an_unrecorded_activation_still_does)
{
  // --- ETD: refuse, and burn nothing ---------------------------------------------------------------
  {
    const int precursor_id = 27;
    Config cfg{trackerConfig("exhaustive", 10, 0.0, 0)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedWinnerScan(tracker, precursor_id, 507, SEVEN_MASSES, TWO_FRAGS, "ETD");

    std::vector<Ms3Target> plan = tracker.planNextScans(precursor_id);
    TEST_EQUAL(static_cast<int>(plan.size()), 0)

    const ProteoformModel* mdl = tracker.getModel(precursor_id);
    TEST_TRUE(mdl != nullptr)
    ABORT_IF(mdl == nullptr)
    // A refusal must cost nothing: no nominal mass may be stamped by a plan that emitted no target.
    TEST_EQUAL(static_cast<int>(mdl->dispatched_nominal_masses.size()), 0)
    // (The companion assertion -- that the objective stays UNMET here so an MS2 rung can still rescue
    // this Precursor -- needs objectiveUnmet, which belongs to the escalation ladder and is not in this
    // branch. It lands with the ladder.)
  }

  // --- "" (no activation recorded): CAPABLE ---------------------------------------------------------
  {
    const int precursor_id = 28;
    Config cfg{trackerConfig("exhaustive", 10, 0.0, 0)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedWinnerScan(tracker, precursor_id, 508, SEVEN_MASSES, TWO_FRAGS, "");

    std::vector<Ms3Target> plan = tracker.planNextScans(precursor_id);
    TEST_EQUAL(static_cast<int>(plan.size()), 7)
  }

  // --- CID: the other classified-capable activation --------------------------------------------------
  {
    const int precursor_id = 29;
    Config cfg{trackerConfig("exhaustive", 10, 0.0, 0)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedWinnerScan(tracker, precursor_id, 509, SEVEN_MASSES, TWO_FRAGS, "CID");

    std::vector<Ms3Target> plan = tracker.planNextScans(precursor_id);
    TEST_EQUAL(static_cast<int>(plan.size()), 7)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// 7. THE REGRESSION GUARD FOR THE 21 EXISTING GOLDEN MODE DIRECTORIES.
//
// The identical fixture -- same peaks, same identification, same feed sequence -- under ambiguity and
// under coverage must still return exactly what it returns on main: the two MAPPED fragments and
// nothing else. The five unassigned masses sit right there in the winner scan and must stay invisible
// to both objectives.
//
// FAILS UNDER: a fork placed above the objective test, an objective that defaults to Exhaustive, or a
// pool build that leaked the raw peak list into the mapped-fragment path. Any of those turns 2 targets
// into 7 and moves every ms3_* golden.
/////////////////////////////////////////////////////////////
START_SECTION(ambiguity_and_coverage_are_untouched)
{
  // --- ambiguity: both fragments CONTAIN the still-ambiguous modification --------------------------
  {
    const int precursor_id = 30;
    Config cfg{trackerConfig("ambiguity", 10, 0.0, 0)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedWinnerScan(tracker, precursor_id, 510, SEVEN_MASSES, TWO_FRAGS, "HCD", /*with_ambiguous_ptm=*/true);

    const ProteoformModel* mdl = tracker.getModel(precursor_id);
    TEST_TRUE(mdl != nullptr)
    ABORT_IF(mdl == nullptr)
    TEST_EQUAL(static_cast<int>(mdl->modifications.size()), 1)

    std::vector<Ms3Target> plan = tracker.planNextScans(precursor_id);
    TEST_EQUAL(static_cast<int>(plan.size()), 2)
    TEST_TRUE(targetAtMass(plan, MASS_B6) != nullptr)
    TEST_TRUE(targetAtMass(plan, MASS_Y6) != nullptr)
    TEST_EQUAL(countUnassigned(plan), 0)
    for (const MassCharge& mc : SEVEN_MASSES)
    {
      if (mc.mass == MASS_B6 || mc.mass == MASS_Y6) continue;
      TEST_TRUE(targetAtMass(plan, mc.mass) == nullptr)
    }
    // The exhaustive-only state stays untouched in the other modes: they are self-limiting by
    // construction and never consult it.
    TEST_EQUAL(static_cast<int>(mdl->dispatched_nominal_masses.size()), 0)
  }

  // --- coverage: both fragments' interiors add unwitnessed backbone bonds ---------------------------
  // b6 witnesses bond 6 and y6 bond 2 (L = 8), leaving {1,3,4,5,7} uncovered; b6's interior is bonds
  // 1..5 and y6's is 3..7, so each adds at least one gap whichever order the tie-broken sort produces.
  {
    const int precursor_id = 31;
    Config cfg{trackerConfig("coverage", 10, 0.0, 0)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedWinnerScan(tracker, precursor_id, 511, SEVEN_MASSES, TWO_FRAGS, "HCD");

    std::vector<Ms3Target> plan = tracker.planNextScans(precursor_id);
    TEST_EQUAL(static_cast<int>(plan.size()), 2)
    TEST_TRUE(targetAtMass(plan, MASS_B6) != nullptr)
    TEST_TRUE(targetAtMass(plan, MASS_Y6) != nullptr)
    TEST_EQUAL(countUnassigned(plan), 0)
    for (const MassCharge& mc : SEVEN_MASSES)
    {
      if (mc.mass == MASS_B6 || mc.mass == MASS_Y6) continue;
      TEST_TRUE(targetAtMass(plan, mc.mass) == nullptr)
    }

    const ProteoformModel* mdl = tracker.getModel(precursor_id);
    TEST_TRUE(mdl != nullptr)
    ABORT_IF(mdl == nullptr)
    TEST_EQUAL(static_cast<int>(mdl->dispatched_nominal_masses.size()), 0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// 8. An unassigned target's MS3 CE sweep is scored by remaining_precursor, whatever the config asked
//    for (ADR-0023 decision 11) -- and a MAPPED target in the same mode is not.
//
// The override is per-TARGET, not per-mode: what disqualifies a reading metric is the missing ion
// frame, which only an unassigned target lacks.
//
// FAILS UNDER: no override at all -- the configured fragment_count would then score every variant by
// matching it, and decision 5's class guard makes the matcher refuse a 'u' outright, so every variant
// would score 0, the winner would be whichever the tie-break reached first, and because a reading
// metric earns no ADR-0020 close-out the sweep would end with nothing acquired at the winning CE.
// Also fails under a per-MODE override, which the converse half below catches.
/////////////////////////////////////////////////////////////
START_SECTION(unassigned_ms3_sweep_is_scored_by_remaining_precursor)
{
  // MS3 exploration ON with the READING metric configured, so the override has something to beat.
  // No `overrides` map anywhere: that is the OTHER reason a production scan is dispatched, and
  // leaving it empty is what makes the close-out assertion below attributable to the metric alone.
  const std::string cfg_json = R"JSON({
    "deconvolution": {"tol": [10, 10, 10]},
    "precursor_selection": {"rank_by": "qscore", "max_precursors": 3},
    "characterization": {
      "mode": "exhaustive",
      "protein_sequence": "PEPTIDEK",
      "max_targets": 3,
      "exploration": {"metric": "fragment_count", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0}
    },
    "ms_settings": {
      "ms1": {"analyzer": "Orbitrap", "resolution": 120000},
      "ms2": {"analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29},
      "ms3": {"analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25}
    }
  })JSON";

  Config cfg{cfg_json};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  TEST_EQUAL(static_cast<int>(cfg.level(3).exploration), static_cast<int>(ExplorationMetric::FragmentCount))

  ScanCommand ms2_ctx = makeMs2Ctx();
  PeakGroup unassigned_pg = makeSyntheticPeakGroup(MASS_U1, 2);

  // 'u'/0 is exactly what planExhaustive_ labels a mass that matched no theoretical fragment.
  std::vector<ScanCommand> cmds = exploration.initiate(3, unassigned_pg, 2, queue, nullptr, &ms2_ctx, 'u', 0);
  TEST_EQUAL(static_cast<int>(cmds.size()), 4)   // CE-0 baseline + CE 20 / 25 / 30
  ABORT_IF(cmds.size() != 4)

  auto grp = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(static_cast<int>(grp.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))
  TEST_EQUAL(grp.msn_level, 3)

  // The sentinel drives the in-engine decision and nothing else: with ion_index 0 every variant takes
  // buildMS3's no-ion branch, so no 'u' appears after the 3-char tracking id (ADR-0023 5b / D-f). The
  // id itself is base-94 and may legitimately contain a 'u', hence the substr.
  for (const ScanCommand& c : cmds)
  {
    const std::string desc(c.scan_description);
    TEST_TRUE(desc.size() > 4)
    TEST_EQUAL(desc.substr(3).find('u') == std::string::npos, true)
  }

  // Drive the sweep with raw arrays. RemainingPrecursor scores from isolation-window intensity alone,
  // so one peak inside the window and two far outside it make the in-window sum exactly controllable.
  const double in_mz = grp.precursor_mz;                 // window is in_mz +- isolation_width/2 (>= 1.0 Th)
  std::vector<double> mzs{in_mz - 50.0, in_mz, in_mz + 50.0};

  std::vector<double> baseline_ints{300.0, 1000.0, 700.0};
  auto info_base = exploration.feedResult(cmds[0].scan_id, mzs.data(), baseline_ints.data(),
                                          static_cast<int>(mzs.size()), 0.5, queue);
  TEST_REAL_SIMILAR(info_base.metric.score, 0.0)   // the baseline never competes

  auto grp_after_baseline = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(grp_after_baseline.has_baseline, true)
  TEST_REAL_SIMILAR(grp_after_baseline.baseline_intensity, 1000.0)
  TEST_EQUAL(grp_after_baseline.baseline_failed, false)

  // score = 1 - |ratio - remaining_precursor_target|, target = 0.1 (the level default).
  // ratios 0.9 / 0.4 / 0.1 -> scores 0.2 / 0.7 / 1.0, so the WINNER is the last-fed variant rather
  // than the first -- a tie-break-order winner would land on index 1 instead.
  std::vector<double> ce20_ints{300.0, 900.0, 700.0};
  std::vector<double> ce25_ints{300.0, 400.0, 700.0};
  std::vector<double> ce30_ints{300.0, 100.0, 700.0};

  auto info20 = exploration.feedResult(cmds[1].scan_id, mzs.data(), ce20_ints.data(),
                                       static_cast<int>(mzs.size()), 1.0, queue);
  auto info25 = exploration.feedResult(cmds[2].scan_id, mzs.data(), ce25_ints.data(),
                                       static_cast<int>(mzs.size()), 1.5, queue);
  auto info30 = exploration.feedResult(cmds[3].scan_id, mzs.data(), ce30_ints.data(),
                                       static_cast<int>(mzs.size()), 2.0, queue);

  TEST_EQUAL(info30.metric.exploration_metric, static_cast<int>(ExplorationMetric::RemainingPrecursor))
  TEST_REAL_SIMILAR(info20.metric.score, 0.2)
  TEST_REAL_SIMILAR(info25.metric.score, 0.7)
  TEST_REAL_SIMILAR(info30.metric.score, 1.0)
  // Distinct and non-zero: the whole point. Under fragment_count these three are 0, 0, 0.
  TEST_EQUAL(info20.metric.score > 0.0 && info25.metric.score > info20.metric.score
             && info30.metric.score > info25.metric.score, true)
  TEST_REAL_SIMILAR(info20.metric.remaining_ratio, 0.9)
  TEST_REAL_SIMILAR(info30.metric.remaining_ratio, 0.1)

  // The winner is the best-scoring variant, and the group closes on the completing feed.
  TEST_EQUAL(info30.group.winner_tracking_id, ScanCommandQueue::encode(cmds[3].scan_id))
  TEST_EQUAL(exploration.activeGroupCount(), 0)

  // ADR-0020's close-out, earned with no edit at the dispatch site: RemainingPrecursor is a MEASURING
  // metric, so a completed MS3 sweep is followed by ONE production scan at the winning CE. This config
  // has no overrides, so this command can only have come from the metric being measuring -- which is
  // the second half of why decision 11 chooses this metric rather than, say, refusing to sweep.
  TEST_EQUAL(static_cast<int>(info30.commands.size()), 1)
  ABORT_IF(info30.commands.size() != 1)
  TEST_EQUAL(info30.commands[0].msn_level, 3)
  TEST_EQUAL(info30.commands[0].num_stages, 2)
  TEST_REAL_SIMILAR(info30.commands[0].stages[1].collision_energy, 30.0)   // the WINNING CE, not the config's 25
  TEST_EQUAL(std::string(info30.commands[0].scan_description).substr(3, 1), std::string("R"))  // production, not 'E'

  // --- the converse: a MAPPED target in the SAME exhaustive run keeps the configured metric --------
  // It has the ion frame a reading metric needs, so nothing is wrong with fragment_count for it. A
  // per-mode override would drag this group along with the unassigned one and silently stop matching
  // fragments for targets that can be matched.
  PeakGroup mapped_pg = makeSyntheticPeakGroup(MASS_B6, 2);
  std::vector<ScanCommand> mapped_cmds = exploration.initiate(3, mapped_pg, 2, queue, nullptr, &ms2_ctx, 'b', 6);
  TEST_EQUAL(static_cast<int>(mapped_cmds.size()), 4)
  ABORT_IF(mapped_cmds.empty())

  auto mapped_grp = ExplorationTestAccess::group(exploration, 2);   // group ids are monotone; 1 has closed
  TEST_EQUAL(static_cast<int>(mapped_grp.exploration_metric), static_cast<int>(ExplorationMetric::FragmentCount))
  TEST_EQUAL(mapped_grp.fragment_ion_type, 'b')
  TEST_EQUAL(mapped_grp.fragment_ion_index, 6)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 9. THE REGRESSION GUARD FOR FOUR EXISTING GOLDEN MODES: exploration_hcd, exploration_etd,
//    exploration_followup, exploration_multiplexed.
//
// Exploration::initiate's ion_type parameter defaults to '\0', which is NOT a known ion class, and the
// MS2 call site (FLASHIda.cpp:200) passes nothing. So the class test alone reads "unknown" on EVERY
// MS2 exploration group; only the msn_level >= 3 gate keeps decision 11 off them.
//
// FAILS UNDER: dropping that gate. Every MS2 group would come out RemainingPrecursor, which changes
// the scores, the winner and the logged exploration_metric column -- with nothing to compensate,
// because a measuring metric at MS2 is not a close-out case. The four goldens above move.
/////////////////////////////////////////////////////////////
START_SECTION(ms2_exploration_metric_is_never_overridden)
{
  // Both configurable metrics, because the guard must be metric-agnostic: the four golden modes do not
  // all sweep on the same one, and an override that only spared mass_count would still move the rest.
  auto ms2ExplorationConfig = [](const std::string& metric) {
    return R"JSON({
      "deconvolution": {"tol": [10, 10, 10]},
      "precursor_selection": {
        "rank_by": "qscore", "max_precursors": 3,
        "exploration": {"metric": ")JSON" + metric + R"JSON(", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0}
      },
      "characterization": {"mode": "off", "protein_sequence": "PEPTIDEK"},
      "ms_settings": {
        "ms1": {"analyzer": "Orbitrap", "resolution": 120000},
        "ms2": {"analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29},
        "ms3": {"analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25}
      }
    })JSON";
  };

  const std::vector<ExplorationMetric> metrics{ExplorationMetric::MassCount, ExplorationMetric::FragmentCount};
  const std::vector<std::string> metric_names{"mass_count", "fragment_count"};

  for (size_t mi = 0; mi < metrics.size(); ++mi)
  {
    Config cfg{ms2ExplorationConfig(metric_names[mi])};
    ScanCommandQueue queue(cfg);
    FragmentAnalysis fragments(cfg);
    Exploration exploration(cfg, fragments);

    TEST_EQUAL(static_cast<int>(cfg.level(2).exploration), static_cast<int>(metrics[mi]))

    PeakGroup pg = makeSyntheticPeakGroup(MASS_U1, 2);

    // No ion_type argument, exactly as the MS2 call site FLASHIda.cpp:200 leaves it -- so the group is
    // formed with ion_type at its '\0' default, the shape the four goldens were captured with.
    std::vector<ScanCommand> cmds = exploration.initiate(2, pg, 2, queue);
    TEST_EQUAL(static_cast<int>(cmds.size()), 4)
    auto grp = ExplorationTestAccess::group(exploration, 1);
    TEST_EQUAL(grp.msn_level, 2)
    TEST_EQUAL(static_cast<int>(grp.exploration_metric), static_cast<int>(metrics[mi]))

    // And passing the default sentinel explicitly changes nothing -- the gate is on the LEVEL, so it
    // cannot be defeated by an unknown class arriving at MS2 through any route.
    std::vector<ScanCommand> explicit_cmds = exploration.initiate(2, pg, 2, queue, nullptr, nullptr, '\0', 0);
    TEST_EQUAL(static_cast<int>(explicit_cmds.size()), 4)
    auto explicit_grp = ExplorationTestAccess::group(exploration, 2);
    TEST_EQUAL(static_cast<int>(explicit_grp.exploration_metric), static_cast<int>(metrics[mi]))
  }
}
END_SECTION

// The case the other three metric sections leave open, and the one that broke CI: 'u' at MS3 is
// pinned (-> RemainingPrecursor), 'b' at MS3 is pinned (-> the configured metric), and '\0' at MS2 is
// pinned by the level gate -- but NOTHING pinned '\0' at MS3.
//
// It is not a hypothetical hole. The override first shipped as `!isKnownIonClass(ion_type)`, which is
// true for '\0', so every MS3 group formed without an ion identity flipped to RemainingPrecursor.
// FLASHIda_exploration_test builds exactly such groups; with its harness feeding no raw arrays the
// baseline window measured 0, the group aborted, and the failure surfaced once as a score of 0 and
// once as a SegFault. '\0' means "no ion identity RECORDED" -- not "an identity the matcher refuses"
// -- and only the latter may change the metric.
START_SECTION(ms3_unrecorded_ion_type_keeps_the_configured_metric)
{
  const std::string cfg_json = R"JSON({
    "deconvolution": {"tol": [10, 10, 10]},
    "precursor_selection": {"rank_by": "qscore", "max_precursors": 3},
    "characterization": {
      "mode": "exhaustive",
      "protein_sequence": "PEPTIDEK",
      "max_targets": 3,
      "exploration": {"metric": "fragment_count", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0}
    },
    "ms_settings": {
      "ms1": {"analyzer": "Orbitrap", "resolution": 120000},
      "ms2": {"analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29},
      "ms3": {"analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25}
    }
  })JSON";

  Config cfg{cfg_json};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  TEST_EQUAL(static_cast<int>(cfg.level(3).exploration), static_cast<int>(ExplorationMetric::FragmentCount))

  ScanCommand ms2_ctx = makeMs2Ctx();
  PeakGroup pg = makeSyntheticPeakGroup(MASS_U1, 2);

  // Level 3, no ion identity recorded -- the initiate() default, reached explicitly here.
  std::vector<ScanCommand> cmds = exploration.initiate(3, pg, 2, queue, nullptr, &ms2_ctx, '\0', 0);
  TEST_EQUAL(static_cast<int>(cmds.size()), 4)   // CE-0 baseline + CE 20 / 25 / 30
  ABORT_IF(cmds.size() != 4)

  auto grp = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(grp.msn_level, 3)
  // The CONFIGURED metric survives. Under the shipped-and-reverted gate this read RemainingPrecursor.
  TEST_EQUAL(static_cast<int>(grp.exploration_metric), static_cast<int>(ExplorationMetric::FragmentCount))
}
END_SECTION

END_TEST
