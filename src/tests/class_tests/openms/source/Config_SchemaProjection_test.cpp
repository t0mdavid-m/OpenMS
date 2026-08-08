// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>

#include <string>

using namespace OpenMS;

// Guards the two-decision-section config schema: the projection that derives every level's
// selection state from characterization.mode, name resolution for ms_settings.additional_ms2, and
// the validation that replaced four different silent guesses.
//
// One binary with four sections rather than four binaries, deliberately: a C++ test runs in CI only
// if it is listed in BOTH the build --target list and the ctest -R alternation in
// .github/workflows/flashida-ci.yml. Four registrations are four chances to build a test that never
// executes.

namespace
{
  // Minimal loadable config. `extra` is spliced in at the top level so a section can add or
  // override whatever it is probing.
  std::string cfgJson(const std::string& precursor_extra = "",
                      const std::string& characterization = R"("characterization": { "mode": "off" },)",
                      const std::string& ms_settings_extra = "")
  {
    return std::string(R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": "qscore", "max_precursors": 5)") + precursor_extra + R"( },
      )" + characterization + R"(
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
        "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 })"
        + ms_settings_extra + R"(
      }
    })";
  }
}

START_TEST(Config_SchemaProjection, "$Id$")

// ---------------------------------------------------------------------------------------------
START_SECTION(projection_covers_all_three_levels)
{
  // THE regression this file exists for.
  //
  // MSLevelConfig::selection defaults to None (it defaulted to Intensity, which meant merely
  // DEFINING a scan config switched a level on). With that default, a projection that assigns only
  // levels 2 and 3 leaves level 1 at None -- and FLASHIda.cpp:168 then short-circuits before MS1
  // selection, so the instrument acquires NOTHING AT ALL, silently. No golden would show a wrong
  // value; the run is simply empty. Assert level 1 explicitly.
  {
    Config cfg(cfgJson());  // mode: off
    TEST_EQUAL((int)cfg.level(1).selection, (int)SelectionMetric::QScore)
    TEST_EQUAL((int)cfg.level(2).selection, (int)SelectionMetric::None)
    TEST_EQUAL((int)cfg.level(3).selection, (int)SelectionMetric::None)
    TEST_EQUAL(cfg.level(1).max_targets, 5)   // max_precursors == the MS2 count
  }

  // rank_by drives level 1 and is INDEPENDENT of mode: MS3 off must not disable MS1 selection.
  {
    Config cfg(cfgJson(R"(, "rank_by": "intensity")"));
    TEST_EQUAL((int)cfg.level(1).selection, (int)SelectionMetric::Intensity)
    TEST_EQUAL((int)cfg.level(2).selection, (int)SelectionMetric::None)
  }

  // mode != off turns BOTH downstream gates on together. MS3 requires level 2 AND level 3 non-None
  // (FLASHIda.cpp:366, Exploration.cpp:728 and :730); driving both from one enum is what makes
  // "MS3 on with MS2 off" unrepresentable instead of merely discouraged.
  {
    Config cfg(cfgJson("",
      R"("characterization": { "mode": "ambiguity", "protein_sequence": "PEPTIDEK", "max_targets": 7, "min_fragment_charge": 2 },)",
      R"(, "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25 })"));
    TEST_EQUAL((int)cfg.level(1).selection, (int)SelectionMetric::QScore)
    TEST_EQUAL(cfg.level(2).selection != SelectionMetric::None, true)
    TEST_EQUAL(cfg.level(3).selection != SelectionMetric::None, true)
    // Authored in characterization, READ off level 2 -- the budget and fragment-charge floor are
    // consumed at ProteoformTracker.cpp:354 and Exploration.cpp:800.
    TEST_EQUAL(cfg.level(2).max_targets, 7)
    TEST_EQUAL(cfg.level(2).min_charge, 2)
    TEST_EQUAL((int)cfg.characterization().objective, (int)CharacterizationObjective::Ambiguity)
  }

  // coverage derives the other objective from the same key.
  {
    Config cfg(cfgJson("",
      R"("characterization": { "mode": "coverage", "protein_sequence": "PEPTIDEK" },)",
      R"(, "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25 })"));
    TEST_EQUAL((int)cfg.characterization().objective, (int)CharacterizationObjective::Coverage)
  }

  // levels_ is always {1,2,3}. toleranceList() walks it POSITIONALLY to build the DoubleList that
  // constructs Deconvolution, so a level going missing would shift every tols_[ms_level-1] index.
  {
    Config cfg(cfgJson());
    TEST_EQUAL(cfg.levels().size(), 3)
  }
}
END_SECTION

// ---------------------------------------------------------------------------------------------
START_SECTION(scan_name_resolution)
{
  const std::string two_extra_scans =
    R"(, "additional_ms2": {
         "etd_extra": { "analyzer": "Orbitrap", "activation": "ETD", "reaction_time": 10 },
         "aaa_first": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 20 } })";

  // ORDER COMES FROM THE REFERENCE ARRAY, never from iterating the map. nlohmann's object_t is a
  // std::map, so map iteration would sort these alphabetically and put aaa_first before etd_extra
  // -- silently reordering MS2 dispatch, and with it scan_commands row order and child_ids.
  {
    Config cfg(cfgJson(R"(, "additional_scans": ["etd_extra", "aaa_first"])",
                       R"("characterization": { "mode": "off" },)", two_extra_scans));
    TEST_EQUAL(cfg.level(2).scans.size(), 3)             // ms2 + two extras
    TEST_EQUAL(cfg.level(2).scans[0].activation, "HCD")  // ms_settings.ms2 always first
    TEST_EQUAL(cfg.level(2).scans[1].activation, "ETD")  // array order, NOT alphabetical
    TEST_EQUAL(cfg.level(2).scans[2].collision_energy, 20)
  }

  // A definition that nobody references is NOT in the roster, so it never fires. This is the whole
  // mechanism keeping a follow-up-backing block from becoming an unconditional MS2.
  {
    Config cfg(cfgJson("", R"("characterization": { "mode": "off" },)", two_extra_scans));
    TEST_EQUAL(cfg.level(2).scans.size(), 1)
  }

  // Dangling reference: hard error, and the message lists what IS defined.
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson(R"(, "additional_scans": ["nope"])",
                   R"("characterization": { "mode": "off" },)", two_extra_scans)))

  // Duplicate reference.
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson(R"(, "additional_scans": ["etd_extra", "etd_extra"])",
                   R"("characterization": { "mode": "off" },)", two_extra_scans)))

  // Name grammar: outer map keys are user-authored so they cannot be allowlisted; they are
  // validated as snake_case identifiers instead.
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson("", R"("characterization": { "mode": "off" },)",
      R"(, "additional_ms2": { "ETD Extra": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 20 } })")))

  // Reserved words would make a reference ambiguous with a level name or with "no scans".
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson("", R"("characterization": { "mode": "off" },)",
      R"(, "additional_ms2": { "none": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 20 } })")))

  // Inner objects are still validated against the 17-key scan allowlist. Returning early for any
  // Dictionary (the old C#-side behaviour) would let this load clean and drop the key.
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson("", R"("characterization": { "mode": "off" },)",
      R"(, "additional_ms2": { "bad": { "IsolationMode": "Quadrupole" } })")))

  // DOUBLE DUTY. additional_ms2 is one flat namespace serving two mutually exclusive roles:
  // additional_scans dispatches unconditionally per precursor, follow_up_scan fires conditionally
  // off a returning MS2. A name in both does BOTH -- the same config acquired twice per precursor,
  // at two different priorities, with no diagnostic. Nothing else catches it: the reference
  // resolves, so the dangling-name check passes, and the block IS referenced, so the unreferenced
  // warning stays quiet. Note conditional_ms2 must be true or validate() throws for the other
  // reason first, which would make this section pass vacuously.
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson(R"(, "additional_scans": ["etd_extra"])",
                   R"("characterization": { "mode": "off" },
                      "conditional_ms2": true,
                      "tagging": { "follow_up_scan": "etd_extra" },)",
                   two_extra_scans)))

  // The same block used ONLY as a follow-up is fine -- that is the intended arrangement, and it is
  // what keeps a follow-up out of the unconditional roster.
  {
    Config cfg(cfgJson("", R"("characterization": { "mode": "off" },
                              "conditional_ms2": true,
                              "tagging": { "follow_up_scan": "etd_extra" },)",
                       two_extra_scans));
    TEST_EQUAL(cfg.level(2).scans.size(), 1)  // ms2 only; the follow-up is NOT in the roster
  }

  // Array form is a migration error, not a generic unknown-key one.
  TEST_EXCEPTION(std::invalid_argument,
    Config(std::string(R"({ "deconvolution": { "tol": [10,10,10] },
      "ms_settings": { "ms1": {}, "ms2": [ { "activation": "HCD", "collision_energy": 29 } ] } })")))

  // So is a surviving selection_strategy block.
  TEST_EXCEPTION(std::invalid_argument,
    Config(std::string(R"({ "deconvolution": { "tol": [10,10,10] },
      "selection_strategy": { "ms1": { "selection": "qscore" } } })")))
}
END_SECTION

// ---------------------------------------------------------------------------------------------
START_SECTION(enum_strictness)
{
  // Four parsers used to guess on an unrecognised value, and they guessed in DIFFERENT directions:
  //   selection  -> Intensity  (a typo ENABLED a level)
  //   metric     -> None       (a typo collapsed an N-variant sweep to one scan)
  //   objective  -> Ambiguity  (a typo silently picked an objective)
  // With mode carrying the on/off bit, a typo'd "Off" would silently enable MS3. All throw now.

  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson("", R"("characterization": { "mode": "Off" },)")))       // wrong case
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson("", R"("characterization": { "mode": "covrage" },)")))   // typo
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson(R"(, "rank_by": "QScore")")))                            // wrong case
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson(R"(, "targeting": "exclusion")")))                       // near-miss: it is exclusion_masses
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson(R"(, "exploration": { "metric": "fragmentcount" })")))

  // The legal values still load.
  {
    Config a(cfgJson(R"(, "targeting": "in_depth")"));
    TEST_EQUAL(a.targeting().mode, 2)          // 2 is IN-DEPTH per PrecursorSelection.cpp:138-141,
    Config b(cfgJson(R"(, "targeting": "exclusion_masses")"));
    TEST_EQUAL(b.targeting().mode, 3)          // ...and 3 is exclusion. The old doc comments had
  }                                            // these two the wrong way round in three files.

  // tolerance_ppm is first-class; leaving it in the overrides map is an explicit error rather than
  // the old silent drop (it used to be extracted and ERASED before Exploration.cpp:605 tested that
  // same map for emptiness to decide whether to acquire the production scan).
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson(R"(, "exploration": { "metric": "mass_count", "overrides": { "tolerance_ppm": "8" } })")))
}
END_SECTION

// ---------------------------------------------------------------------------------------------
START_SECTION(sweep_step_and_required_scan_validation)
{
  // A non-positive step is not a no-op: Exploration.cpp's
  //     for (ce = ce_min; ce <= ce_max + 1e-9; ce += ce_step)
  // never advances, so it spins forever INSIDE processScan -- on the C# ActionBlock thread, with
  // the instrument still waiting for commands. A hang, not an error. Validated on NEITHER side
  // before this change.
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson(R"(, "exploration": { "metric": "mass_count", "ce_min": 15, "ce_max": 50, "ce_step": 0 })")))
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson(R"(, "exploration": { "metric": "mass_count", "ce_min": 15, "ce_max": 50, "ce_step": -1 })")))
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson(R"(, "exploration": { "metric": "mass_count", "reaction_time_min": 1, "reaction_time_max": 9, "reaction_time_step": 0 })")))

  // A valid sweep loads.
  {
    Config cfg(cfgJson(R"(, "exploration": { "metric": "mass_count", "ce_min": 15, "ce_max": 50, "ce_step": 1 })"));
    TEST_EQUAL((int)cfg.level(2).exploration, (int)ExplorationMetric::MassCount)
    // Unset exploration tolerance falls back to the level's base tol, matching what the old
    // overrides path produced.
    TEST_REAL_SIMILAR(cfg.level(2).exploration_tolerance_ppm, 10.0)
  }

  // mode != off REQUIRES ms_settings.ms3. This is the converse of "mode: off does not FORBID an
  // ms3 block", and it is the direction that segfaults: Exploration::initiateNextLevel reads
  // next_cfg.scans[0] unguarded, so a reachable MS3 with no level-3 scan config is an OOB read.
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson("", R"("characterization": { "mode": "ambiguity", "protein_sequence": "PEPTIDEK" },)")))

  // mode != off also requires a protein sequence -- re-keyed onto mode. It used to fire off the
  // UPSTREAM gate, which is why 18 configs that ran no MS3 still carried a placeholder sequence.
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson("", R"("characterization": { "mode": "ambiguity" },)",
      R"(, "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25 })")))

  // ...and mode: off does NOT forbid either of them, so toggling MS3 off is a one-word edit.
  {
    Config cfg(cfgJson("",
      R"("characterization": { "mode": "off", "protein_sequence": "PEPTIDEK" },)",
      R"(, "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25 })"));
    TEST_EQUAL((int)cfg.level(3).selection, (int)SelectionMetric::None)
  }
}
END_SECTION

END_TEST
