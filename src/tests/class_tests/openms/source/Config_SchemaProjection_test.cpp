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
// The fourth mode (ADR-0023). It needed no edit in applyCharacterizationMode_ -- levels 2 and 3
// derive from `mode != Off`, an inequality rather than an enumeration of on-values -- and this
// section is what pins that. If anyone rewrites the projection as a list of on-values, level 1 or
// the level-2/3 pair silently drops out for exhaustive, and the failure is an empty run rather than
// a wrong value.
START_SECTION(exhaustive_projects_every_level)
{
  const std::string ms3 =
    R"(, "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25 })";

  {
    Config cfg(cfgJson("",
      R"("characterization": { "mode": "exhaustive", "protein_sequence": "PEPTIDEK",
                               "max_targets": 3, "min_target_mass": 700 },)", ms3));

    // Level 1 MUST be assigned, and from rank_by rather than from the mode. MSLevelConfig::selection
    // defaults to None, and an unassigned level 1 makes FLASHIda.cpp short-circuit EVERY MS1 -- the
    // instrument acquires nothing at all, silently, with no wrong value anywhere to notice.
    TEST_EQUAL((int)cfg.level(1).selection, (int)SelectionMetric::QScore)
    TEST_EQUAL(cfg.level(2).selection != SelectionMetric::None, true)
    TEST_EQUAL(cfg.level(3).selection != SelectionMetric::None, true)

    // Authored in characterization, read off level 2. The level's own default is 10, so this
    // distinguishes "projected" from "left at the default".
    TEST_EQUAL(cfg.level(2).max_targets, 3)

    // The mode carries its OWN objective. objective defaults to Ambiguity, so a mode that assigned
    // only `mode` would be byte-identical to "ambiguity" everywhere the engine actually branches.
    TEST_EQUAL((int)cfg.characterization().objective, (int)CharacterizationObjective::Exhaustive)
    TEST_REAL_SIMILAR(cfg.characterization().min_target_mass, 700.0)
  }

  // A fourth legal value did not loosen the strictness: unknown and wrong-case mode values still
  // throw, and so does an unknown key beside the new one.
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson("", R"("characterization": { "mode": "Exhaustive", "protein_sequence": "PEPTIDEK" },)", ms3)))
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson("", R"("characterization": { "mode": "exhaustve", "protein_sequence": "PEPTIDEK" },)", ms3)))
  TEST_EXCEPTION(std::invalid_argument,
    Config(cfgJson("", R"("characterization": { "mode": "exhaustive", "protein_sequence": "P", "min_targt_mass": 700 },)", ms3)))

  // min_target_mass is a plain optional key, not an exhaustive-only one: it parses under any mode
  // and simply goes unread. Gating it on the mode would make toggling the mode off a two-key edit.
  {
    Config cfg(cfgJson("",
      R"("characterization": { "mode": "ambiguity", "protein_sequence": "PEPTIDEK", "min_target_mass": 900 },)", ms3));
    TEST_REAL_SIMILAR(cfg.characterization().min_target_mass, 900.0)
    TEST_EQUAL((int)cfg.characterization().objective, (int)CharacterizationObjective::Ambiguity)
  }

  // Absent => 0, i.e. off. This is what keeps every existing config's pool unchanged.
  {
    Config cfg(cfgJson());
    TEST_REAL_SIMILAR(cfg.characterization().min_target_mass, 0.0)
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

// ---------------------------------------------------------------------------------------------
// The two charge-acquisition keys (ADR-0016). Both take the same three values through one shared
// parse, and both default to Single so an existing config is unaffected.
START_SECTION(charge_acquisition_modes_parse)
{
  // Absent => Single at both levels. This is what keeps the feature byte-identical when unused.
  {
    Config cfg(cfgJson());
    TEST_EQUAL((int)cfg.targeting().precursor_charges, (int)ChargeAcquisitionMode::Single)
    TEST_EQUAL((int)cfg.characterization().fragment_charges, (int)ChargeAcquisitionMode::Single)
  }

  {
    Config cfg(cfgJson(R"(, "precursor_charges": "separate")"));
    TEST_EQUAL((int)cfg.targeting().precursor_charges, (int)ChargeAcquisitionMode::Separate)
  }
  {
    Config cfg(cfgJson(R"(, "precursor_charges": "multiplexed")"));
    TEST_EQUAL((int)cfg.targeting().precursor_charges, (int)ChargeAcquisitionMode::Multiplexed)
  }
  {
    Config cfg(cfgJson("",
        R"("characterization": { "mode": "ambiguity", "protein_sequence": "PEPTIDEK", "fragment_charges": "multiplexed" },)",
        R"(, "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25 })"));
    TEST_EQUAL((int)cfg.characterization().fragment_charges, (int)ChargeAcquisitionMode::Multiplexed)
  }
}
END_SECTION

// An unknown VALUE throws rather than defaulting. The whole reason characterization.mode was made to
// throw is that its predecessor silently mapped anything-but-"coverage" to ambiguity, so a typo'd
// mode meant the opposite of what it said. Same failure shape here: a typo'd "Multiplexed" that
// quietly meant "single" would look like the feature simply not working.
START_SECTION(charge_acquisition_unknown_value_throws)
{
  TEST_EXCEPTION(std::invalid_argument, Config(cfgJson(R"(, "precursor_charges": "Multiplexed")")))
  TEST_EXCEPTION(std::invalid_argument, Config(cfgJson(R"(, "precursor_charges": "all")")))
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("",
          R"("characterization": { "mode": "ambiguity", "protein_sequence": "PEPTIDEK", "fragment_charges": "Separate" },)",
          R"(, "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25 })")))
}
END_SECTION

// The retired bool gets a migration error naming its replacement, not the generic unknown-key
// message -- checked before the allowlist precisely so the message can be specific.
START_SECTION(ms3_all_charges_is_a_migration_error)
{
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("",
          R"("characterization": { "mode": "ambiguity", "protein_sequence": "PEPTIDEK", "ms3_all_charges": true },)",
          R"(, "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25 })")))
}
END_SECTION

// Same treatment for the retired developer flag (ADR-0021). It earns a specific message rather than
// the generic unknown-key one for a reason peculiar to this key: its documented job was exclusion
// KEYING, but it was also the only thing that made precursor_charges: "separate" fan out. A reader
// told merely "unknown key" would remove it and silently lose the multi-charge acquisition they were
// relying on, with nothing pointing at "separate" as the replacement. Both values must throw --
// `false` was the default, so a config carrying it is the common case, not the exotic one.
START_SECTION(charge_based_exclusion_is_a_migration_error)
{
  TEST_EXCEPTION(std::invalid_argument, Config(cfgJson(R"(, "charge_based_exclusion": true)")))
  TEST_EXCEPTION(std::invalid_argument, Config(cfgJson(R"(, "charge_based_exclusion": false)")))
}
END_SECTION

// Drift guard, C++ half (ADR-0030). ScanFactory.NeedsReactionTime is the other half, pinned by
// ScanFactoryTests.NeedsReactionTime_MatchesEngineActivationSet.
//
// These predicates now answer TWO questions that used to be answered separately, which is what makes
// pinning them worthwhile. In Config they said which parameters an activation must arrive with; in
// Exploration::buildVariants_ they now also decide which axis of a sweep varies and which axis that
// activation's baseline zeroes; and in C# the reaction-time half decides whether the ReactionTime key
// reaches the instrument at all. A set that drifts therefore changes what is acquired and what is
// commanded, without changing anything the engine logs.
//
// Exact sets, not spot checks: an over-broad predicate is the failure that matters. Adding "CID" to
// needsReactionTime would start emitting a reaction time on every CID scan, and no existing assertion
// anywhere would notice.
START_SECTION(activation_coupling_predicates_are_the_declared_set)
{
  TEST_EQUAL(needsReactionTime("ETD"), true)
  TEST_EQUAL(needsReactionTime("EThcD"), true)
  TEST_EQUAL(needsReactionTime("HCD"), false)
  TEST_EQUAL(needsReactionTime("CID"), false)
  TEST_EQUAL(needsReactionTime("UVPD"), false)
  TEST_EQUAL(needsReactionTime(""), false)
  TEST_EQUAL(needsReactionTime("etd"), false)     // ordinal, case-sensitive -- see the C# mirror

  TEST_EQUAL(needsCollisionEnergy("HCD"), true)
  TEST_EQUAL(needsCollisionEnergy("CID"), true)
  TEST_EQUAL(needsCollisionEnergy("EThcD"), true)
  TEST_EQUAL(needsCollisionEnergy("ETD"), false)
  TEST_EQUAL(needsCollisionEnergy("UVPD"), false)
  TEST_EQUAL(needsCollisionEnergy(""), false)
  TEST_EQUAL(needsCollisionEnergy("hcd"), false)

  // EThcD is the only activation in both sets, and that is load-bearing: it is what makes
  // buildVariants_ take the cross-product branch and turn off BOTH axes for an EThcD baseline --
  // CE at 0, reaction time at MIN_REACTION_TIME_MS, since the two axes turn off at different values.
  TEST_EQUAL(needsCollisionEnergy("EThcD") && needsReactionTime("EThcD"), true)
}
END_SECTION

// An authored scan config pairing an activation with a zero on its coupled axis LOADS (ADR-0030).
// This used to throw. The guard's stated purpose was to stop ScanFactory silently dropping a zero
// reaction time so the instrument fell back to its own method default; that key is now gated on the
// activation rather than the value, so the value reaches the instrument and the silence is gone.
// Both branches were dropped together, so ETD and HCD stay symmetric.
// The sweep GRID is authored, not floored (ADR-0029): Exploration raises only its own synthesized
// baseline to MIN_REACTION_TIME_MS. So a grid starting below the floor would put a scan the
// instrument rejects into every sweep, AND -- because the baseline would no longer coincide with the
// grid's first point -- resurrect the duplicate reference scan the suppression rule removes. Both
// failures are visible only on the hardware, so this rejects the config at load instead.
//
// Both arms matter. The ETD arm is the guard; the HCD arm is what stops it being written over-broad,
// since reaction_time_min defaults to 0 and every CE-only config leaves it there.
START_SECTION(etd_sweep_rejects_sub_floor_reaction_time_min)
{
  auto sweep = [](const std::string& acts, const std::string& rt_min) {
    return R"(, "exploration": { "metric": "mass_count", "ce_min": 15, "ce_max": 50, "ce_step": 5,
                                 "reaction_time_min": )" + rt_min + R"(, "reaction_time_max": 20,
                                 "reaction_time_step": 5, "activations": [)" + acts + R"(] })";
  };

  // ETD swept + a sub-floor reaction_time_min -> rejected.
  TEST_EXCEPTION(std::invalid_argument, Config(cfgJson(sweep(R"("ETD")", "0"))))
  TEST_EXCEPTION(std::invalid_argument, Config(cfgJson(sweep(R"("ETD")", "0.01"))))
  // EThcD sweeps BOTH axes, so it needs the reaction time too.
  TEST_EXCEPTION(std::invalid_argument, Config(cfgJson(sweep(R"("EThcD")", "0"))))
  // One offending activation in a list is enough.
  TEST_EXCEPTION(std::invalid_argument, Config(cfgJson(sweep(R"("HCD", "ETD")", "0"))))

  // At the floor exactly, and above it -> accepted.
  Config at_floor{cfgJson(sweep(R"("ETD")", "0.03"))};
  TEST_REAL_SIMILAR(at_floor.level(2).rt_min, MIN_REACTION_TIME_MS)
  Config above{cfgJson(sweep(R"("ETD")", "5"))};
  TEST_REAL_SIMILAR(above.level(2).rt_min, 5.0)

  // NOT over-broad: a CE-only sweep leaves reaction_time_min at its 0 default and must still load.
  // Getting this wrong would reject every committed exploration config except the ETD one.
  Config hcd_only{cfgJson(sweep(R"("HCD")", "0"))};
  TEST_REAL_SIMILAR(hcd_only.level(2).rt_min, 0.0)
  Config cid_only{cfgJson(sweep(R"("CID")", "0"))};
  TEST_REAL_SIMILAR(cid_only.level(2).rt_min, 0.0)
}
END_SECTION

START_SECTION(zero_on_a_coupled_axis_is_accepted)
{
  // ms_settings.ms2 itself, NOT an additional_ms2 block: the check only ever ran over the dispatch
  // ROSTER (an unreferenced definition never fires), so probing it through an unreferenced block
  // would pass whether or not the guard is still there.
  auto ms2Json = [](const std::string& scan_body) {
    return std::string(R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": "qscore", "max_precursors": 5 },
      "characterization": { "mode": "off" },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
        "ms2": )") + scan_body + R"(
      }
    })";
  };

  // ETD at reaction_time 0 loads. This threw before ADR-0030.
  //
  // KNOWN, ACCEPTED GAP: unlike the exploration sweep, an AUTHORED scan config is not floored to
  // MIN_REACTION_TIME_MS and is not rejected for sitting below it, so this config commands a
  // reaction time the instrument refuses -- on every MS2 of the run. It fails loudly at the device
  // rather than silently inheriting the method default, which is the property ADR-0030 bought, but
  // it is not caught at load. Deliberate: only the sweep path, where a zero floor is the natural way
  // to ask for "no reaction", is guarded (see etd_sweep_rejects_sub_floor_reaction_time_min).
  Config etd{ms2Json(R"({ "analyzer": "Orbitrap", "activation": "ETD", "reaction_time": 0 })")};
  TEST_EQUAL(etd.level(2).scans.empty(), false)
  TEST_REAL_SIMILAR(etd.level(2).scans[0].reaction_time, 0.0)

  // HCD at collision_energy 0 -- the symmetric case, dropped in the same change so the two coupled
  // axes keep the same rule. It is also the CE an exploration baseline has always commanded.
  Config hcd{ms2Json(R"({ "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 0 })")};
  TEST_EQUAL(hcd.level(2).scans.empty(), false)
  TEST_EQUAL(hcd.level(2).scans[0].collision_energy, 0)
}
END_SECTION

// ---------------------------------------------------------------------------------------------
// Quantification (ADR-0038). These live here rather than in ConfigSchemaParity_test because the
// generated reference fixture keeps quantification DISABLED -- enabling it would invert the level-2
// roster, so that test's 17-key ms_settings.ms2 block (which reads level(2).scans[0]) would compare
// the wrong scan. The inversion is behaviour, so it is pinned here.

// cfgJson splices its `characterization` argument verbatim at the top level, so it doubles as the
// slot for a quantification section; `ms_settings_extra` is spliced inside ms_settings.
#define QUANT_ON(body) R"("characterization": { "mode": "off" }, "quantification": )" body R"(,)"
#define MS2_QUANT R"(, "ms2_quant": { "analyzer": "Orbitrap", "activation": "ETD", "reaction_time": 12 })"

START_SECTION(quantification_inverts_the_level2_roster)
{
  // THE structural claim of ADR-0038. With quantification on, the QUANTIFICATION scan takes the
  // roster's primary slot -- it is the screen, so it must be acquired to decide anything -- and
  // ms_settings.ms2 becomes the identification scan a differential verdict buys, held on the quant
  // config rather than dispatched. Without this, ms2 would fire unconditionally and the screen
  // would never be acquired at all.
  Config on(cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex",
        "conditions": [ { "name": "a", "channels": ["126","127","128"] },
                        { "name": "b", "channels": ["129","130","131"] } ] })"), MS2_QUANT));
  TEST_EQUAL(on.level(2).scans.size(), 1)
  // scans[0] is ms2_quant (ETD), NOT ms_settings.ms2 (HCD). Asserting the activation rather than
  // "is non-empty" is the point: both configs are present, so only the value distinguishes them.
  TEST_STRING_EQUAL(on.level(2).scans[0].activation, "ETD")
  TEST_REAL_SIMILAR(on.level(2).scans[0].reaction_time, 12.0)
  TEST_STRING_EQUAL(on.quantification().identification_scan.activation, "HCD")
  TEST_EQUAL(on.quantification().identification_scan.collision_energy, 29)

  // Quantification OFF puts ms_settings.ms2 back on the roster -- i.e. every other mode is
  // untouched, which is what keeps 'R' meaning "identification MS2" everywhere.
  Config off(cfgJson("", R"("characterization": { "mode": "off" },)", MS2_QUANT));
  TEST_EQUAL(off.level(2).scans.size(), 1)
  TEST_STRING_EQUAL(off.level(2).scans[0].activation, "HCD")
  // ...and ms2_quant is still PARSED when off, so the schema fixture can see its keys.
  TEST_STRING_EQUAL(off.quantification().quantification_scan.activation, "ETD")
}
END_SECTION

START_SECTION(quantification_conditions_resolve_by_channel_name)
{
  // Channel NAMES become ordinals into the scheme's getChannelInformation() at load, so nothing
  // downstream needs the name->m/z table and a typo cannot silently read the wrong intensity.
  Config cfg(cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt10plex",
        "conditions": [ { "name": "a", "channels": ["126","127C"] },
                        { "name": "b", "channels": ["130N","131"] } ] })"), MS2_QUANT));
  const auto& q = cfg.quantification();
  TEST_EQUAL(q.conditions.size(), 2)
  TEST_EQUAL(q.conditions[0].channels.size(), 2)
  // tmt10plex order: 126,127N,127C,128N,128C,129N,129C,130N,130C,131
  TEST_EQUAL(q.conditions[0].channels[0], 0)  // 126
  TEST_EQUAL(q.conditions[0].channels[1], 2)  // 127C -- NOT 127N, which is the 6.32 mDa neighbour
  TEST_EQUAL(q.conditions[1].channels[0], 7)  // 130N
  TEST_EQUAL(q.conditions[1].channels[1], 9)  // 131

  // A name that is not a channel of the SELECTED scheme is a load error. "127C" exists in
  // tmt10plex and not in tmt6plex, so this also proves the check is scheme-aware rather than a
  // static list.
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex",
        "conditions": [ { "name": "a", "channels": ["127C"] },
                        { "name": "b", "channels": ["131"] } ] })"), MS2_QUANT)))
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "not_a_scheme",
        "conditions": [ { "name": "a", "channels": ["126"] },
                        { "name": "b", "channels": ["131"] } ] })"), MS2_QUANT)))
}
END_SECTION

START_SECTION(quantification_structural_rejections)
{
  const std::string conds = R"("conditions": [ { "name": "a", "channels": ["126","127","128"] },
                                               { "name": "b", "channels": ["129","130","131"] } ])";

  // 1. enabled with no ms2_quant: nothing is rostered, ms_settings.ms2 stays unconditional, and the
  //    run is plain DDA that silently quantifies nothing. Note the ms_settings_extra is OMITTED.
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex", )" + conds + R"( })"), "")))

  // 2. conditions absent, and not-exactly-two. There is deliberately no fallback behind these: the
  //    positional 3-vs-3 split they replace was correct only for six-plex.
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex" })"), MS2_QUANT)))
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex",
        "conditions": [ { "name": "a", "channels": ["126"] } ] })"), MS2_QUANT)))
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex",
        "conditions": [ { "name": "a", "channels": ["126"] }, { "name": "b", "channels": ["129"] },
                        { "name": "c", "channels": ["131"] } ] })"), MS2_QUANT)))

  // 3. quantification + level-2 exploration. Exploration replaces the roster with CE-sweep variants,
  //    so the quantification scan is never dispatched and every variant is labelled 'E' -- never
  //    measured, never buying anything. The feature dies completely, and the one guard that might
  //    have caught it (exactly-one-scan-config at a swept level) is SATISFIED, because the inverted
  //    roster has exactly one entry. That is why this needs its own rejection.
  TEST_EXCEPTION(std::invalid_argument,
      // `metric` is what actually switches exploration on -- it defaults to "none", so an
      // exploration block without it leaves hasExploration(2) false and would prove nothing here.
      Config(cfgJson(R"(, "exploration": { "metric": "fragment_count", "activations": ["HCD"],
                                           "ce_min": 20, "ce_max": 30, "ce_step": 5 })",
                     QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex", )" + conds + R"( })"), MS2_QUANT)))

  // Disabled + incomplete is LEGAL: 40 of the 41 committed configs carry a quantification block
  // with neither conditions nor ms2_quant, and they must keep loading untouched.
  Config ok(cfgJson("", QUANT_ON(R"({ "enabled": false })"), ""));
  TEST_EQUAL(ok.quantification().enabled, false)
  TEST_EQUAL(ok.level(2).scans.size(), 1)
}
END_SECTION

START_SECTION(quantification_explicit_nulls_are_treated_as_absent)
{
  // REGRESSION. ToCppJson uses the STOCK JavaScriptSerializer, which EMITS nulls -- so every
  // config that authors neither key sends `"conditions": null, "correction_matrix": null`. A bare
  // contains() check then reads "present but not an array" and threw for ALL 41 committed configs,
  // which is 68 C# failures and an engine that will not construct. `is_null()` is the idiom
  // resolveFollowUp_ already uses in this file for exactly this reason.
  //
  // Asserted with quantification DISABLED because that is the shape 40 of the 41 configs have.
  Config nulls(cfgJson("", QUANT_ON(R"({ "enabled": false, "conditions": null,
                                         "correction_matrix": null })"), ""));
  TEST_EQUAL(nulls.quantification().enabled, false)
  TEST_EQUAL(nulls.quantification().conditions.empty(), true)
  TEST_EQUAL(nulls.quantification().correction_matrix.empty(), true)
  // ...and ms_settings.ms2 is still rostered, i.e. such a config is an ordinary DDA run.
  TEST_EQUAL(nulls.level(2).scans.size(), 1)
  TEST_STRING_EQUAL(nulls.level(2).scans[0].activation, "HCD")

  // A null ms2_quant is likewise absent rather than an empty scan config.
  Config nulls2(cfgJson("", QUANT_ON(R"({ "enabled": false })"), R"(, "ms2_quant": null)"));
  TEST_EQUAL(nulls2.quantification().has_quant_scan, false)
  TEST_EQUAL(nulls2.level(2).scans.size(), 1)
}
END_SECTION

START_SECTION(quantification_retired_keys_are_migration_errors)
{
  // Both earn their own message rather than a bare "unknown key", because deleting the key is the
  // wrong fix for follow_up_scan -- the block it named still has to go somewhere (ms2_quant).
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": false, "follow_up_scan": "whatever" })"), "")))
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": false, "only_one_condition": true })"), "")))
}
END_SECTION

// ---------------------------------------------------------------------------------------------
// ADR-0039: the quantification OBJECTIVE. Until this ADR, what a verdict bought was a hardcoded
// `verdict == Differential` in FLASHIda.cpp with no config surface at all.
// ---------------------------------------------------------------------------------------------

// The two conditions every section below shares. Named a/b so `enriched_in` has something to
// resolve against, and so the "either" collision has a control.
// ONE LINE, no continuation: a raw string literal does not process `\`-newline, so a continuation
// would put a literal backslash inside the JSON.
#define QUANT_CONDS R"("conditions": [ { "name": "a", "channels": ["126","127","128"] }, { "name": "b", "channels": ["129","130","131"] } ])"

START_SECTION(quantification_identify_parses_and_rejects)
{
  auto withIdentify = [](const std::string& body) {
    return cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex", )"
                                + std::string(QUANT_CONDS) + R"(, )" + body + R"( })"), MS2_QUANT);
  };

  // All four values reach their enumerator. Asserting the enumerator rather than "it loaded" is the
  // point: a parse that accepted the string and assigned nothing would ship byte-identical to
  // `differential` -- exactly the inert-mode trap ADR-0023 D-a documents for characterization.mode.
  TEST_EQUAL(Config(withIdentify(R"("identify": "differential")")).quantification().identify
                 == QuantIdentify::Differential, true)
  TEST_EQUAL(Config(withIdentify(R"("identify": "quantified")")).quantification().identify
                 == QuantIdentify::Quantified, true)
  TEST_EQUAL(Config(withIdentify(R"("identify": "all")")).quantification().identify
                 == QuantIdentify::All, true)
  TEST_EQUAL(Config(withIdentify(R"("identify": "none")")).quantification().identify
                 == QuantIdentify::None, true)

  // Unauthored defaults to Differential -- this is what makes ADR-0039 byte-identical for the 41
  // committed configs, none of which authors the key.
  TEST_EQUAL(Config(withIdentify(R"("reporter_mz_tol": 0.002)")).quantification().identify
                 == QuantIdentify::Differential, true)

  // Explicit null is "unauthored", NOT a throw. ToCppJson uses the stock JavaScriptSerializer,
  // which emits nulls; a bare .value() here would throw a type_error for every config that leaves
  // the key alone, which is how `conditions` broke all 41 configs once.
  TEST_EQUAL(Config(withIdentify(R"("identify": null)")).quantification().identify
                 == QuantIdentify::Differential, true)

  // Hard-rejected, not defaulted. Case matters, whitespace matters, and a plausible-but-wrong word
  // must not silently select a different acquisition policy.
  TEST_EXCEPTION(std::invalid_argument, Config(withIdentify(R"("identify": "Differential")")))
  TEST_EXCEPTION(std::invalid_argument, Config(withIdentify(R"("identify": "differential ")")))
  TEST_EXCEPTION(std::invalid_argument, Config(withIdentify(R"("identify": "any")")))
  TEST_EXCEPTION(std::invalid_argument, Config(withIdentify(R"("identify": "")")))
  // "off" is the characterization.mode spelling, and the one an author is most likely to reach for
  // when they mean "none". It must not quietly disable the buy.
  TEST_EXCEPTION(std::invalid_argument, Config(withIdentify(R"("identify": "off")")))
}
END_SECTION

START_SECTION(quantification_enriched_in_resolves_by_condition_name)
{
  auto withEnriched = [](const std::string& body) {
    return cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex", )"
                                + std::string(QUANT_CONDS) + R"(, )" + body + R"( })"), MS2_QUANT);
  };

  // Named by CONDITION, resolved to an ordinal at load. Direction is never authored as up/down:
  // fold_change = mean(conditions[0]) / mean(conditions[1]), so "up" would mean enriched in
  // conditions[0] and would invert silently if the array were reordered.
  TEST_EQUAL(Config(withEnriched(R"("enriched_in": "a")")).quantification().enriched_in, 0)
  TEST_EQUAL(Config(withEnriched(R"("enriched_in": "b")")).quantification().enriched_in, 1)

  // -1 is "either direction" -- ADR-0038's symmetric test, and the default.
  TEST_EQUAL(Config(withEnriched(R"("enriched_in": "either")")).quantification().enriched_in, -1)
  TEST_EQUAL(Config(withEnriched(R"("reporter_mz_tol": 0.002)")).quantification().enriched_in, -1)
  TEST_EQUAL(Config(withEnriched(R"("enriched_in": null)")).quantification().enriched_in, -1)

  // An unknown name fails at LOAD. Silently ignoring it would run the experiment in whichever
  // direction happened to be first, with nothing anywhere saying so.
  TEST_EXCEPTION(std::invalid_argument, Config(withEnriched(R"("enriched_in": "treated")")))
  // Condition names are case-sensitive, like every other enum in this schema.
  TEST_EXCEPTION(std::invalid_argument, Config(withEnriched(R"("enriched_in": "A")")))

  // enriched_in authored with NO conditions at all is rejected rather than ignored. This is why the
  // resolution sits outside the conditions block: inside it, this config would load with
  // enriched_in silently at -1.
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": false, "enriched_in": "a" })"), MS2_QUANT)))

  // "either" is the sentinel, so a condition may not claim the name -- otherwise enriched_in:
  // "either" is ambiguous and the losing reading is silently the wrong direction.
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex",
        "conditions": [ { "name": "either", "channels": ["126","127","128"] },
                        { "name": "b", "channels": ["129","130","131"] } ] })"), MS2_QUANT)))
  // Rejected even with quantification OFF: the name is unusable either way, and a config that
  // loads while disabled and throws when enabled is the worse failure.
  TEST_EXCEPTION(std::invalid_argument,
      Config(cfgJson("", QUANT_ON(R"({ "enabled": false,
        "conditions": [ { "name": "a", "channels": ["126"] },
                        { "name": "either", "channels": ["131"] } ] })"), "")))
}
END_SECTION

START_SECTION(quantification_requires_the_identification_scan)
{
  // ADR-0039's fourth structural rejection, and it closes a LIVE gap rather than guarding a new
  // one. Config assigns quant_.identification_scan only `if (has_primary_ms2)`, so without
  // ms_settings.ms2 the scan a verdict buys is built from a DEFAULT-CONSTRUCTED ScanConfig. That
  // was latent while only a Differential verdict reached it; `identify: "all"` fires it on every
  // precursor, so the guard lands in the same change that widens the gap.
  //
  // cfgJson always emits ms_settings.ms2, so these configs are written out in full.
  auto noMs2 = [](const std::string& quant_body) {
    return std::string(R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": "qscore", "max_precursors": 5 },
      "characterization": { "mode": "off" },
      "quantification": )") + quant_body + R"(,
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
        "ms2_quant": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 30 }
      }
    })";
  };

  // The message is asserted, not just the exception type. Config ALREADY rejects a missing
  // ms_settings.ms2 via `level(1).selection != None && level(2).scans.empty()` -- but that guard
  // reads the ROSTER, and in a quant config the roster is non-empty because the 'Q' scan holds it.
  // So the pre-existing guard cannot fire here, and a bare TEST_EXCEPTION would pass on either
  // one: this whole section would stay green with the ADR-0039 throw deleted.
  auto throwMessage = [](const std::string& json) {
    try { Config c(json); } catch (const std::exception& e) { return std::string(e.what()); }
    return std::string("<no exception>");
  };
  auto isTheNewGuard = [](const std::string& msg) {
    // "is not set" (ADR-0039) vs the pre-existing "is not defined".
    return msg.find("ms_settings.ms2 is not set") != std::string::npos;
  };

  const std::string enabled_body = R"({ "enabled": true, "labelling": "tmt6plex", )"
                                   + std::string(QUANT_CONDS) + R"( })";
  TEST_EXCEPTION(std::invalid_argument, Config(noMs2(enabled_body)))
  TEST_EQUAL(isTheNewGuard(throwMessage(noMs2(enabled_body))), true)

  // Still required under identify: "none", where it is INERT. Required-but-inert rather than
  // optional is deliberate (ADR-0013's ms_settings.ms3 rule): it keeps all four identify values
  // interchangeable, so flipping the objective can never invalidate a config.
  const std::string none_body = R"({ "enabled": true, "labelling": "tmt6plex", "identify": "none", )"
                                + std::string(QUANT_CONDS) + R"( })";
  TEST_EXCEPTION(std::invalid_argument, Config(noMs2(none_body)))
  TEST_EQUAL(isTheNewGuard(throwMessage(noMs2(none_body))), true)

  // The two guards cover DISJOINT states, which is why both are needed. With quantification off the
  // roster is empty, so the pre-existing roster guard fires and the ADR-0039 one -- gated on
  // `enabled` alongside its three ADR-0038 siblings -- does not.
  TEST_EQUAL(isTheNewGuard(throwMessage(noMs2(R"({ "enabled": false })"))), false)
  TEST_EQUAL(throwMessage(noMs2(R"({ "enabled": false })")).find("is not defined")
                 != std::string::npos, true)

  // The control: the same config WITH ms_settings.ms2 loads, so the rejection is about the missing
  // key and not about anything else in the fixture.
  Config ok(cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex", )"
                                 + std::string(QUANT_CONDS) + R"( })"), MS2_QUANT));
  TEST_EQUAL(ok.quantification().has_identification_scan, true)
  TEST_STRING_EQUAL(ok.quantification().identification_scan.activation, "HCD")
}
END_SECTION

START_SECTION(quantification_inert_enriched_in_is_a_warning_not_a_throw)
{
  // enriched_in restricts only the `differential` objective. Under the other three it is inert --
  // and inert is announced with [CONFIG-WARN], NOT rejected.
  //
  // This is a deliberate decision, so it gets a test rather than being left to chance. A throw
  // would force a second edit every time `identify` is flipped and would invalidate a lab template
  // that sets enriched_in once and switches objectives between runs. The distinction from
  // only_one_condition, which ADR-0038 deleted, is that this key is LIVE under `differential`;
  // only_one_condition was unreachable in every possible config. It is ms_settings.ms3-under-
  // mode-off, which ADR-0013 explicitly permits.
  auto withBoth = [](const std::string& identify, const std::string& enriched) {
    return cfgJson("", QUANT_ON(R"({ "enabled": true, "labelling": "tmt6plex", )"
                                + std::string(QUANT_CONDS) + R"(, "identify": ")" + identify
                                + R"(", "enriched_in": ")" + enriched + R"(" })"), MS2_QUANT);
  };

  // Loads, and the resolved ordinal is PRESERVED rather than quietly reset -- so flipping identify
  // back to "differential" restores the intended direction with a one-word edit.
  Config all(withBoth("all", "b"));
  TEST_EQUAL(all.quantification().enriched_in, 1)
  TEST_EQUAL(all.quantification().identify == QuantIdentify::All, true)

  Config none(withBoth("none", "a"));
  TEST_EQUAL(none.quantification().enriched_in, 0)

  Config quantified(withBoth("quantified", "b"));
  TEST_EQUAL(quantified.quantification().enriched_in, 1)

  // The combination it is NOT inert in.
  Config diff(withBoth("differential", "b"));
  TEST_EQUAL(diff.quantification().enriched_in, 1)
  TEST_EQUAL(diff.quantification().identify == QuantIdentify::Differential, true)
}
END_SECTION

#undef QUANT_CONDS
#undef QUANT_ON
#undef MS2_QUANT

END_TEST
