// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>

#include <fstream>
#include <sstream>
#include <string>

using namespace OpenMS;

// ---------------------------------------------------------------------------
// Schema drift guard (C++ side). The shared reference fixture
// FlashIDA/test-data/config_schema_reference.json gives every bridge-schema key a UNIQUE
// sentinel value; this test proves the C++ Config reader binds every key to the right field.
// The C# side (Flash.Tests/ConfigSchemaParityTests) proves ToCppJson emits every reference key.
// Together with the config-schema-drift-reminder.sh hook, the two schemas cannot silently diverge.
//
// If you add/rename/move a config key you MUST add its sentinel to the fixture and its assertion
// here (and ToCppJson emission on the C# side).
// ---------------------------------------------------------------------------

START_TEST(ConfigSchemaParity, "$Id$")

std::string json;
{
  std::ifstream in("../../FlashIDA/test-data/config_schema_reference.json");
  std::stringstream ss;
  ss << in.rdbuf();
  json = ss.str();
}

START_SECTION(reference_fixture_is_readable)
{
  TEST_EQUAL(json.empty(), false)
  TEST_EQUAL(json[0], '{')
}
END_SECTION

START_SECTION(EveryKey_ParsesToSentinel)
{
  Config cfg(json);

  // --- deconvolution ---
  TEST_EQUAL(cfg.deconvolution().min_charge, 4)
  TEST_EQUAL(cfg.deconvolution().max_charge, 47)
  TEST_REAL_SIMILAR(cfg.deconvolution().min_mass, 511.0)
  TEST_REAL_SIMILAR(cfg.deconvolution().max_mass, 49001.0)
  TEST_REAL_SIMILAR(cfg.targeting().qscore_threshold, 0.11)
  TEST_REAL_SIMILAR(cfg.targeting().tqscore_threshold, 0.93)
  TEST_REAL_SIMILAR(cfg.targeting().tag_matching_tolerance_ppm, 12.0)  // tol[1]

  // --- precursor_selection (incl. the formerly-developer-nested flags) ---
  TEST_REAL_SIMILAR(cfg.targeting().rt_window, 181.0)
  TEST_EQUAL(cfg.targeting().mode, 1)
  TEST_EQUAL(cfg.targeting().consider_all_charges, true)
  TEST_EQUAL(cfg.targeting().charge_based_exclusion, true)
  TEST_EQUAL(cfg.targeting().hcd_energy, 27)
  TEST_EQUAL(cfg.targeting().strict_inclusion, true)
  TEST_REAL_SIMILAR(cfg.targeting().tie_threshold, 0.13)

  // --- flashtnt ---
  TEST_EQUAL(cfg.targeting().min_tag_length, 4)
  TEST_EQUAL(cfg.targeting().max_tag_length, 9)
  TEST_EQUAL(cfg.targeting().max_total_ptm_count, 5)
  TEST_REAL_SIMILAR(cfg.targeting().max_flanking_mass_diff, 41001.0)
  TEST_EQUAL(cfg.targeting().allow_gap, true)
  TEST_EQUAL(cfg.targeting().max_aa_in_gap, 3)
  TEST_EQUAL(cfg.targeting().fixed_mod.size(), 1)
  TEST_EQUAL(cfg.targeting().fixed_mod.empty() ? std::string("") : cfg.targeting().fixed_mod[0], std::string("Carbamidomethyl (C)"))
  TEST_EQUAL(cfg.targeting().max_blind_mod_count, 1)
  TEST_REAL_SIMILAR(cfg.targeting().max_mod_mass, 733.0)  // load-bearing: NOT the extender default 500

  // --- tagging.follow_up_scan ---
  TEST_EQUAL(cfg.targeting().tagging_follow_up_scan.analyzer, std::string("Orbitrap"))
  TEST_EQUAL(cfg.targeting().tagging_follow_up_scan.activation, std::string("ETD"))
  TEST_EQUAL(cfg.targeting().tagging_follow_up_scan.collision_energy, 24)
  TEST_EQUAL(cfg.targeting().tagging_follow_up_scan.resolution, 15002)

  // --- conditional_ms2 ---
  TEST_EQUAL(cfg.targeting().conditional_ms2_enabled, false)

  // --- quantification ---
  TEST_EQUAL(cfg.quantification().enabled, false)
  TEST_REAL_SIMILAR(cfg.quantification().reporter_mz_tol, 0.0031)
  TEST_REAL_SIMILAR(cfg.quantification().fold_change_threshold, 1.7)

  // --- faims (incl. the renamed cv_precursor_threshold + formerly-developer max_cv_skip) ---
  TEST_EQUAL(cfg.faims().cv_values.size(), 3)
  TEST_REAL_SIMILAR(cfg.faims().cv_values[0], -41.0)
  TEST_REAL_SIMILAR(cfg.faims().cv_values[1], -52.0)
  TEST_REAL_SIMILAR(cfg.faims().cv_values[2], -63.0)
  TEST_EQUAL(cfg.faims().max_cv_skip, 2)
  TEST_EQUAL(cfg.faims().precursor_threshold, 17)
  TEST_EQUAL(cfg.faims().enabled, true)  // derived: cv_values.size() > 1

  // --- ms_settings.ms1 (snake_case -> struct fields) ---
  const auto& ms1 = cfg.level(1).scans[0];
  TEST_EQUAL(ms1.analyzer, std::string("Orbitrap"))
  TEST_REAL_SIMILAR(ms1.first_mass, 501.0)
  TEST_REAL_SIMILAR(ms1.last_mass, 2001.0)
  TEST_EQUAL(ms1.resolution, 120001)
  TEST_EQUAL(ms1.agc_target, 800001)
  TEST_REAL_SIMILAR(ms1.max_it, 247.0)
  TEST_EQUAL(ms1.microscans, 2)
  TEST_REAL_SIMILAR(ms1.rf_lens, 31.0)
  TEST_REAL_SIMILAR(ms1.source_cid, 16.0)
  TEST_EQUAL(ms1.data_type, std::string("Centroid"))

  // --- ms_settings.ms2 / ms3 ---
  const auto& ms2 = cfg.level(2).scans[0];
  TEST_EQUAL(ms2.activation, std::string("HCD"))
  TEST_EQUAL(ms2.collision_energy, 29)
  TEST_EQUAL(ms2.resolution, 120002)
  TEST_EQUAL(ms2.agc_target, 500001)
  TEST_REAL_SIMILAR(ms2.first_mass, 101.0)
  TEST_REAL_SIMILAR(ms2.last_mass, 2002.0)
  TEST_EQUAL(ms2.microscans, 3)
  const auto& ms3 = cfg.level(3).scans[0];
  TEST_EQUAL(ms3.activation, std::string("CID"))
  TEST_EQUAL(ms3.collision_energy, 26)
  TEST_EQUAL(ms3.resolution, 240001)

  // --- scheduling (nested) ---
  TEST_EQUAL(cfg.scheduling().cycle_time_enabled, true)
  TEST_REAL_SIMILAR(cfg.scheduling().cycle_time_ms, 60001.0)
  TEST_EQUAL(cfg.scheduling().timeout_enabled, true)
  TEST_REAL_SIMILAR(cfg.scheduling().timeout_ms, 30001.0)
  TEST_EQUAL((int)cfg.scheduling().agc_interval_ms, 29000)  // 29 s -> ms

  // --- characterization ---
  TEST_EQUAL((int)cfg.characterization().objective, (int)CharacterizationObjective::Coverage)
  TEST_EQUAL(cfg.characterization().protein_sequence, std::string("MSENTINELPEPTIDESEQ"))
  TEST_EQUAL(cfg.characterization().ms3_all_charges, true)

  // --- selection_strategy ---
  TEST_EQUAL((int)cfg.level(1).selection, (int)SelectionMetric::QScore)
  TEST_EQUAL(cfg.level(1).max_targets, 3)
  TEST_EQUAL(cfg.level(1).min_charge, 2)
  TEST_EQUAL((int)cfg.level(2).selection, (int)SelectionMetric::Intensity)
  TEST_EQUAL(cfg.level(2).max_targets, 4)
  TEST_EQUAL(cfg.level(2).min_charge, 1)
  TEST_EQUAL((int)cfg.level(2).exploration, (int)ExplorationMetric::MassCount)
  TEST_REAL_SIMILAR(cfg.level(2).ce_min, 21.0)
  TEST_REAL_SIMILAR(cfg.level(2).ce_max, 39.0)
  TEST_REAL_SIMILAR(cfg.level(2).ce_step, 5.0)
  TEST_REAL_SIMILAR(cfg.level(2).remaining_precursor_target, 0.12)
  TEST_EQUAL((int)cfg.level(3).selection, (int)SelectionMetric::None)
}
END_SECTION

END_TEST
