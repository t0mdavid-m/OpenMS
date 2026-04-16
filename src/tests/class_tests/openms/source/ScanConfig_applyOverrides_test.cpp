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
#include <unordered_map>

using namespace OpenMS;

START_TEST(ScanConfig_applyOverrides, "$Id$")

// Test that all 17 ScanConfig fields can be overridden
START_SECTION(all_17_fields_overridden)
{
  ScanConfig sc;

  // Verify defaults before override
  TEST_EQUAL(sc.analyzer, "Orbitrap")
  TEST_EQUAL(sc.activation.empty(), true)
  TEST_EQUAL(sc.collision_energy, 0)
  TEST_EQUAL(sc.resolution, 0)
  TEST_EQUAL(sc.agc_target, 0)
  TEST_REAL_SIMILAR(sc.first_mass, 0.0)
  TEST_REAL_SIMILAR(sc.last_mass, 0.0)
  TEST_REAL_SIMILAR(sc.max_it, 0.0)
  TEST_EQUAL(sc.microscans, 0)
  TEST_REAL_SIMILAR(sc.rf_lens, 0.0)
  TEST_REAL_SIMILAR(sc.source_cid, 0.0)
  TEST_REAL_SIMILAR(sc.source_cid_scaling, 0.0)
  TEST_EQUAL(sc.data_type.empty(), true)
  TEST_EQUAL(sc.scan_rate.empty(), true)
  TEST_REAL_SIMILAR(sc.reaction_time, 0.0)
  TEST_REAL_SIMILAR(sc.reagent_max_it, 0.0)
  TEST_EQUAL(sc.reagent_agc_target, 0)

  // Apply all 17 overrides
  std::unordered_map<std::string, std::string> overrides = {
    {"analyzer", "FTMS"},
    {"activation", "ETD"},
    {"collision_energy", "35"},
    {"resolution", "60000"},
    {"agc_target", "500000"},
    {"first_mass", "200.5"},
    {"last_mass", "3000.0"},
    {"max_it", "100.0"},
    {"microscans", "3"},
    {"rf_lens", "45.5"},
    {"source_cid", "10.0"},
    {"source_cid_scaling", "0.8"},
    {"data_type", "Centroid"},
    {"scan_rate", "Turbo"},
    {"reaction_time", "15.0"},
    {"reagent_max_it", "200.0"},
    {"reagent_agc_target", "100000"}
  };

  sc.applyOverrides(overrides);

  // Verify all 17 fields changed
  TEST_EQUAL(sc.analyzer, "FTMS")
  TEST_EQUAL(sc.activation, "ETD")
  TEST_EQUAL(sc.collision_energy, 35)
  TEST_EQUAL(sc.resolution, 60000)
  TEST_EQUAL(sc.agc_target, 500000)
  TEST_REAL_SIMILAR(sc.first_mass, 200.5)
  TEST_REAL_SIMILAR(sc.last_mass, 3000.0)
  TEST_REAL_SIMILAR(sc.max_it, 100.0)
  TEST_EQUAL(sc.microscans, 3)
  TEST_REAL_SIMILAR(sc.rf_lens, 45.5)
  TEST_REAL_SIMILAR(sc.source_cid, 10.0)
  TEST_REAL_SIMILAR(sc.source_cid_scaling, 0.8)
  TEST_EQUAL(sc.data_type, "Centroid")
  TEST_EQUAL(sc.scan_rate, "Turbo")
  TEST_REAL_SIMILAR(sc.reaction_time, 15.0)
  TEST_REAL_SIMILAR(sc.reagent_max_it, 200.0)
  TEST_EQUAL(sc.reagent_agc_target, 100000)
}
END_SECTION

// Test that unknown keys are silently ignored
START_SECTION(unknown_keys_ignored)
{
  ScanConfig sc;
  std::unordered_map<std::string, std::string> overrides = {
    {"nonexistent_key", "42"},
    {"resolution", "120000"}
  };

  sc.applyOverrides(overrides);

  // Known key applied, unknown key ignored, other fields untouched
  TEST_EQUAL(sc.resolution, 120000)
  TEST_EQUAL(sc.analyzer, "Orbitrap")
  TEST_EQUAL(sc.collision_energy, 0)
}
END_SECTION

END_TEST
