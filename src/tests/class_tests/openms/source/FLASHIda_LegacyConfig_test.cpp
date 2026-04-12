// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include <stdexcept>

using namespace OpenMS;

START_TEST(FLASHIda_LegacyConfig, "$Id$")

START_SECTION(legacy_config_rejected_via_Config)
{
  // Config directly rejects non-JSON input
  bool threw = false;
  try
  {
    Config cfg("10 100 1 10 5 0.5 -1");
  }
  catch (const std::invalid_argument&)
  {
    threw = true;
  }
  TEST_EQUAL(threw, true)
  (void)threw; // MSVC C4189
}
END_SECTION

START_SECTION(empty_config_rejected_via_Config)
{
  bool threw = false;
  try
  {
    Config cfg("");
  }
  catch (const std::invalid_argument&)
  {
    threw = true;
  }
  TEST_EQUAL(threw, true)
  (void)threw; // MSVC C4189
}
END_SECTION

START_SECTION(legacy_config_rejected_via_FLASHIda)
{
  // FLASHIda constructor delegates to Config and propagates the exception
  const char* legacy_input = "10 100 1 10 5 0.5 -1";
  bool threw = false;
  try
  {
    FLASHIda ida(const_cast<char*>(legacy_input));
  }
  catch (const std::invalid_argument&)
  {
    threw = true;
  }
  TEST_EQUAL(threw, true)
  (void)threw; // MSVC C4189
}
END_SECTION

START_SECTION(empty_config_rejected_via_FLASHIda)
{
  const char* empty_input = "";
  bool threw = false;
  try
  {
    FLASHIda ida(const_cast<char*>(empty_input));
  }
  catch (const std::invalid_argument&)
  {
    threw = true;
  }
  TEST_EQUAL(threw, true)
  (void)threw; // MSVC C4189
}
END_SECTION

START_SECTION(([EXTRA] Config rejects IDScore + exploration combo))
{
  const char* json = R"({
    "deconvolution": { "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10] },
    "precursor_selection": { "IDScore": true, "HCDEnergy": -1 },
    "tagging": {},
    "quantification": { "enabled": false },
    "faims": {},
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [{ "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }]
    },
    "scheduling": { "cycle_time": { "enabled": false }, "scan_timeout": { "enabled": false }, "agc_interval_seconds": 30 },
    "files": {},
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_precursors": 3 },
      "ms2": { "selection": "intensity", "exploration": { "metric": "mass_count", "ce_min": 20, "ce_max": 40, "ce_step": 5, "activation": "HCD" } }
    }
  })";
  bool threw = false;
  try { Config cfg{std::string(json)}; }
  catch (const std::invalid_argument&) { threw = true; }
  TEST_EQUAL(threw, true)
  (void)threw;
}
END_SECTION

START_SECTION(([EXTRA] Config rejects exploration with multiple scan configs))
{
  const char* json = R"({
    "deconvolution": { "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10] },
    "precursor_selection": { "IDScore": false },
    "tagging": {},
    "quantification": { "enabled": false },
    "faims": {},
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 60000 }
      ]
    },
    "scheduling": { "cycle_time": { "enabled": false }, "scan_timeout": { "enabled": false }, "agc_interval_seconds": 30 },
    "files": {},
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_precursors": 3 },
      "ms2": { "selection": "intensity", "exploration": { "metric": "mass_count", "ce_min": 20, "ce_max": 40, "ce_step": 5, "activation": "HCD" } }
    }
  })";
  bool threw = false;
  try { Config cfg{std::string(json)}; }
  catch (const std::invalid_argument&) { threw = true; }
  TEST_EQUAL(threw, true)
  (void)threw;
}
END_SECTION

START_SECTION(([EXTRA] Config rejects conditional_ms2 without tagging follow_up_scan))
{
  const char* json = R"({
    "deconvolution": { "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10] },
    "precursor_selection": {},
    "tagging": { "min_tag_length": 3 },
    "conditional_ms2": true,
    "quantification": { "enabled": false },
    "faims": {},
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [{ "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }]
    },
    "scheduling": { "cycle_time": { "enabled": false }, "scan_timeout": { "enabled": false }, "agc_interval_seconds": 30 },
    "files": {},
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_precursors": 3 },
      "ms2": { "selection": "intensity" }
    }
  })";
  bool threw = false;
  try { Config cfg{std::string(json)}; }
  catch (const std::invalid_argument&) { threw = true; }
  TEST_EQUAL(threw, true)
  (void)threw;
}
END_SECTION

END_TEST
