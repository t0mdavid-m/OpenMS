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

START_SECTION(([EXTRA] Config rejects exploration with multiple scan configs))
{
  const char* json = R"({
  "deconvolution": {
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
  "faims": {},
  "scheduling": {
    "cycle_time": {
      "enabled": false
    },
    "scan_timeout": {
      "enabled": false
    },
    "agc_interval_seconds": 30
  },
  "files": {},
  "precursor_selection": {
    "targeting": "none",
    "rank_by": "qscore",
    "max_precursors": 3,
    "additional_scans": [
      "secondary"
    ],
    "exploration": {
      "metric": "mass_count",
      "ce_min": 20,
      "ce_max": 40,
      "ce_step": 5,
      "activations": [
        "HCD"
      ]
    }
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
    },
    "additional_ms2": {
      "secondary": {
        "analyzer": "Orbitrap",
        "activation": "ETD",
        "collision_energy": 0,
        "reaction_time": 10.0,
        "resolution": 60000
      }
    }
  },
  "quantification": {
    "enabled": false
  }
}
)";
  bool threw = false;
  try { Config cfg{std::string(json)}; }
  catch (const std::invalid_argument&) { threw = true; }
  TEST_EQUAL(threw, true)
  (void)threw;
}
END_SECTION

START_SECTION(([EXTRA] Config rejects legacy ms3.enabled key))
{
  const char* json = R"({
  "deconvolution": {
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
  "faims": {},
  "scheduling": {
    "cycle_time": {
      "enabled": false
    },
    "scan_timeout": {
      "enabled": false
    },
    "agc_interval_seconds": 30
  },
  "files": {},
  "ms3": {
    "enabled": false,
    "protein_sequence": ""
  },
  "precursor_selection": {
    "targeting": "none",
    "rank_by": "qscore",
    "max_precursors": 3
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
  "quantification": {
    "enabled": false
  }
}
)";
  bool threw = false;
  try { Config cfg{std::string(json)}; }
  catch (const std::invalid_argument&) { threw = true; }
  TEST_EQUAL(threw, true)
  (void)threw;
}
END_SECTION

START_SECTION(([EXTRA] Config rejects legacy ms3.active key))
{
  const char* json = R"({
  "deconvolution": {
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
  "faims": {},
  "scheduling": {
    "cycle_time": {
      "enabled": false
    },
    "scan_timeout": {
      "enabled": false
    },
    "agc_interval_seconds": 30
  },
  "files": {},
  "ms3": {
    "active": true,
    "protein_sequence": ""
  },
  "precursor_selection": {
    "targeting": "none",
    "rank_by": "qscore",
    "max_precursors": 3
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
  "quantification": {
    "enabled": false
  }
}
)";
  bool threw = false;
  try { Config cfg{std::string(json)}; }
  catch (const std::invalid_argument&) { threw = true; }
  TEST_EQUAL(threw, true)
  (void)threw;
}
END_SECTION

START_SECTION(([EXTRA] Config rejects legacy ms3.mode key))
{
  const char* json = R"({
  "deconvolution": {
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
  "faims": {},
  "scheduling": {
    "cycle_time": {
      "enabled": false
    },
    "scan_timeout": {
      "enabled": false
    },
    "agc_interval_seconds": 30
  },
  "files": {},
  "ms3": {
    "mode": 1,
    "protein_sequence": ""
  },
  "precursor_selection": {
    "targeting": "none",
    "rank_by": "qscore",
    "max_precursors": 3
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
  "quantification": {
    "enabled": false
  }
}
)";
  bool threw = false;
  try { Config cfg{std::string(json)}; }
  catch (const std::invalid_argument&) { threw = true; }
  TEST_EQUAL(threw, true)
  (void)threw;
}
END_SECTION

START_SECTION(([EXTRA] Config rejects legacy ms3.max_per_ms2 key))
{
  const char* json = R"({
  "deconvolution": {
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
  "faims": {},
  "scheduling": {
    "cycle_time": {
      "enabled": false
    },
    "scan_timeout": {
      "enabled": false
    },
    "agc_interval_seconds": 30
  },
  "files": {},
  "ms3": {
    "max_per_ms2": 4,
    "protein_sequence": ""
  },
  "precursor_selection": {
    "targeting": "none",
    "rank_by": "qscore",
    "max_precursors": 3
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
  "quantification": {
    "enabled": false
  }
}
)";
  bool threw = false;
  try { Config cfg{std::string(json)}; }
  catch (const std::invalid_argument&) { threw = true; }
  TEST_EQUAL(threw, true)
  (void)threw;
}
END_SECTION

START_SECTION(([EXTRA] Config accepts an identification section with only a reference))
{
  const char* json = R"({
  "deconvolution": {
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
  "faims": {},
  "scheduling": {
    "cycle_time": {
      "enabled": false
    },
    "scan_timeout": {
      "enabled": false
    },
    "agc_interval_seconds": 30
  },
  "files": {},
  "precursor_selection": {
    "targeting": "none",
    "rank_by": "qscore",
    "max_precursors": 3
  },
  "identification": {
    "reference": {
      "sequence": "MKWVTFISLLLLFSSAYSRGVFRR"
    }
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
  "quantification": {
    "enabled": false
  }
}
)";
  Config cfg{std::string(json)};
  TEST_EQUAL(cfg.identification().sequence, "MKWVTFISLLLLFSSAYSRGVFRR")
}
END_SECTION

// The fourth characterization mode (ADR-0023). The objective assertion is the point of this
// section: the engine branches on `objective`, never on `mode`, so a mode that parsed without
// deriving its own objective would inherit the in-class Ambiguity default and run byte-identically
// to "ambiguity" -- accepted, logged as exhaustive, and behaving as something else.
START_SECTION(([EXTRA] characterization mode exhaustive parses and carries its own objective))
{
  const std::string json = R"JSON({
    "deconvolution": {"tol": [10, 10, 10]},
    "precursor_selection": {"rank_by": "qscore", "max_precursors": 2},
    "identification": {
      "reference": {
        "sequence": "PEPTIDEPEPTIDE"
      },
      "characterization": {
        "mode": "exhaustive",
        "max_targets": 3,
        "min_target_mass": 0
      }
    },
    "ms_settings": {"ms1": {}, "ms2": {"activation": "HCD", "collision_energy": 29},
                    "ms3": {"activation": "HCD", "collision_energy": 25}}
  })JSON";
  Config cfg(json);
  TEST_EQUAL((int)cfg.characterization().mode, (int)CharacterizationMode::Exhaustive)
  TEST_EQUAL((int)cfg.characterization().objective, (int)CharacterizationObjective::Exhaustive)
  TEST_REAL_SIMILAR(cfg.characterization().min_target_mass, 0.0)
}
END_SECTION

// A fourth legal value must not have loosened the strictness that made mode the single MS3 switch,
// and the new key must be READ, not merely tolerated by the allowlist -- a key that passes
// validation and is then discarded is this file's oldest failure mode.
START_SECTION(([EXTRA] a mistyped mode still throws, and min_target_mass is honoured))
{
  const std::string bad = R"JSON({
    "deconvolution": {"tol": [10, 10, 10]},
    "identification": { "reference": { "sequence": "PEPTIDE" }, "characterization": { "mode": "exhaustve" } },
    "ms_settings": {"ms1": {}, "ms2": {"activation": "HCD", "collision_energy": 29},
                    "ms3": {"activation": "HCD", "collision_energy": 25}}
  })JSON";
  // Assert on the MESSAGE, not merely on the type. Config.cpp throws std::invalid_argument for
  // EVERY error in the file, so a type-only catch goes green off an unrelated defect -- it was
  // passing off a missing collision_energy in this very fixture, and it would still pass with the
  // whole mode allowlist deleted, which is the one regression this section exists to catch.
  bool threw_for_mode = false;
  try { Config cfg{bad}; }
  catch (const std::invalid_argument& e)
  {
    threw_for_mode = std::string(e.what()).find("characterization.mode") != std::string::npos;
  }
  TEST_EQUAL(threw_for_mode, true)

  const std::string floored = R"JSON({
    "deconvolution": {"tol": [10, 10, 10]},
    "precursor_selection": {"rank_by": "qscore", "max_precursors": 2},
    "identification": {
      "reference": {
        "sequence": "PEPTIDEPEPTIDE"
      },
      "characterization": {
        "mode": "exhaustive",
        "min_target_mass": 1500.5
      }
    },
    "ms_settings": {"ms1": {}, "ms2": {"activation": "HCD", "collision_energy": 29},
                    "ms3": {"activation": "HCD", "collision_energy": 25}}
  })JSON";
  Config cfg2(floored);
  TEST_REAL_SIMILAR(cfg2.characterization().min_target_mass, 1500.5)
}
END_SECTION

END_TEST
