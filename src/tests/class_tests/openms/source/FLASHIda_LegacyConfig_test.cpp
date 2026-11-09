// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteomeSearch.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>

#include <nlohmann/json.hpp>

#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;

namespace
{
  // A loadable config around one authored identification section (no characterization, so no ms3).
  std::string withIdentification(const std::string& identification)
  {
    return std::string(R"JSON({
    "deconvolution": {"tol": [10, 10, 10]},
    "precursor_selection": {"rank_by": "qscore", "max_precursors": 2},
    "identification": )JSON") + identification + R"JSON(,
    "ms_settings": {"ms1": {}, "ms2": {"activation": "HCD", "collision_energy": 29}}
  })JSON";
  }

  // The message Config throws for @p json, or "" when it loads.
  std::string loadError(const std::string& json)
  {
    try { Config cfg{json}; }
    catch (const std::invalid_argument& e) { return e.what(); }
    return "";
  }

  bool contains(const std::string& s, const std::string& needle) { return s.find(needle) != std::string::npos; }

  // FLASHTnT's default for every knob, with ADR-0061's one deliberate difference.
  Param flashtntDefaultsWithSkipInference()
  {
    Param p = FLASHTnTAlgorithm().getDefaults();
    p.setValue("ex:skip_precursor_inference", "true");
    return p;
  }

  // A FLASHTnT Param value against a JSON value as the schema writes it: FLASHTnT's "true"/"false" strings
  // are JSON booleans there.
  bool sameValue(const ParamValue& v, const nlohmann::json& j)
  {
    switch (v.valueType())
    {
      case ParamValue::INT_VALUE: return j.is_number_integer() && j.get<long long>() == (long long)(int)v;
      case ParamValue::DOUBLE_VALUE: return j.is_number() && j.get<double>() == (double)v;
      case ParamValue::STRING_VALUE:
      {
        const std::string s = v.toString();
        if (s == "true" || s == "false") return j.is_boolean() && j.get<bool>() == (s == "true");
        return j.is_string() && j.get<std::string>() == s;
      }
      case ParamValue::STRING_LIST:
        return j.is_array() && j.get<std::vector<std::string>>() == v.toStringVector();
      default: return false;
    }
  }

  // Follows a "#/$defs/<name>" $ref, as method.schema.json writes them.
  const nlohmann::json& resolve(const nlohmann::json& schema, const nlohmann::json& node)
  {
    const nlohmann::json* n = &node;
    for (int depth = 0; n->contains("$ref") && depth < 32; ++depth)
    {
      const std::string ref = n->at("$ref").get<std::string>();
      n = &schema.at("$defs").at(ref.substr(std::string("#/$defs/").size()));
    }
    return *n;
  }
}

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

// ADR-0061 d.13: the reference is a FASTA path XOR one inline sequence. Every refusal names the keys and
// what to author instead; null is unauthored, as the C# serializer emits it.
START_SECTION(reference_is_fasta_xor_sequence)
{
  const std::pair<std::string, std::string> refused[] = {
    {R"({ "reference": { "fasta": "C:/db/ecoli.fasta", "sequence": "PEPTIDEK" } })", "both 'fasta' and 'sequence'"},
    {R"({ "reference": {} })", "neither 'fasta' nor 'sequence'"},
    {R"({ "reference": { "fasta": null, "sequence": null } })", "neither 'fasta' nor 'sequence'"},
    {R"({ "reference": { "fasta": "" } })", "identification.reference.fasta is empty"},
    {R"({ "reference": { "sequence": "" } })", "identification.reference.sequence is empty"},
    {R"({ "reference": { "fasta": "", "sequence": null } })", "identification.reference.fasta is empty"},
  };
  for (const auto& [ident, needle] : refused)
  {
    const std::string msg = loadError(withIdentification(ident));
    TEST_EQUAL(contains(msg, needle), true)
    TEST_EQUAL(contains(msg, "ADR-0061"), true)
  }

  {
    Config cfg{withIdentification(R"({ "reference": { "fasta": "C:/db/ecoli.fasta" } })")};
    TEST_STRING_EQUAL(cfg.identification().fasta, "C:/db/ecoli.fasta")
    TEST_STRING_EQUAL(cfg.identification().sequence, "")
    TEST_EQUAL(cfg.identification().enabled(), true)
  }
  {
    // The C# side emits the unauthored half as null.
    Config cfg{withIdentification(R"({ "reference": { "fasta": "db/relative.fasta", "sequence": null } })")};
    TEST_STRING_EQUAL(cfg.identification().fasta, "db/relative.fasta")  // as authored, never resolved here
    TEST_EQUAL(cfg.identification().enabled(), true)
  }
  {
    Config cfg{withIdentification(R"({ "reference": { "sequence": "PEPTIDEK", "fasta": null } })")};
    TEST_STRING_EQUAL(cfg.identification().sequence, "PEPTIDEK")
    TEST_STRING_EQUAL(cfg.identification().fasta, "")
    TEST_EQUAL(cfg.identification().enabled(), true)
  }
}
END_SECTION

// ADR-0067: identification.flashtnt.fixed_mod is the search's fixed modifications, built into one table at load
// (FixedModifications: a name or Unimod accession with its residue, residue anywhere, one per residue). A refused
// modification stops the load with the table's message; a string is not a list; [] and null give an empty table.
START_SECTION(fixed_mod_builds_the_table)
{
  // Delimited raw strings: a parenthesised mod name contains the default raw-string terminator.
  Config cfg{withIdentification(
      R"JSON({ "reference": { "sequence": "PEPTIDEK" }, "flashtnt": { "fixed_mod": ["iodoTMT6plex (C)"] } })JSON")};
  TEST_EQUAL(cfg.identification().flashtnt.fixed_mod.size(), 1)
  TEST_STRING_EQUAL(cfg.identification().flashtnt.fixed_mod[0], "iodoTMT6plex (C)")
  TEST_REAL_SIMILAR(cfg.identification().fixed_mods.delta('C'), 329.226595)
  TEST_EQUAL(cfg.identification().fixed_mods.delta('M'), 0.0)

  TEST_EQUAL(contains(loadError(withIdentification(
      R"JSON({ "reference": { "sequence": "PEPTIDEK" }, "flashtnt": { "fixed_mod": "Carbamidomethyl (C)" } })JSON")),
      "identification.flashtnt.fixed_mod must be a list of modification names"), true)
  const std::string unknown = loadError(withIdentification(
      R"JSON({ "reference": { "sequence": "PEPTIDEK" }, "flashtnt": { "fixed_mod": ["NoSuchModification (C)"] } })JSON"));
  TEST_EQUAL(contains(unknown, "Config: identification.flashtnt.fixed_mod NoSuchModification (C)"), true)
  TEST_EQUAL(contains(loadError(withIdentification(
      R"JSON({ "reference": { "sequence": "PEPTIDEK" }, "flashtnt": { "fixed_mod": ["Acetyl (Protein N-term)"] } })JSON")),
      "only a modification of a residue anywhere can be fixed"), true)

  for (const std::string ft : {std::string(R"("fixed_mod": [])"), std::string(R"("fixed_mod": null)")})
  {
    Config empty{withIdentification(R"({ "reference": { "sequence": "PEPTIDEK" }, "flashtnt": { )" + ft + " } }")};
    TEST_EQUAL(empty.identification().flashtnt.fixed_mod.empty(), true)
    TEST_EQUAL(empty.identification().fixed_mods.empty(), true)
  }
}
END_SECTION

// ADR-0061 d.6 (push 5): identification.threshold. Absent (or null) = identify and log, never close out; an authored
// finite number > 0 round-trips; zero, a negative, a string and a non-number are refused with a message naming the key.
START_SECTION(threshold_absent_never_closes_out_and_authored_round_trips)
{
  {
    Config cfg{withIdentification(R"({ "reference": { "sequence": "PEPTIDEK" } })")};
    TEST_EQUAL(cfg.identification().threshold.has_value(), false)
  }
  {
    Config cfg{withIdentification(R"({ "reference": { "sequence": "PEPTIDEK" }, "threshold": null })")};
    TEST_EQUAL(cfg.identification().threshold.has_value(), false)
  }
  {
    Config cfg{withIdentification(R"({ "reference": { "sequence": "PEPTIDEK" }, "threshold": 42.5 })")};
    TEST_EQUAL(cfg.identification().threshold.has_value(), true)
    TEST_REAL_SIMILAR(*cfg.identification().threshold, 42.5)
  }
  for (const std::string bad : {std::string("0"), std::string("-1"), std::string("\"high\""), std::string("[500]")})
  {
    const std::string msg = loadError(withIdentification(R"({ "reference": { "sequence": "PEPTIDEK" }, "threshold": )" + bad + " }"));
    TEST_EQUAL(contains(msg, "identification.threshold must be a finite number > 0"), true)
  }
}
END_SECTION

// ADR-0066: identification.characterization.per_proteoform. Absent or null = false (every committed config); true
// and false parse; a non-boolean is refused with a message naming the key.
START_SECTION(per_proteoform_absent_is_off_and_only_a_boolean_parses)
{
  auto withCharacterization = [](const std::string& extra) {
    return std::string(R"JSON({
    "deconvolution": {"tol": [10, 10, 10]},
    "precursor_selection": {"rank_by": "qscore", "max_precursors": 2},
    "identification": { "reference": { "sequence": "PEPTIDEPEPTIDE" },
                        "characterization": { "mode": "fragments", "max_targets": 3)JSON") + extra + R"JSON( } },
    "ms_settings": {"ms1": {}, "ms2": {"activation": "HCD", "collision_energy": 29},
                    "ms3": {"activation": "CID", "collision_energy": 35}}
  })JSON";
  };
  TEST_EQUAL(Config{withCharacterization("")}.characterization().per_proteoform, false)
  TEST_EQUAL(Config{withCharacterization(R"(, "per_proteoform": null)")}.characterization().per_proteoform, false)
  TEST_EQUAL(Config{withCharacterization(R"(, "per_proteoform": false)")}.characterization().per_proteoform, false)
  TEST_EQUAL(Config{withCharacterization(R"(, "per_proteoform": true)")}.characterization().per_proteoform, true)
  for (const std::string bad : {std::string("1"), std::string("\"yes\""), std::string("[true]")})
    TEST_EQUAL(contains(loadError(withCharacterization(", \"per_proteoform\": " + bad)),
                        "identification.characterization.per_proteoform must be true or false"), true)
}
END_SECTION

// ADR-0061 d.13: "Every unset key takes FLASHTnT's default". An identification without flashtnt, with an
// empty one and with a null one gives FLASHTnT's Param for every knob (ex:skip_precursor_inference aside,
// ADR-0061's true), and authored values reach FLASHTnT under its own names.
START_SECTION(flashtnt_unset_keys_take_flashtnt_defaults)
{
  const Param tnt = FLASHTnTAlgorithm().getDefaults();
  const Param realtime = ProteomeSearch::defaultParameters();
  const auto& knobs = ProteomeSearch::knobNames();
  TEST_EQUAL(knobs.size(), 14)

  std::set<std::string> keys, names;
  for (const auto& [key, name] : knobs)
  {
    keys.insert(key);
    names.insert(name);
    TEST_EQUAL(tnt.exists(name), true)  // every mapped name is a FLASHTnT parameter
  }
  TEST_EQUAL(keys.size(), knobs.size())
  TEST_EQUAL(names.size(), knobs.size())
  TEST_EQUAL(keys.count("fixed_mod"), 1)  // ADR-0067: FLASHTnTAlgorithm's own fixed_mod

  for (const std::string unset : {std::string(""), std::string(R"(, "flashtnt": {})"), std::string(R"(, "flashtnt": null)")})
  {
    Config cfg{withIdentification(R"({ "reference": { "sequence": "PEPTIDEK" })" + unset + " }")};
    TEST_EQUAL(cfg.identification().flashtnt.fixed_mod.empty(), true)
    const Param p = ProteomeSearch::parametersFor(cfg.identification().flashtnt);
    for (const auto& [key, name] : knobs)
    {
      TEST_EQUAL(p.getValue(name) == realtime.getValue(name), true)
      if (name == "ex:skip_precursor_inference")
      {
        TEST_STRING_EQUAL(p.getValue(name).toString(), "true")           // ADR-0061's choice
        TEST_STRING_EQUAL(tnt.getValue(name).toString(), "false")        // against FLASHTnT's own
      }
      else
      {
        TEST_EQUAL(p.getValue(name) == tnt.getValue(name), true)
      }
    }
    // Nothing but the knobs is touched: every other parameter stays at defaultParameters()'.
    for (auto it = realtime.begin(); it != realtime.end(); ++it)
    {
      if (names.count(it.getName())) continue;
      TEST_EQUAL(p.getValue(it.getName()) == it->value, true)
    }
  }

  // Every knob authored at a value that is not FLASHTnT's default reaches its FLASHTnT name.
  const std::string authored = R"JSON({ "reference": { "sequence": "PEPTIDEK" }, "flashtnt": {
      "min_length": 4, "max_length": 9, "allow_gap": true, "max_aa_in_gap": 3,
      "max_blind_mod_count": 2, "max_mod_mass": 616.5,
      "var_mod": ["Oxidation (M)", "Acetyl (N-term)"], "max_var_mod_count": 1, "max_protein_hits": 7,
      "gapped_tags": "off", "terminal_ladders": false, "precursor_ladders": false,
      "skip_precursor_inference": false, "fixed_mod": ["Carbamidomethyl (C)"] } })JSON";
  const std::vector<std::pair<std::string, ParamValue>> expected = {
    {"min_length", ParamValue(4)},
    {"max_length", ParamValue(9)},
    {"allow_gap", ParamValue("true")},
    {"max_aa_in_gap", ParamValue(3)},
    {"fixed_mod", ParamValue(std::vector<std::string>{"Carbamidomethyl (C)"})},
    {"max_blind_mod_count", ParamValue(2)},
    {"max_mod_mass", ParamValue(616.5)},
    {"var_mod", ParamValue(std::vector<std::string>{"Oxidation (M)", "Acetyl (N-term)"})},
    {"max_var_mod_count", ParamValue(1)},
    {"max_protein_hits", ParamValue(7)},
    {"gapped_tags", ParamValue("off")},
    {"terminal_ladders", ParamValue("false")},
    {"precursor_ladders", ParamValue("false")},
    {"skip_precursor_inference", ParamValue("false")},
  };
  TEST_EQUAL(expected.size(), knobs.size())
  Config cfg{withIdentification(authored)};
  const Param p = ProteomeSearch::parametersFor(cfg.identification().flashtnt);
  const Param defaults = flashtntDefaultsWithSkipInference();
  for (const auto& [key, name] : knobs)
  {
    const ParamValue* want = nullptr;
    for (const auto& e : expected)
      if (e.first == key) want = &e.second;
    TEST_EQUAL(want != nullptr, true)
    if (want == nullptr) continue;
    TEST_EQUAL(p.getValue(name) == *want, true)
    TEST_EQUAL(p.getValue(name) == defaults.getValue(name), false)  // the value is not the default, so it was read
  }

  // FLASHTnT's gapped-tag modes that need the whole run are refused, a typo too.
  for (const std::string v : {std::string("both"), std::string("only"), std::string("Fallback")})
  {
    const std::string msg = loadError(withIdentification(
        R"({ "reference": { "sequence": "PEPTIDEK" }, "flashtnt": { "gapped_tags": ")" + v + R"(" } })"));
    TEST_EQUAL(contains(msg, "identification.flashtnt.gapped_tags"), true)
  }
}
END_SECTION

// The key documentation (docs/site/method.schema.json, ADR-0058) states every identification.flashtnt default:
// its key set is the knobs plus fixed_mod, in both directions, and each default is FLASHTnT's (with
// skip_precursor_inference true, ADR-0061), so the published defaults cannot drift from what the engine runs.
START_SECTION(flashtnt_schema_defaults_equal_flashtnt)
{
  std::string text;
  {
    std::ifstream in("../../FlashIDA/docs/site/method.schema.json");
    std::stringstream ss;
    ss << in.rdbuf();
    text = ss.str();
  }
  ABORT_IF(text.empty())
  const nlohmann::json schema = nlohmann::json::parse(text);
  const nlohmann::json& ident = resolve(schema, schema.at("properties").at("identification"));
  const nlohmann::json& flashtnt = resolve(schema, ident.at("properties").at("flashtnt"));
  const nlohmann::json& props = flashtnt.at("properties");

  const auto& knobs = ProteomeSearch::knobNames();
  std::set<std::string> want{"fixed_mod"};
  for (const auto& kn : knobs) want.insert(kn.first);
  std::set<std::string> documented;
  for (auto it = props.begin(); it != props.end(); ++it) documented.insert(it.key());
  for (const auto& k : want) TEST_EQUAL(documented.count(k), 1)   // every knob is documented
  for (const auto& k : documented) TEST_EQUAL(want.count(k), 1)   // and nothing else is
  TEST_EQUAL(documented.size(), want.size())

  const Param defaults = flashtntDefaultsWithSkipInference();
  for (const auto& [key, name] : knobs)
  {
    const nlohmann::json& node = resolve(schema, props.at(key));
    TEST_EQUAL(node.contains("default"), true)
    if (!node.contains("default")) continue;
    TEST_EQUAL(sameValue(defaults.getValue(name), node.at("default")), true)
  }
  TEST_EQUAL(props.at("fixed_mod").contains("default") && props.at("fixed_mod").at("default") == nlohmann::json::array(), true)
}
END_SECTION

END_TEST
