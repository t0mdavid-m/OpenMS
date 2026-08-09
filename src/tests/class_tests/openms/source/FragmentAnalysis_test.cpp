// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Exploration.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include "FLASHIda_TestAccess.h"   // ExplorationTestAccess::feedResult (private-state access)

#include <fstream>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  struct ScanData
  {
    std::vector<double> mzs;
    std::vector<double> ints;
    double rt;
    std::string scan_id;
  };

  std::vector<ScanData> loadTsvScans(const std::string& path)
  {
    std::vector<ScanData> scans;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line))
    {
      if (line.substr(0, 4) == "Spec")
      {
        scans.emplace_back();
        auto tab = line.find('\t');
        scans.back().scan_id = line.substr(10, tab - 10);
        scans.back().rt = std::stod(line.substr(tab + 1));
      }
      else if (!scans.empty())
      {
        auto tab = line.find('\t');
        if (tab != std::string::npos)
        {
          scans.back().mzs.push_back(std::stod(line.substr(0, tab)));
          scans.back().ints.push_back(std::stod(line.substr(tab + 1)));
        }
      }
    }
    return scans;
  }

  // Horse heart cytochrome c sequence (104 AA)
  const char* cytochrome_c_seq = "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";

  const std::string ms2_cytc_path = "../../FlashIDA/test-data/spectra/ms2_cytc_scan149.txt";

  // JSON config for fragment analysis tests
  const char* fragment_test_config = R"({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.9,
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
  "flashtnt": {
    "min_length": 3,
    "max_length": 8
  },
  "faims": {
    "cv_values": [],
    "max_cv_skip": 0,
    "cv_precursor_threshold": 15
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": false,
      "value_ms": 30000
    }
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
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
    "mode": "ambiguity",
    "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "max_targets": 3
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
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "HCD",
      "collision_energy": 35,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)";

  // Config with fragment_count exploration metric + protein sequence
  const char* fragment_count_exploration_config = R"({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.9,
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
  "flashtnt": {
    "min_length": 3,
    "max_length": 8
  },
  "faims": {
    "cv_values": [],
    "max_cv_skip": 0,
    "cv_precursor_threshold": 15
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": false,
      "value_ms": 30000
    }
  },
  "files": {
    "target_logs": [],
    "fasta": "cytochrome_c.fasta",
    "inclusion_list": "",
    "ptm_list": ""
  },
  "conditional_ms2": false,
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3,
    "exploration": {
      "metric": "fragment_count",
      "ce_min": 20.0,
      "ce_max": 30.0,
      "ce_step": 10.0
    }
  },
  "characterization": {
    "mode": "off",
    "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "max_targets": 3
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
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "HCD",
      "collision_energy": 35,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)";

  // Config with mass_count exploration metric + protein sequence (same sequence present)
  const char* mass_count_exploration_config = R"({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.9,
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
  "flashtnt": {
    "min_length": 3,
    "max_length": 8
  },
  "faims": {
    "cv_values": [],
    "max_cv_skip": 0,
    "cv_precursor_threshold": 15
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": false,
      "value_ms": 30000
    }
  },
  "files": {
    "target_logs": [],
    "fasta": "cytochrome_c.fasta",
    "inclusion_list": "",
    "ptm_list": ""
  },
  "conditional_ms2": false,
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3,
    "exploration": {
      "metric": "mass_count",
      "ce_min": 20.0,
      "ce_max": 30.0,
      "ce_step": 10.0
    }
  },
  "characterization": {
    "mode": "off",
    "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "max_targets": 3
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
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "HCD",
      "collision_energy": 35,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)";

  // Param-builder fixture: every flashtnt key the tagger builder reads carries a NON-default value
  // (min_length 4 != 3, max_length 9 != 8, allow_gap true != false, max_aa_in_gap 3 != 2,
  // fixed_mod non-empty != {}), so a dropped setValue in buildTaggerParam cannot pass by accident.
  // NOTE the custom JSON delimiter -- do not remove it. This fixture carries a real ModificationsDB
  // modification name, Carbamidomethyl (C). In an UNdelimited raw string, a close-paren followed
  // directly by a double-quote is the terminator -- and that sequence occurs inside any
  // parenthesised mod name. Without the delimiter the literal ended 370 chars in, mid-JSON, and
  // everything after it was compiled as C++ (broke CI 2026-08-09).
  const char* tagger_param_config = R"JSON({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.9,
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
  "flashtnt": {
    "min_length": 4,
    "max_length": 9,
    "allow_gap": true,
    "max_aa_in_gap": 3,
    "fixed_mod": [
      "Carbamidomethyl (C)"
    ]
  },
  "faims": {
    "cv_values": [],
    "max_cv_skip": 0,
    "cv_precursor_threshold": 15
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": false,
      "value_ms": 30000
    }
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
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
    "mode": "ambiguity",
    "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "max_targets": 3
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
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "HCD",
      "collision_energy": 35,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)JSON";

  // Same, but with an EXPLICITLY EMPTY fixed_mod -- the case the deleted
  // `if (!config_.targeting().fixed_mod.empty())` guards used to suppress. max_blind_mod_count is
  // non-default (4 != 2) so the extender builder's carry-through is genuinely proven.
  const char* empty_fixed_mod_config = R"({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.9,
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
  "flashtnt": {
    "min_length": 3,
    "max_length": 8,
    "max_blind_mod_count": 4,
    "fixed_mod": []
  },
  "faims": {
    "cv_values": [],
    "max_cv_skip": 0,
    "cv_precursor_threshold": 15
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": false,
      "value_ms": 30000
    }
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
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
    "mode": "ambiguity",
    "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "max_targets": 3
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
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "HCD",
      "collision_energy": 35,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)";

  PeakGroup makeSyntheticPeakGroup(double mz, double mass, int charge)
  {
    PeakGroup pg(charge, charge, true);
    pg.setMonoisotopicMass(mass);
    FLASHHelperClasses::LogMzPeak lp;
    lp.mz = mz;
    lp.abs_charge = charge;
    pg.push_back(lp);
    return pg;
  }

} // anonymous namespace


START_TEST(FragmentAnalysis, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(getTopFragmentMatches_cytochrome_c)
{
  auto scans = loadTsvScans(ms2_cytc_path);
  ABORT_IF(scans.empty())

  Config cfg{std::string(fragment_test_config)};
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});

  // Deconvolve scan 149 (MS2, no precursor info)
  int pg_count = deconv.deconvolveMSn(scans[0].mzs.data(), scans[0].ints.data(),
                                       (int)scans[0].mzs.size(), scans[0].rt, 0.0, 0);
  TEST_EQUAL(pg_count > 0, true)

  FragmentAnalysis fragments(cfg);

  double masses[20], qscores[20], wstarts[20], wends[20];
  int charges[20], frag_indices[20];
  char ion_types[20];

  FragmentAnalysis::ProteoformMatch match_result;
  int found = fragments.getTopFragmentMatches(cytochrome_c_seq, 20,
                                               masses, qscores, charges,
                                               wstarts, wends, ion_types,
                                               frag_indices, deconv.storedMS2(), match_result);
  TEST_EQUAL(found > 0, true)

  for (int i = 0; i < found; ++i)
  {
    TEST_EQUAL(masses[i] > 0, true)
    TEST_EQUAL(charges[i] > 0, true)
    TEST_EQUAL(frag_indices[i] > 0, true)
    TEST_EQUAL(ion_types[i] == 'b' || ion_types[i] == 'y', true)
  }
}
END_SECTION

START_SECTION(getTerminalFragmentIons_cytochrome_c)
{
  auto scans = loadTsvScans(ms2_cytc_path);
  ABORT_IF(scans.empty())

  Config cfg{std::string(fragment_test_config)};
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});

  deconv.deconvolveMSn(scans[0].mzs.data(), scans[0].ints.data(),
                        (int)scans[0].mzs.size(), scans[0].rt, 0.0, 0);

  FragmentAnalysis fragments(cfg);

  double masses[20], qscores[20], wstarts[20], wends[20];
  int charges[20], frag_indices[20];
  char ion_types[20];

  FragmentAnalysis::ProteoformMatch term_result;
  int found = fragments.getTerminalFragmentIons(cytochrome_c_seq, 10,
                                                 masses, qscores, charges,
                                                 wstarts, wends, ion_types,
                                                 frag_indices, deconv.storedMS2(), term_result);
  TEST_EQUAL(found > 0, true)

  for (int i = 0; i < found; ++i)
  {
    TEST_EQUAL(masses[i] > 0, true)
    TEST_EQUAL(ion_types[i] == 'b' || ion_types[i] == 'y', true)
  }
}
END_SECTION

START_SECTION(getTopFragmentMatches_empty_spectrum_returns_zero)
{
  Config cfg{std::string(fragment_test_config)};
  FragmentAnalysis fragments(cfg);

  DeconvolvedSpectrum ds(1);

  double masses[20], qscores[20], wstarts[20], wends[20];
  int charges[20], frag_indices[20];
  char ion_types[20];

  FragmentAnalysis::ProteoformMatch empty_result;
  int found = fragments.getTopFragmentMatches(cytochrome_c_seq, 20,
                                               masses, qscores, charges,
                                               wstarts, wends, ion_types,
                                               frag_indices, ds, empty_result);
  TEST_EQUAL(found, 0)
}
END_SECTION

START_SECTION(getBestMS2Masses_returns_deconvolved_peaks)
{
  auto scans = loadTsvScans(ms2_cytc_path);
  ABORT_IF(scans.empty())

  Config cfg{std::string(fragment_test_config)};
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});

  deconv.deconvolveMSn(scans[0].mzs.data(), scans[0].ints.data(),
                        (int)scans[0].mzs.size(), scans[0].rt, 0.0, 0);

  FragmentAnalysis fragments(cfg);

  double masses[10], qscores[10], wstarts[10], wends[10];
  int charges[10];

  int found = fragments.getBestMS2Masses(10, masses, qscores, charges,
                                          wstarts, wends, deconv.storedMS2());
  TEST_EQUAL(found > 0, true)

  for (int i = 0; i < found; ++i)
  {
    TEST_EQUAL(masses[i] >= 500.0, true)
    TEST_EQUAL(wstarts[i] < wends[i], true)
  }
}
END_SECTION

START_SECTION(fragment_count_populated_for_fragment_count_metric)
{
  // Deconvolve real cytochrome c MS2 spectrum
  auto scans = loadTsvScans(ms2_cytc_path);
  ABORT_IF(scans.empty())

  Config cfg{std::string(fragment_count_exploration_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  // Deconvolve scan 149 to produce a spectrum with real fragment matches
  deconv.deconvolveMSn(scans[0].mzs.data(), scans[0].ints.data(),
                        static_cast<int>(scans[0].mzs.size()), scans[0].rt, 0.0, 0);

  // Initiate exploration group — CE variants 20 and 30
  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  ABORT_IF(cmds.empty())

  // Feed deconvolved spectrum to first variant
  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info = ExplorationTestAccess::feedResult(exploration,tracking_id, deconv.storedMS2(), 1.0, queue);

  // FragmentCount metric: fragment analysis should have run
  TEST_EQUAL(info.metric.fragment_count > 0, true)
  TEST_EQUAL(info.identification.matched_protein.empty(), false)
  TEST_EQUAL(info.identification.proteoform_sequence.empty(), false)
  TEST_STRING_EQUAL(info.identification.proteoform_sequence, std::string(cytochrome_c_seq))
}
END_SECTION

START_SECTION(fragment_analysis_populated_for_mass_count_metric)
{
  // Same real cytochrome c spectrum that produces fragment matches
  auto scans = loadTsvScans(ms2_cytc_path);
  ABORT_IF(scans.empty())

  Config cfg{std::string(mass_count_exploration_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  // Deconvolve — same data as above, would produce matches if fragment analysis ran
  deconv.deconvolveMSn(scans[0].mzs.data(), scans[0].ints.data(),
                        static_cast<int>(scans[0].mzs.size()), scans[0].rt, 0.0, 0);

  // Initiate exploration group with mass_count metric
  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  ABORT_IF(cmds.empty())

  // Feed same deconvolved spectrum
  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info = ExplorationTestAccess::feedResult(exploration,tracking_id, deconv.storedMS2(), 1.0, queue);

  // MassCount metric: fragment analysis still runs (populates metadata for all metrics)
  TEST_EQUAL(info.metric.fragment_count > 0, true)
  TEST_EQUAL(info.identification.matched_protein.empty(), false)
  TEST_EQUAL(info.identification.proteoform_sequence == cytochrome_c_seq, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// buildTaggerParam / buildExtenderParam -- the shared FLASHTnT Param builders.
//
// The same six tagger Params used to be assembled inline in TWO places (FragmentAnalysis.cpp and
// PrecursorSelection.cpp), each wrapping fixed_mod in `if (!config_.targeting().fixed_mod.empty())`.
// The builders replace both copies and set fixed_mod UNCONDITIONALLY. These sections cover the two
// things a reviewer cannot see by diffing two deleted blocks against one new function:
//   T1 -- nothing was dropped in the extraction (every key set, from non-default config values);
//   T2/T3 -- the removed guard stays removed. fixed_mod is inert inside FLASHTnT today, so NO other
//            test in either project can notice if someone reintroduces the conditional.
/////////////////////////////////////////////////////////////
START_SECTION(buildTaggerParam_sets_every_key)
{
  Config cfg{std::string(tagger_param_config)};

  // {"c","z"}, NOT {"b","y"}: {"b","y"} is FLASHTaggerAlgorithm's OWN declared default, so asserting
  // it would pass even if the builder never set ion_type at all.
  Param p = FragmentAnalysis::buildTaggerParam(cfg, {"c", "z"});

  TEST_EQUAL((int)p.getValue("min_length"), 4)
  TEST_EQUAL((int)p.getValue("max_length"), 9)
  TEST_STRING_EQUAL(p.getValue("allow_gap").toString(), std::string("true"))
  TEST_EQUAL((int)p.getValue("max_aa_in_gap"), 3)

  auto ion_types = p.getValue("ion_type").toStringVector();
  TEST_EQUAL((int)ion_types.size(), 2)
  TEST_STRING_EQUAL(ion_types[0], std::string("c"))
  TEST_STRING_EQUAL(ion_types[1], std::string("z"))

  auto fixed_mod = p.getValue("fixed_mod").toStringVector();
  TEST_EQUAL((int)fixed_mod.size(), 1)
  TEST_STRING_EQUAL(fixed_mod[0], std::string("Carbamidomethyl (C)"))
}
END_SECTION

START_SECTION(buildTaggerParam_sets_fixed_mod_when_empty)
{
  // THE regression guard for the removed `if (!fixed_mod.empty())`. With the guard back, the builder
  // would leave FLASHTaggerAlgorithm's declared {""} default in place and the size assertion fails.
  Config cfg{std::string(empty_fixed_mod_config)};

  Param p = FragmentAnalysis::buildTaggerParam(cfg, {"b", "y"});

  // exists() is true EITHER WAY -- getDefaults() already declares fixed_mod, so the guard never
  // removed the key, it only left the {""} placeholder in place. The size assertion below is the
  // load-bearing one; do not "simplify" it away.
  TEST_EQUAL(p.exists("fixed_mod"), true)
  TEST_EQUAL((int)p.getValue("fixed_mod").toStringVector().size(), 0)  // [] verbatim, NOT the declared {""}
}
END_SECTION

START_SECTION(buildExtenderParam_sets_fixed_mod_when_empty)
{
  Config cfg{std::string(empty_fixed_mod_config)};

  Param p = FragmentAnalysis::buildExtenderParam(cfg, {"c", "z"}, 555.0);

  // See buildTaggerParam_sets_fixed_mod_when_empty: exists() passes either way, the size is the guard.
  TEST_EQUAL(p.exists("fixed_mod"), true)
  TEST_EQUAL((int)p.getValue("fixed_mod").toStringVector().size(), 0)  // [] verbatim, NOT the declared {""}

  // The two values the extender builder must carry through: one from config, one from the argument.
  TEST_EQUAL((int)p.getValue("max_blind_mod_count"), 4)   // flashtnt.max_blind_mod_count (non-default)
  TEST_REAL_SIMILAR((double)p.getValue("max_mod_mass"), 555.0)
}
END_SECTION

START_SECTION(buildExtenderParam_sets_every_key)
{
  // The extender counterpart of buildTaggerParam_sets_every_key, and it exists for one key in
  // particular: skip_precursor_inference. The builder sets it "true" against the extender's own
  // declared default of "false", and NOTHING else in either project asserts it -- so if that
  // setValue were ever dropped, the extender would silently resume inferring precursor masses from
  // complementary fragment pairs during MS2 identification, with no test anywhere noticing.
  Config cfg{std::string(empty_fixed_mod_config)};

  // {"c","z"} again, not the declared {"b","y"} default -- otherwise a dropped ion_type set passes.
  Param p = FragmentAnalysis::buildExtenderParam(cfg, {"c", "z"}, 555.0);

  TEST_STRING_EQUAL(p.getValue("skip_precursor_inference").toString(), std::string("true"))

  auto ion_types = p.getValue("ion_type").toStringVector();
  TEST_EQUAL((int)ion_types.size(), 2)
  TEST_STRING_EQUAL(ion_types[0], std::string("c"))
  TEST_STRING_EQUAL(ion_types[1], std::string("z"))

  TEST_EQUAL((int)p.getValue("max_blind_mod_count"), 4)
  TEST_REAL_SIMILAR((double)p.getValue("max_mod_mass"), 555.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// narrowFragmentPTMSites -- per-scan MS3 identification-LEAF localizer (Change L: equiv-frame gap-partition).
//
// REGRESSION SUITE for the MS3 leaf flip/mislocalization bug. Every case INJECTS its matched-fragment set
// directly (no deconvolution -> deterministic, immune to the OpenMS build-nondeterminism that flaps score
// goldens) and drives the PRODUCTION narrowFragmentPTMSites + toProForma. Cases A/B are the two owner-validated
// ground-truth proformas (kept in sync with FlashIDA/test-data/reference/ms3_leaf_expected.tsv); their real
// matched ion sets come from the ms3_cytc capture as raw MS3 ion -> equivalent full-protein ion:
//   !!& (b80): b22,b26,b28 direct + yb67 -> b13 (complement flip)
//   !!' (b70): b26,b34,b44,b48,b61 direct + yb69 -> b1 (complement flip; b1 = 42.01 = M-89)
//
// Leg 1 (ions->leaf):     the injected equiv-frame ions localize each mod exactly as ground truth expects.
// Leg 2 (leaf<->pooled):  the leaf's include/exclude verdict == the pooled predicate
//                         MS3FragmentMatcher::coveredAmbiguousInEquivFrame for the same fragment+mod
//                         (the two localization paths cannot drift apart).
// The retired raw-frame test read fm.ion_type/ion_index; Change L reads fm.equiv_type/equiv_index, so it is
// rewritten here to the equiv contract.
/////////////////////////////////////////////////////////////
START_SECTION(narrowFragmentPTMSites)
{
  using PTM = FragmentAnalysis::PTMSite;
  using FM  = FragmentAnalysis::ProteoformMatch::FragmentMatch;

  // Inject one equiv-frame fragment: raw sub-ion (ion_type/ion_index) + its projected full-protein equivalent
  // (equiv_type/equiv_index) + the complement-flip flag + the sub-frame includes_ptm verdict. narrowFragmentPTMSites
  // reads ONLY the equiv fields + flip + verdict (the raw fields document provenance / the flip source).
  auto mkEq = [](const char* raw_t, int raw_i, const char* eq_t, int eq_i, bool flip, bool inc) {
    FM f; f.ion_type = raw_t; f.ion_index = raw_i; f.equiv_type = eq_t; f.equiv_index = eq_i;
    f.is_complement_flip = flip; f.includes_ptm = inc; return f; };

  // The 105-aa M-start cytC reference the ms3_cytc goldens render against (NOT the 104-aa des-Met
  // `cytochrome_c_seq` above). A b<k> MS3 precursor renders over the k-residue prefix CYTC105[0:k].
  const std::string CYTC105 =
    "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";

  // Ground truth (== test-data/reference/ms3_leaf_expected.tsv). -89.0302 = N-terminal net loss (Met-excision +
  // N-alpha-acetyl) anchored at residue 1; +615.2498 = heme; +526.2196 = -89.0302 + 615.2498 (the summed shift
  // when the two are co-observed and no fragment separates them).
  const std::string EXP_AMP  =  // !!& (b80): -89 over (1-13), +615 over (14-22)
    "(MGDVEKGKKIFVQ)[-89.0302](KCAQCHTVE)[+615.2498]KGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTK";
  const std::string EXP_APOS =  // !!' (b70): -89 localized to M1, +615 over (2-26)
    "M[-89.0302](GDVEKGKKIFVQKCAQCHTVEKGGK)[+615.2498]HKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLE";

  // The a-priori (wide) MS2 model base the two real scans narrow FROM.
  std::vector<PTM> wideAB = { PTM{4, 1, 8, -89.0302}, PTM{20, 15, 26, 615.2498} };

  // ---- Case A (!!&, b80): wide-seed to (1-13)/(14-22) ------------------------------------------------
  // b13 (=flipped yb67) fully-covers -89 (caps hi 13) and no-overlaps +615 (excludes -> lo 14); b22 straddle-
  // includes +615 (hi 22). yields -89 (1-13), +615 (14-22), rendered over the 80-residue b80 prefix.
  auto sitesA = FragmentAnalysis::narrowFragmentPTMSites(
      wideAB, 80, { mkEq("yb", 67, "b", 13, true, true), mkEq("b", 22, "b", 22, false, true),
                    mkEq("b", 26, "b", 26, false, true), mkEq("b", 28, "b", 28, false, true) });
  std::string pfA = FragmentAnalysis::toProForma(CYTC105.substr(0, 80), sitesA);
  TEST_STRING_EQUAL(pfA, EXP_AMP)

  // ---- Case B (!!', b70): THE flip case -------------------------------------------------------------
  // b1 (=flipped yb69 = M-89). The N-terminal-loss exception keeps -89 on M1 DESPITE the flip; b26 caps +615
  // at 26 and b1 excludes it from residue 1 (lo 2). yields -89 localized to M1, +615 (2-26).
  auto sitesB = FragmentAnalysis::narrowFragmentPTMSites(
      wideAB, 70, { mkEq("yb", 69, "b", 1, true, true), mkEq("b", 26, "b", 26, false, true),
                    mkEq("b", 34, "b", 34, false, true), mkEq("b", 44, "b", 44, false, true),
                    mkEq("b", 48, "b", 48, false, true), mkEq("b", 61, "b", 61, false, true) });
  std::string pfB = FragmentAnalysis::toProForma(CYTC105.substr(0, 70), sitesB);
  TEST_STRING_EQUAL(pfB, EXP_APOS)
  // ISSUE(old raw-frame narrower): the flipped b1 (raw yb69) was read as a suffix over [2..70] -> -89 dragged
  // off M1 to (2-8). Reading the EQUIV b1 + the nterm-loss exception pins it to M1.

  // ---- Case C: flip-INVERT (non-N-terminal) ---------------------------------------------------------
  // Controlled scenario for the flip-invert branch: a complement-flipped includer must INVERT its sub-verdict.
  // b15 excludes the heme (lo 16); the flipped b22 (sub_includes_ptm=false) inverts to an includer (hi 22),
  // keeping the heme in the His-region.
  auto sitesC = FragmentAnalysis::narrowFragmentPTMSites(
      { PTM{20, 15, 26, 615.2498} }, 80,
      { mkEq("b", 15, "b", 15, false, false), mkEq("yb", 58, "b", 22, true, false) });
  TEST_EQUAL((int)sitesC.size(), 1)
  TEST_EQUAL(sitesC[0].start_position, 16)
  TEST_EQUAL(sitesC[0].end_position, 22)   // ISSUE(no-invert): a non-inverted flipped b22 would EXCLUDE -> lo 23 -> wrong region

  // ---- Case D: co-observed MERGE --------------------------------------------------------------------
  // Drop the separating b13: nothing splits -89 from +615 this scan, so they MERGE into one summed shift.
  auto sitesD = FragmentAnalysis::narrowFragmentPTMSites(
      wideAB, 80, { mkEq("b", 22, "b", 22, false, true), mkEq("b", 26, "b", 26, false, true),
                    mkEq("b", 28, "b", 28, false, true) });
  TEST_EQUAL((int)sitesD.size(), 1)
  TEST_EQUAL(sitesD[0].start_position, 1)
  TEST_EQUAL(sitesD[0].end_position, 22)
  TEST_REAL_SIMILAR(sitesD[0].mass_shift, 526.2196)   // -89.0302 + 615.2498

  // ---- Case E: separated STAYS split ----------------------------------------------------------------
  // Case A's set (b13 present) separates the mods -> they must NOT merge: 2 sites with the gap at 13/14.
  TEST_EQUAL((int)sitesA.size(), 2)
  TEST_EQUAL(sitesA[0].end_position, 13)     // -89 ends at 13
  TEST_EQUAL(sitesA[1].start_position, 14)   // +615 starts at 14 (b13 keeps them apart)

  // ---- Case F: symmetric SUFFIX (y-precursor equiv) -------------------------------------------------
  // Suffix-equiv y<k> covers [L-k+1, L]. Change L seeds [1,L] and tightens ONLY from fragments (the a-priori
  // mod range is used to CLASSIFY, never as an output bound), so a SINGLE suffix ion bounds one side and leaves
  // the other at the seed edge: y21 includer -> lower bound 60, upper stays L=80 -> [60,80]; y20 excluder ->
  // upper bound 60, lower stays 1 -> [1,60]; the two together bound both sides -> localized [60,60].
  auto sitesF1 = FragmentAnalysis::narrowFragmentPTMSites({ PTM{60, 50, 70, 42.0106} }, 80, { mkEq("y", 21, "y", 21, false, true) });
  TEST_EQUAL(sitesF1[0].start_position, 60)
  TEST_EQUAL(sitesF1[0].end_position, 80)    // single includer bounds only the lower side; upper stays seed L
  auto sitesF2 = FragmentAnalysis::narrowFragmentPTMSites({ PTM{60, 50, 70, 42.0106} }, 80, { mkEq("y", 20, "y", 20, false, false) });
  TEST_EQUAL(sitesF2[0].start_position, 1)   // single excluder bounds only the upper side; lower stays seed 1
  TEST_EQUAL(sitesF2[0].end_position, 60)
  auto sitesF3 = FragmentAnalysis::narrowFragmentPTMSites({ PTM{60, 50, 70, 42.0106} }, 80, { mkEq("y", 21, "y", 21, false, true), mkEq("y", 20, "y", 20, false, false) });
  TEST_EQUAL(sitesF3[0].start_position, 60)
  TEST_EQUAL(sitesF3[0].end_position, 60)

  // ---- No-op guards ---------------------------------------------------------------------------------
  PTM heme_np{20, 15, 26, 615.2498};
  auto rL1 = FragmentAnalysis::narrowFragmentPTMSites({ heme_np }, 1, { mkEq("b", 1, "b", 1, false, true) });
  TEST_EQUAL(rL1[0].start_position, 15)      // L<=1: input returned untouched
  TEST_EQUAL(rL1[0].end_position, 26)
  auto rEmpty = FragmentAnalysis::narrowFragmentPTMSites({ heme_np }, 80, {});
  TEST_EQUAL(rEmpty[0].start_position, 1)     // no fragment evidence -> maximal ambiguity over [1,L]
  TEST_EQUAL(rEmpty[0].end_position, 80)

  // ---- Leg 2: leaf verdict == pooled predicate (coveredAmbiguousInEquivFrame) -----------------------
  // (a) flip + N-terminal-loss (b1 / -89): the pooled predicate keeps -89 on the N-terminus (the nterm
  //     exception ignores the flip) -- matching the leaf localizing -89 to [1,1] (case B).
  MS3FragmentMatcher::ProteoformContext ctxNterm;
  ctxNterm.region_start = 0; ctxNterm.region_end = -1;
  ctxNterm.ptm_sites = { PTM{4, 1, 8, -89.0302} };
  bool pooled_b1 = MS3FragmentMatcher::coveredAmbiguousInEquivFrame(ctxNterm, "b", 1, 70, /*sub_includes_ptm=*/true, /*is_flip=*/true) != 0.0;
  auto leafB1 = FragmentAnalysis::narrowFragmentPTMSites({ PTM{4, 1, 8, -89.0302} }, 70, { mkEq("yb", 69, "b", 1, true, true) });
  bool leaf_b1 = (leafB1[0].start_position == 1 && leafB1[0].end_position == 1);
  TEST_EQUAL(leaf_b1, pooled_b1)   // both true: the flipped b1 pins -89 to M1

  // (b) flip-INVERT (b22 / +615, non-nterm): the pooled predicate inverts the flipped sub-verdict so the
  //     equivalent ion carries the mod -- matching the leaf including b22 (case C).
  MS3FragmentMatcher::ProteoformContext ctxHeme;
  ctxHeme.region_start = 0; ctxHeme.region_end = -1;
  ctxHeme.ptm_sites = { PTM{20, 15, 26, 615.2498} };
  bool pooled_b22 = MS3FragmentMatcher::coveredAmbiguousInEquivFrame(ctxHeme, "b", 22, 80, /*sub_includes_ptm=*/false, /*is_flip=*/true) != 0.0;
  TEST_TRUE(pooled_b22)   // flip inverts sub=false -> equiv carries the mod (== leaf includer, case C)

  // Flip-frame classification (the projection property the old leaf ignored): a raw suffix-family ion ("yb")
  // whose equivalent is a prefix ("b") is a complement flip.
  TEST_TRUE(MS3FragmentMatcher::isPrefixIonType("b") != MS3FragmentMatcher::isPrefixIonType("yb"))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
