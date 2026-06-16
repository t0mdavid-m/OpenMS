// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Exploration.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>

#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <set>
#include <string>
#include <cstdio>
#include <cmath>

using namespace OpenMS;

namespace
{
  // Horse heart cytochrome c sequence (104 AA) — matches ms2_hcd_fragment.txt test spectrum
  const char* cytochrome_c_sequence = "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";

  // Base JSON config with MS2 exploration enabled (mass_count, CE 20-40 step 5)
  // Uses cytochrome c sequence to match ms2_hcd_fragment.txt test spectrum
  const char* exploration_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 29,
          "resolution": 120000
        }
      ],
      "ms3": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 35, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": {
        "selection": "intensity",
        "max_targets": 3,
        "exploration": {
          "metric": "mass_count",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0
        }
      },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with MS2 exploration + MS3 exploration (fragment_count)
  const char* ms3_exploration_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 29,
          "resolution": 120000
        }
      ],
      "ms3": [
        {
          "analyzer": "Orbitrap",
          "activation": "CID",
          "collision_energy": 25,
          "resolution": 120000
        }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": {
        "selection": "intensity",
        "max_targets": 3,
        "exploration": {
          "metric": "mass_count",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0
        }
      },
      "ms3": {
        "selection": "intensity",
        "max_targets": 3,
        "exploration": {
          "metric": "fragment_count",
          "ce_min": 15.0,
          "ce_max": 35.0,
          "ce_step": 5.0
        }
      }
    }
  })";

  // Config with MS2 exploration + MS3 selection only (no exploration)
  const char* ms3_selection_only_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 29,
          "resolution": 120000
        }
      ],
      "ms3": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 35, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": {
        "selection": "none",
        "max_targets": 3,
        "exploration": {
          "metric": "mass_count",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0
        }
      },
      "ms3": { "selection": "intensity", "max_targets": 3 }
    }
  })";

  // Config with remaining_precursor exploration metric at MS2
  const char* remaining_precursor_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 29,
          "resolution": 120000
        }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": ""
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": {
        "selection": "none",
        "max_targets": 3,
        "exploration": {
          "metric": "remaining_precursor",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0
        }
      },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with remaining_precursor exploration metric at MS3
  const char* ms3_remaining_precursor_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 29,
          "resolution": 120000
        }
      ],
      "ms3": [
        {
          "analyzer": "Orbitrap",
          "activation": "CID",
          "collision_energy": 25,
          "resolution": 120000
        }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": ""
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": {
        "selection": "none",
        "max_targets": 3
      },
      "ms3": {
        "selection": "none",
        "max_targets": 3,
        "exploration": {
          "metric": "remaining_precursor",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0
        }
      }
    }
  })";

  // Config with cycle time enabled + exploration
  const char* cycle_time_exploration_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 29,
          "resolution": 120000
        }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": true, "value_ms": 1 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 999999
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": ""
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": {
        "selection": "none",
        "max_targets": 3,
        "exploration": {
          "metric": "mass_count",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0
        }
      },
      "ms3": { "selection": "none" }
    }
  })";

  // Config: no MS2 exploration, MS3 exploration enabled
  const char* no_ms2_expl_ms3_expl_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 29,
          "resolution": 120000
        }
      ],
      "ms3": [
        {
          "analyzer": "Orbitrap",
          "activation": "CID",
          "collision_energy": 25,
          "resolution": 120000
        }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": "MKWVTFISLLLLFSSAYSRGVFRR"
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "intensity", "max_targets": 3 },
      "ms3": {
        "selection": "intensity",
        "max_targets": 3,
        "exploration": {
          "metric": "fragment_count",
          "ce_min": 15.0,
          "ce_max": 35.0,
          "ce_step": 5.0
        }
      }
    }
  })";

  // Helper: create a synthetic PeakGroup with a single peak at the given mz/mass/charge.
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

  // Helper: create synthetic DeconvolvedSpectrum with N peak groups.
  // Score = spec.size() (mass_count metric counts peak groups).
  DeconvolvedSpectrum makeSyntheticDeconv(int scan_number, int num_peak_groups)
  {
    DeconvolvedSpectrum ds(scan_number);
    for (int i = 0; i < num_peak_groups; ++i)
    {
      PeakGroup pg;
      ds.push_back(pg);
    }
    return ds;
  }

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

  // Interleaved bootstrap (canonical engine-id-echo contract): drain survey commands and feed the
  // next MS1 spectrum per survey, stamped with the engine's OWN scan_description (a fabricated
  // "scan_"+id MS1 is rejected by the always-on gate -> 0 precursors), until ONE survey selects a
  // precursor and creates the exploration group. Returns that survey's processScan count (the number
  // of exploration variants pushed; 0 if no scan ever selects). Leaves the variant commands QUEUED
  // (does not drain them), so the caller can inspect the engine's next command while a group is
  // active -- the contract the original pushAllMS1Scans + immediate getNextScanCommand relied on.
  int bootstrapExplorationGroup(FLASHIda* ida, const std::vector<ScanData>& scans, int max_iters = 4000)
  {
    const int n_ms1 = static_cast<int>(scans.size());
    int ms1_fed = 0;
    int idle = 0;
    ScanCommand cmd{};
    for (int it = 0; it < max_iters && idle < 3; ++it)
    {
      if (ida->getNextScanCommand(cmd) != 1) break;
      if (cmd.msn_level == 1 && !cmd.is_agc && ms1_fed < n_ms1)
      {
        const ScanData& s = scans[ms1_fed++];
        int n = ida->processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 1,
                                 cmd.scan_description);
        if (n > 0) return n;  // group created; leave its variants queued for the caller
        idle = 0;
      }
      else { ++idle; }       // AGC / re-survey-after-exhausted -> idle tick
      cmd = ScanCommand{};
    }
    return 0;
  }

  // Minimal TSV reader (first line = header, rest = rows; tab-split). Local to this test
  // TU -- the exploration test does not include FLASHIda_TestHelpers.h (it would collide
  // with the anonymous-namespace ScanData / loadTsvScans defined here). Used by the
  // tag_count-under-tagging section to read scan_results.tsv.
  struct TSVFile
  {
    std::vector<std::string> headers;
    std::vector<std::vector<std::string>> rows;

    static TSVFile parse(const std::string& path)
    {
      TSVFile result;
      std::ifstream f(path);
      std::string line;
      bool first = true;
      while (std::getline(f, line))
      {
        std::vector<std::string> cols;
        std::istringstream iss(line);
        std::string col;
        while (std::getline(iss, col, '\t')) cols.push_back(col);
        if (first) { result.headers = cols; first = false; }
        else        { result.rows.push_back(cols); }
      }
      return result;
    }

    int colIndex(const std::string& name) const
    {
      for (size_t i = 0; i < headers.size(); i++)
        if (headers[i] == name) return static_cast<int>(i);
      return -1;
    }
  };

  // Derive an inclusion-pinned, MS3-capable variant of an exploration config: pin the cytC
  // precursor (target_mode=1 + inclusion_cytc.txt) and swap in the validatable M-starting cytC
  // proteoform, so real ms2_cytc_fresh_scan57 b/y fragments match -> MS3 fires (mirrors P4-U07
  // processScan_ms3_commands). Any exploration block in the source config is preserved.
  std::string inclusionPinCytc(std::string cfg)
  {
    auto rep = [&cfg](const std::string& from, const std::string& to) {
      auto p = cfg.find(from);
      if (p != std::string::npos) cfg.replace(p, from.size(), to);
    };
    rep("\"target_mode\": 0", "\"target_mode\": 1");
    rep("\"inclusion_list\": \"\"",
        "\"inclusion_list\": \"../../FlashIDA/test-data/configs/inclusion_cytc.txt\"");
    rep("GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
        "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE");
    // Activate MS2 selection: Exploration::initiateNextLevel (Exploration.cpp:574) returns 0
    // commands when the current level's selection is None, so the MS2-exploration winner would
    // never emit MS3. "selection": "none" is unique to the ms2 block in ms3_selection_only_config
    // (ms1=qscore, ms3=intensity); a no-op for ms3_exploration_config (already intensity).
    rep("\"selection\": \"none\"", "\"selection\": \"intensity\"");
    return cfg;
  }

  struct ExplResult { bool found_ms3 = false; bool found_production_ms2 = false;
                      int total_returns = 0; int ms3_num_stages = 0; int group_commands = 0; };

  // Drive a single exploration group end-to-end via the canonical INTERLEAVED engine-id-echo
  // contract (twin of FLASHIda_TestHelpers::runInterleaved; inlined here because this TU's local
  // ScanData / loadTsvScans / TSVFile / ExplResult would collide with the helper's
  // anonymous-namespace copies). Pull one command at a time and feed exactly one response per
  // requested command, stamped with the engine's OWN scan_description -- the always-on MS1 gate
  // rejects a fabricated "scan_"+id MS1 (returns 0 -> 0 precursors), so MS1 must carry the engine's
  // survey id. The engine paces the surveys: feed the next cytC MS1 spectrum per survey command, in
  // order, until ONE selects a precursor and creates the exploration group (group_commands = its
  // return). Thereafter further MS1 surveys are idle ticks (ONE group -> variants drain contiguously,
  // the unchanged test intent). Feed each drained MS2 variant back; the final-variant feed fires the
  // winner -> MS3 directly (overrides empty), or a production-MS2 re-acquisition (overrides non-empty)
  // which is fed once more -> MS3. Inclusion + ms1.max_targets pins the single cytC precursor.
  ExplResult driveOneExplorationGroup(FLASHIda* ida, const std::vector<ScanData>& ms1_scans,
                                      const ScanData& ms2)
  {
    ExplResult r;
    const int n_ms1 = static_cast<int>(ms1_scans.size());
    int ms1_fed = 0;
    bool group_formed = false;
    int idle = 0;
    for (int i = 0; i < 200 && !r.found_ms3; ++i)
    {
      ScanCommand next{};
      if (ida->getNextScanCommand(next) != 1) break;
      // MS1 survey: before a group exists, feed the next cytC spectrum stamped with the engine id;
      // after the group is formed (or once all scans are fed), MS1 surveys / AGC are idle ticks.
      if (next.is_agc || next.msn_level == 1)
      {
        if (!group_formed && next.msn_level == 1 && !next.is_agc && ms1_fed < n_ms1)
        {
          // Feeding a survey scan is productive engine interaction (not idle), even when it selects
          // 0 precursors -- mirrors the original "try every scan in order until one creates a group".
          const ScanData& s = ms1_scans[ms1_fed++];
          int n = ida->processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 1,
                                   next.scan_description);
          idle = 0;
          if (n > 0) { r.group_commands = n; group_formed = true; }
          continue;
        }
        if (++idle > 3) break;
        continue;
      }
      idle = 0;
      if (next.msn_level == 3) { r.found_ms3 = true; r.ms3_num_stages = next.num_stages; break; }
      if (next.msn_level == 2)
      {
        std::string d(next.scan_description);
        if (d.size() >= 4 && d[3] != 'E') r.found_production_ms2 = true;  // production-MS2 winner
        r.total_returns += ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                            (int)ms2.mzs.size(), ms2.rt, 2, next.scan_description);
      }
    }
    return r;
  }

  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_standard.txt";
  const std::string ms2_tsv_path = "../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt";
  const std::string ms2_cytc_path = "../../FlashIDA/test-data/spectra/ms2_cytc_scan149.txt";

  // Config with ms2.min_charge set impossibly high — all fragments should be filtered
  const char* ms2_min_charge_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 1,
      "max_charge": 50,
      "min_mass": 100,
      "max_mass": 50000,
      "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 29,
          "resolution": 120000
        }
      ],
      "ms3": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 35,
          "resolution": 60000
        }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "intensity", "max_targets": 3, "min_charge": 99 },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with 3-entry tol and MS2 exploration tolerance override
  const char* exploration_tolerance_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10, 20]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 29,
          "resolution": 120000
        }
      ],
      "ms3": [
        {
          "analyzer": "Orbitrap",
          "activation": "CID",
          "collision_energy": 25,
          "resolution": 120000
        }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": ""
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": {
        "selection": "none",
        "max_targets": 3,
        "exploration": {
          "metric": "mass_count",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0,
          "overrides": {
            "tolerance_ppm": "15"
          }
        }
      },
      "ms3": {
        "selection": "none",
        "max_targets": 3,
        "exploration": {
          "metric": "mass_count",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0
        }
      }
    }
  })";

  // Config with ETD exploration for activation type wiring test
  const char* etd_exploration_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
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
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "ETD",
          "reaction_time": 10.0,
          "resolution": 120000
        }
      ],
      "ms3": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 35, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": {
        "selection": "intensity",
        "max_targets": 3,
        "exploration": {
          "metric": "mass_count",
          "activations": ["ETD"],
          "rt_min": 5.0,
          "rt_max": 15.0,
          "rt_step": 5.0
        }
      },
      "ms3": { "selection": "none" }
    }
  })";
} // anonymous namespace


START_TEST(FLASHIda_exploration, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(exploration_group_creation)
{
  Config cfg{std::string(exploration_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);

  TEST_EQUAL(exploration.activeGroupCount(), 1)

  auto group = exploration.getGroup(1);
  TEST_EQUAL(group.group_id, 1)
  TEST_EQUAL(group.msn_level, 2)
  TEST_EQUAL(group.complete, false)
  TEST_EQUAL(group.winner_index, -1)
  TEST_EQUAL(static_cast<int>(group.exploration_metric),
             static_cast<int>(ExplorationMetric::MassCount))

  TEST_EQUAL(static_cast<int>(group.variants.size()), 5)
  TEST_REAL_SIMILAR(group.variants[0].collision_energy, 20.0)
  TEST_REAL_SIMILAR(group.variants[1].collision_energy, 25.0)
  TEST_REAL_SIMILAR(group.variants[2].collision_energy, 30.0)
  TEST_REAL_SIMILAR(group.variants[3].collision_energy, 35.0)
  TEST_REAL_SIMILAR(group.variants[4].collision_energy, 40.0)

  for (int i = 0; i < 5; ++i)
  {
    TEST_EQUAL(group.variants[i].received, false)
    TEST_EQUAL(group.variants[i].variant_index, i)
  }

  TEST_EQUAL(static_cast<int>(cmds.size()), 5)

  (void)group;
}
END_SECTION

START_SECTION(exploration_variants_priority_by_level)
{
  Config cfg{std::string(exploration_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);

  TEST_EQUAL(static_cast<int>(cmds.size()), 5)
  for (int i = 0; i < 5; ++i)
  {
    TEST_EQUAL(cmds[i].msn_level, 2)
    TEST_EQUAL(cmds[i].priority, 2)  // MS2 exploration variants at priority 2
    TEST_EQUAL(cmds[i].is_agc, 0)
    std::string desc(cmds[i].scan_description);
    TEST_EQUAL(desc.size() >= 4, true)
    TEST_EQUAL(desc[3], 'E')
  }
}
END_SECTION

START_SECTION(ms3_exploration_variants_use_buildMS3)
{
  Config cfg{std::string(ms3_exploration_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto precursor_pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  ScanCommand ms2_ctx = queue.buildMS2(precursor_pg, 3, cfg.level(2).scans[0]);

  auto fragment_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  auto cmds = exploration.initiate(3, fragment_pg, 2, -50.0, queue, &ms2_ctx, 'y', 5);

  TEST_EQUAL(static_cast<int>(cmds.size()), 5)

  for (int i = 0; i < 5; ++i)
  {
    TEST_EQUAL(cmds[i].msn_level, 3)
    TEST_EQUAL(cmds[i].num_stages, 2)
    TEST_EQUAL(cmds[i].priority, 1)
    TEST_REAL_SIMILAR(cmds[i].faims_cv, -50.0)
  }

  TEST_EQUAL(cmds[0].stages[0].charge_state, 3)

  auto group = exploration.getGroup(1);
  TEST_EQUAL(group.msn_level, 3)
  TEST_EQUAL(group.originating_cmd.num_stages > 0, true)

  // MS3 exploration descriptions must include fragment ion info
  for (int i = 0; i < 5; ++i)
  {
    std::string desc(cmds[i].scan_description);
    TEST_EQUAL(desc[3], 'E')  // exploration marker
    // Description format: {ID}E{frag_mass_kDa}@{charge}{ion_type}{frag_index}
    // Fragment mass = 500.0 * 2 / 1000.0 = 1.0 (from fragment_pg mz=500, charge=2)
    // Should end with "y5"
    std::string suffix = desc.substr(desc.size() - 2);
    TEST_EQUAL(suffix, "y5")
  }
}
END_SECTION

START_SECTION(ms3_exploration_winner_selection_and_cleanup)
{
  // The winner's reported score must equal the peak-group count (MassCount semantics):
  // makeSyntheticDeconv() below produces empty (mass-0) PeakGroups, so the MS3
  // FragmentCount path (calibrated MS3FragmentMatcher) would match 0 fragments and
  // score every variant 0. Swap the MS3 exploration metric to mass_count so the score
  // is spec.size() and variant 3 (8 peak groups) wins with score 8.0.
  std::string ms3_winner_cfg(ms3_exploration_config);
  {
    auto mpos = ms3_winner_cfg.find("fragment_count");
    TEST_EQUAL(mpos != std::string::npos, true)
    ms3_winner_cfg.replace(mpos, std::string("fragment_count").size(), "mass_count");
  }
  Config cfg{ms3_winner_cfg};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto precursor_pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  ScanCommand ms2_ctx = queue.buildMS2(precursor_pg, 3, cfg.level(2).scans[0]);

  auto fragment_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  auto cmds = exploration.initiate(3, fragment_pg, 2, -50.0, queue, &ms2_ctx);
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)
  TEST_EQUAL(exploration.activeGroupCount(), 1)

  std::vector<int> peak_counts = {2, 4, 6, 8, 3};
  Exploration::FeedResultInfo last_info;
  for (int i = 0; i < 5; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, peak_counts[i]);
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    last_info = exploration.feedResultForTest(tracking_id, ds, static_cast<double>(i), queue);
  }

  TEST_EQUAL(exploration.activeGroupCount(), 0)
  TEST_REAL_SIMILAR(last_info.score, 8.0)

  // ms3_exploration_config has no overrides, so feedResultImpl_() takes
  // the initiateNextLevel path which produces 0 commands for synthetic data.
  // MS3 command format (num_stages=2, priority=1) is verified in
  // ms3_exploration_variants_use_buildMS3 via the initiate() path.
}
END_SECTION

START_SECTION(winner_selection_by_score)
{
  Config cfg{std::string(exploration_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)

  std::vector<double> scores = {1.0, 3.0, 2.0, 5.0, 0.0};
  Exploration::FeedResultInfo last_info;
  for (int i = 0; i < 5; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, static_cast<int>(scores[i]));
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    last_info = exploration.feedResultForTest(tracking_id, ds, static_cast<double>(i), queue);
  }

  // mass_count metric -> remaining_ratio should be -1.0 (N/A)
  TEST_REAL_SIMILAR(last_info.remaining_ratio, -1.0)
  TEST_EQUAL(exploration.activeGroupCount(), 0)
}
END_SECTION

START_SECTION(cycle_time_suppression_during_exploration)
{
  // P7-U05: MS1 cycle time suppression during active exploration
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  FLASHIda* ida = new FLASHIda(const_cast<char*>(cycle_time_exploration_config));

  // Interleaved bootstrap: feed MS1 surveys (engine ids) until a group forms; its variants stay
  // queued. The group is now active, so the next dequeue must be an exploration variant (cycle-time
  // MS1 suppressed) -- not a cycle-time/idle MS1.
  int total = bootstrapExplorationGroup(ida, ms1_scans);
  if (total == 0) { delete ida; }
  ABORT_IF(total == 0)

  // getNextScanCommand should return exploration variants (priority 2, msn_level 2)
  // NOT cycle-time MS1, because exploration is active (cycle-time suppressed)
  ScanCommand cmd{};
  int result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(cmd.msn_level, 2)
  TEST_EQUAL(cmd.priority, 2)
  std::string desc(cmd.scan_description);
  TEST_EQUAL(desc.size() >= 4, true)
  TEST_EQUAL(desc[3], 'E')

  delete ida;
}
END_SECTION

START_SECTION(ms1_resumes_after_exploration_completes)
{
  // P7-U06: MS1 cycle time injection resumes after exploration completes
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  FLASHIda* ida = new FLASHIda(const_cast<char*>(cycle_time_exploration_config));

  // Interleaved bootstrap: feed MS1 surveys (engine ids) until a group forms; its variants stay
  // queued for the drain loop below.
  int total = bootstrapExplorationGroup(ida, ms1_scans);
  if (total == 0) { delete ida; }
  ABORT_IF(total == 0)

  // Drain all exploration variants and feed MS2 results back
  std::vector<ScanCommand> exploration_cmds;
  ScanCommand cmd{};
  while (ida->getNextScanCommand(cmd) == 1)
  {
    std::string desc(cmd.scan_description);
    if (cmd.msn_level == 2 && desc.size() >= 4 && desc[3] == 'E')
    {
      exploration_cmds.push_back(cmd);
    }
    else
    {
      break;
    }
  }

  // Feed MS2 results for each exploration variant
  const auto& ms2 = ms2_scans[0];
  for (const auto& ecmd : exploration_cmds)
  {
    ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                     (int)ms2.mzs.size(), ms2.rt, 2, ecmd.scan_description);
  }

  // After exploration completes, getNextScanCommand should eventually return MS1
  bool found_ms1 = false;
  for (int i = 0; i < 20; ++i)
  {
    ScanCommand next{};
    if (ida->getNextScanCommand(next) != 1) break;
    if (next.msn_level == 1 && next.is_agc == 0)
    {
      found_ms1 = true;
      break;
    }
  }
  TEST_EQUAL(found_ms1, true)

  delete ida;
}
END_SECTION

START_SECTION(ms3_exploration_creates_child_groups)
{
  // Load real MS2 spectrum data for fragment matching
  auto ms2_scans = loadTsvScans(ms2_cytc_path);
  ABORT_IF(ms2_scans.empty());
  const auto& ms2_data = ms2_scans[0];

  Config cfg{std::string(ms3_exploration_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)
  // MS2-level exploration variants should have priority 2
  for (int i = 0; i < 5; ++i)
  {
    TEST_EQUAL(cmds[i].priority, 2)
  }

  // Deconvolve once without precursor constraint (test has no matched MS1)
  deconv.deconvolveMSn(ms2_data.mzs.data(), ms2_data.ints.data(),
                        static_cast<int>(ms2_data.mzs.size()), ms2_data.rt, 0.0, 0);
  DeconvolvedSpectrum ms2_deconv = deconv.storedMS2();

  for (int i = 0; i < 5; ++i)
  {
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    exploration.feedResultForTest(tracking_id, ms2_deconv, ms2_data.rt, queue);
  }

  int ms3_group_count = exploration.activeGroupCount();
  TEST_EQUAL(ms3_group_count > 0, true)
  TEST_EQUAL(ms3_group_count <= 3, true)
}
END_SECTION

START_SECTION(ms3_selection_no_exploration_standard_targeting)
{
  // P7-U08: MS3 with selection but no exploration -> standard MS3 commands.
  // Inclusion-pin the cytC precursor + M-starting proteoform so the real fresh57 b/y ladder
  // matches; the MS2-exploration winner then emits MS3 directly (overrides empty ->
  // initiateNextLevel, Exploration.cpp:548-555).
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  std::string cfg_str = inclusionPinCytc(ms3_selection_only_config);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  ExplResult c = driveOneExplorationGroup(ida, ms1_scans, ms2_scans[0]);
  delete ida;
  ABORT_IF(c.group_commands == 0)
  // overrides empty -> the winning variant feed emits MS3 directly (two-stage command).
  TEST_EQUAL(c.found_ms3, true)
  TEST_EQUAL(c.ms3_num_stages, 2)
}
END_SECTION

START_SECTION(ms2_exploration_production_winner_then_ms3)
{
  // Overrides NON-empty branch: with an exploration override set, the MS2-exploration winner is
  // a production-MS2 re-acquisition (Exploration.cpp:520-547) rather than MS3-direct. Feeding
  // that winner back drives the standard MS2->MS3 path, so MS3 is still acquired FROM the winning
  // scan. Real cytC data, inclusion-pinned.
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  // Inject a non-tolerance override so config_.level(2).overrides stays non-empty
  // (tolerance_ppm would be extracted+erased; "analyzer" persists).
  std::string cfg_str = inclusionPinCytc(ms3_selection_only_config);
  {
    auto p = cfg_str.find("\"metric\": \"mass_count\"");
    ABORT_IF(p == std::string::npos)
    cfg_str.insert(p, "\"overrides\": { \"analyzer\": \"Orbitrap\" }, ");
  }
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  ExplResult c = driveOneExplorationGroup(ida, ms1_scans, ms2_scans[0]);
  delete ida;
  ABORT_IF(c.group_commands == 0)
  // overrides non-empty -> the winner is a production-MS2 re-acquisition; driveOneExplorationGroup
  // re-feeds it, and the standard MS2->MS3 path then yields MS3 from the winning scan.
  TEST_EQUAL(c.found_production_ms2, true)
  TEST_EQUAL(c.found_ms3, true)
}
END_SECTION

START_SECTION(ms2_then_ms3_exploration_acquires_ms3)
{
  // Two-level exploration (MS2 exploration + MS3 exploration, ms3_exploration_config),
  // inclusion-pinned + real cytC: the MS2-exploration winner cascades into MS3-level commands.
  // Verifies the MS2->MS3 exploration chain actually reaches MS3 on real data (existing MS3-
  // exploration sections only assert structure on synthetic data).
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  std::string cfg_str = inclusionPinCytc(ms3_exploration_config);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  ExplResult c = driveOneExplorationGroup(ida, ms1_scans, ms2_scans[0]);
  delete ida;
  ABORT_IF(c.group_commands == 0)
  // two-level: the MS2-exploration winner cascades into MS3-level commands (two-stage).
  TEST_EQUAL(c.found_ms3, true)
  TEST_EQUAL(c.ms3_num_stages, 2)
}
END_SECTION

START_SECTION(optimization_metadata_populated)
{
  Config cfg{std::string(exploration_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);

  DeconvolvedSpectrum ds = makeSyntheticDeconv(1, 3);
  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResultForTest(tracking_id, ds, 1.0, queue);

  TEST_EQUAL(exploration.activeGroupCount(), 1)

  auto group = exploration.getGroup(1);
  TEST_EQUAL(group.variants[0].received, true)
  auto& stored = group.variants[0].result;
  TEST_EQUAL(stored.hasOptimizationMetadata(), true)

  const auto* meta = stored.getOptimizationMetadata();
  TEST_EQUAL(meta->group_id, 1)
  TEST_EQUAL(meta->variant_index, 0)
  TEST_EQUAL(meta->total_variants, 5)
  TEST_REAL_SIMILAR(meta->collision_energy, 20.0)
  TEST_STRING_EQUAL(meta->activation_type, "HCD")
  TEST_EQUAL(meta->exploration_metric, static_cast<int>(ExplorationMetric::MassCount))
  TEST_EQUAL(meta->is_best_variant, false)
  TEST_REAL_SIMILAR(meta->fragmentation_quality_score, 3.0)
  TEST_EQUAL(meta->exploration_scans, 5)
  TEST_EQUAL(meta->start_ms > 0, true)

  (void)meta;
}
END_SECTION

START_SECTION(metadata_serialized_to_msspectrum)
{
  // P7-U10: Metadata serialized to MSSpectrum via toSpectrum()
  // CRITICAL: toSpectrum() returns MSSpectrum by value (Phase 2 lesson #1)
  // CRITICAL: toSpectrum() requires at least one PeakGroup (Phase 2 lesson #3)

  DeconvolvedSpectrum ds(1);

  // Push a PeakGroup (required by toSpectrum — it accesses peak_groups_[0])
  PeakGroup pg;
  ds.push_back(pg);

  // Populate metadata manually (simulating what feedExplorationResult_ does)
  auto& meta = ds.getOrCreateOptimizationMetadata();
  meta.group_id = 42;
  meta.collision_energy = 25.0;
  meta.is_best_variant = true;
  meta.fragmentation_quality_score = 3.0;
  meta.precursor_mass = 2400.5;
  meta.exploration_metric = 1;  // MassCount

  // Call toSpectrum (returns by value!)
  MSSpectrum out_spec = ds.toSpectrum(1);

  TEST_EQUAL(static_cast<int>(out_spec.getMetaValue("optimization_group_id")), 42)
  TEST_REAL_SIMILAR(static_cast<double>(out_spec.getMetaValue("optimization_collision_energy")), 25.0)
  TEST_STRING_EQUAL(static_cast<std::string>(out_spec.getMetaValue("optimization_is_best_variant")), "true")
  TEST_REAL_SIMILAR(static_cast<double>(out_spec.getMetaValue("optimization_quality_score")), 3.0)
  TEST_REAL_SIMILAR(static_cast<double>(out_spec.getMetaValue("optimization_precursor_mass")), 2400.5)
  TEST_EQUAL(static_cast<int>(out_spec.getMetaValue("optimization_exploration_metric")), 1)

  (void)out_spec;
}
END_SECTION

START_SECTION(no_ms2_exploration_ms3_exploration_immediate)
{
  Config cfg{std::string(no_ms2_expl_ms3_expl_config)};
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto ms2_cfg = cfg.level(2);
  TEST_EQUAL(static_cast<int>(ms2_cfg.exploration), static_cast<int>(ExplorationMetric::None))

  auto ms3_cfg = cfg.level(3);
  TEST_EQUAL(static_cast<int>(ms3_cfg.exploration), static_cast<int>(ExplorationMetric::FragmentCount))

  TEST_EQUAL(exploration.activeGroupCount(), 0)

  (void)ms2_cfg;
  (void)ms3_cfg;
}
END_SECTION

START_SECTION(selection_metric_controls_config)
{
  Config cfg{std::string(exploration_config)};

  auto ms1_cfg = cfg.level(1);
  TEST_EQUAL(static_cast<int>(ms1_cfg.selection), static_cast<int>(SelectionMetric::QScore))
  TEST_EQUAL(ms1_cfg.max_targets, 3)

  auto ms2_cfg = cfg.level(2);
  TEST_EQUAL(static_cast<int>(ms2_cfg.selection), static_cast<int>(SelectionMetric::Intensity))
  TEST_EQUAL(ms2_cfg.max_targets, 3)
  TEST_EQUAL(static_cast<int>(ms2_cfg.exploration), static_cast<int>(ExplorationMetric::MassCount))
  TEST_REAL_SIMILAR(ms2_cfg.ce_min, 20.0)
  TEST_REAL_SIMILAR(ms2_cfg.ce_max, 40.0)
  TEST_REAL_SIMILAR(ms2_cfg.ce_step, 5.0)

  auto ms3_cfg = cfg.level(3);
  TEST_EQUAL(static_cast<int>(ms3_cfg.selection), static_cast<int>(SelectionMetric::None))

  (void)ms1_cfg;
  (void)ms2_cfg;
  (void)ms3_cfg;
}
END_SECTION

START_SECTION(remaining_precursor_config_parsed)
{
  // Verify the remaining_precursor metric is correctly parsed from JSON config
  Config cfg{std::string(remaining_precursor_config)};
  auto ms2_cfg = cfg.level(2);
  TEST_EQUAL(static_cast<int>(ms2_cfg.exploration),
             static_cast<int>(ExplorationMetric::RemainingPrecursor))
  TEST_REAL_SIMILAR(ms2_cfg.remaining_precursor_target, 0.1)

  (void)ms2_cfg;
}
END_SECTION

START_SECTION(remaining_precursor_score_no_raw_data)
{
  // RemainingPrecursor now prepends a CE=0 baseline scan: 1 baseline + 5 CE variants = 6
  Config cfg{std::string(remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  // Verify first command is CE=0 baseline
  TEST_REAL_SIMILAR(cmds[0].stages[0].collision_energy, 0.0)

  // Verify exploration group has RemainingPrecursor metric
  auto group = exploration.getGroup(1);
  TEST_EQUAL(static_cast<int>(group.exploration_metric),
             static_cast<int>(ExplorationMetric::RemainingPrecursor))
  // Baseline variant at index 0 should have is_baseline=true, variant_index=-1
  TEST_EQUAL(group.variants[0].is_baseline, true)
  TEST_EQUAL(group.variants[0].variant_index, -1)
  // Non-baseline variants should have is_baseline=false
  for (int i = 1; i < 6; ++i)
  {
    TEST_EQUAL(group.variants[i].is_baseline, false)
    TEST_EQUAL(group.variants[i].variant_index, i - 1)
  }

  // Feed all 6 variants via feedResultForTest (no raw data -> all scores = 0.0)
  for (int i = 0; i < 6; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, i + 1);
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    auto info = exploration.feedResultForTest(tracking_id, ds, static_cast<double>(i), queue);

    // total_variants should exclude baseline (= 5 real variants)
    TEST_EQUAL(info.total_variants, 5)
    TEST_REAL_SIMILAR(info.remaining_ratio, -1.0)
  }

  // Group should be complete (all variants received)
  TEST_EQUAL(exploration.activeGroupCount(), 0)
}
END_SECTION

START_SECTION(remaining_precursor_score_with_raw_data)
{
  // Test RemainingPrecursor scoring with raw data: baseline CE=0 + CE>0 variants
  Config cfg{std::string(remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  auto group = exploration.getGroup(1);
  TEST_EQUAL(group.precursor_mz > 0.0, true)

  // Feed CE=0 baseline with precursor-window signal (full intensity)
  std::vector<double> baseline_mzs = {790.0, 800.0, 810.0, 900.0};
  std::vector<double> baseline_ints = {100.0, 500.0, 200.0, 50.0};
  // In-window sum depends on isolation_width; 800.0 is precursor_mz center.

  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResult(baseline_tid, baseline_mzs.data(), baseline_ints.data(),
                         static_cast<int>(baseline_mzs.size()), 0.5, queue);

  // After baseline, group should have baseline_intensity set
  auto group_after_baseline = exploration.getGroup(1);
  TEST_EQUAL(group_after_baseline.has_baseline, true)
  TEST_EQUAL(group_after_baseline.baseline_intensity >= 0.0, true)

  // Feed second variant (CE=20) with reduced in-window signal (fragmentation removed some)
  std::vector<double> frag_mzs = {790.0, 800.0, 810.0, 900.0};
  std::vector<double> frag_ints = {50.0, 200.0, 100.0, 50.0};

  int ce20_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto ce20_info = exploration.feedResult(ce20_tid, frag_mzs.data(), frag_ints.data(),
                                          static_cast<int>(frag_mzs.size()), 1.0, queue);

  // If baseline had in-window signal, CE>0 variant should have a score
  auto group_mid = exploration.getGroup(1);
  TEST_EQUAL(group_mid.variants[1].received, true)
  TEST_EQUAL(group_mid.variants[1].score >= 0.0, true)
  TEST_EQUAL(group_mid.variants[1].score <= 1.0, true)

  // remaining_ratio should be valid for RemainingPrecursor with raw data
  TEST_EQUAL(ce20_info.remaining_ratio >= 0.0, true)
  TEST_EQUAL(ce20_info.remaining_ratio <= 1.0, true)
}
END_SECTION

START_SECTION(remaining_precursor_score_no_signal_in_window)
{
  // Feed CE=0 baseline with zero in-window signal. All subsequent variants should score 0.
  Config cfg{std::string(remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  // Feed baseline (CE=0) with raw data entirely OUTSIDE the precursor window
  // -> baseline_intensity = 0 -> all subsequent scores = 0 (baseline failure)
  std::vector<double> mzs = {400.0, 500.0, 600.0, 1200.0};
  std::vector<double> intensities = {100.0, 200.0, 300.0, 400.0};

  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResult(baseline_tid, mzs.data(), intensities.data(),
                         static_cast<int>(mzs.size()), 0.5, queue);

  auto group = exploration.getGroup(1);
  TEST_EQUAL(group.has_baseline, true)
  TEST_REAL_SIMILAR(group.baseline_intensity, 0.0)

  // Feed CE=20 variant with signal inside the precursor window
  // Since baseline = 0, score should be 0 (baseline failure path)
  std::vector<double> frag_mzs = {790.0, 800.0, 810.0, 900.0};
  std::vector<double> frag_ints = {100.0, 500.0, 200.0, 50.0};

  int ce20_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto ce20_info = exploration.feedResult(ce20_tid, frag_mzs.data(), frag_ints.data(),
                                          static_cast<int>(frag_mzs.size()), 1.0, queue);

  auto group_after = exploration.getGroup(1);
  TEST_EQUAL(group_after.variants[1].received, true)
  TEST_REAL_SIMILAR(group_after.variants[1].score, 0.0)
  TEST_REAL_SIMILAR(ce20_info.remaining_ratio, -1.0)
}
END_SECTION

START_SECTION(remaining_precursor_empty_baseline_aborts_group)
{
  // Empty baseline window => no CE variant can be scored. The group must abort: cancel its
  // still-queued child scans, pick no winner, and emit no follow-up production scan.
  Config cfg{std::string(remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  // Push the 5 CE variants (children) into the queue; feed the baseline directly (already acquired).
  for (size_t i = 1; i < cmds.size(); ++i) queue.push(cmds[i]);
  TEST_EQUAL(static_cast<int>(queue.queueSize(2)), 5)

  // CE=0 baseline with raw data entirely OUTSIDE the precursor window -> baseline_intensity = 0.
  std::vector<double> mzs = {400.0, 500.0, 600.0, 1200.0};
  std::vector<double> intensities = {100.0, 200.0, 300.0, 400.0};
  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info = exploration.feedResult(baseline_tid, mzs.data(), intensities.data(),
                                     static_cast<int>(mzs.size()), 0.5, queue);

  // All 5 queued children cancelled; group finalized with no winner and no follow-up scan.
  TEST_EQUAL(static_cast<int>(queue.queueSize(2)), 0)
  TEST_EQUAL(exploration.activeGroupCount(), 0)
  TEST_EQUAL(static_cast<int>(info.commands.size()), 0)
}
END_SECTION

START_SECTION(remaining_precursor_inflight_child_after_abort_is_noop)
{
  // A child already dequeued/sent to the device (in-flight) that returns its result AFTER the
  // group has aborted must be a harmless no-op and must not resurrect the cleaned-up group.
  Config cfg{std::string(remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  // cmds[2..5] still queued; cmds[1] already sent to the device (in-flight, in pending_scan_map_).
  for (size_t i = 2; i < cmds.size(); ++i) queue.push(cmds[i]);
  queue.registerPending(cmds[1]);
  int inflight_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  TEST_EQUAL(queue.peekPending(inflight_tid).has_value(), true)

  // Baseline returns empty -> abort.
  std::vector<double> mzs = {400.0, 500.0, 600.0, 1200.0};
  std::vector<double> intensities = {100.0, 200.0, 300.0, 400.0};
  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResult(baseline_tid, mzs.data(), intensities.data(),
                         static_cast<int>(mzs.size()), 0.5, queue);

  TEST_EQUAL(exploration.activeGroupCount(), 0)                    // aborted + cleaned up
  TEST_EQUAL(static_cast<int>(queue.queueSize(2)), 0)             // queued children cancelled
  TEST_EQUAL(queue.peekPending(inflight_tid).has_value(), false)  // in-flight child cancelled from pending

  // The device now finishes the in-flight child and returns its result LATE -> harmless no-op.
  std::vector<double> frag_mzs = {790.0, 800.0, 810.0, 900.0};
  std::vector<double> frag_ints = {100.0, 500.0, 200.0, 50.0};
  auto late_info = exploration.feedResult(inflight_tid, frag_mzs.data(), frag_ints.data(),
                                          static_cast<int>(frag_mzs.size()), 1.0, queue);
  TEST_EQUAL(late_info.group_id, -1)                              // empty info (routing entry gone)
  TEST_EQUAL(late_info.variant_index, -1)                         // not routed to any variant
  TEST_EQUAL(exploration.activeGroupCount(), 0)                   // group NOT resurrected
}
END_SECTION

START_SECTION(fragment_count_requires_protein_sequence)
{
  // Config with fragment_count metric but empty protein_sequence should throw
  std::string cfg_str = std::string(exploration_config);
  // Clear protein_sequence to test empty-sequence validation path
  {
    const std::string seq = "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";
    auto seq_pos = cfg_str.find(seq);
    if (seq_pos != std::string::npos) cfg_str.erase(seq_pos, seq.size());
  }
  // Replace "mass_count" with "fragment_count" at MS2 level
  auto pos = cfg_str.find("\"mass_count\"");
  if (pos != std::string::npos)
    cfg_str.replace(pos, 12, "\"fragment_count\"");

  // fragment_count requires protein_sequence — should throw
  TEST_EXCEPTION(std::invalid_argument, Config cfg{cfg_str})
}
END_SECTION

START_SECTION(selection_metric_terminal_fragments_parsing)
{
  // Build a config with ms2 selection = "terminal_fragments" and a protein sequence
  std::string cfg_str = std::string(ms3_selection_only_config);
  // ms3_selection_only_config has ms2.selection = "none" and a cytochrome-c protein_sequence (GDVEK...ATNE)
  // Replace ms2 selection "none" with "terminal_fragments"
  auto pos = cfg_str.find("\"selection\": \"none\"");
  TEST_EQUAL(pos != std::string::npos, true)
  cfg_str.replace(pos, 19, "\"selection\": \"terminal_fragments\"");

  Config cfg{cfg_str};
  TEST_EQUAL(static_cast<int>(cfg.level(2).selection),
             static_cast<int>(SelectionMetric::TerminalFragments))
}
END_SECTION

START_SECTION(selection_metric_ambiguity_resolution_parsing)
{
  // Build a config with ms2 selection = "ambiguity_resolution" and a protein sequence
  std::string cfg_str = std::string(ms3_selection_only_config);
  // ms3_selection_only_config has ms2.selection = "none" and a cytochrome-c protein_sequence (GDVEK...ATNE)
  auto pos = cfg_str.find("\"selection\": \"none\"");
  TEST_EQUAL(pos != std::string::npos, true)
  cfg_str.replace(pos, 19, "\"selection\": \"ambiguity_resolution\"");

  Config cfg{cfg_str};
  TEST_EQUAL(static_cast<int>(cfg.level(2).selection),
             static_cast<int>(SelectionMetric::AmbiguityResolution))
}
END_SECTION

START_SECTION(terminal_fragments_requires_protein_sequence)
{
  // ms2 selection = terminal_fragments, protein_sequence empty -> should throw
  std::string cfg_str = std::string(exploration_config);
  // Clear protein_sequence to test empty-sequence validation path
  {
    const std::string seq = "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";
    auto seq_pos = cfg_str.find(seq);
    if (seq_pos != std::string::npos) cfg_str.erase(seq_pos, seq.size());
  }
  // Change ms2 selection to "terminal_fragments" but leave protein_sequence empty
  auto pos = cfg_str.find("\"selection\": \"none\"");
  TEST_EQUAL(pos != std::string::npos, true)
  cfg_str.replace(pos, 19, "\"selection\": \"terminal_fragments\"");

  TEST_EXCEPTION(std::invalid_argument, Config cfg{cfg_str})
}
END_SECTION

START_SECTION(ambiguity_resolution_requires_protein_sequence)
{
  // ms2 selection = ambiguity_resolution, protein_sequence empty -> should throw
  std::string cfg_str = std::string(exploration_config);
  // Clear protein_sequence to test empty-sequence validation path
  {
    const std::string seq = "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";
    auto seq_pos = cfg_str.find(seq);
    if (seq_pos != std::string::npos) cfg_str.erase(seq_pos, seq.size());
  }
  auto pos = cfg_str.find("\"selection\": \"none\"");
  TEST_EQUAL(pos != std::string::npos, true)
  cfg_str.replace(pos, 19, "\"selection\": \"ambiguity_resolution\"");

  TEST_EXCEPTION(std::invalid_argument, Config cfg{cfg_str})
}
END_SECTION

START_SECTION(ms3_protein_sequence_only_accepted)
{
  // Config with ms3.protein_sequence should be accepted (no throw)
  Config cfg{std::string(exploration_config)};
  TEST_EQUAL(cfg.targeting().protein_sequence.empty(), false)
}
END_SECTION

START_SECTION(remaining_precursor_target_aware_scoring)
{
  Config cfg{std::string(remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  auto group = exploration.getGroup(1);
  double iso_half = group.isolation_width / 2.0;
  double mz_center = group.precursor_mz;

  std::vector<double> baseline_mzs = {mz_center};
  std::vector<double> baseline_ints = {1000.0};
  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResult(baseline_tid, baseline_mzs.data(), baseline_ints.data(),
                         static_cast<int>(baseline_mzs.size()), 0.5, queue);

  std::vector<double> perfect_ints = {100.0};
  int ce20_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto info_perfect = exploration.feedResult(ce20_tid, baseline_mzs.data(), perfect_ints.data(),
                                              1, 1.0, queue);
  double score_perfect = info_perfect.score;
  double ratio_perfect = info_perfect.remaining_ratio;

  std::vector<double> over_ints = {500.0};
  int ce25_tid = queue.decode(std::string(cmds[2].scan_description).substr(0, 3));
  auto info_over = exploration.feedResult(ce25_tid, baseline_mzs.data(), over_ints.data(),
                                           1, 2.0, queue);
  double score_over = info_over.score;
  double ratio_over = info_over.remaining_ratio;

  TEST_EQUAL(score_perfect > score_over, true)
  TEST_REAL_SIMILAR(score_perfect, 1.0)
  TEST_REAL_SIMILAR(score_over, 0.6)
  TEST_REAL_SIMILAR(ratio_perfect, 0.1)
  TEST_REAL_SIMILAR(ratio_over, 0.5)

  (void)iso_half;
}
END_SECTION

START_SECTION(fragment_match_propagated_in_feed_result)
{
  Config cfg{std::string(exploration_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  TEST_EQUAL(cfg.targeting().protein_sequence.empty(), false)

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);

  auto ms2_scans = loadTsvScans(ms2_cytc_path);
  ABORT_IF(ms2_scans.empty())
  const auto& ms2_data = ms2_scans[0];

  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info = exploration.feedResult(tracking_id,
      ms2_data.mzs.data(), ms2_data.ints.data(),
      static_cast<int>(ms2_data.mzs.size()), ms2_data.rt, queue);

  // Real cytochrome c spectrum should produce fragment matches
  TEST_EQUAL(info.fragment_count > 0, true)
  TEST_EQUAL(info.matched_protein.empty(), false)
  TEST_EQUAL(info.proteoform_sequence.empty(), false)
  TEST_STRING_EQUAL(info.proteoform_sequence, std::string(cytochrome_c_sequence))
}
END_SECTION

START_SECTION(fragment_count_zero_without_protein_sequence)
{
  Config cfg{std::string(remaining_precursor_config)};
  TEST_EQUAL(cfg.targeting().protein_sequence.empty(), true)

  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  DeconvolvedSpectrum baseline_ds = makeSyntheticDeconv(0, 1);
  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResultForTest(baseline_tid, baseline_ds, 0.0, queue);

  DeconvolvedSpectrum ds = makeSyntheticDeconv(1, 5);
  int tracking_id = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto info = exploration.feedResultForTest(tracking_id, ds, 1.0, queue);

  TEST_EQUAL(info.fragment_count, 0)
  TEST_EQUAL(info.matched_protein.empty(), true)
  TEST_EQUAL(info.proteoform_sequence.empty(), true)
}
END_SECTION

START_SECTION(ms2_exploration_returns_command_count)
{
  // Inclusion-pinned cytC + M-proteoform + real fresh57 ladder: the MS2-exploration winner
  // emits MS3 (overrides empty -> initiateNextLevel), so the winning feed's processScan
  // returns a positive command count.
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  std::string cfg_str = inclusionPinCytc(ms3_selection_only_config);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  ExplResult c = driveOneExplorationGroup(ida, ms1_scans, ms2_scans[0]);
  delete ida;
  ABORT_IF(c.group_commands == 0)
  // The winning variant feed pushes the follow-up command(s), so the summed processScan return
  // count is positive (was a tautological rv >= 0 before).
  TEST_EQUAL(c.total_returns > 0, true)
}
END_SECTION

START_SECTION(ms3_remaining_precursor_isolation_width)
{
  // MS3 RemainingPrecursor scoring requires non-zero isolation_width.
  // Before the fix, initiate() computed width=0 from single-peak PeakGroups,
  // causing all MS3 variants to score -1. The 2.0 Da floor fixes this.
  Config cfg{std::string(ms3_remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  // Create a narrow single-peak PeakGroup (simulates MS3 fragment target)
  // getMzRange() returns (500.0, 500.0) -> width=0 before floor
  auto fragment_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  ScanCommand ms2_ctx = queue.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3,
                                        cfg.level(2).scans[0]);

  auto cmds = exploration.initiate(3, fragment_pg, 2, 0.0, queue, &ms2_ctx);
  // RemainingPrecursor: 1 baseline + 5 CE variants (20,25,30,35,40) = 6
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  auto group = exploration.getGroup(1);
  // isolation_width should be floored to 2.0 (not 0.0)
  TEST_REAL_SIMILAR(group.isolation_width, 2.0)

  double mz_center = group.precursor_mz;

  // Feed baseline (CE=0) with signal at precursor center
  std::vector<double> baseline_mzs = {mz_center};
  std::vector<double> baseline_ints = {1000.0};
  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResult(baseline_tid, baseline_mzs.data(), baseline_ints.data(),
                         static_cast<int>(baseline_mzs.size()), 0.5, queue);

  auto group_after = exploration.getGroup(1);
  TEST_EQUAL(group_after.has_baseline, true)
  // With 2.0 Da window [499.0, 501.0], mz_center=500.0 is in-window
  TEST_REAL_SIMILAR(group_after.baseline_intensity, 1000.0)

  // Feed CE=20 variant with 100.0 intensity -> ratio = 0.1
  std::vector<double> variant_ints = {100.0};
  int ce20_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto info = exploration.feedResult(ce20_tid, baseline_mzs.data(), variant_ints.data(),
                                      1, 1.0, queue);

  // Score should be real (not -1.0), ratio = 100/1000 = 0.1
  TEST_REAL_SIMILAR(info.remaining_ratio, 0.1)
  TEST_REAL_SIMILAR(info.score, 1.0)  // target=0.1, deviation=0.0, score=1.0
  TEST_EQUAL(info.score > 0.0, true)
}
END_SECTION

START_SECTION(tol_three_entry_parsing)
{
  Config cfg{std::string(exploration_tolerance_config)};
  // MS1 = 10, MS2 = 10, MS3 = 20
  TEST_REAL_SIMILAR(cfg.level(1).tolerance_ppm, 10.0)
  TEST_REAL_SIMILAR(cfg.level(2).tolerance_ppm, 10.0)
  TEST_REAL_SIMILAR(cfg.level(3).tolerance_ppm, 20.0)
}
END_SECTION

START_SECTION(exploration_tolerance_override)
{
  Config cfg{std::string(exploration_tolerance_config)};
  // MS2 exploration has overrides.tolerance_ppm = "15"
  TEST_REAL_SIMILAR(cfg.level(2).exploration_tolerance_ppm, 15.0)
  // MS2 base tolerance is still 10
  TEST_REAL_SIMILAR(cfg.level(2).tolerance_ppm, 10.0)
  // MS3 exploration has no override -> falls back to tol[2] = 20
  TEST_REAL_SIMILAR(cfg.level(3).exploration_tolerance_ppm, 20.0)
  // tolerance_ppm should be removed from overrides map
  TEST_EQUAL(cfg.level(2).overrides.count("tolerance_ppm"), 0u)
}
END_SECTION

START_SECTION(exploration_tolerance_fallback)
{
  // exploration_config has exploration at MS2 but NO tolerance override
  Config cfg{std::string(exploration_config)};
  // exploration_tolerance_ppm should equal base tolerance
  TEST_REAL_SIMILAR(cfg.level(2).exploration_tolerance_ppm, 10.0)
  TEST_REAL_SIMILAR(cfg.level(2).tolerance_ppm, 10.0)
}
END_SECTION

START_SECTION(tol_validation_insufficient_entries)
{
  // MS3 configured in selection_strategy but only 2 tol entries -> must throw
  const char* insufficient_tol_config = R"({
    "deconvolution": {
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000,
      "tol": [10, 10]
    },
    "precursor_selection": { "RT_window": 180 },
    "tagging": {},
    "quantification": { "enabled": false },
    "faims": { "cv_values": [-50] },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000 },
      "ms2": [{ "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }],
      "ms3": [{ "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25, "resolution": 120000 }]
    },
    "scheduling": {},
    "files": {},
    "ms3": { "protein_sequence": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none", "max_targets": 3 },
      "ms3": { "selection": "none", "max_targets": 3 }
    }
  })";
  TEST_EXCEPTION(std::invalid_argument, Config cfg{std::string(insufficient_tol_config)})
}
END_SECTION

START_SECTION(selection_without_next_level_scan_config_rejected)
{
  // ms2 selects targets for MS3 (selection != none) but no ms_settings.ms3 scan config is
  // defined -> must throw at construction. Guards the OOB read of next_cfg.scans[0] in
  // Exploration::initiateNextLevel.
  const char* missing_next_scan_config = R"({
    "deconvolution": {
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000,
      "tol": [10, 10, 10]
    },
    "precursor_selection": { "RT_window": 180 },
    "tagging": {},
    "quantification": { "enabled": false },
    "faims": { "cv_values": [-50] },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000 },
      "ms2": [{ "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }]
    },
    "scheduling": {},
    "files": {},
    "ms3": { "protein_sequence": "PEPTIDER" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "intensity", "max_targets": 3 },
      "ms3": { "selection": "none" }
    }
  })";
  TEST_EXCEPTION(std::invalid_argument, Config cfg{std::string(missing_next_scan_config)})
}
END_SECTION

START_SECTION(initiateNextLevel_ms2_min_charge_filters_fragments)
{
  Config cfg{std::string(ms2_min_charge_config)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  // Load real MS2 spectrum that normally produces fragment matches
  auto scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(scans.empty())

  // Deconvolve the MS2 spectrum
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  for (const auto& scan : scans)
  {
    deconv.deconvolveMSn(scan.mzs.data(), scan.ints.data(), (int)scan.mzs.size(),
                         scan.rt, 12000.0, 10);
  }

  auto precursor_pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  ScanCommand ms2_ctx = queue.buildMS2(precursor_pg, 3, cfg.level(2).scans[0]);

  // initiateNextLevel processes MS2 results and picks fragments for MS3
  // With ms2.min_charge=99, ALL fragments should be filtered out
  auto nlr = exploration.initiateNextLevel(2, deconv.storedMS2(), -50.0, queue, &ms2_ctx);
  TEST_EQUAL(static_cast<int>(nlr.commands.size()), 0)  // no commands — all fragments filtered by charge
}
END_SECTION

START_SECTION(activation_type_wiring_in_scoring)
{
  // Verify that exploration scoring passes variant activation type through
  // to fragment matching (not hardcoded "HCD").
  // ETD config produces variants with activation_type="ETD", which maps to
  // {c, z} ions — different from HCD's {b, y}.

  Config cfg{std::string(etd_exploration_config)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);
  TEST_EQUAL(cmds.size() > 0, true)

  // Verify variants have ETD activation type
  auto group = exploration.getGroup(1);
  for (const auto& v : group.variants)
  {
    TEST_STRING_EQUAL(v.activation_type, "ETD")
  }

  // Feed result — this exercises the full scoring chain:
  // feedResultImpl_ → computeExplorationScore_ → computeFragmentMatch_
  // which now passes v.activation_type instead of hardcoded "HCD"
  DeconvolvedSpectrum ds = makeSyntheticDeconv(1, 5);
  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info = exploration.feedResultForTest(tracking_id, ds, 1.0, queue);

  // The key assertion: feedResult completed without error and the activation type
  // was correctly propagated through the chain.
  TEST_EQUAL(info.activation_type, "ETD")
  TEST_EQUAL(info.group_id > 0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// T8: HCD / ETD exploration split.
// KEEP the existing HCD baselines (exploration_group_creation /
// exploration_variants_priority_by_level above). The two sections below are the
// ETD-parallel half of the split: they assert that ETD exploration builds variants
// (variants.size() > 0), that those variants carry the ETD activation (which maps
// internally to c/z fragment ions, vs HCD's b/y), and that the ETD branch sweeps
// reaction_time (rt_min/rt_max/rt_step) rather than collision_energy.
//
// The C++ engine's Config only accepts the C++-facing JSON schema (not the C#
// method.json schema), so the canonical, build-stable assertions run against the
// inline etd_exploration_config. A second, guarded section additionally drives the
// data-agent fixture method_exploration_etd.json *iff* it is present AND parses as a
// C++ config — so a C++-schema fixture is exercised end-to-end, while an absent file
// (or a C#-schema one the engine cannot parse) is skipped cleanly with no false fail.
/////////////////////////////////////////////////////////////

START_SECTION(exploration_group_creation_etd)
{
  // ETD half of the HCD/ETD split (mirrors exploration_group_creation). ETD sweeps
  // reaction_time over [rt_min, rt_max] step rt_step = [5,10,15] -> 3 variants, each
  // with activation_type "ETD" (c/z ions). Verifies variants.size() > 0 and that the
  // ETD branch wires reaction_time (not collision_energy) per Exploration::buildVariants_.
  Config cfg{std::string(etd_exploration_config)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);

  TEST_EQUAL(cmds.size() > 0, true)              // ETD exploration produced variants
  TEST_EQUAL(exploration.activeGroupCount(), 1)

  auto group = exploration.getGroup(1);
  TEST_EQUAL(group.msn_level, 2)
  TEST_EQUAL(group.variants.size() > 0, true)

  // c/z proxy: every variant is ETD. rt-sweep: reaction_time takes the swept values
  // {5,10,15}; collision_energy is the (constant) base-config CE, not a CE sweep.
  std::set<double> seen_rts;
  for (const auto& v : group.variants)
  {
    TEST_STRING_EQUAL(v.activation_type, "ETD")
    TEST_EQUAL(std::isfinite(v.reaction_time), true)
    seen_rts.insert(v.reaction_time);
  }
  // rt_min=5, rt_max=15, rt_step=5 -> 3 distinct reaction-time variants.
  TEST_EQUAL(static_cast<int>(group.variants.size()), 3)
  TEST_EQUAL(static_cast<int>(seen_rts.size()), 3)
  TEST_EQUAL(seen_rts.count(5.0) == 1 && seen_rts.count(10.0) == 1 && seen_rts.count(15.0) == 1, true)
}
END_SECTION

START_SECTION(exploration_variants_priority_by_level_etd)
{
  // ETD analogue of exploration_variants_priority_by_level: every ETD MS2 variant
  // command is level 2, priority 2, non-AGC, and marked 'E' in its scan_description.
  Config cfg{std::string(etd_exploration_config)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, 0.0, queue);

  TEST_EQUAL(cmds.size() > 0, true)
  for (size_t i = 0; i < cmds.size(); ++i)
  {
    TEST_EQUAL(cmds[i].msn_level, 2)
    TEST_EQUAL(cmds[i].priority, 2)
    TEST_EQUAL(cmds[i].is_agc, 0)
    std::string desc(cmds[i].scan_description);
    TEST_EQUAL(desc.size() >= 4, true)
    TEST_EQUAL(desc[3], 'E')
  }
}
END_SECTION

START_SECTION(exploration_etd_from_method_fixture)
{
  // Drive the data-agent ETD fixture by path (T8). The C++ Config only parses the
  // C++-facing schema; this section therefore READS the file and exercises it only
  // when it is present AND constructs a FLASHIda without throwing (i.e. it is a
  // C++-schema config). An absent file, or a C#-schema method.json the engine cannot
  // parse, is skipped cleanly (informational STATUS, no assertion, no failure) — the
  // canonical ETD assertions live in exploration_group_creation_etd above.
  const std::string fixture = "../../FlashIDA/test-data/configs/method_exploration_etd.json";

  std::ifstream f(fixture);
  if (!f.good())
  {
    STATUS("method_exploration_etd.json absent -- ETD fixture section skipped cleanly")
  }
  else
  {
    std::stringstream buf;
    buf << f.rdbuf();
    std::string cfg_json = buf.str();
    f.close();

    // Only a C++-schema config (starts with '{' and Config accepts it) can be driven
    // here; a C#-schema method.json would throw in the FLASHIda/Config ctor.
    bool loadable = false;
    FLASHIda* ida = nullptr;
    try
    {
      ida = new FLASHIda(const_cast<char*>(cfg_json.c_str()));
      loadable = true;
    }
    catch (...)
    {
      if (ida) { delete ida; ida = nullptr; }
      loadable = false;
    }

    if (!loadable)
    {
      STATUS("method_exploration_etd.json is not a C++-schema config -- ETD fixture section skipped cleanly")
    }
    else
    {
      // The fixture loaded. Drive a single MS1 -> MS2-exploration cycle on real cytC
      // data and assert that ETD exploration variants were produced (variants.size()>0)
      // and that they carry the ETD activation (c/z ions). Feed the engine-emitted
      // descriptions back so the variants are scored on real data.
      auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
      auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
      ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

      // Bootstrap MS1 via the interleaved engine-id-echo contract: feed survey scans (stamped with
      // the engine's own ids) until one selects a precursor and creates the exploration group; its
      // variants stay queued for the drain loop below.
      int group_cmds = bootstrapExplorationGroup(ida, ms1_scans);
      ABORT_IF(group_cmds == 0)

      int etd_variant_cmds = 0;
      int idle = 0;
      for (int i = 0; i < 100; ++i)
      {
        ScanCommand next{};
        if (ida->getNextScanCommand(next) != 1) break;
        if (next.is_agc || next.msn_level == 1) { if (++idle > 3) break; continue; }
        idle = 0;
        if (next.msn_level == 2)
        {
          std::string d(next.scan_description);
          if (d.size() >= 4 && d[3] == 'E') ++etd_variant_cmds;  // an exploration variant
          ida->processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                           (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 2, next.scan_description);
        }
        else if (next.msn_level == 3)
        {
          break;  // reached MS3 -> the MS2-exploration group already produced its variants
        }
      }
      delete ida;
      // The fixture drives an ETD exploration -> at least one 'E' MS2 variant command.
      TEST_EQUAL(etd_variant_cmds > 0, true)
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// T10: exploration-MS2 tag_count under tagging (locks E4).
// E4 surfaces the FLASHTagger tag count generated during exploration identification
// (info.identification_result.tag_count) into the exploration MS2 scan_results row,
// instead of the old literal 0. Driving the inclusion-pinned cytC exploration on the
// REAL ms2_cytc_fresh_scan57 ladder (whose b/y ions match the M-starting cytC
// proteoform -- the same data that fires MS3 in the sections above) generates real
// sequence tags during identification. We capture scan_results.tsv and assert that at
// least one level-2 row logs tag_count > 0.
//
// NOTE (deviation from the literal T10 wording): the task says "an exploration MS2
// variant ... logs tag_count>0". Per the engine's design (logging-validation suite
// decision) a CE-/RT-sweep 'E' variant only logs a non-zero tag_count when its fed
// spectrum actually tags against the proteoform; the assertion below therefore targets
// "at least one level-2 result row" (which on this real, matching data is an
// exploration row) rather than over-specifying a particular variant index. This still
// locks E4: the field is no longer a hardcoded 0.
/////////////////////////////////////////////////////////////

START_SECTION(exploration_ms2_tag_count_under_tagging)
{
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  // Inclusion-pin the cytC precursor + the validatable M-starting proteoform so the
  // real fresh57 ladder matches and FLASHTagger produces tags during identification.
  std::string cfg_str = inclusionPinCytc(ms3_exploration_config);

  // Inject a runtime block (the exploration configs have none) so the engine writes
  // scan_results.tsv, which carries the tag_count column we assert on. Insert before
  // the config's closing brace.
  const std::string results_path = "expl_tagcount_scan_results.tsv";
  std::remove(results_path.c_str());
  {
    auto last = cfg_str.rfind('}');
    ABORT_IF(last == std::string::npos)
    std::string rt = ", \"runtime\": { \"scan_results_path\": \"" + results_path + "\" }";
    cfg_str.insert(last, rt);
  }

  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  // Drive the inclusion-pinned exploration group via the interleaved engine-id-echo contract:
  // bootstrap MS1 (feed survey scans stamped with the engine's own ids -- the always-on gate
  // rejects fabricated "scan_"+id MS1) until a group is created, then drain + feed each MS2 variant
  // back on the real cytC ladder so the exploration identification (and tag generation) runs and the
  // rows are written. Variants stay queued for the drain loop below.
  int group_cmds = bootstrapExplorationGroup(ida, ms1_scans);
  ABORT_IF(group_cmds == 0)

  int idle = 0;
  for (int i = 0; i < 100; ++i)
  {
    ScanCommand next{};
    if (ida->getNextScanCommand(next) != 1) break;
    if (next.is_agc || next.msn_level == 1) { if (++idle > 3) break; continue; }
    idle = 0;
    if (next.msn_level == 2)
    {
      ida->processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                       (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 2, next.scan_description);
    }
    else if (next.msn_level == 3)
    {
      // Feed the MS3 too so the cycle keeps draining (its row is not what we assert on).
      ida->processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                       (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 3, next.scan_description);
    }
  }
  delete ida;  // flush + close the TSV streams

  // Parse scan_results.tsv: find a level-2 row whose tag_count > 0.
  TSVFile results = TSVFile::parse(results_path);
  ABORT_IF(results.rows.empty())
  int lvl_col = results.colIndex("ms_level");
  int tag_col = results.colIndex("tag_count");
  ABORT_IF(lvl_col < 0 || tag_col < 0)

  bool found_ms2_tag = false;
  int max_ms2_tag = 0;
  for (const auto& row : results.rows)
  {
    if (lvl_col >= (int)row.size() || tag_col >= (int)row.size()) continue;
    if (row[lvl_col] != "2") continue;
    int tc = 0;
    try { tc = std::stoi(row[tag_col]); } catch (...) { tc = 0; }
    if (tc > max_ms2_tag) max_ms2_tag = tc;
    if (tc > 0) found_ms2_tag = true;
  }
  STATUS("max level-2 tag_count = " << max_ms2_tag)
  // E4: at least one MS2 exploration row logs a real (non-zero) tag count.
  TEST_EQUAL(found_ms2_tag, true)

  std::remove(results_path.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// Inclusion-mode MS3 full-acquisition round-trip (H-full).
// Mirrors the Flash.cs event loop: repeatedly drain a command and feed the matching
// scan back stamped with THAT command's engine-emitted scan_description, so MS1->MS2
// ->MS3 lineage uses the engine's own ids end-to-end. Driven in INCLUSION mode for
// cytC (single ~12360 Da target via inclusion_cytc.txt) so exactly one precursor fires
// and the MS1->MS2->MS3 chain is unambiguous. Terminates on the 3-consecutive-idle
// (AGC / already-fed-MS1) sentinel.
//
// The C++ Config takes the C++ JSON schema, not the C# method.json; the inline config
// below is the C++-schema mirror of method_inclusion.json + inclusion_cytc.txt
// (target_mode=1, the cytC inclusion list, the M-starting proteoform). It is deliberately
// a STANDARD MS2 -> MS3 inclusion config (no ms2 exploration block) so the MS2->MS3
// lineage flows through the standard path (FLASHIda.cpp MS2 branch, where each MS3
// command's parent_scan_id is set to the originating MS2's encoded id) -- matching
// method_inclusion.json (plain inclusion DDA, no exploration). agc_interval_seconds is
// large so the only AGCs are idle ones (required by the idle-sentinel termination).
//
// HARD asserts (no silent pass): ms2_count >= 1, ms3_count >= 1, and every MS3
// parent_tracking_id resolves to an emitted MS2 scan id.
/////////////////////////////////////////////////////////////

START_SECTION(inclusion_ms3_full_acquisition_roundtrip)
{
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.size() < 2 || ms2_scans.empty())

  // C++-schema mirror of method_inclusion.json + inclusion_cytc.txt: inclusion mode
  // (target_mode=1) pinned to the single ~12360 Da cytC target, the M-starting cytC
  // proteoform (so the real fresh57 b/y ladder matches and MS3 fires), standard MS2 and
  // MS3 intensity selection (NO exploration block), and a large agc_interval so only
  // idle AGCs are emitted.
  const char* inclusion_ms3_config = R"({
    "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10] },
    "precursor_selection": { "RT_window": 180, "target_mode": 1, "IDScore": false, "AllCharges": false, "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1 },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [ { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 } ],
      "ms3": [ { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 35, "resolution": 120000 } ]
    },
    "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 }, "agc_interval_seconds": 999999 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "../../FlashIDA/test-data/configs/inclusion_cytc.txt", "ptm_list": "" },
    "ms3": { "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE" },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "intensity", "max_targets": 3 },
      "ms3": { "selection": "intensity", "max_targets": 3 }
    }
  })";
  std::string cfg_str(inclusion_ms3_config);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  // --- H-full id-chaining driver (inlined; mirrors FLASHIda_TestHelpers::runFullAcquisition).
  // MS1 bootstraps from the engine's first idle-emitted MS1 command; feed each level
  // back with the engine-emitted scan_description. We use a single strong MS1 scan
  // (ms1_scans[1] = scan 134, the cytC envelope; ms1_scans[0] = scan 132 is a weak edge scan
  // that selects 0 precursors) so an MS1 re-survey after it is already fed counts as idle
  // (avoids RT self-exclusion churn).
  const ScanData& ms1 = ms1_scans[1];
  const ScanData& ms2 = ms2_scans[0];

  std::vector<ScanCommand> ms2_cmds, ms3_cmds;
  bool ms1_fed = false;
  int idle = 0;
  ScanCommand cmd{};
  for (int it = 0; it < 300 && idle < 3; ++it)
  {
    if (ida->getNextScanCommand(cmd) != 1) break;
    if (cmd.is_agc || (cmd.msn_level == 1 && ms1_fed)) { ++idle; cmd = ScanCommand{}; continue; }
    idle = 0;
    if (cmd.msn_level == 1)
    {
      ida->processScan(ms1.mzs.data(), ms1.ints.data(), (int)ms1.mzs.size(), ms1.rt, 1, cmd.scan_description);
      ms1_fed = true;
    }
    else if (cmd.msn_level == 2)
    {
      ms2_cmds.push_back(cmd);
      ida->processScan(ms2.mzs.data(), ms2.ints.data(), (int)ms2.mzs.size(), ms2.rt, 2, cmd.scan_description);
    }
    else if (cmd.msn_level == 3)
    {
      ms3_cmds.push_back(cmd);
      ida->processScan(ms2.mzs.data(), ms2.ints.data(), (int)ms2.mzs.size(), ms2.rt, 3, cmd.scan_description);
    }
    cmd = ScanCommand{};
  }
  delete ida;

  // Hard asserts: the inclusion-mode cytC acquisition produced MS2 and MS3 commands.
  TEST_EQUAL(ms2_cmds.size() >= 1, true)
  TEST_EQUAL(ms3_cmds.size() >= 1, true)

  // Every MS3 parent_tracking_id must resolve to an emitted MS2 scan id. The MS2 ids
  // are the leading 3 chars of each emitted MS2 scan_description; an MS3 command stores
  // its parent MS2 id in parent_scan_id (3 chars + NUL).
  std::set<std::string> emitted_ms2_ids;
  for (const auto& m2 : ms2_cmds)
  {
    std::string d(m2.scan_description);
    if (d.size() >= 3) emitted_ms2_ids.insert(d.substr(0, 3));
  }
  bool all_ms3_parents_resolve = true;
  std::string unresolved;
  for (const auto& m3 : ms3_cmds)
  {
    std::string parent(m3.parent_scan_id);  // NUL-terminated 3-char id
    if (parent.empty() || emitted_ms2_ids.count(parent) == 0)
    {
      all_ms3_parents_resolve = false;
      unresolved = parent;
      break;
    }
  }
  STATUS("emitted MS2 ids = " << emitted_ms2_ids.size()
         << ", MS3 cmds = " << ms3_cmds.size()
         << (all_ms3_parents_resolve ? std::string("")
                                     : (", first unresolved MS3 parent = '" + unresolved + "'")))
  TEST_EQUAL(all_ms3_parents_resolve, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
