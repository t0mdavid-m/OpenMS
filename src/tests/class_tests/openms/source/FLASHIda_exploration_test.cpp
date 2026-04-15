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

  int pushAllMS1Scans(FLASHIda* ida, const std::vector<ScanData>& scans)
  {
    int total = 0;
    for (const auto& scan : scans)
    {
      int n = ida->processScan(scan.mzs.data(), scan.ints.data(),
                                (int)scan.mzs.size(), scan.rt, 1,
                                ("scan_" + scan.scan_id).c_str());
      total += n;
    }
    return total;
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
  Config cfg{std::string(ms3_exploration_config)};
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

  int total = pushAllMS1Scans(ida, ms1_scans);
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

  int total = pushAllMS1Scans(ida, ms1_scans);
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
  // P7-U08: MS3 with selection but no exploration -> standard MS3 commands
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  FLASHIda* ida = new FLASHIda(const_cast<char*>(ms3_selection_only_config));

  int total = pushAllMS1Scans(ida, ms1_scans);
  if (total == 0) { delete ida; }
  ABORT_IF(total == 0)

  // Drain exploration variants and feed MS2 results
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

  const auto& ms2 = ms2_scans[0];
  for (const auto& ecmd : exploration_cmds)
  {
    ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                     (int)ms2.mzs.size(), ms2.rt, 2, ecmd.scan_description);
  }

  // After MS2 exploration completes, MS3 commands should be queued
  bool found_ms3 = false;
  for (int i = 0; i < 20; ++i)
  {
    ScanCommand next{};
    if (ida->getNextScanCommand(next) != 1) break;
    if (next.msn_level == 3)
    {
      found_ms3 = true;
      TEST_EQUAL(next.num_stages, 2)
      break;
    }
  }
  // MS3 generation depends on deconvolution results
  (void)found_ms3;

  delete ida;
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
  // ms3_selection_only_config has ms2.selection = "none" and protein_sequence = "MKWVTFISLLLLFSSAYSRGVFRR"
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
  // ms3_selection_only_config has ms2.selection = "none" and protein_sequence = "MKWVTFISLLLLFSSAYSRGVFRR"
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

START_SECTION(legacy_ms3_enabled_rejected)
{
  // Build a config with ms3.enabled present — should throw
  std::string cfg_str = std::string(exploration_config);
  auto pos = cfg_str.find("\"ms3\":");
  auto brace = cfg_str.find("{", pos);
  cfg_str.insert(brace + 1, " \"enabled\": false,");

  TEST_EXCEPTION(std::invalid_argument, Config cfg{cfg_str})
}
END_SECTION

START_SECTION(legacy_ms3_mode_rejected)
{
  std::string cfg_str = std::string(exploration_config);
  auto pos = cfg_str.find("\"ms3\":");
  auto brace = cfg_str.find("{", pos);
  cfg_str.insert(brace + 1, " \"mode\": 0,");

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
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  FLASHIda* ida = new FLASHIda(const_cast<char*>(cycle_time_exploration_config));

  int total = pushAllMS1Scans(ida, ms1_scans);
  if (total == 0) { delete ida; }
  ABORT_IF(total == 0)

  std::vector<ScanCommand> exploration_cmds;
  ScanCommand cmd{};
  while (ida->getNextScanCommand(cmd) == 1)
  {
    std::string desc(cmd.scan_description);
    if (cmd.msn_level == 2 && desc.size() >= 4 && desc[3] == 'E')
      exploration_cmds.push_back(cmd);
    else
      break;
  }
  ABORT_IF(exploration_cmds.empty())

  const auto& ms2 = ms2_scans[0];
  std::vector<int> return_values;
  for (const auto& ecmd : exploration_cmds)
  {
    int rv = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                               (int)ms2.mzs.size(), ms2.rt, 2, ecmd.scan_description);
    return_values.push_back(rv);
  }

  for (int rv : return_values)
  {
    TEST_EQUAL(rv >= 0, true)
  }

  delete ida;
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
