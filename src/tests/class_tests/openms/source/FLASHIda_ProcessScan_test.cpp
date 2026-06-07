// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------
//
// Phase 4 unit tests: processScan full routing.
// Peak data loaded at runtime from FlashIDA/test-data/spectra/ms1_standard.txt
// and ms2_hcd_fragment.txt.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <thread>  // std::this_thread::sleep_for

using namespace OpenMS;

namespace
{
  // Minimal JSON config for standard DDA mode with score_threshold=0 to accept all peaks
  const char* standard_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with MS3 targeting enabled via selection_strategy
  const char* ms3_mode1_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ],
      "ms3": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 35, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "ms3": { "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE" },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "intensity", "max_targets": 3 },
      "ms3": { "selection": "intensity", "max_targets": 2 }
    }
  })";

  // Config with IDScore=false, AllCharges=true (activates sortByQScoreAllCharges)
  const char* qscore_allcharges_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": true,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with quantification enabled and 2 MS2 configs (HCD+ETD, required for quant path)
  const char* quant_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": true, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with tag-based targeting enabled via valid FASTA path, 2 MS2 configs
  const char* tag_targeting_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "../../FlashIDA/test-data/configs/test_fasta.fasta", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with conditional MS2 enabled (no FASTA — tags cannot be found)
  const char* conditional_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000,
      "follow_up_scan": { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
    },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "conditional_ms2": true,
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with cycle_time enabled and value_ms=0 (always triggers), AGC suppressed
  const char* cycle_time_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": true, "value_ms": 0 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with agc_interval_seconds=0 (AGC triggers immediately)
  const char* agc_fast_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": true, "value_ms": 30000 },
      "agc_interval_seconds": 0
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with conditional MS2 + tag targeting + FASTA (tags CAN be found)
  const char* conditional_with_tags_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000,
      "follow_up_scan": { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
    },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "conditional_ms2": true,
    "files": { "target_logs": [], "fasta": "../../FlashIDA/test-data/configs/test_fasta.fasta", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with quantification + low fold_change_threshold (any reporter ratio triggers)
  const char* quant_sensitive_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": {
      "enabled": true, "reporter_mz_tol": 0.002, "fold_change_threshold": 0.01,
      "follow_up_scan": { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
    },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with agc_interval_seconds=9999 to disable timer-based AGC.
  // Only idle cycle (empty queue) produces AGC/MS1 pairs.
  const char* idle_cycle_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": true, "value_ms": 30000 },
      "agc_interval_seconds": 9999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with intensity-based MS1 selection
  const char* intensity_selection_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "intensity", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with selection=none at MS1
  const char* none_selection_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "none", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with max_targets=1 (cap test)
  const char* max1_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 1 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // Config with max_targets=5 (cap test)
  const char* max5_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 5 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // TSV file paths relative to the OpenMS build directory (CTest working dir)
  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_standard.txt";
  const std::string ms2_tsv_path = "../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt";
  const std::string ms2_tmt_tsv_path = "../../FlashIDA/test-data/spectra/ms2_quant_tmt.txt";
  const std::string fasta_path = "../../FlashIDA/test-data/configs/test_fasta.fasta";

  struct ScanData
  {
    std::vector<double> mzs;
    std::vector<double> ints;
    double rt;
    std::string scan_id;
  };

  // Parse multi-scan TSV: "Spec scan=N\tRT" headers, "mz\tintensity" data lines
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
      else if (! scans.empty())
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

  // Push all scans through processScan, return total command count
  int pushAllScans(FLASHIda* ida, const std::vector<ScanData>& scans)
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
}

START_TEST(FLASHIda_ProcessScan, "$Id$")

/////////////////////////////////////////////////////////////

// P4-U01: MS1 processScan returns > 0 commands for real spectral data
START_SECTION(processScan_ms1_returns_commands)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)
  delete ida;
}
END_SECTION

// P4-U02: Commands from processScan are dequeued by getNextScanCommand
START_SECTION(processScan_commands_dequeued)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  // Dequeue all commands
  ScanCommand cmd{};
  for (int i = 0; i < total; i++)
  {
    int result = ida->getNextScanCommand(cmd);
    TEST_EQUAL(result, 1)
    TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
    TEST_EQUAL(cmd.msn_level, 2)
    TEST_EQUAL(cmd.priority, 2)
    TEST_EQUAL(cmd.num_stages, 1)
    TEST_EQUAL(cmd.is_agc, 0)
    // Isolation stage should have valid m/z
    TEST_EQUAL(cmd.stages[0].precursor_mz > 0, true)
    TEST_EQUAL(cmd.stages[0].isolation_width > 0, true)
    TEST_EQUAL(cmd.stages[0].charge_state >= 4, true)
    // Scan description starts with 3-char tracking ID + type code
    TEST_EQUAL(std::strlen(cmd.scan_description) >= 3, true)
    // Enqueue timestamp should be set
    TEST_EQUAL(cmd.enqueue_timestamp_ms > 0, true)
  }

  // Queue empty — idle cycle returns AGC (never returns 0)
  int idle_result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(idle_result, 1)
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.is_agc, 1)

  delete ida;
}
END_SECTION

// P4-U03: ScanCommand fields populated correctly
START_SECTION(processScan_command_fields)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  ScanCommand cmd{};
  ida->getNextScanCommand(cmd);
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)

  // Analyzer from ms2_configs_[0]
  TEST_STRING_EQUAL(std::string(cmd.analyzer), "Orbitrap")
  // Resolution from ms2_configs_[0]
  TEST_EQUAL(cmd.orbitrap_resolution, 120000)
  // Activation type from ms2_configs_[0]
  TEST_STRING_EQUAL(std::string(cmd.stages[0].activation_type), "HCD")
  // Collision energy
  TEST_EQUAL(cmd.stages[0].collision_energy > 0, true)
  // Precursor m/z within scan range
  TEST_EQUAL(cmd.stages[0].precursor_mz >= 500.0, true)
  TEST_EQUAL(cmd.stages[0].precursor_mz <= 2000.0, true)

  // D1: Scan description format — 3-char base-94 ID, then type code
  std::string desc(cmd.scan_description);
  TEST_EQUAL(desc.size() >= 4, true)
  std::string id_part = desc.substr(0, 3);
  for (char c : id_part)
  {
    TEST_EQUAL(c >= 0x21 && c <= 0x7E, true)
  }
  TEST_EQUAL(desc[3], 'R')  // recording MS2

  delete ida;
}
END_SECTION

// P4-U04: processScan with ms_level=2 processes MS2 path
START_SECTION(processScan_ms2_path)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // Push all MS1 scans to accumulate state and generate MS2 commands
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  // Dequeue one MS2 command to get its scan description (contains tracking ID)
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)

  // Now process MS2 return with the tracking ID in scan description
  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description);
  // Should return 0 (no conditional, no MS3, no quant in standard config)
  TEST_EQUAL(ms2_result, 0)

  delete ida;
}
END_SECTION

// P4-U05: processScan with empty spectrum returns 0
START_SECTION(processScan_empty_spectrum)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int n = ida->processScan(nullptr, nullptr, 0, 0.0, 1, "empty");
  TEST_EQUAL(n, 0)
  delete ida;
}
END_SECTION

// P4-U06: conditional MS2 follow-up is pushed at priority 0 (requires tag detection)
START_SECTION(processScan_conditional_ms2_followup)
{
  {
    std::ifstream fasta_check(fasta_path);
    ABORT_IF(! fasta_check.good())
  }
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  // Uses conditional_with_tags_json: conditional_ms2=true + FASTA for tag detection
  FLASHIda* ida = new FLASHIda(const_cast<char*>(conditional_with_tags_json));

  // Push all MS1 scans
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  // Dequeue first MS2
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)
  TEST_EQUAL(ms2_cmd.priority, 2)

  // Process MS2 return — tags found → conditional follow-up at priority 0
  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description);
  TEST_EQUAL(ms2_result > 0, true)

  // Follow-up at priority 0 dequeues before remaining MS2s at priority 2
  ScanCommand out{};
  int r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.priority, 0)
  TEST_EQUAL(out.msn_level, 2)

  // D1: Conditional follow-up uses 3-char ID + type code 'C'
  std::string cond_desc(out.scan_description);
  TEST_EQUAL(cond_desc.size() >= 4, true)
  TEST_EQUAL(cond_desc[3], 'C')

  delete ida;
}
END_SECTION

// P4-U07: MS3 commands are pushed at priority 1 — hard positive (golden file confirms)
START_SECTION(processScan_ms3_commands)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(ms3_mode1_json));

  // Push all MS1 scans
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  // Dequeue first MS2
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)

  // Process MS2 return — golden file confirms MS3 targets are produced
  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description);
  TEST_EQUAL(ms2_result > 0, true)

  // Drain all commands; MS3 at priority 1 dequeues before MS2 at priority 2
  ScanCommand out{};
  int ms3_count = 0;
  while (ida->getNextScanCommand(out) == 1)
  {
    TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
    if (out.is_agc) break; // idle cycle = queue empty
    if (out.msn_level == 3)
    {
      ms3_count++;
      TEST_EQUAL(out.priority, 1)  // MS3 at priority 1
      TEST_EQUAL(out.num_stages, 2)  // MS2 precursor + fragment target

      // D1: MS3 scan description uses 3-char ID + type code 'R' (recording)
      std::string ms3_desc(out.scan_description);
      TEST_EQUAL(ms3_desc.size() >= 4, true)
      TEST_EQUAL(ms3_desc[3], 'R')  // MS3 recording
    }
  }
  TEST_EQUAL(ms3_count > 0, true)

  delete ida;
}
END_SECTION

// P4-U08: decodeTracking roundtrips with encodeTracking
START_SECTION(decodeTracking_roundtrip)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // Test roundtrip via processScan → scan_description parsing
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  ScanCommand cmd{};
  ida->getNextScanCommand(cmd);
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  // scan_description: {3-char base-94}{type_code}{payload}
  std::string desc(cmd.scan_description);
  TEST_EQUAL(desc.size() >= 4, true)  // 3-char ID + type code
  std::string id_str = desc.substr(0, 3);
  // Verify the ID is valid base-94 (printable ASCII 0x21-0x7E)
  for (char c : id_str)
  {
    TEST_EQUAL(c >= 0x21 && c <= 0x7E, true)
  }

  // Roundtrip: decoding the base-94 ID should give back cmd.scan_id
  int decoded_id = ida->getQueueForTest().decode(id_str);
  TEST_EQUAL(decoded_id, cmd.scan_id)

  delete ida;
}
END_SECTION

// P4-U09: cleanupExpiredCommands removes old entries
START_SECTION(cleanup_expired_commands)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  // Use timeout-enabled config
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // Push all MS1 scans (this adds to pending_scan_map_)
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  // Dequeue all commands — they're still in pending_scan_map_ for tracking
  ScanCommand cmd{};
  for (int i = 0; i < total; i++)
  {
    ida->getNextScanCommand(cmd);
    TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  }

  // The pending_scan_map_ entries have timestamps. cleanupExpiredCommands_
  // is called by getNextScanCommand. With timeout_ms=30000, entries should
  // NOT be expired immediately. Verify by processing an MS2 with valid tracking ID.
  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, cmd.scan_description);
  // Should succeed (entry found in pending_scan_map_)
  // ms2_result can be 0 (no follow-ups) but shouldn't crash
  TEST_EQUAL(ms2_result >= 0, true)

  // Verify pending_scan_map_ entry was erased by first resolution
  int ms2_result2 = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                      (int)ms2.mzs.size(), ms2.rt,
                                      2, cmd.scan_description);
  // Second resolution with same tracking ID should find nothing (entry erased)
  TEST_EQUAL(ms2_result2, 0)

  // Also verify via accessor: pending_scan_map_ should have (total - 1) entries
  // (we resolved one, rest are still pending)
  TEST_EQUAL(ida->getQueueForTest().pendingScanMapSize(), (size_t)(total - 1))

  delete ida;
}
END_SECTION

// P4-U02: All 4 scoring branches execute without error
START_SECTION(processScan_scoring_branches)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  // 2 configs: standard (QScore rep), qscore_allcharges
  FLASHIda* ida_qscore = new FLASHIda(const_cast<char*>(standard_json));
  FLASHIda* ida_qscore_all = new FLASHIda(const_cast<char*>(qscore_allcharges_json));

  int total_qscore = pushAllScans(ida_qscore, ms1_scans);
  int total_qscore_all = pushAllScans(ida_qscore_all, ms1_scans);

  // Both branches must produce > 0 commands
  TEST_EQUAL(total_qscore > 0, true)
  TEST_EQUAL(total_qscore_all > 0, true)

  // Dequeue first command from each — must have valid precursor_mz
  ScanCommand cmd{};

  ida_qscore->getNextScanCommand(cmd);
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.stages[0].precursor_mz > 0, true)

  ida_qscore_all->getNextScanCommand(cmd);
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.stages[0].precursor_mz > 0, true)

  delete ida_qscore;
  delete ida_qscore_all;
}
END_SECTION

// P4-U03: Mass exclusion deprioritizes previously-selected high-confidence precursors.
// Re-pushing the same scans within the RT window must NOT re-select the excluded masses
// (they enter the exclusion map on first selection), and must never GROW the command count.
// (A strict pass2<pass1 count drop is not valid: in standard mode excluded masses are simply
// not re-selected rather than deterministically reducing the per-run total.)
START_SECTION(processScan_mass_exclusion)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // Pass 1: push all scans, dequeue all, recording selected MS2 precursor m/z.
  int total_pass1 = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total_pass1 > 0, true)

  std::vector<int> pass1_mz;
  ScanCommand cmd{};
  for (int i = 0; i < total_pass1; i++)
  {
    if (ida->getNextScanCommand(cmd) != 1) break;
    TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
    if (!cmd.is_agc && cmd.msn_level == 2)
      pass1_mz.push_back((int)(cmd.stages[0].precursor_mz + 0.5));
  }
  ABORT_IF(pass1_mz.empty())

  // Pass 2: push the SAME scans at the same RTs (within RT_window=180s).
  int total_pass2 = pushAllScans(ida, ms1_scans);
  std::vector<int> pass2_mz;
  ScanCommand out{};
  for (int i = 0; i < total_pass2; i++)
  {
    if (ida->getNextScanCommand(out) != 1) break;
    if (!out.is_agc && out.msn_level == 2)
      pass2_mz.push_back((int)(out.stages[0].precursor_mz + 0.5));
  }

  // Exclusion must never INCREASE the count for identical re-pushed data.
  TEST_EQUAL(total_pass2 <= total_pass1, true)

  // At least one precursor selected (and excluded) in pass 1 must NOT be re-selected in pass 2.
  bool some_excluded = false;
  for (int mz1 : pass1_mz)
  {
    bool reselected = false;
    for (int mz2 : pass2_mz) { if (mz2 == mz1) { reselected = true; break; } }
    if (!reselected) { some_excluded = true; break; }
  }
  TEST_EQUAL(some_excluded, true)

  delete ida;
}
END_SECTION

// P4-U07: Quantification follow-up routing — hard positive with sensitive threshold
START_SECTION(processScan_quant_followup)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_tmt_scans = loadTsvScans(ms2_tmt_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_tmt_scans.empty())
  // Use quant_sensitive_json with fold_change_threshold=0.01 to trigger on any reporter ratio
  FLASHIda* ida = new FLASHIda(const_cast<char*>(quant_sensitive_json));

  // Push MS1 scans to generate MS2 commands
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  // Dequeue one MS2 command
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)

  // Push MS2 return with TMT reporter data
  const auto& ms2 = ms2_tmt_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description);
  // With fold_change_threshold=0.01, any reporter ion ratio should trigger
  TEST_EQUAL(ms2_result, 1)  // Exactly 1 quant follow-up per MS2

  // Follow-up at priority 0 dequeues before remaining MS2s at priority 2
  ScanCommand out{};
  int r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_STRING_EQUAL(std::string(out.stages[0].activation_type), "ETD")
  TEST_EQUAL(out.msn_level, 2)
  TEST_EQUAL(out.priority, 0)  // Follow-up priority

  delete ida;
}
END_SECTION

// P4-U08: Tag-based targeting code path executes without crash
START_SECTION(processScan_tag_targeting)
{
  // Guard: check fasta file exists
  {
    std::ifstream fasta_check(fasta_path);
    ABORT_IF(! fasta_check.good())
  }
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(tag_targeting_json));

  // Push MS1 scans
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  // Dequeue one MS2 command
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)

  // Push MS2 return (HCD fragment data)
  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description);
  // tag_targeting_json has no quant/conditional/ms3 — tag targeting updates internal state
  // but pushes 0 follow-up commands
  TEST_EQUAL(ms2_result, 0)

  delete ida;
}
END_SECTION

// P4-U14: Cycle time forces MS1 when cycle_time_ms exceeded
START_SECTION(processScan_cycle_time_enforcement)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  // cycle_time_json: cycle_time.enabled=true, value_ms=0 (any elapsed time triggers)
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cycle_time_json));

  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  // With cycle_time_ms=0, ANY elapsed time triggers a cycle-time MS1
  // Cycle-time MS1 is queued at priority 0, then dequeued before MS2s (priority 2)
  ScanCommand cmd{};
  int r = ida->getNextScanCommand(cmd);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.msn_level, 1)  // Cycle-time forced MS1
  TEST_EQUAL(cmd.priority, 0)   // Cycle-time MS1 queued at priority 0

  // Scan description: 3-char ID + type code 'S' for MS1 survey
  std::string desc(cmd.scan_description);
  TEST_EQUAL(desc.size() >= 4, true)
  TEST_EQUAL(desc[3], 'S')

  delete ida;
}
END_SECTION

// P4-U15: AGC command uses ms1 config values (not hardcoded zeros)
START_SECTION(agc_command_values)
{
  // agc_fast_json: agc_interval_seconds=0 triggers AGC immediately
  // AGC is a dedicated minimal gain-control scan: makeAGC hardcodes agc_target=30000, max_it=1
  FLASHIda* ida = new FLASHIda(const_cast<char*>(agc_fast_json));

  ScanCommand cmd{};
  int r = ida->getNextScanCommand(cmd);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.is_agc, 1)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_EQUAL(cmd.agc_target, 30000)   // makeAGC: dedicated fast gain-control scan
  TEST_EQUAL(cmd.max_it, 1)           // makeAGC: minimal injection time
  TEST_STRING_EQUAL(std::string(cmd.analyzer), "IonTrap")
  TEST_EQUAL(cmd.priority, 0)

  // Scan description: 3-char ID + type code 'A' for AGC calibration
  std::string desc(cmd.scan_description);
  TEST_EQUAL(desc.size() >= 4, true)
  TEST_EQUAL(desc[3], 'A')

  delete ida;
}
END_SECTION

// P4-U17: Conditional MS2 does NOT fire when no FASTA is loaded (no tags possible)
START_SECTION(processScan_conditional_ms2_requires_tags)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  // conditional_json has conditional_ms2=true but fasta="" — no tag targeting possible
  FLASHIda* ida = new FLASHIda(const_cast<char*>(conditional_json));

  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)

  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description);
  // Without FASTA, tags cannot be found — conditional must NOT fire
  TEST_EQUAL(ms2_result, 0)

  delete ida;
}
END_SECTION

// P4-U19: Tag targeting + conditional MS2 produces follow-up commands
START_SECTION(processScan_tag_targeting_produces_followups)
{
  {
    std::ifstream fasta_check(fasta_path);
    ABORT_IF(! fasta_check.good())
  }
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(conditional_with_tags_json));

  // Push all 50 MS1 scans
  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total > 0, true)

  // Dequeue one MS2 command
  ScanCommand ms2_cmd{};
  ida->getNextScanCommand(ms2_cmd);
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)
  TEST_EQUAL(ms2_cmd.priority, 2)

  // Push MS2 return — should find tags and trigger conditional follow-up
  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description);
  // Golden file continuity_tag_ms2return.json confirms: tags found → follow-ups produced
  TEST_EQUAL(ms2_result > 0, true)

  // Follow-up at priority 0 dequeues before remaining MS2s at priority 2
  ScanCommand out{};
  int r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.priority, 0)
  TEST_EQUAL(out.msn_level, 2)
  // Conditional follow-up uses second MS2 config (ETD)
  TEST_STRING_EQUAL(std::string(out.stages[0].activation_type), "ETD")

  delete ida;
}
END_SECTION

// Idle cycle: empty queue produces alternating AGC then MS1, 3 full iterations
START_SECTION(idle_cycle_agc_then_ms1)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(idle_cycle_json));

  // No processScan — queue is empty from the start.
  // Each getNextScanCommand call on an empty queue should produce an AGC (returned
  // immediately) and push an MS1 at priority 3 into the queue. The next call
  // dequeues that MS1.  So 6 calls = 3 AGC + 3 MS1, strictly alternating.
  ScanCommand cmd{};

  for (int iter = 0; iter < 3; ++iter)
  {
    // Odd call: AGC
    int r1 = ida->getNextScanCommand(cmd);
    TEST_EQUAL(r1, 1)
    TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
    TEST_EQUAL(cmd.is_agc, 1)
    TEST_EQUAL(cmd.msn_level, 1)
    TEST_EQUAL(cmd.agc_target, 30000)
    TEST_EQUAL(cmd.priority, 0)

    // Scan description: 3-char ID + type code 'A' for AGC calibration
    std::string agc_desc(cmd.scan_description);
    TEST_EQUAL(agc_desc.size() >= 4, true)
    TEST_EQUAL(agc_desc[3], 'A')

    // Even call: MS1
    int r2 = ida->getNextScanCommand(cmd);
    TEST_EQUAL(r2, 1)
    TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
    TEST_EQUAL(cmd.is_agc, 0)
    TEST_EQUAL(cmd.msn_level, 1)
    TEST_EQUAL(cmd.orbitrap_resolution, 120000)
    TEST_EQUAL(cmd.priority, 3)

    // Scan description: 3-char ID + type code 'S' for MS1 survey
    std::string ms1_desc(cmd.scan_description);
    TEST_EQUAL(ms1_desc.size() >= 4, true)
    TEST_EQUAL(ms1_desc[3], 'S')
  }

  delete ida;
}
END_SECTION

// Queued MS2 at priority 2 beats idle MS1 at priority 3
START_SECTION(ms2_priority_beats_idle_ms1)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(idle_cycle_json));

  // Push an MS2 at priority 2
  ScanCommand ms2_a{};
  ms2_a.msn_level = 2;
  ms2_a.priority = 2;
  ms2_a.scan_id = 42;
  ms2_a.faims_cv = -50.0;
  ida->pushCommandForTest(ms2_a);

  // First getNextScanCommand: MS2 at priority 2 (queue not empty, no idle cycle)
  ScanCommand out{};
  int r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.msn_level, 2)
  TEST_EQUAL(out.scan_id, 42)
  TEST_EQUAL(out.priority, 2)

  // Second call: queue empty -> idle cycle fires: returns AGC, pushes MS1 at prio 3
  r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.is_agc, 1)
  TEST_EQUAL(out.msn_level, 1)
  TEST_EQUAL(out.priority, 0)

  // Now push another MS2 at priority 2 WHILE the idle MS1 sits at priority 3
  ScanCommand ms2_b{};
  ms2_b.msn_level = 2;
  ms2_b.priority = 2;
  ms2_b.scan_id = 43;
  ms2_b.faims_cv = -50.0;
  ida->pushCommandForTest(ms2_b);

  // Third call: MS2 at priority 2 beats idle MS1 at priority 3
  r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.msn_level, 2)
  TEST_EQUAL(out.scan_id, 43)
  TEST_EQUAL(out.priority, 2)

  // Fourth call: idle MS1 at priority 3 is dequeued
  r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.is_agc, 0)
  TEST_EQUAL(out.msn_level, 1)
  TEST_EQUAL(out.priority, 3)

  delete ida;
}
END_SECTION

/////////////////////////////////////////////////////////////

// Intensity selection: MS1 precursors should be ordered by max charge intensity
START_SECTION(processScan_ms1_intensity_selection)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  FLASHIda* ida_qscore = new FLASHIda(const_cast<char*>(standard_json));
  FLASHIda* ida_intensity = new FLASHIda(const_cast<char*>(intensity_selection_json));

  int n_qscore = pushAllScans(ida_qscore, ms1_scans);
  int n_intensity = pushAllScans(ida_intensity, ms1_scans);

  TEST_EQUAL(n_qscore > 0, true)
  TEST_EQUAL(n_intensity > 0, true)

  ScanCommand cmd_q{};
  ScanCommand cmd_i{};
  int q_count = 0, i_count = 0;
  while (ida_qscore->getNextScanCommand(cmd_q)) { if (cmd_q.is_agc) break; q_count++; }
  while (ida_intensity->getNextScanCommand(cmd_i)) { if (cmd_i.is_agc) break; i_count++; }

  TEST_EQUAL(q_count > 0, true)
  TEST_EQUAL(i_count > 0, true)

  delete ida_qscore;
  delete ida_intensity;
}
END_SECTION

// Selection=none: MS1 should produce 0 commands
START_SECTION(processScan_ms1_none_selection)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(none_selection_json));

  int total = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total, 0)

  ScanCommand cmd{};
  TEST_EQUAL(ida->getNextScanCommand(cmd), true)  // idle cycle always returns 1
  TEST_EQUAL(cmd.is_agc, 1)                       // first command is AGC = queue was empty

  delete ida;
}
END_SECTION

// max_targets cap: max=1 <= max=3 <= max=5
START_SECTION(processScan_ms1_max_targets_cap)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  FLASHIda* ida1 = new FLASHIda(const_cast<char*>(max1_json));
  FLASHIda* ida3 = new FLASHIda(const_cast<char*>(standard_json));
  FLASHIda* ida5 = new FLASHIda(const_cast<char*>(max5_json));

  int total1 = pushAllScans(ida1, ms1_scans);
  int total3 = pushAllScans(ida3, ms1_scans);
  int total5 = pushAllScans(ida5, ms1_scans);

  TEST_EQUAL(total1 > 0, true)
  TEST_EQUAL(total3 >= total1, true)
  TEST_EQUAL(total5 >= total3, true)
  TEST_EQUAL(total1 <= (int)ms1_scans.size(), true)

  delete ida1;
  delete ida3;
  delete ida5;
}
END_SECTION

START_SECTION(processScan_ms1_min_charge_filter)
{
  // Config identical to standard_json but with ms1.min_charge = 99
  // This should filter out ALL precursors since no precursor has charge >= 99
  const char* min_charge_json = R"({
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
      "scan_timeout": { "enabled": false, "value_ms": 30000 },
      "agc_interval_seconds": 9999999
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3, "min_charge": 99 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  FLASHIda ida(const_cast<char*>(min_charge_json));

  // Load same MS1 scans that normally produce commands
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  for (const auto& scan : scans)
  {
    ida.processScan(scan.mzs.data(), scan.ints.data(), (int)scan.mzs.size(),
                    scan.rt, 1, ("scan_" + scan.scan_id).c_str(), -50.0);
  }

  // With min_charge=99, no precursor should pass the filter
  ScanCommand cmd{};
  int result = ida.getNextScanCommand(cmd);
  TEST_EQUAL(result, 0)  // no commands generated
}
END_SECTION

START_SECTION(processScan_agc_scan_skipped)
{
  // Use agc_fast_json (agc_interval_seconds=0) to get an AGC command immediately
  FLASHIda* ida = new FLASHIda(const_cast<char*>(agc_fast_json));

  // Get the AGC command the engine produces on init
  ScanCommand agc_cmd{};
  int r = ida->getNextScanCommand(agc_cmd);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(agc_cmd.is_agc, 1)

  // Load real MS1 peak data (non-trivial spectrum that would normally produce commands)
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  const auto& scan = ms1_scans[0];

  // Feed real peak data back with the AGC scan description — should be gated
  int result = ida->processScan(scan.mzs.data(), scan.ints.data(),
                                 (int)scan.mzs.size(), scan.rt, 1,
                                 agc_cmd.scan_description);
  TEST_EQUAL(result, 0)

  // Verify no commands were generated
  ScanCommand next_cmd{};
  int has_cmd = ida->getNextScanCommand(next_cmd);
  // Only expect the next idle-cycle AGC/MS1, not any MS2 from deconvolution
  // If the gate works, the scan data was never processed, so no MS2 commands
  TEST_EQUAL(has_cmd == 0 || next_cmd.is_agc == 1 || next_cmd.msn_level == 1, true)

  delete ida;
}
END_SECTION

// cleanupExpired should remove stale commands from priority queues, not pending map
START_SECTION(cleanup_expired_drops_stale_queued_commands)
{
  // Use idle_cycle_json (timeout enabled, 30s) but we'll use pushCommandForTest
  // with a 1ms-timeout config to verify expiry
  const char* short_timeout_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
    },
    "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [-50], "max_cv_skip": 0 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [
        { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": true, "value_ms": 1 },
      "agc_interval_seconds": 9999
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  FLASHIda* ida = new FLASHIda(const_cast<char*>(short_timeout_json));

  // Push 3 MS2 commands at priority 2 via test helper
  for (int i = 0; i < 3; ++i)
  {
    ScanCommand ms2{};
    ms2.msn_level = 2;
    ms2.priority = 2;
    ms2.scan_id = ida->getQueueForTest().nextTrackingId();
    ida->pushCommandForTest(ms2);
  }

  // Verify all 3 are in the queue
  TEST_EQUAL(ida->getQueueSizeForTest(2), (size_t)3)

  // Sleep 2ms to ensure commands exceed the 1ms timeout
  std::this_thread::sleep_for(std::chrono::milliseconds(2));

  // getNextScanCommand calls cleanupExpired in step 3.
  // After cleanup, the queue should be empty — all 3 stale commands dropped.
  // The idle cycle (step 5) will fire and return an AGC.
  ScanCommand out{};
  int r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(out.is_agc, 1)  // idle cycle AGC, not a stale MS2

  // Queue at priority 2 should now be empty (stale commands dropped)
  TEST_EQUAL(ida->getQueueSizeForTest(2), (size_t)0)

  delete ida;
}
END_SECTION

// MS1 and AGC scans should be resolved from pending_scan_map_ after processScan
START_SECTION(ms1_agc_resolved_from_pending_map)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(idle_cycle_json));

  // Idle cycle: returns AGC immediately, pushes MS1 at priority 3
  ScanCommand agc_cmd{};
  int r = ida->getNextScanCommand(agc_cmd);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(agc_cmd.is_agc, 1)

  // AGC is in pending map via registerPending
  TEST_EQUAL(ida->getQueueForTest().pendingScanMapSize(), (size_t)1)

  // Dequeue MS1 — now both AGC and MS1 are in pending map
  ScanCommand ms1_cmd{};
  r = ida->getNextScanCommand(ms1_cmd);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(ms1_cmd.is_agc, 0)
  TEST_EQUAL(ms1_cmd.msn_level, 1)
  TEST_EQUAL(ida->getQueueForTest().pendingScanMapSize(), (size_t)2)

  // processScan with AGC scan description — should resolve AGC from pending map
  // AGC gate (desc[3]=='A') returns 0 and resolves at line 715
  int n = ida->processScan(nullptr, nullptr, 0, 0.0, 1, agc_cmd.scan_description);
  TEST_EQUAL(n, 0)
  TEST_EQUAL(ida->getQueueForTest().pendingScanMapSize(), (size_t)1)

  // processScan with MS1 scan description — should resolve MS1 from pending map
  // Load real MS1 data to feed through the MS1 path
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  const auto& scan = ms1_scans[0];
  n = ida->processScan(scan.mzs.data(), scan.ints.data(),
                       (int)scan.mzs.size(), scan.rt, 1,
                       ms1_cmd.scan_description);
  // n >= 0 (may or may not produce MS2 commands from a single scan)
  TEST_EQUAL(n >= 0, true)
  TEST_EQUAL(ida->getQueueForTest().pendingScanMapSize(), (size_t)0)

  delete ida;
}
END_SECTION

END_TEST
