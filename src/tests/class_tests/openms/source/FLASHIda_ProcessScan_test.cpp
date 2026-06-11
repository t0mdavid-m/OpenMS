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
#include <map>
#include <set>
#include <thread>  // std::this_thread::sleep_for
#include <algorithm>  // std::max
#include <sstream>    // std::istringstream / std::ostringstream
#include <cstdio>     // std::remove

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

  // Config with MS3 targeting enabled: inclusion mode (target_mode=1) + cytC inclusion list
  // pins the precursor, plus MS3 selection via selection_strategy. Mirrors C# CT35.
  const char* ms3_mode1_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 1,
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
    "ms3": { "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE" },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "../../FlashIDA/test-data/configs/inclusion_cytc.txt", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 1 },
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

  // TSV file paths relative to the OpenMS build directory (CTest working dir)
  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_standard.txt";
  // 8 rich, co-eluting MS1 scans (raw 326-350) from a complex E. coli run, each with many distinct
  // selectable precursors (charge>=4, ChargeSNR>=1, mass 500-50k) -> the max_targets cap actually bites.
  const std::string ms1_ecoli_rich_tsv_path = "../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt";
  const std::string ms2_tsv_path = "../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt";
  const std::string ms2_tmt_tsv_path = "../../FlashIDA/test-data/spectra/ms2_quant_tmt.txt";
  const std::string ms1_cytc_tsv_path = "../../FlashIDA/test-data/spectra/ms1_cytc.txt";
  const std::string ms2_cytc_tsv_path = "../../FlashIDA/test-data/spectra/ms2_cytc_scan149.txt";
  // Fresh real Mode-2 MS3 CytC run (20250121): intact-cytC MS2 (HCD CE40, scan 57) with a rich b/y
  // ladder that FLASHExtender can validate -> real MS3. Extracted via test-scripts/prepare-test-data.py.
  const std::string ms2_cytc_fresh_tsv_path = "../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt";
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

  // Build a cap-test config from the max1_json template: set ms1.max_targets, optionally extend the
  // single HCD MS2 scan config to HCD+ETD, and enable per-scan results logging (runtime.scan_results_path).
  std::string buildCapConfig(int max_targets, bool etd, const std::string& results_path)
  {
    std::string j(max1_json);
    {
      const std::string key = "\"max_targets\": 1";  // sole occurrence (ms1 selection_strategy)
      auto p = j.find(key);
      if (p != std::string::npos) j.replace(p, key.size(), "\"max_targets\": " + std::to_string(max_targets));
    }
    if (etd)
    {
      const std::string hcd = R"({ "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 })";
      auto p = j.find(hcd);
      if (p != std::string::npos)
        j.replace(p, hcd.size(),
                  hcd + R"(,
        { "analyzer": "Orbitrap", "activation": "ETD", "collision_energy": 0, "resolution": 120000 })");
    }
    {
      const std::string files_key = "\"files\":";
      auto p = j.find(files_key);
      if (p != std::string::npos)
        j.insert(p, "\"runtime\": { \"scan_results_path\": \"" + results_path + "\" }, ");
    }
    return j;
  }

  // Parse a FLASHIda scan_results.tsv. Only MS1 scans are pushed in the cap test, so every data row is
  // an MS1 result row; return {commands_pushed, count(child_ids)} per row, in file order.
  std::vector<std::pair<int, int>> readMs1ResultRows(const std::string& path)
  {
    std::vector<std::pair<int, int>> rows;
    std::ifstream f(path);
    std::string header;
    if (! std::getline(f, header)) return rows;
    std::vector<std::string> cols;
    {
      std::istringstream hs(header);
      std::string c;
      while (std::getline(hs, c, '\t')) cols.push_back(c);
    }
    int ci_pushed = -1, ci_child = -1;
    for (int i = 0; i < (int)cols.size(); i++)
    {
      if (cols[i] == "commands_pushed") ci_pushed = i;
      else if (cols[i] == "child_ids") ci_child = i;
    }
    if (ci_pushed < 0 || ci_child < 0) return rows;
    std::string line;
    while (std::getline(f, line))
    {
      if (line.empty()) continue;
      std::vector<std::string> fields;
      {
        std::istringstream ls(line);
        std::string fld;
        while (std::getline(ls, fld, '\t')) fields.push_back(fld);
      }
      int pushed = (ci_pushed < (int)fields.size() && ! fields[ci_pushed].empty()) ? std::stoi(fields[ci_pushed]) : 0;
      int nchild = 0;
      if (ci_child < (int)fields.size() && ! fields[ci_child].empty())
      {
        // child_ids are space-separated (separator is outside the tracking-id alphabet, which
        // includes ';'); split on ' ' to match writeScanResultRow_'s child_str join.
        std::istringstream cs(fields[ci_child]);
        std::string cid;
        while (std::getline(cs, cid, ' ')) { if (! cid.empty()) nchild++; }
      }
      rows.emplace_back(pushed, nchild);
    }
    return rows;
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

// P4-U07: MS3 commands are pushed at priority 1 — hard positive.
// Drives the real processScan -> FragmentAnalysis/FLASHExtender -> buildMS3 chain on a FRESH real
// Mode-2 MS3 CytC run (20250121): the intact-cytC MS2 (HCD CE40, scan 57, ms2_cytc_fresh_scan57.txt)
// has a rich b/y ladder, and ms3_mode1_json now carries the matching M-starting cytC proteoform.
// Inclusion mode (target_mode=1 + inclusion_cytc.txt 12360 Da, ms1.max_targets=1 over ms1_cytc.txt)
// pins the MS2 command's precursor; FLASHExtender validates the proteoform -> MS3. The source run's
// MS3_Mappings.tsv (300+ sub-ppm b/y matches) is ground truth that the engine finds hits from this MS2.
START_SECTION(processScan_ms3_commands)
{
  auto ms1_scans = loadTsvScans(ms1_cytc_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_cytc_fresh_tsv_path);
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

  // Process MS2 return — inclusion-pinned cytC precursor -> b/y fragment matches -> MS3 targets
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

// P4-U03a: Mass exclusion — previously-selected precursors are NOT re-selected on an identical
// re-push (so they are never picked first). With tqscore_threshold=0.0 every selected mass (qscore
// is always > 0) arms the within-run exclusion map; in target_mode=0 only phase 0 runs and it SKIPS
// armed masses (PrecursorSelection.cpp:609-616), so the re-push selects strictly fewer commands
// (here zero — every selected mass armed). It does NOT re-pick them. ("Picked last, not first" is a
// different code path — phase 1 / mode 2 — out of scope here.)
START_SECTION(processScan_mass_exclusion)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  // standard_json uses tqscore_threshold=0.9; on ms1_standard.txt no selected mass reaches
  // qscore>0.9, so the within-run dynamic exclusion never arms. Lower it to 0.0 so EVERY
  // pass-1-selected mass arms the exclusion map (the gate at PrecursorSelection.cpp:654 is a
  // strict '>', and qscore is always > 0). Local copy only; shared standard_json is untouched.
  std::string excl_json(standard_json);
  {
    auto p = excl_json.find("\"tqscore_threshold\": 0.9");
    TEST_EQUAL(p != std::string::npos, true)
    excl_json.replace(p, std::string("\"tqscore_threshold\": 0.9").size(),
                      "\"tqscore_threshold\": 0.0");
  }
  FLASHIda* ida = new FLASHIda(const_cast<char*>(excl_json.c_str()));

  // Pass 1: every selected MS2 precursor m/z becomes "excluded" (qscore > 0 > threshold 0.0).
  int total_pass1 = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total_pass1 > 0, true)
  std::set<int> excluded;
  ScanCommand cmd{};
  for (int i = 0; i < total_pass1; i++)
  {
    if (ida->getNextScanCommand(cmd) != 1 || cmd.is_agc) break;
    if (cmd.msn_level == 2) excluded.insert((int)(cmd.stages[0].precursor_mz + 0.5));
  }
  ABORT_IF(excluded.empty())

  // Pass 2: re-push the SAME scans at the same RTs (within RT_window=180s), in dequeue order.
  int total_pass2 = pushAllScans(ida, ms1_scans);
  std::vector<int> pass2_order;
  ScanCommand out{};
  for (int i = 0; i < total_pass2; i++)
  {
    if (ida->getNextScanCommand(out) != 1 || out.is_agc) break;
    if (out.msn_level == 2) pass2_order.push_back((int)(out.stages[0].precursor_mz + 0.5));
  }
  // Dynamic exclusion SKIPS armed masses on the re-push, so the command count strictly drops
  // (every selected mass armed at threshold 0.0, so total_pass2 == 0 here).
  TEST_EQUAL(total_pass2 < total_pass1, true)
  // No previously-excluded mass is re-selected (hence none is picked first). Vacuous-safe: if pass2
  // is empty the loop simply runs zero times.
  for (int mz : pass2_order) TEST_EQUAL(excluded.count(mz) == 0, true)

  delete ida;
}
END_SECTION

// P4-U03b: Thresholded mass exclusion — ONLY masses selected with qscore > 0.1 arm the exclusion
// map and are therefore not re-selected on the re-push; lower-confidence masses stay eligible.
// Exercises the tqscore_threshold gate; contrast with processScan_mass_exclusion above (threshold
// 0.0, every selected mass armed).
START_SECTION(processScan_mass_exclusion_thresholded)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  std::string excl_json(standard_json);
  {
    auto p = excl_json.find("\"tqscore_threshold\": 0.9");
    TEST_EQUAL(p != std::string::npos, true)
    excl_json.replace(p, std::string("\"tqscore_threshold\": 0.9").size(),
                      "\"tqscore_threshold\": 0.1");
  }
  FLASHIda* ida = new FLASHIda(const_cast<char*>(excl_json.c_str()));

  // Pass 1: record, per integer-m/z, the MAX selection qscore (the value the engine gates on at
  // PrecursorSelection.cpp:654; cmd.qscore == pg.getQscore() since AllCharges=false).
  int total_pass1 = pushAllScans(ida, ms1_scans);
  TEST_EQUAL(total_pass1 > 0, true)
  std::map<int, float> max_q;
  ScanCommand cmd{};
  for (int i = 0; i < total_pass1; i++)
  {
    if (ida->getNextScanCommand(cmd) != 1 || cmd.is_agc) break;
    if (cmd.msn_level == 2)
    {
      int mz = (int)(cmd.stages[0].precursor_mz + 0.5);
      auto it = max_q.find(mz);
      if (it == max_q.end() || cmd.qscore > it->second) max_q[mz] = cmd.qscore;
    }
  }
  // Only masses whose selection qscore exceeded 0.1 arm the exclusion map.
  std::set<int> excluded;
  for (const auto& kv : max_q) if (kv.second > 0.1f) excluded.insert(kv.first);
  TEST_EQUAL(excluded.size() > 0, true)  // non-vacuous: there ARE qscore>0.1 masses to exclude

  // Pass 2: re-push the SAME scans, in dequeue order.
  int total_pass2 = pushAllScans(ida, ms1_scans);
  std::vector<int> pass2_order;
  ScanCommand out{};
  for (int i = 0; i < total_pass2; i++)
  {
    if (ida->getNextScanCommand(out) != 1 || out.is_agc) break;
    if (out.msn_level == 2) pass2_order.push_back((int)(out.stages[0].precursor_mz + 0.5));
  }
  // Exclusion strictly reduces the count; no qscore>0.1 mass is re-selected (hence not picked first).
  // Vacuous-safe over an empty pass2.
  TEST_EQUAL(total_pass2 < total_pass1, true)
  for (int mz : pass2_order) TEST_EQUAL(excluded.count(mz) == 0, true)

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

// max_targets caps precursor selection PER MS1 scan: children (MS2 commands) per scan rise with the cap,
// and equal (precursors selected) x (MS2 activations per precursor). Verified via the engine's own
// per-scan logging (runtime.scan_results_path -> commands_pushed / child_ids), across a 2 (activation) x
// 3 (cap) matrix. A *cumulative* command total across all scans is NOT a valid cap check: selection
// couples to persistent cross-scan state (mass_qscore_map_ score-drop skip, PrecursorSelection.cpp:647-659)
// that depends on the cap, so the running total is non-monotonic in max_targets. We therefore compare per
// scan, and across caps only on scan 0 (all engines start from an empty exclusion state). Precursor
// selection is independent of the ms2 array, so an HCD+ETD engine selects the same precursors as its HCD
// twin and logs exactly 2x the children (asserted on every scan).
START_SECTION(processScan_ms1_max_targets_cap)
{
  auto ms1_scans = loadTsvScans(ms1_ecoli_rich_tsv_path);
  ABORT_IF(ms1_scans.empty())
  const int n = (int)ms1_scans.size();

  // Run all MS1 scans through a fresh engine (max_targets, HCD or HCD+ETD) with per-scan results logging;
  // return children-per-MS1-scan and assert the logging join contract (commands_pushed == #child_ids).
  auto runCase = [&](int max_targets, bool etd) -> std::vector<int> {
    std::string path = std::string("p4u_maxcap_") + (etd ? "etd" : "hcd") + "_"
                     + std::to_string(max_targets) + "_results.tsv";
    std::remove(path.c_str());  // results stream opens in append mode -> ensure a fresh file
    std::string cfg = buildCapConfig(max_targets, etd, path);
    FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg.c_str()));
    for (const auto& s : ms1_scans)
      ida->processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 1,
                       ("scan_" + s.scan_id).c_str());
    delete ida;  // flush + close the results TSV
    std::vector<int> children;
    for (const auto& r : readMs1ResultRows(path)) { TEST_EQUAL(r.first, r.second) children.push_back(r.first); }
    std::remove(path.c_str());
    return children;
  };

  std::vector<int> A1 = runCase(1, false), A3 = runCase(3, false), A5 = runCase(5, false);
  std::vector<int> B1 = runCase(1, true),  B3 = runCase(3, true),  B5 = runCase(5, true);

  // 1) exactly one MS1 result row per scan
  TEST_EQUAL((int)A1.size(), n) TEST_EQUAL((int)A3.size(), n) TEST_EQUAL((int)A5.size(), n)
  TEST_EQUAL((int)B1.size(), n) TEST_EQUAL((int)B3.size(), n) TEST_EQUAL((int)B5.size(), n)
  ABORT_IF((int)A1.size() != n || (int)A3.size() != n || (int)A5.size() != n
           || (int)B1.size() != n || (int)B3.size() != n || (int)B5.size() != n)

  // 3) per-scan cap bound (with activation factor): children never exceed max_targets x activations
  for (int v : A1) TEST_EQUAL(v <= 1, true)
  for (int v : A3) TEST_EQUAL(v <= 3, true)
  for (int v : A5) TEST_EQUAL(v <= 5, true)
  for (int v : B1) TEST_EQUAL(v <= 2, true)
  for (int v : B3) TEST_EQUAL(v <= 6, true)
  for (int v : B5) TEST_EQUAL(v <= 10, true)

  // 4) more max_targets -> more children, on scan 0 (clean empty exclusion state for all engines)
  TEST_EQUAL(A1[0] <= A3[0], true) TEST_EQUAL(A3[0] <= A5[0], true)
  TEST_EQUAL(B1[0] <= B3[0], true) TEST_EQUAL(B3[0] <= B5[0], true)

  // 5) the cap bites at scan 0 (rich co-eluting E. coli precursors): exactly 1 at cap=1, strictly more
  //    as the cap grows (these rich scans offer >=9 selectable precursors each)
  TEST_EQUAL(A1[0], 1)
  TEST_EQUAL(A3[0] >= 2, true)
  TEST_EQUAL(A5[0] > A1[0], true)

  // 6) HCD+ETD logs exactly 2x the children of its HCD twin, on EVERY scan (the activation multiplier
  //    that made the old cumulative comparison invalid)
  for (int i = 0; i < n; i++)
  {
    TEST_EQUAL(B1[i], 2 * A1[i])
    TEST_EQUAL(B3[i], 2 * A3[i])
    TEST_EQUAL(B5[i], 2 * A5[i])
  }

  // 7) non-vacuous: a larger cap really selects > 1 precursor on some scan
  int maxA5 = 0; for (int v : A5) maxA5 = std::max(maxA5, v);
  TEST_EQUAL(maxA5 > 1, true)
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

  // With min_charge=99, no precursor passes the filter, so no deconvolution-derived
  // MS2 is produced. getNextScanCommand has no return-0 path: when the queue is empty
  // it emits an idle-cycle AGC and returns 1 (cf. processScan_commands_dequeued).
  ScanCommand cmd{};
  int result = ida.getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(cmd.is_agc, 1)
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
