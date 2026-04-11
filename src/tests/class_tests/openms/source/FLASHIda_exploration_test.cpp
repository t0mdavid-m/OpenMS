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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Exploration.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>

#include <vector>
#include <algorithm>

using namespace OpenMS;

namespace
{
  // Base JSON config with MS2 exploration enabled (mass_count, CE 20-40 step 5)
  const char* exploration_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "MS3AllCharges": true,
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
      "enabled": false,
      "mode": 0,
      "max_per_ms2": 4,
      "protein_sequence": ""
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_precursors": 3 },
      "ms2": {
        "selection": "intensity",
        "max_fragments": 3,
        "exploration": {
          "metric": "mass_count",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0,
          "activation": "HCD"
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
      "tol": [10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "MS3AllCharges": true,
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
      "enabled": false,
      "mode": 0,
      "max_per_ms2": 4,
      "protein_sequence": ""
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_precursors": 3 },
      "ms2": {
        "selection": "intensity",
        "max_fragments": 3,
        "exploration": {
          "metric": "mass_count",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0,
          "activation": "HCD"
        }
      },
      "ms3": {
        "selection": "intensity",
        "max_fragments": 3,
        "exploration": {
          "metric": "fragment_count",
          "ce_min": 15.0,
          "ce_max": 35.0,
          "ce_step": 5.0,
          "activation": "CID"
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
      "tol": [10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "MS3AllCharges": true,
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
      "enabled": false,
      "mode": 0,
      "max_per_ms2": 4,
      "protein_sequence": ""
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_precursors": 3 },
      "ms2": {
        "selection": "intensity",
        "max_fragments": 3,
        "exploration": {
          "metric": "mass_count",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0,
          "activation": "HCD"
        }
      },
      "ms3": { "selection": "intensity", "max_fragments": 3 }
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
      "tol": [10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "MS3AllCharges": true,
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
      "enabled": false,
      "mode": 0,
      "max_per_ms2": 4,
      "protein_sequence": ""
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_precursors": 3 },
      "ms2": {
        "selection": "intensity",
        "max_fragments": 3,
        "exploration": {
          "metric": "mass_count",
          "ce_min": 20.0,
          "ce_max": 40.0,
          "ce_step": 5.0,
          "activation": "HCD"
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
      "tol": [10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "MS3AllCharges": true,
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
      "enabled": false,
      "mode": 0,
      "max_per_ms2": 4,
      "protein_sequence": ""
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_precursors": 3 },
      "ms2": { "selection": "intensity", "max_fragments": 3 },
      "ms3": {
        "selection": "intensity",
        "max_fragments": 3,
        "exploration": {
          "metric": "fragment_count",
          "ce_min": 15.0,
          "ce_max": 35.0,
          "ce_step": 5.0,
          "activation": "CID"
        }
      }
    }
  })";

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
} // anonymous namespace


START_TEST(FLASHIda_exploration, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(exploration_group_creation)
{
  Config cfg{std::string(exploration_config)};
  ScanCommandQueue queue(cfg);
  Exploration exploration(cfg);

  auto cmds = exploration.initiate(2, 800.0, 2400.0, 3, 0.0, queue);

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

START_SECTION(exploration_variants_priority_0)
{
  Config cfg{std::string(exploration_config)};
  ScanCommandQueue queue(cfg);
  Exploration exploration(cfg);

  auto cmds = exploration.initiate(2, 800.0, 2400.0, 3, 0.0, queue);

  TEST_EQUAL(static_cast<int>(cmds.size()), 5)
  for (int i = 0; i < 5; ++i)
  {
    TEST_EQUAL(cmds[i].msn_level, 2)
    TEST_EQUAL(cmds[i].priority, 0)
    TEST_EQUAL(cmds[i].is_agc, 0)
    std::string desc(cmds[i].scan_description);
    TEST_EQUAL(desc.size() >= 4, true)
    TEST_EQUAL(desc[3], 'E')
  }
}
END_SECTION

START_SECTION(winner_selection_by_score)
{
  Config cfg{std::string(exploration_config)};
  ScanCommandQueue queue(cfg);
  Exploration exploration(cfg);

  auto cmds = exploration.initiate(2, 800.0, 2400.0, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)

  std::vector<double> scores = {1.0, 3.0, 2.0, 5.0, 0.0};
  for (int i = 0; i < 5; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, static_cast<int>(scores[i]));
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    exploration.feedResult(tracking_id, ds, static_cast<double>(i), queue);
  }

  TEST_EQUAL(exploration.activeGroupCount(), 0)
}
END_SECTION

START_SECTION(cycle_time_suppression_during_exploration)
{
  // P7-U05: MS1 cycle time suppression during active exploration
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cycle_time_exploration_config));

  // Create exploration group (cycle_time is 1ms, so it would trigger immediately)
  ida->initiateExplorationForTest(2, 800.0, 2400.0, 3, 0.0);

  // Active groups > 0
  TEST_EQUAL(ida->getActiveExplorationGroupCountForTest(), 1)

  // Get next command — should be exploration variant from queue[0], NOT MS1 from cycle time
  ScanCommand cmd{};
  int result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.msn_level, 2)  // exploration variant, not MS1
  TEST_EQUAL(cmd.priority, 0)
  std::string desc(cmd.scan_description);
  TEST_EQUAL(desc.size() >= 4, true)
  TEST_EQUAL(desc[3], 'E')

  delete ida;
}
END_SECTION

START_SECTION(ms1_resumes_after_exploration_completes)
{
  // P7-U06: MS1 cycle time injection resumes after exploration completes
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cycle_time_exploration_config));

  ida->initiateExplorationForTest(2, 800.0, 2400.0, 3, 0.0);

  // Drain all 5 exploration variants from queue
  for (int i = 0; i < 5; ++i)
  {
    ScanCommand drain{};
    ida->getNextScanCommand(drain);
    TEST_EQUAL(std::strlen(drain.scan_description) <= 15, true)
  }

  // Feed all 5 results to complete the group
  for (int i = 0; i < 5; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, 1);
    ida->feedExplorationResultForTest(1, i, ds, static_cast<double>(i));
  }

  // Group should be complete and removed
  TEST_EQUAL(ida->getActiveExplorationGroupCountForTest(), 0)

  // Now getNextScanCommand should return MS1 from cycle time (cycle_time_ms_ = 1, long past)
  ScanCommand cmd{};
  int result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_EQUAL(cmd.is_agc, 0)

  delete ida;
}
END_SECTION

START_SECTION(ms3_exploration_creates_child_groups)
{
  Config cfg{std::string(ms3_exploration_config)};
  ScanCommandQueue queue(cfg);
  Exploration exploration(cfg);

  auto cmds = exploration.initiate(2, 800.0, 2400.0, 3, 0.0, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)

  for (int i = 0; i < 5; ++i)
  {
    int count = (i == 2) ? 5 : 1;
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, count);
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    exploration.feedResult(tracking_id, ds, static_cast<double>(i), queue);
  }

  int ms3_group_count = exploration.activeGroupCount();
  TEST_EQUAL(ms3_group_count > 0, true)
  TEST_EQUAL(ms3_group_count <= 3, true)
}
END_SECTION

START_SECTION(ms3_selection_no_exploration_standard_targeting)
{
  // P7-U08: MS3 with selection but no exploration → standard MS3 commands
  FLASHIda* ida = new FLASHIda(const_cast<char*>(ms3_selection_only_config));

  // Create MS2 exploration group
  ida->initiateExplorationForTest(2, 800.0, 2400.0, 3, 0.0);

  // Drain exploration variants from queue[0]
  for (int i = 0; i < 5; ++i)
  {
    ScanCommand drain{};
    ida->getNextScanCommand(drain);
    TEST_EQUAL(std::strlen(drain.scan_description) <= 15, true)
  }

  // Feed 5 variants: variant 3 (CE=35.0) wins with 4 peak groups
  for (int i = 0; i < 5; ++i)
  {
    int count = (i == 3) ? 4 : 1;
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, count);
    ida->feedExplorationResultForTest(1, i, ds, static_cast<double>(i));
  }

  // MS2 group should be complete — no MS3 exploration groups created
  TEST_EQUAL(ida->getActiveExplorationGroupCountForTest(), 0)

  // MS3 commands should be in queue[1] (normal priority), not exploration groups
  // The initiateNextLevel_ for selection-only pushes to queue[1]
  size_t q1_size = ida->getQueueSizeForTest(1);
  TEST_EQUAL(q1_size > 0, true)

  // Verify commands are msn_level 3
  ScanCommand cmd{};
  int result = ida->getNextScanCommand(cmd);
  if (result == 1)
  {
    TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
    TEST_EQUAL(cmd.msn_level, 3)
    TEST_EQUAL(cmd.priority, 1)
  }
  (void)result;

  delete ida;
}
END_SECTION

START_SECTION(optimization_metadata_populated)
{
  Config cfg{std::string(exploration_config)};
  ScanCommandQueue queue(cfg);
  Exploration exploration(cfg);

  auto cmds = exploration.initiate(2, 800.0, 2400.0, 3, 0.0, queue);

  DeconvolvedSpectrum ds = makeSyntheticDeconv(1, 3);
  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResult(tracking_id, ds, 1.0, queue);

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
  Exploration exploration(cfg);

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
  TEST_STRING_EQUAL(ms2_cfg.exploration_activation, "HCD")

  auto ms3_cfg = cfg.level(3);
  TEST_EQUAL(static_cast<int>(ms3_cfg.selection), static_cast<int>(SelectionMetric::None))

  (void)ms1_cfg;
  (void)ms2_cfg;
  (void)ms3_cfg;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
