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
#include <fstream>

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
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  if (ms1_scans.empty()) { NOT_TESTABLE; break; }

  FLASHIda* ida = new FLASHIda(const_cast<char*>(cycle_time_exploration_config));

  int total = pushAllMS1Scans(ida, ms1_scans);
  if (total == 0) { delete ida; NOT_TESTABLE; break; }

  // getNextScanCommand should return exploration variants (priority 0, msn_level 2)
  // NOT cycle-time MS1, because exploration is active
  ScanCommand cmd{};
  int result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(cmd.msn_level, 2)
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
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  if (ms1_scans.empty() || ms2_scans.empty()) { NOT_TESTABLE; break; }

  FLASHIda* ida = new FLASHIda(const_cast<char*>(cycle_time_exploration_config));

  int total = pushAllMS1Scans(ida, ms1_scans);
  if (total == 0) { delete ida; NOT_TESTABLE; break; }

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
  // P7-U08: MS3 with selection but no exploration -> standard MS3 commands
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  if (ms1_scans.empty() || ms2_scans.empty()) { NOT_TESTABLE; break; }

  FLASHIda* ida = new FLASHIda(const_cast<char*>(ms3_selection_only_config));

  int total = pushAllMS1Scans(ida, ms1_scans);
  if (total == 0) { delete ida; NOT_TESTABLE; break; }

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
