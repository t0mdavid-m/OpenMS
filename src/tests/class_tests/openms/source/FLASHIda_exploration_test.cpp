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
// optimal_window_margin_ (ADR-0026 decision 2). Named DIRECTLY rather than taken on ScanCommandQueue.h's
// transitive pull: optimal_window_margin_has_one_definition guards this header's copy, so the guard has to
// name the header that owns it.
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/NotchSelection.h>

#include "FLASHIda_TestHelpers.h"  // ground-truth harness: ScanData/loadTsvScans/TSVFile/AcqResult/runInterleaved/runFullCycle/ExplResult/driveOneExplorationGroup/inclusionPinCytc
#include "FLASHIda_TestAccess.h"   // ExplorationTestAccess::feedResult/group (private-state access)

#include <vector>
#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <set>
#include <string>
#include <cstdio>
#include <cmath>
#include <tuple>   // PeakGroup::getMzRange's return type, read explicitly in single_isotope_charge_yields_non_degenerate_scan_range

using namespace OpenMS;

namespace
{
  // Carbonic anhydrase 2 (CA2) sequence (259 AA) — matches ms2_ca_hcd25_scan181.txt / ms2_ca_hcd45_scan185.txt
  // (real CA MS2; intact ~29006 Da / m/z 968.4916 / z30 precursor)
  const char* ca_sequence = "SHHWGYGKHNGPEHWHKDFPIANGERQSPVDIDTKAVVQDPALKPLALVYGEATSRRMVNNGHSFNVEYDDSQDKAVLKDGPLTGTYRLVQFHFHWGSSDDQGSEHTVDRKKYAAELHLVHWNTKYGDFGTAAQQPDGLAVVGVFLKVGDANPALQKVLDALDSIKTKGKSTDFPNFDPGSLLPNVLDYWTYPGSLTTPPLLESVTWIVLKEPISVSSQQMLKFRTLNFNAEGEPELLMLANWRPAQPLKNRQVRGFPK";
  const std::string ms2_ca_path = "../../FlashIDA/test-data/spectra/ms2_ca_hcd25_scan181.txt";

  // Minimal loadable config; `mode` and `rank_by` are the two spots a removed metric could be typed.
  std::string metricProbe(const std::string& mode, const std::string& rank_by)
  {
    return R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": ")" + rank_by + R"(", "max_precursors": 3 },
      "characterization": { "mode": ")" + mode + R"(", "protein_sequence": "PEPTIDER" },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
        "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 },
        "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25 }
      }
    })";
  }

  // ---- remaining_precursor rejection probes (ADR-0026) ---------------------------------------
  //
  // ADR-0026 binds a remaining_precursor pre-scan's scan range to the isolation window it reads.
  // That narrowing is what buys the sweep its speed, and it is sound only because the winner is
  // ALWAYS re-acquired by a production scan built from the un-overridden config. Two config shapes
  // break the reasoning, and Config::validate() rejects both -- the four sections that sit beside
  // ms3_protein_sequence_only_accepted pin them.
  //
  // Spliced, never find/replaced. The sections commented at :2108-2110 died on `invalid string
  // position` doing surgery on a finished literal; a splice has no position to compute, so the only
  // way the probes below can fail is the rejection itself. Every knob a section varies is a
  // parameter here, so no section reaches into another section's JSON.

  // A non-empty exploration.overrides map. Rejection A (Config.cpp:924) fires on EMPTINESS and is
  // checked BEFORE the multiplexing pair-checks, so any section aimed at a different throw has to
  // carry one -- otherwise Rejection A throws first and the section passes for the wrong reason.
  const std::string sweep_overrides = R"(, "overrides": { "analyzer": "Orbitrap" })";

  // A remaining_precursor exploration block, ready to splice into either decision section. Pass ""
  // for `overrides` to probe Rejection A and `sweep_overrides` otherwise. One CE range serves both
  // levels: needsCollisionEnergy is true for HCD and CID alike, and Config.cpp:1019-1022 asks only
  // for ce_min < ce_max.
  std::string rpSweep(const std::string& overrides)
  {
    return R"("exploration": { "metric": "remaining_precursor", "ce_min": 20.0, "ce_max": 40.0, "ce_step": 5.0)"
           + overrides + R"( })";
  }

  // The characterization body that makes characterization.exploration reachable at level 3: the MS3
  // gate (mode != off) plus the sequence Config.cpp:968-972 then requires. Trailing comma included,
  // so a caller appends its own keys directly.
  const std::string ms3_on = R"("mode": "ambiguity", "protein_sequence": "PEPTIDER", )";

  // Minimal loadable config. `ps_extra` is appended inside precursor_selection (lead with a comma),
  // `charact` REPLACES the characterization body, and `with_ms3` adds the level-3 scan config that a
  // non-off mode requires (Config.cpp:978 -- without it Exploration::initiateNextLevel OOB-reads
  // scans[0], so its absence is a crash rather than a no-op).
  std::string sweepProbe(const std::string& ps_extra = "",
                         const std::string& charact = R"("mode": "off")",
                         bool with_ms3 = false)
  {
    return std::string(R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": "qscore", "max_precursors": 3)") + ps_extra + R"( },
      "characterization": { )" + charact + R"( },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
        "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 })"
      + (with_ms3 ? R"(,
        "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25 })" : "")
      + R"(
      }
    })";
  }

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
    "max_precursors": 3,
    "exploration": {
      "metric": "mass_count",
      "ce_min": 20.0,
      "ce_max": 40.0,
      "ce_step": 5.0
    }
  },
  "characterization": {
    "mode": "off",
    "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
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

  // Config with MS2 exploration + MS3 exploration (fragment_count)
  const char* ms3_exploration_config = R"({
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
    "max_precursors": 3,
    "exploration": {
      "metric": "mass_count",
      "ce_min": 20.0,
      "ce_max": 40.0,
      "ce_step": 5.0
    }
  },
  "characterization": {
    "mode": "ambiguity",
    "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "max_targets": 3,
    "exploration": {
      "metric": "fragment_count",
      "ce_min": 15.0,
      "ce_max": 35.0,
      "ce_step": 5.0
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
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "CID",
      "collision_energy": 25,
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

  // Config with MS2 exploration + MS3 selection only (no exploration)
  const char* ms3_selection_only_config = R"({
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
    "max_precursors": 3,
    "exploration": {
      "metric": "mass_count",
      "ce_min": 20.0,
      "ce_max": 40.0,
      "ce_step": 5.0
    }
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

  // Config with remaining_precursor exploration metric at MS2.
  //
  // The exploration.overrides map below MUST be non-empty: ADR-0026 rejects "remaining_precursor" with
  // an empty map at config load (Config.cpp:924-930, at BOTH levels). Such a sweep's pre-scans are
  // narrowed to the isolation window and then thrown away, so the schema now forces the author to
  // declare the settings they run at, which is what guarantees the ADR-0020 gate-#1 re-acquisition.
  // "Orbitrap" is deliberately a NO-OP: it is what ms_settings.ms2.analyzer already says, so
  // Exploration::initiate's base_config.applyOverrides (Exploration.cpp:128) writes back the value the
  // base config already carries and every fixture behaviour is unchanged. Same inert-override trick as
  // the injection at :1659.
  //
  // THE ONE BRANCH A NON-EMPTY MAP OPENS is ADR-0020 gate #1 -- `!level_config.overrides.empty() ||
  // measuring_ms3_sweep` at Exploration.cpp:748 -- and it is unreached here. THREE sections driven off
  // this fixture complete their group, and all three reach completion through the baseline_failed abort
  // at Exploration.cpp:434-437, which takes the `group.baseline_failed || best_idx < 0` arm at :683 and
  // returns before the gate is ever evaluated:
  //   * remaining_precursor_empty_baseline_aborts_group and
  //     remaining_precursor_inflight_child_after_abort_is_noop -- aborts by construction, on a raw
  //     baseline whose peaks all sit outside the isolation window.
  //   * remaining_precursor_score_no_raw_data -- NOT an abort section by name or by intent. It feeds all
  //     SIX variants and asserts activeGroupCount() == 0, and it lands on the same arm only because
  //     ExplorationTestAccess::feedResult hands feedResultImpl_ mzs == ints == nullptr, so
  //     precursorWindowIntensity_ short-circuits to 0.0 (Exploration.cpp:1214-1218) and the baseline
  //     reads as empty.
  // Every other section stops feeding partway, so no group on this fixture has ever selected a winner.
  //
  // LIVE HAZARD, not a historical note. Give remaining_precursor_score_no_raw_data -- or any section
  // added here -- a raw-array baseline with in-window signal and a winner IS selected, and gate #1 is
  // then not merely reached but TRUE, because this fixture's overrides map is non-empty by
  // construction. Gate #1 and the MS2 cascade `else if (group.msn_level < 3)` (Exploration.cpp:816) are
  // mutually exclusive arms of ONE if/else-if chain, so the cascade would be REPLACED by a production-
  // MS2 re-acquisition, never added to. With characterization.mode "off" the cascade yields nothing
  // here anyway (level 3 selection is None, so initiateNextLevel returns empty), which is exactly what
  // makes the swap quiet: today its only symptom is an extra command in info.commands, and the lost
  // cascade surfaces only once someone gives this fixture an ms_settings.ms3 block.
  const char* remaining_precursor_config = R"({
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
    "max_precursors": 3,
    "exploration": {
      "metric": "remaining_precursor",
      "ce_min": 20.0,
      "ce_max": 40.0,
      "ce_step": 5.0,
      "overrides": { "analyzer": "Orbitrap" }
    }
  },
  "characterization": {
    "mode": "off",
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

  // Config with remaining_precursor exploration metric at MS3.
  //
  // Same ADR-0026 requirement as the MS2 fixture above -- the rejection at Config.cpp:924-930 loops over
  // every level, so characterization.exploration.overrides must be non-empty too. "Orbitrap" is again a
  // deliberate no-op (ms_settings.ms3.analyzer already says Orbitrap), so this fixture drives exactly the
  // behaviour it did before. Note the isolation-window narrowing itself is NOT gated on this map -- it is
  // suppressed only by an explicit first_mass/last_mass override (Exploration.cpp:255-256) -- so adding an
  // analyzer key does not perturb ms3_remaining_precursor_isolation_width's window arithmetic.
  const char* ms3_remaining_precursor_config = R"({
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
    "mode": "off",
    "max_targets": 3,
    "exploration": {
      "metric": "remaining_precursor",
      "ce_min": 20.0,
      "ce_max": 40.0,
      "ce_step": 5.0,
      "overrides": { "analyzer": "Orbitrap" }
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
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "CID",
      "collision_energy": 25,
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

  // Config with mass_count exploration at MS3 and NO overrides key -- the fixture that isolates ADR-0020
  // gate #2 now that ADR-0026 has made "remaining_precursor + empty overrides" unrepresentable.
  //
  // Both of ADR-0020's re-acquisition reasons are MEASURING-metric-driven, but only gate #2 keys on the
  // metric; gate #1 keys on a non-empty overrides map (Exploration.cpp:748). Since ADR-0026 forces every
  // authored remaining_precursor sweep to carry overrides, remaining_precursor can no longer exercise
  // gate #2 in isolation -- gate #1 would fire first and the assertion would prove nothing. mass_count is
  // the other measuring metric (isMeasuringMetric, Config.h:111-114), is untouched by either ADR-0026
  // rejection, and so is now the only way to reach gate #2 with an empty map.
  //
  // Written as a SELF-CONTAINED literal rather than a find/replace over ms3_exploration_config: the note
  // at :2108-2110 records four sections that died on `invalid string position` when the token they
  // searched for was renamed out from under them.
  const char* ms3_mass_count_config = R"({
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
    "mode": "off",
    "max_targets": 3,
    "exploration": {
      "metric": "mass_count",
      "ce_min": 20.0,
      "ce_max": 40.0,
      "ce_step": 5.0
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
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "CID",
      "collision_energy": 25,
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

  // Config with cycle time enabled + exploration
  const char* cycle_time_exploration_config = R"({
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
      "enabled": true,
      "value_ms": 1
    },
    "scan_timeout": {
      "enabled": false,
      "value_ms": 30000
    },
    "agc_interval_seconds": 999999
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
    "max_precursors": 3,
    "exploration": {
      "metric": "mass_count",
      "ce_min": 20.0,
      "ce_max": 40.0,
      "ce_step": 5.0
    }
  },
  "characterization": {
    "mode": "off",
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

  // Config: no MS2 exploration, MS3 exploration enabled
  const char* no_ms2_expl_ms3_expl_config = R"({
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
    "protein_sequence": "MKWVTFISLLLLFSSAYSRGVFRR",
    "max_targets": 3,
    "exploration": {
      "metric": "fragment_count",
      "ce_min": 15.0,
      "ce_max": 35.0,
      "ce_step": 5.0
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
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "CID",
      "collision_energy": 25,
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


  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_standard.txt";
  const std::string ms2_tsv_path = "../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt";

  // Config with 3-entry tol and MS2 exploration tolerance override
  const char* exploration_tolerance_config = R"({
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
      20
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
    "max_precursors": 3,
    "exploration": {
      "metric": "mass_count",
      "ce_min": 20.0,
      "ce_max": 40.0,
      "ce_step": 5.0,
      "tolerance_ppm": 15.0
    }
  },
  "characterization": {
    "mode": "off",
    "max_targets": 3,
    "exploration": {
      "metric": "mass_count",
      "ce_min": 20.0,
      "ce_max": 40.0,
      "ce_step": 5.0
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
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "CID",
      "collision_energy": 25,
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

  // Config with ETD exploration for activation type wiring test
  const char* etd_exploration_config = R"({
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
    "max_precursors": 3,
    "exploration": {
      "metric": "mass_count",
      "reaction_time_min": 5.0,
      "reaction_time_max": 15.0,
      "reaction_time_step": 5.0,
      "activations": [
        "ETD"
      ]
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
      "activation": "ETD",
      "reaction_time": 10.0,
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
  auto cmds = exploration.initiate(2, pg, 3, queue);

  TEST_EQUAL(exploration.activeGroupCount(), 1)

  auto group = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group.group_id, 1)
  TEST_EQUAL(group.msn_level, 2)
  TEST_EQUAL(group.complete, false)
  TEST_EQUAL(group.winner_index, -1)
  TEST_EQUAL(static_cast<int>(group.exploration_metric),
             static_cast<int>(ExplorationMetric::MassCount))

  // baseline-on-all (#18) prepends one CE-0 baseline variant (variant_index -1, is_baseline) at [0].
  // Assert the 5 real CE-sweep variants by filtering the baseline out — robust to the baseline rule.
  TEST_EQUAL(static_cast<int>(group.variants.size()), 6)
  TEST_EQUAL(group.variants[0].is_baseline, true)
  TEST_EQUAL(group.variants[0].variant_index, -1)
  TEST_REAL_SIMILAR(group.variants[0].collision_energy, 0.0)

  const double expected_ce[5] = {20.0, 25.0, 30.0, 35.0, 40.0};
  int real_count = 0;
  for (const auto& v : group.variants)
  {
    if (v.is_baseline) continue;
    TEST_EQUAL(v.received, false)
    TEST_EQUAL(v.variant_index, real_count)
    TEST_REAL_SIMILAR(v.collision_energy, expected_ce[real_count])
    ++real_count;
  }
  TEST_EQUAL(real_count, 5)

  TEST_EQUAL(static_cast<int>(cmds.size()), 6)  // baseline + 5 CE-sweep variants

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
  auto cmds = exploration.initiate(2, pg, 3, queue);

  TEST_EQUAL(static_cast<int>(cmds.size()), 6)  // baseline + 5 CE-sweep variants (all priority-2 'E')
  for (int i = 0; i < 6; ++i)
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
  ScanCommand ms2_ctx = queue.buildMS2(precursor_pg, 3, cfg.level(2).scans[0], 2, 0);

  auto fragment_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  ms2_ctx.faims_cv = -50.0;  // Item 1: CV travels via the context
  auto cmds = exploration.initiate(3, fragment_pg, 2, queue, nullptr, &ms2_ctx, 'y', 5);

  TEST_EQUAL(static_cast<int>(cmds.size()), 6)  // baseline + 5 CE-sweep variants (all MS3 buildMS3)

  for (int i = 0; i < 6; ++i)
  {
    TEST_EQUAL(cmds[i].msn_level, 3)
    TEST_EQUAL(cmds[i].num_stages, 2)
    TEST_EQUAL(cmds[i].priority, 1)
    TEST_REAL_SIMILAR(cmds[i].faims_cv, -50.0)
  }

  TEST_EQUAL(cmds[0].stages[0].charge_state, 3)

  auto group = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group.msn_level, 3)
  TEST_EQUAL(group.originating_cmd.num_stages > 0, true)

  // MS3 exploration descriptions must include fragment ion info
  for (int i = 0; i < 6; ++i)
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
  // makeSyntheticDeconv() below produces empty (mass-0) PeakGroups, so the MS3 FragmentCount path
  // (calibrated MS3FragmentMatcher) would match 0 fragments and score every variant 0. Swap the MS3
  // exploration metric to mass_count so the score is spec.size(): variant 3 (8 peak groups) WINS.
  // F5: the completing feedResult (variant 4, the LAST fed, 3 peak groups) now reports ITS OWN metrics
  // (score 3.0, variant_index 4) — it is no longer overwritten with the winner's. The winner (variant 3)
  // is identified by the new winner_tracking_id, which equals variant 3's encoded command id.
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
  ScanCommand ms2_ctx = queue.buildMS2(precursor_pg, 3, cfg.level(2).scans[0], 2, 0);

  auto fragment_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  ms2_ctx.faims_cv = -50.0;  // Item 1: CV travels via the context
  auto cmds = exploration.initiate(3, fragment_pg, 2, queue, nullptr, &ms2_ctx);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)  // baseline (cmds[0]) + 5 CE-sweep variants
  TEST_EQUAL(exploration.activeGroupCount(), 1)

  // cmds[0] is the CE-0 baseline (variant_index -1, skipped in winner selection); cmds[1..5] are the
  // 5 real variants (variant_index 0..4) with peak counts {2,4,6,8,3}. Feeding ALL 6 completes the group;
  // the winner is the real variant with 8 peaks (variant_index 3 == cmds[4]).
  std::vector<int> peak_counts = {0, 2, 4, 6, 8, 3};  // baseline + 5 reals
  Exploration::FeedResultInfo last_info;
  for (int i = 0; i < 6; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, peak_counts[i]);
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    last_info = ExplorationTestAccess::feedResult(exploration,tracking_id, ds, static_cast<double>(i), queue);
  }

  TEST_EQUAL(exploration.activeGroupCount(), 0)
  // F5: last_info is the COMPLETING variant (variant_index 4, 3 peak groups) reporting its OWN metrics,
  // not the winner's (variant_index 3, score 8.0) — the winner-overwrite was removed.
  TEST_REAL_SIMILAR(last_info.metric.score, 3.0)            // completing variant's own score (was 8.0)
  TEST_EQUAL(last_info.group.variant_index, 4)             // its own index, not the winner's 3
  // The winner stays identifiable via winner_tracking_id == variant_index 3's encoded command id (cmds[4]).
  // (getGroup throws here — the group was erased on completion — so read it from last_info.)
  TEST_EQUAL(last_info.group.winner_tracking_id.empty(), false)
  TEST_EQUAL(last_info.group.winner_tracking_id, std::string(cmds[4].scan_description).substr(0, 3))

  // MS3 command format (num_stages=2, priority=1) is verified in
  // ms3_exploration_variants_use_buildMS3 via the initiate() path.
  //
  // ADR-0020: mass_count is a MEASURING metric, so a completed MS3 group re-acquires even though
  // overrides is empty. This used to assert nothing and merely comment "returns 0 commands" — which
  // was true then, and is the exact defect ADR-0020 fixes (a measuring MS3 sweep left no evidence at
  // all). The cascade branch is NOT what produced this: an MS3 group is terminal (msn_level < 3 is
  // false), so the one command here can only be the production re-acquisition.
  TEST_EQUAL(static_cast<int>(last_info.commands.size()), 1)
  // Release build: no _ITERATOR_DEBUG_LEVEL and no OPENMS_PRECONDITION. Indexing an EMPTY
  // `commands` here is a read at address 0x4 that kills the binary and takes every later section
  // with it -- which is exactly what an empty-baseline group abort produced: one wrong metric
  // surfaced as a score of 0 in this section and a SegFault that hid 44 of the file's 48.
  // One failure must report as one failed line.
  ABORT_IF(last_info.commands.empty())
  TEST_EQUAL(last_info.commands[0].msn_level, 3)
  TEST_EQUAL(last_info.commands[0].num_stages, 2)
  TEST_EQUAL(last_info.commands[0].priority, 1)
  // Production, not a variant: it must return on the REGULAR MS3 path (marker 'R', not 'E'), which is
  // the whole point — that path runs the calibrated matcher against the live winner.
  TEST_EQUAL(std::string(last_info.commands[0].scan_description)[3], 'R')
  // ...and the regular path resolves its parent MS2 context out of the cache, so the caller must be
  // handed a seed for it. A miss there silently skips identification (the FLASHIda.cpp comment on the
  // MS2->MS3 cascade records exactly that failure), which would reintroduce the bug one layer down.
  TEST_EQUAL(static_cast<int>(last_info.ms2_context_cache.size()), 1)
  ABORT_IF(last_info.ms2_context_cache.empty())
  TEST_EQUAL(last_info.ms2_context_cache[0].first, last_info.commands[0].scan_id)
  // The production CE is the WINNER's, not the level default (25) and not the completing variant's
  // (cmds[5], CE 35). mass_count scores spec.size(), so variant_index 3 — 8 peak groups, cmds[4],
  // CE 30 on the 15/20/25/30/35 grid — won. Reading it off cmds[4] rather than hard-coding 30 keeps
  // the assertion true if the fixture's sweep grid is ever retuned.
  TEST_REAL_SIMILAR(last_info.commands[0].stages[1].collision_energy,
                    cmds[4].stages[1].collision_energy)
}
END_SECTION

START_SECTION(ms3_measuring_metric_always_reacquires_without_overrides)
{
  // ADR-0020, the empty-overrides half. A MEASURING metric scores from bulk signal and never
  // matches fragments, so its pre-scans leave NO identification behind: identification_result stays
  // default-constructed (its single write site is the FragmentCount batch re-score), the inline fold
  // at group completion is gated on that field being non-empty, and so a completed MS3 sweep used to
  // contribute nothing whatsoever. Real-run evidence: 66 MS3 scans acquired, 3 CE winners chosen,
  // zero MS3 rows in identification.tsv and zero MS3 evidence in pooled_identification.tsv.
  //
  // The fix is a production re-acquisition at the winning CE, which returns on the regular MS3 path
  // and IS identified. Asserted here at the seam that decides it: commands emitted at completion.
  //
  // THE METRIC IS mass_count, NOT remaining_precursor, AND THAT IS FORCED BY ADR-0026. The gate is
  // `!level_config.overrides.empty() || measuring_ms3_sweep` (Exploration.cpp:748) -- a disjunction, so
  // proving the SECOND term requires the first to be false, i.e. an empty overrides map. ADR-0026 now
  // rejects exactly that pairing for remaining_precursor at config load (Config.cpp:924-930), so this
  // section's old fixture no longer loads, and adding overrides to it would have satisfied gate #1 and
  // left the assertion below passing for the wrong reason. mass_count is the other measuring metric
  // (Config.h:111-114) and neither ADR-0026 rejection touches it, so it is the only remaining way to
  // reach gate #2 alone. The metrics are interchangeable HERE precisely because gate #2 keys on
  // isMeasuringMetric and not on a specific enumerator.
  //
  // Feeding pre-deconvolved spectra (ExplorationTestAccess) rather than raw peaks is the matching
  // change: mass_count scores spec.size(), so the peak-count vector below IS the score vector, whereas
  // remaining_precursor had to be driven through raw isolation-window intensities. It also drops the
  // empty-baseline hazard entirely -- the baseline_failed abort is gated on RemainingPrecursor
  // (Exploration.cpp:434-435), so a zero-intensity baseline is simply a score-0 variant here.
  Config cfg{std::string(ms3_mass_count_config)};
  TEST_EQUAL(cfg.level(3).overrides.empty(), true)  // the precondition: no overrides to trigger the old gate
  TEST_EQUAL(isMeasuringMetric(cfg.level(3).exploration), true)

  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto fragment_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  ScanCommand ms2_ctx = queue.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3,
                                       cfg.level(2).scans[0], 2, 0);
  auto cmds = exploration.initiate(3, fragment_pg, 2, queue, nullptr, &ms2_ctx);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)  // baseline + CE 20,25,30,35,40

  // cmds[0] is the CE-0 baseline (variant_index -1, skipped in winner selection); cmds[1..5] are the
  // 5 real variants. mass_count scores spec.size(), so the counts below are the scores directly and
  // the winner is the 8-peak variant at cmds[4].
  std::vector<int> peak_counts = {0, 2, 4, 6, 8, 3};  // baseline + 5 reals
  Exploration::FeedResultInfo last_info;
  for (int i = 0; i < 6; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, peak_counts[i]);
    int tid = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    last_info = ExplorationTestAccess::feedResult(exploration, tid, ds, static_cast<double>(i), queue);
  }

  TEST_EQUAL(exploration.activeGroupCount(), 0)                      // group completed
  TEST_EQUAL(last_info.group.winner_tracking_id.empty(), false)      // a winner was chosen
  TEST_EQUAL(last_info.group.winner_tracking_id, std::string(cmds[4].scan_description).substr(0, 3))

  // THE ASSERTION THIS TEST EXISTS FOR: exactly one production MS3, despite empty overrides.
  TEST_EQUAL(static_cast<int>(last_info.commands.size()), 1)
  // Release build: no _ITERATOR_DEBUG_LEVEL and no OPENMS_PRECONDITION. Indexing an EMPTY
  // `commands` here is a read at address 0x4 that kills the binary and takes every later section
  // with it -- which is exactly what an empty-baseline group abort produced: one wrong metric
  // surfaced as a score of 0 in this section and a SegFault that hid 44 of the file's 48.
  // One failure must report as one failed line.
  ABORT_IF(last_info.commands.empty())
  TEST_EQUAL(last_info.commands[0].msn_level, 3)
  TEST_EQUAL(last_info.commands[0].num_stages, 2)
  TEST_EQUAL(std::string(last_info.commands[0].scan_description)[3], 'R')
  TEST_EQUAL(static_cast<int>(last_info.ms2_context_cache.size()), 1)
  ABORT_IF(last_info.ms2_context_cache.empty())
  TEST_EQUAL(last_info.ms2_context_cache[0].first, last_info.commands[0].scan_id)
}
END_SECTION

START_SECTION(ms3_reading_metric_does_not_reacquire_without_overrides)
{
  // The other side of ADR-0020, and the reason the rule is keyed on the METRIC rather than made
  // unconditional: FragmentCount is a READING metric. It matches every variant in order to score it,
  // so the batch re-score populates identification_result and the inline fold at group completion
  // already carries the evidence. Re-acquiring would spend one extra MS3 per target for data the
  // group already holds.
  //
  // Pinning this asymmetry is what stops a later "simplify: always re-acquire" from silently
  // doubling the MS3 cost of every fragment_count method.
  Config cfg{std::string(ms3_exploration_config)};
  TEST_EQUAL(cfg.level(3).overrides.empty(), true)
  TEST_EQUAL(isMeasuringMetric(cfg.level(3).exploration), false)  // FragmentCount == reading

  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto fragment_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  ScanCommand ms2_ctx = queue.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3,
                                       cfg.level(2).scans[0], 2, 0);
  auto cmds = exploration.initiate(3, fragment_pg, 2, queue, nullptr, &ms2_ctx);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  std::vector<int> peak_counts = {0, 2, 4, 6, 8, 3};
  Exploration::FeedResultInfo last_info;
  for (int i = 0; i < 6; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, peak_counts[i]);
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    last_info = ExplorationTestAccess::feedResult(exploration, tracking_id, ds,
                                                  static_cast<double>(i), queue);
  }

  TEST_EQUAL(exploration.activeGroupCount(), 0)
  // No production scan: an MS3 group is terminal so the cascade branch cannot fire either, leaving
  // the inline fold as the sole consumer — exactly as before ADR-0020.
  TEST_EQUAL(static_cast<int>(last_info.commands.size()), 0)
  TEST_EQUAL(static_cast<int>(last_info.ms2_context_cache.size()), 0)
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
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)  // baseline (cmds[0]) + 5 CE-sweep variants

  std::vector<double> scores = {0.0, 1.0, 3.0, 2.0, 5.0, 0.0};  // baseline + 5 reals
  Exploration::FeedResultInfo last_info;
  for (int i = 0; i < 6; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, static_cast<int>(scores[i]));
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    last_info = ExplorationTestAccess::feedResult(exploration,tracking_id, ds, static_cast<double>(i), queue);
  }

  // mass_count metric -> remaining_ratio should be -1.0 (N/A)
  TEST_REAL_SIMILAR(last_info.metric.remaining_ratio, -1.0)
  TEST_EQUAL(exploration.activeGroupCount(), 0)
}
END_SECTION

START_SECTION(cycle_time_suppression_during_exploration)
{
  // P7-U05: MS1 cycle-time injection is suppressed WHILE an exploration group is active. Drive the FULL
  // acquisition through the general harness (runInterleaved) and scope the assertion to the exploration-ACTIVE
  // window via AcqResult.all_active (the engine's exploration_active_ flag captured at each dequeue). A
  // cycle-time MS1 (msn_level 1, priority 0, is_agc 0) injected while exploration is INACTIVE — before a group
  // forms or between groups — is CORRECT engine behavior, so suppression is asserted ONLY where all_active is true.
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  FLASHIda* ida = new FLASHIda(const_cast<char*>(cycle_time_exploration_config));
  AcqResult a = runInterleaved(ida, ms1_scans, std::vector<ScanData>{ms2_scans[0]});

  // The only prio-0 non-AGC MS1 the engine emits is the cycle-time injection (idle survey MS1 is priority 3;
  // AGC is is_agc==1). FAIMS is off in this config, so there is no CV-transition prio-0 MS1 to confuse it.
  auto isCycleTimeMS1 = [](const ScanCommand& c) {
    return c.msn_level == 1 && c.priority == 0 && c.is_agc == 0;
  };

  // Anti-vacuous: a group must have formed (>=1 active-window command) and a variant must have fired.
  int first_active = -1;
  bool any_e = false;
  for (int i = 0; i < (int)a.all_cmds.size(); ++i)
  {
    if (first_active < 0 && a.all_active[i]) first_active = i;
    std::string d(a.all_cmds[i].scan_description);
    if (d.size() >= 4 && d[3] == 'E') any_e = true;
  }
  if (first_active < 0 || !any_e) { delete ida; }
  ABORT_IF(first_active < 0 || !any_e)

  // Suppression invariant, scoped to the active window: no cycle-time MS1 is dequeued while exploration is active.
  for (int i = 0; i < (int)a.all_cmds.size(); ++i)
    if (a.all_active[i]) TEST_EQUAL(isCycleTimeMS1(a.all_cmds[i]), false);

  // Original local intent recovered (now correctly scoped): the FIRST command dequeued while the group is
  // active is the suppressed exploration variant (priority 2, msn_level 2, scan_description[3]=='E').
  TEST_EQUAL(a.all_cmds[first_active].msn_level, 2)
  TEST_EQUAL(a.all_cmds[first_active].priority, 2)
  std::string first_desc(a.all_cmds[first_active].scan_description);
  TEST_EQUAL(first_desc.size() >= 4, true)
  TEST_EQUAL(first_desc[3], 'E')

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

  // Full interleaved drive over the one canonical contract (runInterleaved), replacing the
  // hand-rolled bootstrap + `while (getNextScanCommand==1) { ... else break; }` variant-drain +
  // MS2-feed + trailing MS1-probe loop. Feed the SAME MS1 survey list the old bootstrap consumed
  // (engine-emitted ids); one survey forms the exploration group, runInterleaved drains its
  // variants + feeds MS2 back so the group completes, and -- because cycle-time MS1 injection is
  // NOT suppressed once exploration is done -- the engine resumes MS1, which the driver feeds as a
  // further recorded MS1 command (idle>=3 then bounds the drive).
  //   - exploration variants drained + MS2 fed back  -> recorded in a.ms2_cmds (the 'E' variants)
  //   - cycle-time MS1 resumes after the group done   -> >= 2 recorded MS1 cmds in a.ms1_cmds
  AcqResult a = runInterleaved(ida, ms1_scans, std::vector<ScanData>{ms2_scans[0]});
  delete ida;

  // The exploration group must have formed (a forming MS1) and produced MS2 exploration variants.
  ABORT_IF(a.ms1_cmds.empty())
  TEST_EQUAL(a.ms2_cmds.size() >= 1, true)
  bool drained_exploration_variant = false;
  for (const auto& c : a.ms2_cmds)
  {
    std::string desc(c.scan_description);
    if (desc.size() >= 4 && desc[3] == 'E') { drained_exploration_variant = true; break; }
  }
  TEST_EQUAL(drained_exploration_variant, true)

  // MS1 cycle-time injection resumes after exploration completes: more than one MS1 survey was
  // emitted + fed (the post-exploration MS1 was injected only because exploration was finished).
  TEST_EQUAL(a.ms1_cmds.size() >= 2, true)
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

START_SECTION(ms2_exploration_fixed_ce_ms3_carries_parent_context)
{
  // Regression pin for the dropped-context defect. Same config as the section above (MS2 exploration,
  // overrides EMPTY, MS3 selection on, MS3 exploration off), but asserting the step past dispatch: each
  // fixed-CE MS3 the winner dispatches must arrive with its parent MS2 context registered.
  //
  // Why the section above did not catch it: it asserts an MS3 command was EMITTED. The commands were
  // always emitted — feedResult just dropped NextLevelResult::ms3_contexts on the floor, so the MS3s
  // returned on the REGULAR path, missed ms2_context_cache_, and skipped the whole identification block
  // (no identification.tsv row, no tracker feedScan/foldMs3, nothing in pooled_identification). The
  // regular MS2->MS3 path seeded the same cache correctly, which is why only the exploration route broke.
  //
  // The drive passes a non-null but EMPTY ms3_ion_map: MS3 commands are dequeued and recorded but never
  // fed back (the harness never fabricates MS3 data), so the cache entries are still outstanding when we
  // read them — a returning MS3 erases its own entry.
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  std::string cfg_str = inclusionPinCytc(ms3_selection_only_config);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  const std::map<std::string, std::vector<ScanData>> no_ms3_fixtures;  // non-null + empty => record, never feed
  const int budget = 256 + 64 * static_cast<int>(ms1_scans.size() + 1);
  AcqResult a = runInterleaved(ida, ms1_scans, std::vector<ScanData>{ms2_scans[0]}, &no_ms3_fixtures,
                               budget, /*single_group_only=*/true);
  ABORT_IF(a.first_group_commands == 0)
  ABORT_IF(a.ms3_cmds.empty())

  // Every dispatched MS3 carries exactly one parent-MS2 context. Structural, not numeric: if engine drift
  // changes how many MS3s are selected, both sides move together.
  TEST_EQUAL(FLASHIdaTestAccess::ms2ContextCacheSize(*ida), a.ms3_cmds.size())
  delete ida;
}
END_SECTION

START_SECTION(ms3_exploration_variants_do_not_retain_parent_contexts)
{
  // Companion to the section above: contexts must be RETIRED as well as seeded. An MS3 exploration
  // variant is dispatched through the same initiateNextLevel path as a regular MS3, so it gets a
  // parent-MS2 context — but it returns on the EXPLORATION branch, which never reads the cache, and
  // only the regular MS3 branch erased. Every variant therefore left a permanent entry, growing for
  // the lifetime of the engine. Unobservable in the logs (a tracking id is resolved out of the pending
  // map exactly once, so an erased entry can never be re-read), which is why nothing caught it.
  //
  // Here all MS3 variants are fed back, so a correct engine ends with an empty cache.
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  std::string cfg_str = inclusionPinCytc(ms3_exploration_config);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  const int budget = 256 + 64 * static_cast<int>(ms1_scans.size() + 1);
  AcqResult a = runInterleaved(ida, ms1_scans, std::vector<ScanData>{ms2_scans[0]}, nullptr,
                               budget, /*single_group_only=*/true);
  ABORT_IF(a.ms3_cmds.empty())   // nothing to retain if nothing was dispatched

  // runInterleaved feeds every MS3 command it dequeues, and it terminates on idle>=3 (all workload
  // queues drained), so every dispatched MS3 has returned and consumed its context by now.
  TEST_EQUAL(FLASHIdaTestAccess::ms2ContextCacheSize(*ida), 0)
  delete ida;
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
  auto cmds = exploration.initiate(2, pg, 3, queue);

  // cmds[0] is the CE-0 baseline (variant_index -1); feed the FIRST REAL variant (cmds[1], variant_index 0,
  // CE 20) so the metadata under test belongs to a real CE-sweep variant, not the baseline.
  DeconvolvedSpectrum ds = makeSyntheticDeconv(1, 3);
  int tracking_id = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  ExplorationTestAccess::feedResult(exploration,tracking_id, ds, 1.0, queue);

  TEST_EQUAL(exploration.activeGroupCount(), 1)

  auto group = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group.variants[1].received, true)
  auto& stored = group.variants[1].result;
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
  // Levels 2 and 3 are no longer selected independently: applyCharacterizationMode_ projects BOTH
  // from characterization.mode, and this fixture is mode:"off", so both are None. That is the
  // point of ADR-0013 -- MS3 needs level 2 AND level 3 non-None (FLASHIda.cpp:366,
  // Exploration.cpp:728 and :730), so driving them from one enum makes "MS3 on with MS2 off"
  // unrepresentable rather than merely discouraged.
  //
  // No behaviour moved. level(2).selection is a pure MS3 gate, never a ranking metric --
  // filterAndRank early-returns for ms_level != 1 and its sole caller passes literal 1 -- so with
  // MS3 off, Intensity and None were always indistinguishable at runtime. Only the value an
  // observer reads back off the Config changed.
  TEST_EQUAL(static_cast<int>(ms2_cfg.selection), static_cast<int>(SelectionMetric::None))
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
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  // Verify first command is CE=0 baseline
  TEST_REAL_SIMILAR(cmds[0].stages[0].collision_energy, 0.0)

  // Verify exploration group has RemainingPrecursor metric
  auto group = ExplorationTestAccess::group(exploration,1);
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
    auto info = ExplorationTestAccess::feedResult(exploration,tracking_id, ds, static_cast<double>(i), queue);

    // total_variants should exclude baseline (= 5 real variants)
    TEST_EQUAL(info.group.total_variants, 5)
    TEST_REAL_SIMILAR(info.metric.remaining_ratio, -1.0)
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
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  auto group = ExplorationTestAccess::group(exploration,1);
  // Exact, not `> 0.0`: initiate() sets precursor_mz = midpoint of pg.getMzRange(charge), and the
  // synthetic group holds ONE peak at 800.0 at charge 3, so the midpoint is that peak (800+800)/2.
  // The old `> 0.0` could not fail -- getMzRange returns (-1, -10) only when the charge is outside
  // [min_abs_charge_, max_abs_charge_], which makeSyntheticPeakGroup(mz, mass, c) makes impossible
  // by constructing PeakGroup(c, c, true) and feeding the peak at that same charge.
  TEST_REAL_SIMILAR(group.precursor_mz, 800.0)

  // Feed CE=0 baseline with precursor-window signal (full intensity)
  std::vector<double> baseline_mzs = {790.0, 800.0, 810.0, 900.0};
  std::vector<double> baseline_ints = {100.0, 500.0, 200.0, 50.0};
  // In-window sum depends on isolation_width; 800.0 is precursor_mz center.

  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResult(baseline_tid, baseline_mzs.data(), baseline_ints.data(),
                         static_cast<int>(baseline_mzs.size()), 0.5, queue);

  // After baseline, group should have baseline_intensity set.
  //
  // The exact value, not `>= 0.0` -- precursorWindowIntensity_ sums non-negative intensities, so the
  // old form was a tautology that would have held just as well for the empty window this section
  // exists to rule out. At MS2 there is NO 2.0 Th floor (Exploration.cpp:165-167 applies that to
  // msn_level >= 3 only) and the synthetic group is a single peak, so mz2 - mz1 == 0 and the window is
  // the ADR-0026 decision-2 margin alone: 2 * optimal_window_margin_ == 0.8 Th, i.e. [799.6, 800.4].
  // precursorWindowIntensity_ tests `>= mz_low && <= mz_high`, so that window admits the 800.0 peak and
  // rejects 790/810/900: the sum is that peak's 500.0 alone.
  auto group_after_baseline = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group_after_baseline.has_baseline, true)
  TEST_REAL_SIMILAR(group_after_baseline.baseline_intensity, 500.0)

  // Feed second variant (CE=20) with reduced in-window signal (fragmentation removed some)
  std::vector<double> frag_mzs = {790.0, 800.0, 810.0, 900.0};
  std::vector<double> frag_ints = {50.0, 200.0, 100.0, 50.0};

  int ce20_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto ce20_info = exploration.feedResult(ce20_tid, frag_mzs.data(), frag_ints.data(),
                                          static_cast<int>(frag_mzs.size()), 1.0, queue);

  // The CE>0 variant's exact score, not a [0, 1] range check. computeRemainingPrecursorScore_ ends in
  // `score = 1.0 - deviation; if (score < 0.0) score = 0.0;` with deviation = |ratio - target| >= 0, so
  // both bounds held BY CONSTRUCTION and the pair could not fail for any input at all -- they would have
  // passed just as happily on the score-0 baseline-failure path the next section exercises.
  // Derivation: the same 0.8 Th window admits only the 800.0 peak, so the variant's
  // in-window intensity is 200.0 against the baseline's 500.0 -> ratio 0.4; remaining_precursor_target
  // is unset in this fixture and defaults to 0.1 (Config.h:174, asserted in remaining_precursor_config_parsed),
  // so deviation = 0.3 and score = 0.7.
  auto group_mid = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group_mid.variants[1].received, true)
  TEST_REAL_SIMILAR(group_mid.variants[1].score, 0.7)

  // remaining_ratio: 200.0 in-window / 500.0 baseline = 0.4. Unlike the score above, the two bounds this
  // replaces were WEAK rather than unfailable -- the ratio is a raw quotient, so -1.0 (the "N/A" sentinel
  // for no baseline or an empty one) and values above 1.0 (a variant brighter than the un-fragmented
  // reference) are both reachable and both would have failed. But [0, 1] still admits every plausible
  // depletion, so it never pinned WHICH scan the metric-independent block (Exploration.cpp:468-473)
  // ratioed against -- the whole point of that block being keyed on the baseline.
  TEST_REAL_SIMILAR(ce20_info.metric.remaining_ratio, 0.4)
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
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  // Feed baseline (CE=0) with raw data entirely OUTSIDE the precursor window
  // -> baseline_intensity = 0 -> all subsequent scores = 0 (baseline failure)
  std::vector<double> mzs = {400.0, 500.0, 600.0, 1200.0};
  std::vector<double> intensities = {100.0, 200.0, 300.0, 400.0};

  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResult(baseline_tid, mzs.data(), intensities.data(),
                         static_cast<int>(mzs.size()), 0.5, queue);

  auto group = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group.has_baseline, true)
  TEST_REAL_SIMILAR(group.baseline_intensity, 0.0)

  // Feed CE=20 variant with signal inside the precursor window
  // Since baseline = 0, score should be 0 (baseline failure path)
  std::vector<double> frag_mzs = {790.0, 800.0, 810.0, 900.0};
  std::vector<double> frag_ints = {100.0, 500.0, 200.0, 50.0};

  int ce20_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto ce20_info = exploration.feedResult(ce20_tid, frag_mzs.data(), frag_ints.data(),
                                          static_cast<int>(frag_mzs.size()), 1.0, queue);

  auto group_after = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group_after.variants[1].received, true)
  TEST_REAL_SIMILAR(group_after.variants[1].score, 0.0)
  TEST_REAL_SIMILAR(ce20_info.metric.remaining_ratio, -1.0)
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
  auto cmds = exploration.initiate(2, pg, 3, queue);
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
  auto cmds = exploration.initiate(2, pg, 3, queue);
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
  TEST_EQUAL(late_info.group.group_id, -1)                              // empty info (routing entry gone)
  TEST_EQUAL(late_info.group.variant_index, -1)                         // not routed to any variant
  TEST_EQUAL(exploration.activeGroupCount(), 0)                   // group NOT resurrected
}
END_SECTION

START_SECTION(fragment_count_requires_protein_sequence)
{
  // Config with fragment_count metric but empty protein_sequence should throw
  std::string cfg_str = std::string(exploration_config);
  // Clear protein_sequence to test empty-sequence validation path
  {
    const std::string seq = "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";
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

// The two removed SELECTION METRICS.
//
// selection_strategy.ms3.selection used to be both the MS3 on/off gate and the MS2-matcher
// chooser: intensity/qscore -> getTopFragmentMatches, terminal_fragments -> getTerminalFragmentIons,
// ambiguity_resolution -> getAmbiguityEnclosingIons. ADR-0013/0014 deleted that key. The gate is now
// characterization.mode and the matcher is hardcoded to getTopFragmentMatches -- which is what every
// MS3-enabled config selected anyway, so the other two SelectionMetric values survive in the engine
// but NO CONFIG CAN REACH THEM. That is a deliberate capability removal, recorded in ADR-0013.
//
// These sections used to assert the values PARSED. They now assert they are REJECTED, which pins the
// removal: re-introducing a config path to either matcher fails here rather than silently reviving a
// code path nothing else covers. Both surviving enum sites are checked, because either would be the
// natural place for someone to re-add them.
//
// Written as self-contained JSON rather than find/replace surgery on a shared template. The previous
// versions searched for "selection": "none" -- a key the reshape deleted -- so all four died on
// `invalid string position` instead of testing anything.
START_SECTION(terminal_fragments_is_unreachable_from_config)
{
  // Sanity: the probe itself loads, so a rejection below is the metric and not a malformed fixture.
  Config ok{metricProbe("ambiguity", "qscore")};
  TEST_EQUAL(static_cast<int>(ok.level(3).selection), static_cast<int>(SelectionMetric::Intensity))

  TEST_EXCEPTION(std::invalid_argument, Config cfg{metricProbe("terminal_fragments", "qscore")})
  TEST_EXCEPTION(std::invalid_argument, Config cfg{metricProbe("ambiguity", "terminal_fragments")})
}
END_SECTION

START_SECTION(ambiguity_resolution_is_unreachable_from_config)
{
  // Note "ambiguity_resolution" is NOT "ambiguity". The former is the deleted MATCHER; the latter is
  // a live characterization.mode value. Conflating the two is the single easiest mistake here, and it
  // is why mode hard-rejects instead of falling through the way the old objective parse did.
  TEST_EXCEPTION(std::invalid_argument, Config cfg{metricProbe("ambiguity_resolution", "qscore")})
  TEST_EXCEPTION(std::invalid_argument, Config cfg{metricProbe("ambiguity", "ambiguity_resolution")})

  Config live{metricProbe("ambiguity", "qscore")};
  TEST_EQUAL(static_cast<int>(live.characterization().objective), static_cast<int>(CharacterizationObjective::Ambiguity))
}
END_SECTION

START_SECTION(ms3_protein_sequence_only_accepted)
{
  // Config with ms3.protein_sequence should be accepted (no throw)
  Config cfg{std::string(exploration_config)};
  TEST_EQUAL(cfg.characterization().protein_sequence.empty(), false)
}
END_SECTION

// The two remaining_precursor config rejections (ADR-0026).
//
// A remaining_precursor pre-scan is scanned over its isolation window ALONE. That narrowing is what
// makes the sweep cheap, and it is sound only because the winner is always re-acquired by a
// production scan built from the un-overridden config. Two config shapes break that reasoning, and
// Config::validate() rejects both at load rather than letting a run produce plausible-looking
// emptiness:
//
//   Rejection A (Config.cpp:922-931) -- metric "remaining_precursor" with an EMPTY overrides map.
//   Rejection B (Config.cpp:946-953 and :955-962) -- metric "remaining_precursor" at a level whose
//                OWN charge mode is "multiplexed". Two pair-specific clauses, not one
//                "multiplexed anywhere" test.
//
// Three sections pin the throws; the fourth pins what must keep LOADING, because a rejection
// written one clause too broad is invisible to every positive test.
START_SECTION(remaining_precursor_without_overrides_throws)
{
  // Rejection A is LEVEL-AGNOSTIC: the reasoning is about the metric, not about which decision
  // section authored it, so validate() loops over levels_ instead of testing level 2. What makes
  // this section fail is that check being absent -- or scoped to a single level, which is the shape
  // a "fix it where it was reported" edit produces and which the second throw below exists to catch.
  //
  // The defect it stands in for is silent, not loud. An MS2 sweep with empty overrides cascades the
  // winning variant's window-narrowed spectrum into Exploration::initiateNextLevel, which finds zero
  // MS3 targets and logs `[MS3-PLAN] no_containing_fragment`. The run reads as a protein that did
  // not fragment, and nothing in it points back at the config.
  TEST_EXCEPTION(std::invalid_argument, Config(sweepProbe(", " + rpSweep(""))))
  TEST_EXCEPTION(std::invalid_argument, Config(sweepProbe("", ms3_on + rpSweep(""), true)))

  // Sanity: both probes load once the overrides map is non-empty, so the rejections above are the
  // empty map and not a malformed fixture.
  Config ms2_ok{sweepProbe(", " + rpSweep(sweep_overrides))};
  TEST_EQUAL(static_cast<int>(ms2_ok.level(2).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))

  Config ms3_ok{sweepProbe("", ms3_on + rpSweep(sweep_overrides), true)};
  TEST_EQUAL(static_cast<int>(ms3_ok.level(3).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))
}
END_SECTION

START_SECTION(remaining_precursor_multiplexed_ms2_throws)
{
  // Rejection B at level 2. [first_mass, last_mass] is ONE interval and a multiplexed readout is a
  // notch SET: the charge states of a ~12 kDa protein scatter their 2 Th windows over hundreds of
  // Th. Binding the pre-scan to the anchor window alone would isolate the whole set and read one
  // member of it; spanning the set would give back most of what the narrowing bought. Neither is
  // the scan the author asked for, so the pair is rejected instead of silently resolved one way.
  //
  // LOAD-BEARING: the probe carries a non-empty overrides map. Rejection A is checked first and
  // fires on an empty one, so without those overrides this section would throw for A's reason and
  // stay green with Rejection B deleted outright.
  TEST_EXCEPTION(std::invalid_argument,
      Config(sweepProbe(R"(, "precursor_charges": "multiplexed", )" + rpSweep(sweep_overrides))))

  // Positive control: the identical config at "single" loads, so what is rejected is the charge
  // mode and not the sweep.
  Config ok{sweepProbe(R"(, "precursor_charges": "single", )" + rpSweep(sweep_overrides))};
  TEST_EQUAL(static_cast<int>(ok.targeting().precursor_charges), static_cast<int>(ChargeAcquisitionMode::Single))
  TEST_EQUAL(static_cast<int>(ok.level(2).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))
}
END_SECTION

START_SECTION(remaining_precursor_multiplexed_ms3_throws)
{
  // Rejection B at level 3 -- a separate clause from level 2's, deliberately, because the two
  // cross-level pairs stay legal (see the negative section below). An MS3 sweep spreads its readout
  // when characterization.fragment_charges multiplexes the SUB-fragments it is reading, which is
  // the level-3 mirror of the MS2 case and is not something the level-2 clause can see.
  //
  // LOAD-BEARING, same as the MS2 section: the non-empty overrides map keeps Rejection A from
  // firing first and masking whether Rejection B exists at all.
  TEST_EXCEPTION(std::invalid_argument,
      Config(sweepProbe("", ms3_on + R"("fragment_charges": "multiplexed", )" + rpSweep(sweep_overrides), true)))

  Config ok{sweepProbe("", ms3_on + R"("fragment_charges": "single", )" + rpSweep(sweep_overrides), true)};
  TEST_EQUAL(static_cast<int>(ok.characterization().fragment_charges), static_cast<int>(ChargeAcquisitionMode::Single))
  TEST_EQUAL(static_cast<int>(ok.level(3).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))
}
END_SECTION

START_SECTION(remaining_precursor_separate_and_crosslevel_are_legal)
{
  // The NEGATIVE half of Rejection B: four shapes that must keep loading. This section is what
  // fails when the check is written one clause too broad -- "remaining_precursor plus multiplexed
  // ANYWHERE", or "remaining_precursor plus any non-single charge mode" -- neither of which the
  // three throwing sections above can tell apart from the correct pair-specific form.
  //
  // (a) and (b): `separate` FANS OUT to one scan per charge state instead of co-isolating them.
  // buildMS2's notch guard (ScanCommandQueue.cpp:314) tests `== Multiplexed` only, so a `separate`
  // scan carries no notches at all and every readout is a single contiguous interval -- exactly
  // what ADR-0026's binding describes. Rejecting `separate` would cost the mode for no reason.
  //
  // (c) and (d): the CROSS-LEVEL pairs, where the multiplexed level is not the level that sweeps.
  // Stage-0 notches change WHICH precursors are fragmented, not where a stage-1 readout sits, so an
  // MS3 sub-fragment scan under a multiplexed MS2 is still one contiguous window; and a multiplexed
  // MS3 does not touch the MS2 pre-scan above it. Only the sweeping level's OWN charge mode can
  // spread the readout that the sweep is bound to.
  {
    // (a) level-2 sweep, precursor_charges "separate"
    Config cfg{sweepProbe(R"(, "precursor_charges": "separate", )" + rpSweep(sweep_overrides))};
    TEST_EQUAL(static_cast<int>(cfg.level(2).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))
    TEST_EQUAL(static_cast<int>(cfg.targeting().precursor_charges), static_cast<int>(ChargeAcquisitionMode::Separate))
  }
  {
    // (b) level-3 sweep, fragment_charges "separate"
    Config cfg{sweepProbe("", ms3_on + R"("fragment_charges": "separate", )" + rpSweep(sweep_overrides), true)};
    TEST_EQUAL(static_cast<int>(cfg.level(3).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))
    TEST_EQUAL(static_cast<int>(cfg.characterization().fragment_charges), static_cast<int>(ChargeAcquisitionMode::Separate))
  }
  {
    // (c) level-3 sweep under a MULTIPLEXED level 2. Level 2 does not sweep, so nothing binds its
    // scan range; the MS3 readout that IS bound sits a stage below the notches.
    Config cfg{sweepProbe(R"(, "precursor_charges": "multiplexed")", ms3_on + rpSweep(sweep_overrides), true)};
    TEST_EQUAL(static_cast<int>(cfg.level(3).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))
    TEST_EQUAL(static_cast<int>(cfg.targeting().precursor_charges), static_cast<int>(ChargeAcquisitionMode::Multiplexed))
    TEST_EQUAL(static_cast<int>(cfg.level(2).exploration), static_cast<int>(ExplorationMetric::None))
  }
  {
    // (d) the mirror: level-2 sweep under MULTIPLEXED fragment charges. mode is on, so level 3 is
    // live and its charge mode is real -- it simply is not the level being swept.
    Config cfg{sweepProbe(", " + rpSweep(sweep_overrides), ms3_on + R"("fragment_charges": "multiplexed")", true)};
    TEST_EQUAL(static_cast<int>(cfg.level(2).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))
    TEST_EQUAL(static_cast<int>(cfg.characterization().fragment_charges), static_cast<int>(ChargeAcquisitionMode::Multiplexed))
    TEST_EQUAL(static_cast<int>(cfg.level(3).exploration), static_cast<int>(ExplorationMetric::None))
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// The ADR-0026 BINDING itself (Exploration.cpp:254-261), where the four sections above pin only the
// config shapes that make it representable. A RemainingPrecursor pre-scan's first_mass/last_mass are
// set to precursor_mz -/+ isolation_width/2 -- at MS2 and MS3 alike, and for EVERY variant including
// the CE-0 baseline, because a full-range trap fill and a narrow one are not comparable denominators
// for the ratio. Five sections: the two levels, the production scan the binding must NOT reach, the
// explicit-range escape hatch, and the ADR-0023 FORCED metric that no config declares.
//
// EVERY RANGE ASSERTION BELOW IS RELATIONAL -- `group.precursor_mz -/+ group.isolation_width / 2.0`
// read back through ExplorationTestAccess, never a literal, and that is load-bearing rather than
// stylistic. ADR-0026 decision 2 ships in THIS push: initiate() no longer takes isolation_width as the
// bare PeakGroup m/z span, because buildMS2 pads the commanded window by optimal_window_margin_ on
// each side (ScanCommandQueue.cpp:293-295, .4 Th => 0.8 Th total) and the interval the metric SUMS
// must be the interval the instrument was told to ISOLATE. So MS2 is now (mz2 - mz1) + 0.8 while MS3
// keeps the 2.0 Th floor, the two arms of Exploration.cpp:165-167. Every window in these five sections
// moved when that landed and not one line here needed touching, because the group and the command are
// read from the same computation; hard-coded 499.0/501.0 bounds would have had to be re-derived in
// five places, which is how a golden-shaped assertion turns into a reason not to make the change.
//
// The only literals below are the two MS3 isolation_width pins, and they are there because the
// relational pair cannot fail on its own once std::max(..., 2.0) is in play.
/////////////////////////////////////////////////////////////

START_SECTION(remaining_precursor_binds_scan_range_ms2)
{
  // FAILS WHEN the binding is absent at level 2: every command then carries ms_settings.ms2's
  // 300/1800 and the sweep pays a 1500 Th fill to read the 2 Th it scores from.
  //
  // ALSO FAILS under the suppression test written against base_config instead of cfg.overrides --
  // the shape Exploration.cpp:244-249 warns about. applyOverrides has already run at :128 and an
  // authored ms_settings range lands in the SAME ScanConfig field with the SAME 0 default
  // (Config.cpp:130-131, Config.h:124), so `base_config.first_mass != 0` reads true here purely
  // because ms2 names a range, and the binding would be skipped for a config that asked for no
  // range override at all. THE NON-ZERO 300/1800 IS WHAT MAKES THAT DETECTABLE: leave ms2 at the 0
  // default and both forms agree, and this section passes under either.
  const std::string cfg_json = R"JSON({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": {
      "rank_by": "qscore", "max_precursors": 3,
      "exploration": { "metric": "remaining_precursor", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0,
                       "overrides": { "analyzer": "IonTrap" } }
    },
    "characterization": { "mode": "off" },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29,
               "first_mass": 300, "last_mass": 1800 }
    }
  })JSON";

  Config cfg{cfg_json};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  // The control the loop below is read against: what an UNBOUND command carries. Asserted off the
  // parsed config rather than restated as a number, so a fixture typo cannot make the section pass.
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].first_mass, 300.0)
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].last_mass, 1800.0)

  // A SECOND peak at the same charge, so the group's m/z span is non-degenerate. The shared helper
  // makes ONE peak and getMzRange then returns (mz, mz); at MS2 there is no 2.0 Th floor
  // (Exploration.cpp:165-167 applies that to msn_level >= 3 only), so a one-peak group would bind the
  // decision-2 margin alone -- precursor_mz -/+ 0.4, symmetric about the centre, which lo and hi would
  // then be two readings of. Two peaks make them independent claims about the measured SPAN.
  PeakGroup pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  FLASHHelperClasses::LogMzPeak lp2;
  lp2.mz = 802.0;
  lp2.abs_charge = 3;
  pg.push_back(lp2);

  std::vector<ScanCommand> cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 4)   // CE-0 baseline + CE 20 / 25 / 30
  ABORT_IF(cmds.size() != 4)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(group.msn_level, 2)
  TEST_EQUAL(static_cast<int>(group.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))

  // cmds[0] IS THE CE-0 BASELINE, asserted rather than assumed. It is variant_params[0] (inserted at
  // Exploration.cpp:138) and the binding reaches it only by being written onto base_config ABOVE the
  // loop, so that the baseline takes the same unconditional `variant_config = base_config` at :273 as
  // every other variant. A binding moved inside an `is_baseline == false` branch would leave the CE-0
  // reference scanning 300-1800 while its variants scan ~2 Th, and every ratio in the group would then
  // divide a narrow fill by a full-range one.
  TEST_EQUAL(group.variants[0].is_baseline, true)
  TEST_EQUAL(group.variants[0].variant_index, -1)
  TEST_REAL_SIMILAR(group.variants[0].collision_energy, 0.0)
  TEST_EQUAL(group.variants[0].cmd.scan_id, cmds[0].scan_id)

  const double lo = group.precursor_mz - group.isolation_width / 2.0;
  const double hi = group.precursor_mz + group.isolation_width / 2.0;
  TEST_EQUAL(hi > lo, true)   // non-degenerate, so the pair below cannot be satisfied by one value

  for (const ScanCommand& c : cmds)
  {
    TEST_EQUAL(c.msn_level, 2)
    TEST_REAL_SIMILAR(c.first_mass, lo)
    TEST_REAL_SIMILAR(c.last_mass, hi)
  }
}
END_SECTION

START_SECTION(remaining_precursor_binds_scan_range_ms3)
{
  // FAILS WHEN the binding is absent at level 3 -- or present but written on the MS2 branch only.
  // The binding sits ABOVE the variant loop and before the buildMS2/buildMS3 fork
  // (Exploration.cpp:254-261 vs :280-290), which is exactly why one edit covers both levels; a
  // version placed inside that fork is the plausible way to get MS2 right and MS3 wrong, and the
  // MS2 section alone cannot see it.
  //
  // ms_settings.ms3 names 350/1750 for the same two reasons ms2 named 300/1800 above: it is the
  // control for "unbound", and being non-zero it is also what makes a base_config-keyed suppression
  // test fail here rather than silently agree.
  const std::string cfg_json = R"JSON({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": { "rank_by": "qscore", "max_precursors": 3 },
    "characterization": {
      "mode": "ambiguity", "protein_sequence": "PEPTIDER", "max_targets": 3,
      "exploration": { "metric": "remaining_precursor", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0,
                       "overrides": { "analyzer": "IonTrap" } }
    },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 },
      "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25,
               "first_mass": 350, "last_mass": 1750 }
    }
  })JSON";

  Config cfg{cfg_json};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  TEST_REAL_SIMILAR(cfg.level(3).scans[0].first_mass, 350.0)
  TEST_REAL_SIMILAR(cfg.level(3).scans[0].last_mass, 1750.0)

  PeakGroup fragment_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  ScanCommand ms2_ctx = queue.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3,
                                       cfg.level(2).scans[0], 2, 0);

  std::vector<ScanCommand> cmds = exploration.initiate(3, fragment_pg, 2, queue, nullptr, &ms2_ctx);
  TEST_EQUAL(static_cast<int>(cmds.size()), 4)   // CE-0 baseline + CE 20 / 25 / 30
  ABORT_IF(cmds.size() != 4)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(group.msn_level, 3)
  TEST_EQUAL(static_cast<int>(group.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))

  // Same baseline pin as the MS2 section: the CE-0 scan is in the loop below, not exempt from it.
  TEST_EQUAL(group.variants[0].is_baseline, true)
  TEST_EQUAL(group.variants[0].variant_index, -1)
  TEST_EQUAL(group.variants[0].cmd.scan_id, cmds[0].scan_id)

  const double lo = group.precursor_mz - group.isolation_width / 2.0;
  const double hi = group.precursor_mz + group.isolation_width / 2.0;
  // Exact, not `hi > lo`: at MS3 isolation_width is std::max(mz2 - mz1, 2.0) (Exploration.cpp:165-167),
  // so hi - lo >= 2.0 for EVERY input the floor can see -- including the (-1, -10) sentinel getMzRange
  // returns for a charge with no peaks. A comparison against zero was documentation, not a test. The
  // fixture is a single peak at 500.0, so mz2 - mz1 == 0 and the floor IS the window: this fails if the
  // floor is dropped, re-tuned, or if the MS2 margin branch is taken at level 3.
  TEST_REAL_SIMILAR(group.isolation_width, 2.0)

  // The window bound is the FRAGMENT stage's, not the inherited MS2 precursor's. An MS3 command
  // carries both (stage[0] ~ 800 Th from ms2_ctx, stage[1] = the ~500 Th fragment), and
  // precursorWindowIntensity_ sums around group.precursor_mz -- so a binding that reached for
  // stages[0] would command one window and score another.
  TEST_REAL_SIMILAR(cmds[0].stages[1].precursor_mz, group.precursor_mz)
  TEST_EQUAL(std::abs(cmds[0].stages[0].precursor_mz - group.precursor_mz) > 100.0, true)

  for (const ScanCommand& c : cmds)
  {
    TEST_EQUAL(c.msn_level, 3)
    TEST_EQUAL(c.num_stages, 2)
    TEST_REAL_SIMILAR(c.first_mass, lo)
    TEST_REAL_SIMILAR(c.last_mass, hi)
  }
}
END_SECTION

START_SECTION(production_scan_keeps_configured_range)
{
  // FAILS WHEN the binding leaks past the pre-scans into the post-winner close-out.
  //
  // The production scan is the ONE acquisition of the group that is meant to be identified, and it
  // is rebuilt from level_config.scans[0] (Exploration.cpp:750) rather than from base_config, which
  // is why it keeps its configured 200-2000 while the sweep that chose its CE ran at ~2 Th. That
  // separation is structural -- no flag defends it -- so an edit that "hoisted" the narrowing onto
  // level_config, or pushed it down into buildMS3 where every MS3 command would inherit it, would
  // narrow the close-out too. Nothing would throw: ADR-0026's entire safety argument (decision 3,
  // the winner is always re-acquired at full range) would simply stop being true, and the run would
  // look like a working narrow sweep that identified nothing.
  //
  // The pre-scan assertion below is the positive control. It proves the binding DID fire in this
  // very group, so the production scan's 200/2000 is a contrast rather than a binding that never ran.
  const std::string cfg_json = R"JSON({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": { "rank_by": "qscore", "max_precursors": 3 },
    "characterization": {
      "mode": "ambiguity", "protein_sequence": "PEPTIDER", "max_targets": 3,
      "exploration": { "metric": "remaining_precursor", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0,
                       "overrides": { "analyzer": "IonTrap" } }
    },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 },
      "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25,
               "first_mass": 200, "last_mass": 2000 }
    }
  })JSON";

  Config cfg{cfg_json};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  PeakGroup fragment_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  ScanCommand ms2_ctx = queue.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3,
                                       cfg.level(2).scans[0], 2, 0);

  std::vector<ScanCommand> cmds = exploration.initiate(3, fragment_pg, 2, queue, nullptr, &ms2_ctx);
  TEST_EQUAL(static_cast<int>(cmds.size()), 4)
  ABORT_IF(cmds.size() != 4)

  // Read the group BEFORE the sweep completes -- completion erases it from active_groups_.
  auto group = ExplorationTestAccess::group(exploration, 1);
  const double lo = group.precursor_mz - group.isolation_width / 2.0;
  const double hi = group.precursor_mz + group.isolation_width / 2.0;
  TEST_REAL_SIMILAR(cmds[0].first_mass, lo)   // positive control: the pre-scans ARE narrowed here
  TEST_REAL_SIMILAR(cmds[3].last_mass, hi)

  // Drive the sweep to completion with RAW ARRAYS. RemainingPrecursor scores from isolation-window
  // intensity alone, so one peak at the centre and two 50 Th away make the in-window sum exactly
  // controllable; the ExplorationTestAccess deconvolution bypass the mass_count sections use cannot
  // drive this metric at all. A baseline with no in-window signal aborts the group on
  // baseline_failed (set at Exploration.cpp:434-437), which suppresses winner selection at :664 and
  // takes the no-follow-up arm at :683, so no production scan is emitted -- the assertion would then
  // read as a missing command rather than a narrowed one.
  const double mz_c = group.precursor_mz;
  std::vector<double> mzs{mz_c - 50.0, mz_c, mz_c + 50.0};
  // score = 1 - |ratio - remaining_precursor_target|, target = 0.1 (the level default), so the
  // ratios 0.9 / 0.4 / 0.1 give 0.2 / 0.7 / 1.0 and the LAST-fed variant wins. A tie-break-order
  // winner would land on cmds[1] and the CE assertion below would catch it.
  std::vector<double> baseline_ints{300.0, 1000.0, 700.0};
  std::vector<double> ce20_ints{300.0, 900.0, 700.0};
  std::vector<double> ce25_ints{300.0, 400.0, 700.0};
  std::vector<double> ce30_ints{300.0, 100.0, 700.0};

  exploration.feedResult(cmds[0].scan_id, mzs.data(), baseline_ints.data(),
                         static_cast<int>(mzs.size()), 0.5, queue);
  exploration.feedResult(cmds[1].scan_id, mzs.data(), ce20_ints.data(),
                         static_cast<int>(mzs.size()), 1.0, queue);
  exploration.feedResult(cmds[2].scan_id, mzs.data(), ce25_ints.data(),
                         static_cast<int>(mzs.size()), 1.5, queue);
  auto info = exploration.feedResult(cmds[3].scan_id, mzs.data(), ce30_ints.data(),
                                     static_cast<int>(mzs.size()), 2.0, queue);

  TEST_EQUAL(exploration.activeGroupCount(), 0)
  TEST_REAL_SIMILAR(info.metric.score, 1.0)
  TEST_EQUAL(info.group.winner_tracking_id, ScanCommandQueue::encode(cmds[3].scan_id))

  // Exactly one close-out command. Both ADR-0020 gates hold here (non-empty overrides AND a
  // measuring MS3 metric), which is not a weakness of this section: WHICH gate fired is
  // ms3_measuring_metric_always_reacquires_without_overrides' subject, and this one only needs a
  // production scan to exist so that its range can be read.
  TEST_EQUAL(static_cast<int>(info.commands.size()), 1)
  // Release build: no _ITERATOR_DEBUG_LEVEL and no OPENMS_PRECONDITION, so indexing an empty
  // `commands` is a read near address 0 that kills the binary and takes every later section with it.
  ABORT_IF(info.commands.size() != 1)
  TEST_EQUAL(info.commands[0].msn_level, 3)
  TEST_REAL_SIMILAR(info.commands[0].stages[1].collision_energy, 30.0)  // rebuilt at the WINNING CE

  // THE ASSERTIONS THIS SECTION EXISTS FOR: the configured MS3 range, not the ~2 Th window its own
  // pre-scans were bound to. Stated twice on purpose -- against the literals so the intent is
  // readable, and against the parsed config so the pair cannot both drift to some third value.
  TEST_REAL_SIMILAR(info.commands[0].first_mass, 200.0)
  TEST_REAL_SIMILAR(info.commands[0].last_mass, 2000.0)
  TEST_REAL_SIMILAR(info.commands[0].first_mass, cfg.level(3).scans[0].first_mass)
  TEST_REAL_SIMILAR(info.commands[0].last_mass, cfg.level(3).scans[0].last_mass)
  TEST_EQUAL(info.commands[0].first_mass < lo && info.commands[0].last_mass > hi, true)
}
END_SECTION

START_SECTION(explicit_range_override_suppresses_binding)
{
  // FAILS WHEN the escape hatch (ADR-0026 decision 5) is missing: the binding would overwrite the
  // author's 600/1600 with the ~2 Th window, and the one knob that lets someone tune a sweep against
  // real hardware without a rebuild would be inert -- quietly, because the resulting scans still look
  // like a working narrow sweep and only scan_commands.tsv's range columns (decision 6) can tell the
  // two apart.
  //
  // WHAT THIS SECTION CANNOT DO, AND WHERE THAT IS COVERED. It does NOT separate a suppression test
  // written against cfg.overrides from one written against base_config: applyOverrides has already
  // folded 600/1600 into base_config by the time the binding is reached (Exploration.cpp:128), so
  // both forms suppress here and both emit 600/1600. The discrimination lives in the two binding
  // sections above, whose ms_settings ranges are non-zero and must STILL be narrowed -- a
  // base_config test reads those as "a range override is present" and skips the binding. The three
  // sections only cover ADR-0026 decision 5 together.
  //
  // ms_settings.ms2 carries its OWN 100/2000, deliberately different from the override, so the
  // 600/1600 below is attributable to the overrides map: set ms2 to 600/1600 and the assertions
  // would pass even if applyOverrides never ran at all.
  //
  // The override VALUES ARE JSON STRINGS. Config.cpp:741 reads the map with `.get<std::string>()`,
  // so a bare `"first_mass": 600` throws a nlohmann type error out of the constructor and the
  // section fails as an unloadable fixture rather than as a suppression defect.
  const std::string cfg_json = R"JSON({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": {
      "rank_by": "qscore", "max_precursors": 3,
      "exploration": { "metric": "remaining_precursor", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0,
                       "overrides": { "analyzer": "Orbitrap", "first_mass": "600", "last_mass": "1600" } }
    },
    "characterization": { "mode": "off" },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29,
               "first_mass": 100, "last_mass": 2000 }
    }
  })JSON";

  Config cfg{cfg_json};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  // The map survived parsing under both keys -- the escape hatch's actual input. tolerance_ppm was
  // erased from this map once; that is now only history in a comment (Config.cpp:721-724) and the live
  // code throws on the key instead (Config.cpp:736-740), so "an exploration key reached overrides" is
  // not something to take on trust.
  TEST_EQUAL(cfg.level(2).overrides.count("first_mass"), 1u)
  TEST_EQUAL(cfg.level(2).overrides.count("last_mass"), 1u)
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].first_mass, 100.0)   // the AUTHORED scan range, untouched
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].last_mass, 2000.0)

  PeakGroup pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  std::vector<ScanCommand> cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 4)
  ABORT_IF(cmds.size() != 4)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(static_cast<int>(group.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))

  // What the binding WOULD have written, computed rather than asserted, so the section states the
  // collision it is protecting against instead of merely avoiding it. ~800 Th against 600 -- and it
  // stays ~800 after the decision-2 margin push, so the 1.0 Th slack is not a tolerance in disguise.
  const double bound_lo = group.precursor_mz - group.isolation_width / 2.0;
  TEST_EQUAL(std::abs(bound_lo - 600.0) > 1.0, true)

  for (const ScanCommand& c : cmds)
  {
    TEST_REAL_SIMILAR(c.first_mass, 600.0)
    TEST_REAL_SIMILAR(c.last_mass, 1600.0)
  }
}
END_SECTION

START_SECTION(forced_remaining_precursor_binds_scan_range)
{
  // FAILS WHEN the binding is gated on cfg.exploration rather than on group.exploration_metric.
  // The two agree everywhere except here: ADR-0023 decision 11 DRAGS an exhaustive-mode unassigned
  // mass -- ion class 'u', which MS3FragmentMatcher::isKnownIonClass rejects -- onto
  // RemainingPrecursor whatever the config asked for (Exploration.cpp:196-200). This config asks for
  // fragment_count, so a cfg-keyed binding emits ms_settings.ms3's 200/2000 and the forced sweep
  // pays a full-range fill for every pre-scan it then reads ~2 Th of. Mirrors
  // ProteoformTracker_Exhaustive_test.cpp:682-719, which pins the force itself; this section pins
  // what the force implies for the scan range.
  //
  // WHY THE FORCED PATH NEEDS NO CONFIG GUARANTEE, and why no third Config::validate() rejection
  // should be added for it. Narrowing is safe because the winner is always re-acquired at full
  // range; for a CONFIGURED sweep that is guaranteed by Rejection A (ADR-0026 decision 3), and here
  // it is structural instead:
  //     force fires => msn_level >= 3 AND metric == RemainingPrecursor (measuring by definition)
  //                 => measuring_ms3_sweep (Exploration.cpp:747-748)
  //                 => ADR-0020 gate #2 fires => a full-range production scan re-acquires.
  // The force's own precondition IS gate #2's condition. Nor is there a cascade to strand: MS3 is
  // terminal at the `< 3` MS4 wall (:816), so "a window-only spectrum yields zero next-level
  // targets" -- the hazard Rejection A exists for -- is an MS2-only failure.
  //
  // The overrides map is non-empty and deliberately carries NO range key. Unlike the four sections
  // above it is NOT there to keep Rejection A from firing first: that rejection keys on the
  // CONFIGURED metric, which is fragment_count here, so it cannot fire and the forced metric never
  // reaches validate() at all -- which is the asymmetry this section documents. It is there so the
  // narrowing below is attributable to group.exploration_metric alone, with the escape hatch of
  // ADR-0026 decision 5 provably not engaged.
  const std::string cfg_json = R"JSON({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": { "rank_by": "qscore", "max_precursors": 3 },
    "characterization": {
      "mode": "exhaustive", "protein_sequence": "PEPTIDEK", "max_targets": 3,
      "exploration": { "metric": "fragment_count", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0,
                       "overrides": { "analyzer": "IonTrap" } }
    },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 },
      "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25,
               "first_mass": 200, "last_mass": 2000 }
    }
  })JSON";

  Config cfg{cfg_json};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  // The CONFIG says reading. The GROUP will say measuring. That gap is the whole section, so both
  // halves are asserted rather than one inferred from the other.
  TEST_EQUAL(static_cast<int>(cfg.level(3).exploration), static_cast<int>(ExplorationMetric::FragmentCount))
  TEST_EQUAL(cfg.level(3).overrides.count("first_mass"), 0u)
  TEST_EQUAL(cfg.level(3).overrides.count("last_mass"), 0u)
  TEST_REAL_SIMILAR(cfg.level(3).scans[0].first_mass, 200.0)   // what a cfg-keyed binding would emit
  TEST_REAL_SIMILAR(cfg.level(3).scans[0].last_mass, 2000.0)

  PeakGroup unassigned_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  ScanCommand ms2_ctx = queue.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3,
                                       cfg.level(2).scans[0], 2, 0);

  // 'u'/0 is exactly what planExhaustive_ labels a mass that matched no theoretical fragment.
  std::vector<ScanCommand> cmds = exploration.initiate(3, unassigned_pg, 2, queue, nullptr, &ms2_ctx, 'u', 0);
  TEST_EQUAL(static_cast<int>(cmds.size()), 4)   // CE-0 baseline + CE 20 / 25 / 30
  ABORT_IF(cmds.size() != 4)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(group.msn_level, 3)
  TEST_EQUAL(static_cast<int>(group.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))
  TEST_EQUAL(group.variants[0].is_baseline, true)
  TEST_EQUAL(group.variants[0].cmd.scan_id, cmds[0].scan_id)

  const double lo = group.precursor_mz - group.isolation_width / 2.0;
  const double hi = group.precursor_mz + group.isolation_width / 2.0;
  // The FORCED sweep takes the same MS3 floor a configured one does (Exploration.cpp:165-167), so pin
  // the width rather than `hi > lo`, which std::max(..., 2.0) makes true for every possible input --
  // the (-1, -10) sentinel getMzRange returns for a charge with no peaks included. 2.0 also says the
  // width was computed from the FRAGMENT group and not from ms2_ctx, whose buildMS2 window is 0.8 Th.
  TEST_REAL_SIMILAR(group.isolation_width, 2.0)
  TEST_EQUAL(lo > cfg.level(3).scans[0].first_mass, true)   // strictly inside the configured range,
  TEST_EQUAL(hi < cfg.level(3).scans[0].last_mass, true)    // so "bound" and "unbound" cannot coincide

  for (const ScanCommand& c : cmds)
  {
    TEST_EQUAL(c.msn_level, 3)
    TEST_REAL_SIMILAR(c.first_mass, lo)
    TEST_REAL_SIMILAR(c.last_mass, hi)
  }
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
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  auto group = ExplorationTestAccess::group(exploration,1);
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
  double score_perfect = info_perfect.metric.score;
  double ratio_perfect = info_perfect.metric.remaining_ratio;

  std::vector<double> over_ints = {500.0};
  int ce25_tid = queue.decode(std::string(cmds[2].scan_description).substr(0, 3));
  auto info_over = exploration.feedResult(ce25_tid, baseline_mzs.data(), over_ints.data(),
                                           1, 2.0, queue);
  double score_over = info_over.metric.score;
  double ratio_over = info_over.metric.remaining_ratio;

  // (The former `score_perfect > score_over` ordering check is gone: the two exact values below
  // already imply it, and an ordering assertion that survives a shift of BOTH scores is the weaker
  // of the two claims, never an independent one.)
  TEST_REAL_SIMILAR(score_perfect, 1.0)
  TEST_REAL_SIMILAR(score_over, 0.6)
  TEST_REAL_SIMILAR(ratio_perfect, 0.1)
  TEST_REAL_SIMILAR(ratio_over, 0.5)

  (void)iso_half;
}
END_SECTION

START_SECTION(fragment_match_propagated_in_feed_result)
{
  // CA2 real MS2 (HCD25, scan 181, ~20 fragments) against the CA2 sequence: a confident full-length ID.
  // Build a CA-sequence config from the shared exploration_config WITHOUT editing the shared literal
  // (mirror the ms3_winner_cfg.replace pattern): swap the cytC protein_sequence for CA2.
  std::string ca_cfg(exploration_config);
  const std::string cytc_seq = "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";
  ca_cfg.replace(ca_cfg.find(cytc_seq), cytc_seq.size(), std::string(ca_sequence));
  Config cfg{ca_cfg};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  TEST_EQUAL(cfg.characterization().protein_sequence.empty(), false)

  auto pg = makeSyntheticPeakGroup(968.4916, 29006.3881, 30);   // INTACT CA precursor z30
  auto cmds = exploration.initiate(2, pg, 30, queue);

  auto ms2_scans = loadTsvScans(ms2_ca_path);
  ABORT_IF(ms2_scans.empty())
  const auto& ms2_data = ms2_scans[0];

  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info = exploration.feedResult(tracking_id,
      ms2_data.mzs.data(), ms2_data.ints.data(),
      static_cast<int>(ms2_data.mzs.size()), ms2_data.rt, queue);

  // Real CA2 spectrum should produce a confident full-length fragment match
  TEST_EQUAL(info.metric.fragment_count > 0, true)
  TEST_EQUAL(info.identification.matched_protein.empty(), false)
  TEST_STRING_EQUAL(info.identification.proteoform_sequence, std::string(ca_sequence))
  // Full-length: FLASHExtender leaves region at -1 ("full sequence") or sets 0/259; assert NOT truncated
  // (rejects an N-terminal start>0 or a C-terminal end in [0,259)).
  TEST_EQUAL(info.identification.result.region_start <= 0, true)
  TEST_EQUAL(info.identification.result.region_end < 0 || info.identification.result.region_end == 259, true)
}
END_SECTION

START_SECTION(full_length_ca_fragment_match_high_count)
{
  // CA2 real MS2 at higher collision energy (HCD45, scan 185, ~36 fragments) against the CA2 sequence:
  // a confident full-length ID with a high fragment count. Same shared-config CA swap as above.
  const std::string ms2_ca_hcd45_path = "../../FlashIDA/test-data/spectra/ms2_ca_hcd45_scan185.txt";

  std::string ca_cfg(exploration_config);
  const std::string cytc_seq = "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";
  ca_cfg.replace(ca_cfg.find(cytc_seq), cytc_seq.size(), std::string(ca_sequence));
  Config cfg{ca_cfg};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  TEST_EQUAL(cfg.characterization().protein_sequence.empty(), false)

  auto pg = makeSyntheticPeakGroup(968.4916, 29006.3881, 30);   // INTACT CA precursor z30
  auto cmds = exploration.initiate(2, pg, 30, queue);

  auto ms2_scans = loadTsvScans(ms2_ca_hcd45_path);
  ABORT_IF(ms2_scans.empty())
  const auto& ms2_data = ms2_scans[0];

  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info = exploration.feedResult(tracking_id,
      ms2_data.mzs.data(), ms2_data.ints.data(),
      static_cast<int>(ms2_data.mzs.size()), ms2_data.rt, queue);

  // High-energy CA2 spectrum should produce a high-count full-length fragment match
  TEST_EQUAL(info.metric.fragment_count >= 20, true)
  TEST_EQUAL(info.identification.matched_protein.empty(), false)
  TEST_STRING_EQUAL(info.identification.proteoform_sequence, std::string(ca_sequence))
  // Full-length: FLASHExtender leaves region at -1 ("full sequence") or sets 0/259; assert NOT truncated
  // (rejects an N-terminal start>0 or a C-terminal end in [0,259)).
  TEST_EQUAL(info.identification.result.region_start <= 0, true)
  TEST_EQUAL(info.identification.result.region_end < 0 || info.identification.result.region_end == 259, true)
}
END_SECTION

START_SECTION(fragment_count_zero_without_protein_sequence)
{
  Config cfg{std::string(remaining_precursor_config)};
  TEST_EQUAL(cfg.characterization().protein_sequence.empty(), true)

  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  DeconvolvedSpectrum baseline_ds = makeSyntheticDeconv(0, 1);
  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  ExplorationTestAccess::feedResult(exploration,baseline_tid, baseline_ds, 0.0, queue);

  DeconvolvedSpectrum ds = makeSyntheticDeconv(1, 5);
  int tracking_id = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto info = ExplorationTestAccess::feedResult(exploration,tracking_id, ds, 1.0, queue);

  TEST_EQUAL(info.metric.fragment_count, 0)
  TEST_EQUAL(info.identification.matched_protein.empty(), true)
  TEST_EQUAL(info.identification.proteoform_sequence.empty(), true)
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
                                        cfg.level(2).scans[0], 2, 0);

  auto cmds = exploration.initiate(3, fragment_pg, 2, queue, nullptr, &ms2_ctx);
  // RemainingPrecursor: 1 baseline + 5 CE variants (20,25,30,35,40) = 6
  TEST_EQUAL(static_cast<int>(cmds.size()), 6)

  auto group = ExplorationTestAccess::group(exploration,1);
  // isolation_width should be floored to 2.0 (not 0.0)
  TEST_REAL_SIMILAR(group.isolation_width, 2.0)

  double mz_center = group.precursor_mz;

  // Feed baseline (CE=0) with signal at precursor center
  std::vector<double> baseline_mzs = {mz_center};
  std::vector<double> baseline_ints = {1000.0};
  int baseline_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  exploration.feedResult(baseline_tid, baseline_mzs.data(), baseline_ints.data(),
                         static_cast<int>(baseline_mzs.size()), 0.5, queue);

  auto group_after = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group_after.has_baseline, true)
  // With 2.0 Da window [499.0, 501.0], mz_center=500.0 is in-window
  TEST_REAL_SIMILAR(group_after.baseline_intensity, 1000.0)

  // Feed CE=20 variant with 100.0 intensity -> ratio = 0.1
  std::vector<double> variant_ints = {100.0};
  int ce20_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto info = exploration.feedResult(ce20_tid, baseline_mzs.data(), variant_ints.data(),
                                      1, 1.0, queue);

  // Score should be real (not -1.0), ratio = 100/1000 = 0.1
  TEST_REAL_SIMILAR(info.metric.remaining_ratio, 0.1)
  // target=0.1, deviation=0.0, score=1.0. (The former trailing `score > 0.0` is gone -- it was implied
  // by this line and unfailable on its own: computeRemainingPrecursorScore_ floors at 0 and caps at 1,
  // so it could only have caught the exact-zero case this exact assertion already excludes.)
  TEST_REAL_SIMILAR(info.metric.score, 1.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0026 DECISION 2: THE MARGIN CORRECTION. The interval a metric SUMS must be the interval the
// instrument was told to ISOLATE.
//
// buildMS2 has always padded the measured [mz1, mz2] by optimal_window_margin_ per side before
// commanding it (ScanCommandQueue.cpp:291-295), while Exploration::initiate took the BARE span. So a
// RemainingPrecursor sweep scored a window 0.8 Th narrower than the one it acquired: every peak in
// those two flanks was paid for, transmitted, and then dropped on the floor by
// precursorWindowIntensity_. Nothing read wrong -- the ratio was simply computed over a different
// interval than the one the fill was spent on.
//
// The correction is Exploration.cpp:165-167's MS2 arm, `(mz2 - mz1) + 2.0 * optimal_window_margin_`,
// and the three sections below are the whole case for shipping it in the SAME push as the ADR-0026
// binding rather than a later one:
//   * ms2_scoring_window_matches_commanded_isolation -- the ONLY section in this file that goes red
//     without it, because it is the only one whose fixture puts signal in the 0.4 Th flanks.
//   * optimal_window_margin_has_one_definition -- the drift guard for the constant the correction
//     had to be collapsed onto (three TU-local copies -> NotchSelection.h:123).
//   * single_isotope_charge_yields_non_degenerate_scan_range -- the defect that makes "later"
//     impossible: without the margin the binding emits a zero-width scan range.
/////////////////////////////////////////////////////////////

START_SECTION(ms2_scoring_window_matches_commanded_isolation)
{
  // FAILS WHEN Exploration.cpp:167 loses `+ 2.0 * optimal_window_margin_` and reverts to the bare
  // `(mz2 - mz1)`. That is the ONE edit this file can detect numerically, and only from here: every
  // other remaining_precursor section feeds peaks that are 10 Th or more off centre (790/810/900) or
  // sits at MS3, where std::max(..., 2.0) swallows a 0.8 Th difference whole. This fixture is built
  // the other way round -- two of its four baseline peaks land INSIDE the commanded window and
  // OUTSIDE the un-margined one, so the two definitions of "the window" give two different sums.
  //
  // THE ARITHMETIC, derived from the fixture and from precursorWindowIntensity_'s actual predicate
  // (`mzs[i] >= mz_low && mzs[i] <= mz_high`, Exploration.cpp:1227), not guessed:
  //   getMzRange(3) = (800.0, 802.0)          two peaks at one charge -> a NON-degenerate span
  //   precursor_mz  = 801.0                   (mz1 + mz2) / 2, computed before any margin
  //   raw_half      = (802.0 - 800.0) / 2 = 1.0
  //   OLD window    = [800.0, 802.0]          precursor_mz +/- raw_half
  //   NEW window    = [799.6, 802.4]          precursor_mz +/- (raw_half + optimal_window_margin_)
  // The two flank peaks sit at precursor_mz +/- (raw_half + 0.2), i.e. 0.2 Th inside the new edge and
  // 0.2 Th outside the old one -- far enough from both that the `>=`/`<=` boundary cannot decide the
  // result, which a peak placed exactly on an edge would.
  //   NEW sum = 300 + 1000 + 200 = 1500.0     (810.0's 700.0 is outside BOTH, so this is not merely
  //   OLD sum =       1000.0                   "the sum of everything fed")
  const std::string cfg_json = R"JSON({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": {
      "rank_by": "qscore", "max_precursors": 3,
      "exploration": { "metric": "remaining_precursor", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0,
                       "overrides": { "analyzer": "IonTrap" } }
    },
    "characterization": { "mode": "off" },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 }
    }
  })JSON";

  Config cfg{cfg_json};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  // A SECOND peak at the same charge. makeSyntheticPeakGroup makes one, and a one-peak group has
  // mz2 == mz1, which would make the commanded window the margin ALONE -- symmetric about the centre,
  // so "inside the new window but outside the old" would have no room to exist. The single-peak case
  // is the third section's subject, not this one's.
  PeakGroup pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  FLASHHelperClasses::LogMzPeak lp2;
  lp2.mz = 802.0;
  lp2.abs_charge = 3;
  pg.push_back(lp2);

  std::vector<ScanCommand> cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 4)   // CE-0 baseline + CE 20 / 25 / 30
  ABORT_IF(cmds.size() != 4)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(group.msn_level, 2)
  TEST_REAL_SIMILAR(group.precursor_mz, 801.0)

  // THE INVARIANT, stated rather than inferred from the sum below: an MS2 exploration window is the
  // measured span plus the margin on each side (2.8 Th here). Expressed through the shared symbol, so
  // it says "the same margin buildMS2 uses" and not "0.8" -- a retuned constant should move this
  // section through optimal_window_margin_has_one_definition, not silently agree with a stale literal.
  TEST_REAL_SIMILAR(group.isolation_width, (802.0 - 800.0) + 2.0 * optimal_window_margin_)

  // ...and the COMMANDED side of the same claim, read off the baseline command the group actually
  // emitted. This is ADR-0026 decision 2 in one line: the scored interval IS the isolated interval. It
  // is a second, independent way for a reverted MS2 arm to go red (group 2.0 against command 2.8), and
  // unlike the sum below it needs no fed peaks at all.
  TEST_REAL_SIMILAR(cmds[0].stages[0].precursor_mz, group.precursor_mz)
  TEST_REAL_SIMILAR(cmds[0].stages[0].isolation_width, group.isolation_width)

  std::vector<double> baseline_mzs  = {799.8, 801.0, 802.2, 810.0};
  std::vector<double> baseline_ints = {300.0, 1000.0, 200.0, 700.0};
  exploration.feedResult(cmds[0].scan_id, baseline_mzs.data(), baseline_ints.data(),
                         static_cast<int>(baseline_mzs.size()), 0.5, queue);

  // THE ASSERTION THIS SECTION EXISTS FOR. 1500.0 with the margin, 1000.0 without it -- so a revert
  // does not merely widen a tolerance, it changes the number. Read back off the group rather than out
  // of the returned info because baseline_intensity is where the denominator of every later ratio in
  // this group is stored (Exploration.cpp:427).
  auto after = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(after.has_baseline, true)
  TEST_REAL_SIMILAR(after.baseline_intensity, 1500.0)
  TEST_EQUAL(after.baseline_failed, false)   // 1500 > 0, so this is the live path and not the abort one
}
END_SECTION

START_SECTION(optimal_window_margin_has_one_definition)
{
  // DRIFT GUARD for the constant this push collapsed from three TU-local copies -- ScanCommandQueue.cpp,
  // PrecursorSelection.cpp, and FragmentAnalysis.cpp's anonymous namespace, all `.4`, none of them aware
  // of the others -- onto the single inline constexpr at NotchSelection.h:123. Collapsing two
  // definitions into one earns a permanent bidirectional guard; this is it.
  //
  // FAILS WHEN a second definition reappears and drifts. The literal check below catches a retune of
  // the shared symbol, but on its own it cannot see a private copy: a `static const double
  // optimal_window_margin_ = .5;` back inside ScanCommandQueue.cpp shadows the header for that TU
  // alone, and every assertion that reads the header would keep passing. The buildMS2 assertion is the
  // half that catches it -- it compares the width ScanCommandQueue ACTUALLY computed against one
  // derived from the shared symbol, so the two sides have to agree numerically or the section is red.
  TEST_REAL_SIMILAR(optimal_window_margin_, 0.4)

  // mass_count, deliberately, and not remaining_precursor: this section is about GEOMETRY, and
  // mass_count is untouched by either ADR-0026 config rejection (Config.cpp:940-946 and :962-983), so the fixture
  // needs no overrides map to keep Rejection A quiet and nothing here can be misread as a test of the
  // sweep contract.
  const std::string cfg_json = R"JSON({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": {
      "rank_by": "qscore", "max_precursors": 3,
      "exploration": { "metric": "mass_count", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0 }
    },
    "characterization": { "mode": "off" },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 }
    }
  })JSON";

  Config cfg{cfg_json};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  // getMzRange(4) == (600.0, 603.0): a span of 3.0, deliberately different from the other two sections'
  // so a fixture copy-pasted between them cannot pass by coincidence.
  PeakGroup pg = makeSyntheticPeakGroup(600.0, 2396.0, 4);
  FLASHHelperClasses::LogMzPeak lp2;
  lp2.mz = 603.0;
  lp2.abs_charge = 4;
  pg.push_back(lp2);

  // Computed from the SHARED symbol, which is the whole point: 3.8 written as a literal would agree
  // with a divergent private copy exactly as happily as with the header.
  const double expected_width = (603.0 - 600.0) + 2.0 * optimal_window_margin_;

  // THE ScanCommandQueue SIDE. buildMS2 pads mz1/mz2 (ScanCommandQueue.cpp:293-294) AFTER taking the
  // centre (:292), so the centre is the un-padded midpoint and only the width carries the margin --
  // asserted separately, because a margin applied before the centre would leave the width right and the
  // centre wrong.
  ScanCommand ms2 = queue.buildMS2(pg, 4, cfg.level(2).scans[0], 2, 0);
  TEST_EQUAL(ms2.msn_level, 2)
  TEST_EQUAL(ms2.num_stages, 1)
  TEST_REAL_SIMILAR(ms2.stages[0].precursor_mz, 601.5)
  TEST_REAL_SIMILAR(ms2.stages[0].isolation_width, expected_width)

  // THE Exploration SIDE -- the invariant ADR-0026 decision 2 establishes, on the same peak group.
  // Before the correction this read 3.0 against the command's 3.8.
  std::vector<ScanCommand> cmds = exploration.initiate(2, pg, 4, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 4)   // CE-0 baseline + CE 20 / 25 / 30
  ABORT_IF(cmds.size() != 4)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_REAL_SIMILAR(group.precursor_mz, 601.5)
  TEST_REAL_SIMILAR(group.isolation_width, expected_width)
  // The two sides against each other, with no literal in between: whatever value the margin takes, the
  // scored window and the commanded window are the same interval.
  TEST_REAL_SIMILAR(group.isolation_width, cmds[0].stages[0].isolation_width)
  TEST_REAL_SIMILAR(group.precursor_mz, cmds[0].stages[0].precursor_mz)
}
END_SECTION

START_SECTION(single_isotope_charge_yields_non_degenerate_scan_range)
{
  // THE REASON THE MARGIN CORRECTION SHIPS IN THIS PUSH RATHER THAN A LATER ONE, and the section that
  // pins it. A charge state the deconvolution resolved at ONE isotope gives mz2 == mz1, so the bare
  // span is 0. Without the correction Exploration.cpp:165-167 hands the ADR-0026 binding
  // isolation_width == 0, the binding at :254-261 writes first_mass = last_mass = precursor_mz, and
  // BOTH are positive -- which is exactly the shape ScanFactory.cs:245-247 forwards:
  //     if (cmd.FirstMass > 0) p.FirstMass = ...;  if (cmd.LastMass > 0) p.LastMass = ...;
  //     if (cmd.FirstMass > 0 && cmd.LastMass > 0) p.ScanRangeMode = "DefineMZRange";
  // A zero-width DefineMZRange reaches the instrument. Nothing on either side of the bridge rejects it,
  // and no other section in this file would have noticed, so the binding could not land alone.
  //
  // THE MARGIN IS THE FLOOR HERE: MS2 has no 2.0 Th minimum of its own. That arm is msn_level >= 3 only
  // (Exploration.cpp:165-166), mirrored by buildMS3's stage[1] re-floor at ScanCommandQueue.cpp:451, so
  // an MS3 sweep on the same degenerate span is protected by something this section cannot see.
  //
  // FAILS WHEN the MS2 margin is dropped: last_mass - first_mass collapses to exactly 0, which the
  // `> 0.0` half reports as a degenerate range and the `2 * margin` half reports as the wrong width.
  const std::string cfg_json = R"JSON({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": {
      "rank_by": "qscore", "max_precursors": 3,
      "exploration": { "metric": "remaining_precursor", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0,
                       "overrides": { "analyzer": "IonTrap" } }
    },
    "characterization": { "mode": "off" },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 }
    }
  })JSON";

  Config cfg{cfg_json};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  // ONE peak, and it stays one peak. Adding a second is the workaround
  // remaining_precursor_binds_scan_range_ms2 uses to get an independent lo/hi pair, and applying it
  // here would delete the case: a single-isotope charge is the input, not an inconvenience.
  // makeSyntheticPeakGroup builds precisely that, so no helper needed changing.
  PeakGroup pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  const std::tuple<double, double> mz_range = pg.getMzRange(3);
  const double mz1 = std::get<0>(mz_range);
  const double mz2 = std::get<1>(mz_range);
  TEST_REAL_SIMILAR(mz2 - mz1, 0.0)      // the degenerate span, asserted rather than assumed
  TEST_REAL_SIMILAR(mz1, 800.0)          // and not the (-1, -10) no-such-charge sentinel

  // ms_settings.ms2 names NO range, so its 0/0 default is what a command carries when the binding does
  // not fire -- which keeps the loop below honest in the other direction too: a missing binding gives
  // 0 - 0 == 0 and fails on the same two assertions as a missing margin.
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].first_mass, 0.0)
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].last_mass, 0.0)
  TEST_EQUAL(cfg.level(2).overrides.count("first_mass"), 0u)   // the decision-5 escape hatch is NOT engaged
  TEST_EQUAL(cfg.level(2).overrides.count("last_mass"), 0u)

  std::vector<ScanCommand> cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 4)   // CE-0 baseline + CE 20 / 25 / 30
  ABORT_IF(cmds.size() != 4)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(static_cast<int>(group.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))
  TEST_REAL_SIMILAR(group.precursor_mz, 800.0)
  // Span 0 plus the margin on each side. Nothing else contributes at MS2, so this is the floor in its
  // entirety: 0.8 Th, neither 2.0 nor 0.
  TEST_REAL_SIMILAR(group.isolation_width, 2.0 * optimal_window_margin_)

  // EVERY command, the CE-0 baseline included: the binding is written onto base_config above the
  // variant loop (Exploration.cpp:258-260), so a degenerate range would reach the reference scan too.
  for (const ScanCommand& c : cmds)
  {
    TEST_EQUAL(c.msn_level, 2)
    TEST_EQUAL(c.last_mass - c.first_mass > 0.0, true)
    TEST_REAL_SIMILAR(c.last_mass - c.first_mass, 2.0 * optimal_window_margin_)
    TEST_REAL_SIMILAR(c.first_mass, group.precursor_mz - group.isolation_width / 2.0)
    TEST_REAL_SIMILAR(c.last_mass, group.precursor_mz + group.isolation_width / 2.0)
  }
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
  // DELIBERATELY 2 tol ENTRIES -- do not "fix" this to 3. levels_ materialises {1,2,3}
  // unconditionally and requires tol.size() >= 3, and proving that throw is the whole point
  // of this section. A bulk sweep that extended every 2-entry tol array in the tree defanged
  // it once already, and the assertion then passes vacuously because nothing throws.
  const char* insufficient_tol_config = R"({
  "deconvolution": {
    "min_charge": 4,
    "max_charge": 50,
    "min_mass": 500,
    "max_mass": 50000,
    "tol": [
      10,
      10
    ]
  },
  "faims": {
    "cv_values": []
  },
  "scheduling": {},
  "files": {},
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "rank_by": "qscore",
    "max_precursors": 3
  },
  "characterization": {
    "mode": "off",
    "max_targets": 3
  },
  "ms_settings": {
    "ms1": {
      "analyzer": "Orbitrap",
      "first_mass": 500,
      "last_mass": 2000,
      "resolution": 120000
    },
    "ms2": {
      "analyzer": "Orbitrap",
      "activation": "HCD",
      "collision_energy": 29,
      "resolution": 120000
    },
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "CID",
      "collision_energy": 25,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false
  }
}
)";
  TEST_EXCEPTION(std::invalid_argument, Config cfg{std::string(insufficient_tol_config)})
}
END_SECTION

START_SECTION(selection_without_next_level_scan_config_rejected)
{
  // MS3 is ON but no ms_settings.ms3 scan config is defined -> must throw at construction.
  // Guards the OOB read of next_cfg.scans[0] in Exploration::initiateNextLevel.
  //
  // `mode` MUST NOT be "off" here. The schema migration mechanically rewrote this fixture's
  // ms3.selection:"none" into mode:"off", which silently DEFANGED the test: with MS3 off no
  // level-3 scan config is required, so the expected throw correctly stopped happening and the
  // section failed with "no exception thrown". The premise is "MS3 reachable, nowhere to put it".
  const char* missing_next_scan_config = R"({
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
  "faims": {
    "cv_values": []
  },
  "scheduling": {},
  "files": {},
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "rank_by": "qscore",
    "max_precursors": 3
  },
  "characterization": {
    "mode": "ambiguity",
    "protein_sequence": "PEPTIDER",
    "max_targets": 3
  },
  "ms_settings": {
    "ms1": {
      "analyzer": "Orbitrap",
      "first_mass": 500,
      "last_mass": 2000,
      "resolution": 120000
    },
    "ms2": {
      "analyzer": "Orbitrap",
      "activation": "HCD",
      "collision_energy": 29,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false
  }
}
)";
  TEST_EXCEPTION(std::invalid_argument, Config cfg{std::string(missing_next_scan_config)})
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
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(cmds.size() > 0, true)

  // Verify variants have ETD activation type
  auto group = ExplorationTestAccess::group(exploration,1);
  for (const auto& v : group.variants)
  {
    TEST_STRING_EQUAL(v.activation_type, "ETD")
  }

  // Feed result — this exercises the full scoring chain:
  // feedResultImpl_ → computeExplorationScore_ → computeFragmentMatch_
  // which now passes v.activation_type instead of hardcoded "HCD"
  DeconvolvedSpectrum ds = makeSyntheticDeconv(1, 5);
  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info = ExplorationTestAccess::feedResult(exploration,tracking_id, ds, 1.0, queue);

  // The key assertion: feedResult completed without error and the activation type
  // was correctly propagated through the chain.
  TEST_EQUAL(info.fragmentation.activation_type, "ETD")
  TEST_EQUAL(info.group.group_id > 0, true)
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
  auto cmds = exploration.initiate(2, pg, 3, queue);

  TEST_EQUAL(cmds.size() > 0, true)              // ETD exploration produced variants
  TEST_EQUAL(exploration.activeGroupCount(), 1)

  auto group = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group.msn_level, 2)
  TEST_EQUAL(group.variants.size() > 0, true)

  // c/z proxy: every variant is ETD. rt-sweep: reaction_time takes the swept values
  // {5,10,15}; collision_energy is the (constant) base-config CE, not a CE sweep.
  std::set<double> seen_rts;
  int real_variants = 0;
  for (const auto& v : group.variants)
  {
    if (v.is_baseline) continue;  // CE-0 baseline (RT 0): a reference, not an RT-sweep point
    ++real_variants;
    TEST_STRING_EQUAL(v.activation_type, "ETD")
    TEST_EQUAL(std::isfinite(v.reaction_time), true)
    seen_rts.insert(v.reaction_time);
  }
  // rt_min=5, rt_max=15, rt_step=5 -> 3 distinct reaction-time variants (+ the prepended CE-0 baseline).
  TEST_EQUAL(static_cast<int>(group.variants.size()), 4)  // baseline + 3 RT variants
  TEST_EQUAL(real_variants, 3)
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
  auto cmds = exploration.initiate(2, pg, 3, queue);

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

      // Drive the full MS1 -> MS2-exploration cycle over the one canonical interleaved engine-id-echo
      // contract (runInterleaved): feed survey scans stamped with the engine's own ids; one selects a
      // precursor and forms the ETD exploration group, whose CE-/RT-sweep variants are drained and fed
      // back on the real cytC ladder. ms2_scans[0] is fed for every MS2 command.
      AcqResult a = runInterleaved(ida, ms1_scans, std::vector<ScanData>{ms2_scans[0]});
      delete ida;

      // Re-express etd_variant_cmds>0 over the harness result: count MS2 commands that are ETD exploration
      // variants (scan_description[3]=='E'); the ETD activation is wired by the fixture's exploration block.
      int etd_variant_cmds = 0;
      for (const auto& c : a.ms2_cmds)
      {
        std::string d(c.scan_description);
        if (d.size() >= 4 && d[3] == 'E') ++etd_variant_cmds;
      }
      // The fixture drives an ETD exploration -> at least one 'E' MS2 variant command.
      TEST_EQUAL(etd_variant_cmds > 0, true)
    }
  }
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
    },
    "agc_interval_seconds": 999999
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "../../FlashIDA/test-data/configs/inclusion_cytc.txt",
    "ptm_list": ""
  },
  "conditional_ms2": false,
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "inclusion",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3
  },
  "characterization": {
    "mode": "ambiguity",
    "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
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
  std::string cfg_str(inclusion_ms3_config);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  // H-full id-chaining drive over the one canonical contract (runInterleaved), replacing the
  // inlined `for (it<300; idle<3)` driver that re-implemented it. We feed a SINGLE strong MS1 scan
  // (ms1_scans[1] = scan 134, the cytC envelope; ms1_scans[0] = scan 132 is a weak edge scan that
  // selects 0 precursors), so an MS1 re-survey after it is fed counts as idle (avoids RT
  // self-exclusion churn). MS2 = ms2_scans[0]; the MS3 feed uses the MS2-as-MS3 shortcut
  // (ms3_ion_map == nullptr), matching the original inline driver's MS3 feed.
  AcqResult a = runInterleaved(ida, std::vector<ScanData>{ms1_scans[1]},
                               std::vector<ScanData>{ms2_scans[0]}, nullptr, 300);
  delete ida;

  // Hard asserts: the inclusion-mode cytC acquisition produced MS2 and MS3 commands.
  TEST_EQUAL(a.ms2_cmds.size() >= 1, true)
  TEST_EQUAL(a.ms3_cmds.size() >= 1, true)

  // Every MS3 parent_tracking_id must resolve to an emitted MS2 scan id. The MS2 ids
  // are the leading 3 chars of each emitted MS2 scan_description; an MS3 command stores
  // its parent MS2 id in parent_scan_id (3 chars + NUL).
  std::set<std::string> emitted_ms2_ids;
  for (const auto& m2 : a.ms2_cmds)
  {
    std::string d(m2.scan_description);
    if (d.size() >= 3) emitted_ms2_ids.insert(d.substr(0, 3));
  }
  bool all_ms3_parents_resolve = true;
  std::string unresolved;
  for (const auto& m3 : a.ms3_cmds)
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
         << ", MS3 cmds = " << a.ms3_cmds.size()
         << (all_ms3_parents_resolve ? std::string("")
                                     : (", first unresolved MS3 parent = '" + unresolved + "'")))
  TEST_EQUAL(all_ms3_parents_resolve, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Exploration x co-isolation. This combination had NO coverage on either CI path: every charge-mode
// golden runs without exploration, and every exploration test leaves precursor_charges at "single".
//
// What it guards: an MS3 cascaded from an MS2-EXPLORATION winner takes its stage 0 from
// group.variants[best_idx].cmd (Exploration.cpp:687 -> initiateNextLevel -> buildMS3), which is a
// different command object than the regular path's resolved parent_ctx. If that object ever stops
// carrying the notch block -- a variant rebuilt from descriptors, a snapshot taken too early, a
// stage-0 re-selection -- the MS3 replays a FRACTION of the precursor it was supposed to regenerate
// the fragment from, and nothing else in the tree notices.
//
// Two drives, identical but for the two charge keys, because the notch assertion below would pass
// vacuously on data with single-charge envelopes: the "single" drive must acquire the same way and
// carry ZERO notches. That is what makes the multiplexed drive's counts evidence.
/////////////////////////////////////////////////////////////

START_SECTION(exploration_multiplexed_ms3_inherits_the_variant_notch_set)
{
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.size() < 2 || ms2_scans.empty())

  // Same cytC MS2->MS3 acquisition as the section above, with targeting off and an MS2 CE sweep added,
  // so the MS3s arrive via the exploration cascade rather than the regular path. @PC@/@FC@ are spliced
  // per drive. R"JSON(...)JSON" is delimited deliberately: an undelimited R"(...)" ends at the first
  // ")" inside the payload, which has silently truncated a config in this tree before.
  const std::string cfg_template = R"JSON({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.9,
    "min_charge": 4,
    "max_charge": 50,
    "min_mass": 500,
    "max_mass": 50000,
    "tol": [10, 10, 10]
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
    },
    "agc_interval_seconds": 999999
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
    "max_precursors": 1,
    "precursor_charges": "@PC@",
    "exploration": {
      "metric": "mass_count",
      "ce_min": 20,
      "ce_max": 40,
      "ce_step": 5,
      "activations": ["HCD"]
    }
  },
  "characterization": {
    "mode": "ambiguity",
    "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "max_targets": 3,
    "fragment_charges": "@FC@"
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

  // inclusionPinCytc is the canonical cytC recipe every other cytC section here uses -- it swaps
  // targeting none -> inclusion and points inclusion_list at inclusion_cytc.txt. Required, not
  // decorative: under targeting "none" the top-qscore precursor of this survey is not reliably cytC,
  // and if it is not, the fresh57 ladder matches nothing and no MS3 fires at all. Its other two
  // replacements (M-start sequence, mode) are already satisfied above, so they no-op.
  auto instantiate = [&cfg_template](const std::string& pc, const std::string& fc) {
    std::string s = cfg_template;
    s.replace(s.find("@PC@"), 4, pc);
    s.replace(s.find("@FC@"), 4, fc);
    return inclusionPinCytc(s);
  };

  // Notch charge list of one cascade stage, as ints -- structural, so no float tolerance is involved.
  auto notch_charges = [](const ScanCommand& c, int k) {
    std::vector<int> out;
    auto r = notchesForStage(c, k);
    for (int i = 0; i < r.second; ++i) out.push_back(r.first[i].charge_state);
    return out;
  };

  auto drive = [&](const std::string& cfg) {
    FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg.c_str()));
    AcqResult a = runInterleaved(ida, std::vector<ScanData>{ms1_scans[1]},
                                 std::vector<ScanData>{ms2_scans[0]}, nullptr, 300);
    delete ida;
    return a;
  };

  // ---- Drive 1: "single". The control. Same acquisition, no notch anywhere. ----
  AcqResult s = drive(instantiate("single", "single"));
  TEST_EQUAL(s.ms2_cmds.size() >= 1, true)
  TEST_EQUAL(s.ms3_cmds.size() >= 1, true)

  int single_notches = 0;
  for (const auto& c : s.ms2_cmds) single_notches += c.stage0_notch_count + c.stage1_notch_count;
  for (const auto& c : s.ms3_cmds) single_notches += c.stage0_notch_count + c.stage1_notch_count;
  STATUS("single: ms2=" << s.ms2_cmds.size() << " ms3=" << s.ms3_cmds.size()
         << " notches=" << single_notches)
  TEST_EQUAL(single_notches, 0)

  // ---- Drive 2: "multiplexed" at both levels. ----
  AcqResult m = drive(instantiate("multiplexed", "multiplexed"));
  TEST_EQUAL(m.ms2_cmds.size() >= 1, true)
  TEST_EQUAL(m.ms3_cmds.size() >= 1, true)

  // The exploration path co-isolates at all: at least one MS2 variant carries notches. Without this the
  // inheritance check below would be satisfied by every MS3 inheriting an empty set.
  int notched_ms2 = 0, widest_ms2 = 0;
  std::map<std::string, std::vector<int>> ms2_stage0_by_id;   // 3-char id -> its stage-0 notch charges
  for (const auto& c : m.ms2_cmds)
  {
    if (c.stage0_notch_count > 0) ++notched_ms2;
    if (c.stage0_notch_count > widest_ms2) widest_ms2 = c.stage0_notch_count;
    // No stage can exceed its OWN block; a larger value means the two stages are sharing again.
    TEST_EQUAL(c.stage0_notch_count <= MAX_NOTCHES_PER_STAGE, true)
    std::string d(c.scan_description);
    if (d.size() >= 3) ms2_stage0_by_id[d.substr(0, 3)] = notch_charges(c, 0);
  }
  STATUS("multiplexed: ms2=" << m.ms2_cmds.size() << " notched_ms2=" << notched_ms2
         << " widest_stage0=" << widest_ms2 << " ms3=" << m.ms3_cmds.size())
  TEST_EQUAL(notched_ms2 >= 1, true)

  // Every MS3 replays its parent's FULL precursor isolation, not just the anchor.
  int checked = 0, mismatched = 0;
  for (const auto& c : m.ms3_cmds)
  {
    std::string parent(c.parent_scan_id);
    auto pit = ms2_stage0_by_id.find(parent);
    if (pit == ms2_stage0_by_id.end()) continue;   // parent resolution is the section above's job
    ++checked;
    if (notch_charges(c, 0) != pit->second) ++mismatched;
    TEST_EQUAL(c.stage0_notch_count <= MAX_NOTCHES_PER_STAGE, true)
    TEST_EQUAL(c.stage1_notch_count <= MAX_NOTCHES_PER_STAGE, true)
  }
  STATUS("MS3 stage-0 inheritance: checked=" << checked << " mismatched=" << mismatched)
  TEST_EQUAL(checked >= 1, true)
  TEST_EQUAL(mismatched, 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// B1: model-driven MS3 charge floor (selection_strategy.ms2.min_charge).
// Real cytC inclusion MS2->MS3 (M-start proteoform), driven through the standard model-driven path
// (M3: tracker fed in initiateNextLevel -> planNextScans). The charge floor is applied at
// Exploration.cpp:802 (std::abs(target.frag_charge) < charge_floor). Two sibling tests so the 0 in
// the floor-filters case is NOT vacuous: the no-floor sibling proves the model DOES emit MS3 here.
/////////////////////////////////////////////////////////////

START_SECTION(ms3_min_charge_floor_filters_all_targets)
{
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.size() < 2 || ms2_scans.empty())

  // Clone of inclusion_ms3_config (A5) with selection_strategy.ms2.min_charge = 99 —
  // impossible floor that filters every charge-1-3 MS3 target at Exploration.cpp:802.
  const char* cfg99 = R"({
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
    },
    "agc_interval_seconds": 999999
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "../../FlashIDA/test-data/configs/inclusion_cytc.txt",
    "ptm_list": ""
  },
  "conditional_ms2": false,
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "inclusion",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3
  },
  "characterization": {
    "mode": "ambiguity",
    "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "max_targets": 3,
    "min_fragment_charge": 99
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
  std::string s99(cfg99);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(s99.c_str()));
  AcqResult a = runInterleaved(ida, std::vector<ScanData>{ms1_scans[1]},
                               std::vector<ScanData>{ms2_scans[0]}, nullptr, 300);
  delete ida;

  // The precursor still fires MS2 (inclusion selects it), but every model MS3 target has
  // charge 1-3 which is < 99, so the charge floor at Exploration.cpp:802 filters them all
  // -> zero MS3 commands.
  TEST_EQUAL(a.ms2_cmds.size() >= 1, true)
  TEST_EQUAL(a.ms3_cmds.size(), 0)
}
END_SECTION

START_SECTION(ms3_min_charge_default_emits_ms3)
{
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.size() < 2 || ms2_scans.empty())

  // Clone of inclusion_ms3_config (A5) with selection_strategy.ms2.min_charge = 1 —
  // a permissive floor that lets all fragment charges through; positive control proving
  // the model-driven path DOES emit MS3 so the zero in ms3_min_charge_floor_filters_all_targets
  // is real (not a structural dead-end).
  const char* cfg1 = R"({
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
    },
    "agc_interval_seconds": 999999
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "../../FlashIDA/test-data/configs/inclusion_cytc.txt",
    "ptm_list": ""
  },
  "conditional_ms2": false,
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "inclusion",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3
  },
  "characterization": {
    "mode": "ambiguity",
    "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "max_targets": 3,
    "min_fragment_charge": 1
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
  std::string s1(cfg1);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(s1.c_str()));
  AcqResult a = runInterleaved(ida, std::vector<ScanData>{ms1_scans[1]},
                               std::vector<ScanData>{ms2_scans[0]}, nullptr, 300);
  delete ida;

  // Positive control: with a permissive floor the model-driven path emits MS3 commands,
  // proving the zero in ms3_min_charge_floor_filters_all_targets is caused by the charge
  // floor at Exploration.cpp:802, not by a structural dead-end.
  TEST_EQUAL(a.ms2_cmds.size() >= 1, true)
  TEST_EQUAL(a.ms3_cmds.size() >= 1, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
