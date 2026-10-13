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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>          // a_trap_pre_scan_is_measured_never_identified
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>  // builds a real tracker to prove nothing was pooled
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

  // ---- remaining_precursor config probes (ADR-0026 as amended by ADR-0044) --------------------
  //
  // ADR-0026 once bound a remaining_precursor pre-scan's scan range to the isolation window it reads
  // and, as that narrowing's interlock, refused the metric with an empty overrides map. ADR-0044
  // withdrew both: a pre-scan reads out its level's configured range and overrides are optional
  // under every metric. What survives is the level-matched multiplexed rejection, on a ground of its
  // own -- the remaining-precursor ratio sums the ANCHOR's window only -- and the sections beside
  // ms3_protein_sequence_only_accepted pin that, plus what must keep loading.
  //
  // Spliced, never find/replaced. The sections commented at :2113-2115 died on `invalid string
  // position` doing surgery on a finished literal; a splice has no position to compute, so the only
  // way the probes below can fail is the rejection itself. Every knob a section varies is a
  // parameter here, so no section reaches into another section's JSON.

  // A non-empty exploration.overrides map on the Orbitrap. Overrides are optional (ADR-0044), so this
  // is now only a way to say "degraded pre-scans" without naming the ion trap -- which ADR-0045 lets
  // a remaining_precursor sweep do, but which would change what the section is about.
  const std::string sweep_overrides = R"(, "overrides": { "analyzer": "Orbitrap" })";

  // A remaining_precursor exploration block, ready to splice into either decision section. Pass ""
  // for `overrides` to probe the empty-map shape and `sweep_overrides` otherwise. One CE range serves
  // both levels: needsCollisionEnergy is true for HCD and CID alike, and Config::validate asks only
  // for ce_min < ce_max.
  std::string rpSweep(const std::string& overrides)
  {
    return R"("exploration": { "metric": "remaining_precursor", "ce_min": 20.0, "ce_max": 40.0, "ce_step": 5.0)"
           + overrides + R"( })";
  }

  // The characterization body that makes characterization.exploration reachable at level 3: the MS3
  // gate (mode != off) plus the sequence Config.cpp:1012-1016 then requires. Trailing comma included,
  // so a caller appends its own keys directly.
  const std::string ms3_on = R"("mode": "ambiguity", "protein_sequence": "PEPTIDER", )";

  // Minimal loadable config. `ps_extra` is appended inside precursor_selection (lead with a comma),
  // `charact` REPLACES the characterization body, and `with_ms3` adds the level-3 scan config that a
  // non-off mode requires (Config.cpp:1022-1027 -- without it Exploration::initiateNextLevel OOB-reads
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
    },
    "agc_interval_seconds": 9999999
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
    },
    "agc_interval_seconds": 9999999
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
    },
    "agc_interval_seconds": 9999999
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

  // ADR-0055 fixtures: ms3_selection_only_config, cytC-pinned and switched to exhaustive, but naming UBIQUITIN as
  // the protein -- so the real cytC MS2 identifies nothing, exactly like a contaminant MS2 under a cytC method.
  // Without @p ms2_exploration the mass_count sweep is removed and the precursor gets one regular MS2. Every
  // replace THROWS when it misses: a silent no-op would leave the fixture identifying and the sections vacuous.
  std::string unidentifiedExhaustiveCytc(bool ms2_exploration)
  {
    std::string cfg = inclusionPinCytc(ms3_selection_only_config);
    auto rep = [&cfg](const std::string& from, const std::string& to) {
      auto p = cfg.find(from);
      if (p == std::string::npos) throw std::runtime_error("unidentifiedExhaustiveCytc: no '" + from + "'");
      cfg.replace(p, from.size(), to);
    };
    rep("\"mode\": \"ambiguity\"", "\"mode\": \"exhaustive\"");
    rep("MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG");
    if (!ms2_exploration)
      rep(",\n    \"exploration\": {\n      \"metric\": \"mass_count\",\n      \"ce_min\": 20.0,\n"
          "      \"ce_max\": 40.0,\n      \"ce_step\": 5.0\n    }", "");
    return cfg;
  }

  // An MS3 description in the no-ion form -- nothing after "@{fragment charge}" -- which is what a 'u'/0 target
  // becomes on the wire (buildMS3's ion_index == 0 branch). A mapped b/y target would carry its ion here, so
  // "every MS3 is no-ion" is also the proof that the precursor really went unidentified.
  bool isNoIonDescription(const char* desc)
  {
    const std::string d(desc);
    const auto at = d.rfind('@');
    if (at == std::string::npos) return false;
    size_t i = at + 1;
    while (i < d.size() && std::isdigit(static_cast<unsigned char>(d[i]))) ++i;
    return i > at + 1 && i == d.size();
  }

  // Config with remaining_precursor exploration metric at MS2.
  //
  // The exploration.overrides map below MUST be non-empty: ADR-0026 rejects "remaining_precursor" with
  // an empty map at config load (Config.cpp:963-969, at BOTH levels). Such a sweep's pre-scans are
  // narrowed to the isolation window and then thrown away, so the schema now forces the author to
  // declare the settings they run at, which is what guarantees the ADR-0020 gate-#1 re-acquisition.
  // "Orbitrap" is deliberately a NO-OP: it is what ms_settings.ms2.analyzer already says, so
  // Exploration::initiate's base_config.applyOverrides (Exploration.cpp:129) writes back the value the
  // base config already carries and every fixture behaviour is unchanged. Same inert-override trick as
  // the injection at :1664.
  //
  // THE ONE BRANCH A NON-EMPTY MAP OPENS is ADR-0020 gate #1 -- `!level_config.overrides.empty() ||
  // measuring_ms3_sweep` at Exploration.cpp:749 -- and it is unreached here. The sections driven off
  // this fixture that complete their group all reach completion with NO WINNER, taking the
  // `best_idx < 0` arm and returning before the gate is ever evaluated. Since ADR-0029 they get there
  // by scoring rather than by aborting: an activation whose reference returned empty scores -1.0, and
  // winner selection seeds best_score at -1.0 and takes `score > best_score`, so nothing is eligible.
  //   * remaining_precursor_empty_baseline_bars_activation_from_winning and
  //     remaining_precursor_late_child_after_group_cleanup_is_noop -- de-referenced by construction,
  //     on a raw baseline whose peaks all sit outside the isolation window.
  //   * remaining_precursor_score_no_raw_data -- NOT a de-referenced section by name or by intent. It
  //     feeds all SIX variants and asserts activeGroupCount() == 0, and it lands on the same arm only
  //     because ExplorationTestAccess::feedResult hands feedResultImpl_ mzs == ints == nullptr, so
  //     precursorWindowIntensity_ short-circuits to 0.0 (Exploration.cpp:1215-1219) and the reference
  //     reads as empty.
  // Every other section stops feeding partway, so no group on this fixture has ever selected a winner.
  //
  // LIVE HAZARD, not a historical note. Give remaining_precursor_score_no_raw_data -- or any section
  // added here -- a raw-array baseline with in-window signal and a winner IS selected, and gate #1 is
  // then not merely reached but TRUE, because this fixture's overrides map is non-empty by
  // construction. Gate #1 and the MS2 cascade `else if (group.msn_level < 3)` (Exploration.cpp:817) are
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
    },
    "agc_interval_seconds": 9999999
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
  // Same ADR-0026 requirement as the MS2 fixture above -- the rejection at Config.cpp:963-969 loops over
  // every level, so characterization.exploration.overrides must be non-empty too. "Orbitrap" is again a
  // deliberate no-op (ms_settings.ms3.analyzer already says Orbitrap), so this fixture drives exactly the
  // behaviour it did before. Note the isolation-window narrowing itself is NOT gated on this map -- it is
  // suppressed only by an explicit first_mass/last_mass override (Exploration.cpp:256-257) -- so adding an
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
    },
    "agc_interval_seconds": 9999999
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

  // Config with mass_count exploration at MS3 and NO overrides key -- one of the two fixtures that
  // isolate ADR-0020 gate #2. Both of ADR-0020's re-acquisition reasons are MEASURING-metric-driven,
  // but only gate #2 keys on the metric; gate #1 keys on a non-empty overrides map, so proving gate #2
  // needs an empty one. mass_count drives through pre-deconvolved spectra (its score IS spec.size());
  // the remaining_precursor arm of the same section derives its fixture from
  // ms3_remaining_precursor_config with the overrides line erased, and drives through raw arrays.
  //
  // Written as a SELF-CONTAINED literal rather than a find/replace over ms3_exploration_config: the note
  // at :2113-2115 records four sections that died on `invalid string position` when the token they
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
    },
    "agc_interval_seconds": 9999999
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
    },
    "agc_interval_seconds": 9999999
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
    },
    "agc_interval_seconds": 9999999
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
    },
    "agc_interval_seconds": 9999999
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

  // --- ADR-0029 helpers: one un-fragmented reference per swept activation -------------------
  //
  // These live HERE, not beside the sections that use them: START_TEST opens `int main()`, so
  // everything past it is function-local and a namespace-scope helper cannot be declared there.

  // A sweep config whose every ADR-0029-relevant knob is a parameter: which activations are swept,
  // where the CE and RT sweeps start, and what the base scan config carries on the axes that are
  // NOT swept (which is what that activation's baseline must copy).
  //
  // ce_max/rt_max are derived as min + 20 because Config::validate rejects an empty range on a
  // swept axis; the steps are coarse (10) so a block stays countable by hand: 3 sweep points.
  //
  // `overrides` defaults to empty and MUST be passed as `sweep_overrides` for a
  // remaining_precursor sweep: Config::validate rejects that metric with an empty overrides map
  // (ADR-0026 -- such a sweep never keeps its pre-scans, so it has to declare the settings they run
  // at). Omitting it does not fail an assertion, it throws out of the Config constructor, which
  // ClassTest reports as a failed section with an EMPTY "Failed lines:" list.
  std::string sweepCfg(const std::string& activations, double ce_min, double rt_min,
                       const std::string& metric = "mass_count",
                       const std::string& ms2_activation = "HCD",
                       int ms2_collision_energy = 35,
                       const std::string& overrides = "")
  {
    std::ostringstream os;
    os << R"({
  "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 1,
                     "max_charge": 50, "min_mass": 600, "max_mass": 50000, "tol": [10, 10, 10] },
  "precursor_selection": {
    "rt_window": 180, "targeting": "none", "rank_by": "qscore", "max_precursors": 3,
    "exploration": {
      "metric": ")" << metric << R"(",
      "ce_min": )" << ce_min << R"(, "ce_max": )" << (ce_min + 20.0) << R"(, "ce_step": 10.0,
      "reaction_time_min": )" << rt_min << R"(, "reaction_time_max": )" << (rt_min + 20.0)
       << R"(, "reaction_time_step": 10.0,
      "activations": [)" << activations << R"(])" << overrides << R"(
    }
  },
  "characterization": { "mode": "off",
    "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
    "max_targets": 3 },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
    "ms2": { "analyzer": "Orbitrap", "activation": ")" << ms2_activation
       << R"(", "collision_energy": )" << ms2_collision_energy << R"(, "resolution": 120000 }
  }
})";
    return os.str();
  }

  /// Index of the first variant belonging to @p act, or -1. Each activation's variants are
  /// contiguous and the blocks appear in the order the activations were authored.
  int blockStart(const std::vector<Exploration::ExplorationVariant>& vs, const std::string& act)
  {
    for (size_t i = 0; i < vs.size(); ++i)
      if (vs[i].activation_type == act) return static_cast<int>(i);
    return -1;
  }

  /// An MS2 kinetic_model sweep (ADR-0047): HCD ce_min..ce_max/ce_step with the authored k2 JSON object and an
  /// optional overrides fragment (", \"overrides\": {...}").
  std::string kineticCfg(double ce_min, double ce_max, double ce_step,
                         const std::string& k2_json = R"({"HCD": 0.039046})", const std::string& overrides = "")
  {
    return std::string(R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": "qscore", "max_precursors": 3,
        "exploration": { "metric": "kinetic_model", "ce_min": )") + std::to_string(ce_min)
        + R"(, "ce_max": )" + std::to_string(ce_max) + R"(, "ce_step": )" + std::to_string(ce_step)
        + R"(, "activations": ["HCD"], "kinetic_model": { "k2": )" + k2_json + R"( })" + overrides + R"( } },
      "characterization": { "mode": "off" },
      "ms_settings": { "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
                       "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29,
                                "first_mass": 150, "last_mass": 2000 } } })";
  }

  /// An MS3 kinetic_model sweep: CID ce_min..ce_max/ce_step on characterization.exploration with an authored
  /// CID k2; mode off (the sweep is driven through initiate(3, ...) directly, which reads the level config).
  std::string kineticMs3Cfg(double ce_min, double ce_max, double ce_step, double k2 = 0.0388)
  {
    return std::string(R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": "qscore", "max_precursors": 3 },
      "characterization": { "mode": "off", "max_targets": 3,
        "exploration": { "metric": "kinetic_model", "ce_min": )") + std::to_string(ce_min)
        + R"(, "ce_max": )" + std::to_string(ce_max) + R"(, "ce_step": )" + std::to_string(ce_step)
        + R"(, "activations": ["CID"], "kinetic_model": { "k2": { "CID": )" + std::to_string(k2) + R"( } } } },
      "ms_settings": { "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
                       "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 },
                       "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25,
                                "first_mass": 200, "last_mass": 2000 } } })";
  }

  /// Feed one variant a two-peak spectrum whose precursor share is exactly `share`: that fraction of 1000
  /// at the anchor m/z, the rest 300 Th away. share < 0 feeds an EMPTY spectrum (no measurement).
  Exploration::FeedResultInfo feedShare(Exploration& e, ScanCommandQueue& q, const ScanCommand& cmd,
                                        double anchor_mz, double share)
  {
    const int tid = q.decode(std::string(cmd.scan_description).substr(0, 3));
    if (share < 0.0) return e.feedResult(tid, nullptr, nullptr, 0, 1.0, q);
    std::vector<double> mzs{anchor_mz, anchor_mz + 300.0}, ints{1000.0 * share, 1000.0 * (1.0 - share)};
    return e.feedResult(tid, mzs.data(), ints.data(), 2, 1.0, q);
  }

  /// A CE-only sweep config with an explicit range: ce_min..ce_max in steps of ce_step (ADR-0047 d.8).
  std::string sweepCfgRange(const std::string& activations, double ce_min, double ce_max, double ce_step)
  {
    return std::string(R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": "qscore", "max_precursors": 3,
        "exploration": { "metric": "mass_count", "ce_min": )") + std::to_string(ce_min)
        + R"(, "ce_max": )" + std::to_string(ce_max) + R"(, "ce_step": )" + std::to_string(ce_step)
        + R"(, "activations": [)" + activations + R"(] } },
      "characterization": { "mode": "off" },
      "ms_settings": { "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
                       "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 } } })";
  }
} // anonymous namespace



/////////////////////////////////////////////////////////////
// ADR-0042 -- the MONITOR scan: an MS1 that observes the source and decides nothing.
/////////////////////////////////////////////////////////////

namespace
{
  /// Config for the monitor-scan tests. Level 2 sweeps HCD 20-40 over ms1_standard -- the same shape
  /// cycle_time_suppression_during_exploration drives. The two monitor_ms1 blocks are passed verbatim
  /// (empty string = key absent) so one helper covers on/off, per-level, and the invalid-interval case.
  std::string monitorScanConfig(const std::string& ms2_monitor,
                                const std::string& ms3_monitor = "",
                                bool cycle_time_enabled = false,
                                const std::string& ms3_metric = "none")
  {
    const std::string ct = cycle_time_enabled ? "true" : "false";
    return std::string("{"
      "\"deconvolution\": { \"score_threshold\": 0.0, \"tqscore_threshold\": 0.9, \"min_charge\": 4,"
      "                     \"max_charge\": 50, \"min_mass\": 500, \"max_mass\": 50000, \"tol\": [10, 10, 10] },"
      "\"faims\": { \"cv_values\": [] },"
      "\"scheduling\": { \"cycle_time\": { \"enabled\": ") + ct + ", \"value_ms\": 1 },"
      "                  \"agc_interval_seconds\": 999999 },"
      "\"precursor_selection\": {"
      "  \"rank_by\": \"qscore\", \"max_precursors\": 3,"
      "  \"exploration\": { \"metric\": \"mass_count\", \"ce_min\": 20.0, \"ce_max\": 40.0, \"ce_step\": 5.0"
      + ms2_monitor + " }"
      "},"
      // protein_sequence is set unconditionally: an ms3_metric of "fragment_count" makes level 3 a
      // READING metric, and Config::validate requires a sequence for one. Harmless when the metric
      // is "none", which is this helper's default.
      "\"characterization\": { \"mode\": \"off\", \"max_targets\": 3,"
      "  \"protein_sequence\": \"PEPTIDERPEPTIDER\","
      "  \"exploration\": { \"metric\": \"" + ms3_metric + "\"" + ms3_monitor + " }"
      "},"
      "\"ms_settings\": {"
      "  \"ms1\": { \"analyzer\": \"Orbitrap\", \"first_mass\": 500, \"last_mass\": 2000,"
      "             \"resolution\": 120000, \"agc_target\": 800000, \"max_it\": 246 },"
      "  \"ms2\": { \"analyzer\": \"Orbitrap\", \"activation\": \"HCD\", \"collision_energy\": 30 },"
      "  \"ms3\": { \"analyzer\": \"Orbitrap\", \"activation\": \"CID\", \"collision_energy\": 25 }"
      "}"
      "}";
  }

  const char* kMonitorOn  = ", \"monitor_ms1\": { \"enabled\": true, \"interval_ms\": 30000 }";

  /// Pull commands off the drain until one carries @p want, or the bound is hit. Returns whether it
  /// was found; @p out holds it. Non-matching commands are DISCARDED, which is fine below because
  /// every caller snapshots engine state AFTER this returns.
  bool drainUntilRole(FLASHIda* ida, ScanRole want, ScanCommand& out, int max_pulls = 12)
  {
    for (int i = 0; i < max_pulls; ++i)
    {
      out = ScanCommand{};
      if (ida->getNextScanCommand(out) != 1) return false;
      if (roleOf(out) == want) return true;
    }
    return false;
  }
}

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

  // A sweep is exactly its grid points (ADR-0048 d.3): no CE-0 reference scan under any metric.
  TEST_EQUAL(static_cast<int>(group.variants.size()), 5)

  const double expected_ce[5] = {20.0, 25.0, 30.0, 35.0, 40.0};
  int real_count = 0;
  for (const auto& v : group.variants)
  {
    TEST_EQUAL(v.received, false)
    TEST_EQUAL(v.variant_index, real_count)
    TEST_REAL_SIMILAR(v.collision_energy, expected_ce[real_count])
    ++real_count;
  }
  TEST_EQUAL(real_count, 5)

  TEST_EQUAL(static_cast<int>(cmds.size()), 5)  // the 5 CE-sweep variants, nothing else

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

  TEST_EQUAL(static_cast<int>(cmds.size()), 5)  // 5 CE-sweep variants (all priority-2 'E'); no baseline (ADR-0048)
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
  ScanCommand ms2_ctx = queue.buildMS2(precursor_pg, 3, cfg.level(2).scans[0], 2, 0);

  auto fragment_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  ms2_ctx.faims_cv = -50.0;  // Item 1: CV travels via the context
  auto cmds = exploration.initiate(3, fragment_pg, 2, queue, nullptr, &ms2_ctx, 'y', 5);

  TEST_EQUAL(static_cast<int>(cmds.size()), 5)  // 5 CE-sweep variants (all MS3 buildMS3); no baseline (ADR-0048)

  for (int i = 0; i < 5; ++i)
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

// ADR-0054 d.5: an MS3 sweep on a c ion (an ETD winner's fragment) keeps its CONFIGURED metric, here fragment_count --
// every variant then scores 0 and ce_min wins, accepted knowingly. Only a target with NO frame ('u') is forced onto
// remaining_precursor (ADR-0023 d.11), and that stays so.
//
// FAILS UNDER: annotatability folded into MS3FragmentMatcher::isKnownIonClass -- the c sweep would be forced too.
START_SECTION(ms3_sweep_on_a_c_ion_keeps_its_configured_metric)
{
  Config cfg{std::string(ms3_exploration_config)};
  TEST_EQUAL(static_cast<int>(cfg.level(3).exploration), static_cast<int>(ExplorationMetric::FragmentCount))

  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);
  ScanCommand ms2_ctx = queue.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3, cfg.level(2).scans[0], 2, 0);
  TEST_EQUAL(exploration.initiate(3, makeSyntheticPeakGroup(500.0, 1000.0, 2), 2, queue, nullptr, &ms2_ctx, 'c', 5).empty(), false)
  TEST_EQUAL(static_cast<int>(ExplorationTestAccess::group(exploration, 1).exploration_metric),
             static_cast<int>(ExplorationMetric::FragmentCount))

  // the control: an unassigned target in the same config is still forced
  ScanCommandQueue queue_u(cfg);
  FragmentAnalysis fragments_u(cfg);
  Exploration exploration_u(cfg, fragments_u);
  ScanCommand ms2_ctx_u = queue_u.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3, cfg.level(2).scans[0], 2, 0);
  TEST_EQUAL(exploration_u.initiate(3, makeSyntheticPeakGroup(500.0, 1000.0, 2), 2, queue_u, nullptr, &ms2_ctx_u, 'u', 0).empty(), false)
  TEST_EQUAL(static_cast<int>(ExplorationTestAccess::group(exploration_u, 1).exploration_metric),
             static_cast<int>(ExplorationMetric::RemainingPrecursor))
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
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)  // 5 CE-sweep variants; no baseline (ADR-0048)
  TEST_EQUAL(exploration.activeGroupCount(), 1)

  // cmds[0..4] are the 5 variants (variant_index 0..4) with peak counts {2,4,6,8,3}. Feeding ALL 5
  // completes the group; the winner is the variant with 8 peaks (variant_index 3 == cmds[3]).
  std::vector<int> peak_counts = {2, 4, 6, 8, 3};
  Exploration::FeedResultInfo last_info;
  for (int i = 0; i < 5; ++i)
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
  // The winner stays identifiable via winner_tracking_id == variant_index 3's encoded command id (cmds[3]).
  // (getGroup throws here — the group was erased on completion — so read it from last_info.)
  TEST_EQUAL(last_info.group.winner_tracking_id.empty(), false)
  TEST_EQUAL(last_info.group.winner_tracking_id, std::string(cmds[3].scan_description).substr(0, 3))

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
  // (cmds[4], CE 35). mass_count scores spec.size(), so variant_index 3 — 8 peak groups, cmds[3],
  // CE 30 on the 15/20/25/30/35 grid — won. Reading it off cmds[3] rather than hard-coding 30 keeps
  // the assertion true if the fixture's sweep grid is ever retuned.
  TEST_REAL_SIMILAR(last_info.commands[0].stages[1].collision_energy,
                    cmds[3].stages[1].collision_energy)
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
  // BOTH MEASURING METRICS, EACH WITH AN EMPTY OVERRIDES MAP. The gate is
  // `!level_config.overrides.empty() || measuring_ms3_sweep` -- a disjunction, so proving the SECOND
  // term requires the first to be false. ADR-0026 once rejected "remaining_precursor + empty
  // overrides" at load, which left mass_count as the only way to reach gate #2; ADR-0044 lifted that
  // rejection, so the remaining_precursor arm at the end of this section is reachable again and pins
  // that gate #2 keys on isMeasuringMetric rather than on one enumerator.
  //
  // The mass_count arm feeds pre-deconvolved spectra (ExplorationTestAccess): mass_count scores
  // spec.size(), so the peak-count vector below IS the score vector, whereas remaining_precursor has
  // to be driven through raw isolation-window intensities (second arm). It also drops the
  // empty-reference hazard entirely -- only RemainingPrecursor divides by a reference, so under
  // mass_count a zero-intensity baseline is simply a score-0 variant here.
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
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)  // CE 20,25,30,35,40; no baseline (ADR-0048)

  // cmds[0..4] are the 5 variants. mass_count scores spec.size(), so the counts below are the scores
  // directly and the winner is the 8-peak variant at cmds[3].
  std::vector<int> peak_counts = {2, 4, 6, 8, 3};
  Exploration::FeedResultInfo last_info;
  for (int i = 0; i < 5; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, peak_counts[i]);
    int tid = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    last_info = ExplorationTestAccess::feedResult(exploration, tid, ds, static_cast<double>(i), queue);
  }

  TEST_EQUAL(exploration.activeGroupCount(), 0)                      // group completed
  TEST_EQUAL(last_info.group.winner_tracking_id.empty(), false)      // a winner was chosen
  TEST_EQUAL(last_info.group.winner_tracking_id, std::string(cmds[3].scan_description).substr(0, 3))

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

  // ---- the remaining_precursor arm (ADR-0044 decision 3): same gate, other measuring metric ----
  // FAILS WHEN the overrides-mandatory rejection returns (the config below no longer loads), or when
  // gate #2 is narrowed to mass_count. Fixture: ms3_remaining_precursor_config with its one overrides
  // line erased -- a splice on a located token, not a find/replace over the whole literal.
  {
    std::string rp_json = ms3_remaining_precursor_config;
    const std::string ov_line = R"("overrides": { "analyzer": "Orbitrap" })";
    const size_t ov_pos = rp_json.find(ov_line);
    ABORT_IF(ov_pos == std::string::npos)
    const size_t comma = rp_json.rfind(',', ov_pos);   // the comma that ended "ce_step": 5.0
    ABORT_IF(comma == std::string::npos)
    rp_json.erase(comma, ov_pos + ov_line.size() - comma);

    Config rp_cfg{rp_json};
    TEST_EQUAL(rp_cfg.level(3).overrides.empty(), true)
    TEST_EQUAL(static_cast<int>(rp_cfg.level(3).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))

    ScanCommandQueue rp_queue(rp_cfg);
    FragmentAnalysis rp_fragments(rp_cfg);
    Exploration rp_expl(rp_cfg, rp_fragments);
    auto rp_frag_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
    ScanCommand rp_ms2_ctx = rp_queue.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3,
                                               rp_cfg.level(2).scans[0], 2, 0);
    auto rp_cmds = rp_expl.initiate(3, rp_frag_pg, 2, rp_queue, nullptr, &rp_ms2_ctx);
    TEST_EQUAL(static_cast<int>(rp_cmds.size()), 5)   // CE 20,25,30,35,40; no baseline (ADR-0048)
    ABORT_IF(rp_cmds.size() != 5)
    for (const auto& c : rp_cmds) TEST_STRING_EQUAL(c.analyzer, "Orbitrap")   // no override: the level's own analyzer

    // Raw arrays: one peak at the window centre, two 50 Th outside it. The precursor SHARE is the centre
    // peak over the whole scan (ADR-0048): 900/1900, 700/1700, 400/1400, 100/1100, 50/1050 = 0.474 /
    // 0.412 / 0.286 / 0.091 / 0.048 against target 0.1 give scores 0.626 / 0.688 / 0.814 / 0.991 /
    // 0.948, so the CE-35 variant (rp_cmds[3]) wins.
    auto rp_group = ExplorationTestAccess::group(rp_expl, 1);
    const double c_mz = rp_group.precursor_mz;
    std::vector<double> mzs{c_mz - 50.0, c_mz, c_mz + 50.0};
    std::vector<std::vector<double>> ints{{300.0, 900.0, 700.0}, {300.0, 700.0, 700.0},
                                          {300.0, 400.0, 700.0}, {300.0, 100.0, 700.0}, {300.0, 50.0, 700.0}};
    Exploration::FeedResultInfo rp_last;
    for (int i = 0; i < 5; ++i)
      rp_last = rp_expl.feedResult(rp_cmds[i].scan_id, mzs.data(), ints[i].data(), 3, 0.5 * (i + 1), rp_queue);

    TEST_EQUAL(rp_expl.activeGroupCount(), 0)
    TEST_EQUAL(rp_last.group.winner_tracking_id, ScanCommandQueue::encode(rp_cmds[3].scan_id))
    TEST_EQUAL(static_cast<int>(rp_last.commands.size()), 1)   // gate #2, with gate #1 provably false
    ABORT_IF(rp_last.commands.empty())
    TEST_EQUAL(rp_last.commands[0].msn_level, 3)
    TEST_EQUAL(std::string(rp_last.commands[0].scan_description)[3], 'R')
    TEST_REAL_SIMILAR(rp_last.commands[0].stages[1].collision_energy, 35.0)
    TEST_EQUAL(static_cast<int>(rp_last.ms2_context_cache.size()), 1)
  }
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
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)   // no baseline (ADR-0048)

  std::vector<int> peak_counts = {2, 4, 6, 8, 3};
  Exploration::FeedResultInfo last_info;
  for (int i = 0; i < 5; ++i)
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
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)  // 5 CE-sweep variants; no baseline (ADR-0048)

  std::vector<double> scores = {1.0, 3.0, 2.0, 5.0, 0.0};
  Exploration::FeedResultInfo last_info;
  for (int i = 0; i < 5; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, static_cast<int>(scores[i]));
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    last_info = ExplorationTestAccess::feedResult(exploration,tracking_id, ds, static_cast<double>(i), queue);
  }

  // No raw arrays reached the engine (the deconvolution bypass), so no share could be read: -1 = N/A.
  TEST_REAL_SIMILAR(last_info.metric.precursor_share, -1.0)
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
  // initiateNextLevel, Exploration.cpp:817-822).
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
  // a production-MS2 re-acquisition (Exploration.cpp:749-810) rather than MS3-direct. Feeding
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

START_SECTION(unidentified_regular_ms2_dispatches_exhaustive_ms3)
{
  // ADR-0055, no exploration at either level -- the shape of the only committed exhaustive method. The pinned
  // cytC precursor is acquired but the method names ubiquitin, so its one MS2 identifies nothing. Before
  // ADR-0055 planNextScans answered reason=unidentified_precursor and no MS3 left; now the MS2's own masses are
  // the pool, every target is 'u'/0, and every MS3's stage 0 replays that MS2 (ms_settings.ms2: HCD at 29).
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  std::string cfg_str = unidentifiedExhaustiveCytc(/*ms2_exploration=*/false);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  const std::map<std::string, std::vector<ScanData>> no_ms3_fixtures;  // non-null + empty => record, never feed
  const int budget = 256 + 64 * static_cast<int>(ms1_scans.size() + 1);
  AcqResult a = runInterleaved(ida, ms1_scans, std::vector<ScanData>{ms2_scans[0]}, &no_ms3_fixtures,
                               budget, /*single_group_only=*/true);
  delete ida;

  TEST_EQUAL(a.ms3_cmds.empty(), false)
  for (const ScanCommand& c : a.ms3_cmds)
  {
    TEST_EQUAL(isNoIonDescription(c.scan_description), true)
    TEST_REAL_SIMILAR(c.stages[0].collision_energy, 29.0)
  }
}
END_SECTION

START_SECTION(unidentified_ms2_sweep_dispatches_exhaustive_ms3_from_the_metric_winner)
{
  // ADR-0055 decision 2 on the cascade: a mass_count sweep 20/25/30/35/40 with no overrides, nothing identified.
  // The pool must come from the METRIC's winner, so the fixture makes that winner something other than the
  // first variant: every CE is fed the unfragmented CE-0 spectrum except 30, which gets the rich scan57.
  // mass_count crowns 30; a "first staged scan" rule would put CE 20 on every MS3's stage 0 instead.
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto rich = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  auto poor = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_ce0.txt");
  ABORT_IF(ms1_scans.empty() || rich.empty() || poor.empty())
  std::map<int, ScanData> ce_map{{20, poor[0]}, {25, poor[0]}, {30, rich[0]}, {35, poor[0]}, {40, poor[0]}};

  std::string cfg_str = unidentifiedExhaustiveCytc(/*ms2_exploration=*/true);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));

  const std::map<std::string, std::vector<ScanData>> no_ms3_fixtures;
  const int budget = 256 + 64 * static_cast<int>(ms1_scans.size() + 1);
  AcqResult a = runInterleaved(ida, ms1_scans, std::vector<ScanData>{rich[0]}, &no_ms3_fixtures, budget,
                               /*single_group_only=*/true, 0.0, &ce_map);
  delete ida;

  TEST_EQUAL(a.ms3_cmds.empty(), false)
  for (const ScanCommand& c : a.ms3_cmds)
  {
    TEST_EQUAL(isNoIonDescription(c.scan_description), true)
    TEST_REAL_SIMILAR(c.stages[0].collision_energy, 30.0)
  }
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

  // Feed the first variant (cmds[0], variant_index 0, CE 20); there is no baseline (ADR-0048).
  DeconvolvedSpectrum ds = makeSyntheticDeconv(1, 3);
  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  ExplorationTestAccess::feedResult(exploration,tracking_id, ds, 1.0, queue);

  TEST_EQUAL(exploration.activeGroupCount(), 1)

  auto group = ExplorationTestAccess::group(exploration,1);
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
  // Levels 2 and 3 are no longer selected independently: applyCharacterizationMode_ projects BOTH
  // from characterization.mode, and this fixture is mode:"off", so both are None. That is the
  // point of ADR-0013 -- MS3 needs level 2 AND level 3 non-None (FLASHIda.cpp:368,
  // Exploration.cpp:903 and :905), so driving them from one enum makes "MS3 on with MS2 off"
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
  // 5 CE variants and no reference scan (ADR-0048): the share needs none.
  Config cfg{std::string(remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)

  // The first command is the first GRID point, not a CE-0 reference.
  TEST_REAL_SIMILAR(cmds[0].stages[0].collision_energy, 20.0)

  auto group = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(static_cast<int>(group.exploration_metric),
             static_cast<int>(ExplorationMetric::RemainingPrecursor))
  for (int i = 0; i < 5; ++i) TEST_EQUAL(group.variants[i].variant_index, i)

  // Feed all 5 variants via the deconvolution bypass (no raw arrays -> no share -> not scored).
  for (int i = 0; i < 5; ++i)
  {
    DeconvolvedSpectrum ds = makeSyntheticDeconv(i + 1, i + 1);
    int tracking_id = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    auto info = ExplorationTestAccess::feedResult(exploration,tracking_id, ds, static_cast<double>(i), queue);

    TEST_EQUAL(info.group.total_variants, 5)
    TEST_REAL_SIMILAR(info.metric.precursor_share, -1.0)   // no arrays reached the engine
    TEST_REAL_SIMILAR(info.metric.score, -1.0)             // a share of -1 is "not scored", never 0
  }

  // Group completes with no winner (nothing scored) and is erased.
  TEST_EQUAL(exploration.activeGroupCount(), 0)
}
END_SECTION

START_SECTION(remaining_precursor_score_with_raw_data)
{
  // RemainingPrecursor scores the precursor SHARE of each pre-scan (ADR-0048): the anchor window's raw
  // intensity over the raw intensity of the same scan. No reference scan exists any more.
  Config cfg{std::string(remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)

  auto group = ExplorationTestAccess::group(exploration,1);
  // Exact, not `> 0.0`: initiate() sets precursor_mz = midpoint of pg.getMzRange(charge), and the
  // synthetic group holds ONE peak at 800.0 at charge 3, so the midpoint is that peak (800+800)/2.
  TEST_REAL_SIMILAR(group.precursor_mz, 800.0)

  // At MS2 there is NO 2.0 Th floor (that arm is msn_level >= 3 only) and the synthetic group is a
  // single peak, so mz2 - mz1 == 0 and the window is the ADR-0026 decision-2 margin alone:
  // 2 * optimal_window_margin_ == 0.8 Th, i.e. [799.6, 800.4]. KineticModel::precursorShare tests
  // `>= mz_low && <= mz_high`, so that window admits the 800.0 peak and rejects 790/810/900.
  // Share = 200 / (50 + 200 + 100 + 50) = 0.5; target 0.1 -> score = 1 - 0.4 = 0.6.
  std::vector<double> frag_mzs = {790.0, 800.0, 810.0, 900.0};
  std::vector<double> frag_ints = {50.0, 200.0, 100.0, 50.0};
  int ce20_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto ce20_info = exploration.feedResult(ce20_tid, frag_mzs.data(), frag_ints.data(),
                                          static_cast<int>(frag_mzs.size()), 1.0, queue);

  auto group_mid = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group_mid.variants[0].received, true)
  TEST_REAL_SIMILAR(group_mid.variants[0].precursor_share, 0.5)
  TEST_REAL_SIMILAR(group_mid.variants[0].score, 0.6)
  TEST_REAL_SIMILAR(ce20_info.metric.precursor_share, 0.5)
  TEST_REAL_SIMILAR(ce20_info.metric.score, 0.6)

  // A second variant whose window carries a different fraction: 450 / 900 = 0.5 as well -- the share
  // is a RATIO within the scan, so a scan twice as bright reads the same.
  std::vector<double> bright_ints = {100.0, 450.0, 250.0, 100.0};
  int ce25_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto ce25_info = exploration.feedResult(ce25_tid, frag_mzs.data(), bright_ints.data(),
                                          static_cast<int>(frag_mzs.size()), 1.5, queue);
  TEST_REAL_SIMILAR(ce25_info.metric.precursor_share, 0.5)
  TEST_REAL_SIMILAR(ce25_info.metric.score, 0.6)
}
END_SECTION

START_SECTION(remaining_precursor_score_no_signal_in_window)
{
  // Two distinct "nothing in the window" cases, and they must stay distinct (ADR-0048 d.5, Review
  // Focus 2): a scan with signal elsewhere but none in the window has share 0 -- a fully depleted
  // precursor, which is CLOSER to a 0.1 target than an intact one and scores 0.9 -- while a scan with
  // no signal at all has NO share (-1) and is not scored (-1, so it can never win).
  Config cfg{std::string(remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)

  // Raw data entirely OUTSIDE the precursor window, total > 0 -> share 0 -> score 1 - |0 - 0.1| = 0.9.
  std::vector<double> mzs = {400.0, 500.0, 600.0, 1200.0};
  std::vector<double> intensities = {100.0, 200.0, 300.0, 400.0};
  int ce20_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto ce20_info = exploration.feedResult(ce20_tid, mzs.data(), intensities.data(),
                                          static_cast<int>(mzs.size()), 0.5, queue);
  TEST_REAL_SIMILAR(ce20_info.metric.precursor_share, 0.0)
  TEST_REAL_SIMILAR(ce20_info.metric.score, 0.9)

  // A scan with NO signal anywhere (an empty trap read): no share, not scored. Acquired and received
  // all the same -- nothing is cancelled.
  std::vector<double> zeros = {0.0, 0.0, 0.0, 0.0};
  int ce25_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto ce25_info = exploration.feedResult(ce25_tid, mzs.data(), zeros.data(),
                                          static_cast<int>(mzs.size()), 1.0, queue);
  auto group_after = ExplorationTestAccess::group(exploration,1);
  TEST_EQUAL(group_after.variants[1].received, true)
  TEST_REAL_SIMILAR(ce25_info.metric.precursor_share, -1.0)
  TEST_REAL_SIMILAR(ce25_info.metric.score, -1.0)
  TEST_REAL_SIMILAR(group_after.variants[1].score, -1.0)
}
END_SECTION

START_SECTION(remaining_precursor_late_child_after_group_cleanup_is_noop)
{
  // A result arriving for a variant whose group has already been cleaned up must be a harmless
  // no-op and must not resurrect the group.
  //
  // RETARGETED, not weakened. This hazard used to be reached through the empty-baseline abort,
  // which cancelled the in-flight child and tore the group down while it was still outstanding.
  // Nothing is cancelled any more (ADR-0029), so a group is only ever cleaned up by COMPLETING --
  // and the late/duplicate result the device can still deliver afterwards is exactly the same
  // hazard, guarded by the same routing-map erase at Exploration.cpp:457.
  Config cfg{std::string(remaining_precursor_config)};
  ScanCommandQueue queue(cfg);
  Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)   // no baseline (ADR-0048)
  TEST_EQUAL(exploration.activeGroupCount(), 1)

  // Every child with in-window signal -> the group completes and is erased.
  std::vector<double> base_mzs  = {790.0, 800.0, 810.0, 900.0};
  std::vector<double> frag_ints = {100.0, 500.0, 200.0, 50.0};
  int replay_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  for (size_t i = 0; i < cmds.size(); ++i)
  {
    int tid = queue.decode(std::string(cmds[i].scan_description).substr(0, 3));
    exploration.feedResult(tid, base_mzs.data(), frag_ints.data(),
                           static_cast<int>(base_mzs.size()), 1.0, queue);
  }
  TEST_EQUAL(exploration.activeGroupCount(), 0)   // completed + cleaned up

  // The device now delivers a LATE duplicate for a child of that erased group -> harmless no-op.
  auto late_info = exploration.feedResult(replay_tid, base_mzs.data(), frag_ints.data(),
                                          static_cast<int>(base_mzs.size()), 1.0, queue);
  TEST_EQUAL(late_info.group.group_id, -1)        // empty info (routing entry gone)
  TEST_EQUAL(late_info.group.variant_index, -1)   // not routed to any variant
  TEST_EQUAL(exploration.activeGroupCount(), 0)   // group NOT resurrected
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
// used to make the sweep cheap by narrowing its pre-scans to the window it reads, and refused an
// empty overrides map as that narrowing's interlock. ADR-0044 withdrew both (the narrowing bought
// nothing on the instrument), so the config surface is now:
//
//   metric "remaining_precursor" with an EMPTY overrides map -- LOADS. Its pre-scans run at
//                production settings; at MS3 ADR-0020 gate #2 still re-acquires the winner.
//   metric "remaining_precursor" at a level whose OWN charge mode is "multiplexed" -- THROWS,
//                because the ratio sums the anchor's window only (ADR-0044 decision 4). Two
//                pair-specific clauses, not one "multiplexed anywhere" test.
//
// One section pins the load, two pin the throws, and the fourth pins what must keep LOADING under
// the multiplexed rule, because a rejection written one clause too broad is invisible to every
// positive test.
START_SECTION(remaining_precursor_loads_without_overrides)
{
  // FAILS WHEN ADR-0026 decision 3 returns: the metric with an empty overrides map was rejected at
  // load, at both levels, and the rejection was the interlock for a scan-range binding that no
  // longer exists. Level-agnostic on purpose, like the throw was -- a rule resurrected "where it was
  // reported" would trip one of the two lines below and not the other.
  Config ms2{sweepProbe(", " + rpSweep(""))};
  TEST_EQUAL(ms2.level(2).overrides.empty(), true)
  TEST_EQUAL(static_cast<int>(ms2.level(2).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))

  Config ms3{sweepProbe("", ms3_on + rpSweep(""), true)};
  TEST_EQUAL(ms3.level(3).overrides.empty(), true)
  TEST_EQUAL(static_cast<int>(ms3.level(3).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))

  // And the non-empty shape still loads too -- overrides are optional, not forbidden.
  Config ms2_ov{sweepProbe(", " + rpSweep(sweep_overrides))};
  TEST_EQUAL(ms2_ov.level(2).overrides.empty(), false)
}
END_SECTION

START_SECTION(remaining_precursor_multiplexed_ms2_throws)
{
  // The multiplexed rejection at level 2 (ADR-0026 decision 4, re-grounded by ADR-0044). The
  // remaining-precursor ratio sums the ANCHOR's isolation window and nothing else, so under
  // co-isolation it reports one charge state's depletion while the same collision energy depletes
  // its siblings at other rates -- a target met by the anchor decides nothing about the rest. The
  // engine already RECORDS that anchor-only number for every multiplexed variant; this stops it
  // being DECIDED on.
  //
  // The probe carries a non-empty overrides map only because rpSweep takes one; overrides are
  // optional now (ADR-0044), so the map has no bearing on which rule fires.
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
  // The multiplexed rejection at level 3 -- a separate clause from level 2's, deliberately, because
  // the two cross-level pairs stay legal (see the negative section below). An MS3 sweep co-isolates
  // when characterization.fragment_charges multiplexes the SUB-fragments it is reading, which is
  // the level-3 mirror of the MS2 case and is not something the level-2 clause can see.
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
  // scan carries no notches at all and the anchor window IS the whole isolation -- the ratio reads
  // everything that was isolated. Rejecting `separate` would cost the mode for no reason.
  //
  // (c) and (d): the CROSS-LEVEL pairs, where the multiplexed level is not the level that sweeps.
  // Stage-0 notches change WHICH precursors are fragmented, not where a stage-1 readout sits, so an
  // MS3 sub-fragment scan under a multiplexed MS2 still isolates one window; and a multiplexed
  // MS3 does not touch the MS2 pre-scan above it. Only the sweeping level's OWN charge mode can
  // put more in the isolation than the anchor window the ratio reads.
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
// THE SCAN RANGE OF A PRE-SCAN (ADR-0044, withdrawing ADR-0026 decision 1). A pre-scan has no range
// of its own: it reads out its level's configured first_mass/last_mass, patched by overrides like any
// other key -- at MS2 and MS3 alike, for EVERY variant including the CE-0 baseline, and for the
// ADR-0023 FORCED metric that no config declares. ADR-0026 used to bind a remaining_precursor
// pre-scan to precursor_mz -/+ isolation_width/2 instead; measured on the instrument, the narrowing
// bought nothing, and it discarded the rest of every spectrum.
//
// Every range assertion below is read off the PARSED CONFIG (`cfg.level(n).scans[0].first_mass`),
// never restated as a literal, so a fixture typo cannot make a section pass; and every fixture names
// a NON-ZERO range, so "the configured range" and "a ~2 Th window" cannot coincide. The
// isolation_width pins survive from the old sections: ADR-0026 decision 2 -- the margin at MS2, the
// 2.0 Th floor at MS3 -- still governs what the metric SUMS, only no longer what the scan READS.
/////////////////////////////////////////////////////////////

START_SECTION(pre_scans_read_out_the_levels_scan_range_ms2)
{
  // FAILS WHEN a per-target scan-range binding returns at level 2 (ADR-0026 decision 1, withdrawn by
  // ADR-0044): every command carries ms_settings.ms2's 300/1800 like any other scan built from
  // scans[0] + overrides, the CE-0 baseline included. 300/1800 rather than the 0 default so "the
  // configured range" and "a ~2 Th window" cannot coincide.
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

  // The control the loop below is read against, asserted off the parsed config rather than restated
  // as a number, so a fixture typo cannot make the section pass.
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].first_mass, 300.0)
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].last_mass, 1800.0)

  PeakGroup pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  std::vector<ScanCommand> cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 3)   // CE 20 / 25 / 30; no baseline (ADR-0048)
  ABORT_IF(cmds.size() != 3)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(group.msn_level, 2)
  TEST_EQUAL(static_cast<int>(group.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))
  TEST_EQUAL(group.variants[0].variant_index, 0)
  TEST_REAL_SIMILAR(group.variants[0].collision_energy, 20.0)
  TEST_EQUAL(group.variants[0].cmd.scan_id, cmds[0].scan_id)

  // What a binding WOULD write, so the loop is a contrast and not a tautology: ~800 Th against 300.
  const double bound_lo = group.precursor_mz - group.isolation_width / 2.0;
  TEST_EQUAL(std::abs(bound_lo - cfg.level(2).scans[0].first_mass) > 100.0, true)

  for (const ScanCommand& c : cmds)
  {
    TEST_EQUAL(c.msn_level, 2)
    TEST_REAL_SIMILAR(c.first_mass, cfg.level(2).scans[0].first_mass)
    TEST_REAL_SIMILAR(c.last_mass, cfg.level(2).scans[0].last_mass)
    TEST_STRING_EQUAL(c.analyzer, "IonTrap")   // the override reached the command
  }
}
END_SECTION

START_SECTION(pre_scans_read_out_the_levels_scan_range_ms3)
{
  // FAILS WHEN a per-target scan-range binding returns at level 3 -- the MS3 twin of the section
  // above, because buildMS3 is a different builder and an edit could narrow one level and not the
  // other. ms_settings.ms3 names 350/1750 as the non-zero control for the same reason ms2 named
  // 300/1800.
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
  TEST_EQUAL(static_cast<int>(cmds.size()), 3)   // CE 20 / 25 / 30; no baseline (ADR-0048)
  ABORT_IF(cmds.size() != 3)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(group.msn_level, 3)
  TEST_EQUAL(static_cast<int>(group.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))
  TEST_EQUAL(group.variants[0].variant_index, 0)
  TEST_EQUAL(group.variants[0].cmd.scan_id, cmds[0].scan_id)

  // ADR-0026 decision 2 STANDS: at MS3 isolation_width is std::max(mz2 - mz1, 2.0), and the fixture
  // is a single peak at 500.0, so the floor IS the window the metric sums. Fails if the floor is
  // dropped, re-tuned, or if the MS2 margin branch is taken at level 3.
  TEST_REAL_SIMILAR(group.isolation_width, 2.0)

  // The window the metric sums is the FRAGMENT stage's, not the inherited MS2 precursor's. An MS3
  // command carries both (stage[0] ~ 800 Th from ms2_ctx, stage[1] = the ~500 Th fragment), and
  // precursorWindowIntensity_ sums around group.precursor_mz.
  TEST_REAL_SIMILAR(cmds[0].stages[1].precursor_mz, group.precursor_mz)
  TEST_EQUAL(std::abs(cmds[0].stages[0].precursor_mz - group.precursor_mz) > 100.0, true)

  for (const ScanCommand& c : cmds)
  {
    TEST_EQUAL(c.msn_level, 3)
    TEST_EQUAL(c.num_stages, 2)
    TEST_REAL_SIMILAR(c.first_mass, cfg.level(3).scans[0].first_mass)
    TEST_REAL_SIMILAR(c.last_mass, cfg.level(3).scans[0].last_mass)
    TEST_STRING_EQUAL(c.analyzer, "IonTrap")
  }
}
END_SECTION

START_SECTION(production_scan_is_built_from_the_unoverridden_config)
{
  // FAILS WHEN the post-winner close-out is built from base_config (scans[0] + overrides) instead of
  // from level_config.scans[0]. The production scan is the ONE acquisition of the group that is meant
  // to be identified, so it must carry the level's own analyzer and range while the pre-scans that
  // chose its CE ran under the overrides -- here on the ION TRAP (ADR-0045). That separation is
  // structural, no flag defends it, and this is the ctest that pins the analyzer on both sides
  // (scan_commands.tsv has no analyzer column, so no golden can).
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
  TEST_EQUAL(static_cast<int>(cmds.size()), 3)   // CE 20 / 25 / 30; no baseline (ADR-0048)
  ABORT_IF(cmds.size() != 3)

  // Read the group BEFORE the sweep completes -- completion erases it from active_groups_.
  auto group = ExplorationTestAccess::group(exploration, 1);
  // The pre-scans ran under the override (ion trap) at the level's full range -- the contrast the
  // production-scan assertions at the end are read against.
  for (const ScanCommand& c : cmds)
  {
    TEST_STRING_EQUAL(c.analyzer, "IonTrap")
    TEST_REAL_SIMILAR(c.first_mass, cfg.level(3).scans[0].first_mass)
    TEST_REAL_SIMILAR(c.last_mass, cfg.level(3).scans[0].last_mass)
  }

  // Drive the sweep to completion with RAW ARRAYS. RemainingPrecursor scores the precursor SHARE --
  // the centre peak over the whole scan (ADR-0048) -- so one peak at the centre and two 50 Th away
  // make the share exactly controllable; the ExplorationTestAccess deconvolution bypass the mass_count
  // sections use cannot drive this metric at all.
  const double mz_c = group.precursor_mz;
  std::vector<double> mzs{mz_c - 50.0, mz_c, mz_c + 50.0};
  // score = 1 - |share - remaining_precursor_target|, target = 0.1 (the level default). Shares
  // 900/1900 = 0.474, 400/1400 = 0.286, 100/1100 = 0.0909 give 0.626 / 0.814 / 0.991 and the LAST-fed
  // variant wins. A tie-break-order winner would land on cmds[0] and the CE assertion below would catch it.
  std::vector<double> ce20_ints{300.0, 900.0, 700.0};
  std::vector<double> ce25_ints{300.0, 400.0, 700.0};
  std::vector<double> ce30_ints{300.0, 100.0, 700.0};

  exploration.feedResult(cmds[0].scan_id, mzs.data(), ce20_ints.data(),
                         static_cast<int>(mzs.size()), 1.0, queue);
  exploration.feedResult(cmds[1].scan_id, mzs.data(), ce25_ints.data(),
                         static_cast<int>(mzs.size()), 1.5, queue);
  auto info = exploration.feedResult(cmds[2].scan_id, mzs.data(), ce30_ints.data(),
                                     static_cast<int>(mzs.size()), 2.0, queue);

  TEST_EQUAL(exploration.activeGroupCount(), 0)
  TEST_REAL_SIMILAR(info.metric.precursor_share, 100.0 / 1100.0)
  TEST_REAL_SIMILAR(info.metric.score, 1.0 - std::abs(100.0 / 1100.0 - 0.1))
  TEST_EQUAL(info.group.winner_tracking_id, ScanCommandQueue::encode(cmds[2].scan_id))

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

  // THE ASSERTIONS THIS SECTION EXISTS FOR: the un-overridden config's analyzer and range. Stated
  // twice on purpose -- against the literals so the intent is readable, and against the parsed config
  // so the pair cannot both drift to some third value.
  TEST_STRING_EQUAL(info.commands[0].analyzer, "Orbitrap")
  TEST_REAL_SIMILAR(info.commands[0].first_mass, 200.0)
  TEST_REAL_SIMILAR(info.commands[0].last_mass, 2000.0)
  TEST_REAL_SIMILAR(info.commands[0].first_mass, cfg.level(3).scans[0].first_mass)
  TEST_REAL_SIMILAR(info.commands[0].last_mass, cfg.level(3).scans[0].last_mass)
}
END_SECTION

START_SECTION(range_overrides_patch_pre_scans_like_any_key)
{
  // ADR-0044 decision 2: with no per-target binding there is nothing for a range override to
  // suppress -- first_mass/last_mass in exploration.overrides patch the pre-scan exactly as
  // `analyzer` or `resolution` do. This is how a STATIC narrow range stays expressible; what went
  // away is the dynamic, per-target one. FAILS WHEN applyOverrides loses its two range branches.
  //
  // ms_settings.ms2 carries its OWN 100/2000, deliberately different from the override, so the
  // 600/1600 below is attributable to the overrides map: set ms2 to 600/1600 and the assertions
  // would pass even if applyOverrides never ran at all.
  //
  // The override VALUES ARE JSON STRINGS. The overrides loop reads the map with `.get<std::string>()`,
  // so a bare `"first_mass": 600` throws a nlohmann type error out of the constructor and the
  // section fails as an unloadable fixture rather than as an override defect.
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

  // The map survived parsing under both keys. tolerance_ppm was erased from this map once; that is
  // now only history in a comment and the live code throws on the key instead, so "an exploration
  // key reached overrides" is not something to take on trust.
  TEST_EQUAL(cfg.level(2).overrides.count("first_mass"), 1u)
  TEST_EQUAL(cfg.level(2).overrides.count("last_mass"), 1u)
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].first_mass, 100.0)   // the AUTHORED scan range, untouched
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].last_mass, 2000.0)

  PeakGroup pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  std::vector<ScanCommand> cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 3)   // no baseline (ADR-0048)
  ABORT_IF(cmds.size() != 3)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(static_cast<int>(group.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))

  // Neither the authored 100/2000 nor a ~800 Th window: the override, on every command.
  const double window_lo = group.precursor_mz - group.isolation_width / 2.0;
  TEST_EQUAL(std::abs(window_lo - 600.0) > 1.0, true)

  for (const ScanCommand& c : cmds)
  {
    TEST_REAL_SIMILAR(c.first_mass, 600.0)
    TEST_REAL_SIMILAR(c.last_mass, 1600.0)
  }
}
END_SECTION

START_SECTION(forced_remaining_precursor_reads_the_configured_range)
{
  // ADR-0023 decision 11 DRAGS an exhaustive-mode unassigned mass -- ion class 'u', which
  // MS3FragmentMatcher::isKnownIonClass rejects -- onto RemainingPrecursor whatever the config asked
  // for. This config asks for fragment_count. Under ADR-0026 the forced sweep was the one path that
  // was narrowed with no config entry to say so; under ADR-0044 it reads out ms_settings.ms3's
  // 200/2000 like every other pre-scan. FAILS WHEN a binding keyed on group.exploration_metric
  // returns. Mirrors ProteoformTracker_Exhaustive_test.cpp, which pins the force itself.
  //
  // The override is the ORBITRAP, deliberately: ADR-0045 refuses an ion-trap override under a
  // counting metric at load, and the CONFIGURED metric here is fragment_count. (The forced sweep would
  // be compatible with the trap -- it forces the one metric that can score a trap pre-scan -- but the
  // load rule reads the config, not the force, and that is the right side to err on.)
  const std::string cfg_json = R"JSON({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": { "rank_by": "qscore", "max_precursors": 3 },
    "characterization": {
      "mode": "exhaustive", "protein_sequence": "PEPTIDEK", "max_targets": 3,
      "exploration": { "metric": "fragment_count", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0,
                       "overrides": { "analyzer": "Orbitrap" } }
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
  TEST_REAL_SIMILAR(cfg.level(3).scans[0].first_mass, 200.0)   // the range every pre-scan must carry
  TEST_REAL_SIMILAR(cfg.level(3).scans[0].last_mass, 2000.0)

  PeakGroup unassigned_pg = makeSyntheticPeakGroup(500.0, 1000.0, 2);
  ScanCommand ms2_ctx = queue.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3,
                                       cfg.level(2).scans[0], 2, 0);

  // 'u'/0 is exactly what planExhaustive_ labels a mass that matched no theoretical fragment.
  std::vector<ScanCommand> cmds = exploration.initiate(3, unassigned_pg, 2, queue, nullptr, &ms2_ctx, 'u', 0);
  TEST_EQUAL(static_cast<int>(cmds.size()), 3)   // CE 20 / 25 / 30; no baseline (ADR-0048)
  ABORT_IF(cmds.size() != 3)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(group.msn_level, 3)
  TEST_EQUAL(static_cast<int>(group.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))
  TEST_EQUAL(group.variants[0].cmd.scan_id, cmds[0].scan_id)

  const double lo = group.precursor_mz - group.isolation_width / 2.0;
  const double hi = group.precursor_mz + group.isolation_width / 2.0;
  // The FORCED sweep takes the same MS3 floor a configured one does, so pin the width: 2.0 also says
  // it was computed from the FRAGMENT group and not from ms2_ctx, whose buildMS2 window is 0.8 Th.
  TEST_REAL_SIMILAR(group.isolation_width, 2.0)
  TEST_EQUAL(lo > cfg.level(3).scans[0].first_mass, true)   // the window sits strictly inside the range,
  TEST_EQUAL(hi < cfg.level(3).scans[0].last_mass, true)    // so a returning binding cannot pass the loop

  for (const ScanCommand& c : cmds)
  {
    TEST_EQUAL(c.msn_level, 3)
    TEST_REAL_SIMILAR(c.first_mass, cfg.level(3).scans[0].first_mass)
    TEST_REAL_SIMILAR(c.last_mass, cfg.level(3).scans[0].last_mass)
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
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)   // no baseline (ADR-0048)

  auto group = ExplorationTestAccess::group(exploration,1);
  double mz_center = group.precursor_mz;

  // Two peaks: the anchor at the window centre, the rest of the scan 100 Th away. The share is the
  // anchor's fraction of the two (ADR-0048), so 100 / (100 + 900) = 0.1 hits the target exactly.
  std::vector<double> two_mzs = {mz_center, mz_center + 100.0};
  std::vector<double> perfect_ints = {100.0, 900.0};
  int ce20_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info_perfect = exploration.feedResult(ce20_tid, two_mzs.data(), perfect_ints.data(),
                                              2, 1.0, queue);
  double score_perfect = info_perfect.metric.score;
  double share_perfect = info_perfect.metric.precursor_share;

  // 500 / (500 + 500) = 0.5 -> |0.5 - 0.1| = 0.4 -> score 0.6.
  std::vector<double> over_ints = {500.0, 500.0};
  int ce25_tid = queue.decode(std::string(cmds[1].scan_description).substr(0, 3));
  auto info_over = exploration.feedResult(ce25_tid, two_mzs.data(), over_ints.data(),
                                           2, 2.0, queue);
  double score_over = info_over.metric.score;
  double share_over = info_over.metric.precursor_share;

  // (The former `score_perfect > score_over` ordering check is gone: the two exact values below
  // already imply it, and an ordering assertion that survives a shift of BOTH scores is the weaker
  // of the two claims, never an independent one.)
  TEST_REAL_SIMILAR(score_perfect, 1.0)
  TEST_REAL_SIMILAR(score_over, 0.6)
  TEST_REAL_SIMILAR(share_perfect, 0.1)
  TEST_REAL_SIMILAR(share_over, 0.5)
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
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)   // no baseline (ADR-0048)

  DeconvolvedSpectrum ds = makeSyntheticDeconv(1, 5);
  int tracking_id = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
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
  // RemainingPrecursor: 5 CE variants (20,25,30,35,40); no baseline (ADR-0048)
  TEST_EQUAL(static_cast<int>(cmds.size()), 5)

  auto group = ExplorationTestAccess::group(exploration,1);
  // isolation_width should be floored to 2.0 (not 0.0)
  TEST_REAL_SIMILAR(group.isolation_width, 2.0)

  double mz_center = group.precursor_mz;

  // With the 2.0 Da window [499.0, 501.0], mz_center=500.0 is in-window and a peak 50 Th away is not.
  // Share = 100 / (100 + 900) = 0.1 -> target 0.1, deviation 0, score 1.0. Under a zero width the
  // window would admit nothing, the share would read 0 and the score 0.9.
  std::vector<double> mzs = {mz_center, mz_center + 50.0};
  std::vector<double> variant_ints = {100.0, 900.0};
  int ce20_tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info = exploration.feedResult(ce20_tid, mzs.data(), variant_ints.data(), 2, 1.0, queue);

  TEST_REAL_SIMILAR(info.metric.precursor_share, 0.1)
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
// The correction is Exploration.cpp:166-168's MS2 arm, `(mz2 - mz1) + 2.0 * optimal_window_margin_`,
// and the three sections below are the whole case for shipping it in the SAME push as the ADR-0026
// binding rather than a later one. ALL THREE go red on a revert -- there is no one guard section to
// point at -- but each fails in a different quantity, and only the third makes "later" impossible:
//   * ms2_scoring_window_matches_commanded_isolation -- the only one that fails on a SUM
//     (baseline_intensity 1500 -> 1000) rather than on a width, because it is the only one whose
//     fixture puts signal in the 0.4 Th flanks. So it is the only place here that shows the margin
//     changing what a metric MEASURES, and not merely what the command carries.
//   * optimal_window_margin_has_one_definition -- the drift guard for the constant the correction
//     had to be collapsed onto (three TU-local copies -> NotchSelection.h:129). Its widths go
//     3.8 -> 3.0 on a revert.
//   * single_isotope_charge_yields_non_degenerate_scan_range -- the defect that makes "later"
//     impossible: without the margin the binding emits a zero-width scan range (0.8 -> 0.0).
/////////////////////////////////////////////////////////////

START_SECTION(ms2_scoring_window_matches_commanded_isolation)
{
  // FAILS WHEN Exploration.cpp:168 loses `+ 2.0 * optimal_window_margin_` and reverts to the bare
  // `(mz2 - mz1)`. The two sections after this one go red on that same edit, both of them on a WIDTH;
  // this is the only one that goes red on a SUM, which is the only way the file shows that the margin
  // changes what a metric MEASURES and not just what the command carries. Every other
  // remaining_precursor section feeds peaks that are 10 Th or more off centre (790/810/900) or sits at
  // MS3, where std::max(..., 2.0) swallows a 0.8 Th difference whole. This fixture is built the other
  // way round -- two of its four baseline peaks land INSIDE the commanded window and OUTSIDE the
  // un-margined one, so the two definitions of "the window" give two different sums.
  //
  // THE ARITHMETIC, derived from the fixture and from precursorWindowIntensity_'s actual predicate
  // (`mzs[i] >= mz_low && mzs[i] <= mz_high`, Exploration.cpp:1228), not guessed:
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
  TEST_EQUAL(static_cast<int>(cmds.size()), 3)   // CE 20 / 25 / 30; no baseline (ADR-0048)
  ABORT_IF(cmds.size() != 3)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(group.msn_level, 2)
  TEST_REAL_SIMILAR(group.precursor_mz, 801.0)

  // THE INVARIANT, stated rather than inferred from the sum below: an MS2 exploration window is the
  // measured span plus the margin on each side (2.8 Th here). Expressed through the shared symbol, so
  // it says "the same margin buildMS2 uses" and not "0.8" -- a retuned constant should move this
  // section through optimal_window_margin_has_one_definition, not silently agree with a stale literal.
  TEST_REAL_SIMILAR(group.isolation_width, (802.0 - 800.0) + 2.0 * optimal_window_margin_)

  // ...and the COMMANDED side of the same claim, read off the first command the group actually
  // emitted. This is ADR-0026 decision 2 in one line: the scored interval IS the isolated interval. It
  // is a second, independent way for a reverted MS2 arm to go red (group 2.0 against command 2.8), and
  // unlike the sum below it needs no fed peaks at all.
  TEST_REAL_SIMILAR(cmds[0].stages[0].precursor_mz, group.precursor_mz)
  TEST_REAL_SIMILAR(cmds[0].stages[0].isolation_width, group.isolation_width)

  std::vector<double> peak_mzs  = {799.8, 801.0, 802.2, 810.0};
  std::vector<double> peak_ints = {300.0, 1000.0, 200.0, 700.0};
  auto info = exploration.feedResult(cmds[0].scan_id, peak_mzs.data(), peak_ints.data(),
                                     static_cast<int>(peak_mzs.size()), 0.5, queue);

  // THE ASSERTION THIS SECTION EXISTS FOR. The precursor SHARE (ADR-0048) is the window sum over the
  // whole scan's 2200: 1500/2200 with the margin, 1000/2200 without it -- so a revert does not merely
  // widen a tolerance, it changes the number.
  TEST_REAL_SIMILAR(info.metric.precursor_share, 1500.0 / 2200.0)
  TEST_EQUAL(std::abs(info.metric.precursor_share - 1000.0 / 2200.0) > 0.1, true)
}
END_SECTION

START_SECTION(optimal_window_margin_has_one_definition)
{
  // DRIFT GUARD for the constant this push collapsed from three TU-local copies -- ScanCommandQueue.cpp,
  // PrecursorSelection.cpp, and FragmentAnalysis.cpp's anonymous namespace, all `.4`, none of them aware
  // of the others -- onto the single inline constexpr at NotchSelection.h:129. Collapsing two
  // definitions into one earns a permanent bidirectional guard; this is it.
  //
  // WHAT IT COVERS, WHICH IS NOT ALL THREE TUs. The cross-TU half below calls queue.buildMS2, so
  // ScanCommandQueue is the only one of the three collapsed copies this section actually executes;
  // Exploration -- which held no copy and was handed the symbol by the same decision -- is covered by
  // the initiate() half. PrecursorSelection.cpp and FragmentAnalysis.cpp are NOT: a FragmentAnalysis
  // IS constructed below (Exploration's ctor requires one) but this section calls no method on it
  // that pads a window -- its four margin uses are in getBestMS2Masses, getTopFragmentMatches,
  // getAmbiguityEnclosingIons and getTerminalFragmentIons, none of which initiate() reaches -- and
  // PrecursorSelection is never constructed at all. The only other C++ test that sees their padded
  // windows checks ordering
  // alone (FragmentAnalysis_test's `wstarts[i] < wends[i]`), which any non-zero margin satisfies. So a
  // private copy or a stale hard-coded .4 in those two files stays green everywhere. Widening this
  // guard means driving those two consumers, not adding assertions to the two already covered.
  //
  // FAILS WHEN the shared symbol and a consumer disagree. The literal check below catches a retune of
  // the symbol itself, so an intentional 0.4 -> 0.5 has to come through here and be acknowledged; on
  // its own it cannot see a consumer that stopped reading the symbol. That is the buildMS2 half's job:
  // it compares the width ScanCommandQueue ACTUALLY computed against one derived from the shared
  // symbol, so the two sides have to agree numerically or the section is red. The reachable shapes are
  // a block-scope `const double optimal_window_margin_ = .5;` inside buildMS2, which hides the
  // namespace-scope symbol for that function alone, and a retuned header with a hard-coded 0.4 left
  // behind at the use site. NOT a namespace-scope re-declaration in ScanCommandQueue.cpp: its code
  // sits in `namespace OpenMS`, the same namespace the header declares the symbol in, so re-adding it
  // there is a redefinition -- the TU fails to compile rather than shadowing anything.
  TEST_REAL_SIMILAR(optimal_window_margin_, 0.4)

  // mass_count, deliberately, and not remaining_precursor: this section is about GEOMETRY, and
  // mass_count is untouched by either ADR-0026 config rejection (Config.cpp:963-969 and :985-1006), so the fixture
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
  TEST_EQUAL(static_cast<int>(cmds.size()), 3)   // CE 20 / 25 / 30; no baseline (ADR-0048)
  ABORT_IF(cmds.size() != 3)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_REAL_SIMILAR(group.precursor_mz, 601.5)
  TEST_REAL_SIMILAR(group.isolation_width, expected_width)
  // The two sides against each other, with no literal in between: whatever value the margin takes, the
  // scored window and the commanded window are the same interval.
  TEST_REAL_SIMILAR(group.isolation_width, cmds[0].stages[0].isolation_width)
  TEST_REAL_SIMILAR(group.precursor_mz, cmds[0].stages[0].precursor_mz)
}
END_SECTION

START_SECTION(single_isotope_charge_keeps_the_margin_window)
{
  // ADR-0026 decision 2 STANDS. A charge state the deconvolution resolved at ONE isotope gives
  // mz2 == mz1, so the bare span is 0; the margin turns that into a 0.8 Th window for the metric to
  // sum over, and MS2 has no 2.0 Th floor of its own (that arm is msn_level >= 3 only). Under
  // ADR-0026 this width also became the scan range, and a dropped margin would have sent a zero-width
  // "DefineMZRange" to the instrument; under ADR-0044 nothing turns the width into a range, so the
  // commands carry ms2's unset 0/0 whatever the window is -- asserted below, so a returning binding
  // is caught here as a non-zero range.
  //
  // FAILS WHEN the MS2 margin is dropped (isolation_width collapses to 0) or a binding returns.
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

  // ONE peak, and it stays one peak: a single-isotope charge is the input, not an inconvenience.
  // makeSyntheticPeakGroup builds precisely that, so no helper needed changing.
  PeakGroup pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  const std::tuple<double, double> mz_range = pg.getMzRange(3);
  const double mz1 = std::get<0>(mz_range);
  const double mz2 = std::get<1>(mz_range);
  TEST_REAL_SIMILAR(mz2 - mz1, 0.0)      // the degenerate span, asserted rather than assumed
  TEST_REAL_SIMILAR(mz1, 800.0)          // and not the (-1, -10) no-such-charge sentinel

  // ms_settings.ms2 names NO range, so 0/0 (= unset, the instrument method default) is what every
  // command must carry.
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].first_mass, 0.0)
  TEST_REAL_SIMILAR(cfg.level(2).scans[0].last_mass, 0.0)
  TEST_EQUAL(cfg.level(2).overrides.count("first_mass"), 0u)   // no range override either
  TEST_EQUAL(cfg.level(2).overrides.count("last_mass"), 0u)

  std::vector<ScanCommand> cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL(static_cast<int>(cmds.size()), 3)   // CE 20 / 25 / 30; no baseline (ADR-0048)
  ABORT_IF(cmds.size() != 3)

  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(static_cast<int>(group.exploration_metric), static_cast<int>(ExplorationMetric::RemainingPrecursor))
  TEST_REAL_SIMILAR(group.precursor_mz, 800.0)
  // Span 0 plus the margin on each side. Nothing else contributes at MS2, so this is the floor in its
  // entirety: 0.8 Th, neither 2.0 nor 0.
  TEST_REAL_SIMILAR(group.isolation_width, 2.0 * optimal_window_margin_)

  // EVERY command, the CE-0 baseline included, carries the unset range: the window above is what the
  // metric SUMS, never what the scan READS (ADR-0044).
  for (const ScanCommand& c : cmds)
  {
    TEST_EQUAL(c.msn_level, 2)
    TEST_REAL_SIMILAR(c.first_mass, 0.0)
    TEST_REAL_SIMILAR(c.last_mass, 0.0)
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
  "scheduling": { "agc_interval_seconds": 9999999 },
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
  "scheduling": { "agc_interval_seconds": 9999999 },
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
    ++real_variants;
    TEST_STRING_EQUAL(v.activation_type, "ETD")
    TEST_EQUAL(std::isfinite(v.reaction_time), true)
    seen_rts.insert(v.reaction_time);
  }
  // rt_min=5, rt_max=15, rt_step=5 -> 3 distinct reaction-time variants and nothing else: no
  // reference scan under any metric (ADR-0048).
  TEST_EQUAL(static_cast<int>(group.variants.size()), 3)
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
// group.variants[best_idx].cmd (Exploration.cpp:819-821 -> initiateNextLevel -> buildMS3), which is a
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
// Exploration.cpp:1037-1038 (std::abs(target.frag_charge) < charge_floor). Two sibling tests so the 0 in
// the floor-filters case is NOT vacuous: the no-floor sibling proves the model DOES emit MS3 here.
/////////////////////////////////////////////////////////////

START_SECTION(ms3_min_charge_floor_filters_all_targets)
{
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.size() < 2 || ms2_scans.empty())

  // Clone of inclusion_ms3_config (A5) with selection_strategy.ms2.min_charge = 99 —
  // impossible floor that filters every charge-1-3 MS3 target at Exploration.cpp:1038.
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
  // charge 1-3 which is < 99, so the charge floor at Exploration.cpp:1038 filters them all
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
  // floor at Exploration.cpp:1038, not by a structural dead-end.
  TEST_EQUAL(a.ms2_cmds.size() >= 1, true)
  TEST_EQUAL(a.ms3_cmds.size() >= 1, true)
}
END_SECTION

START_SECTION(a_sweep_always_acquires_both_edges)
{
  // Summary: 15..40/2 does not divide. The accumulating loop stopped at 39 and the top edge was never
  // acquired. Purpose: the decimation study's result rests on "both edges acquired" (ADR-0047 d.8); the
  // production method's MS3 grid is 15..40/2. In practice: no committed golden config has a
  // non-dividing grid, so only etc/method.json's MS3 sweep gains a scan.
  Config cfg{sweepCfgRange(R"("HCD")", 15.0, 40.0, 2.0)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);
  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  auto group = ExplorationTestAccess::group(exploration, 1);
  std::vector<double> ces;
  for (const auto& v : group.variants) ces.push_back(v.collision_energy);
  TEST_EQUAL((int)ces.size(), 14)                 // 15,17,...,39 (13) + 40
  ABORT_IF(ces.size() != 14)
  TEST_REAL_SIMILAR(ces.front(), 15.0)
  TEST_REAL_SIMILAR(ces[12], 39.0)
  TEST_REAL_SIMILAR(ces.back(), 40.0)
  TEST_EQUAL((int)cmds.size(), (int)group.variants.size())

  // A dividing grid is unchanged: 20..40/5 is still exactly 5 points.
  Config cfg2{sweepCfgRange(R"("HCD")", 20.0, 40.0, 5.0)};
  ScanCommandQueue q2(cfg2);
  FragmentAnalysis f2(cfg2);
  Exploration e2(cfg2, f2);
  e2.initiate(2, pg, 3, q2);
  TEST_EQUAL((int)ExplorationTestAccess::group(e2, 1).variants.size(), 5)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0048: a sweep is exactly its grid points -- no reference scan under any metric
/////////////////////////////////////////////////////////////

START_SECTION(a_sweep_acquires_exactly_its_grid_points)
{
  // Summary: CID + HCD + ETD swept together; every activation contributes its grid points and nothing
  // else. Purpose: ADR-0048 d.3 supersedes ADR-0029 -- the precursor share is referenced against the
  // scan's own spectrum, so the per-activation CE-0 / RT-floor reference scan is gone, under every
  // metric. In practice: a sweep costs one scan less per swept activation.
  Config cfg{sweepCfg(R"("CID", "HCD", "ETD")", 15.0, 5.0)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  auto group = ExplorationTestAccess::group(exploration, 1);
  const auto& vs = group.variants;

  // CID 15/25/35 + HCD 15/25/35 + ETD RT 5/15/25: 9 variants, 9 commands, no head-of-block reference.
  TEST_EQUAL(static_cast<int>(vs.size()), 9)
  TEST_EQUAL(static_cast<int>(cmds.size()), 9)
  TEST_EQUAL(blockStart(vs, "CID"), 0)
  TEST_EQUAL(blockStart(vs, "HCD"), 3)
  TEST_EQUAL(blockStart(vs, "ETD"), 6)
  for (int i = 0; i < 9; ++i) TEST_EQUAL(vs[i].variant_index, i)      // one running counter, no -1
  // No CE-0 and no RT-floor point anywhere: the grid starts at 15 / 5.
  for (const auto& v : vs)
  {
    if (v.activation_type == "ETD") TEST_EQUAL(v.reaction_time >= 5.0, true)
    else TEST_EQUAL(v.collision_energy >= 15.0, true)
  }
}
END_SECTION

START_SECTION(an_authored_ce_zero_is_an_ordinary_pre_scan)
{
  // A grid authored to start at CE 0 acquires that point as an ordinary first variant (variant_index 0),
  // scored like any other; nothing is flagged, nothing is skipped (ADR-0048 d.3).
  Config cfg{sweepCfg(R"("HCD")", 0.0, 5.0)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(static_cast<int>(group.variants.size()), 3)             // 0 / 10 / 20
  TEST_EQUAL(static_cast<int>(cmds.size()), 3)
  TEST_REAL_SIMILAR(group.variants[0].collision_energy, 0.0)
  TEST_EQUAL(group.variants[0].variant_index, 0)
  // Fed a scorable result it competes: mass_count scores spec.size(), so it wins outright here.
  DeconvolvedSpectrum ds = makeSyntheticDeconv(1, 4);
  int tid = queue.decode(std::string(cmds[0].scan_description).substr(0, 3));
  auto info = ExplorationTestAccess::feedResult(exploration, tid, ds, 0.5, queue);
  TEST_REAL_SIMILAR(info.metric.score, 4.0)
  TEST_EQUAL(info.group.variant_index, 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0047: the Kinetic Model decides by the modeled maximum and always follows up
/////////////////////////////////////////////////////////////

START_SECTION(kinetic_model_sweep_emits_the_follow_up_at_the_modeled_maximum)
{
  // Summary: six synthetic pre-scans carrying ubiquitin 9+'s real shares at the study's offset-0 grid
  // (1, 11, 21, 31, 41, 50); the group must end with ONE follow-up at the modeled maximum the vendored
  // pipeline computed for that sweep (kinetic_expected.tsv, offset 0), no winner, no cascade.
  // Purpose: ADR-0047 d.4 / d.5, and "the correct energy is found" on real shares (layer B).
  const auto shares = ubiquitinOffsetZeroShares();
  ABORT_IF(shares.size() != 6)
  const auto expected = kineticExpected("Ubiquitin", "9", "1", 0);
  ABORT_IF(expected.first < 1)

  Config cfg{kineticCfg(1.0, 50.0, 10.0)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);
  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  TEST_EQUAL((int)cmds.size(), 6)                        // 1, 11, 21, 31, 41, 50: top edge acquired, no baseline
  ABORT_IF(cmds.size() != 6)
  for (size_t i = 0; i < 6; ++i) TEST_REAL_SIMILAR(cmds[i].stages[0].collision_energy, (double)shares[i].first)

  Exploration::FeedResultInfo last;
  for (size_t i = 0; i < cmds.size(); ++i)
  {
    last = feedShare(exploration, queue, cmds[i], 800.0, shares[i].second);
    TEST_REAL_SIMILAR(last.metric.precursor_share, shares[i].second)   // the share the engine read IS the share fed
    TEST_REAL_SIMILAR(last.metric.score, -1.0)                          // no pre-scan is scored under this metric
  }
  // The completing row carries the fit and the decision.
  TEST_EQUAL(last.kinetic.fitted, true)
  TEST_EQUAL(last.kinetic.fit_points, 6)
  TEST_EQUAL(last.kinetic.energy_source, std::string("modeled"))
  TEST_REAL_SIMILAR(last.kinetic.k2, expected.second)                   // the authored k2 is the fold k2 of the table
  TEST_EQUAL(last.kinetic.modeled_ce, expected.first)                   // identical: the engine IS the study's lattice fit
  TEST_EQUAL(last.kinetic.p0 > 0.0 && last.kinetic.width > 0.0 && last.kinetic.modeled_b > 0.0, true)
  TEST_EQUAL(last.group.winner_tracking_id.empty(), true)               // no pre-scan wins under this metric
  // Exactly one follow-up: an MS2 'R' at the modeled energy, built from the un-overridden config.
  TEST_EQUAL((int)last.commands.size(), 1)
  ABORT_IF(last.commands.empty())
  TEST_EQUAL(last.commands[0].msn_level, 2)
  TEST_EQUAL(std::string(last.commands[0].scan_description)[3], 'R')
  TEST_REAL_SIMILAR(last.commands[0].stages[0].collision_energy, (double)last.kinetic.modeled_ce)
  TEST_STRING_EQUAL(last.commands[0].analyzer, "Orbitrap")
  TEST_EQUAL(exploration.activeGroupCount(), 0)
}
END_SECTION

START_SECTION(kinetic_model_follows_up_with_and_without_overrides)
{
  // ADR-0047 d.5: the follow-up does not depend on overrides. Same shares, trap overrides -> pre-scans on the
  // IonTrap, follow-up on the Orbitrap at the same modeled energy as the no-override run.
  const auto shares = ubiquitinOffsetZeroShares();
  ABORT_IF(shares.size() != 6)
  int modeled_plain = 0, modeled_trap = 0;
  for (int arm = 0; arm < 2; ++arm)
  {
    const std::string ov = arm == 0 ? "" : R"(, "overrides": { "analyzer": "IonTrap" })";
    Config cfg{kineticCfg(1.0, 50.0, 10.0, R"({"HCD": 0.039046})", ov)};
    ScanCommandQueue queue(cfg);
    FragmentAnalysis fragments(cfg);
    Exploration exploration(cfg, fragments);
    auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
    auto cmds = exploration.initiate(2, pg, 3, queue);
    ABORT_IF(cmds.size() != 6)
    for (const auto& c : cmds) TEST_STRING_EQUAL(c.analyzer, arm == 0 ? "Orbitrap" : "IonTrap")
    Exploration::FeedResultInfo last;
    for (size_t i = 0; i < cmds.size(); ++i) last = feedShare(exploration, queue, cmds[i], 800.0, shares[i].second);
    TEST_EQUAL((int)last.commands.size(), 1)
    ABORT_IF(last.commands.empty())
    TEST_STRING_EQUAL(last.commands[0].analyzer, "Orbitrap")          // the follow-up is production-grade either way
    (arm == 0 ? modeled_plain : modeled_trap) = last.kinetic.modeled_ce;
  }
  TEST_EQUAL(modeled_plain > 0, true)
  TEST_EQUAL(modeled_plain, modeled_trap)                              // the share is analyzer-independent in this feed
}
END_SECTION

START_SECTION(kinetic_model_no_fit_aborts_without_a_follow_up)
{
  // Two finite shares (four scans empty) -> fewer than kMinPoints -> no fit, no follow-up, group erased
  // (ADR-0047 d.7). No fallback to the authored collision energy: the owner chose to lose the target.
  Config cfg{kineticCfg(1.0, 50.0, 10.0)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);
  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  ABORT_IF(cmds.size() != 6)
  Exploration::FeedResultInfo last;
  last = feedShare(exploration, queue, cmds[0], 800.0, 0.9);
  last = feedShare(exploration, queue, cmds[1], 800.0, 0.8);
  for (size_t i = 2; i < cmds.size(); ++i) last = feedShare(exploration, queue, cmds[i], 800.0, -1.0);   // empty scans
  TEST_EQUAL(last.kinetic.fitted, false)
  TEST_EQUAL(last.kinetic.fit_points, 2)
  TEST_EQUAL(last.kinetic.energy_source, std::string("none"))
  TEST_EQUAL(last.kinetic.modeled_ce, 0)
  TEST_EQUAL((int)last.commands.size(), 0)
  TEST_EQUAL(exploration.activeGroupCount(), 0)
}
END_SECTION

START_SECTION(kinetic_model_excludes_non_finite_shares)
{
  // Review Focus 2: an empty scan in the middle is excluded (share -1), the remaining five still fit.
  const auto shares = ubiquitinOffsetZeroShares();
  ABORT_IF(shares.size() != 6)
  Config cfg{kineticCfg(1.0, 50.0, 10.0)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);
  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  ABORT_IF(cmds.size() != 6)
  Exploration::FeedResultInfo last;
  for (size_t i = 0; i < cmds.size(); ++i)
    last = feedShare(exploration, queue, cmds[i], 800.0, i == 2 ? -1.0 : shares[i].second);   // CE 21 empty
  TEST_EQUAL(last.kinetic.fitted, true)
  TEST_EQUAL(last.kinetic.fit_points, 5)
  TEST_EQUAL((int)last.commands.size(), 1)
}
END_SECTION

START_SECTION(kinetic_model_flat_shares_choose_the_upper_edge_and_flag_it)
{
  // Review Focus 1: shares 0.80, 0.80, 0.79, 0.79, 0.78, 0.78 -> the curve's fall lies beyond the range. An
  // edge maximum is still acquired, and the completing row says which edge (ADR-0047 d.4).
  Config cfg{kineticCfg(1.0, 50.0, 10.0)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);
  auto pg = makeSyntheticPeakGroup(800.0, 2400.0, 3);
  auto cmds = exploration.initiate(2, pg, 3, queue);
  ABORT_IF(cmds.size() != 6)
  const double flat[6] = {0.80, 0.80, 0.79, 0.79, 0.78, 0.78};
  Exploration::FeedResultInfo last;
  for (size_t i = 0; i < cmds.size(); ++i) last = feedShare(exploration, queue, cmds[i], 800.0, flat[i]);
  TEST_EQUAL(last.kinetic.fitted, true)
  TEST_EQUAL(last.kinetic.modeled_ce, 50)
  TEST_EQUAL(last.kinetic.at_upper_edge, true)
  TEST_EQUAL(last.kinetic.energy_source, std::string("modeled_edge_upper"))   // one column carries both facts
  TEST_EQUAL((int)last.commands.size(), 1)                                      // an edge maximum is still acquired
  ABORT_IF(last.commands.empty())
  TEST_REAL_SIMILAR(last.commands[0].stages[0].collision_energy, 50.0)
}
END_SECTION

START_SECTION(kinetic_model_ms3_reads_the_stage1_window_and_follows_up)
{
  // Review Focus 3 + ADR-0047 d.2 at MS3, on the REAL b32 CID sweep of 2026-07-08 (test-data/spectra/
  // FIXTURES.md): stage 0 is cytC 12+ at m/z 1027.31 (the ms2 context), stage 1 is b32 (4+) at m/z 1027.63.
  // The 13 pre-scans at CID 15..39 step 2 are the cut fixtures; the follow-up is an MS3 at the modeled energy
  // with stage 0 untouched.
  const std::map<int, std::string> b32 = {
    {15, "450"}, {17, "451"}, {19, "452"}, {21, "453"}, {23, "455"}, {25, "456"}, {27, "457"},
    {29, "458"}, {31, "459"}, {33, "460"}, {35, "461"}, {37, "462"}, {39, "463"}};
  Config cfg{kineticMs3Cfg(15.0, 39.0, 2.0)};
  ScanCommandQueue queue(cfg);
  FragmentAnalysis fragments(cfg);
  Exploration exploration(cfg, fragments);

  ScanCommand ms2_ctx = queue.buildMS2(makeSyntheticPeakGroup(1027.3131, 12315.67, 12), 12, cfg.level(2).scans[0], 2, 0);
  PeakGroup frag_pg = makeSyntheticPeakGroup(1027.6311, 4106.49, 4);
  auto cmds = exploration.initiate(3, frag_pg, 4, queue, nullptr, &ms2_ctx, 'b', 32);
  TEST_EQUAL((int)cmds.size(), 13)                                              // 15, 17, ..., 39; no baseline
  ABORT_IF(cmds.size() != 13)
  auto group = ExplorationTestAccess::group(exploration, 1);
  TEST_EQUAL(group.msn_level, 3)
  TEST_REAL_SIMILAR(group.isolation_width, 2.0)                                 // the MS3 floor: [1026.63, 1028.63]

  Exploration::FeedResultInfo last;
  int fed = 0;
  for (size_t i = 0; i < cmds.size(); ++i)
  {
    const int ce = static_cast<int>(std::lround(cmds[i].stages[1].collision_energy));
    ABORT_IF(b32.find(ce) == b32.end())
    // The directory literal stands alone so the container's fixture-path pre-check (a grep for the
    // relative fixture-path literals in the test sources) sees a path that exists; the file name is
    // assembled from it.
    auto scans = loadTsvScans(std::string("../../FlashIDA/test-data/spectra/")
                              + "ms3_cytc3_b32_ce" + std::to_string(ce) + "_scan" + b32.at(ce) + ".txt");
    ABORT_IF(scans.empty())
    const auto& d = scans[0];
    last = exploration.feedResult(cmds[i].scan_id, d.mzs.data(), d.ints.data(), static_cast<int>(d.mzs.size()), d.rt, queue);
    TEST_EQUAL(last.metric.precursor_share >= 0.0 && last.metric.precursor_share <= 1.0, true)
    ++fed;
  }
  TEST_EQUAL(fed, 13)
  TEST_EQUAL(last.kinetic.fitted, true)
  TEST_EQUAL(last.kinetic.fit_points, 13)
  TEST_EQUAL(last.kinetic.modeled_ce >= 1 && last.kinetic.modeled_ce <= 39, true)
  TEST_REAL_SIMILAR(last.kinetic.k2, 0.0388)
  TEST_EQUAL((int)last.commands.size(), 1)
  ABORT_IF(last.commands.empty())
  TEST_EQUAL(last.commands[0].msn_level, 3)
  TEST_EQUAL(last.commands[0].num_stages, 2)
  TEST_EQUAL(std::string(last.commands[0].scan_description)[3], 'R')
  TEST_REAL_SIMILAR(last.commands[0].stages[1].collision_energy, (double)last.kinetic.modeled_ce)
  TEST_REAL_SIMILAR(last.commands[0].stages[0].collision_energy, ms2_ctx.stages[0].collision_energy)   // stage 0 untouched
  TEST_EQUAL((int)last.ms2_context_cache.size(), 1)                            // the regular MS3 path can identify it
  TEST_EQUAL(exploration.activeGroupCount(), 0)

  // The stage-1 window, on a group whose two stages sit far apart (800 vs 500 Th): a spectrum with 1000 at the
  // STAGE-0 m/z and 100 at the STAGE-1 m/z reads share 100/1100 -- reading stage 0's window would give 1000/1100.
  Config cfg2{kineticMs3Cfg(15.0, 39.0, 2.0)};
  ScanCommandQueue q2(cfg2);
  FragmentAnalysis f2(cfg2);
  Exploration e2(cfg2, f2);
  ScanCommand ctx2 = q2.buildMS2(makeSyntheticPeakGroup(800.0, 2400.0, 3), 3, cfg2.level(2).scans[0], 2, 0);
  auto cmds2 = e2.initiate(3, makeSyntheticPeakGroup(500.0, 1000.0, 2), 2, q2, nullptr, &ctx2, 'b', 5);
  ABORT_IF(cmds2.empty())
  std::vector<double> mzs{500.0, 800.0}, ints{100.0, 1000.0};
  auto info2 = e2.feedResult(cmds2[0].scan_id, mzs.data(), ints.data(), 2, 1.0, q2);
  TEST_REAL_SIMILAR(info2.metric.precursor_share, 100.0 / 1100.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


START_SECTION(a_monitor_scan_leaves_the_decision_state_bit_identical)
{
  // THE ADR, as one assertion. A monitor scan is DEFINED by what it does not write, so this compares
  // the engine's WHOLE acquisition memory across one rather than a hand-picked subset of maps.
  // SurveyMemory's operator== is compiler-generated, so a container added to that struct years from
  // now is covered here with no edit to this test.
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  std::string cfg_str = monitorScanConfig(kMonitorOn);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));
  TEST_EQUAL(FLASHIdaTestAccess::hasMonitorDeconv(*ida), true)   // the feature is actually on

  int pushed = bootstrapExplorationGroup(ida, ms1_scans);
  ABORT_IF(pushed <= 0)                                          // a sweep really is running

  // Force one to be due, then take it off the drain.
  FLASHIdaTestAccess::backdateMonitorClock(*ida, 0, 60000);
  ScanCommand mon{};
  TEST_EQUAL(drainUntilRole(ida, ScanRole::Monitor, mon), true)

  // Snapshot AFTER the pull, so the commands discarded above cannot be mistaken for the effect.
  auto before      = FLASHIdaTestAccess::surveyMemory(*ida);
  int  pid_before  = FLASHIdaTestAccess::nextPrecursorId(*ida);
  double cv_before = FLASHIdaTestAccess::faimsCurrentCV(*ida);
  size_t q_before[4];
  for (int p = 0; p < 4; ++p) q_before[p] = FLASHIdaTestAccess::queueSize(*ida, p);

  const ScanData& s = ms1_scans[0];
  ida->processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 1,
                   mon.scan_description, 0.0, instrumentScanNumberOf(s));

  TEST_EQUAL(FLASHIdaTestAccess::surveyMemory(*ida) == before, true)      // no acquisition memory written
  TEST_EQUAL(FLASHIdaTestAccess::nextPrecursorId(*ida), pid_before)       // no precursor id minted
  TEST_REAL_SIMILAR(FLASHIdaTestAccess::faimsCurrentCV(*ida), cv_before)  // the FAIMS wheel did not turn
  for (int p = 0; p < 4; ++p)
    TEST_EQUAL(FLASHIdaTestAccess::queueSize(*ida, p), q_before[p])       // nothing pushed

  // ANTI-VACUOUS MIRROR. Without this the test passes just as well if 'M' scans are dropped on the
  // floor entirely -- the fingerprint would be unchanged for the wrong reason. Feed the SAME spectrum
  // as a SURVEY and require the fingerprint to MOVE.
  // ANTI-VACUOUS MIRROR, on a FRESH engine driven until a survey genuinely SELECTS.
  //
  // Both details are load-bearing. Re-feeding this spectrum to `ida` moves nothing for a real engine
  // reason rather than a bug -- the species is already acquired and the qscore ratchet refuses one
  // whose score does not IMPROVE, so the maps rebuild identically at any rt. And a single survey of
  // ms1_scans[0] moves nothing either, because acquisition memory records what was ACQUIRED: a
  // spectrum with nothing selectable runs the whole survey arm and writes none of it.
  // bootstrapExplorationGroup returns > 0 only once a survey has actually selected and pushed, which
  // is exactly the precondition this mirror needs.
  FLASHIda* fresh = new FLASHIda(const_cast<char*>(cfg_str.c_str()));
  auto fresh_before = FLASHIdaTestAccess::surveyMemory(*fresh);
  ABORT_IF(bootstrapExplorationGroup(fresh, ms1_scans) <= 0)
  TEST_EQUAL(FLASHIdaTestAccess::surveyMemory(*fresh) == fresh_before, false)

  delete fresh;
  delete ida;
}
END_SECTION

START_SECTION(a_monitor_scan_is_never_priority_three)
{
  // Priority 3 is ADR-0031's idle-survey sentinel, and FIVE role-blind drain loops break on it
  // (FLASHIdaWrapper.cs:459/:477, ContinuityTestHarness.cs:123, plus the two fixture classifiers).
  // A monitor scan emitted at 3 would silently truncate every offline drive -- the one way this
  // feature could break something with no connection to it.
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  std::string cfg_str = monitorScanConfig(kMonitorOn);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));
  AcqResult a = runInterleaved(ida, ms1_scans, std::vector<ScanData>{});

  int monitors = 0;
  for (const auto& c : a.all_cmds)
  {
    if (roleOf(c) != ScanRole::Monitor) continue;
    ++monitors;
    TEST_EQUAL(c.priority, 0)          // the lane cycle-time and CV-transition MS1s already use
    TEST_EQUAL(c.msn_level, 1)
    TEST_EQUAL(c.is_agc, 0)
  }
  TEST_EQUAL(monitors > 0, true)       // anti-vacuous: some actually fired
  delete ida;
}
END_SECTION

START_SECTION(a_monitor_scan_does_not_reset_the_survey_clock)
{
  // Pins the roleDecides conjunct on the drain tail's recordMS1Time. A monitor scan satisfies
  // `msn_level == 1 && is_agc == 0` exactly as a survey does, so without that conjunct it stands in
  // for one -- character-for-character the ADR-0031 bug, where the idle path called recordAGCTime()
  // and the authored AGC interval stopped governing the cadence. cycle_time is the only config under
  // which the survey clock is observable at all, which is why it is enabled here and nowhere else.
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  std::string cfg_str = monitorScanConfig(kMonitorOn, "", true);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));
  ABORT_IF(bootstrapExplorationGroup(ida, ms1_scans) <= 0)

  FLASHIdaTestAccess::backdateMonitorClock(*ida, 0, 60000);
  ScanCommand mon{};
  TEST_EQUAL(drainUntilRole(ida, ScanRole::Monitor, mon), true)

  // The clock the cycle-time gate reads must not have been reset by that dequeue. With value_ms = 1
  // it is already overdue, so "still overdue" is the observable.
  uint64_t since = FLASHIdaTestAccess::queue(*ida).msSinceLastMS1();
  TEST_EQUAL(since > 1, true)
  delete ida;
}
END_SECTION

START_SECTION(monitor_cadence_anchors_each_episode_and_then_honours_the_interval)
{
  // The cadence, exactly as authored: one ANCHOR when a level starts sweeping, then one per
  // interval_ms while it keeps sweeping, and RE-ARMED once it goes quiet. Exercised on the queue
  // primitive directly, because the alternative is waiting 30 real seconds -- which is precisely why
  // agc_interval_seconds' equivalent branch has never had a test (ADR-0031).
  std::string cfg_str = monitorScanConfig(kMonitorOn);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));
  ScanCommandQueue& q = FLASHIdaTestAccess::queue(*ida);

  const double interval = 30000.0;

  // Not sweeping: never due -- and this call is also what keeps the level re-armed.
  TEST_EQUAL(q.takeMonitorDue(0, false, interval), false)

  // Sweep starts -> ANCHOR immediately, independent of the interval.
  TEST_EQUAL(q.takeMonitorDue(0, true, interval), true)
  // Still sweeping, no time passed -> not due.
  TEST_EQUAL(q.takeMonitorDue(0, true, interval), false)

  // Interval elapses -> due once, and the clock is CONSUMED (a second ask returns false).
  q.backdateMonitorClock(0, interval + 1);
  TEST_EQUAL(q.takeMonitorDue(0, true, interval), true)
  TEST_EQUAL(q.takeMonitorDue(0, true, interval), false)

  // Sweep ends, then a new one starts -> RE-ARMED, so it owes a fresh anchor. Drop the re-arm and the
  // next episode silently begins with no data point at all.
  TEST_EQUAL(q.takeMonitorDue(0, false, interval), false)
  TEST_EQUAL(q.takeMonitorDue(0, true, interval), true)

  // PER LEVEL: level 3 has its own arm and clock, untouched by everything above.
  TEST_EQUAL(q.takeMonitorDue(1, true, interval), true)
  TEST_EQUAL(q.takeMonitorDue(1, true, interval), false)
  TEST_EQUAL(q.takeMonitorDue(0, true, interval), false)   // level 2 still not due
  delete ida;
}
END_SECTION

START_SECTION(a_monitor_scan_does_not_disturb_its_sweep)
{
  // monitor_deconv_ is a SEPARATE Deconvolution instance, and this is why. deconv_ holds ONE
  // SpectralDeconvolution shared between MS1 and MSn, and deconvolveMS1 restores the global charge
  // window on it (Deconvolution.cpp, flagged LOAD-BEARING). Routed through deconv_, a monitor scan
  // would reset that window under the running sweep and land in the MS1 slot the survey arm reads.
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  std::string cfg_str = monitorScanConfig(kMonitorOn);
  FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));
  ABORT_IF(bootstrapExplorationGroup(ida, ms1_scans) <= 0)

  FLASHIdaTestAccess::backdateMonitorClock(*ida, 0, 60000);
  ScanCommand mon{};
  TEST_EQUAL(drainUntilRole(ida, ScanRole::Monitor, mon), true)
  TEST_EQUAL(FLASHIdaTestAccess::hasMonitorDeconv(*ida), true)

  const ScanData& s = ms1_scans[0];
  ida->processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 1,
                   mon.scan_description, 0.0, instrumentScanNumberOf(s));

  // The group is still active -- a monitor scan cancels nothing and completes nothing.
  TEST_EQUAL(FLASHIdaTestAccess::explorationActive(*ida), true)
  delete ida;
}
END_SECTION

START_SECTION(the_engine_record_decides_the_role_not_the_echo)
{
  // ADR-0008/0035 separate a scan's identity channels but never said which copy of its ROLE wins.
  // It is the ENGINE's: cmd_role comes from the queued command via pending_scan_map_, not from the
  // description the instrument echoed back. Reading the echo instead would let an altered trailer
  // flip a survey into an observation -- and a survey silently observed acquires nothing.
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())

  std::string cfg_str = monitorScanConfig(kMonitorOn);
  const ScanData& s = ms1_scans[0];

  // DIRECTION 1 -- engine record says 'S', the echo claims 'M'. The SURVEY arm must run.
  // A fresh engine, so "memory moved" is unconditional (see the mirror in
  // a_monitor_scan_leaves_the_decision_state_bit_identical for why re-feeding an established
  // engine proves nothing).
  {
    FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));
    ScanCommand surv{};
    TEST_EQUAL(drainUntilRole(ida, ScanRole::Survey, surv, 8), true)

    std::string spoofed(surv.scan_description);
    ABORT_IF(spoofed.size() < 4)
    spoofed[3] = 'M';

    // WHICH DECONVOLUTION ENGINE CARRIES THE FED RT is the discriminator: the survey arm deconvolves
    // into deconv_, the monitor arm into its own monitor_deconv_. Acquisition memory would not do --
    // it moves only if the survey also SELECTED, and this spectrum need not.
    const double probe_rt = 777.25;
    ida->processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), probe_rt, 1,
                     spoofed.c_str(), 0.0, instrumentScanNumberOf(s));
    TEST_REAL_SIMILAR(FLASHIdaTestAccess::surveyDeconvRT(*ida), probe_rt)      // SURVEY arm ran
    TEST_NOT_EQUAL(FLASHIdaTestAccess::monitorDeconvRT(*ida), probe_rt)        // observing arm did not
    delete ida;
  }

  // DIRECTION 2 -- engine record says 'M', the echo claims 'S'. The MONITOR arm must run.
  // Without both directions this test is half a test: direction 1 alone also passes for an engine
  // that simply ignores the marker and always runs the survey arm.
  {
    FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg_str.c_str()));
    ABORT_IF(bootstrapExplorationGroup(ida, ms1_scans) <= 0)
    FLASHIdaTestAccess::backdateMonitorClock(*ida, 0, 60000);

    ScanCommand mon{};
    TEST_EQUAL(drainUntilRole(ida, ScanRole::Monitor, mon), true)

    std::string spoofed(mon.scan_description);
    ABORT_IF(spoofed.size() < 4)
    spoofed[3] = 'S';

    auto before = FLASHIdaTestAccess::surveyMemory(*ida);
    const double probe_rt = 888.75;
    ida->processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), probe_rt, 1,
                     spoofed.c_str(), 0.0, instrumentScanNumberOf(s));

    // The MONITOR arm ran: its private engine carries the fed rt, the survey's does not, and no
    // acquisition memory moved -- despite the echoed description claiming 'S'.
    TEST_REAL_SIMILAR(FLASHIdaTestAccess::monitorDeconvRT(*ida), probe_rt)
    TEST_NOT_EQUAL(FLASHIdaTestAccess::surveyDeconvRT(*ida), probe_rt)
    TEST_EQUAL(FLASHIdaTestAccess::surveyMemory(*ida) == before, true)
    delete ida;
  }
}
END_SECTION

START_SECTION(roleName_matches_the_legacy_scan_type_strings)
{
  // scan_commands.tsv's `scan_type` column is a GOLDEN column and ScanRole.h's table is now its only
  // source. These six pairs are exactly what IdaLogger::scanTypeFromDescription_'s switch returned
  // before ADR-0042 moved it: the byte-identity proof for 28 golden files, as a test rather than as
  // a diff someone has to remember to run.
  TEST_EQUAL(std::string(roleName(ScanRole::Survey)),             "survey")
  TEST_EQUAL(std::string(roleName(ScanRole::Agc)),                "agc")
  TEST_EQUAL(std::string(roleName(ScanRole::Identification)),     "recording")
  TEST_EQUAL(std::string(roleName(ScanRole::Quantification)),     "quantification")
  TEST_EQUAL(std::string(roleName(ScanRole::FollowUp)),           "conditional")
  TEST_EQUAL(std::string(roleName(ScanRole::ExplorationVariant)), "exploration")
  TEST_EQUAL(std::string(roleName(ScanRole::Monitor)),            "monitor")   // the only new value

  TEST_EQUAL(roleMarker(ScanRole::Survey),             'S')
  TEST_EQUAL(roleMarker(ScanRole::Agc),                'A')
  TEST_EQUAL(roleMarker(ScanRole::Identification),     'R')
  TEST_EQUAL(roleMarker(ScanRole::Quantification),     'Q')
  TEST_EQUAL(roleMarker(ScanRole::FollowUp),           'C')
  TEST_EQUAL(roleMarker(ScanRole::ExplorationVariant), 'E')
  TEST_EQUAL(roleMarker(ScanRole::Monitor),            'M')

  // AGC is the only role whose peaks the engine does not read; Monitor is the only other that drives
  // no decision. Both are asserted as EXACT sets, so widening either has to be deliberate.
  for (size_t i = 0; i < kScanRoleTraits.size(); ++i)
  {
    TEST_EQUAL(kScanRoleTraits[i].observes, kScanRoleTraits[i].role != ScanRole::Agc)
    TEST_EQUAL(kScanRoleTraits[i].decides,
               !(kScanRoleTraits[i].role == ScanRole::Agc || kScanRoleTraits[i].role == ScanRole::Monitor))
  }

  // An UNDECODABLE description must resolve as a survey does -- observes and decides both true.
  // Failing closed here would silently stop processing every scan the engine did not mint, and it is
  // what keeps the migrated read sites byte-identical with the char comparisons they replaced.
  ScanCommand blank{};
  blank.scan_description[0] = '\0';
  TEST_EQUAL(roleOf(blank).has_value(), false)
  TEST_EQUAL(roleObserves(roleOf(blank)), true)
  TEST_EQUAL(roleDecides(roleOf(blank)), true)
  TEST_EQUAL(std::string(roleNameOf(blank)), "unknown")
}
END_SECTION

START_SECTION(monitor_ms1_config_is_strict_on_both_ends)
{
  // THROW on a non-positive interval when enabled: a monitor scan would then be due on EVERY drain,
  // minting a priority-0 MS1 ahead of the sweep indefinitely so the sweep could never progress. That
  // is a hang, not a degraded run, and it must not load.
  TEST_EXCEPTION(std::invalid_argument,
      Config cfg{monitorScanConfig(", \"monitor_ms1\": { \"enabled\": true, \"interval_ms\": 0 }")})
  TEST_EXCEPTION(std::invalid_argument,
      Config cfg{monitorScanConfig(", \"monitor_ms1\": { \"enabled\": true, \"interval_ms\": -1 }")})

  // An unknown key inside the block is rejected like every other section (ADR-0007).
  TEST_EXCEPTION(std::invalid_argument,
      Config cfg{monitorScanConfig(", \"monitor_ms1\": { \"enabled\": true, \"interval_sec\": 30 }")})

  // A non-positive interval while DISABLED is fine: the key is inert, and a config must stay loadable
  // when the feature is switched off.
  {
    Config cfg{monitorScanConfig(", \"monitor_ms1\": { \"enabled\": false, \"interval_ms\": 0 }")};
    TEST_EQUAL(cfg.level(2).monitor_ms1_enabled, false)
  }

  // Values reach the level they were authored on, and the two levels are INDEPENDENT.
  {
    Config cfg{monitorScanConfig(", \"monitor_ms1\": { \"enabled\": true, \"interval_ms\": 1111 }",
                                 ", \"monitor_ms1\": { \"enabled\": false, \"interval_ms\": 2222 }",
                                 false, "fragment_count")};
    TEST_EQUAL(cfg.level(2).monitor_ms1_enabled, true)
    TEST_REAL_SIMILAR(cfg.level(2).monitor_ms1_interval_ms, 1111.0)
    TEST_EQUAL(cfg.level(3).monitor_ms1_enabled, false)
    TEST_REAL_SIMILAR(cfg.level(3).monitor_ms1_interval_ms, 2222.0)
  }

  // Absent block => off, at the documented default. This is the state of every committed config.
  {
    Config cfg{monitorScanConfig("")};
    TEST_EQUAL(cfg.level(2).monitor_ms1_enabled, false)
    TEST_REAL_SIMILAR(cfg.level(2).monitor_ms1_interval_ms, 30000.0)
  }

  // Enabled on a level that never sweeps is INERT, not invalid: it must load (and the engine prints a
  // [CONFIG-WARN]). Turning exploration off must never invalidate a config -- ADR-0039's call for an
  // enriched_in that no objective reads, taken again here.
  {
    Config cfg{monitorScanConfig("", kMonitorOn, false, "none")};
    TEST_EQUAL(cfg.level(3).monitor_ms1_enabled, true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// A TRAP PRE-SCAN IS MEASURED, NEVER IDENTIFIED (ADR-0045). With the ADR-0026 narrowing gone a trap
// pre-scan is a real spectrum, and at MS2 every variant used to be deconvolved, matched AND pooled
// under every metric. Unit-resolution peaks have no business in a ppm-tolerance identification, so a
// variant whose OWN command names the ion trap is read for its window sum and nothing else.
/////////////////////////////////////////////////////////////

START_SECTION(a_trap_pre_scan_is_measured_never_identified)
{
  // The same real CA MS2 spectrum is fed to two sweeps that differ ONLY in the override analyzer.
  // Orbitrap: deconvolved, matched (fragment_count > 0), fed to the tracker. IonTrap: the window sum
  // is read and nothing else -- no match, no model, and the deconv accessors FLASHIda.cpp reads for
  // the scan_results row report nothing rather than a STALE spectrum.
  auto make = [](const std::string& analyzer) {
    return std::string(R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": "qscore", "max_precursors": 3,
        "exploration": { "metric": "remaining_precursor", "ce_min": 20.0, "ce_max": 30.0, "ce_step": 5.0,
                         "overrides": { "analyzer": ")" + analyzer + R"(" } } },
      "characterization": { "mode": "off", "protein_sequence": ")" + std::string(ca_sequence) + R"(" },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
        "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "first_mass": 150, "last_mass": 2000 }
      }
    })");
  };
  auto ms2_scans = loadTsvScans(ms2_ca_path);
  ABORT_IF(ms2_scans.empty())
  const auto& d = ms2_scans[0];
  const int precursor_id = 7;

  struct Arm { Exploration::FeedResultInfo variant; bool has_model = false; int deconv_count = -1; bool deconv_null = false; };
  auto run = [&](const std::string& analyzer) {
    Arm a;
    Config cfg{make(analyzer)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    ScanCommandQueue queue(cfg);
    FragmentAnalysis fragments(cfg);
    Exploration exploration(cfg, fragments);
    auto pg = makeSyntheticPeakGroup(968.4916, 29006.3881, 30);   // intact CA z30 -- the fixture's precursor
    auto cmds = exploration.initiate(2, pg, 30, queue);
    TEST_EQUAL(static_cast<int>(cmds.size()), 3)   // CE 20 / 25 / 30; no baseline (ADR-0048)
    if (cmds.empty()) return a;   // ABORT_IF is a `break`, unusable inside a lambda; the assertions below then fail on the defaults
    for (const auto& c : cmds) TEST_STRING_EQUAL(c.analyzer, analyzer)
    a.variant = exploration.feedResult(cmds[0].scan_id, d.mzs.data(), d.ints.data(), static_cast<int>(d.mzs.size()),
                                       d.rt, queue, &tracker, precursor_id);
    a.deconv_count = exploration.explorationDeconvMassCount();
    a.deconv_null  = (exploration.explorationDeconvSpectrum() == nullptr);
    a.has_model    = (tracker.getModel(precursor_id) != nullptr);
    return a;
  };
  Arm orbi = run("Orbitrap");
  Arm trap = run("IonTrap");

  // Positive control: the Orbitrap arm IS identified, so the trap arm's zeros are the skip, not the data.
  TEST_EQUAL(orbi.variant.metric.fragment_count > 0, true)
  TEST_EQUAL(orbi.deconv_count > 0, true)
  TEST_EQUAL(orbi.deconv_null, false)
  TEST_EQUAL(orbi.has_model, true)

  // Measured: the precursor share is read from raw peaks on both arms, so it is analyzer-independent
  // and a real fraction of the scan (ADR-0048).
  TEST_EQUAL(trap.variant.metric.precursor_share > 0.0 && trap.variant.metric.precursor_share <= 1.0, true)
  TEST_REAL_SIMILAR(trap.variant.metric.precursor_share, orbi.variant.metric.precursor_share)
  // Never identified.
  TEST_EQUAL(trap.variant.metric.fragment_count, 0)
  TEST_REAL_SIMILAR(trap.variant.metric.tic_coverage, 0.0)
  TEST_EQUAL(trap.variant.identification.matched_protein.empty(), true)
  TEST_EQUAL(trap.deconv_count, 0)     // FAILS WHEN deconvolveMSn is skipped but the stored spectrum is left STALE
  TEST_EQUAL(trap.deconv_null, true)
  TEST_EQUAL(trap.has_model, false)    // FAILS WHEN the tracker feed is not gated
}
END_SECTION

START_SECTION(a_trap_sweep_must_use_remaining_precursor)
{
  // ADR-0045 decision 2. A trap pre-scan is never deconvolved, so mass_count / fragment_count score every
  // variant 0 and winner selection (seeded -1.0, strictly greater) crowns ce_min every time. Refused at load,
  // at both levels; the level-3 arm is guarded on characterization.mode like the multiplexed pair.
  const std::string trap = R"(, "overrides": { "analyzer": "IonTrap" })";
  auto sweep = [](const std::string& metric, const std::string& ov) {
    return R"("exploration": { "metric": ")" + metric + R"(", "ce_min": 20.0, "ce_max": 40.0, "ce_step": 5.0)" + ov + R"( })";
  };
  const std::string seq_off = R"("mode": "off", "protein_sequence": "PEPTIDER")";   // fragment_count needs a sequence at any level

  Config ok_rp{sweepProbe(", " + sweep("remaining_precursor", trap))};
  TEST_EQUAL(static_cast<int>(ok_rp.level(2).exploration), static_cast<int>(ExplorationMetric::RemainingPrecursor))
  Config ok_orbi{sweepProbe(", " + sweep("mass_count", sweep_overrides))};   // Orbitrap override + counting metric: fine
  TEST_EQUAL(static_cast<int>(ok_orbi.level(2).exploration), static_cast<int>(ExplorationMetric::MassCount))

  TEST_EXCEPTION(std::invalid_argument, Config(sweepProbe(", " + sweep("mass_count", trap))))
  TEST_EXCEPTION(std::invalid_argument, Config(sweepProbe(", " + sweep("fragment_count", trap), seq_off)))
  TEST_EXCEPTION(std::invalid_argument, Config(sweepProbe("", ms3_on + sweep("mass_count", trap), true)))

  // mode off: no MS3 sweep is ever emitted, so no throw (the same guard the multiplexed pair carries).
  Config off_l3{sweepProbe("", R"("mode": "off", )" + sweep("mass_count", trap), true)};
  TEST_EQUAL(static_cast<int>(off_l3.level(3).exploration), static_cast<int>(ExplorationMetric::MassCount))
}
END_SECTION
END_TEST
