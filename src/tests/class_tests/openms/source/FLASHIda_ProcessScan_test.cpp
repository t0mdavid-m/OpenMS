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

#include "FLASHIda_TestHelpers.h"  // ground-truth harness: ScanData/loadTsvScans/AcqResult/runInterleaved/runFullCycle
#include "FLASHIda_TestAccess.h"   // FLASHIdaTestAccess: push/queue/queueSize/explorationActive (private-state access)

#include <nlohmann/json.hpp>  // buildCapConfig edits the parsed document, not the text

#include <fstream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <set>
#include <thread>  // std::this_thread::sleep_for
#include <chrono>  // std::chrono::milliseconds
#include <algorithm>  // std::max
#include <sstream>    // std::istringstream / std::ostringstream
#include <cstdio>     // std::remove
#include <future>     // std::async / std::future_status -- the drain-blocking pins
#include <mutex>      // std::unique_lock over FLASHIdaTestAccess::analysisMutex
#include <iostream>   // the drain pins report which path they covered

using namespace OpenMS;

namespace
{
  // Minimal JSON config for standard DDA mode with score_threshold=0 to accept all peaks
  const char* standard_json = R"({
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
    "max_cv_skip": 0
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
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3,
    "additional_scans": [
      "secondary"
    ]
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
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
        "resolution": 120000
      }
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

  // Config with MS3 targeting enabled: inclusion mode (target_mode=1) + cytC inclusion list
  // pins the precursor, plus MS3 selection via selection_strategy. Mirrors C# CT35.
  const char* ms3_mode1_json = R"({
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
    "max_cv_skip": 0
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
    "inclusion_list": "../../FlashIDA/test-data/configs/inclusion_cytc.txt",
    "ptm_list": ""
  },
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "inclusion",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 1
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

  // Config with IDScore=false, AllCharges=true (activates sortByQScoreAllCharges)
  const char* qscore_allcharges_json = R"({
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
    "max_cv_skip": 0
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
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": true,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3,
    "additional_scans": [
      "secondary"
    ]
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
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
        "resolution": 120000
      }
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

  // Config with quantification enabled and 2 MS2 configs (HCD+ETD, required for quant path)
  const char* quant_json = R"({
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
    "max_cv_skip": 0
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
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3,
    "additional_scans": [
      "secondary"
    ]
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
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
        "resolution": 120000
      }
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": true,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)";

  // Config with tag-based targeting enabled via valid FASTA path, 2 MS2 configs
  const char* tag_targeting_json = R"({
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
    "max_cv_skip": 0
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
    "fasta": "../../FlashIDA/test-data/configs/test_fasta.fasta",
    "inclusion_list": "",
    "ptm_list": ""
  },
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3,
    "additional_scans": [
      "secondary"
    ]
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
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
        "resolution": 120000
      }
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

  // Config with conditional MS2 enabled (no FASTA — tags cannot be found)
  const char* conditional_json = R"({
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
    "max_cv_skip": 0
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
  "conditional_ms2": true,
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
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
    "max_targets": 10
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
      "tagging_follow_up": {
        "analyzer": "Orbitrap",
        "activation": "ETD",
        "collision_energy": 0,
        "reaction_time": 10.0,
        "resolution": 120000
      }
    }
  },
  "tagging": {
    "follow_up_scan": "tagging_follow_up"
  },
  "quantification": {
    "enabled": false,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)";

  // Config with cycle_time enabled and value_ms=0 (always triggers), AGC suppressed
  const char* cycle_time_json = R"({
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
    "max_cv_skip": 0
  },
  "scheduling": {
    "cycle_time": {
      "enabled": true,
      "value_ms": 0
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
    "max_targets": 10
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

  // Config with agc_interval_seconds=0 (AGC triggers immediately)
  const char* agc_fast_json = R"({
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
    "max_cv_skip": 0
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": true,
      "value_ms": 30000
    },
    "agc_interval_seconds": 0
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
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
    "max_targets": 10
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

  // Config with conditional MS2 + tag targeting + FASTA (tags CAN be found)
  const char* conditional_with_tags_json = R"({
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
    "max_cv_skip": 0
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
  "conditional_ms2": true,
  "files": {
    "target_logs": [],
    "fasta": "../../FlashIDA/test-data/configs/test_fasta.fasta",
    "inclusion_list": "",
    "ptm_list": ""
  },
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
    "max_targets": 10
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
      "tagging_follow_up": {
        "analyzer": "Orbitrap",
        "activation": "ETD",
        "collision_energy": 0,
        "reaction_time": 10.0,
        "resolution": 120000
      }
    }
  },
  "tagging": {
    "follow_up_scan": "tagging_follow_up"
  },
  "quantification": {
    "enabled": false,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)";

  // Follow-up parameter-ownership configs. The MS2 and the follow_up_scan differ in EVERY
  // activation-coupled parameter on purpose: if the follow-up inherited them from the triggering
  // MS2 (the defect), each assertion below would read the MS2's value instead of the follow-up's.
  // A config where the two agree would make the test vacuous.
  const char* followup_owns_params_conditional_json = R"({
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
    "max_cv_skip": 0
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
  "conditional_ms2": true,
  "files": {
    "target_logs": [],
    "fasta": "../../FlashIDA/test-data/configs/test_fasta.fasta",
    "inclusion_list": "",
    "ptm_list": ""
  },
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
    "max_targets": 10
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
      "resolution": 120000,
      "reaction_time": 10.0,
      "reagent_max_it": 200.0,
      "reagent_agc_target": 700000
    },
    "additional_ms2": {
      "tagging_follow_up": {
        "analyzer": "Orbitrap",
        "activation": "EThcD",
        "collision_energy": 7,
        "reaction_time": 5.0,
        "reagent_max_it": 111.0,
        "reagent_agc_target": 222,
        "resolution": 120000
      }
    }
  },
  "tagging": {
    "follow_up_scan": "tagging_follow_up"
  },
  "quantification": {
    "enabled": false,
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 1.4
  }
}
)";

  const char* followup_owns_params_quant_json = R"({
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
    "max_cv_skip": 0
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
    "max_targets": 10
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
      "resolution": 120000,
      "reaction_time": 10.0,
      "reagent_max_it": 200.0,
      "reagent_agc_target": 700000
    },
    "ms2_quant": {
      "analyzer": "Orbitrap",
      "activation": "EThcD",
      "collision_energy": 7,
      "reaction_time": 5.0,
      "reagent_max_it": 111.0,
      "reagent_agc_target": 222,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": true,
    "labelling": "tmt6plex",
    "reporter_mz_tol": 0.01,
    "fold_change_threshold": 0.01,
    "conditions": [
      { "name": "a", "channels": ["126", "127", "128"] },
      { "name": "b", "channels": ["129", "130", "131"] }
    ]
  }
}
)";

  // Config with quantification + low fold_change_threshold (any reporter ratio triggers)
  const char* quant_sensitive_json = R"({
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
    "max_cv_skip": 0
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
    "max_targets": 10
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
    "ms2_quant": {
      "analyzer": "Orbitrap",
      "activation": "ETD",
      "collision_energy": 0,
      "reaction_time": 10.0,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": true,
    "labelling": "tmt6plex",
    "reporter_mz_tol": 0.002,
    "fold_change_threshold": 0.01,
    "conditions": [
      { "name": "a", "channels": ["126", "127", "128"] },
      { "name": "b", "channels": ["129", "130", "131"] }
    ]
  }
}
)";

  // Config with agc_interval_seconds=9999 to disable timer-based AGC.
  // Only idle cycle (empty queue) produces AGC/MS1 pairs.
  const char* idle_cycle_json = R"({
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
    "max_cv_skip": 0
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": true,
      "value_ms": 30000
    },
    "agc_interval_seconds": 9999
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
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
    "max_targets": 10
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

  // Config with intensity-based MS1 selection
  const char* intensity_selection_json = R"({
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
    "max_cv_skip": 0
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
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "intensity",
    "max_precursors": 3,
    "additional_scans": [
      "secondary"
    ]
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
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
        "resolution": 120000
      }
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

  // Config with selection=none at MS1
  const char* none_selection_json = R"({
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
    "max_cv_skip": 0
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
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "none",
    "max_precursors": 3
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
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

  // Config with max_targets=1 (cap test)
  const char* max1_json = R"({
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
    "max_cv_skip": 0
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
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 1
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
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

  // TSV file paths relative to the OpenMS build directory (CTest working dir)
  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_standard.txt";
  // 8 rich, co-eluting MS1 scans (raw 326-350) from a complex E. coli run, each with many distinct
  // selectable precursors (charge>=4, ChargeSNR>=1, mass 500-50k) -> the max_targets cap actually bites.
  const std::string ms1_ecoli_rich_tsv_path = "../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt";
  const std::string ms2_tsv_path = "../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt";
  const std::string ms2_tmt_tsv_path = "../../FlashIDA/test-data/spectra/ms2_quant_tmt.txt";
  // ADR-0039. ms2_quant_tmt.txt with the three CONTROL reporter peaks (126/127/128) deleted and
  // nothing else changed -- a species present under treatment and absent from control. That is the
  // one shape the six-channel fixture cannot produce, and the only one that distinguishes an
  // enriched_in implemented on condition_means from one implemented on fold_change.
  const std::string ms2_quant_absent_tsv_path =
      "../../FlashIDA/test-data/spectra/ms2_quant_tmt_absent.txt";
  // The MIRROR: 129/130/131 removed instead, so `treated` is the absent condition and
  // condition_means is [X, 0] rather than [0, X]. Both directions are needed because a sentinel bug
  // spelled `fold_change < 1.0` happens to AGREE with the correct answer when condition[0] is the
  // absent one -- it is only distinguishable when the absent condition is the numerator.
  const std::string ms2_quant_treated_absent_tsv_path =
      "../../FlashIDA/test-data/spectra/ms2_quant_tmt_treated_absent.txt";
  const std::string ms1_cytc_tsv_path = "../../FlashIDA/test-data/spectra/ms1_cytc.txt";
  const std::string ms2_cytc_tsv_path = "../../FlashIDA/test-data/spectra/ms2_cytc_scan149.txt";
  // Fresh real Mode-2 MS3 CytC run (20250121): intact-cytC MS2 (HCD CE40, scan 57) with a rich b/y
  // ladder that FLASHExtender can validate -> real MS3. Extracted via test-scripts/prepare-test-data.py.
  const std::string ms2_cytc_fresh_tsv_path = "../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt";
  const std::string fasta_path = "../../FlashIDA/test-data/configs/test_fasta.fasta";


  // Build a cap-test config from the max1_json template: set the precursor cap, optionally add a
  // second (ETD) MS2 scan config, and enable per-scan results logging.
  //
  // EDITS THE PARSED DOCUMENT, NOT THE TEXT. This used to be three std::string::find calls, each
  // guarded by `if (p != npos)` -- so when the schema reshape renamed the key and re-indented the
  // template, all three edits silently became no-ops and this test ran SIX IDENTICAL CONFIGS and
  // compared them to each other. It stayed green until both the cap edit and the ETD edit died at
  // once; with only one broken it would have kept passing while testing half of nothing.
  //
  // Two specific traps the textual version walked into:
  //   - it searched for "max_targets": 1, but the cap is precursor_selection.max_precursors now.
  //     The only surviving max_targets is characterization.max_targets -- the MS3 budget -- so a
  //     near-miss would have silently changed the wrong knob.
  //   - it appended `, {...ETD...}` after the HCD object to make an array. ms_settings.ms2 is a
  //     bare object now; a second MS2 must be DEFINED in additional_ms2 and then REFERENCED from
  //     precursor_selection.additional_scans, or it parses fine and never dispatches.
  std::string buildCapConfig(int max_targets, bool etd, const std::string& log_dir)
  {
    nlohmann::json j = nlohmann::json::parse(max1_json);

    j["precursor_selection"]["max_precursors"] = max_targets;

    if (etd)
    {
      j["ms_settings"]["additional_ms2"]["etd"] = {
        {"analyzer", "Orbitrap"}, {"activation", "ETD"}, {"collision_energy", 0},
        {"reaction_time", 10.0}, {"resolution", 120000}};
      // Definition alone is inert -- the roster is built from this array.
      j["precursor_selection"]["additional_scans"] = nlohmann::json::array({"etd"});
    }

    // One folder, five fixed basenames. This test only reads scan_results.tsv, but the engine
    // now opens all five streams under log_dir -- there is no way to ask for just one.
    j["runtime"]["log_dir"] = log_dir;
    return j.dump();
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

  // ADR-0039. One quantification config, parameterized on the three knobs the objective sections
  // vary. Deliberately carries NO fasta and conditional_ms2: false, so nothing but the
  // quantification screen can ever push a command -- which is what lets those sections assert on
  // `pushed` directly rather than having to classify what came back.
  //
  // The two conditions are named control/treated to match method_quant.json, and their channel
  // split is the same 3-vs-3. fold_change = mean(control) / mean(treated), so on the TMT fixture
  // (control mean 5088.06, treated mean 9925.51) the species is enriched in TREATED and the raw
  // ratio is < 1 -- which is exactly why direction is authored by condition NAME and never as
  // "up": "up" here would have meant enriched in control.
  std::string quantObjectiveJson(const std::string& fold_change_threshold,
                                 const std::string& identify_key_or_empty,
                                 const std::string& enriched_in_key_or_empty)
  {
    return std::string(R"({
  "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 4,
                     "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10] },
  "flashtnt": { "min_length": 3, "max_length": 8 },
  "faims": { "cv_values": [], "max_cv_skip": 0 },
  "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 },
                  "scan_timeout": { "enabled": false, "value_ms": 30000 },
                  "agc_interval_seconds": 9999999 },
  "conditional_ms2": false,
  "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
  "precursor_selection": { "rt_window": 180, "targeting": "none", "consider_all_charges": false,
                           "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore",
                           "max_precursors": 1 },
  "characterization": { "mode": "off", "max_targets": 3 },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000,
             "agc_target": 800000, "max_it": 246 },
    "ms2": { "analyzer": "Orbitrap", "activation": "ETD", "reaction_time": 7, "resolution": 120000 },
    "ms2_quant": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 30,
                   "resolution": 120000 }
  },
  "quantification": { "enabled": true, "labelling": "tmt6plex", "reporter_mz_tol": 0.002,
                      "fold_change_threshold": )") + fold_change_threshold + R"(,)"
           + identify_key_or_empty + enriched_in_key_or_empty + R"(
                      "conditions": [ { "name": "control", "channels": ["126", "127", "128"] },
                                      { "name": "treated", "channels": ["129", "130", "131"] } ] }
})";
  }

  /// Drive one MS1, feed the rostered 'Q' scan the given MS2 spectrum, and return how many commands
  /// its return pushed. 0 = the verdict bought nothing; 1 = it bought the identification scan.
  /// Asserts the rostered scan really is the 'Q', so a roster regression cannot make this vacuous.
  int quantPushedFor(const char* cfg_json, const std::vector<ScanData>& ms1_scans,
                     const std::vector<ScanData>& ms2_scans)
  {
    FLASHIda ida(const_cast<char*>(cfg_json));
    AcqResult acq = runInterleaved(&ida, ms1_scans, {});
    if (acq.ms2_cmds.empty()) return -1;
    const ScanCommand& q = acq.ms2_cmds[0];
    if (std::strlen(q.scan_description) < 4 || q.scan_description[3] != 'Q') return -2;
    const auto& ms2 = ms2_scans[0];
    return ida.processScan(ms2.mzs.data(), ms2.ints.data(), (int)ms2.mzs.size(), ms2.rt, 2,
                           q.scan_description, 0.0, instrumentScanNumberOf(ms2));
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
  // Drive MS1 by echoing the engine's own survey id (the always-on MS1 gate rejects fabricated ids).
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  delete ida;
}
END_SECTION

// P4-U02: Commands from processScan are dequeued by getNextScanCommand
START_SECTION(processScan_commands_dequeued)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  // Drive MS1 via the engine-id-echo contract; runInterleaved collects every dequeued MS2 command, so the
  // per-command field checks below run over acq.ms2_cmds (each was returned by getNextScanCommand with
  // result==1 inside runInterleaved). The interleaving does not change the spectra/config -> same commands.
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)

  for (const auto& cmd : acq.ms2_cmds)
  {
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

  // After a full drive the queue is drained, so every further call emits a fresh idle survey MS1 --
  // no AGC prescan alternates with them any more (ADR-0031), and this fixture pins
  // agc_interval_seconds high so the scheduled prescan cannot fire either. Two consecutive calls
  // therefore both yield priority-3 surveys, with DIFFERENT tracking ids: each is minted on the
  // call that returns it, rather than one being left queued by a preceding prescan.
  ScanCommand cmd{};
  int idle_result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(idle_result, 1)
  TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_EQUAL(cmd.is_agc, 0)
  TEST_EQUAL(cmd.priority, 3)
  const int first_idle_id = cmd.scan_id;

  int second_result = ida->getNextScanCommand(cmd);
  TEST_EQUAL(second_result, 1)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_EQUAL(cmd.is_agc, 0)
  TEST_EQUAL(cmd.priority, 3)
  TEST_NOT_EQUAL(cmd.scan_id, first_idle_id)

  delete ida;
}
END_SECTION

// P4-U03: ScanCommand fields populated correctly
START_SECTION(processScan_command_fields)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  // Drive MS1 via engine-id echo; inspect the first dequeued MS2 command's fields.
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  const ScanCommand& cmd = acq.ms2_cmds[0];
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

  // Drive MS1 via engine-id echo with NO MS2 feed, so the emitted MS2 commands stay in pending_scan_map_
  // (unfed) and we can feed the first one ourselves below to inspect the MS2-path return value.
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  // The first emitted MS2 command carries the engine's tracking ID in its scan description.
  const ScanCommand& ms2_cmd = acq.ms2_cmds[0];
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)

  // Now process MS2 return with the tracking ID in scan description
  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
  // Should return 0 (no conditional, no MS3, no quant in standard config). Interleaving does not change
  // the standard-config MS2-path follow-up count: this re-derives to the same 0 as before.
  TEST_EQUAL(ms2_result, 0)

  delete ida;
}
END_SECTION

// P4-U04b: unified context-support gate (TODO #8/#9). processScan validates ONCE, before the level
// dispatch, that the resolved queued command supports the inbound scan: ms_level in {1,2,3},
// ctx.msn_level == ms_level, and ctx.num_stages >= required(ms_level). Any violation -> return 0.
START_SECTION(processScan_context_gate_rejects_unsupported_levels)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  // Drive MS1 with NO MS2 feed so an emitted MS2 command (msn_level=2, num_stages=1) stays pending, then
  // feed it back at the given inbound ms_level. Fresh engine per call so resolvePending isn't re-consumed.
  auto feed_ms2_cmd_at_level = [&](int inbound_level) -> int
  {
    FLASHIda ida(const_cast<char*>(standard_json));
    AcqResult acq = runInterleaved(&ida, ms1_scans, {});
    if (acq.ms2_cmds.empty()) return -999;  // setup failure -> fails the assertion below, visibly
    const ScanCommand& ms2_cmd = acq.ms2_cmds[0];
    const auto& ms2 = ms2_scans[0];
    return ida.processScan(ms2.mzs.data(), ms2.ints.data(), (int)ms2.mzs.size(), ms2.rt,
                           inbound_level, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
  };

  // (a) MS4+ unsupported (required_stages < 0): would previously have run the MS3 path. (TODO #9)
  TEST_EQUAL(feed_ms2_cmd_at_level(4), 0)
  // (b) Inbound ms_level=3 but the queued command is msn_level=2 / num_stages=1 (< required 2): the
  //     context cannot support the inbound scan level -> rejected. (TODO #8 + msn_level==ms_level)
  TEST_EQUAL(feed_ms2_cmd_at_level(3), 0)
  // (c) Inbound ms_level=1 but the queued command is msn_level=2 (level mismatch) -> rejected.
  TEST_EQUAL(feed_ms2_cmd_at_level(1), 0)
}
END_SECTION

// P4-U05: processScan with empty spectrum returns 0
START_SECTION(processScan_empty_spectrum)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  int n = ida->processScan(nullptr, nullptr, 0, 0.0, 1, "empty", 0.0, 0);
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

  // Drive MS1 via engine-id echo with NO MS2 feed (so the emitted MS2 commands stay pending and we feed
  // the first one ourselves below to exercise the conditional follow-up).
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  const ScanCommand& ms2_cmd = acq.ms2_cmds[0];
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)
  TEST_EQUAL(ms2_cmd.priority, 2)

  // Process MS2 return — tags found → conditional follow-up at priority 0
  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
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

// P4-U06b (F6): a follow-up MS2 ('C' conditional / 'F' quant) must NOT spawn a SECOND-generation follow-up.
// A follow-up is itself processed through the MS2 path, so without a single-generation guard its (same)
// tag-bearing spectrum re-satisfies the predicate and recurses. Feeds a fired 'C' follow-up's own response
// back and asserts zero new follow-ups. Set/count relationship only -> drift-stable.
START_SECTION(processScan_followup_does_not_recurse)
{
  {
    std::ifstream fasta_check(fasta_path);
    ABORT_IF(! fasta_check.good())
  }
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(conditional_with_tags_json));

  // Drive the FULL acquisition through the canonical harness with a NON-empty ms2_scans: runInterleaved feeds
  // each dequeued MS2 command's own response back BY LEVEL (the engine paces it; bounded by idle>=3/max_iters,
  // so NO raw getNextScanCommand drain loop and no hang). The 'R' MS2 is fed -> tags -> a first-generation 'C'
  // follow-up; that 'C' is itself an MS2-level command, so the harness dequeues it, collects it in acq.ms2_cmds,
  // AND feeds ITS response back too. A second-generation follow-up (pre-fix) would likewise be dequeued + collected.
  AcqResult acq = runInterleaved(ida, ms1_scans, ms2_scans);
  ABORT_IF(acq.ms2_cmds.empty())

  // conditional_with_tags_json has quant disabled, so each original ('R') MS2 yields at most ONE 'C' follow-up.
  // With the single-generation guard (FLASHIda.cpp:984) every follow-up's parent is an 'R', never another
  // follow-up, so followups <= r_count.
  // ISSUE(F6): pre-fix feeding a follow-up re-triggered buildFollowUp -> a 'C' spawned another 'C' -> followups
  // ISSUE      grows past r_count (bounded only by max_iters). followups <= r_count is the no-recursion invariant.
  int r_count = 0, followups = 0;
  for (const auto& c : acq.ms2_cmds)
  {
    std::string d(c.scan_description);
    if (d.size() < 4) continue;
    if (d[3] == 'R')      ++r_count;     // original production MS2
    else if (d[3] == 'C') ++followups;   // conditional follow-up
    // 'F' is gone (ADR-0038) and is deliberately NOT replaced by 'Q' here: a 'Q' is a ROSTERED
    // quantification scan, not a follow-up, so counting it would break the invariant this test
    // asserts. This fixture disables quantification anyway, so no 'Q' is emitted.
  }
  TEST_TRUE(followups >= 1)          // a follow-up WAS produced (the test is meaningful)
  TEST_TRUE(followups <= r_count)    // no follow-up spawned another (single-generation guard holds)

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

  // Drive the WHOLE MS1 -> MS2 -> MS3 chain through the canonical harness. runFullCycle pulls each command,
  // feeds the inclusion-pinned cytC MS1 survey by engine id, then feeds ms2_cytc (ms2_scans[0]) back for the
  // emitted MS2 command -> the b/y fragment matches fire MS3 targets; MS3 commands are collected in
  // result.ms3_cmds (with no ms3_ion_map the MS2 spectrum is fed back for each MS3 via the legacy shortcut,
  // which the engine accepts). Bounded by idle>=3 / auto-sized budget -> no raw drain loop, no hang.
  CycleResult result = runFullCycle(ida, ms1_scans, ms2_scans);
  TEST_EQUAL(result.ms2_cmds.size() > 0, true)
  ABORT_IF(result.ms2_cmds.empty())

  // MS3 at priority 1 dequeues before MS2 at priority 2; the harness collected every MS3 command it drained.
  // Re-express the old per-command checks over result.ms3_cmds (each was returned by getNextScanCommand==1
  // inside runFullCycle). The ms3_count>0 hard positive becomes "at least one MS3 command emitted".
  int ms3_count = 0;
  for (const ScanCommand& out : result.ms3_cmds)
  {
    TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
    ms3_count++;
    TEST_EQUAL(out.priority, 1)    // MS3 at priority 1
    TEST_EQUAL(out.num_stages, 2)  // MS2 precursor + fragment target

    // D1: MS3 scan description uses 3-char ID + type code 'R' (recording)
    std::string ms3_desc(out.scan_description);
    TEST_EQUAL(ms3_desc.size() >= 4, true)
    TEST_EQUAL(ms3_desc[3], 'R')  // MS3 recording
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

  // Test roundtrip via the engine-id-echo drive -> scan_description parsing on the first emitted MS2 command.
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  const ScanCommand& cmd = acq.ms2_cmds[0];
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
  int decoded_id = FLASHIdaTestAccess::queue(*ida).decode(id_str);
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

  // Drive MS1 via engine-id echo with NO MS2 feed: every dequeued command is registered in
  // pending_scan_map_; the fed MS1 surveys resolve (erase) themselves, but the dequeued MS2 commands and
  // the idle AGCs are NOT fed back, so they stay pending and are available for the resolution check below.
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())
  for (const auto& c : acq.ms2_cmds) TEST_EQUAL(std::strlen(c.scan_description) <= 15, true)

  // The pending_scan_map_ entries have timestamps. cleanupExpired is called by getNextScanCommand. With
  // timeout_ms=30000, entries should NOT be expired immediately. Verify by resolving one MS2 tracking ID.
  // Snapshot the pending size first so the post-resolution check is exact (one entry erased), independent of
  // how many AGC/MS2 entries the interleaved drive left pending.
  const size_t pending_before = FLASHIdaTestAccess::queue(*ida).pendingScanMapSize();
  ABORT_IF(pending_before == 0)
  const ScanCommand& ms2_cmd = acq.ms2_cmds[0];

  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
  // Should succeed (entry found in pending_scan_map_)
  // ms2_result can be 0 (no follow-ups) but shouldn't crash
  TEST_EQUAL(ms2_result >= 0, true)

  // Verify pending_scan_map_ entry was erased by first resolution
  int ms2_result2 = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                      (int)ms2.mzs.size(), ms2.rt,
                                      2, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
  // Second resolution with same tracking ID should find nothing (entry erased)
  TEST_EQUAL(ms2_result2, 0)

  // Also verify via accessor: the first resolution erased EXACTLY one pending entry (pending_before - 1).
  // (Re-derived from interleaving: the snapshot baseline replaces the old `total` MS2-only baseline; the
  // engine still erases exactly one entry per resolved MS2 tracking id -- see FLASHIda.cpp:918.)
  TEST_EQUAL(FLASHIdaTestAccess::queue(*ida).pendingScanMapSize(), (size_t)(pending_before - 1))

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

  // Drive MS1 via engine-id echo for each scoring branch.
  AcqResult acq_qscore = runInterleaved(ida_qscore, ms1_scans, {});
  AcqResult acq_qscore_all = runInterleaved(ida_qscore_all, ms1_scans, {});

  // Both branches must produce > 0 commands
  TEST_EQUAL(acq_qscore.ms2_cmds.size() > 0, true)
  TEST_EQUAL(acq_qscore_all.ms2_cmds.size() > 0, true)
  ABORT_IF(acq_qscore.ms2_cmds.empty() || acq_qscore_all.ms2_cmds.empty())

  // First command from each — must have valid precursor_mz
  TEST_EQUAL(std::strlen(acq_qscore.ms2_cmds[0].scan_description) <= 15, true)
  TEST_EQUAL(acq_qscore.ms2_cmds[0].stages[0].precursor_mz > 0, true)

  TEST_EQUAL(std::strlen(acq_qscore_all.ms2_cmds[0].scan_description) <= 15, true)
  TEST_EQUAL(acq_qscore_all.ms2_cmds[0].stages[0].precursor_mz > 0, true)

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

  // Pass 1: drive MS1 via engine-id echo (NO MS2 feed). Every selected MS2 precursor m/z becomes
  // "excluded" (qscore > 0 > threshold 0.0). Exclusion is armed at MS1 selection time (not on the MS2
  // feed), so {} for ms2_scans is correct; the exclusion state persists in the engine into pass 2.
  AcqResult p1 = runInterleaved(ida, ms1_scans, {});
  int total_pass1 = (int)p1.ms2_cmds.size();
  TEST_EQUAL(total_pass1 > 0, true)
  std::set<int> excluded;
  for (const auto& cmd : p1.ms2_cmds) excluded.insert((int)(cmd.stages[0].precursor_mz + 0.5));
  ABORT_IF(excluded.empty())

  // Pass 2: re-drive the SAME scans at the same RTs (within RT_window=180s) on the SAME engine.
  AcqResult p2 = runInterleaved(ida, ms1_scans, {});
  int total_pass2 = (int)p2.ms2_cmds.size();
  std::vector<int> pass2_order;
  for (const auto& out : p2.ms2_cmds) pass2_order.push_back((int)(out.stages[0].precursor_mz + 0.5));
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

  // Pass 1: drive MS1 via engine-id echo (NO MS2 feed) and record, per integer-m/z, the MAX selection
  // qscore (the value the engine gates on at PrecursorSelection.cpp:654; cmd.qscore == pg.getQscore()
  // since AllCharges=false). The exclusion map is armed at MS1 selection time and persists into pass 2.
  AcqResult p1 = runInterleaved(ida, ms1_scans, {});
  int total_pass1 = (int)p1.ms2_cmds.size();
  TEST_EQUAL(total_pass1 > 0, true)
  std::map<int, float> max_q;
  for (const auto& cmd : p1.ms2_cmds)
  {
    int mz = (int)(cmd.stages[0].precursor_mz + 0.5);
    auto it = max_q.find(mz);
    if (it == max_q.end() || cmd.qscore > it->second) max_q[mz] = cmd.qscore;
  }
  // Only masses whose selection qscore exceeded 0.1 arm the exclusion map.
  std::set<int> excluded;
  for (const auto& kv : max_q) if (kv.second > 0.1f) excluded.insert(kv.first);
  TEST_EQUAL(excluded.size() > 0, true)  // non-vacuous: there ARE qscore>0.1 masses to exclude

  // Pass 2: re-drive the SAME scans on the SAME engine.
  AcqResult p2 = runInterleaved(ida, ms1_scans, {});
  int total_pass2 = (int)p2.ms2_cmds.size();
  std::vector<int> pass2_order;
  for (const auto& out : p2.ms2_cmds) pass2_order.push_back((int)(out.stages[0].precursor_mz + 0.5));
  // Exclusion strictly reduces the count; no qscore>0.1 mass is re-selected (hence not picked first).
  // Vacuous-safe over an empty pass2.
  TEST_EQUAL(total_pass2 < total_pass1, true)
  for (int mz : pass2_order) TEST_EQUAL(excluded.count(mz) == 0, true)

  delete ida;
}
END_SECTION

// P4-U07: Quantification follow-up routing — hard positive with sensitive threshold
// ADR-0038: a quantification scan is MEASURED and BUYS the identification scan -- and does nothing
// else. Its activation and energy were chosen to release a reporter ion, not to fragment
// informatively, so tags read off it would come from the wrong spectrum and a 'C' follow-up bought
// here would be a third scan on the precursor before the identification scan has even run.
//
// The control is what makes this a test rather than an observation: tags_with_quant_json is
// conditional_with_tags_json with quantification ADDED and nothing else changed -- same fasta, same
// tagging.follow_up_scan, same spectra. processScan_followup_does_not_recurse already asserts the
// control produces at least one 'C'. So if the is_quant_scan gate on the tagging branch were
// removed, this section would see those same 'C' commands and fail.
START_SECTION(processScan_quant_scan_raises_no_tagging_or_ms3)
{
  const char* tags_with_quant_json = R"({
  "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 4,
                     "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10] },
  "flashtnt": { "min_length": 3, "max_length": 8 },
  "faims": { "cv_values": [], "max_cv_skip": 0 },
  "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 },
                  "scan_timeout": { "enabled": false, "value_ms": 30000 },
                  "agc_interval_seconds": 9999999 },
  "conditional_ms2": true,
  "files": { "target_logs": [], "fasta": "../../FlashIDA/test-data/configs/test_fasta.fasta",
             "inclusion_list": "", "ptm_list": "" },
  "precursor_selection": { "rt_window": 180, "targeting": "none", "consider_all_charges": false,
                           "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore",
                           "max_precursors": 3 },
  "characterization": { "mode": "off", "max_targets": 10 },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000,
             "agc_target": 800000, "max_it": 246 },
    "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
    "ms2_quant": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 30, "resolution": 120000 },
    "additional_ms2": { "tagging_follow_up": { "analyzer": "Orbitrap", "activation": "ETD",
                                               "collision_energy": 0, "reaction_time": 10.0,
                                               "resolution": 120000 } }
  },
  "tagging": { "follow_up_scan": "tagging_follow_up" },
  "quantification": { "enabled": true, "labelling": "tmt6plex", "reporter_mz_tol": 0.002,
                      "fold_change_threshold": 0.01,
                      "conditions": [ { "name": "a", "channels": ["126", "127", "128"] },
                                      { "name": "b", "channels": ["129", "130", "131"] } ] }
})";

  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_tmt_scans = loadTsvScans(ms2_tmt_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_tmt_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(tags_with_quant_json));

  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  ABORT_IF(acq.ms2_cmds.empty())

  // The rostered scan is the 'Q'. Feed its return the TMT spectrum, which the fasta ALSO yields
  // tags from -- that is what the control demonstrates.
  const ScanCommand& q = acq.ms2_cmds[0];
  TEST_EQUAL(q.scan_description[3], 'Q')
  const auto& ms2 = ms2_tmt_scans[0];
  int pushed = ida->processScan(ms2.mzs.data(), ms2.ints.data(), (int)ms2.mzs.size(), ms2.rt,
                                2, q.scan_description, 0.0, instrumentScanNumberOf(ms2));

  // EXACTLY ONE command: the identification scan it bought. Not two, which is what a 'C' alongside
  // it would make, and not more.
  TEST_EQUAL(pushed, 1)

  // Drain everything and confirm what came out: no 'C' anywhere, and no MS3 at any point.
  int c_count = 0, r_count = 0, q_count = 0;
  ScanCommand out{};
  for (int i = 0; i < 40; ++i)
  {
    if (ida->getNextScanCommand(out) != 1) break;
    if (out.msn_level == 1 && out.priority == 3) break;   // idle survey == drained
    if (out.msn_level == 3) { TEST_EQUAL(out.msn_level, 2) }  // fails loudly naming the level
    if (std::strlen(out.scan_description) >= 4)
    {
      if (out.scan_description[3] == 'C') ++c_count;
      if (out.scan_description[3] == 'R') ++r_count;
      if (out.scan_description[3] == 'Q') ++q_count;
    }
    out = ScanCommand{};
  }
  TEST_EQUAL(c_count, 0)          // the whole point: tagging is suppressed on a 'Q'
  TEST_TRUE(r_count >= 1)         // ...but the identification scan WAS bought, so this is not
                                  //    "nothing happened", which would pass vacuously

  // The quantification scans are counted from the DRIVE, not from the drain above: runInterleaved
  // already dequeued them into acq.ms2_cmds, so the post-return drain only sees what the return
  // pushed. Counting them there found 0 and said nothing about the roster.
  int q_rostered = 0;
  for (const auto& c : acq.ms2_cmds)
    if (std::strlen(c.scan_description) >= 4 && c.scan_description[3] == 'Q') ++q_rostered;
  TEST_TRUE(q_rostered >= 1)      // the roster really is emitting quantification scans
  (void)q_count;

  delete ida;
}
END_SECTION

// ---------------------------------------------------------------------------------------------
// ADR-0039: the quantification OBJECTIVE. Until this ADR the decision was a hardcoded
// `verdict == Differential` with no config surface, so "identify everything I could quantify",
// "quantify only" and "only the ones enriched in treated" were all unauthorable.
//
// Thresholds are chosen far from any plausible ratio rather than close to it: isotope correction is
// ON by default (ADR-0038), so the fixture's corrected fold change is NOT its raw 0.512627 and a
// threshold picked to sit just beside that number would be pinning a value this test does not own.
// 0.01 makes every ratio differential; 100.0 makes none.
// ---------------------------------------------------------------------------------------------

START_SECTION(processScan_quant_identify_selects_what_the_verdict_buys)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_tmt_scans = loadTsvScans(ms2_tmt_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_tmt_scans.empty())

  // A Differential verdict under the DEFAULT objective buys the identification scan. This is
  // ADR-0038's behaviour, and it is asserted here so the sections below are read as differences
  // from a known-good baseline rather than in isolation.
  const std::string diff_default = quantObjectiveJson("0.01", "", "");
  TEST_EQUAL(quantPushedFor(diff_default.c_str(), ms1_scans, ms2_tmt_scans), 1)

  // ...and authoring the default explicitly changes nothing. That is what makes ADR-0039
  // byte-identical for every committed config.
  const std::string diff_explicit = quantObjectiveJson("0.01", R"( "identify": "differential",)", "");
  TEST_EQUAL(quantPushedFor(diff_explicit.c_str(), ms1_scans, ms2_tmt_scans), 1)

  // NotDifferential under `differential` buys NOTHING -- the precursor is screened and dropped.
  const std::string notdiff_default = quantObjectiveJson("100.0", "", "");
  TEST_EQUAL(quantPushedFor(notdiff_default.c_str(), ms1_scans, ms2_tmt_scans), 0)

  // ...and that is the exact case `quantified` exists to change: a cleanly measured species now
  // earns identification whether or not it moved.
  const std::string notdiff_quantified =
      quantObjectiveJson("100.0", R"( "identify": "quantified",)", "");
  TEST_EQUAL(quantPushedFor(notdiff_quantified.c_str(), ms1_scans, ms2_tmt_scans), 1)

  // `all` buys regardless of verdict.
  const std::string notdiff_all = quantObjectiveJson("100.0", R"( "identify": "all",)", "");
  TEST_EQUAL(quantPushedFor(notdiff_all.c_str(), ms1_scans, ms2_tmt_scans), 1)

  // `none` buys nothing even on a Differential verdict -- a quantify-only run. Paired with the
  // FIRST assertion in this section, which is the same config differing only in this key, so a
  // 0 here is the objective and not a screen that failed to fire.
  const std::string diff_none = quantObjectiveJson("0.01", R"( "identify": "none",)", "");
  TEST_EQUAL(quantPushedFor(diff_none.c_str(), ms1_scans, ms2_tmt_scans), 0)
}
END_SECTION

START_SECTION(processScan_enriched_in_restricts_by_condition)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_tmt_scans = loadTsvScans(ms2_tmt_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_tmt_scans.empty())

  // The fixture's species is enriched in TREATED (means 5088.06 control vs 9925.51 treated), so
  // naming treated buys and naming control does not. Both configs are identical but for this one
  // string, which is what makes the pair a direction test rather than two unrelated observations.
  const std::string in_treated =
      quantObjectiveJson("0.01", "", R"( "enriched_in": "treated",)");
  TEST_EQUAL(quantPushedFor(in_treated.c_str(), ms1_scans, ms2_tmt_scans), 1)

  const std::string in_control =
      quantObjectiveJson("0.01", "", R"( "enriched_in": "control",)");
  TEST_EQUAL(quantPushedFor(in_control.c_str(), ms1_scans, ms2_tmt_scans), 0)

  // "either" is the default and restores the symmetric test ADR-0038 shipped.
  const std::string either = quantObjectiveJson("0.01", "", R"( "enriched_in": "either",)");
  TEST_EQUAL(quantPushedFor(either.c_str(), ms1_scans, ms2_tmt_scans), 1)

  // Inert under a non-differential objective, and inert means the scan is still bought -- the
  // direction does not quietly survive as a filter. Emits [CONFIG-WARN] at load.
  const std::string inert_all =
      quantObjectiveJson("0.01", R"( "identify": "all",)", R"( "enriched_in": "control",)");
  TEST_EQUAL(quantPushedFor(inert_all.c_str(), ms1_scans, ms2_tmt_scans), 1)
}
END_SECTION

START_SECTION(processScan_enriched_in_survives_a_wholly_absent_condition)
{
  // THE trap test. A species present under treatment and absent from control is the strongest
  // result the experiment can produce: Quantification reports it Differential with
  // fold_change == -1.0, a SENTINEL rather than a ratio (condition_means carries the truth).
  //
  // An `enriched_in` implemented on fold_change therefore reads -1, finds it below any threshold,
  // and silently DROPS exactly the species the experiment exists to find. Nothing else in this
  // suite would catch that: every other quant fixture has all six channels populated, so
  // fold_change is a real number and a wrong implementation agrees with a right one everywhere.
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_absent_scans = loadTsvScans(ms2_quant_absent_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_absent_scans.empty())

  // Control is wholly absent -> enriched in treated. Buys.
  const std::string in_treated =
      quantObjectiveJson("0.01", "", R"( "enriched_in": "treated",)");
  TEST_EQUAL(quantPushedFor(in_treated.c_str(), ms1_scans, ms2_absent_scans), 1)

  // ...and the other direction is correctly refused, so the assertion above is not just "the
  // absent-condition path always buys".
  const std::string in_control =
      quantObjectiveJson("0.01", "", R"( "enriched_in": "control",)");
  TEST_EQUAL(quantPushedFor(in_control.c_str(), ms1_scans, ms2_absent_scans), 0)

  // Control: with no direction restriction the same spectrum buys, confirming the verdict really
  // is Differential and the 0 above is the direction test rather than a rejected measurement.
  const std::string either = quantObjectiveJson("0.01", "", "");
  TEST_EQUAL(quantPushedFor(either.c_str(), ms1_scans, ms2_absent_scans), 1)

  // BOTH directions of absence, because one is not enough. With control absent, condition_means is
  // [0, X] and a bug spelled `fold_change < 1.0` gives the right answer by coincidence: -1 is
  // indeed below 1, so "enriched in conditions[1]" comes out true. Only the mirror -- treated
  // absent, condition_means [X, 0] -- separates them, because there the correct answer is
  // "enriched in conditions[0]" and every fold_change-based spelling of that reads the sentinel as
  // a tiny ratio and refuses.
  auto ms2_treated_absent_scans = loadTsvScans(ms2_quant_treated_absent_tsv_path);
  ABORT_IF(ms2_treated_absent_scans.empty())

  const std::string mirror_control =
      quantObjectiveJson("0.01", "", R"( "enriched_in": "control",)");
  TEST_EQUAL(quantPushedFor(mirror_control.c_str(), ms1_scans, ms2_treated_absent_scans), 1)

  const std::string mirror_treated =
      quantObjectiveJson("0.01", "", R"( "enriched_in": "treated",)");
  TEST_EQUAL(quantPushedFor(mirror_treated.c_str(), ms1_scans, ms2_treated_absent_scans), 0)
}
END_SECTION

START_SECTION(processScan_quant_followup)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_tmt_scans = loadTsvScans(ms2_tmt_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_tmt_scans.empty())
  // Use quant_sensitive_json with fold_change_threshold=0.01 to trigger on any reporter ratio
  FLASHIda* ida = new FLASHIda(const_cast<char*>(quant_sensitive_json));

  // Drive MS1 via engine-id echo with NO MS2 feed (the emitted MS2 commands stay pending; we feed the TMT
  // MS2 ourselves below to exercise the quant follow-up).
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  // ADR-0038: the ROSTERED MS2 is now the quantification scan -- marked 'Q', carrying
  // ms_settings.ms2_quant's parameters (ETD here). It is the scan that gets measured.
  const ScanCommand& ms2_cmd = acq.ms2_cmds[0];
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)
  TEST_EQUAL(ms2_cmd.scan_description[3], 'Q')
  TEST_STRING_EQUAL(std::string(ms2_cmd.stages[0].activation_type), "ETD")
  TEST_EQUAL(ms2_cmd.priority, 2)  // rostered, not a follow-up

  // Push its return with TMT reporter data
  const auto& ms2 = ms2_tmt_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
  // With fold_change_threshold=0.01, any reporter ion ratio is differential
  TEST_EQUAL(ms2_result, 1)  // exactly 1 identification scan bought per quantification scan

  // What a differential verdict BUYS is the identification scan: ms_settings.ms2, marked 'R', at
  // priority 0 so it dequeues ahead of the remaining rostered MS2s. Asserting HCD here is the whole
  // point of the fixture -- the two configs carry different activations, so this distinguishes
  // "bought ms_settings.ms2" from "re-emitted the screen".
  ScanCommand out{};
  int r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.scan_description[3], 'R')
  TEST_STRING_EQUAL(std::string(out.stages[0].activation_type), "HCD")
  TEST_EQUAL(out.stages[0].collision_energy, 29)
  TEST_EQUAL(out.msn_level, 2)
  TEST_EQUAL(out.priority, 0)  // bought, not rostered

  delete ida;
}
END_SECTION

// A conditional ('C') follow-up takes its fragmentation parameters from tagging.follow_up_scan,
// never from the MS2 that triggered it. Regression guard for the defect where buildFollowUp
// overrode the activation but left the activation-coupled parameters inherited, emitting a scan
// whose activation and reaction settings came from different configs.
START_SECTION(processScan_conditional_followup_owns_its_scan_parameters)
{
  {
    std::ifstream fasta_check(fasta_path);
    ABORT_IF(! fasta_check.good())
  }
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(followup_owns_params_conditional_json));

  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  // The triggering MS2 carries the OTHER set of values -- assert that first, so a failure
  // distinguishes "the follow-up inherited" from "the MS2 itself was misconfigured".
  const ScanCommand& ms2_cmd = acq.ms2_cmds[0];
  TEST_STRING_EQUAL(std::string(ms2_cmd.stages[0].activation_type), "HCD")
  TEST_REAL_SIMILAR(ms2_cmd.stages[0].reaction_time, 10.0)

  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                    (int)ms2.mzs.size(), ms2.rt,
                                    2, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
  TEST_EQUAL(ms2_result > 0, true)
  ABORT_IF(ms2_result <= 0)

  ScanCommand out{};
  int r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(out.priority, 0)
  TEST_STRING_EQUAL(std::string(out.stages[0].activation_type), "EThcD")

  // Every one of these differs from the triggering MS2's value, so each fails under the defect.
  TEST_REAL_SIMILAR(out.stages[0].reaction_time, 5.0)        // MS2 had 10.0
  TEST_REAL_SIMILAR(out.stages[0].reagent_max_it, 111.0)     // MS2 had 200.0
  TEST_EQUAL(out.stages[0].reagent_agc_target, 222)          // MS2 had 700000
  TEST_REAL_SIMILAR(out.stages[0].collision_energy, 7.0)     // MS2 had 29
  TEST_EQUAL(out.hcd_energy, 7)                              // MS2 had 29

  // Precursor context still comes from the triggering MS2 -- a follow-up of THAT precursor.
  TEST_REAL_SIMILAR(out.mono_mass, ms2_cmd.mono_mass)
  TEST_EQUAL(out.stages[0].charge_state, ms2_cmd.stages[0].charge_state)
  TEST_REAL_SIMILAR(out.faims_cv, ms2_cmd.faims_cv)

  delete ida;
}
END_SECTION

// Same contract on the quantification ('F') path: the two follow-up kinds share buildFollowUp but
// have independent config blocks, so both must be pinned.
START_SECTION(processScan_quant_followup_owns_its_scan_parameters)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_tmt_scans = loadTsvScans(ms2_tmt_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_tmt_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(followup_owns_params_quant_json));

  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  // ADR-0038 swapped which config is rostered. The rostered scan is now the QUANTIFICATION scan,
  // carrying ms_settings.ms2_quant (EThcD / rt 5.0 / reagent 111 / 222 / CE 7).
  const ScanCommand& ms2_cmd = acq.ms2_cmds[0];
  TEST_EQUAL(ms2_cmd.scan_description[3], 'Q')
  TEST_STRING_EQUAL(std::string(ms2_cmd.stages[0].activation_type), "EThcD")
  TEST_REAL_SIMILAR(ms2_cmd.stages[0].reaction_time, 5.0)

  const auto& ms2 = ms2_tmt_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                    (int)ms2.mzs.size(), ms2.rt,
                                    2, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
  TEST_EQUAL(ms2_result > 0, true)
  ABORT_IF(ms2_result <= 0)

  // The scan it buys is the IDENTIFICATION scan and takes every parameter from ms_settings.ms2,
  // never from the quantification scan that triggered it. The two configs disagree on all five
  // fields below, which is what makes this a test rather than a coincidence.
  ScanCommand out{};
  int r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(out.priority, 0)
  TEST_EQUAL(out.scan_description[3], 'R')
  TEST_STRING_EQUAL(std::string(out.stages[0].activation_type), "HCD")

  TEST_REAL_SIMILAR(out.stages[0].reaction_time, 10.0)       // the Q scan had 5.0
  TEST_REAL_SIMILAR(out.stages[0].reagent_max_it, 200.0)     // the Q scan had 111.0
  TEST_EQUAL(out.stages[0].reagent_agc_target, 700000)       // the Q scan had 222
  TEST_EQUAL(out.stages[0].collision_energy, 29)             // the Q scan had 7

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

  // Drive MS1 via engine-id echo with NO MS2 feed (the emitted MS2 commands stay pending; we feed one back
  // ourselves below to exercise the tag-targeting MS2 path).
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  const ScanCommand& ms2_cmd = acq.ms2_cmds[0];
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)

  // Push MS2 return (HCD fragment data)
  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
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

  // Bootstrap the survey + MS2-command seeding through the canonical harness instead of a hand-rolled
  // pull->echo->feed loop. runInterleaved drives MS1 by echoing each engine survey id (the always-on MS1
  // gate rejects fabricated ids) with NO MS2 feed (ms2_scans={}), so the selectable MS1's emitted MS2
  // commands stay QUEUED at priority 2 -- exactly the state the next cycle-time MS1 (priority 0) must beat.
  // acq.ms2_cmds non-empty == the old `fed_ms1` (a selectable MS1 produced MS2 commands).
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  // Guarantee a strictly positive gap since the last survey so cycle_time_ms=0 ("any elapsed time triggers")
  // fires deterministically rather than racing the sub-millisecond idle path.
  std::this_thread::sleep_for(std::chrono::milliseconds(2));

  // Next command: the elapsed time re-triggers a cycle-time MS1, queued at priority 0, dequeued before the
  // MS2s (priority 2). Same fields/values as before the migration -- only the MS1 feed mechanism changed.
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
  // agc_fast_json: agc_interval_seconds=0, so a scheduled prescan is due as soon as any measurable
  // time has passed. Step 1 is now the ONLY source of an AGC prescan (ADR-0031).
  // AGC is a dedicated minimal gain-control scan: makeAGC hardcodes agc_target=30000, max_it=1
  FLASHIda* ida = new FLASHIda(const_cast<char*>(agc_fast_json));

  // needsAGC() is `elapsed > agc_interval_ms`, so with the interval at 0 any non-zero elapsed
  // qualifies. Sleep past the millisecond boundary the steady_clock difference is measured in.
  // Before ADR-0031 this section passed either way -- the drained-queue path fabricated a prescan
  // regardless of whether the timer had fired, so it could not distinguish the two and never
  // proved Step 1 ran at all.
  std::this_thread::sleep_for(std::chrono::milliseconds(5));

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

  // Drive MS1 via engine-id echo with NO MS2 feed; feed one emitted MS2 back ourselves below.
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  const ScanCommand& ms2_cmd = acq.ms2_cmds[0];
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)

  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
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

  // Drive MS1 via engine-id echo with NO MS2 feed; feed one emitted MS2 back ourselves below.
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size() > 0, true)
  ABORT_IF(acq.ms2_cmds.empty())

  const ScanCommand& ms2_cmd = acq.ms2_cmds[0];
  TEST_EQUAL(std::strlen(ms2_cmd.scan_description) <= 15, true)
  TEST_EQUAL(ms2_cmd.msn_level, 2)
  TEST_EQUAL(ms2_cmd.priority, 2)

  // Push MS2 return — should find tags and trigger conditional follow-up
  const auto& ms2 = ms2_scans[0];
  int ms2_result = ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                     (int)ms2.mzs.size(), ms2.rt,
                                     2, ms2_cmd.scan_description, 0.0, instrumentScanNumberOf(ms2));
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

// Idle cycle: an empty queue produces a survey MS1 on EVERY call, never an AGC prescan.
//
// This section used to assert a strictly alternating AGC/MS1 pair, which was the drained-queue path
// fabricating a prescan as filler and pushing the survey behind it. Prescans are now scheduled by
// scheduling.agc_interval_seconds alone (ADR-0031), and idle_cycle_json pins that at 9999 s, so no
// prescan can appear here at all -- which is exactly what makes `is_agc == 0` a real assertion
// rather than a description of call parity.
START_SECTION(idle_cycle_returns_survey_ms1)
{
  FLASHIda* ida = new FLASHIda(const_cast<char*>(idle_cycle_json));

  // No processScan — the queue is empty from the start, so every call takes Step 5.
  ScanCommand cmd{};
  std::set<int> ids;

  for (int iter = 0; iter < 6; ++iter)
  {
    int r = ida->getNextScanCommand(cmd);
    TEST_EQUAL(r, 1)
    TEST_EQUAL(std::strlen(cmd.scan_description) <= 15, true)
    TEST_EQUAL(cmd.is_agc, 0)
    TEST_EQUAL(cmd.msn_level, 1)
    TEST_EQUAL(cmd.orbitrap_resolution, 120000)
    TEST_EQUAL(cmd.priority, 3)

    // Scan description: 3-char ID + type code 'S' for MS1 survey
    std::string ms1_desc(cmd.scan_description);
    TEST_EQUAL(ms1_desc.size() >= 4, true)
    TEST_EQUAL(ms1_desc[3], 'S')

    ids.insert(cmd.scan_id);
  }

  // Each idle survey is minted on the call that returns it, so six calls yield six distinct ids.
  // A regression that returned the same queued command twice would still satisfy every assertion
  // above.
  TEST_EQUAL(ids.size(), (size_t)6)

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
  FLASHIdaTestAccess::push(*ida,ms2_a);

  // First getNextScanCommand: MS2 at priority 2 (queue not empty, no idle cycle)
  ScanCommand out{};
  int r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.msn_level, 2)
  TEST_EQUAL(out.scan_id, 42)
  TEST_EQUAL(out.priority, 2)

  // The priority claim needs an MS1 that is genuinely WAITING at priority 3 while an MS2 sits at
  // priority 2. The idle survey no longer provides one: Step 5 pushes it and dequeues it within the
  // same call (ADR-0031), so it is never observably queued. Inject the pair directly instead --
  // which is a sharper test of the ladder anyway, since it does not depend on how the survey got
  // there.
  ScanCommand ms1_waiting{};
  ms1_waiting.msn_level = 1;
  ms1_waiting.priority = 3;
  ms1_waiting.scan_id = 44;
  FLASHIdaTestAccess::push(*ida, ms1_waiting);

  ScanCommand ms2_b{};
  ms2_b.msn_level = 2;
  ms2_b.priority = 2;
  ms2_b.scan_id = 43;
  ms2_b.faims_cv = -50.0;
  FLASHIdaTestAccess::push(*ida,ms2_b);

  // Second call: MS2 at priority 2 beats the MS1 waiting at priority 3, despite being pushed later.
  r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.msn_level, 2)
  TEST_EQUAL(out.scan_id, 43)
  TEST_EQUAL(out.priority, 2)

  // Third call: the MS1 at priority 3 is dequeued.
  r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(std::strlen(out.scan_description) <= 15, true)
  TEST_EQUAL(out.is_agc, 0)
  TEST_EQUAL(out.msn_level, 1)
  TEST_EQUAL(out.scan_id, 44)
  TEST_EQUAL(out.priority, 3)

  // Fourth call: queue drained -> a fresh idle survey, minted here rather than left over.
  r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(out.is_agc, 0)
  TEST_EQUAL(out.msn_level, 1)
  TEST_EQUAL(out.priority, 3)
  TEST_NOT_EQUAL(out.scan_id, 44)

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

  // Drive MS1 via engine-id echo for each selection metric; the collected MS2 commands are exactly the
  // non-AGC commands the old drain-until-first-AGC loop counted.
  AcqResult acq_q = runInterleaved(ida_qscore, ms1_scans, {});
  AcqResult acq_i = runInterleaved(ida_intensity, ms1_scans, {});

  TEST_EQUAL(acq_q.ms2_cmds.size() > 0, true)
  TEST_EQUAL(acq_i.ms2_cmds.size() > 0, true)

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

  // Drive MS1 via engine-id echo. selection=none -> the MS1 path selects nothing and pushes 0 MS2 commands
  // (FLASHIda.cpp:789-792), so no MS2 commands are ever emitted.
  AcqResult acq = runInterleaved(ida, ms1_scans, {});
  TEST_EQUAL(acq.ms2_cmds.size(), (size_t)0)

  // After the drive the queue is drained, so every call yields a fresh idle survey MS1 -- no AGC
  // prescan is interleaved with them any more (ADR-0031).
  ScanCommand cmd{};
  TEST_EQUAL(ida->getNextScanCommand(cmd), true)  // idle cycle always returns 1
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_EQUAL(cmd.is_agc, 0)
  TEST_EQUAL(cmd.priority, 3)
  TEST_EQUAL(ida->getNextScanCommand(cmd), true)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_EQUAL(cmd.is_agc, 0)
  TEST_EQUAL(cmd.priority, 3)

  delete ida;
}
END_SECTION

// max_targets caps precursor selection PER MS1 scan: children (MS2 commands) per scan rise with the cap,
// and equal (precursors selected) x (MS2 activations per precursor). Verified via the engine's own
// per-scan logging (runtime.log_dir -> scan_results.tsv commands_pushed / child_ids), across a 2 (activation) x
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
  // MS1 is driven via the engine-id-echo contract (runInterleaved): the always-on MS1 gate rejects the old
  // fabricated "scan_"+id feed, so we feed each survey with the engine's own scan_description. Only MS1
  // scans are fed (ms2_scans={} -> MS2 commands are collected but never fed back), so every data row in the
  // results TSV is an MS1 result row, in scan order -- exactly as readMs1ResultRows assumes. Per-scan
  // selection is unchanged by interleaving (same spectra, same order, same accumulated exclusion state).
  auto runCase = [&](int max_targets, bool etd) -> std::vector<int> {
    // Six engines run through this lambda, so each needs its OWN folder: the basenames are fixed
    // now, and the streams still open in append mode, so a shared folder would concatenate runs
    // and inject a second header row mid-file. freshLogDir wipes and creates -- the engine never
    // creates a directory, and a missing one would leave every stream closed and this test
    // reading zero rows without any error.
    const std::string dir = freshLogDir(std::string("ps_maxcap_") + (etd ? "etd" : "hcd") + "_"
                                        + std::to_string(max_targets));
    std::string cfg = buildCapConfig(max_targets, etd, dir);
    FLASHIda* ida = new FLASHIda(const_cast<char*>(cfg.c_str()));
    runInterleaved(ida, ms1_scans, {});
    delete ida;  // flush + close the results TSV
    std::vector<int> children;
    for (const auto& r : readMs1ResultRows(dir + "/scan_results.tsv"))
    { TEST_EQUAL(r.first, r.second) children.push_back(r.first); }
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
    "min_precursor_charge": 99
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
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

  FLASHIda ida(const_cast<char*>(min_charge_json));

  // Load same MS1 scans that normally produce commands
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  // Drive MS1 via engine-id echo (the always-on MS1 gate rejects fabricated "scan_"+id ids; feeding the
  // engine's own survey id lets the spectrum actually reach deconvolution/selection so the min_charge=99
  // FILTER -- not the gate -- is what produces zero MS2 commands, preserving this test's intent).
  AcqResult acq = runInterleaved(&ida, scans, {});

  // With min_charge=99, no precursor passes the filter, so no deconvolution-derived MS2 is produced.
  TEST_EQUAL(acq.ms2_cmds.size(), (size_t)0)

  // getNextScanCommand has no return-0 path: after the drive the idle cycle's pending priority-3 survey MS1 is
  // dequeued first (is_agc==0), then the idle AGC (cf. processScan_commands_dequeued).
  ScanCommand cmd{};
  int result = ida.getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_EQUAL(cmd.is_agc, 0)
  TEST_EQUAL(cmd.priority, 3)
  result = ida.getNextScanCommand(cmd);
  TEST_EQUAL(result, 1)
  TEST_EQUAL(cmd.msn_level, 1)
  TEST_EQUAL(cmd.is_agc, 0)   // a drained queue yields surveys, never a prescan (ADR-0031)
  TEST_EQUAL(cmd.priority, 3)
}
END_SECTION

START_SECTION(processScan_agc_scan_skipped)
{
  // Use agc_fast_json (agc_interval_seconds=0) to get an AGC command immediately
  FLASHIda* ida = new FLASHIda(const_cast<char*>(agc_fast_json));

  // Get a scheduled AGC prescan from the engine. agc_fast_json sets agc_interval_seconds=0, and the
  // sleep carries `elapsed` past the millisecond boundary needsAGC() measures in -- Step 1 is the
  // only path that emits a prescan now (ADR-0031).
  std::this_thread::sleep_for(std::chrono::milliseconds(5));

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
                                 agc_cmd.scan_description, 0.0, instrumentScanNumberOf(scan));
  TEST_EQUAL(result, 0)

  // Verify no commands were generated
  ScanCommand next_cmd{};
  int has_cmd = ida->getNextScanCommand(next_cmd);
  // If the gate works, the scan data was never deconvolved, so nothing MS2 was queued and the next
  // drain can only be an idle survey or another scheduled prescan -- both MS1, both stage-less.
  //
  // The old form was `has_cmd == 0 || is_agc == 1 || msn_level == 1`, which is vacuous: an AGC is
  // itself msn_level 1, so the disjunction reduced to "msn_level == 1" and getNextScanCommand never
  // returns 0. Assert the absence of MS2 directly instead.
  TEST_EQUAL(has_cmd, 1)
  TEST_EQUAL(next_cmd.msn_level, 1)
  TEST_EQUAL(next_cmd.num_stages, 0)

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
    "max_cv_skip": 0
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": true,
      "value_ms": 1
    },
    "agc_interval_seconds": 9999
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
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
    "max_targets": 10
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

  FLASHIda* ida = new FLASHIda(const_cast<char*>(short_timeout_json));

  // Push 3 MS2 commands at priority 2 via test helper
  for (int i = 0; i < 3; ++i)
  {
    ScanCommand ms2{};
    ms2.msn_level = 2;
    ms2.priority = 2;
    ms2.scan_id = FLASHIdaTestAccess::queue(*ida).nextTrackingId();
    FLASHIdaTestAccess::push(*ida,ms2);
  }

  // Verify all 3 are in the queue
  TEST_EQUAL(FLASHIdaTestAccess::queueSize(*ida,2), (size_t)3)

  // Sleep 2ms to ensure commands exceed the 1ms timeout
  std::this_thread::sleep_for(std::chrono::milliseconds(2));

  // getNextScanCommand calls cleanupExpired in step 3.
  // After cleanup, the queue should be empty — all 3 stale commands dropped.
  // The idle cycle (step 5) will fire and return a fresh survey MS1.
  ScanCommand out{};
  int r = ida->getNextScanCommand(out);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(out.msn_level, 1)  // idle survey, not a stale MS2
  TEST_EQUAL(out.is_agc, 0)
  TEST_EQUAL(out.priority, 3)

  // Queue at priority 2 should now be empty (stale commands dropped)
  TEST_EQUAL(FLASHIdaTestAccess::queueSize(*ida,2), (size_t)0)

  delete ida;
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0031: the scheduled prescan is the ONLY prescan, and priority 3 is the drained-queue signal.
// Both facts are now load-bearing and neither was pinned before.
/////////////////////////////////////////////////////////////

// The interval AGC must fire on its schedule -- roughly once per interval, not never and not on
// every drain.
//
// Why this exists: with the drained-queue prescan deleted, agc_interval_seconds is the only thing
// that emits a prescan at all. A broken gate means the instrument acquires an entire run with no
// flux estimate, and every scan is gain-corrected from whatever the startup handshake left behind.
// Nothing else would notice -- the same silent shape as the faims cv_values `size() > 1` bug and
// the source_cid_scaling `> 0` guard.
START_SECTION(interval_agc_fires_once_per_interval)
{
  nlohmann::json j = nlohmann::json::parse(standard_json);
  j["scheduling"]["agc_interval_seconds"] = 0.05;   // 50 ms
  std::string agc_periodic_json = j.dump();

  FLASHIda ida(const_cast<char*>(agc_periodic_json.c_str()));

  int agc = 0, other = 0;
  for (int i = 0; i < 20; ++i)
  {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));
    ScanCommand c {};
    int r = ida.getNextScanCommand(c);
    TEST_EQUAL(r, 1)
    if (c.is_agc) ++agc; else ++other;
  }

  std::cout << "[AGC-CADENCE] interval_ms=50 drains=20 agc=" << agc << " other=" << other << std::endl;

  // ~200 ms of sleep against a 50 ms interval. The bounds are deliberately loose and asymmetric: a
  // loaded runner only ever sleeps LONGER, which raises `agc`, so the upper bound carries the
  // slack and there is no flake direction that a slow machine can push this into.
  TEST_EQUAL(agc >= 2, true)
  TEST_EQUAL(agc <= 12, true)
  // Anti-vacuous: idle surveys still flow between prescans. A gate stuck open would give agc == 20.
  TEST_EQUAL(other > 0, true)
}
END_SECTION

// Nothing but the idle survey is emitted at priority 3.
//
// Why this exists: three C# drain loops (FLASHIdaWrapper's offline harness x2, ContinuityTestHarness
// .PushScan) terminate on `msn_level == 1 && priority == 3`. If any other command were ever emitted
// at priority 3 those loops would stop early, and the failure is invisible -- a truncated drain
// looks exactly like "the engine had nothing more to do".
//
// Scoped to ENGINE-EMITTED commands via runInterleaved. Tests may legitimately inject a priority-3
// command through FLASHIdaTestAccess::push (queue_priority_dequeue_order does); that bypasses the
// builders and says nothing about what production emits.
START_SECTION(only_the_idle_survey_is_emitted_at_priority_3)
{
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  auto ms2_scans = loadTsvScans(ms2_tsv_path);
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())

  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));
  AcqResult a = runInterleaved(ida, ms1_scans, std::vector<ScanData>{ms2_scans[0]});
  delete ida;

  int p3 = 0, ms2_seen = 0;
  for (const auto& c : a.all_cmds)
  {
    if (c.msn_level == 2) ++ms2_seen;
    if (c.priority != 3) continue;
    ++p3;
    TEST_EQUAL(c.msn_level, 1)
    TEST_EQUAL(c.is_agc, 0)
  }

  std::cout << "[P3-CONTRACT] commands=" << a.all_cmds.size() << " p3=" << p3
            << " ms2=" << ms2_seen << std::endl;

  // Anti-vacuous twice over: the drive must have produced real MS2 workload (so the loop had
  // non-p3 commands to reject) AND at least one priority-3 command (so the assertions above ran).
  TEST_EQUAL(ms2_seen > 0, true)
  TEST_EQUAL(p3 > 0, true)
}
END_SECTION

// Tracking id 0 is RESERVED and must never be issued.
//
// Why this exists: 0 doubles as buildMS2's "no parent / root" sentinel (`if (parent_scan_id > 0)`,
// ScanCommandQueue.cpp), and FLASHIda.cpp passes a literal 0 for a root MS2. So a survey that
// actually HOLDS id 0 is indistinguishable from one with no parent -- every MS2 it parents ships
// with an empty parent_tracking_id and the lineage graph loses its first generation, on the
// instrument as much as in tests.
//
// This was latent for the whole life of the code. The drained-queue AGC prescan was the first
// command minted on every fresh engine, so it absorbed id 0 and no survey ever held it. Deleting
// that prescan (ADR-0031) handed 0 straight to the first survey, and
// FLASHIda_LoggingFields_test::parent_tracking_id_resolution caught it. A drift guard belongs here
// rather than there, because the defect is in the id allocator, not in logging.
START_SECTION(tracking_id_zero_is_never_issued)
{
  FLASHIda ida(const_cast<char*>(standard_json));

  // Fresh engine: the very first id handed out must not be the sentinel.
  const int first = ida.getNextTrackingId();
  std::cout << "[ID-RESERVE] first_id=" << first << std::endl;
  TEST_NOT_EQUAL(first, 0)

  // ...and across the wrap. The counter resets when it passes 94^3 - 1; resetting it to 0 rather
  // than 1 would reintroduce the defect once per 830k scans -- rare enough to never be noticed and
  // permanent once it happened, since the id would then recur every wrap.
  // No early exit: issuing max_id + 5 ids necessarily crosses the wrap, so the count below is
  // anti-vacuous by construction and one defect produces exactly one failure.
  const int max_id = 94 * 94 * 94 - 1;
  int zeros = (first == 0) ? 1 : 0;
  for (int i = 0; i < max_id + 4; ++i)
    if (ida.getNextTrackingId() == 0) ++zeros;

  std::cout << "[ID-RESERVE] issued=" << (max_id + 5) << " zeros=" << zeros << std::endl;
  TEST_EQUAL(zeros, 0)
}
END_SECTION

// MS1 and AGC scans should be resolved from pending_scan_map_ after processScan
START_SECTION(ms1_agc_resolved_from_pending_map)
{
  // agc_fast_json, not idle_cycle_json: this section needs an AGC prescan in the pending map, and
  // Step 1 is now the only thing that puts one there (ADR-0031). idle_cycle_json pins the interval
  // at 9999 s, so under it this section would never see a prescan at all.
  FLASHIda* ida = new FLASHIda(const_cast<char*>(agc_fast_json));

  // Interval 0 + a sleep past the millisecond boundary -> the first drain takes Step 1.
  std::this_thread::sleep_for(std::chrono::milliseconds(5));

  ScanCommand agc_cmd{};
  int r = ida->getNextScanCommand(agc_cmd);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(agc_cmd.is_agc, 1)

  // AGC is in pending map via registerPending
  TEST_EQUAL(FLASHIdaTestAccess::queue(*ida).pendingScanMapSize(), (size_t)1)

  // Next drain: Step 1 just called recordAGCTime(), so `elapsed > 0` is false again and the drained
  // queue takes Step 5 -- an idle survey MS1, registered pending by dequeue(). Both are now held.
  ScanCommand ms1_cmd{};
  r = ida->getNextScanCommand(ms1_cmd);
  TEST_EQUAL(r, 1)
  TEST_EQUAL(ms1_cmd.is_agc, 0)
  TEST_EQUAL(ms1_cmd.msn_level, 1)
  TEST_EQUAL(ms1_cmd.priority, 3)
  TEST_EQUAL(FLASHIdaTestAccess::queue(*ida).pendingScanMapSize(), (size_t)2)

  // processScan with AGC scan description — should resolve AGC from pending map
  // AGC gate (desc[3]=='A') returns 0 and resolves the pending entry.
  int n = ida->processScan(nullptr, nullptr, 0, 0.0, 1, agc_cmd.scan_description, 0.0, 0);
  TEST_EQUAL(n, 0)
  TEST_EQUAL(FLASHIdaTestAccess::queue(*ida).pendingScanMapSize(), (size_t)1)

  // processScan with MS1 scan description — should resolve MS1 from pending map
  // Load real MS1 data to feed through the MS1 path
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  const auto& scan = ms1_scans[0];
  n = ida->processScan(scan.mzs.data(), scan.ints.data(),
                       (int)scan.mzs.size(), scan.rt, 1,
                       ms1_cmd.scan_description, 0.0, instrumentScanNumberOf(scan));
  // n >= 0 (may or may not produce MS2 commands from a single scan)
  TEST_EQUAL(n >= 0, true)
  TEST_EQUAL(FLASHIdaTestAccess::queue(*ida).pendingScanMapSize(), (size_t)0)

  delete ida;
}
END_SECTION

// ADDITIVE: pins the always-on MS1 gate (FLASHIda.cpp:775) — an MS1 fed with a tracking id that was never
// emitted as a survey command is rejected (return 0), an MS1 fed with the engine's OWN survey id is accepted
// (selectable -> >=1 MS2), and that resolved id cannot be reused. This is the engine-side enforcement that
// FORCES the interleaved engine-id-echo harness contract (docs/kb/test-harness/README.md, invariant 6): a
// test harness cannot fabricate MS1 ids; it must pull the survey command first, then echo its scan_description.
START_SECTION(processScan_ms1_gate_rejects_unrequested_id)
{
  // Same selectable MS1 fixture + selecting config as processScan_ms1_returns_commands.
  auto ms1_scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(ms1_scans.empty())
  FLASHIda* ida = new FLASHIda(const_cast<char*>(standard_json));

  // ---- NEGATIVE: a never-emitted tracking id is rejected by the gate ----
  // "ZZZ" is a syntactically valid 3-char base-94 id, 'S' is the MS1-survey type code, but this id was never
  // minted/emitted as a command, so it is not in pending_scan_map_ -> the gate returns 0 before any selection.
  const size_t pending_before = FLASHIdaTestAccess::queue(*ida).pendingScanMapSize();
  const auto& neg = ms1_scans[0];
  int ret_neg = ida->processScan(neg.mzs.data(), neg.ints.data(), (int)neg.mzs.size(), neg.rt, 1, "ZZZS", 0.0, instrumentScanNumberOf(neg));
  TEST_EQUAL(ret_neg, 0)
  // No MS2 was queued and no pending entry was added: the spectrum never reached deconvolution/selection.
  TEST_EQUAL(FLASHIdaTestAccess::queue(*ida).pendingScanMapSize(), pending_before)

  // ---- POSITIVE: an engine-emitted survey id + a SELECTABLE spectrum is accepted ----
  // Drain a FRESH survey per attempt and feed successive scans until one selects a precursor (scan 0 of
  // ms1_standard is not selectable; a feed resolves its survey id, so each attempt needs its own survey). The
  // idle cycle alternates AGC (skip) then a priority-3 survey MS1 whose minted id IS in pending_scan_map_ at
  // dequeue, so the gate passes and selection runs.
  int ret_pos = 0;
  ScanCommand used{};
  for (int si = 0; si < (int)ms1_scans.size() && ret_pos == 0; ++si)
  {
    ScanCommand survey{};
    bool got_survey = false;
    for (int k = 0; k < 8 && !got_survey; ++k)
    {
      int rs = ida->getNextScanCommand(survey);
      TEST_EQUAL(rs, 1)
      if (survey.is_agc || survey.scan_description[0] == '\0' || survey.msn_level != 1) continue;
      got_survey = true;
    }
    ABORT_IF(!got_survey)
    TEST_EQUAL(std::strlen(survey.scan_description) <= 15, true)
    TEST_EQUAL(survey.is_agc, 0)
    const auto& pos = ms1_scans[si];
    int r = ida->processScan(pos.mzs.data(), pos.ints.data(), (int)pos.mzs.size(), pos.rt, 1,
                             survey.scan_description, 0.0, instrumentScanNumberOf(pos));
    if (r >= 1) { ret_pos = r; used = survey; }
  }
  TEST_EQUAL(ret_pos >= 1, true)  // engine-emitted id + selectable MS1 -> at least one MS2 command pushed

  // ---- NO-REUSE: the now-resolved survey id is rejected on a second feed ----
  // processScan resolved (erased) the survey id from pending_scan_map_, so a repeat feed with the SAME id is
  // no longer "emitted" -> the gate rejects it (return 0), symmetric with the MS2/MS3 resolve-once gate.
  const auto& reuse_scan = ms1_scans[0];
  int ret_reuse = ida->processScan(reuse_scan.mzs.data(), reuse_scan.ints.data(), (int)reuse_scan.mzs.size(),
                                   reuse_scan.rt, 1, used.scan_description, 0.0, instrumentScanNumberOf(reuse_scan));
  TEST_EQUAL(ret_reuse, 0)

  delete ida;
}
END_SECTION

/////////////////////////////////////////////////////////////
// Drain-blocking pins.
//
// These two are a matched pair and only mean something together. The first says the drain must NOT
// wait on analysis_mutex_; the second says processScan must STILL take it. Either one alone can be
// satisfied by deleting the wrong lock.
/////////////////////////////////////////////////////////////

// The drain must not acquire analysis_mutex_ on ANY of its three emitting paths.
START_SECTION(drain_completes_while_analysis_mutex_held)
{
  // processScan takes analysis_mutex_ function-scoped (FLASHIda.cpp:84) across its whole body,
  // including the MS1 deconvolution. Holding that same mutex here is an exact stand-in for "a
  // deconvolution is in flight" -- there is no timing threshold and no scheduling assumption to
  // flake on: either the drain wants the lock or it does not.
  //
  // THREE PASSES, all load-bearing. getNextScanCommand has three emitting paths and each acquired
  // the mutex at its own site, so covering one proves nothing about the other two:
  //   Step 1 (:667) time-triggered AGC -- reachable only when needsAGC() is true
  //   Step 4 (:708) priority dequeue   -- reachable only with something queued
  //   Step 5 (:758) idle fallback      -- what a fresh engine does
  // needsAGC() is FALSE on a fresh engine: last_agc_time_ is stamped in the ScanCommandQueue ctor
  // and standard_json pins agc_interval_seconds at 9999999. So Step 1 needs its own config with the
  // interval at zero, or that guard ships unpinned.
  // The helper deliberately returns a status instead of asserting: the ClassTest macros belong at
  // section scope, and keeping them out of the lambda also keeps the lock's lifetime obvious.
  auto drainWhileHeld = [](FLASHIda& ida, int& rc_out) -> std::future_status {
    std::unique_lock<std::mutex> lk(FLASHIdaTestAccess::analysisMutex(ida));

    ScanCommand cmd {};
    // std::launch::async EXPLICITLY. The default policy is allowed to defer, and a deferred future
    // reports future_status::deferred rather than ready -- which would fail this test against
    // perfectly correct code.
    std::future<int> fut = std::async(std::launch::async, [&ida, &cmd]() { return ida.getNextScanCommand(cmd); });

    std::future_status st = fut.wait_for(std::chrono::seconds(2));

    // Release BEFORE touching the future. If the drain is still blocked (the pre-fix state) then
    // get() -- or just the future's destructor -- would wait on it while we still hold the very lock
    // it is waiting for. That is a real deadlock, and ctest would only break it at its 1500 s
    // default timeout. Unlocking first turns the pre-fix state into a clean assertion failure.
    lk.unlock();

    rc_out = fut.get();
    return st;
  };

  // --- Pass 1: fresh engine, empty queue -> Step 5, the idle fallback ---
  {
    FLASHIda ida(const_cast<char*>(standard_json));
    int rc = 0;
    std::future_status st = drainWhileHeld(ida, rc);
    std::cout << "[DRAIN-PIN] path=step5_idle completed_while_held="
              << (st == std::future_status::ready ? 1 : 0) << std::endl;
    TEST_EQUAL(st == std::future_status::ready, true)
    TEST_EQUAL(rc, 1)
  }

  // --- Pass 2: something queued -> Step 4, the priority dequeue (the one real map read) ---
  {
    FLASHIda ida(const_cast<char*>(standard_json));
    ScanCommand queued {};
    queued.msn_level = 2;
    queued.priority = 2;
    queued.scan_id = 7;
    FLASHIdaTestAccess::push(ida, queued);

    int rc = 0;
    std::future_status st = drainWhileHeld(ida, rc);
    std::cout << "[DRAIN-PIN] path=step4_dequeue completed_while_held="
              << (st == std::future_status::ready ? 1 : 0) << std::endl;
    TEST_EQUAL(st == std::future_status::ready, true)
    TEST_EQUAL(rc, 1)
  }

  // --- Pass 3: AGC due -> Step 1, which returns before the dequeue is even reached ---
  {
    nlohmann::json j = nlohmann::json::parse(standard_json);
    j["scheduling"]["agc_interval_seconds"] = 0;
    std::string agc_due_json = j.dump();

    FLASHIda ida(const_cast<char*>(agc_due_json.c_str()));
    // needsAGC() is `elapsed > agc_interval_ms`, so with the interval at 0 any non-zero elapsed
    // qualifies. Sleep past the millisecond boundary the steady_clock difference is measured in.
    std::this_thread::sleep_for(std::chrono::milliseconds(5));

    int rc = 0;
    std::future_status st = drainWhileHeld(ida, rc);
    std::cout << "[DRAIN-PIN] path=step1_agc completed_while_held="
              << (st == std::future_status::ready ? 1 : 0) << std::endl;
    TEST_EQUAL(st == std::future_status::ready, true)
    TEST_EQUAL(rc, 1)
  }
}
END_SECTION

// processScan must STILL hold analysis_mutex_ -- the inverse pin.
START_SECTION(process_scan_still_blocks_while_analysis_mutex_held)
{
  // Why this exists: the cheapest way to make the test above pass is to delete the lock at
  // FLASHIda.cpp:84 rather than the three drain-side guards. That "fix" leaves deconv_, selection_,
  // exploration_, tracker_ and the rest unsynchronised against a concurrent drain, and every other
  // test in this suite is single-threaded so none of them would notice. This one fails.
  //
  // One-sided by construction: it asserts something does NOT happen inside a window, so a slow or
  // loaded runner can only make it more likely to pass. There is no flake direction.
  FLASHIda ida(const_cast<char*>(standard_json));

  std::unique_lock<std::mutex> lk(FLASHIdaTestAccess::analysisMutex(ida));

  // The lock is the FIRST statement of processScan, so the arguments are irrelevant -- this call
  // parks before it can look at any of them. Empty peak arrays keep it that way if it ever does run.
  std::vector<double> no_mzs;
  std::vector<double> no_ints;
  std::future<int> fut = std::async(std::launch::async, [&ida, &no_mzs, &no_ints]() {
    return ida.processScan(no_mzs.data(), no_ints.data(), 0, 0.0, 1, "", 0.0, 0);
  });

  std::future_status st = fut.wait_for(std::chrono::milliseconds(500));

  // Release first, for the same deadlock reason as above -- here the future is EXPECTED to be
  // outstanding, so this unlock is what lets it finish.
  lk.unlock();

  std::cout << "[DRAIN-PIN] path=process_scan blocked_while_held="
            << (st == std::future_status::timeout ? 1 : 0) << std::endl;
  TEST_EQUAL(st == std::future_status::timeout, true)
  TEST_EQUAL(fut.get(), 0)  // empty description -> the size<3 gate, once it finally runs
}
END_SECTION

END_TEST
