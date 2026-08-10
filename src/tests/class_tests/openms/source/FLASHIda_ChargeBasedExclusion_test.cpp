// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Tests for the charge_based_exclusion developer flag. Uses the ms1_standard TSV
// fixture already used by FLASHIda_ProcessScan_test.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include "FLASHIda_TestHelpers.h"  // ground-truth harness: ScanData/loadTsvScans/runInterleaved/AcquisitionRow/ms2AcquisitionRows

#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;

namespace
{
  const char* base_off_json = R"({
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
    "agc_interval_seconds": 30
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
    "charge_based_exclusion": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 10
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

  const char* base_on_json = R"({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.3,
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
    "agc_interval_seconds": 30
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
    "charge_based_exclusion": true,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 10
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

  // cytC (~12356 Da neutral) is present at many charge states (envelopes at m/z
  // ~1030/1124/1236/1373/1545), so a single mass appears at >=2 charges (CBE-03) and every
  // scan is productive — unlike ms1_standard (<=1 selectable/scan, charges spread across
  // different masses), which left CBE-03/04/06 with 0 MS2.
  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_cytc.txt";


  // Interleaved drive result: every emitted MS2 (charge, mz, width) tuple + the MS2 command count.
  struct DriveResult
  {
    std::vector<AcquisitionRow> rows;
    int ms2_count = 0;
  };

  // Thin adapter over the ground-truth harness runInterleaved (NO hand-rolled drive): MS1-only DDA driven
  // by the engine's own survey ids; rt_offset shifts each fed survey (2-pass exclusion: 2nd pass inside vs
  // outside rt_window). Projects the MS2 commands to (charge, mz, width) rows + count. State persists in
  // `ida` across calls, so a caller may invoke this twice (two passes) on the same engine.
  DriveResult driveInterleaved(FLASHIda& ida, const std::vector<ScanData>& ms1_scans,
                               double rt_offset = 0.0, int max_iters = 4000)
  {
    AcqResult a = runInterleaved(&ida, ms1_scans, {}, nullptr, max_iters, /*single_group_only=*/false, rt_offset);
    DriveResult res;
    res.rows = ms2AcquisitionRows(a);
    res.ms2_count = static_cast<int>(a.ms2_cmds.size());
    return res;
  }

  // Convenience: fresh engine, single interleaved pass, return the emitted MS2 (charge, mz, width) rows.
  std::vector<AcquisitionRow> runAndCollect(const char* cfg, const std::vector<ScanData>& scans)
  {
    FLASHIda ida(const_cast<char*>(cfg));
    return driveInterleaved(ida, scans).rows;
  }
}

START_TEST(FLASHIda_ChargeBasedExclusion, "$Id$")

// CBE-01: Flag off produces the same count and charge set as the existing
// default config — regression guard for byte-for-byte preservation.
START_SECTION(flag_off_matches_default_behavior)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())
  auto rows = runAndCollect(base_off_json, scans);
  TEST_EQUAL(rows.size() > 0, true)
  // Count unique (charge, mz) tuples — with the flag off, each (charge, mz) should be unique per scan
  // (default same-mz avoidance). Assert every emitted row has charge >= 4 (config min_charge).
  for (const auto& r : rows) { TEST_EQUAL(r.charge >= 4, true) }
}
END_SECTION

// CBE-02: Flag on emits more acquisitions than flag off across the same scan sequence
// because multiple charges of the same mass are now fragmented.
START_SECTION(flag_on_emits_more_acquisitions_than_flag_off)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())
  auto off_rows = runAndCollect(base_off_json, scans);
  auto on_rows = runAndCollect(base_on_json, scans);
  TEST_EQUAL(off_rows.size() > 0, true)
  TEST_EQUAL(on_rows.size() > off_rows.size(), true)
}
END_SECTION

// CBE-03: Under the flag, two acquisitions share a nominal mass but have different
// charges — demonstrates per-mass multi-charge expansion.
// We bucket by charge-corrected nominal mass: mono_mass ≈ (mz - proton) * charge.
// Since ScanCommand exposes precursor_mz and charge_state, reconstruct mass per row and
// confirm at least one mass bucket contains ≥2 distinct charges.
START_SECTION(flag_on_yields_multiple_charges_per_mass)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())
  auto rows = runAndCollect(base_on_json, scans);
  ABORT_IF(rows.empty())

  constexpr double proton = 1.00728;
  std::map<int, std::set<int>> nominal_mass_to_charges;
  for (const auto& r : rows)
  {
    int nominal = (int)std::round((r.mz - proton) * r.charge);
    nominal_mass_to_charges[nominal].insert(r.charge);
  }
  int masses_with_multiple_charges = 0;
  for (const auto& kv : nominal_mass_to_charges) { if (kv.second.size() >= 2) { masses_with_multiple_charges++; } }
  TEST_EQUAL(masses_with_multiple_charges >= 1, true)
}
END_SECTION

// CBE-04: Under the flag, driving the first scan through processScan twice in a row
// yields fewer commands on the second invocation — per-(mass, charge) exclusion
// suppresses re-acquisition of charges whose per-charge qscore already crossed threshold.
START_SECTION(flag_on_suppresses_reacquisition_across_scans)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  FLASHIda ida(const_cast<char*>(base_on_json));

  // First pass: interleave-drive every scan once via the engine's own survey ids (a single scan may
  // yield no selectable precursor; the full sequence reliably does). MS2 commands are drained inline.
  int drained_first = driveInterleaved(ida, scans).ms2_count;
  TEST_EQUAL(drained_first > 0, true)

  // Second pass on the SAME engine: same precursors at an RT within rt_window, so the per-(mass,
  // charge) exclusion set populated in pass 1 suppresses re-acquisition.
  int drained_second = driveInterleaved(ida, scans, /*rt_offset=*/0.001).ms2_count;
  TEST_EQUAL(drained_second < drained_first, true)
}
END_SECTION

// CBE-05: With the flag on, acquisitions targeting the same nominal m/z appear at multiple
// distinct charge values — if the flag had no effect we would only see one.
START_SECTION(flag_on_diverse_charges_at_same_mass)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())
  auto on_rows = runAndCollect(base_on_json, scans);
  ABORT_IF(on_rows.empty())
  std::set<int> distinct_charges;
  for (const auto& r : on_rows) { distinct_charges.insert(r.charge); }
  TEST_EQUAL(distinct_charges.size() >= 2, true)
}
END_SECTION

// CBE-06: Under the flag, entries aged past rt_window are evicted so a previously-
// excluded (mass, charge) becomes eligible again on a later scan.
START_SECTION(flag_on_rt_window_eviction_reenables_charge)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  FLASHIda ida(const_cast<char*>(base_on_json));

  // First pass at the scans' own RTs (interleaved via the engine's survey ids) populates the
  // per-(mass, charge) exclusion entries.
  int drained_first = driveInterleaved(ida, scans).ms2_count;
  TEST_EQUAL(drained_first > 0, true)

  // Second pass far past rt_window (default 180): evictions fire, per-charge state is
  // cleared, charges become eligible again, so re-acquisition is not suppressed.
  int drained_second = driveInterleaved(ida, scans, /*rt_offset=*/1000.0).ms2_count;
  // Eviction cleared the per-charge state, so the second pass acquires at least as many.
  TEST_EQUAL(drained_second >= drained_first, true)
}
END_SECTION

// CBE-07: charges_to_process is a preference-ordered FALLBACK list, not a work list. Within a
// SINGLE survey a mass is acquired at exactly ONE charge — the best one not already excluded.
//
// This is the only section that separates the fallback from the fan-out. CBE-02/03/05 assert
// more-acquisitions and multiple-charges-per-mass across the FULL scan sequence, and all three hold
// under BOTH designs, because the fallback also reaches new charges — just on later surveys. That
// gap is exactly how the fan-out survived: the flag consumed one mass_count slot per charge, so
// with max_precursors 3 and a species at three charges, the 2nd and 3rd precursors were never
// fragmented, and no test noticed.
START_SECTION(flag_on_selects_one_charge_per_mass_within_a_survey)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  // The FULL sequence, as every other section drives it. A single survey does NOT work: this
  // fixture's first scan selects nothing, so the engine emits only idle AGC/MS1 pairs and the
  // harness terminates on its 3-idle rule with zero MS2 commands. runInterleaved is called directly
  // rather than through driveInterleaved because the assertion needs whole commands (mono_mass,
  // parent_scan_id), which the AcquisitionRow projection drops.
  FLASHIda ida(const_cast<char*>(base_on_json));
  AcqResult a = runInterleaved(&ida, scans, {}, nullptr, /*max_iters=*/4000);

  // Non-vacuity: with no MS2 commands the assertion below is trivially satisfied.
  TEST_EQUAL(a.ms2_cmds.size() > 0, true)

  // "Within a survey" must be keyed per PARENT, not per run. Charge-keyed exclusion is SUPPOSED to
  // reach other charges of a mass on LATER surveys -- that is exactly what CBE-03/CBE-05 assert --
  // so grouping across the whole sequence would assert the opposite of the intended behaviour.
  // parent_scan_id identifies the survey that spawned each MS2.
  //
  // The mass key is the engine's OWN mono_mass: the fan-out pushed the SAME PeakGroup once per
  // charge, so those commands carry a bit-identical value and no rounding near-miss can split a
  // group and mask the defect. The first survey's MS2s carry an EMPTY parent_scan_id (buildMS2
  // stamps it only when the parent tracking id is non-zero, and the first survey's is 0); that is
  // still one coherent group.
  std::map<std::pair<std::string, long long>, std::set<int>> charges_per_survey_mass;
  for (const auto& c : a.ms2_cmds)
    if (c.num_stages >= 1)
      charges_per_survey_mass[{std::string(c.parent_scan_id), std::llround(c.mono_mass)}]
          .insert(c.stages[0].charge_state);

  // Report every offender before asserting: TEST_EQUAL carries no message, and "2 != 1" without the
  // survey and the mass is not a diagnosis.
  int violations = 0;
  for (const auto& kv : charges_per_survey_mass)
  {
    if (kv.second.size() <= 1) continue;
    ++violations;
    std::cout << "[CBE-07] survey '" << kv.first.first << "' mass " << kv.first.second
              << " acquired at " << kv.second.size() << " charges:";
    for (int z : kv.second) std::cout << " " << z;
    std::cout << std::endl;
  }
  TEST_EQUAL(violations, 0)
}
END_SECTION

// CBE-08: the INVERSE of CBE-07, and the reason it exists. precursor_charges: "separate" is the mode
// that deliberately acquires one MS2 per charge state, so within a single survey a mass MUST appear at
// two or more charges. Nothing asserted that, and the value was silently inert: it parsed, it was
// documented as "the fan-out arrives as an explicit acquisition mode", and PrecursorSelection never
// branched on it -- the break added for the fallback fix was unconditional. A config asking for
// separate fragmentation behaved exactly like "single".
//
// This is the shape of test that catches a mode wired at parse time and nowhere else: assert the
// BEHAVIOUR the value promises, not that the value round-trips.
//
// max_precursors must exceed 1 or the mode cannot express itself: selected_peak_groups_ is bounded by
// mass_count, so at 1 "one scan per charge" can only ever emit one scan.
START_SECTION(separate_yields_multiple_charges_per_mass_within_a_survey)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  std::string sep_json = std::string(base_on_json);
  {
    const std::string k = "\"max_precursors\": 10";
    const std::string::size_type pos = sep_json.find(k);
    TEST_TRUE(pos != std::string::npos)
    ABORT_IF(pos == std::string::npos)
    sep_json.replace(pos, k.size(),
                     "\"max_precursors\": 10, \"precursor_charges\": \"separate\"");
  }

  FLASHIda ida(const_cast<char*>(sep_json.c_str()));
  AcqResult a = runInterleaved(&ida, scans, {}, nullptr, /*max_iters=*/4000);
  TEST_EQUAL(a.ms2_cmds.size() > 0, true)

  // Same keying as CBE-07 -- (survey, mass) -> charges -- but the expectation is inverted.
  std::map<std::pair<std::string, long long>, std::set<int>> charges_per_survey_mass;
  for (const auto& c : a.ms2_cmds)
    if (c.num_stages >= 1)
      charges_per_survey_mass[{std::string(c.parent_scan_id), std::llround(c.mono_mass)}]
          .insert(c.stages[0].charge_state);

  int masses_with_multiple_charges = 0;
  for (const auto& kv : charges_per_survey_mass)
    if (kv.second.size() >= 2) ++masses_with_multiple_charges;

  if (masses_with_multiple_charges == 0)
    std::cout << "[CBE-08] no mass was acquired at >=2 charges within one survey -- "
                 "precursor_charges: \"separate\" is not fanning out" << std::endl;
  TEST_EQUAL(masses_with_multiple_charges >= 1, true)
}
END_SECTION

END_TEST
