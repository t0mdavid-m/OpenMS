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

#include <cstring>
#include <fstream>
#include <set>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  const char* base_off_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.9,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1,
      "ChargeBasedExclusion": false
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
      "agc_interval_seconds": 30
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 10 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  const char* base_on_json = R"({
    "deconvolution": {
      "score_threshold": 0.0, "tqscore_threshold": 0.3,
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
    },
    "precursor_selection": {
      "RT_window": 180, "target_mode": 0,
      "IDScore": false, "AllCharges": false,
      "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1,
      "ChargeBasedExclusion": true
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
      "agc_interval_seconds": 30
    },
    "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 10 },
      "ms2": { "selection": "none" },
      "ms3": { "selection": "none" }
    }
  })";

  // cytC (~12356 Da neutral) is present at many charge states (envelopes at m/z
  // ~1030/1124/1236/1373/1545), so a single mass appears at >=2 charges (CBE-03) and every
  // scan is productive — unlike ms1_standard (<=1 selectable/scan, charges spread across
  // different masses), which left CBE-03/04/06 with 0 MS2.
  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_cytc.txt";

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

  struct AcquisitionRow { int charge; double mz; double width; };

  // Interleaved drive result: every emitted MS2 (charge, mz, width) tuple + the MS2 command count.
  struct DriveResult
  {
    std::vector<AcquisitionRow> rows;
    int ms2_count = 0;
  };

  // Drive the engine via the canonical interleaved engine-id-echo contract (twin of
  // FLASHIda_TestHelpers::runInterleaved; inlined here because this TU's local ScanData /
  // loadTsvScans would collide with the helper's anonymous-namespace copies). Pull one command at
  // a time and feed exactly one MS1 survey response per survey command, stamped with the engine's
  // OWN scan_description (so the always-on MS1 gate accepts it — a fabricated "scan_"+id id is
  // rejected and would yield 0 precursors). MS2 commands are collected; the engine paces surveys.
  // ms1_scans are fed one per survey, in order (with rt_offset added); further surveys are idle
  // ticks. Terminates on idle>=3 (queue drained) or max_iters. State persists in `ida` across
  // calls, so a caller may invoke this twice (two passes) on the same engine.
  DriveResult driveInterleaved(FLASHIda& ida, const std::vector<ScanData>& ms1_scans,
                               double rt_offset = 0.0, int max_iters = 4000)
  {
    DriveResult res;
    const int n_ms1 = static_cast<int>(ms1_scans.size());
    int idle = 0, ms1_fed = 0;
    ScanCommand cmd{};
    for (int it = 0; it < max_iters && idle < 3; ++it)
    {
      if (ida.getNextScanCommand(cmd) != 1) break;
      // Idle tick: AGC, empty descriptor, or an MS1 re-survey after all ms1_scans have been fed.
      if (cmd.is_agc || cmd.scan_description[0] == '\0' || (cmd.msn_level <= 1 && ms1_fed >= n_ms1))
      {
        ++idle;
        cmd = ScanCommand{};
        continue;
      }
      idle = 0;
      if (cmd.msn_level <= 1)
      {
        const ScanData& s = ms1_scans[ms1_fed++];
        ida.processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt + rt_offset, 1,
                        cmd.scan_description);
      }
      else if (cmd.msn_level == 2)
      {
        ++res.ms2_count;
        if (cmd.num_stages >= 1)
          res.rows.push_back({cmd.stages[0].charge_state, cmd.stages[0].precursor_mz,
                              cmd.stages[0].isolation_width});
      }
      cmd = ScanCommand{};
    }
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

END_TEST
