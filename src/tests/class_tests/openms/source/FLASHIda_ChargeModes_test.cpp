// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Behavioural tests for precursor_selection.precursor_charges -- the MS2 acquisition GEOMETRY axis
// (ADR-0016): single | separate | multiplexed.
//
// These assert what each value PROMISES, under a config that sets nothing else unusual. That
// distinction is the whole reason the file exists. "separate" was previously covered by exactly one
// section, and that section configured charge_based_exclusion: true -- an unrelated exclusion-KEYING
// developer flag which, as a side effect, was the only thing that made the mode's charge list
// multi-valued. With the flag at its default the mode emitted one scan, i.e. behaved exactly like
// "single", in every shipped config, and the suite was green throughout.
//
// That flag is now gone (ADR-0021), so geometry has exactly one source and these sections exercise
// it directly. The lesson generalises past this key: a mode wired at parse time and asserted only
// for round-trip is indistinguishable from a mode that does nothing.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PlanUnit.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include "FLASHIda_TestAccess.h"   // PrecursorSelection's acquisition memory, read without a public getter
#include "FLASHIda_TestHelpers.h"  // ground-truth harness: ScanData/loadTsvScans/runInterleaved

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;

namespace
{
  // cytC (~12356 Da neutral) resolves at 9-10 charge states per species in this fixture, so the
  // difference between "one charge" and "the envelope" is large and unambiguous. ms1_standard would
  // not do: it yields <=1 selectable mass per scan.
  const std::string ms1_tsv_path = "../../FlashIDA/test-data/spectra/ms1_cytc.txt";

  // One config, one varying value. max_precursors is 10 so the budget never masks the geometry:
  // the species guard bounds selection by SPECIES, and a whole envelope costs ONE slot (ADR-0016),
  // but a cramped budget would still make "did it fan out?" ambiguous.
  // The AUTHORED CHARGE SET fixtures (ADR-0028). cytC resolves at EVERY charge z=8..20 in this
  // spectrum, so naming three of them is a large, visible restriction rather than a no-op, and the
  // envelope's most intense charge (z=15) is deliberately NOT among them -- an anchor that ignored
  // the set would land outside it.
  const std::string charge_set_list      = "../../FlashIDA/test-data/configs/inclusion_cytc_charges.txt";
  const std::string charge_set_rows_list = "../../FlashIDA/test-data/configs/inclusion_cytc_charges_rows.txt";
  const std::string unrestricted_list    = "../../FlashIDA/test-data/configs/inclusion_cytc.txt";
  const long long   authored_mass        = 12351;   // llround of the authored monoisotopic mass
  const std::set<int> authored_charges{10, 13, 16};

  // The WIDE authored row (ADR-0036), naming a contiguous span the winning PeakGroup often cannot
  // cover on its own. `10;13;16` above cannot exercise split-envelope completion at all: its winner
  // carries all three in 18 of the 25 productive surveys, and in the 7 where it does not, the
  // missing charges were already spent by an earlier survey -- so no sibling ever has anything to
  // complete. Against `12..18` the same spectrum behaves differently: the winner carries all seven
  // in 12 surveys, carries a strict non-empty subset in several more (survey 6 is the first), and
  // the species' PeakGroups collectively carry all seven in 23 of 25.
  const std::string wide_set_list      = "../../FlashIDA/test-data/configs/inclusion_cytc_charges_wide.txt";
  const std::string wide_set_rows_list = "../../FlashIDA/test-data/configs/inclusion_cytc_charges_wide_rows.txt";
  const std::set<int> wide_charges{12, 13, 14, 15, 16, 17, 18};

  bool isSubsetOfAuthored(const std::set<int>& charges)
  {
    for (int z : charges)
      if (authored_charges.count(z) == 0) return false;
    return true;
  }

  // @p max_precursors and @p strict_inclusion default to what every CM-01..CM-08 case already used,
  // so adding them changes no existing config text. CM-09 needs both: the split-envelope completion
  // it exercises is reachable only while the candidate loop is still running (ADR-0036 leaves the
  // cap hard), and under strict inclusion a non-target never spends a slot, which is what makes a
  // budget of 2 enough for one target and its siblings.
  std::string cfgFor(const std::string& precursor_charges,
                     const std::string& targeting = "none",
                     const std::string& inclusion_list = "",
                     int max_precursors = 10,
                     bool strict_inclusion = false,
                     // ADR-0040. Emitted unconditionally at their INERT defaults, so every existing
                     // caller's config is unchanged in behaviour while the keys stay exercised.
                     int max_charge_states = 0,
                     int min_charge_states = 1,
                     // ADR-0059. Emitted unconditionally at its default ("snr").
                     const std::string& charge_selection = "snr")
  {
    return std::string(R"JSON({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.3,
    "min_charge": 4,
    "max_charge": 50,
    "min_mass": 500,
    "max_mass": 50000,
    "tol": [10, 10, 10]
  },
  "flashtnt": { "min_length": 3, "max_length": 8 },
  "faims": { "cv_values": [], "max_cv_skip": 0 },
  "scheduling": {
    "cycle_time": { "enabled": false, "value_ms": 60000 },
    "scan_timeout": { "enabled": true, "value_ms": 30000 },
    "agc_interval_seconds": 30
  },
  "files": { "target_logs": [], "fasta": "", "inclusion_list": ")JSON")
           + inclusion_list + R"JSON(", "ptm_list": "" },
  "precursor_selection": {
    "rt_window": 180,
    "targeting": ")JSON" + targeting + R"JSON(",
    "consider_all_charges": false,
    "strict_inclusion": )JSON" + (strict_inclusion ? "true" : "false") + R"JSON(,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": )JSON" + std::to_string(max_precursors) + R"JSON(,
    "precursor_charges": ")JSON" + precursor_charges + R"JSON(",
    "max_charge_states": )JSON" + std::to_string(max_charge_states) + R"JSON(,
    "min_charge_states": )JSON" + std::to_string(min_charge_states) + R"JSON(,
    "charge_selection": ")JSON" + charge_selection + R"JSON("
  },
  "characterization": { "mode": "off", "max_targets": 10 },
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
  "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 }
}
)JSON";
  }

  // Every charge the emitted MS2s isolate, keyed by (survey, mass).
  //
  // Keying per SURVEY, not per run, is the point: charge-keyed fallback legitimately reaches other
  // charges of a mass on LATER surveys, so a run-wide grouping would conflate "acquired the envelope
  // at once" with "drifted onto another charge eventually" -- the exact conflation that let the
  // defect survive. parent_scan_id identifies the spawning survey.
  //
  // @p with_notches folds a multiplexed command's co-isolation notches in alongside its anchor, so
  // the two geometries become comparable as SETS.
  std::map<std::pair<std::string, long long>, std::set<int>>
  chargesPerSurveyMass(const AcqResult& a, bool with_notches)
  {
    std::map<std::pair<std::string, long long>, std::set<int>> out;
    for (const auto& c : a.ms2_cmds)
    {
      if (c.num_stages < 1) continue;
      auto& charges = out[{std::string(c.parent_scan_id), std::llround(c.mono_mass)}];
      charges.insert(c.stages[0].charge_state);
      if (with_notches)
      {
        auto [first, count] = notchesForStage(c, 0);
        for (int i = 0; i < count; ++i) charges.insert(first[i].charge_state);
      }
    }
    return out;
  }

  // The FIRST survey's acquisitions only, as mass -> charges.
  //
  // "First survey" is read from the run itself -- the parent of the first MS2 command emitted -- and
  // NOT assumed to be the empty string. buildMS2 stamps parent_scan_id only when the parent tracking
  // id is non-zero (ScanCommandQueue.cpp:345), and by the time the first productive survey runs, the
  // idle AGC/MS1 ticks have already consumed ids, so every MS2 here carries a real parent. Assuming
  // "" matched nothing at all and compared zero groups.
  //
  // Reading it per-run is not incidental: the two runs mint different ids once one emits several
  // times as many commands, so the id is only meaningful within its own run. The first survey itself
  // is comparable across runs because both start from the same engine state and the same MS1.
  std::map<long long, std::set<int>> firstSurveyCharges(const AcqResult& a, bool with_notches)
  {
    std::map<long long, std::set<int>> out;
    if (a.ms2_cmds.empty()) return out;

    const std::string first_survey(a.ms2_cmds.front().parent_scan_id);
    for (const auto& c : a.ms2_cmds)
    {
      if (c.num_stages < 1) continue;
      if (std::string(c.parent_scan_id) != first_survey) continue;
      auto& charges = out[std::llround(c.mono_mass)];
      charges.insert(c.stages[0].charge_state);
      if (with_notches)
      {
        auto [first, count] = notchesForStage(c, 0);
        for (int i = 0; i < count; ++i) charges.insert(first[i].charge_state);
      }
    }
    return out;
  }

  AcqResult runMode(FLASHIda& ida, const std::vector<ScanData>& scans)
  {
    return runInterleaved(&ida, scans, {}, nullptr, /*max_iters=*/4000);
  }

  // ---- ADR-0059 oracle -------------------------------------------------------------------------
  //
  // The charge_selection sections below check the engine against an INDEPENDENT reading of the
  // same PeakGroup: raw accessors and a purity summed here over the survey spectrum, never the
  // engine's peakGroupNotchCandidates / windowPurity / rankCharges. A synthetic PeakGroup cannot
  // carry per-charge intensity or qscores (only PeakGroup's private deconvolution path fills them),
  // so these run on real ms1_cytc surveys.
  struct OracleCharge { int z; double snr, intensity, qscore, purity; };

  std::vector<OracleCharge> oracleCharges(const PeakGroup& pg, const MSSpectrum& raw, int min_charge)
  {
    std::vector<OracleCharge> out;
    const auto qs = pg.getAllQscores();
    auto [min_c, max_c] = pg.getAbsChargeRange();
    for (int z = min_c; z <= max_c; ++z)
    {
      auto [lo, hi] = pg.getMzRange(z);
      if (hi <= lo) continue;
      if (min_charge > 0 && z < min_charge) continue;
      const double s = pg.getChargeIntensity(z);
      double total = 0.0;
      for (const auto& p : raw)
        if (p.getMZ() >= lo - 0.4 && p.getMZ() <= hi + 0.4) total += p.getIntensity();
      const double purity = (s <= 0.0 || total <= 0.0) ? -1.0 : std::min(1.0, s / total);
      const auto q = qs.find(z);
      out.push_back({z, pg.getChargeSNR(z), s, q != qs.end() ? (double)q->second : pg.getQscore(), purity});
    }
    return out;
  }

  double oracleKey(const OracleCharge& c, const std::string& rule)
  {
    if (rule == "snr") return c.snr;
    if (rule == "qscore") return c.qscore;
    if (rule == "intensity") return c.intensity;
    return c.purity;
  }

  // The rule's order: key descending, then intensity descending, then the lower charge.
  std::vector<int> oracleOrder(std::vector<OracleCharge> cands, const std::string& rule)
  {
    std::sort(cands.begin(), cands.end(), [&rule](const OracleCharge& a, const OracleCharge& b) {
      const double ka = oracleKey(a, rule), kb = oracleKey(b, rule);
      if (ka != kb) return ka > kb;
      if (a.intensity != b.intensity) return a.intensity > b.intensity;
      return a.z < b.z;
    });
    std::vector<int> out;
    for (const auto& c : cands) out.push_back(c.z);
    return out;
  }

  // Drive PrecursorSelection survey by survey -- the decision whose charges buildMS2 is handed.
  template <typename F>
  void forEachSurvey(const std::string& cfg_json, const std::vector<ScanData>& scans, F fn)
  {
    Config cfg{cfg_json};
    Deconvolution deconv(cfg, cfg.toleranceList());
    PrecursorSelection ps(cfg, deconv);
    for (size_t si = 0; si < scans.size(); ++si)
    {
      const auto& s = scans[si];
      const int n = ps.filterAndRank(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 1, 0.0);
      fn(cfg, ps, n, (int)si);
    }
  }

  const std::vector<std::string> kRules{"snr", "qscore", "intensity", "purity"};

  // ---- Survey re-planning: the selection journal -----------------------------------------------
  //
  // filterAndRank still writes acquisition memory in place; it ALSO journals every such write as a
  // MemoryStamp, so a later push can roll the in-place writes back and apply a species' stamps only
  // when it commits. These helpers serve the sections that pin the journal against the writes.

  const std::string ms1_standard_path = "../../FlashIDA/test-data/spectra/ms1_standard.txt";
  const std::string ms1_ecoli_path    = "../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt";
  const std::string target_log_path   = "../../FlashIDA/test-data/configs/test_target_log.log";

  /// Every stamp one survey journaled, merged back into WRITE order -- the order the in-place writes
  /// were made in, which replay must reproduce because mass_qscore_map_ is an overwrite.
  std::vector<MemoryStamp> allStampsBySeq(const SelectionJournal& j)
  {
    std::vector<MemoryStamp> all = j.survey_level;
    for (const auto& kv : j.species)
    {
      const std::vector<MemoryStamp> f = PrecursorSelection::flatten(kv.second);
      all.insert(all.end(), f.begin(), f.end());
    }
    std::sort(all.begin(), all.end(), [](const MemoryStamp& a, const MemoryStamp& b) { return a.seq < b.seq; });
    return all;
  }

  /// cfgFor with one target log -- the input of the two modes that read one (in_depth, legacy inclusion).
  std::string withTargetLog(std::string cfg, const std::string& log_path)
  {
    const std::string none = "\"target_logs\": []";
    cfg.replace(cfg.find(none), none.size(), "\"target_logs\": [\"" + log_path + "\"]");
    return cfg;
  }

  // The ida.log target behaviour before FDmerge 1fe5ee0b3e, under its two switches: a target is matched only at a
  // charge the reference acquired it at and each acquisition uses that charge up (the erase the journal stamps as
  // LegacyTargetChargeErase), at every retention time (these logs' RTs do not match the fixtures' surveys).
  std::string withReferenceLogTargets(std::string cfg)
  {
    const std::string rt = "\"rt_window\": 180,";
    cfg.replace(cfg.find(rt), rt.size(), rt + " \"log_target_rt_window\": false, \"log_target_charges\": \"reference\",");
    return cfg;
  }

  // ---- A survey with more species than slots (CM-12, CM-13) -------------------------------------
  //
  // No committed fixture has one: ms1_standard yields <= 1 selectable mass per scan and ms1_cytc is
  // one species. So these build the deconvolved survey directly and run the engine's own selection
  // on it (FLASHIdaTestAccess::selectFromDeconvolved), which needs no fixture at all.

  /// A synthetic MS1 PeakGroup resolving every charge in [z_lo, z_hi]. Three peaks per charge, as in
  /// FLASHIda_CandidateAdmission_test's makeEnvelope: with one peak getMzRange returns {mz, mz},
  /// which every consumer reads as "not resolved", and the fan-out would silently shrink to the anchor.
  PeakGroup syntheticSpecies(double mono_mass, int z_lo, int z_hi, double qscore)
  {
    const double proton = 1.00727646677;
    const double iso_gap = 1.00235;
    PeakGroup pg(z_lo, z_hi, true);
    pg.setMonoisotopicMass(mono_mass);
    for (int z = z_lo; z <= z_hi; ++z)
    {
      for (int iso = 0; iso < 3; ++iso)
      {
        FLASHHelperClasses::LogMzPeak lp;
        lp.mz = (mono_mass + z * proton + iso * iso_gap) / z;
        lp.abs_charge = z;
        pg.push_back(lp);
      }
      pg.setChargeSNR(z, 10.0f);
      pg.setChargeIsotopeCosine(z, 0.99f);
    }
    pg.setQscore(qscore);
    pg.setRepAbsCharge(z_lo);
    return pg;
  }

  /// Selection time (seconds) for splitEnvelope (FLASHIda_TestHelpers.h), inside the wide row's window.
  const double kSplitRt = 60.0;

  /// selectFrom on a fresh selection built from @p cfg, with stdout captured into @p log.
  int selectSplit(const std::string& cfg, std::unique_ptr<Config>& config, std::unique_ptr<Deconvolution>& deconv,
                  std::unique_ptr<PrecursorSelection>& ps, std::string& log)
  {
    config = std::make_unique<Config>(cfg);
    deconv = std::make_unique<Deconvolution>(*config, config->toleranceList());
    ps = std::make_unique<PrecursorSelection>(*config, *deconv);
    std::ostringstream captured;
    std::streambuf* const real_cout = std::cout.rdbuf(captured.rdbuf());
    int n = 0;
    try { n = FLASHIdaTestAccess::selectFrom(*ps, splitEnvelope(), kSplitRt); }
    catch (...) { std::cout.rdbuf(real_cout); throw; }
    std::cout.rdbuf(real_cout);
    log = captured.str();
    return n;
  }

  /// Six well-separated species, each resolved at charges 8..11, with strictly decreasing qscores so
  /// the ranking is total: rank k is mass 10000 + 1500 k.
  const int syn_species = 6;
  double synMass(int k) { return 10000.0 + 1500.0 * k; }

  DeconvolvedSpectrum syntheticSurvey()
  {
    DeconvolvedSpectrum ds(1);
    for (int k = 0; k < syn_species; ++k)
      ds.push_back(syntheticSpecies(synMass(k), 8, 11, 0.95 - 0.01 * k));
    return ds;
  }

  /// The species (rounded masses) one survey acquired.
  std::set<long long> speciesOf(const std::vector<std::pair<double, int>>& acquisitions)
  {
    std::set<long long> out;
    for (const auto& a : acquisitions) out.insert(std::llround(a.first));
    return out;
  }

  /// One survey of a fresh engine on the synthetic spectrum. @p targets as in selectFromDeconvolved.
  std::vector<std::pair<double, int>> selectSynthetic(const std::string& cfg, std::vector<double> targets = {})
  {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    return FLASHIdaTestAccess::selectFromDeconvolved(ida, syntheticSurvey(), /*rt=*/60.0, std::move(targets));
  }
}

START_TEST(FLASHIda_ChargeModes, "$Id$")

// CM-01: the assertion that was missing. "separate" means one MS2 per SNR-positive charge state, so
// within a SINGLE survey a mass must appear at two or more charges -- under a config that does
// nothing but ask for it.
//
// Under the defect every group has exactly one charge and this fails, which is the point.
START_SECTION(separate_fans_out_under_a_default_config)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("separate");
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)  // non-vacuity: with no MS2 the assertion is trivial

  const auto groups = chargesPerSurveyMass(a, /*with_notches=*/false);
  int masses_with_multiple_charges = 0;
  for (const auto& kv : groups)
    if (kv.second.size() >= 2) ++masses_with_multiple_charges;

  if (masses_with_multiple_charges == 0)
    std::cout << "[CM-01] no mass acquired at >=2 charges within one survey -- "
                 "precursor_charges: \"separate\" is not fanning out" << std::endl;
  TEST_EQUAL(masses_with_multiple_charges >= 1, true)
}
END_SECTION

// CM-02: ADR-0016's "the two modes differ only in SCAN COUNT, never in which charges" as an
// executable invariant. Both must resolve the same acquisition charge set from the same
// peakGroupNotchCandidates + selectNotches pair; sourcing them differently is precisely how
// "separate" ended up reading an exclusion flag while "multiplexed" read the PeakGroup.
//
// Scoped to the FIRST survey, which both runs reach from the same engine state and the same MS1.
// Later surveys are not comparable across runs: "separate" emits several times as many commands, so
// it consumes more tracking ids and the two runs' survey ids -- and then their exclusion state --
// diverge. Comparing those would test id allocation, not geometry.
START_SECTION(separate_and_multiplexed_acquire_the_same_charge_set)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string sep_cfg = cfgFor("separate");
  const std::string mux_cfg = cfgFor("multiplexed");

  FLASHIda sep_ida(const_cast<char*>(sep_cfg.c_str()));
  FLASHIda mux_ida(const_cast<char*>(mux_cfg.c_str()));
  AcqResult sep = runMode(sep_ida, scans);
  AcqResult mux = runMode(mux_ida, scans);

  TEST_EQUAL(sep.ms2_cmds.size() > 0, true)
  TEST_EQUAL(mux.ms2_cmds.size() > 0, true)

  // separate: N commands, anchor only. multiplexed: 1 command, anchor + notches.
  const auto sep_groups = firstSurveyCharges(sep, /*with_notches=*/false);
  const auto mux_groups = firstSurveyCharges(mux, /*with_notches=*/true);

  int compared = 0, mismatches = 0;
  for (const auto& kv : mux_groups)
  {
    auto it = sep_groups.find(kv.first);
    if (it == sep_groups.end()) continue;
    ++compared;
    if (it->second == kv.second) continue;
    ++mismatches;
    std::cout << "[CM-02] mass " << kv.first << " separate={";
    for (int z : it->second) std::cout << " " << z;
    std::cout << " } multiplexed={";
    for (int z : kv.second) std::cout << " " << z;
    std::cout << " }" << std::endl;
  }

  TEST_EQUAL(compared > 0, true)  // non-vacuity: an empty intersection would pass trivially
  TEST_EQUAL(mismatches, 0)

  // The modes must differ in scan count, or "separate" collapsed into "multiplexed" instead of
  // fanning out and the set comparison above would agree for the wrong reason.
  TEST_EQUAL(sep.ms2_cmds.size() > mux.ms2_cmds.size(), true)
}
END_SECTION

// CM-03: the default mode's within-survey invariant -- a mass is acquired at exactly ONE charge.
// Guards the other direction: a fan-out that fires regardless of mode would spend the whole budget
// on one proteoform's charge states and never fragment P2 or P3.
START_SECTION(single_acquires_one_charge_per_mass_per_survey)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("single");
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)

  const auto groups = chargesPerSurveyMass(a, /*with_notches=*/false);
  int violations = 0;
  for (const auto& kv : groups)
  {
    if (kv.second.size() <= 1) continue;
    ++violations;
    std::cout << "[CM-03] survey '" << kv.first.first << "' mass " << kv.first.second
              << " acquired at " << kv.second.size() << " charges:";
    for (int z : kv.second) std::cout << " " << z;
    std::cout << std::endl;
  }
  TEST_EQUAL(violations, 0)

  // "single" must also command no notches -- one charge means one isolation window.
  int notched = 0;
  for (const auto& c : a.ms2_cmds)
    if (notchesForStage(c, 0).second > 0) ++notched;
  TEST_EQUAL(notched, 0)
}
END_SECTION

// CM-04: an AUTHORED CHARGE SET is a RESTRICTION (ADR-0028 decision 1). Under multiplexed the scan
// co-isolates the envelope, so if the set were ignored the command would carry ~10 charges, most of
// them unnamed. It must carry only named ones -- and still more than one, or "restrict" quietly
// became "acquire a single charge" and the co-isolation was thrown away with the rest.
START_SECTION(authored_charge_set_restricts_the_acquisition)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("multiplexed", "inclusion", charge_set_list);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)  // non-vacuity

  const auto groups = chargesPerSurveyMass(a, /*with_notches=*/true);
  int checked = 0, violations = 0, co_isolated = 0;
  for (const auto& kv : groups)
  {
    if (kv.first.second != authored_mass) continue;
    ++checked;
    if (kv.second.size() >= 2) ++co_isolated;
    if (isSubsetOfAuthored(kv.second)) continue;
    ++violations;
    std::cout << "[CM-04] survey '" << kv.first.first << "' isolated unauthored charges:";
    for (int z : kv.second)
      if (authored_charges.count(z) == 0) std::cout << " " << z;
    std::cout << " (authored 10;13;16)" << std::endl;
  }

  TEST_EQUAL(checked > 0, true)      // the authored mass was selected at all
  TEST_EQUAL(violations, 0)          // FAILS if the notch filter is skipped
  TEST_EQUAL(co_isolated > 0, true)  // FAILS if the filter emptied the notch set
}
END_SECTION

// CM-05: the anchor is drawn FROM the set (ADR-0028 decision 2). cytC resolves at 13 charges here
// and three are named, so an anchor picked from the envelope rather than the set lands outside it
// almost surely -- which is exactly the regression this guards: the old code reassigned the charge
// inside the inclusion matcher, and that reassignment is now gone.
START_SECTION(authored_charge_set_moves_the_anchor)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("single", "inclusion", charge_set_list);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)

  const auto first = firstSurveyCharges(a, /*with_notches=*/false);
  auto it = first.find(authored_mass);
  TEST_EQUAL(it != first.end(), true)  // non-vacuity: the authored mass must have been acquired
  ABORT_IF(it == first.end())

  // "single" acquires exactly one charge per species per survey, and it must be a named one.
  TEST_EQUAL(it->second.size(), 1)
  if (!isSubsetOfAuthored(it->second))
    std::cout << "[CM-05] anchor " << *it->second.begin()
              << " is not one of the authored charges 10;13;16 -- the set did not reach the anchor"
              << std::endl;
  TEST_EQUAL(isSubsetOfAuthored(it->second), true)
}
END_SECTION

// CM-06: per-charge exclusion walks the set across surveys (ADR-0028 decision 3b) -- the load-
// bearing one. Under mass-keyed exclusion the species is retired by its own first acquisition and
// the second and third named charges are unreachable, so this fails with exactly one anchor.
START_SECTION(authored_charge_set_is_walked_across_surveys)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("single", "inclusion", charge_set_list);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)

  // Anchors for the authored mass, one entry per survey that acquired it.
  std::vector<int> anchors;
  std::set<std::string> surveys_seen;
  for (const auto& c : a.ms2_cmds)
  {
    if (c.num_stages < 1) continue;
    if (std::llround(c.mono_mass) != authored_mass) continue;
    const std::string survey(c.parent_scan_id);
    if (!surveys_seen.insert(survey).second) continue;
    anchors.push_back(c.stages[0].charge_state);
  }

  std::cout << "[CM-06] authored mass acquired at charges:";
  for (int z : anchors) std::cout << " " << z;
  std::cout << " (authored 10;13;16)" << std::endl;

  const std::set<int> distinct(anchors.begin(), anchors.end());
  TEST_EQUAL(anchors.size() >= 2, true)           // FAILS under mass-keyed exclusion
  TEST_EQUAL(distinct.size(), anchors.size())     // FAILS if a spent charge is re-acquired
  TEST_EQUAL(isSubsetOfAuthored(distinct), true)  // FAILS if the walk escapes the set
  TEST_EQUAL(anchors.size() <= 3, true)           // FAILS if the set never terminates
}
END_SECTION

// CM-07: several rows naming one mass union their sets (ADR-0028 decision 5), so the one-row and
// three-row spellings are the same method. It is also the regression guard for the matcher that
// took the first ACTIVE row regardless of its mass: under that bug the three-row file hands over
// the charge of whichever row std::sort happened to leave first, and the two runs diverge.
START_SECTION(rows_naming_one_mass_union_their_charge_sets)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string one_cfg  = cfgFor("multiplexed", "inclusion", charge_set_list);
  const std::string rows_cfg = cfgFor("multiplexed", "inclusion", charge_set_rows_list);

  FLASHIda one_ida(const_cast<char*>(one_cfg.c_str()));
  FLASHIda rows_ida(const_cast<char*>(rows_cfg.c_str()));
  AcqResult one = runMode(one_ida, scans);
  AcqResult rows = runMode(rows_ida, scans);

  TEST_EQUAL(one.ms2_cmds.size() > 0, true)
  TEST_EQUAL(rows.ms2_cmds.size() > 0, true)

  // Scoped to the first survey for the same reason CM-02 is: later surveys are only comparable
  // while both runs have minted the same tracking ids.
  const auto one_first  = firstSurveyCharges(one, /*with_notches=*/true);
  const auto rows_first = firstSurveyCharges(rows, /*with_notches=*/true);

  auto a_it = one_first.find(authored_mass);
  auto b_it = rows_first.find(authored_mass);
  TEST_EQUAL(a_it != one_first.end(), true)  // non-vacuity on both sides
  TEST_EQUAL(b_it != rows_first.end(), true)
  ABORT_IF(a_it == one_first.end() || b_it == rows_first.end())

  if (a_it->second != b_it->second)
  {
    std::cout << "[CM-07] one-row={";
    for (int z : a_it->second) std::cout << " " << z;
    std::cout << " } three-row={";
    for (int z : b_it->second) std::cout << " " << z;
    std::cout << " }" << std::endl;
  }
  TEST_EQUAL(a_it->second == b_it->second, true)
}
END_SECTION

// CM-08: a row naming NO charge stays unrestricted. This is the guard behind the claim that every
// committed config -- all of which write -1 -- is byte-identical. If authoredChargesFor_ ever
// returned a non-empty set for such a row, the whole envelope would collapse onto three charges
// here and every existing inclusion golden would move with it.
START_SECTION(a_row_naming_no_charge_is_unrestricted)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string cfg = cfgFor("multiplexed", "inclusion", unrestricted_list);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runMode(ida, scans);

  TEST_EQUAL(a.ms2_cmds.size() > 0, true)

  const auto first = firstSurveyCharges(a, /*with_notches=*/true);
  auto it = first.find(authored_mass);
  TEST_EQUAL(it != first.end(), true)
  ABORT_IF(it == first.end())

  // cytC co-isolates ~10 of its 13 resolved charges here, so an unrestricted row must reach well
  // outside {10,13,16}.
  TEST_EQUAL(isSubsetOfAuthored(it->second), false)
}
END_SECTION

// CM-09: a SPLIT ENVELOPE is completed within one survey (ADR-0036), and the budget cap is what
// bounds it.
//
// The discriminator has to be "two charges of one mass in ONE survey", not "the named set completed
// eventually": within a single survey under `single`, a mass reaching two charges is only possible if a
// second PeakGroup of it was admitted. The split is synthetic (splitEnvelope): FDmerge's deconvolution
// no longer leaves one in ms1_cytc.
//
// The max_precursors 1 arm pins the other half of the decision: the cap stays HARD, so a sibling is
// not examined once the loop has ended, and the run says so rather than silently acquiring less.
START_SECTION(a_split_authored_envelope_is_completed_within_one_survey)
{
  std::map<int, std::set<int>> taken;   // max_precursors -> charges of the authored mass in the one survey
  std::map<int, std::string> logs;
  for (int cap : {2, 1})
  {
    std::unique_ptr<Config> config;
    std::unique_ptr<Deconvolution> deconv;
    std::unique_ptr<PrecursorSelection> ps;
    const int n = selectSplit(cfgFor("single", "inclusion", wide_set_list, cap, /*strict_inclusion=*/true), config, deconv, ps,
                              logs[cap]);
    for (int i = 0; i < n; ++i)
      if (std::llround(ps->selectedPeakGroups()[i].getMonoMass()) == authored_mass) taken[cap].insert(ps->triggerCharges()[i]);
  }

  // Both PeakGroups contributed a charge: the sibling was admitted and completed part of the named set.
  TEST_EQUAL(taken[2].size(), 2)
  TEST_EQUAL(taken[2].count(12) + taken[2].count(13) + taken[2].count(14), 1)
  TEST_EQUAL(taken[2].count(15) + taken[2].count(16), 1)
  // Everything acquired is a NAMED charge: completion must not become "acquire more".
  for (int z : taken[2]) TEST_EQUAL(wide_charges.count(z) > 0, true)

  // The hard cap: at max_precursors 1 the loop ends on the first PeakGroup and no sibling is reached.
  TEST_EQUAL(taken[1].size(), 1)

  // ...and it is not silent about it. This is the whole reason the cap was left hard rather than
  // special-cased: the configuration remains legal and unchanged, and the run reports what it cost.
  TEST_EQUAL(logs[1].find("reason=budget_exhausted") != std::string::npos, true)
  TEST_EQUAL(logs[2].find("reason=budget_exhausted") == std::string::npos, true)
}
END_SECTION

// CM-10: the union rule (ADR-0028) survives a SPLIT. CM-07 already pins that one row naming
// `10;13;16` and three rows naming one charge each acquire identically -- but its winner resolves
// the whole named set every survey, so no sibling is ever admitted and the equivalence never
// crosses the completion path. Against the wide row it does.
START_SECTION(rows_naming_one_mass_union_their_charge_sets_under_a_split)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  const std::string one_cfg  = cfgFor("single", "inclusion", wide_set_list, 2, true);
  const std::string rows_cfg = cfgFor("single", "inclusion", wide_set_rows_list, 2, true);

  FLASHIda one_ida(const_cast<char*>(one_cfg.c_str()));
  FLASHIda rows_ida(const_cast<char*>(rows_cfg.c_str()));
  AcqResult one  = runMode(one_ida, scans);
  AcqResult rows = runMode(rows_ida, scans);

  TEST_EQUAL(one.ms2_cmds.size() > 0, true)
  TEST_EQUAL(one.ms2_cmds.size(), rows.ms2_cmds.size())
  ABORT_IF(one.ms2_cmds.size() != rows.ms2_cmds.size())

  const auto one_groups  = chargesPerSurveyMass(one, /*with_notches=*/true);
  const auto rows_groups = chargesPerSurveyMass(rows, /*with_notches=*/true);
  TEST_EQUAL(one_groups.size(), rows_groups.size())
  TEST_EQUAL(one_groups == rows_groups, true)
}
END_SECTION

// CM-11: an UNRESTRICTED row is untouched by any of this (ADR-0036 as amended; ADR-0037 withdrawn).
//
// This is the regression guard for the scope decision, and it is the test whose absence let a
// change aimed at charge completeness move thirteen log goldens -- every one of them a `-1` row,
// while the two authored-charge modes stayed still. Two things must hold for a `-1` row: no sibling
// is ever admitted (so a mass reaches at most one anchor per survey under `single`), and matching
// an inclusion row does not lift dynamic exclusion (so the mass is acquired once per rt window
// rather than on every survey that beats the stored qscore by any margin).
START_SECTION(an_unrestricted_row_is_unaffected_by_split_envelope_completion)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  // Under `single`: one anchor per mass per survey, exactly as CM-03 requires of a non-targeted run.
  {
    const std::string cfg = cfgFor("single", "inclusion", unrestricted_list);
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    AcqResult a = runMode(ida, scans);
    TEST_EQUAL(a.ms2_cmds.size() > 0, true)

    int violations = 0, surveys_with_the_mass = 0;
    for (const auto& kv : chargesPerSurveyMass(a, /*with_notches=*/false))
    {
      if (kv.first.second != authored_mass) continue;
      ++surveys_with_the_mass;
      if (kv.second.size() <= 1) continue;
      ++violations;
      std::cout << "[CM-11] survey '" << kv.first.first << "' acquired the unrestricted mass at "
                << kv.second.size() << " charges -- a sibling was admitted for a `-1` row"
                << std::endl;
    }
    TEST_EQUAL(surveys_with_the_mass > 0, true)   // non-vacuity
    TEST_EQUAL(violations, 0)

    // Dynamic exclusion still governs: the target is acquired ONCE inside the rt window, not on
    // every survey that resolves it a fraction better. rt_window is 180 s and this fixture is a
    // 63 s run, so nothing expires -- a count above one here means the bars were exempted.
    TEST_EQUAL(surveys_with_the_mass, 1)
  }

  // Under `multiplexed`: still one scan for the mass per survey, however many PeakGroups of it the
  // survey produced. The notch set is the winner's own envelope and is NOT extended by a sibling.
  {
    const std::string cfg = cfgFor("multiplexed", "inclusion", unrestricted_list);
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    AcqResult a = runMode(ida, scans);
    TEST_EQUAL(a.ms2_cmds.size() > 0, true)

    int scans_for_the_mass = 0;
    for (const auto& c : a.ms2_cmds)
      if (c.num_stages >= 1 && std::llround(c.mono_mass) == authored_mass) ++scans_for_the_mass;

    TEST_EQUAL(scans_for_the_mass, 1)
  }
}
END_SECTION

// CM-09: precursor_selection.max_charge_states caps the acquired set END TO END (ADR-0040).
//
// This is the test that pins the CALLER's arithmetic. max_charge_states counts the SET (anchor +
// notches) while AdmissionContext::max_notches counts notches ALONE, so PrecursorSelection must
// pass `max_charge_states - 1`. Passing it unconverted yields one charge too many -- a defect that
// is invisible in every log, because "max 3" producing 4 charges reads as entirely plausible.
// CA-B18 pins the +1 at the unit level; only this section pins the subtraction at the call site.
//
// Uses `separate`, where the cap bounds the SCAN COUNT per species and was previously bounded only
// by MAX_NOTCHES_PER_STAGE -- a notch-array capacity that `separate` never writes into.
START_SECTION(max_charge_states_caps_the_acquired_set)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  // Uncapped first, to establish that this fixture DOES resolve more charges than the cap allows.
  // Without this the capped assertion below could pass on a survey that never had 4 charges to cap.
  const std::string uncapped_cfg = cfgFor("separate");
  FLASHIda uncapped(const_cast<char*>(uncapped_cfg.c_str()));
  AcqResult ua = runMode(uncapped, scans);
  const auto ugroups = chargesPerSurveyMass(ua, /*with_notches=*/false);
  size_t widest = 0;
  for (const auto& kv : ugroups) widest = std::max(widest, kv.second.size());
  if (widest <= 3)
    std::cout << "[CM-09] uncapped run never acquired a mass at >3 charges (widest=" << widest
              << ") -- the cap assertion below would be vacuous" << std::endl;
  TEST_EQUAL(widest > 3, true)

  // Capped at 3. No mass may exceed it, and at least one must REACH it.
  const std::string capped_cfg = cfgFor("separate", "none", "", 10, false, /*max_charge_states=*/3);
  FLASHIda capped(const_cast<char*>(capped_cfg.c_str()));
  AcqResult ca = runMode(capped, scans);
  TEST_EQUAL(ca.ms2_cmds.size() > 0, true)   // non-vacuity

  const auto cgroups = chargesPerSurveyMass(ca, /*with_notches=*/false);
  size_t capped_widest = 0;
  for (const auto& kv : cgroups) capped_widest = std::max(capped_widest, kv.second.size());
  TEST_EQUAL(capped_widest <= 3, true)   // 4 here means the -1 was dropped at the call site
  TEST_EQUAL(capped_widest, 3)           // and the cap is REACHED, so it is not merely unreachable
}
END_SECTION

// CM-10 (ADR-0059): under multiplexed the WIRE isolates the set admission recorded. buildMS2 used to
// re-derive the notches with its own cap of nine, so max_charge_states never reached the instrument
// for an unrestricted row (the quant_multiplexed golden: cap 3, ten charges isolated).
START_SECTION(multiplexed_wire_respects_max_charge_states)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  // Uncapped, the fixture co-isolates more than three charges -- else the cap below is vacuous.
  const std::string uncapped_cfg = cfgFor("multiplexed");
  FLASHIda uncapped(const_cast<char*>(uncapped_cfg.c_str()));
  AcqResult ua = runMode(uncapped, scans);
  int uncapped_widest = 0;
  for (const auto& c : ua.ms2_cmds) uncapped_widest = std::max(uncapped_widest, 1 + notchesForStage(c, 0).second);
  TEST_EQUAL(uncapped_widest > 3, true)

  // End to end, capped at 3: no command isolates more, and at least one reaches it.
  const std::string capped_cfg = cfgFor("multiplexed", "none", "", 10, false, /*max_charge_states=*/3);
  FLASHIda capped(const_cast<char*>(capped_cfg.c_str()));
  AcqResult ca = runMode(capped, scans);
  TEST_EQUAL(ca.ms2_cmds.size() > 0, true)
  int capped_widest = 0;
  for (const auto& c : ca.ms2_cmds) capped_widest = std::max(capped_widest, 1 + notchesForStage(c, 0).second);
  TEST_EQUAL(capped_widest, 3)

  // The wire equals the record: anchor then notches, in order, is exactly the acquisition set.
  int compared = 0, mismatches = 0;
  forEachSurvey(capped_cfg, scans, [&](const Config& cfg, PrecursorSelection& ps, int n, int) {
    ScanCommandQueue q(cfg);
    const auto& sel = ps.selectedPeakGroups();
    for (int i = 0; i < n; ++i)
    {
      const std::vector<int>& set = ps.triggerAcquisitionSets()[i];
      ScanCommand cmd = q.buildMS2(sel[i], ps.triggerCharges()[i], cfg.level(2).scans[0], 2, 1, &set);
      std::vector<int> wire{cmd.stages[0].charge_state};
      auto [first, count] = notchesForStage(cmd, 0);
      for (int k = 0; k < count; ++k) wire.push_back(first[k].charge_state);
      ++compared;
      if (wire != set) ++mismatches;
    }
  });
  TEST_EQUAL(compared > 0, true)
  TEST_EQUAL(mismatches, 0)
}
END_SECTION

// CM-11 (ADR-0059): under single, the anchor is the first charge of charge_selection's ranking.
// snr and qscore are defined by the PeakGroup accessors they used before the key existed, so those
// two are checked as "the accessor, and nothing ranks higher"; intensity and purity against the
// oracle's argmax. The non-vacuity guard needs intensity and purity each to disagree with snr on
// at least one species, or every assertion here could pass with the rule ignored.
START_SECTION(anchor_follows_charge_selection)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  std::map<std::string, std::map<std::pair<int, long long>, int>> anchors;   // rule -> (survey, mass) -> z
  for (const std::string& rule : kRules)
  {
    int checked = 0, wrong = 0;
    forEachSurvey(cfgFor("single", "none", "", 10, false, 0, 1, rule), scans,
                  [&](const Config& cfg, PrecursorSelection& ps, int n, int si) {
      const auto& sel = ps.selectedPeakGroups();
      const MSSpectrum& raw = ps.deconvolvedMS1().getOriginalSpectrum();
      for (int i = 0; i < n; ++i)
      {
        const PeakGroup& pg = sel[i];
        const int z = ps.triggerCharges()[i];
        const auto cands = oracleCharges(pg, raw, cfg.level(1).min_charge);
        anchors[rule][{si, std::llround(pg.getMonoMass())}] = z;
        ++checked;
        bool ok;
        if (rule == "snr" || rule == "qscore")
        {
          ok = (rule == "snr") ? z == pg.getRepAbsCharge() : z == pg.getBestQScoreCharge();
          if (rule == "qscore") ok = ok && !pg.getAllQscores().empty();
          double kz = -1e300;
          for (const auto& c : cands) if (c.z == z) kz = oracleKey(c, rule);
          for (const auto& c : cands) if (oracleKey(c, rule) > kz + 1e-9) ok = false;
        }
        else
        {
          ok = !cands.empty() && z == oracleOrder(cands, rule).front();
        }
        if (!ok)
        {
          ++wrong;
          std::cout << "[CM-11] rule=" << rule << " survey=" << si << " mass=" << pg.getMonoMass()
                    << " engine z=" << z << std::endl;
        }
      }
    });
    TEST_EQUAL(checked > 0, true)
    TEST_EQUAL(wrong, 0)
  }

  for (const std::string& rule : {std::string("intensity"), std::string("purity")})
  {
    int differs = 0;
    for (const auto& kv : anchors[rule])
    {
      auto it = anchors["snr"].find(kv.first);
      if (it != anchors["snr"].end() && it->second != kv.second) ++differs;
    }
    if (differs == 0) std::cout << "[CM-11] " << rule << " never disagreed with snr -- vacuous" << std::endl;
    TEST_EQUAL(differs > 0, true)
  }
}
END_SECTION

// CM-12 (ADR-0059): under separate the emitted charges of one species are the anchor followed by
// the rule's next members among the SNR-gated charges, cut at max_charge_states (3).
START_SECTION(separate_set_follows_charge_selection)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  for (const std::string& rule : kRules)
  {
    int groups = 0, wrong = 0;
    forEachSurvey(cfgFor("separate", "none", "", 10, false, /*max_charge_states=*/3, 1, rule), scans,
                  [&](const Config& cfg, PrecursorSelection& ps, int n, int si) {
      const auto& sel = ps.selectedPeakGroups();
      const auto& trig = ps.triggerCharges();
      const MSSpectrum& raw = ps.deconvolvedMS1().getOriginalSpectrum();
      for (int i = 0; i < n;)
      {
        // One species' emissions are consecutive and share its PeakGroup.
        int j = i;
        std::vector<int> emitted;
        while (j < n && sel[j].getMonoMass() == sel[i].getMonoMass()) emitted.push_back(trig[j++]);
        const PeakGroup& pg = sel[i];
        const int anchor = emitted.front();

        std::vector<OracleCharge> gated;
        for (const auto& c : oracleCharges(pg, raw, 0))
          if (c.z != anchor && c.snr >= cfg.targeting().snr_threshold) gated.push_back(c);
        std::vector<int> expected{anchor};
        for (int z : oracleOrder(gated, rule))
          if ((int)expected.size() < 3) expected.push_back(z);

        bool ok = emitted == expected;
        if (rule == "intensity" || rule == "purity")
          ok = ok && anchor == oracleOrder(oracleCharges(pg, raw, cfg.level(1).min_charge), rule).front();
        ++groups;
        if (!ok)
        {
          ++wrong;
          std::cout << "[CM-12] rule=" << rule << " survey=" << si << " mass=" << pg.getMonoMass() << " emitted=";
          for (int z : emitted) std::cout << " " << z;
          std::cout << " expected=";
          for (int z : expected) std::cout << " " << z;
          std::cout << std::endl;
        }
        i = j;
      }
    });
    TEST_EQUAL(groups > 0, true)
    TEST_EQUAL(wrong, 0)
  }
}
END_SECTION

// CM-13 (ADR-0059): an authored charge set only narrows -- the rule chooses among the named charges.
// Checked on the species' FIRST emission of each run, where no named charge is spent yet.
START_SECTION(authored_anchor_follows_charge_selection)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  for (const std::string& rule : kRules)
  {
    bool seen = false;
    int checked = 0, wrong = 0;
    forEachSurvey(cfgFor("single", "inclusion", charge_set_list, 10, false, 0, 1, rule), scans,
                  [&](const Config& cfg, PrecursorSelection& ps, int n, int si) {
      const auto& sel = ps.selectedPeakGroups();
      const MSSpectrum& raw = ps.deconvolvedMS1().getOriginalSpectrum();
      for (int i = 0; i < n && !seen; ++i)
      {
        if (std::llround(sel[i].getMonoMass()) != authored_mass) continue;
        seen = true;
        std::vector<OracleCharge> named;
        for (const auto& c : oracleCharges(sel[i], raw, cfg.level(1).min_charge))
          if (authored_charges.count(c.z)) named.push_back(c);
        const int z = ps.triggerCharges()[i];
        ++checked;
        if (named.empty() || z != oracleOrder(named, rule).front())
        {
          ++wrong;
          std::cout << "[CM-13] rule=" << rule << " survey=" << si << " engine z=" << z << std::endl;
        }
      }
    });
    TEST_EQUAL(checked, 1)
    TEST_EQUAL(wrong, 0)
  }
}
END_SECTION

// ---- Survey re-planning: the selection journal ---------------------------------------------------
//
// The journal is only worth anything if replaying it IS the in-place write sequence. Each survey is
// therefore compared against itself: A is acquisition memory as production left it, B is the survey's
// entry state (rollbackAcquisitionMemory) plus every journaled stamp in write order. Replay restores A
// exactly, so the run continues on production's own trajectory and every later survey is compared too.
//
// The modes cover every write site: single / separate / multiplexed (the mass-keyed writes and the
// species fan-out), an authored charge set (the per-charge record, ADR-0028), in_depth (two iterations
// over one survey, the intra-survey read-after-write the journal must not disturb) and legacy .log
// inclusion (the charge erase in target_mass_charge_map_, the one write outside SurveyMemory).
START_SECTION(journal_replay_reproduces_in_place_memory)
{
  const auto cytc = loadTsvScans(ms1_tsv_path);
  const auto standard = loadTsvScans(ms1_standard_path);
  ABORT_IF(cytc.empty() || standard.empty())

  // test_target_log.log names three small masses ms1_cytc never carries, so on its own the legacy erase
  // could pass by never firing. This log names the fixture's own species at every charge it resolves.
  const std::string cytc_log = "chargemodes_journal_cytc_target.log";
  {
    std::ofstream lg(cytc_log);
    lg << "MS1 RT 0.05 (scan)\n";
    for (int z = 8; z <= 20; ++z) lg << "Mass 12351.3\tScore=0.9\tZ=" << z << "\n";
  }

  struct Case
  {
    std::string name;
    std::string cfg;
    const std::vector<ScanData>* scans;
  };
  const std::vector<Case> cases = {
    {"single", cfgFor("single"), &cytc},
    {"separate", cfgFor("separate"), &cytc},
    {"multiplexed", cfgFor("multiplexed"), &cytc},
    {"authored inclusion", cfgFor("single", "inclusion", charge_set_list, 5), &cytc},
    {"in_depth", withTargetLog(cfgFor("single", "in_depth"), target_log_path), &standard},
    {"in_depth cytc", withTargetLog(cfgFor("single", "in_depth"), target_log_path), &cytc},
    {"legacy inclusion", withReferenceLogTargets(withTargetLog(cfgFor("single", "inclusion"), target_log_path)), &standard},
    {"legacy inclusion cytc", withReferenceLogTargets(withTargetLog(cfgFor("single", "inclusion"), cytc_log)), &cytc},
  };

  int legacy_erases = 0;
  for (const auto& c : cases)
  {
    int stamps_seen = 0, rollbacks_that_undid_something = 0, mismatches = 0;
    forEachSurvey(c.cfg, *c.scans, [&](const Config&, PrecursorSelection& ps, int, int si) {
      const auto in_place = FLASHIdaTestAccess::surveyMemory(ps);
      const auto in_place_targets = FLASHIdaTestAccess::targetMassCharges(ps);
      const std::vector<MemoryStamp> stamps = allStampsBySeq(ps.journal());

      stamps_seen += (int)stamps.size();
      for (size_t i = 1; i < stamps.size(); ++i) TEST_EQUAL(stamps[i - 1].seq < stamps[i].seq, true)
      for (const auto& s : stamps)
        if (s.kind == MemoryStamp::Kind::LegacyTargetChargeErase) ++legacy_erases;

      ps.rollbackAcquisitionMemory();
      if (!(FLASHIdaTestAccess::surveyMemory(ps) == in_place)
          || FLASHIdaTestAccess::targetMassCharges(ps) != in_place_targets)
        ++rollbacks_that_undid_something;

      ps.applyStamps(stamps);
      const bool same = FLASHIdaTestAccess::surveyMemory(ps) == in_place
                        && FLASHIdaTestAccess::targetMassCharges(ps) == in_place_targets;
      if (!same)
      {
        ++mismatches;
        std::cout << "[JOURNAL] " << c.name << " survey=" << si << ": replay does not reproduce the in-place memory"
                  << std::endl;
      }
    });
    std::cout << "[JOURNAL] " << c.name << " stamps=" << stamps_seen << " rollbacks=" << rollbacks_that_undid_something
              << std::endl;
    TEST_EQUAL(c.name + " mismatches=" + std::to_string(mismatches), c.name + " mismatches=0")
    // Non-vacuity: a rollback and a replay that both did nothing would agree with anything.
    TEST_EQUAL(c.name + (stamps_seen > 0 ? " journaled" : " journaled nothing"), c.name + " journaled")
    TEST_EQUAL(c.name + (rollbacks_that_undid_something > 0 ? " rolled back" : " never rolled back"),
               c.name + " rolled back")
  }
  TEST_EQUAL(legacy_erases > 0, true)
}
END_SECTION

// The rollback restores the ACQUISITION containers and nothing else. window_id_ and the id_*_map_
// window bookkeeping must keep advancing -- they are never a decision input -- so rolling them back
// with the rest would hand the next survey ids the last one already issued.
//
// The entry state is read the only way the class offers: a survey whose spectrum selects nothing runs
// the same expiry pass at the same retention time and writes nothing after it.
START_SECTION(rollback_restores_exactly_the_acquisition_containers)
{
  const auto scans = loadTsvScans(ms1_ecoli_path);
  ABORT_IF(scans.size() < 2)

  const std::string cfg_json = cfgFor("single");
  Config cfg{cfg_json};
  Deconvolution deconv(cfg, cfg.toleranceList());
  PrecursorSelection ps(cfg, deconv);

  const auto& s0 = scans[0];
  const auto& s1 = scans[1];
  TEST_EQUAL(ps.filterAndRank(s0.mzs.data(), s0.ints.data(), (int)s0.mzs.size(), s0.rt, 1, 0.0) > 0, true)

  const std::vector<double> lone_mz{500.0}, lone_int{1.0};
  TEST_EQUAL(ps.filterAndRank(lone_mz.data(), lone_int.data(), 1, s1.rt, 1, 0.0), 0)
  const auto entry = FLASHIdaTestAccess::surveyMemory(ps);
  const auto entry_targets = FLASHIdaTestAccess::targetMassCharges(ps);
  TEST_EQUAL(entry.all_mass_rt_map_.empty(), false)  // the first survey's writes survived: a non-trivial entry

  ps.filterAndRank(s1.mzs.data(), s1.ints.data(), (int)s1.mzs.size(), s1.rt, 1, 0.0);
  TEST_EQUAL(ps.journal().species.size() >= 2, true)
  const auto after = FLASHIdaTestAccess::surveyMemory(ps);
  TEST_EQUAL(after == entry, false)
  TEST_EQUAL(after.window_id_ > entry.window_id_, true)

  ps.rollbackAcquisitionMemory();
  const auto rolled = FLASHIdaTestAccess::surveyMemory(ps);

  TEST_EQUAL(rolled.tqscore_exceeding_mz_rt_map_ == entry.tqscore_exceeding_mz_rt_map_, true)
  TEST_EQUAL(rolled.tqscore_exceeding_mass_rt_map_ == entry.tqscore_exceeding_mass_rt_map_, true)
  TEST_EQUAL(rolled.all_mass_rt_map_ == entry.all_mass_rt_map_, true)
  TEST_EQUAL(rolled.mass_qscore_map_ == entry.mass_qscore_map_, true)
  TEST_EQUAL(rolled.authored_acquired_rt_map_ == entry.authored_acquired_rt_map_, true)
  TEST_EQUAL(FLASHIdaTestAccess::targetMassCharges(ps) == entry_targets, true)

  // Bookkeeping is NOT rolled back.
  TEST_EQUAL(rolled.window_id_, after.window_id_)
  TEST_EQUAL(rolled.id_mass_map_ == after.id_mass_map_, true)
  TEST_EQUAL(rolled.id_mz_map_ == after.id_mz_map_, true)
  TEST_EQUAL(rolled.id_qscore_map_ == after.id_qscore_map_, true)
  TEST_EQUAL(rolled.id_charge_map_ == after.id_charge_map_, true)
}
END_SECTION

// A non-MS1 filterAndRank resets the journal and takes no snapshot, so the entry state still on hand belongs to
// the PREVIOUS survey. A rollback after it must therefore do nothing, or it would undo a survey nobody asked to undo.
START_SECTION(rollback_after_a_non_ms1_call_is_a_noop)
{
  const auto scans = loadTsvScans(ms1_ecoli_path);
  ABORT_IF(scans.empty())

  Config cfg{cfgFor("single")};
  Deconvolution deconv(cfg, cfg.toleranceList());
  PrecursorSelection ps(cfg, deconv);

  const auto& s0 = scans[0];
  TEST_EQUAL(ps.filterAndRank(s0.mzs.data(), s0.ints.data(), (int)s0.mzs.size(), s0.rt, 1, 0.0) > 0, true)
  const auto after_survey = FLASHIdaTestAccess::surveyMemory(ps);
  const auto after_survey_targets = FLASHIdaTestAccess::targetMassCharges(ps);
  TEST_EQUAL(after_survey.all_mass_rt_map_.empty(), false)  // the survey wrote: a rollback that ran would show

  TEST_EQUAL(ps.filterAndRank(s0.mzs.data(), s0.ints.data(), (int)s0.mzs.size(), s0.rt, 2, 0.0), 0)
  ps.rollbackAcquisitionMemory();

  TEST_EQUAL(FLASHIdaTestAccess::surveyMemory(ps) == after_survey, true)
  TEST_EQUAL(FLASHIdaTestAccess::targetMassCharges(ps) == after_survey_targets, true)
}
END_SECTION

// ADR-0036: the PeakGroups of one split envelope are ONE Precursor. Their stamps therefore land in one
// SpeciesStamps, keyed by nominal mass, with every charge the survey acquired across all of them --
// which is what lets a later push commit the species once. The split is synthetic (splitEnvelope).
START_SECTION(siblings_share_one_species_stamp_set)
{
  std::unique_ptr<Config> config;
  std::unique_ptr<Deconvolution> deconv;
  std::unique_ptr<PrecursorSelection> ps;
  std::string log;
  const int n = selectSplit(cfgFor("single", "inclusion", wide_set_list, /*max_precursors=*/2, /*strict_inclusion=*/true),
                            config, deconv, ps, log);

  std::map<int, std::set<int>> charges_by_species;
  std::map<int, int> rows_by_species;
  for (int i = 0; i < n; ++i)
  {
    const int nm = SpectralDeconvolution::getNominalMass(ps->selectedPeakGroups()[i].getMonoMass());
    charges_by_species[nm].insert(ps->triggerCharges()[i]);
    ++rows_by_species[nm];
  }
  TEST_EQUAL(rows_by_species.size(), 1)   // one species...
  ABORT_IF(rows_by_species.size() != 1)
  const int nominal = rows_by_species.begin()->first;
  TEST_EQUAL(rows_by_species[nominal], 2) // ...acquired through both of its PeakGroups

  TEST_EQUAL(ps->journal().species.count(nominal), 1)
  const auto it = ps->journal().species.find(nominal);
  ABORT_IF(it == ps->journal().species.end())
  std::set<int> journaled;
  for (const auto& pc : it->second.per_charge) journaled.insert(pc.first);
  TEST_EQUAL(journaled == charges_by_species[nominal], true)
  TEST_EQUAL(it->second.nominal_mass, nominal)
}
END_SECTION

// A candidate that clears exclusion and is then refused still refreshes all_mass_rt_map_ -- that is
// what keeps its stored qscore alive through expiry. No unit owns it, so its stamp is SURVEY-level.
//
// The refusal is ScoreNotBetter: the species was acquired on an earlier survey at a score this one
// does not beat. tqscore_threshold is lifted above every qscore so the mass bars never arm and the
// qscore ledger is what decides.
START_SECTION(refused_candidate_stamp_is_survey_level)
{
  auto scans = loadTsvScans(ms1_tsv_path);
  ABORT_IF(scans.empty())

  std::string cfg = cfgFor("single");
  const std::string tq = "\"tqscore_threshold\": 0.3";
  cfg.replace(cfg.find(tq), tq.size(), "\"tqscore_threshold\": 1.0");

  int refused = 0, admitted = 0;
  // SurveyMemory is private; its type is reached through the test seam, as PlanUnit_test does.
  decltype(FLASHIdaTestAccess::surveyMemory(std::declval<const PrecursorSelection&>())) before;  // as the previous survey left it
  bool have_before = false;
  forEachSurvey(cfg, scans, [&](const Config&, PrecursorSelection& ps, int, int) {
    const SelectionJournal& j = ps.journal();
    for (const auto& s : j.survey_level)
    {
      if (s.kind != MemoryStamp::Kind::AllMassRt) continue;
      ++refused;
      TEST_EQUAL(j.species.count(s.key), 0)
      // ...refused against a bar an EARLIER survey set, i.e. the second candidate of this mass.
      TEST_EQUAL(have_before && before.mass_qscore_map_.count(s.key) > 0, true)
    }
    for (const auto& kv : j.species)
    {
      ++admitted;
      int all_mass = 0;
      for (const auto& s : kv.second.species_level)
        if (s.kind == MemoryStamp::Kind::AllMassRt && s.key == kv.first) ++all_mass;
      TEST_EQUAL(all_mass, 1)  // an admitted candidate's stamp goes with its species, not the survey
    }
    before = FLASHIdaTestAccess::surveyMemory(ps);
    have_before = true;
  });
  TEST_EQUAL(admitted > 0, true)
  TEST_EQUAL(refused > 0, true)
}
END_SECTION

// CM-12: max_precursors is a budget per SURVEY under "separate" too, when the survey walks the
// candidate list twice. In-depth mode does: iteration 0, then the back-fill iteration 1, which may
// only fill what iteration 0 left. "separate" skips the acquisition-count guard that bounds single and
// multiplexed across both passes, and its own species count restarted at 0 on the second pass, so
// with six eligible species and max_precursors 2 it fragmented FOUR species (16 scans) per survey.
//
// Multiplexed is the reference: ADR-0016 says the two modes acquire the same species and charge
// sets and differ only in scan count. No fixture: the survey is built from synthetic PeakGroups.
START_SECTION(separate_respects_the_survey_budget_in_depth_mode)
{
  const std::vector<std::pair<double, int>> sep = selectSynthetic(cfgFor("separate", "in_depth", "", /*max_precursors=*/2));
  const std::vector<std::pair<double, int>> mux = selectSynthetic(cfgFor("multiplexed", "in_depth", "", 2));
  const std::vector<std::pair<double, int>> one = selectSynthetic(cfgFor("single", "in_depth", "", 2));

  // Non-vacuity: "separate" really fans out here (four charges per species), and the budget is
  // really contended (six eligible species against two slots).
  TEST_EQUAL(sep.size() > speciesOf(sep).size(), true)
  TEST_EQUAL((int)speciesOf(one).size(), 2)

  if (speciesOf(sep).size() > 2)
    std::cout << "[CM-12] separate acquired " << speciesOf(sep).size() << " species in one survey "
              << "against max_precursors 2 -- the back-fill iteration spent the budget again" << std::endl;
  TEST_EQUAL((int)speciesOf(sep).size(), 2)
  TEST_EQUAL((int)sep.size(), 8)   // 2 species x charges 8..11
  TEST_EQUAL(speciesOf(sep) == speciesOf(mux), true)
  TEST_EQUAL(speciesOf(sep) == speciesOf(one), true)
  TEST_EQUAL(speciesOf(sep).count(std::llround(synMass(0))) == 1 && speciesOf(sep).count(std::llround(synMass(1))) == 1, true)
}
END_SECTION

// CM-13: the same budget under NON-STRICT inclusion, the other configuration with two passes: phase 0
// takes the targets, phase 1 everything else. One target (rank 4) and max_precursors 2 leave ONE slot
// for phase 1. Under the defect "separate" restarted its count in phase 1 and took two non-targets,
// three species in all.
START_SECTION(separate_respects_the_survey_budget_under_non_strict_inclusion)
{
  const std::vector<double> targets {synMass(4)};
  const std::vector<std::pair<double, int>> sep =
      selectSynthetic(cfgFor("separate", "inclusion", "", /*max_precursors=*/2, /*strict_inclusion=*/false), targets);
  const std::vector<std::pair<double, int>> mux =
      selectSynthetic(cfgFor("multiplexed", "inclusion", "", 2, false), targets);

  TEST_EQUAL(sep.size() > speciesOf(sep).size(), true)   // non-vacuity: it fans out

  if (speciesOf(sep).size() > 2)
    std::cout << "[CM-13] separate acquired " << speciesOf(sep).size() << " species in one survey "
              << "against max_precursors 2 -- the non-target phase spent the budget again" << std::endl;
  TEST_EQUAL((int)speciesOf(sep).size(), 2)
  TEST_EQUAL(speciesOf(sep) == speciesOf(mux), true)
  // The target comes first (phase 0), and the one slot left goes to the best non-target.
  TEST_EQUAL(sep.empty() ? 0LL : std::llround(sep.front().first), std::llround(synMass(4)))
  TEST_EQUAL((int)speciesOf(sep).count(std::llround(synMass(0))), 1)
}
END_SECTION

// CM-14: the single-pass configurations are untouched -- "separate" with no targeting already spent
// exactly max_precursors species, and must still fan out over all of them.
START_SECTION(separate_single_pass_budget_is_unchanged)
{
  const std::vector<std::pair<double, int>> sep = selectSynthetic(cfgFor("separate", "none", "", /*max_precursors=*/2));
  TEST_EQUAL((int)speciesOf(sep).size(), 2)
  TEST_EQUAL((int)sep.size(), 8)
}
END_SECTION

END_TEST
