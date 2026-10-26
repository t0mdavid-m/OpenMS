// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Survey re-planning: a survey's MS2 decisions stay provisional while they wait in FLASHIda's queue,
// every newer survey at the same CV re-plans them, and only the drain's dequeue commits.
//
// The survey hands the queue its whole plan for one CV in one call (replacePlan), one PLAN UNIT per
// species, and prints `[REPLAN] cv= survey= kept= dropped= added= adjusted=`. These sections drive the
// real engine and read that line and the commands the drain hands out.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>

#include "FLASHIda_TestAccess.h"
#include "FLASHIda_TestHelpers.h"  // ground-truth harness: ScanData/loadTsvScans/runInterleaved, TSVFile

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <streambuf>
#include <string>
#include <thread>
#include <vector>

using namespace OpenMS;

namespace
{
  /// One `[REPLAN]` line, as the survey printed it.
  struct ReplanLine
  {
    std::string survey;
    int kept = -1, dropped = -1, added = -1, adjusted = -1;
  };

  /// Every `[REPLAN]` line of @p out, in print order.
  std::vector<ReplanLine> replanLines(const std::string& out)
  {
    std::vector<ReplanLine> lines;
    std::istringstream in(out);
    std::string line;
    while (std::getline(in, line))
    {
      if (line.rfind("[REPLAN] ", 0) != 0) continue;
      ReplanLine r;
      std::istringstream tok(line.substr(9));
      std::string kv;
      while (tok >> kv)
      {
        const auto eq = kv.find('=');
        if (eq == std::string::npos) continue;
        const std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
        if (k == "survey") r.survey = v;
        else if (k == "kept") r.kept = std::atoi(v.c_str());
        else if (k == "dropped") r.dropped = std::atoi(v.c_str());
        else if (k == "added") r.added = std::atoi(v.c_str());
        else if (k == "adjusted") r.adjusted = std::atoi(v.c_str());
      }
      lines.push_back(r);
    }
    return lines;
  }

  /// The C++-schema twin of test-data/configs/method_inclusion_charge_set.json (the C++ Config cannot read
  /// the C# schema): inclusion of cytC at the authored charge set 10;13;16, `single`, five precursors.
  std::string inclusionChargeSetConfig(const std::string& log_dir,
                                       const std::string& list = "inclusion_cytc_charges.txt",
                                       int max_precursors = 5, bool strict_inclusion = false)
  {
    return std::string(R"JSON({
  "deconvolution": {
    "score_threshold": 0.0, "tqscore_threshold": 0.0,
    "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
  },
  "flashtnt": { "min_length": 3, "max_length": 8 },
  "faims": { "cv_values": [], "max_cv_skip": 0 },
  "scheduling": {
    "cycle_time": { "enabled": false, "value_ms": 60000 },
    "scan_timeout": { "enabled": false, "value_ms": 30000 },
    "agc_interval_seconds": 9999999
  },
  "files": { "target_logs": [], "fasta": "", "inclusion_list": "../../FlashIDA/test-data/configs/)JSON") + list + R"JSON(", "ptm_list": "" },
  "precursor_selection": {
    "rt_window": 180, "targeting": "inclusion", "precursor_charges": "single",
    "consider_all_charges": false, "strict_inclusion": )JSON" + (strict_inclusion ? "true" : "false") + R"JSON(,
    "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": )JSON" + std::to_string(max_precursors) + R"JSON(
  },
  "characterization": { "mode": "off", "max_targets": 3 },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
    "ms2": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 500000,
             "max_it": 100, "microscans": 4, "activation": "ETD", "reaction_time": 7,
             "reagent_agc_target": 700000, "reagent_max_it": 200 }
  },
  "tagging": {},
  "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
  "runtime": { "log_dir": ")JSON" + log_dir + R"JSON(" }
})JSON";
  }

  /// Run @p fn with std::cout captured; returns what it printed.
  template <typename Fn>
  std::string captureCout(Fn fn)
  {
    std::ostringstream captured;
    std::streambuf* const real_cout = std::cout.rdbuf(captured.rdbuf());
    try
    {
      fn();
    }
    catch (...)
    {
      std::cout.rdbuf(real_cout);
      throw;
    }
    std::cout.rdbuf(real_cout);
    return captured.str();
  }

  /// One survey's plan as the run shows it: the `[REPLAN]` line's counts, and the nominal masses of the
  /// MS2 commands the survey parented.
  struct SurveyPlan
  {
    int added = -1;
    int dropped = -1;
    size_t commands = 0;
    std::set<int> nominals;
  };

  /// Drive @p cfg over cytC and pair every `[REPLAN]` line with the survey's MS2 commands. The drive is
  /// depth 1, so every MS2 a survey planned is dequeued before the next survey: the survey's plan is
  /// exactly the MS2 commands whose parent it is.
  std::vector<SurveyPlan> plansOf(const std::string& cfg)
  {
    auto scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
    AcqResult a;
    const std::string out = captureCout([&] {
      FLASHIda ida(const_cast<char*>(cfg.c_str()));
      a = runInterleaved(&ida, scans, std::vector<ScanData>{}, nullptr, 4000);
    });
    std::map<std::string, std::vector<const ScanCommand*>> by_survey;
    for (const auto& c : a.ms2_cmds) by_survey[std::string(c.parent_scan_id)].push_back(&c);

    std::vector<SurveyPlan> plans;
    for (const auto& r : replanLines(out))
    {
      SurveyPlan p;
      p.added = r.added;
      p.dropped = r.dropped;
      const auto it = by_survey.find(r.survey);
      if (it != by_survey.end())
      {
        p.commands = it->second.size();
        for (const ScanCommand* c : it->second) p.nominals.insert(SpectralDeconvolution::getNominalMass(c->mono_mass));
      }
      plans.push_back(p);
    }
    return plans;
  }

  /// The decision_survey_id cell of every scan_commands.tsv row in @p log_dir, keyed by the row's tracking_id (both
  /// encoded). The drain consumes the engine's decision record as it writes the row, so after the drain the log is
  /// where the decision survey is read -- and the logged value is the decided one.
  std::map<std::string, std::string> decisionSurveys(const std::string& log_dir)
  {
    const TSVFile tsv = TSVFile::parse(log_dir + "/scan_commands.tsv");
    std::map<std::string, std::string> out;
    for (const auto& row : tsv.rows) out[cell(tsv, row, "tracking_id")] = cell(tsv, row, "decision_survey_id");
    return out;
  }

  /// The logged decision survey of the command with tracking id @p id; empty when it has no row or no survey decided it.
  std::string decidedBy(const std::map<std::string, std::string>& surveys, int id)
  {
    const auto it = surveys.find(ScanCommandQueue::encode(id));
    return it == surveys.end() ? std::string() : it->second;
  }

  // ---- The bought R (Task 10): a quantification verdict buys an identification R that waits ----------------

  /// The C++-schema twin of test-data/configs/method_quant_welch.json (@p charges "separate", the quant_welch
  /// golden's method) or of method_quant_multiplexed.json ("multiplexed"): one species per survey, two to three
  /// charges, a Q per charge or one Q for the set, and the verdict buys an ETD identification R.
  std::string quantTwinConfig(const std::string& charges)
  {
    const bool separate = charges == "separate";
    return std::string(R"JSON({
  "deconvolution": {
    "score_threshold": 0.0, "tqscore_threshold": 0.9,
    "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
  },
  "flashtnt": { "min_length": 3, "max_length": 8 },
  "conditional_ms2": false,
  "faims": { "cv_values": [], "max_cv_skip": 0 },
  "scheduling": {
    "cycle_time": { "enabled": false, "value_ms": 60000 },
    "scan_timeout": { "enabled": false, "value_ms": 30000 },
    "agc_interval_seconds": 9999999
  },
  "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
  "precursor_selection": {
    "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false,
    "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 1,
    "precursor_charges": ")JSON") + charges + R"JSON(", "min_charge_states": 2, "max_charge_states": 3
  },
  "characterization": { "mode": "off", "max_targets": 3 },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
    "ms2": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 500000,
             "max_it": 100, "microscans": 4, "activation": "ETD", "reaction_time": 7,
             "reagent_agc_target": 700000, "reagent_max_it": 200 },
    "ms2_quant": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 500000,
                   "max_it": 100, "microscans": 4, "activation": "HCD", "collision_energy": 30 }
  },
  "tagging": {},
  "quantification": {
    "enabled": true, "labelling": "tmt6plex", "reporter_mz_tol": 0.002,)JSON"
           + (separate ? std::string(R"JSON(
    "test": "welch_and_fold_change", "p_value_threshold": 0.05, "fold_change_threshold": 1.5,
    "empty_condition": "incomplete", "aggregate": "median", "min_usable_scans": 1,)JSON")
                       : std::string(R"JSON(
    "fold_change_threshold": 1.4,)JSON"))
           + R"JSON(
    "conditions": [ { "name": "control", "channels": ["126", "127", "128"] },
                    { "name": "treated", "channels": ["129", "130", "131"] } ]
  }
})JSON";
  }

  /// A bought R, caught while it waits: the command as queued, the charges it was registered with, the
  /// survey that decided its Q, and every charge the Qs fed back to the engine isolated (sorted, unique).
  struct WaitingR
  {
    bool bought = false;
    ScanCommand r{};
    std::vector<int> covered;
    int q_survey = 0;
    std::vector<int> q_charges;
  };

  /// Feed @p a as the engine's first survey and return its 'Q' scans the two real iTMT spectra of one species by
  /// charge (odd -> z23, even -> z22, as the quant_welch golden does) until the verdict buys an R. Stops BEFORE
  /// the drain could take the R: it is read off the queue, never dequeued.
  WaitingR buyR(FLASHIda& ida, const ScanData& a, const ScanData& q_odd, const ScanData& q_even)
  {
    WaitingR out;
    ScanCommandQueue& q = FLASHIdaTestAccess::queue(ida);
    ScanCommand cmd{};
    ida.getNextScanCommand(cmd); // the idle survey
    if (cmd.msn_level != 1) return out;
    out.q_survey = cmd.scan_id;
    ida.processScan(a.mzs.data(), a.ints.data(), (int)a.mzs.size(), a.rt, 1, cmd.scan_description, cmd.faims_cv,
                    instrumentScanNumberOf(a));
    for (int it = 0; it < 20; ++it)
    {
      for (const ScanCommand& c : FLASHIdaTestAccess::queuedCommands(q))
      {
        const auto covered = FLASHIdaTestAccess::adjustableResolvedCharges(q, c.scan_id);
        if (!covered) continue;
        out.bought = true;
        out.r = c;
        out.covered = *covered;
        std::sort(out.q_charges.begin(), out.q_charges.end());
        out.q_charges.erase(std::unique(out.q_charges.begin(), out.q_charges.end()), out.q_charges.end());
        return out;
      }
      cmd = ScanCommand{};
      ida.getNextScanCommand(cmd);
      if (cmd.msn_level != 2 || cmd.scan_description[3] != 'Q') return out; // the plan ran out before a buy
      out.q_charges.push_back(std::abs(cmd.stages[0].charge_state));
      for (int j = 0; j < cmd.stage0_notch_count; ++j) out.q_charges.push_back(std::abs(cmd.notches[j].charge_state));
      const ScanData& s = std::abs(cmd.stages[0].charge_state) % 2 == 1 ? q_odd : q_even;
      ida.processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 2, cmd.scan_description, 0.0,
                      instrumentScanNumberOf(s));
    }
    return out;
  }

  /// A survey that reaches the engine while the R still waits -- a cycle-time or CV-transition MS1 the drain
  /// sent ahead of it. Minted and registered as the drain would; returns its tracking id.
  int feedSurvey(FLASHIda& ida, const ScanData& s, double rt)
  {
    ScanCommandQueue& q = FLASHIdaTestAccess::queue(ida);
    ScanCommand ms1 = q.makeMS1();
    ms1.scan_id = q.nextTrackingId();
    ScanCommandQueue::writeMs1Description(ms1, ScanRole::Survey);
    q.registerPending(ms1);
    ida.processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), rt, 1, ms1.scan_description, 0.0, 0);
    return ms1.scan_id;
  }

  /// @p s without the isotope envelope of the species of mono mass @p mass at each charge of @p charges: every peak
  /// between the monoisotopic m/z and 25 Da above it, padded by 0.3 Th. The other charges' windows lie tens of Th
  /// away, so the species stays in the spectrum at every charge not named.
  ScanData withoutCharges(const ScanData& s, double mass, const std::vector<int>& charges)
  {
    ScanData out = s;
    out.mzs.clear();
    out.ints.clear();
    for (size_t i = 0; i < s.mzs.size(); ++i)
    {
      bool inside = false;
      for (int z : charges)
      {
        const double lo = mass / z + 1.007276 - 0.3;
        const double hi = (mass + 25.0) / z + 1.007276 + 0.3;
        if (s.mzs[i] >= lo && s.mzs[i] <= hi) inside = true;
      }
      if (inside) continue;
      out.mzs.push_back(s.mzs[i]);
      out.ints.push_back(s.ints[i]);
    }
    return out;
  }

  /// The species' charges on deconvolution @p d as charge_selection "snr" (the default) ranks them, computed here
  /// from the PeakGroups rather than through the engine: within ±1 nominal mass, isolatable, charge SNR >= 1,
  /// restricted to @p allowed when non-empty; SNR descending, then intensity, then the lower charge.
  std::vector<int> rankedCharges(const DeconvolvedSpectrum& d, int nominal, const std::vector<int>& allowed)
  {
    struct C { int z; double snr; double intensity; };
    std::vector<C> cs;
    for (const PeakGroup& pg : d)
    {
      if (std::abs(SpectralDeconvolution::getNominalMass(pg.getMonoMass()) - nominal) > 1) continue;
      const auto [min_z, max_z] = pg.getAbsChargeRange();
      for (int z = min_z; z <= max_z; ++z)
      {
        const auto [lo, hi] = pg.getMzRange(z);
        if (hi <= lo || pg.getChargeSNR(z) < 1.0) continue;
        if (!allowed.empty() && std::find(allowed.begin(), allowed.end(), z) == allowed.end()) continue;
        cs.push_back({z, pg.getChargeSNR(z), pg.getChargeIntensity(z)});
      }
    }
    std::stable_sort(cs.begin(), cs.end(), [](const C& a, const C& b) {
      if (a.snr != b.snr) return a.snr > b.snr;
      if (a.intensity != b.intensity) return a.intensity > b.intensity;
      return a.z < b.z;
    });
    std::vector<int> out;
    for (const C& c : cs)
      if (std::find(out.begin(), out.end(), c.z) == out.end()) out.push_back(c.z);
    return out;
  }

  /// The anchor first, then the stage-0 notches: every charge @p c isolates.
  std::vector<int> isolatedCharges(const ScanCommand& c)
  {
    std::vector<int> zs = {std::abs(c.stages[0].charge_state)};
    for (int j = 0; j < c.stage0_notch_count; ++j) zs.push_back(std::abs(c.notches[j].charge_state));
    return zs;
  }

  /// The queued command with tracking id @p id, read without dequeuing; scan_id 0 when it is not queued.
  ScanCommand queued(FLASHIda& ida, int id)
  {
    for (const ScanCommand& c : FLASHIdaTestAccess::queuedCommands(FLASHIdaTestAccess::queue(ida)))
      if (c.scan_id == id) return c;
    return ScanCommand{};
  }

  // ---- Re-plan semantics (Task 9): the methods, as C++-schema twins --------------------------------------------

  /// @p cfg with the one occurrence of @p from replaced by @p to; throws when @p from is absent, so a twin can never
  /// silently keep the base value it meant to change.
  std::string with(std::string cfg, const std::string& from, const std::string& to)
  {
    const auto p = cfg.find(from);
    if (p == std::string::npos) throw std::runtime_error("config twin: '" + from + "' not found");
    return cfg.replace(p, from.size(), to);
  }

  /// The C++-schema twin of test-data/configs/method_default_topn5.json: standard DDA, five precursors, ETD MS2.
  std::string topn5Config(const std::string& log_dir = "")
  {
    return std::string(R"JSON({
  "deconvolution": {
    "score_threshold": 0.0, "tqscore_threshold": 0.9,
    "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
  },
  "flashtnt": { "min_length": 3, "max_length": 8 },
  "conditional_ms2": false,
  "faims": { "cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
  "scheduling": {
    "cycle_time": { "enabled": false, "value_ms": 60000 },
    "scan_timeout": { "enabled": false, "value_ms": 30000 },
    "agc_interval_seconds": 9999999
  },
  "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
  "precursor_selection": {
    "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false,
    "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 5
  },
  "characterization": { "mode": "off", "max_targets": 3 },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
    "ms2": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 500000,
             "max_it": 100, "microscans": 4, "activation": "ETD", "reaction_time": 7,
             "reagent_agc_target": 700000, "reagent_max_it": 200 },
    "ms3": { "analyzer": "Orbitrap", "first_mass": 200, "last_mass": 2000, "resolution": 240000, "agc_target": 50000000,
             "max_it": 500, "microscans": 8, "activation": "CID", "collision_energy": 25 }
  },
  "tagging": {},
  "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
  "runtime": { "log_dir": ")JSON") + log_dir + R"JSON(" }
})JSON";
  }

  /// method_faims_3cv.json's twin, on FLASHIdaFAIMS_test's three CVs: CV cycling over -40/-50/-60, three precursors.
  std::string faims3cvConfig(const std::string& log_dir = "")
  {
    return with(with(topn5Config(log_dir),R"("cv_values": [])", R"("cv_values": [-40, -50, -60])"), R"("max_precursors": 5)",
                R"("max_precursors": 3)");
  }

  /// method_cytc_separate_charges.json's twin with max_charge_states 3: inclusion of cytC, one precursor, `separate`
  /// -- one MS2 per charge -- and ambiguity characterization of the M-starting cytC proteoform.
  std::string cytcSeparateConfig(const std::string& log_dir = "")
  {
    std::string cfg = with(topn5Config(log_dir),R"("tqscore_threshold": 0.9)", R"("tqscore_threshold": 0.0)");
    cfg = with(cfg, R"("inclusion_list": "")", R"("inclusion_list": "../../FlashIDA/test-data/configs/inclusion_cytc.txt")");
    cfg = with(cfg, R"("targeting": "none")", R"("targeting": "inclusion")");
    cfg = with(cfg, R"("max_precursors": 5)", R"("max_precursors": 1, "precursor_charges": "separate", "max_charge_states": 3)");
    cfg = with(cfg, R"("characterization": { "mode": "off", "max_targets": 3 })",
               R"("characterization": { "mode": "ambiguity", "max_targets": 3, "fragment_charges": "separate",
    "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE" })");
    return with(cfg, R"("activation": "ETD", "reaction_time": 7,
             "reagent_agc_target": 700000, "reagent_max_it": 200 })",
                R"("activation": "HCD", "collision_energy": 30 })");
  }

  /// method_tag_targeting.json's twin: inclusion targeting with no list, the FASTA tags are matched against, and an
  /// HCD conditional follow-up ('C') raised by an MS2 whose tags hit.
  std::string tagTargetingConfig(const std::string& log_dir)
  {
    std::string cfg = with(topn5Config(log_dir), R"("conditional_ms2": false)", R"("conditional_ms2": true)");
    cfg = with(cfg, R"("fasta": "")", R"("fasta": "../../FlashIDA/test-data/configs/test_fasta.fasta")");
    cfg = with(cfg, R"("targeting": "none")", R"("targeting": "inclusion")");
    cfg = with(cfg, R"("tagging": {})", R"("tagging": { "follow_up_scan": "tagging_follow_up" })");
    return with(cfg, R"("activation": "CID", "collision_energy": 25 })",
                R"("activation": "CID", "collision_energy": 25 },
    "additional_ms2": {
      "tagging_follow_up": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000,
                             "agc_target": 500000, "max_it": 100, "microscans": 4, "activation": "HCD", "collision_energy": 30 }
    })");
  }

  /// method_exploration.json's selection side: one precursor, an HCD mass_count sweep over CE 20..40 step 5. The
  /// shared explorationConfig() with its precursor budget at the method's 1; characterization stays off, since the
  /// method's protein is not cytC's.
  std::string explorationTwinConfig()
  {
    return with(explorationConfig(), R"("max_precursors": 3)", R"("max_precursors": 1)");
  }

  /// method_quant_separate.json's twin: quantTwinConfig's multiplexed method (the same quantification block, fold
  /// change 1.4) with the charges acquired `separate` -- a Q per charge of the one species. With @p log_dir set, the
  /// five streams are written there (quantTwinConfig itself carries no runtime block).
  std::string quantSeparateConfig(const std::string& log_dir = "")
  {
    std::string cfg = with(quantTwinConfig("multiplexed"), R"("precursor_charges": "multiplexed")", R"("precursor_charges": "separate")");
    if (log_dir.empty()) return cfg;
    return with(cfg, R"("tagging": {},)", R"("tagging": {},
  "runtime": { "log_dir": ")" + log_dir + R"(" },)");
  }

  int nominalOf(const ScanCommand& c) { return SpectralDeconvolution::getNominalMass(c.mono_mass); }

  /// The ms1_cytc.txt scan that serves as the newer survey B for the separate remainder: scan 191 (rt 21.3724), the
  /// nearest real survey after scan 8 -- index 9 is a noise scan that deconvolves to no PeakGroup. It resolves cytC at
  /// a charge ranking that differs from scan 8's, so a re-derived remainder is told from a kept one.
  const size_t kSeparateB = 12;

  /// The `[REPLAN]` line survey @p id printed; counts -1 when it printed none.
  ReplanLine replanOf(const std::string& out, int id)
  {
    for (const auto& r : replanLines(out))
      if (r.survey == ScanCommandQueue::encode(id)) return r;
    return ReplanLine{};
  }

  /// The MS2 commands of @p cmds, in order.
  std::vector<ScanCommand> ms2Of(const std::vector<ScanCommand>& cmds)
  {
    std::vector<ScanCommand> out;
    for (const auto& c : cmds)
      if (c.msn_level == 2) out.push_back(c);
    return out;
  }

  bool hasId(const std::vector<ScanCommand>& cmds, int id)
  {
    return std::any_of(cmds.begin(), cmds.end(), [id](const ScanCommand& c) { return c.scan_id == id; });
  }

  /// The commands of @p cmds whose nominal mass is within 1 Da of @p nominal -- one species' -- in order.
  std::vector<ScanCommand> ofSpecies(const std::vector<ScanCommand>& cmds, int nominal)
  {
    std::vector<ScanCommand> out;
    for (const auto& c : cmds)
      if (std::abs(nominalOf(c) - nominal) <= 1) out.push_back(c);
    return out;
  }

  /// @p s at retention time @p rt (seconds).
  ScanData at(ScanData s, double rt)
  {
    s.rt = rt;
    return s;
  }

  /// A streambuf that appends under a mutex, for capturing std::cout while a drain thread and the analysis thread
  /// both print.
  struct LockedStringBuf : std::streambuf
  {
    std::mutex m;
    std::string s;
    int_type overflow(int_type c) override
    {
      if (c != traits_type::eof())
      {
        std::lock_guard<std::mutex> lk(m);
        s.push_back(static_cast<char>(c));
      }
      return c;
    }
    std::streamsize xsputn(const char* p, std::streamsize n) override
    {
      std::lock_guard<std::mutex> lk(m);
      s.append(p, static_cast<size_t>(n));
      return n;
    }
  };

  /// One race: survey A plans, a drain thread starts dequeuing @p drain_delay into survey B's processing and keeps
  /// dequeuing every @p drain_period, and the seam dequeues one more command between B's plan and its swap. Until the
  /// seam has run the drain thread takes at most |A's plan| - 1 commands (none when @p drain_before_seam is false),
  /// so the seam always commits one of A's still-waiting commands in the gap. The drain is finished on this thread
  /// afterwards, through the idle survey. Returns every command dequeued, in the order each thread took them (the
  /// seam's first), and what the run printed.
  struct RaceRun
  {
    int a_id = 0, b_id = 0;
    std::vector<ScanCommand> a_plan, seam, drained;
    std::string out;
    int units_left = -1;
    size_t queued_left = 0;
  };

  RaceRun race(const std::string& cfg, const ScanData& ms1, std::chrono::microseconds drain_delay,
               std::chrono::microseconds drain_period, bool drain_before_seam)
  {
    RaceRun r;
    LockedStringBuf buf;
    std::streambuf* const real_cout = std::cout.rdbuf(&buf);
    {
      FLASHIda ida(const_cast<char*>(cfg.c_str()));
      ScanCommandQueue& q = FLASHIdaTestAccess::queue(ida);
      ScanCommand a{};
      ida.getNextScanCommand(a);
      r.a_id = a.scan_id;
      ida.processScan(ms1.mzs.data(), ms1.ints.data(), (int)ms1.mzs.size(), ms1.rt, 1, a.scan_description, a.faims_cv, 0);
      r.a_plan = ms2Of(FLASHIdaTestAccess::queuedCommands(q));

      std::atomic<bool> seam_done{false};
      FLASHIdaTestAccess::onBeforeReplan(ida, [&] {
        ScanCommand c{};
        ida.getNextScanCommand(c);
        r.seam.push_back(c);
        seam_done.store(true);
      });
      ScanCommand b = q.makeMS1();
      b.scan_id = q.nextTrackingId();
      ScanCommandQueue::writeMs1Description(b, ScanRole::Survey);
      q.registerPending(b);
      r.b_id = b.scan_id;

      std::atomic<bool> b_done{false};
      std::vector<ScanCommand> thread_drained;
      const size_t budget = drain_before_seam && !r.a_plan.empty() ? r.a_plan.size() - 1 : 0;
      std::thread drain([&] {
        std::this_thread::sleep_for(drain_delay);
        while (!b_done.load())
        {
          if (!seam_done.load() && thread_drained.size() >= budget)
          {
            std::this_thread::sleep_for(drain_period); // leave the seam one of A's commands to commit
            continue;
          }
          ScanCommand c{};
          ida.getNextScanCommand(c);
          thread_drained.push_back(c);
          std::this_thread::sleep_for(drain_period);
        }
      });
      ida.processScan(ms1.mzs.data(), ms1.ints.data(), (int)ms1.mzs.size(), ms1.rt + 2.0, 1, b.scan_description, b.faims_cv, 0);
      b_done.store(true);
      drain.join();
      FLASHIdaTestAccess::onBeforeReplan(ida, {});

      r.drained = std::move(thread_drained);
      for (int i = 0; i < 1000; ++i)
      {
        ScanCommand c{};
        ida.getNextScanCommand(c);
        r.drained.push_back(c);
        if (c.msn_level == 1 && c.priority == 3) break;
      }
      r.units_left = FLASHIdaTestAccess::planUnitCount(q);
      r.queued_left = FLASHIdaTestAccess::queuedCommands(q).size();
    }
    std::cout.rdbuf(real_cout);
    r.out = buf.s;
    return r;
  }

  /// Every call of a queue push site in @p source, as "<enclosing FLASHIda:: function> | <the call's line, trimmed>",
  /// comments stripped. The enclosing function is the last `FLASHIda::name(` definition opened at column 2.
  std::multiset<std::string> pushSites(const std::string& source)
  {
    static const std::vector<std::string> calls = {"stampAndPush_(", "queue_.push(", "queue_.pushAndDequeue(", "queue_.replacePlan(",
                                                   "queue_.registerAdjustable("};
    std::multiset<std::string> sites;
    std::istringstream in(source);
    std::string line, fn;
    bool in_block_comment = false;
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      std::string code;
      for (size_t i = 0; i < line.size(); ++i)
      {
        if (in_block_comment)
        {
          if (line.compare(i, 2, "*/") == 0)
          {
            in_block_comment = false;
            ++i;
          }
          continue;
        }
        if (line.compare(i, 2, "/*") == 0)
        {
          in_block_comment = true;
          ++i;
          continue;
        }
        if (line.compare(i, 2, "//") == 0) break;
        code += line[i];
      }
      if (code.size() > 2 && code[0] == ' ' && code[1] == ' ' && code[2] != ' ')
      {
        const auto p = code.find("FLASHIda::");
        const auto paren = p == std::string::npos ? std::string::npos : code.find('(', p);
        if (paren != std::string::npos) fn = code.substr(p + 10, paren - p - 10);
      }
      for (const auto& call : calls)
      {
        if (code.find(call) == std::string::npos) continue;
        const auto b = code.find_first_not_of(' ');
        const auto e = code.find_last_not_of(' ');
        sites.insert(fn + " | " + code.substr(b, e - b + 1));
      }
    }
    return sites;
  }
} // namespace

START_TEST(FLASHIda_Replan, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Two PeakGroups of one species in one survey (ADR-0036 siblings) are ONE plan unit: the unit is keyed by
// nominal mass, so a survey adds exactly as many units as it selected distinct nominal masses.
START_SECTION(siblings_form_one_unit)
{
  // method_inclusion_charge_set: cytC at 10;13;16. Its winner carries all three named charges in every
  // survey of this fixture that acquires it, so no sibling is ever admitted -- the assertion holds, but it
  // cannot tell one unit per species from one unit per PeakGroup.
  const auto plans = plansOf(inclusionChargeSetConfig(freshLogDir("replan_siblings")));
  TEST_EQUAL(plans.empty(), false)
  for (const auto& p : plans)
  {
    TEST_EQUAL(p.added, static_cast<int>(p.nominals.size()))
    TEST_EQUAL(p.dropped, 0)
  }

  // So the same assertion again on the WIDE row 12..18 (FLASHIda_ChargeModes_test CM-09's config), whose
  // envelope arrives split: there a survey admits a sibling, and two commands of one nominal mass in one
  // survey must still be one unit.
  const auto split = plansOf(inclusionChargeSetConfig(freshLogDir("replan_siblings_split"), "inclusion_cytc_charges_wide.txt",
                                                      /*max_precursors=*/2, /*strict_inclusion=*/true));
  TEST_EQUAL(split.empty(), false)
  int surveys_with_siblings = 0;
  for (const auto& p : split)
  {
    TEST_EQUAL(p.added, static_cast<int>(p.nominals.size()))
    TEST_EQUAL(p.dropped, 0)
    if (p.commands > p.nominals.size()) ++surveys_with_siblings;
  }
  TEST_EQUAL(surveys_with_siblings > 0, true) // non-vacuity: a sibling was admitted
}
END_SECTION

// At depth 1 every MS2 a survey planned is dequeued before the next survey arrives, so no survey ever meets
// a waiting unit of its CV and re-planning must change nothing: no survey drops, keeps or adjusts a command, and
// every survey-decided MS2 goes out decided by the survey that is its parent. The drive is FLASHIda_Logging_test's
// (method_default's C++ twin, ms1_standard + ms2_hcd_fragment). It used to compare against a capture taken at the
// Push 1 tip; FDmerge's deconvolution changed the selection itself, so the claim is now asserted directly.
START_SECTION(depth1_drive_is_unchanged_without_a_stale_survey)
{
  const std::string log_dir = freshLogDir("replan_depth1");
  const std::string cfg = buildJsonWithLogDir(log_dir);
  auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_standard.txt");
  auto ms2 = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt");
  ABORT_IF(ms1.empty() || ms2.empty())
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    runFullCycle(&ida, ms1, ms2);
  });
  const auto lines = replanLines(out);
  TEST_EQUAL(lines.empty(), false)
  for (const auto& r : lines)
  {
    TEST_EQUAL(r.dropped, 0)
    TEST_EQUAL(r.kept, 0)
    TEST_EQUAL(r.adjusted, 0)
  }

  const TSVFile tsv = TSVFile::parse(log_dir + "/scan_commands.tsv");
  std::map<std::string, std::string> level_of;
  for (const auto& row : tsv.rows) level_of[cell(tsv, row, "tracking_id")] = cell(tsv, row, "ms_level");
  int survey_decided = 0;
  for (const auto& row : tsv.rows)
  {
    const std::string parent = cell(tsv, row, "parent_tracking_id");
    const bool from_survey = cell(tsv, row, "ms_level") == "2" && level_of.count(parent) && level_of[parent] == "1";
    if (from_survey)
    {
      ++survey_decided;
      TEST_EQUAL(cell(tsv, row, "decision_survey_id"), parent)
    }
    else
    {
      TEST_EQUAL(cell(tsv, row, "decision_survey_id"), "")
    }
  }
  TEST_EQUAL(survey_decided > 0, true) // non-vacuity: the drive selected something
}
END_SECTION

// ---- A stale survey re-plans its CV's waiting work ------------------------------------------------------------
//
// driveStaleSurvey (FLASHIda_TestHelpers.h): survey A plans, k of its commands are dequeued (sent, so committed),
// survey B reaches the engine while the rest still wait, and the drain then empties the queue. The drain hands
// out every command exactly once, so "dequeued" below is "has a scan_commands row".

// Nothing of A was sent, so all of A's plan is uncommitted: B drops it whole and plans afresh, and everything the
// drain hands out is B's decision.
START_SECTION(newer_survey_replaces_unsent_species)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt");
  ABORT_IF(ms1.empty())
  const std::string log_dir = freshLogDir("replan_newer_survey");
  const std::string cfg = topn5Config(log_dir);
  StaleDrive d;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    d = driveStaleSurvey(ida, ms1[0], at(ms1[0], ms1[0].rt + 2.0), 0);
  });
  const auto decided_by = decisionSurveys(log_dir);
  TEST_EQUAL(d.a_plan.size(), 5) // S1..S5
  TEST_EQUAL(replanOf(out, d.b_id).dropped, 5)

  const auto sent = ms2Of(d.after_b);
  TEST_EQUAL(sent.size(), 5)
  for (const auto& c : sent)
  {
    TEST_EQUAL(std::string(c.parent_scan_id), ScanCommandQueue::encode(d.b_id))
    TEST_EQUAL(decidedBy(decided_by, c.scan_id), ScanCommandQueue::encode(d.b_id))
  }
  for (const auto& a : d.a_plan) TEST_EQUAL(hasId(d.after_b, a.scan_id), false)
}
END_SECTION

// S1 was sent before B arrived, so it committed: B is never handed it back and never plans it again. S2..S5 were
// not sent, so B drops them and plans them again as its own.
START_SECTION(sent_species_is_never_dropped_and_not_reselected)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt");
  ABORT_IF(ms1.empty())
  const std::string cfg = topn5Config();
  StaleDrive d;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    d = driveStaleSurvey(ida, ms1[0], at(ms1[0], ms1[0].rt + 2.0), 1);
  });
  ABORT_IF(d.a_plan.size() != 5 || d.before_b.size() != 1)
  const ScanCommand& s1 = d.before_b[0];
  TEST_EQUAL(s1.scan_id, d.a_plan[0].scan_id)
  TEST_EQUAL(replanOf(out, d.b_id).dropped, 4)

  std::vector<ScanCommand> all = ms2Of(d.before_b);
  for (const auto& c : ms2Of(d.after_b)) all.push_back(c);
  TEST_EQUAL(std::count_if(all.begin(), all.end(), [&](const ScanCommand& c) { return nominalOf(c) == nominalOf(s1); }), 1)
  for (size_t i = 1; i < d.a_plan.size(); ++i)
  {
    const ScanCommand& a = d.a_plan[i];
    TEST_EQUAL(hasId(d.after_b, a.scan_id), false) // dropped ...
    TEST_EQUAL(std::any_of(d.after_b.begin(), d.after_b.end(),
                           [&](const ScanCommand& c) {
                             return c.msn_level == 2 && nominalOf(c) == nominalOf(a)
                                    && std::string(c.parent_scan_id) == ScanCommandQueue::encode(d.b_id);
                           }),
               true) // ... and re-added under B
  }
}
END_SECTION

// A quantification group freezes with its first Q: once one went out, the species is committed, its other Qs are
// sent as A decided them (A's parent, A's decision), its group stays open for them, and B opens no group for it.
START_SECTION(quant_group_frozen_after_first_q)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  ABORT_IF(ms1.size() < 9)
  const std::string log_dir = freshLogDir("replan_quant_frozen");
  const std::string cfg = quantSeparateConfig(log_dir);
  const ScanData& a = ms1[8];
  StaleDrive d;
  size_t open_groups = 0;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    d = driveStaleSurvey(ida, a, at(a, a.rt + 2.0), 1);
    open_groups = FLASHIdaTestAccess::quantOpenGroupCount(ida);
  });
  const auto decided_by = decisionSurveys(log_dir);
  ABORT_IF(d.a_plan.size() < 2 || d.before_b.size() != 1)
  for (const auto& c : d.a_plan) TEST_EQUAL(c.scan_description[3], 'Q')
  const int x = nominalOf(d.before_b[0]);

  for (size_t i = 1; i < d.a_plan.size(); ++i)
  {
    const ScanCommand& q = d.a_plan[i];
    TEST_EQUAL(hasId(d.after_b, q.scan_id), true)
    TEST_EQUAL(std::string(q.parent_scan_id), ScanCommandQueue::encode(d.a_id))
    TEST_EQUAL(decidedBy(decided_by, q.scan_id), ScanCommandQueue::encode(d.a_id))
  }
  std::set<int> new_species;
  for (const auto& c : ms2Of(d.after_b))
  {
    if (std::string(c.parent_scan_id) != ScanCommandQueue::encode(d.b_id)) continue;
    TEST_EQUAL(std::abs(nominalOf(c) - x) <= 1, false) // no new plan for X under B
    if (c.scan_description[3] == 'Q') new_species.insert(nominalOf(c));
  }
  TEST_EQUAL(open_groups, 1 + new_species.size()) // X's group is still open, and B opened none for it
  const ReplanLine r = replanOf(out, d.b_id);
  TEST_EQUAL(r.kept, 1)
  TEST_EQUAL(r.dropped, 0)
}
END_SECTION

// A `separate` species whose first charge went out keeps a remainder, and B re-derives it: B ranks the species'
// charges by charge_selection, drops the one already sent, and keeps min(max_charge_states, |admissible on B|) - 1
// of the rest -- now B's decision, queued afresh.
START_SECTION(separate_r_remainder_is_rederived)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  ABORT_IF(ms1.size() < kSeparateB + 1)
  const std::string log_dir = freshLogDir("replan_separate_rederived");
  const std::string cfg = cytcSeparateConfig(log_dir);
  const ScanData& a = ms1[8];
  StaleDrive d;
  std::vector<int> ranked_b;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    d = driveStaleSurvey(ida, a, at(ms1[kSeparateB], a.rt + 2.0), 1);
    if (!d.before_b.empty()) ranked_b = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominalOf(d.before_b[0]), {});
  });
  const auto decided_by = decisionSurveys(log_dir);
  ABORT_IF(d.before_b.size() != 1)
  const int x = nominalOf(d.before_b[0]);
  const int sent_z = std::abs(d.before_b[0].stages[0].charge_state);
  // The inclusion row admits two species on this survey (12351.3 and 12367.3): X is the one whose first charge went.
  const auto x_plan = ofSpecies(d.a_plan, x);
  ABORT_IF(x_plan.size() != 3)
  TEST_EQUAL(ranked_b.empty(), false) // B resolves X, so the remainder is B's to re-derive
  ABORT_IF(ranked_b.empty())

  std::vector<int> expected;
  for (int z : ranked_b)
    if (z != sent_z) expected.push_back(z);
  const int count = std::min<int>(3, static_cast<int>(ranked_b.size())) - 1;
  if (static_cast<int>(expected.size()) > count) expected.resize(std::max(0, count));
  std::sort(expected.begin(), expected.end());

  std::vector<int> a_rest;
  for (size_t i = 1; i < x_plan.size(); ++i) a_rest.push_back(std::abs(x_plan[i].stages[0].charge_state));
  std::sort(a_rest.begin(), a_rest.end());
  TEST_EQUAL(expected.empty(), false)
  TEST_EQUAL(expected != a_rest, true) // non-vacuity: B chooses another remainder than A did

  std::vector<int> got;
  for (const auto& c : ms2Of(d.after_b))
  {
    if (std::abs(nominalOf(c) - x) > 1) continue;
    got.push_back(std::abs(c.stages[0].charge_state));
    TEST_EQUAL(std::string(c.parent_scan_id), ScanCommandQueue::encode(d.b_id))
    TEST_EQUAL(decidedBy(decided_by, c.scan_id), ScanCommandQueue::encode(d.b_id))
    TEST_EQUAL(c.enqueue_timestamp_ms >= d.b_start_ms, true)
  }
  std::sort(got.begin(), got.end());
  TEST_EQUAL(got == expected, true)
  for (size_t i = 1; i < x_plan.size(); ++i) TEST_EQUAL(hasId(d.after_b, x_plan[i].scan_id), false) // the old remainder went
}
END_SECTION

// A survey that does not resolve the species at all (ms1_standard's first scan deconvolves to nothing) keeps the
// remainder as A decided it: same commands, A's parent, A's decision.
START_SECTION(separate_r_remainder_kept_when_unresolved)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto standard = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_standard.txt");
  ABORT_IF(ms1.size() < 9 || standard.empty())
  const std::string log_dir = freshLogDir("replan_separate_kept");
  const std::string cfg = cytcSeparateConfig(log_dir);
  const ScanData& a = ms1[8];
  StaleDrive d;
  std::vector<int> ranked_b;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    d = driveStaleSurvey(ida, a, at(standard[0], a.rt + 2.0), 1);
    if (!d.before_b.empty()) ranked_b = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominalOf(d.before_b[0]), {});
  });
  const auto decided_by = decisionSurveys(log_dir);
  ABORT_IF(d.before_b.size() != 1)
  const auto x_plan = ofSpecies(d.a_plan, nominalOf(d.before_b[0])); // X: the species whose first charge went
  ABORT_IF(x_plan.size() != 3)
  TEST_EQUAL(ranked_b.empty(), true) // non-vacuity: B does not resolve X
  for (size_t i = 1; i < x_plan.size(); ++i)
  {
    const ScanCommand& rest = x_plan[i];
    TEST_EQUAL(hasId(d.after_b, rest.scan_id), true)
    TEST_EQUAL(std::string(rest.parent_scan_id), ScanCommandQueue::encode(d.a_id))
    TEST_EQUAL(decidedBy(decided_by, rest.scan_id), ScanCommandQueue::encode(d.a_id))
  }
  // B drops only the second species (never sent, unresolved on B); nothing of X's is among the dropped commands
  size_t gone = 0;
  for (const auto& c : d.a_plan)
  {
    if (hasId(d.before_b, c.scan_id) || hasId(d.after_b, c.scan_id)) continue;
    ++gone;
    TEST_EQUAL(std::abs(nominalOf(c) - nominalOf(d.before_b[0])) > 1, true)
  }
  const ReplanLine r = replanOf(out, d.b_id);
  TEST_EQUAL(r.kept, 1)
  TEST_EQUAL(r.dropped, static_cast<int>(d.a_plan.size() - x_plan.size()))
  TEST_EQUAL(gone, d.a_plan.size() - x_plan.size())
}
END_SECTION

// A sweep none of whose variants went out is uncommitted: B drops it, the exploration group goes with it, and
// only B's own sweep is active.
START_SECTION(unsent_sweep_dropped_and_group_forgotten)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  ABORT_IF(ms1.size() < 9)
  const std::string cfg = explorationTwinConfig();
  const ScanData& a = ms1[8];
  StaleDrive d;
  int g = 0, last = 0, active = -1;
  bool flag = false;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    g = FLASHIdaTestAccess::explorationLastGroupId(ida) + 1; // the id A's sweep is minted with
    d = driveStaleSurvey(ida, a, at(a, a.rt + 2.0), 0);
    last = FLASHIdaTestAccess::explorationLastGroupId(ida);
    active = FLASHIdaTestAccess::explorationActiveGroupCount(ida);
    flag = FLASHIdaTestAccess::explorationActive(ida);
  });
  ABORT_IF(d.a_plan.size() < 2)
  for (const auto& c : d.a_plan) TEST_EQUAL(c.scan_description[3], 'E')
  TEST_EQUAL(out.find("[EXPL-DROP] group=" + std::to_string(g) + "\n") != std::string::npos, true)
  TEST_EQUAL(last, g + 1) // B's sweep
  TEST_EQUAL(active, 1)   // ... and only B's (its variants were never answered)
  TEST_EQUAL(flag, true)
  for (const auto& c : d.a_plan) TEST_EQUAL(hasId(d.after_b, c.scan_id), false)
  for (const auto& c : ms2Of(d.after_b)) TEST_EQUAL(std::string(c.parent_scan_id), ScanCommandQueue::encode(d.b_id))
}
END_SECTION

// One variant out commits the sweep: it is never dropped, every variant goes out and is answered, the group
// completes, and B plans no sweep of its own for the species.
START_SECTION(started_sweep_runs_to_completion)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto ms2 = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1.size() < 9 || ms2.empty())
  const std::string cfg = explorationTwinConfig();
  const ScanData& a = ms1[8];
  StaleDrive d;
  int g = 0;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    g = FLASHIdaTestAccess::explorationLastGroupId(ida) + 1;
    d = driveStaleSurvey(ida, a, at(a, a.rt + 2.0), 1, &ms2[0]);
  });
  ABORT_IF(d.a_plan.size() < 2 || d.before_b.size() != 1)
  const int x = nominalOf(d.before_b[0]);
  for (const auto& c : d.a_plan) TEST_EQUAL(hasId(d.before_b, c.scan_id) || hasId(d.after_b, c.scan_id), true)
  TEST_EQUAL(out.find("[EXPL-WINNER] group=" + std::to_string(g) + " ") != std::string::npos, true)
  TEST_EQUAL(out.find("[EXPL-DROP] group=" + std::to_string(g) + "\n") == std::string::npos, true)
  for (const auto& c : ms2Of(d.after_b))
  {
    if (std::string(c.parent_scan_id) != ScanCommandQueue::encode(d.b_id)) continue;
    TEST_EQUAL(std::abs(nominalOf(c) - x) <= 1, false)
  }
}
END_SECTION

// FAIMS cycling: a survey at another CV leaves A's CV alone -- A's plan goes out whole, as A decided it.
START_SECTION(other_cv_survey_does_not_replan)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt");
  ABORT_IF(ms1.empty())
  const std::string log_dir = freshLogDir("replan_other_cv");
  const std::string cfg = faims3cvConfig(log_dir);
  StaleDrive d;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    d = driveStaleSurvey(ida, ms1[0], at(ms1[0], ms1[0].rt + 2.0), 0, nullptr, -60.0);
  });
  const auto decided_by = decisionSurveys(log_dir);
  ABORT_IF(d.a_plan.empty())
  for (const auto& c : d.a_plan)
  {
    TEST_REAL_SIMILAR(c.faims_cv, -40.0)
    TEST_EQUAL(hasId(d.after_b, c.scan_id), true)
    TEST_EQUAL(std::string(c.parent_scan_id), ScanCommandQueue::encode(d.a_id))
    TEST_EQUAL(decidedBy(decided_by, c.scan_id), ScanCommandQueue::encode(d.a_id))
  }
  const ReplanLine r = replanOf(out, d.b_id);
  TEST_EQUAL(r.dropped, 0)
  TEST_EQUAL(r.kept, 0)
  TEST_EQUAL(r.added > 0, true) // B planned at its own CV
}
END_SECTION

// The race: the drain commits while B is being planned. A drain thread starts at a different point of B's
// processing in every round and keeps dequeuing, and the seam commits one more command between B's plan and its
// swap -- the one gap a fold cannot see. Whatever the drain took, no species goes out twice, every unit B added
// goes out, and nothing is left registered.
START_SECTION(race_commit_between_plan_and_swap)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt");
  ABORT_IF(ms1.empty())
  // tqscore_threshold 0 (owner 2026-10-23): every admitted species then carries a hard bar, so a species can go out
  // twice only through the race -- at 0.9 the qscore bar re-admits a sub-0.9 species at equal score, by design.
  const std::string cfg = with(topn5Config(), R"("tqscore_threshold": 0.9)", R"("tqscore_threshold": 0.0)");
  // One survey's processing time, to spread the drain thread's start across B's.
  std::chrono::microseconds t_survey{0};
  captureCout([&] {
    FLASHIda scratch(const_cast<char*>(cfg.c_str()));
    ScanCommand c{};
    scratch.getNextScanCommand(c);
    const auto t0 = std::chrono::steady_clock::now();
    scratch.processScan(ms1[0].mzs.data(), ms1[0].ints.data(), (int)ms1[0].mzs.size(), ms1[0].rt, 1, c.scan_description, 0.0, 0);
    t_survey = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - t0);
  });

  int rounds = 0, twice = 0, lost = 0, leaked = 0, invariant = 0, seam_misfires = 0, committed_in_gap = 0;
  for (int i = 0; i < 200; ++i)
  {
    const auto delay = t_survey * (i % 12) / 10;
    const auto period = std::max(std::chrono::microseconds(50), t_survey / 40);
    const RaceRun r = race(cfg, ms1[0], delay, period, true);
    ++rounds;
    std::vector<ScanCommand> all = r.seam;
    all.insert(all.end(), r.drained.begin(), r.drained.end());
    std::map<int, int> by_id, by_nominal;
    int b_sent = 0;
    for (const auto& c : ms2Of(all))
    {
      ++by_id[c.scan_id];
      ++by_nominal[nominalOf(c)];
      if (std::string(c.parent_scan_id) == ScanCommandQueue::encode(r.b_id)) ++b_sent;
    }
    bool dup = false;
    for (const auto& [k, n] : by_id) dup = dup || n > 1;
    for (const auto& [k, n] : by_nominal) dup = dup || n > 1;
    if (dup) ++twice;
    // b_sent == added holds only because topn5 is `single`: one command per unit. Under `separate` or a sweep a unit
    // sends several commands and the two would differ without anything being lost.
    if (b_sent != replanOf(r.out, r.b_id).added) ++lost;
    if (r.units_left != 0 || r.queued_left != 0) ++leaked;
    if (r.out.find("[REPLAN-INVARIANT]") != std::string::npos) ++invariant;
    if (r.seam.size() != 1) ++seam_misfires;
    if (r.seam.size() == 1 && r.seam[0].msn_level == 2 && std::string(r.seam[0].parent_scan_id) == ScanCommandQueue::encode(r.a_id))
      ++committed_in_gap;
  }
  TEST_EQUAL(rounds, 200)
  TEST_EQUAL(twice, 0)
  TEST_EQUAL(lost, 0)
  TEST_EQUAL(leaked, 0)
  TEST_EQUAL(invariant, 0)
  TEST_EQUAL(seam_misfires, 0)
  TEST_EQUAL(committed_in_gap, rounds) // non-vacuity: every round committed one of A's commands in the gap
}
END_SECTION

// The same race under `separate`, where the seam commits X's first charge: B's fresh unit for X is reconciled away,
// A's remainder is kept, and no charge of X goes out twice or is lost.
START_SECTION(race_separate_commit_between_plan_and_swap)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  ABORT_IF(ms1.size() < 9)
  const std::string cfg = cytcSeparateConfig(); // tqscore_threshold 0 already
  const ScanData& a = ms1[8];
  std::chrono::microseconds t_survey{0};
  captureCout([&] {
    FLASHIda scratch(const_cast<char*>(cfg.c_str()));
    ScanCommand c{};
    scratch.getNextScanCommand(c);
    const auto t0 = std::chrono::steady_clock::now();
    scratch.processScan(a.mzs.data(), a.ints.data(), (int)a.mzs.size(), a.rt, 1, c.scan_description, 0.0, 0);
    t_survey = std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::steady_clock::now() - t0);
  });

  int rounds = 0, twice = 0, lost = 0, leaked = 0, invariant = 0, seam_misfires = 0, x_in_gap = 0;
  for (int i = 0; i < 24; ++i)
  {
    const auto delay = t_survey * (i % 12) / 10;
    const auto period = std::max(std::chrono::microseconds(50), t_survey / 40);
    // The drain thread waits for the seam, so the seam takes the queue's head: X's first charge.
    const RaceRun r = race(cfg, a, delay, period, false);
    ++rounds;
    if (r.a_plan.empty())
    {
      ++lost;
      continue;
    }
    // The inclusion row admits two species on this survey (12351.3 and 12367.3); X is the one queued first.
    const int x = nominalOf(r.a_plan[0]);
    const auto x_plan = ofSpecies(r.a_plan, x);
    std::vector<ScanCommand> all = r.seam;
    all.insert(all.end(), r.drained.begin(), r.drained.end());
    std::map<int, int> by_charge;
    std::map<int, int> by_id;
    for (const auto& c : ms2Of(all))
    {
      ++by_id[c.scan_id];
      if (std::abs(nominalOf(c) - x) <= 1) ++by_charge[std::abs(c.stages[0].charge_state)];
    }
    bool dup = false;
    for (const auto& [k, n] : by_id) dup = dup || n > 1;
    for (const auto& [k, n] : by_charge) dup = dup || n > 1;
    if (dup) ++twice;
    if (by_charge.size() != x_plan.size()) ++lost;
    if (r.units_left != 0 || r.queued_left != 0) ++leaked;
    if (r.out.find("[REPLAN-INVARIANT]") != std::string::npos) ++invariant;
    if (r.seam.size() != 1) ++seam_misfires;
    if (r.seam.size() == 1 && r.seam[0].scan_id == x_plan[0].scan_id) ++x_in_gap;
  }
  TEST_EQUAL(rounds, 24)
  TEST_EQUAL(twice, 0)
  TEST_EQUAL(lost, 0)
  TEST_EQUAL(leaked, 0)
  TEST_EQUAL(invariant, 0)
  TEST_EQUAL(seam_misfires, 0)
  TEST_EQUAL(x_in_gap, rounds) // non-vacuity: every round committed X's first charge in the gap
}
END_SECTION

// A follow-up is derived from a command that went out, so its parent is never a staged (still-registered)
// command: the queue prints [REPLAN-INVARIANT] if one ever is, and every 'C' row's parent has its own row.
START_SECTION(follow_up_parent_is_never_staged)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_standard.txt");
  const auto ms2 = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt");
  ABORT_IF(ms1.empty() || ms2.empty())
  // ms1_standard selects at most one species per scan, and most scans none: take the first that selects.
  int pick = -1;
  captureCout([&] {
    const std::string cfg = tagTargetingConfig("");
    for (size_t i = 0; i < ms1.size() && pick < 0; ++i)
    {
      FLASHIda scratch(const_cast<char*>(cfg.c_str()));
      ScanCommand c{};
      scratch.getNextScanCommand(c);
      if (scratch.processScan(ms1[i].mzs.data(), ms1[i].ints.data(), (int)ms1[i].mzs.size(), ms1[i].rt, 1, c.scan_description, 0.0, 0) > 0)
        pick = static_cast<int>(i);
    }
  });
  ABORT_IF(pick < 0)

  int conditional = 0;
  for (int k = 0; k <= 1; ++k)
  {
    const std::string log_dir = freshLogDir("replan_follow_up_k" + std::to_string(k));
    const std::string cfg = tagTargetingConfig(log_dir);
    const ScanData& s = ms1[pick];
    const std::string out = captureCout([&] {
      FLASHIda ida(const_cast<char*>(cfg.c_str()));
      driveStaleSurvey(ida, s, at(s, s.rt + 2.0), k, &ms2[0]);
    });
    TEST_EQUAL(out.find("[REPLAN-INVARIANT]"), std::string::npos)
    const TSVFile tsv = TSVFile::parse(log_dir + "/scan_commands.tsv");
    std::set<std::string> ids;
    for (const auto& row : tsv.rows) ids.insert(cell(tsv, row, "tracking_id"));
    for (const auto& row : tsv.rows)
    {
      if (cell(tsv, row, "scan_type") != "conditional") continue;
      ++conditional;
      TEST_EQUAL(ids.count(cell(tsv, row, "parent_tracking_id")), 1)
    }
  }
  TEST_EQUAL(conditional > 0, true) // non-vacuity: a follow-up was raised
}
END_SECTION

// Every place FLASHIda.cpp puts a command in the queue, by enclosing function and call text. Re-planning owns the
// survey's MS2s through replacePlan alone; every other site is a permanent push of a command derived from one that
// went out, or an MS1. A new or moved site fails here until this table is updated -- and whoever updates it decides
// whether the new site must go through the plan.
START_SECTION(push_site_inventory)
{
  std::ifstream in("../src/openms/source/ANALYSIS/TOPDOWN/FLASHIda.cpp");
  ABORT_IF(!in.good())
  std::stringstream ss;
  ss << in.rdbuf();
  const std::multiset<std::string> expected = {
    // the survey's plan, in one call
    "runSurveyMS1_ | const ReplanResult replan = queue_.replacePlan(faims_cv, std::move(units), adjust, replace_remainder_of);",
    // the six permanent sites
    "refineRoundClose_ | stampAndPush_(c, precursor_id);",                  // refine re-dispatch
    "processScan | stampAndPush_(c, precursor_id);",                        // exploration MS2 return
    "processScan | stampAndPush_(ident, cons.winner_precursor_id);",        // bought R
    "processScan | stampAndPush_(cond, precursor_id);",                     // tag follow-up
    "processScan | stampAndPush_(c, precursor_id);",                        // MS3 targeting
    "processScan | stampAndPush_(c, precursor_id);",                        // exploration MS3 return
    // the bought R, registered adjustable at the quant verdict
    "processScan | queue_.registerAdjustable(ident.scan_id, cons.q_covered_charges, cons.winner_survey_id,",
    // the MS1 family
    "runSurveyMS1_ | queue_.push(ms1);",                                    // FAIMS CV transition
    "getNextScanCommand | queue_.push(ms1_cmd);",                           // cycle-time MS1
    "getNextScanCommand | queue_.push(mon_cmd);",                           // monitor scan
    "getNextScanCommand | dequeued = queue_.pushAndDequeue(ms1_cmd);",     // idle survey
  };
  const std::multiset<std::string> found = pushSites(ss.str());
  for (const auto& s : found)
    if (expected.count(s) != found.count(s)) std::cerr << "push site: " << s << " x" << found.count(s) << std::endl;
  for (const auto& s : expected)
    if (expected.count(s) != found.count(s)) std::cerr << "expected push site: " << s << " x" << expected.count(s) << std::endl;
  TEST_EQUAL(found == expected, true)
}
END_SECTION

// ---- The bought R is adjustable, never droppable --------------------------------------------------------------
//
// ms1_cytc.txt[8] fed as the first survey: a `separate` quantification group of cytC (three Qs, decided on the
// two real iTMT spectra) buys an identification R at its purest window, z13 -- while charge_selection ("snr")
// ranks z14 first on that survey. The R is caught WAITING, a newer survey reaches the engine ahead of it, and the
// drain then takes it. Every expected charge is ranked by rankedCharges from the PeakGroups, not by the engine.
// The R may move only among the charges its group's Q scans covered (owner 2026-10-23): the anchors and notches of
// the Qs that returned, registered with the R.

// A newer survey that ranks another Q's charge first moves the R there: the R row's charge is the newer survey's
// choice, that survey decided it, its parent is still the winning Q, and its queue time restarts.
START_SECTION(bought_r_adjusted_to_better_charge)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("separate");
  const ScanData& a = ms1[8];

  WaitingR w;
  int b_id = 0, expected = 0;
  ScanCommand after{}, sent{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, a, q_odd[0], q_even[0]);
    if (!w.bought) return;
    std::this_thread::sleep_for(std::chrono::milliseconds(5)); // a queue-time reset is then a strictly later stamp
    b_id = feedSurvey(ida, a, a.rt + 2.0);                      // the same spectrum, two seconds on
    const auto ranked = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), SpectralDeconvolution::getNominalMass(w.r.mono_mass),
                                      w.covered);
    if (!ranked.empty()) expected = ranked.front();
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
    ida.getNextScanCommand(sent);
  });
  ABORT_IF(!w.bought)
  const int bought_z = std::abs(w.r.stages[0].charge_state);   // the winning Q's charge (the verdict's purest window)
  TEST_EQUAL(std::find(w.covered.begin(), w.covered.end(), bought_z) != w.covered.end(), true)
  TEST_EQUAL(w.covered == w.q_charges, true) // registered with exactly the charges the group's Qs isolated
  TEST_EQUAL(w.covered.size() > 1, true)
  TEST_EQUAL(expected, 14) // non-vacuity: the newer survey prefers another charge, and a Q covered it
  TEST_NOT_EQUAL(bought_z, expected) // so the adjust below is a MOVE, not a re-confirm
  TEST_EQUAL(std::find(w.covered.begin(), w.covered.end(), expected) != w.covered.end(), true)

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(std::abs(after.stages[0].charge_state), expected)
  const std::string desc(after.scan_description);
  TEST_EQUAL(desc.substr(desc.rfind('@') + 1), std::to_string(expected))
  TEST_EQUAL(std::string(after.parent_scan_id), std::string(w.r.parent_scan_id)) // lineage: still the winning Q
  TEST_EQUAL(after.enqueue_timestamp_ms > w.r.enqueue_timestamp_ms, true)
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, b_id)
  TEST_EQUAL(decision->isolation_purity >= 0.0 && decision->isolation_purity <= 1.0, true)
  const auto lines = replanLines(out);
  ABORT_IF(lines.size() != 2)
  TEST_EQUAL(lines.back().survey, ScanCommandQueue::encode(b_id))
  TEST_EQUAL(lines.back().adjusted, 1)

  // ...and that is the command the drain hands out.
  TEST_EQUAL(sent.scan_id, w.r.scan_id)
  TEST_EQUAL(std::abs(sent.stages[0].charge_state), expected)
}
END_SECTION

// A newer survey that ranks the R's own charge first re-confirms it (owner 2026-10-22): the charge stays, but the
// newer survey is now the one that decided it, and the queue time restarts. The newer survey is A without whichever
// charges outrank the R's on it, so the R's charge leads among those the Qs covered.
START_SECTION(bought_r_reconfirmed_same_charge_resets_decision)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("separate");
  const ScanData& a = ms1[8];

  WaitingR w;
  int b_id = 0, expected = 0;
  std::vector<int> ranked_on_a;
  ScanCommand after{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, a, q_odd[0], q_even[0]);
    if (!w.bought) return;
    const int nominal = SpectralDeconvolution::getNominalMass(w.r.mono_mass);
    const int own = std::abs(w.r.stages[0].charge_state);
    ranked_on_a = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, w.covered);
    // Take out whatever charge outranks the R's, one at a time, until the R's leads: removing an envelope moves the
    // others' SNR, so the charges to remove are found on a scratch engine's deconvolution of each candidate.
    std::vector<int> removed;
    ScanData b = a;
    for (int i = 0; i < 6; ++i)
    {
      std::vector<int> ranked;
      {
        FLASHIda scratch(const_cast<char*>(cfg.c_str()));
        feedSurvey(scratch, b, a.rt + 2.0);
        ranked = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(scratch), nominal, w.covered);
      }
      if (ranked.empty() || ranked.front() == own) break;
      removed.push_back(ranked.front());
      b = withoutCharges(a, w.r.mono_mass, removed);
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    b_id = feedSurvey(ida, b, a.rt + 2.0);
    const auto ranked = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, w.covered);
    if (!ranked.empty()) expected = ranked.front();
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
  });
  ABORT_IF(!w.bought || ranked_on_a.empty())
  const int bought_z = std::abs(w.r.stages[0].charge_state);
  TEST_EQUAL(ranked_on_a.front() != bought_z, true) // non-vacuity: on the Q's survey another charge led
  TEST_EQUAL(expected, bought_z)                    // ...on the newer survey the R's own leads

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(std::abs(after.stages[0].charge_state), bought_z)
  TEST_EQUAL(std::string(after.scan_description), std::string(w.r.scan_description))
  TEST_EQUAL(after.enqueue_timestamp_ms > w.r.enqueue_timestamp_ms, true)
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, b_id)
  const auto lines = replanLines(out);
  ABORT_IF(lines.empty())
  TEST_EQUAL(lines.back().adjusted, 1)
}
END_SECTION

// A newer survey that resolves the species only at charges no Q covered cannot move the R: the newer survey is A
// without every covered charge, so it still resolves cytC, at charges the group never measured. The R goes out as
// decided, at its charge, still attributed to the Q's survey.
START_SECTION(bought_r_keeps_charge_when_no_covered_charge_resolves)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("separate");
  const ScanData& a = ms1[8];

  WaitingR w;
  std::vector<int> ranked_on_b, ranked_on_both;
  ScanCommand after{}, sent{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, a, q_odd[0], q_even[0]);
    if (!w.bought) return;
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    feedSurvey(ida, withoutCharges(a, w.r.mono_mass, w.covered), a.rt + 2.0);
    const int nominal = SpectralDeconvolution::getNominalMass(w.r.mono_mass);
    ranked_on_b = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, {});
    ranked_on_both = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, w.covered);
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
    ida.getNextScanCommand(sent);
  });
  ABORT_IF(!w.bought || ranked_on_b.empty())
  // Non-vacuity: the newer survey resolves cytC, at no charge a Q covered.
  TEST_EQUAL(std::find(w.covered.begin(), w.covered.end(), ranked_on_b.front()) == w.covered.end(), true)
  TEST_EQUAL(ranked_on_both.empty(), true)

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(std::string(after.scan_description), std::string(w.r.scan_description))
  TEST_EQUAL(after.stages[0].charge_state, w.r.stages[0].charge_state)
  TEST_REAL_SIMILAR(after.stages[0].precursor_mz, w.r.stages[0].precursor_mz)
  TEST_EQUAL(after.enqueue_timestamp_ms, w.r.enqueue_timestamp_ms) // nothing re-chose it
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, w.q_survey)
  const auto lines = replanLines(out);
  ABORT_IF(lines.empty())
  TEST_EQUAL(lines.back().adjusted, 0)
  TEST_EQUAL(sent.scan_id, w.r.scan_id)
  TEST_EQUAL(sent.stages[0].charge_state, w.r.stages[0].charge_state)
}
END_SECTION

// A multiplexed bought R re-chooses its whole SET (owner 2026-10-22) within the charges its Q covered -- the Q's
// anchor and notches (owner 2026-10-23): anchor and notches become those charges in the newer survey's order, and
// its scan_description names the new anchor. The newer survey is A without the last charge of the R's set: it
// still resolves cytC at charges the Q never isolated, and none of them may fill the freed slot.
START_SECTION(bought_r_multiplexed_set_rechosen)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("multiplexed");
  const ScanData& a = ms1[8];

  WaitingR w;
  int b_id = 0;
  std::vector<int> expected, ranked_on_b;
  ScanCommand after{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, a, q_odd[0], q_even[0]);
    if (!w.bought) return;
    b_id = feedSurvey(ida, withoutCharges(a, w.r.mono_mass, {isolatedCharges(w.r).back()}), a.rt + 2.0);
    const int nominal = SpectralDeconvolution::getNominalMass(w.r.mono_mass);
    ranked_on_b = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, {});
    expected = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, w.covered);
    if (expected.size() > 3) expected.resize(3);
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
  });
  ABORT_IF(!w.bought)
  TEST_EQUAL(isolatedCharges(w.r).size(), 3)
  std::vector<int> isolated_sorted = isolatedCharges(w.r);
  std::sort(isolated_sorted.begin(), isolated_sorted.end());
  TEST_EQUAL(w.covered == isolated_sorted, true) // the one Q's anchor and notches, which the R inherited
  TEST_EQUAL(expected.size(), 2)                 // the freed slot stays empty ...
  TEST_EQUAL(ranked_on_b.size() > expected.size(), true) // ... though the newer survey resolves more charges
  TEST_EQUAL(expected != isolatedCharges(w.r), true)     // non-vacuity: the newer survey's set differs

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(isolatedCharges(after) == expected, true)
  const std::string desc(after.scan_description);
  TEST_EQUAL(desc.substr(desc.rfind('@') + 1), std::to_string(expected.front()))
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, b_id)
  const auto lines = replanLines(out);
  ABORT_IF(lines.empty())
  TEST_EQUAL(lines.back().adjusted, 1)
}
END_SECTION

// The covered set is a FILTER, not a preference: a newer survey that ranks a charge no Q covered first and a
// covered one after it moves the R to the covered one. The Q's survey is A without A's own top charge, so that
// charge is not covered; the newer survey is A without the R's charge, so the top charge leads again and the R
// cannot stay where it is.
START_SECTION(bought_r_moves_only_among_covered_charges)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("separate");
  const ScanData& a = ms1[8];
  const double cytc = 12351.3; // the species ms1_cytc.txt[8]'s verdict buys an R for (the sections above)

  int top_on_a = 0;
  {
    FLASHIda scratch(const_cast<char*>(cfg.c_str()));
    captureCout([&] { feedSurvey(scratch, a, a.rt); });
    const auto ranked = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(scratch), SpectralDeconvolution::getNominalMass(cytc), {});
    if (!ranked.empty()) top_on_a = ranked.front();
  }
  ABORT_IF(top_on_a == 0)

  WaitingR w;
  int b_id = 0, expected = 0;
  std::vector<int> ranked_on_b;
  ScanCommand after{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, withoutCharges(a, cytc, {top_on_a}), q_odd[0], q_even[0]);
    if (!w.bought) return;
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    b_id = feedSurvey(ida, withoutCharges(a, w.r.mono_mass, {std::abs(w.r.stages[0].charge_state)}), a.rt + 2.0);
    const int nominal = SpectralDeconvolution::getNominalMass(w.r.mono_mass);
    ranked_on_b = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, {});
    const auto ranked = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, w.covered);
    if (!ranked.empty()) expected = ranked.front();
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
  });
  ABORT_IF(!w.bought || ranked_on_b.empty() || expected == 0)
  const int bought_z = std::abs(w.r.stages[0].charge_state);
  // Non-vacuity: the newer survey's best charge is one no Q covered, and a covered charge other than the R's own
  // follows it.
  TEST_EQUAL(std::find(w.covered.begin(), w.covered.end(), top_on_a) == w.covered.end(), true)
  TEST_EQUAL(std::find(w.covered.begin(), w.covered.end(), ranked_on_b.front()) == w.covered.end(), true)
  TEST_EQUAL(expected != bought_z, true)

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(std::abs(after.stages[0].charge_state), expected)
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, b_id)
  const auto lines = replanLines(out);
  ABORT_IF(lines.empty())
  TEST_EQUAL(lines.back().adjusted, 1)
}
END_SECTION

// A multiplexed R the newer survey could re-isolate at fewer than min_charge_states (2) of its covered charges is
// kept as decided (owner 2026-10-23, ADR-0040's floor): its set, its description and its decision survey stay. The
// newer survey is A without every covered charge but one.
START_SECTION(bought_r_multiplexed_below_floor_is_kept)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("multiplexed");
  const ScanData& a = ms1[8];

  WaitingR w;
  std::vector<int> ranked_on_b;
  ScanCommand after{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, a, q_odd[0], q_even[0]);
    if (!w.bought || w.covered.empty()) return;
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    const std::vector<int> others(w.covered.begin() + 1, w.covered.end());
    feedSurvey(ida, withoutCharges(a, w.r.mono_mass, others), a.rt + 2.0);
    ranked_on_b = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), SpectralDeconvolution::getNominalMass(w.r.mono_mass),
                                w.covered);
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
  });
  ABORT_IF(!w.bought)
  TEST_EQUAL(isolatedCharges(w.r).size(), 3)
  TEST_EQUAL(ranked_on_b.size(), 1) // non-vacuity: one covered charge resolves, so there was something to move to

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(isolatedCharges(after) == isolatedCharges(w.r), true)
  TEST_EQUAL(std::string(after.scan_description), std::string(w.r.scan_description))
  TEST_EQUAL(after.enqueue_timestamp_ms, w.r.enqueue_timestamp_ms) // nothing re-chose it
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, w.q_survey)
  const auto lines = replanLines(out);
  ABORT_IF(lines.empty())
  TEST_EQUAL(lines.back().adjusted, 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
