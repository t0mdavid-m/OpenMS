// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Trajectory pin: pooled_identification.tsv is a per-precursor TRAJECTORY, not a single snapshot.
//
//   ProteoformTracker::finalizeMS2(pid)               -> emits ONE MS2 baseline row (trigger = "MS2").
//   ProteoformTracker::foldMs3(pid, ion, scan_id)  -> ADDITIVELY folds the staged MS3 scan into the
//                                                     EXISTING finalized model (no winner re-pick, no
//                                                     fragments reset), bumps update_index, and emits
//                                                     ONE more row tagged trigger = <ion>,
//                                                     trigger_scan_id = <id> (ProteoformTracker.cpp:296-311).
//
// The old code RESET / dropped MS3 evidence (it re-ran finalize, which clears m.fragments at
// ProteoformTracker.cpp:280 and re-tags the row "MS2"). This test pins BOTH the accumulation and the
// trajectory tagging:
//
//   (A) model accumulation  -- mdl->fragments.size() grows 2 -> 3 -> 4 (a reset would DROP it to the
//                              single MS3 fragment, i.e. 1), update_index 1 -> 2 -> 3, coverage
//                              non-decreasing.
//   (B) pooled trajectory   -- exactly 3 data rows; per-row trigger / trigger_scan_id / update_index /
//                              n_fragments columns are pinned.
//
// ADVERSARIAL (C): if foldMs3 were replaced by the old finalize, finalize would CLEAR m.fragments
// (ProteoformTracker.cpp:280) and re-map only the staged MS3 scan -> fragments.size() would drop to 1
// (the single MS3 equiv fragment), failing assertion (A) (size must be > 2). finalize also tags the row
// "MS2" with trigger_scan_id = encode(winner_scan_id), failing assertion (B) (trigger must be "b3"/"y3"
// and trigger_scan_id "AAB"/"AAC"). Either failure catches the regression.
//
// CONSTRUCTION NOTE: the WINNER MS2 scan's fragments are pooled verbatim (winner path) and MS3 folds map
// by equiv index, so neither depends on a peak match. The NON-winner MS2 scan (102) is re-matched against
// the winner ladder (Change 1: winner-anchored pooling), so its raw peak mass MUST equal a real winner b/y
// ion (computed in-test with FragmentAnalysis::computePTMAdjustedFragmentMasses) or it would be DROPPED.
// Hand-built PeakGroups have getChargeIntensity() == 0, so intensities are 0; mapping is by MASS, so that
// is fine. The accumulation under test depends on the fragments map keying by (ion_type, winner index) --
// distinct ions per fold -> new entries.

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms2Params.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3YieldModel.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/SYSTEM/File.h>

#include "FLASHIda_EventCapture.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

using namespace OpenMS;

namespace
{
  const char* WINNER_SEQ = "PEPTIDEK";   // 8-residue winner proteoform (L = 8)

  // Build a minimal synthetic PeakGroup at a given (mz, mono_mass, charge). getChargeIntensity is 0 for
  // such hand-built groups (no scoring pass); the fed fragment still maps (closest-mass fallback).
  PeakGroup makeSyntheticPeakGroup(double mz, double mono_mass, int charge)
  {
    PeakGroup pg(charge, charge, true);
    pg.setMonoisotopicMass(mono_mass);
    FLASHHelperClasses::LogMzPeak lp;
    lp.mz = mz;
    lp.abs_charge = charge;
    lp.intensity = 1000.0f;
    pg.push_back(lp);
    return pg;
  }

  // A single-fragment MS2 ProteoformMatch for the winner proteoform: score >= 0 + non-empty sequence so
  // finalize accepts it as the winner; full-sequence region (region_start/end = -1) -> identity-frame
  // mapping; carries one MS2 fragment (region-relative ion_type/ion_index/observed_mass).
  FragmentAnalysis::ProteoformMatch makeMatch(const std::string& ion_type, int ion_index, double observed_mass)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.score = 1.0;
    m.identification_score = 1.0;  // finalizeMS2 ranks the winner by it (ADR-0061 d.5)
    m.region_start = -1;
    m.region_end = -1;
    m.matched_protein = "synthetic";
    m.proteoform_sequence = WINNER_SEQ;
    FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
    fm.ion_type = ion_type;      // MS2: region-relative ion type/index
    fm.ion_index = ion_index;
    fm.observed_mass = observed_mass;
    m.fragments.push_back(fm);
    return m;
  }

  // A single-fragment MS3 ProteoformMatch. mapScanOntoModel_'s ms_level == 3 branch keys the fragment on
  // equiv_type / equiv_index, pools adjusted_mass, and finds the fragment's own peak in the scan by
  // observed_mass (the measured sub-fragment mass: what the scan's peak list holds). score = -1 and
  // proteoform_sequence = "" because foldMs3 never re-picks a winner: the model already carries the winner
  // sequence, and mapScanOntoModel_ uses the MODEL's sequence (m.proteoform_sequence), not the match's.
  FragmentAnalysis::ProteoformMatch makeMatchMs3(const std::string& equiv_type, int equiv_index, double adjusted_mass,
                                                 double observed_mass)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.score = -1.0;              // not a winner candidate (foldMs3 folds onto the existing baseline)
    m.region_start = -1;         // full-protein equiv index maps straight into the winner region
    m.region_end = -1;
    m.proteoform_sequence = "";  // model already holds WINNER_SEQ; mapScanOntoModel_ uses the model's
    FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
    fm.equiv_type = equiv_type;       // MS3: full-protein equivalent ion type ("b"/"y")
    fm.equiv_index = equiv_index;     // MS3: full-protein equivalent ion index
    fm.adjusted_mass = adjusted_mass; // MS3: the measured mass re-expressed in the proteoform frame
    fm.observed_mass = observed_mass; // MS3: the measured sub-fragment mass, as it stands in the MS3 scan
    m.fragments.push_back(fm);
    return m;
  }

  // A minimal 2-stage MS2 context command. Only used as the feedScan ms2_ctx argument; on the MS3
  // re-feeds it is ignored because has_ms2_ctx is already true (captured at the first MS2 feedScan).
  ScanCommand makeMs2Ctx()
  {
    ScanCommand ctx{};
    ctx.scan_id = 1;
    ctx.msn_level = 2;
    ctx.num_stages = 1;
    ctx.mono_mass = 900.0;
    ctx.stages[0].precursor_mz = 451.0;
    ctx.stages[0].isolation_width = 2.0;
    ctx.stages[0].charge_state = 2;
    ctx.stages[0].collision_energy = 29.0;
    std::strncpy(ctx.stages[0].activation_type, "HCD", sizeof(ctx.stages[0].activation_type) - 1);
    return ctx;
  }

  // Split a TSV line on '\t'. NOTE: std::getline drops a TRAILING empty field, but in this scenario the
  // last column (protein, the "synthetic" accession of every fed match) is always non-empty, so every data
  // row yields all 20 fields; interior
  // empty fields (localized_mods/ambiguous_mods) sit between tabs and ARE preserved as "".
  std::vector<std::string> splitTab(const std::string& line)
  {
    std::vector<std::string> fields;
    std::istringstream ls(line);
    std::string fld;
    while (std::getline(ls, fld, '\t')) fields.push_back(fld);
    return fields;
  }

  // ---- ADR-0066: the per-proteoform plan gate (Ms3PlanGate) ---------------------------------------------------
  //
  // A hand-built PeakGroup has getChargeIntensity() == 0, so every observation below has intensity 0 and a window
  // SNR of 0 (FragmentAnalysis::windowSnr returns 0 for no signal). The sections that need a charge ADMITTED pass no
  // spectrum (no SNR test) or snr_min -1; the SNR section passes a spectrum and moves only snr_min. The 2x rule is
  // driven through the snapshot: a sent intensity of 0 admits a re-acquisition at 0, a sent intensity of 1 refuses it.

  /// The tracker config with characterization's objective and budget as parameters (no log folder: nothing is read).
  /// @p order writes characterization.order (ADR-0073); empty leaves the key out, which is "ranked".
  std::string gateConfig(const std::string& mode, int max_targets, const std::string& fragment_charges = "", const std::string& order = "")
  {
    const std::string charges = fragment_charges.empty() ? std::string() : ", \"fragment_charges\": \"" + fragment_charges + "\"";
    const std::string send_order = order.empty() ? std::string() : ", \"order\": \"" + order + "\"";
    return std::string(R"JSON({
    "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 1, "max_charge": 50, "min_mass": 100, "max_mass": 50000, "tol": [10, 10, 10] },
    "precursor_selection": { "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 3 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
      "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25, "resolution": 120000 }
    },
    "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 } },
    "files": { "target_logs": [], "inclusion_list": "" },
    "identification": { "reference": { "sequence": "PEPTIDEK" }, "flashtnt": { "max_blind_mod_count": 2, "max_mod_mass": 700 },
      "characterization": { "mode": ")JSON") + mode + "\", \"max_targets\": " + std::to_string(max_targets) + charges + send_order + " } } }";
  }

  const double PROTON_MASS = 1.00728;
  const double ISOTOPE_STEP = 1.00235;

  /// A fragment PeakGroup seen at every charge in @p charges with TWO isotopes each, so PeakRecord::by_charge holds a
  /// real (non-zero-width) record per charge and ms2_by_charge carries all of them. The representative charge is the
  /// lowest (getMaxIntensityAbsCharge with no intensity). Each charge carries a charge SNR of 5, above the notch
  /// selection's snr_threshold, so separate and multiplexed see every charge.
  PeakGroup makeEnvelope(double mono_mass, const std::vector<int>& charges)
  {
    const int zmin = *std::min_element(charges.begin(), charges.end());
    const int zmax = *std::max_element(charges.begin(), charges.end());
    PeakGroup pg(zmin, zmax, true);
    pg.setMonoisotopicMass(mono_mass);
    for (int z : charges)
    {
      pg.setChargeSNR(z, 5.0f);
      for (int iso = 0; iso < 2; ++iso)
      {
        FLASHHelperClasses::LogMzPeak lp;
        lp.mz = (mono_mass + iso * ISOTOPE_STEP) / z + PROTON_MASS;
        lp.abs_charge = z;
        lp.isotopeIndex = iso;
        lp.intensity = 1000.0f;
        pg.push_back(lp);
      }
    }
    return pg;
  }

  struct GateFrag
  {
    std::string ion_type;
    int ion_index;
    double mass;
    std::vector<int> charges;   ///< one entry: a single-peak group (no envelope); more: makeEnvelope
    double iso_cos = -1.0;      ///< ADR-0069: >= 0 sets the group's isotope cosine (the yield model's input); < 0 leaves it
  };

  /// Feed ONE winner MS2 scan carrying @p frags (the winner match claims each, pooled verbatim) and finalize. With
  /// @p ambiguous the match carries one modification on residues 4..5, which b5, b6 and y6 contain and nothing brackets.
  /// ADR-0073: every entry of @p sites is a residue range (first..last) that becomes one more modification, each with
  /// its own mass shift.
  void feedGateWinner(ProteoformTracker& tracker, int precursor_id, const std::vector<GateFrag>& frags, bool ambiguous = false,
                      const std::vector<std::pair<int, int>>& sites = {})
  {
    Ms2Params p;
    p.collision_energy = 29.0;
    p.activation_type = "HCD";
    DeconvolvedSpectrum d(301);
    FragmentAnalysis::ProteoformMatch m;
    m.score = 1.0;
    m.identification_score = 1.0;
    m.region_start = -1;
    m.region_end = -1;
    m.matched_protein = "synthetic";
    m.proteoform_sequence = WINNER_SEQ;
    for (const GateFrag& g : frags)
    {
      PeakGroup pg = (g.charges.size() == 1) ? makeSyntheticPeakGroup(g.mass / g.charges[0] + PROTON_MASS, g.mass, g.charges[0])
                                             : makeEnvelope(g.mass, g.charges);
      if (g.iso_cos >= 0.0) pg.setIsotopeCosine(static_cast<float>(g.iso_cos));
      d.push_back(pg);
      FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
      fm.ion_type = g.ion_type;
      fm.ion_index = g.ion_index;
      fm.observed_mass = g.mass;
      m.fragments.push_back(fm);
    }
    if (ambiguous)
    {
      FragmentAnalysis::PTMSite site;
      site.start_position = 4;
      site.end_position = 5;
      site.position = 4;
      site.mass_shift = 42.0106;
      m.ptm_sites.push_back(site);
    }
    for (size_t i = 0; i < sites.size(); ++i)
    {
      FragmentAnalysis::PTMSite site;
      site.start_position = sites[i].first;
      site.end_position = sites[i].second;
      site.position = sites[i].first;
      site.mass_shift = 42.0106 + 37.9557 * static_cast<double>(i);   // acetyl, then phospho: one shift per modification
      m.ptm_sites.push_back(site);
    }
    tracker.feedScan(precursor_id, 2, p, 301, d, m, 1.0, makeMs2Ctx());
    tracker.finalizeMS2(precursor_id);
  }

  // Four fragments of PEPTIDEK (L = 8), each with a bond inside it, at charge 2; b5, b6 and y6 contain residues 4..5.
  const std::vector<GateFrag> FOUR_FRAGS{{"b", 3, 330.0, {2}}, {"b", 5, 550.0, {2}}, {"b", 6, 700.0, {2}}, {"y", 6, 900.0, {2}}};

  /// A gate that admits every isolatable charge: nothing sent, no SNR test, no charge floor, room for 100 cleavages.
  Ms3PlanGate openGate()
  {
    Ms3PlanGate g;
    g.max_targets = 100;
    return g;
  }

  using TargetKey = std::tuple<std::string, int, int, double>;
  std::multiset<TargetKey> keysOf(const std::vector<Ms3Target>& plan)
  {
    std::multiset<TargetKey> out;
    for (const Ms3Target& t : plan) out.insert(TargetKey{t.ion_type, t.ion_index, t.frag_charge, t.frag_mz});
    return out;
  }

  bool planHas(const std::vector<Ms3Target>& plan, const std::string& ion_type, int ion_index)
  {
    for (const Ms3Target& t : plan)
      if (t.ion_type == ion_type && t.ion_index == ion_index) return true;
    return false;
  }

  // ---- ADR-0069: the yield gate -------------------------------------------------------------------------------
  /// The dispatching MS2 of the yield sections: stage 0 of every MS3 (precursor intensity, charge, mass).
  ScanCommand yieldCtx()
  {
    ScanCommand ctx = makeMs2Ctx();
    ctx.mono_mass = 10000.0;
    ctx.precursor_intensity = 5.0e5;
    ctx.stages[0].charge_state = 10;
    return ctx;
  }

  /// openGate() plus the yield gate @p yg against @p ctx.
  Ms3PlanGate yieldGate(const YieldGateConfig& yg, const ScanCommand& ctx)
  {
    Ms3PlanGate g = openGate();
    g.yield = &yg;
    g.ms_ctx = &ctx;
    return g;
  }

  /// @p frags with the fragments named in @p low at isotope cosine @p low_iso and the rest at 0.99.
  std::vector<GateFrag> withIsoCos(std::vector<GateFrag> frags, const std::set<std::pair<std::string, int>>& low, double low_iso)
  {
    for (GateFrag& g : frags) g.iso_cos = low.count({g.ion_type, g.ion_index}) ? low_iso : 0.99;
    return frags;
  }

  /// The predicted bond yield of @p g as the planner computes it: intensity 0 (fixture groups carry none), its mass,
  /// isotope cosine and ion index, stage 0 from @p ctx.
  double yieldOf(const GateFrag& g, const ScanCommand& ctx)
  {
    FragmentAnalysis::FragmentScores s1;
    s1.mono_mass = g.mass;
    s1.iso_cos = g.iso_cos;
    s1.precursor_intensity = 0.0;
    return predictBondYield(YieldCoefficients{}, yieldInputsOf(ctx, s1, g.ion_index)).e;
  }

  /// A threshold between the yield of @p frags' fragment @p low (set below) and the lowest yield of every other fragment.
  double thresholdBetween(const std::vector<GateFrag>& frags, const std::pair<std::string, int>& low, const ScanCommand& ctx)
  {
    const std::vector<GateFrag> set = withIsoCos(frags, {low}, 0.85);
    double low_e = 0.0, high_e = std::numeric_limits<double>::max();
    for (const GateFrag& g : set)
    {
      if (g.ion_type == low.first && g.ion_index == low.second) low_e = yieldOf(g, ctx);
      else high_e = std::min(high_e, yieldOf(g, ctx));
    }
    if (!(low_e < high_e)) return -1.0;   // the fixture does not separate them: a threshold that admits all fails the sections
    return std::sqrt(low_e * high_e);
  }

  int countOf(const std::string& hay, const std::string& needle)
  {
    int n = 0;
    for (size_t pos = hay.find(needle); pos != std::string::npos; pos = hay.find(needle, pos + needle.size())) ++n;
    return n;
  }

  // ---- ADR-0073: the interleaved send order -------------------------------------------------------------------
  /// @p frags with the masses and isotope cosines that make the list order the ranked order under a yield gate that
  /// ranks by yield: a place down the list is 1500 Da heavier and 0.0375 lower in isotope cosine, and both lower the
  /// predicted yield. (The model clips the cosine to [0.85, 1.0005], so the cosine alone cannot separate five places.)
  std::vector<GateFrag> rankedFrags(std::vector<GateFrag> frags)
  {
    for (size_t i = 0; i < frags.size(); ++i)
    {
      frags[i].mass = 1000.0 + 1500.0 * static_cast<double>(i);
      frags[i].iso_cos = 1.0 - 0.0375 * static_cast<double>(i);
    }
    return frags;
  }

  /// rankedFrags' claim: the predicted yields of @p frags fall strictly down the list and none is zero.
  bool yieldsFall(const std::vector<GateFrag>& frags, const ScanCommand& ctx)
  {
    if (frags.empty() || !(yieldOf(frags.back(), ctx) > 0.0)) return false;
    for (size_t i = 1; i < frags.size(); ++i)
      if (!(yieldOf(frags[i - 1], ctx) > yieldOf(frags[i], ctx))) return false;
    return true;
  }

  /// A send order as text, "0 3 1 4 2", so a failed comparison prints both orders.
  std::string orderText(const std::vector<size_t>& order)
  {
    std::string out;
    for (size_t i : order) out += (out.empty() ? "" : " ") + std::to_string(i);
    return out;
  }

  /// A plan's ions in send order as text, "y6 b3 y5", one entry per target (under separate, one per charge).
  std::string ionsText(const std::vector<Ms3Target>& plan)
  {
    std::string out;
    for (const Ms3Target& t : plan) out += (out.empty() ? "" : " ") + t.ion_type + std::to_string(t.ion_index);
    return out;
  }

  /// One planner fixture of the interleaved sections: what is planned, on which modifications, behind which snapshot.
  struct SendOrderCase
  {
    std::string mode;
    int budget = 10;
    std::string fragment_charges;
    std::vector<std::pair<int, int>> sites;          ///< the ambiguous modifications (residue ranges)
    const Ms3RegistrySnapshot* snapshot = nullptr;   ///< the proteoform's sent fragments (null: not governed)
  };

  /// The plan of @p frags under @p c with characterization.order @p order ("" = ranked), on a fresh tracker, behind a
  /// yield gate every fragment clears (half the lowest yield) that ranks by yield.
  std::vector<Ms3Target> planInSendOrder(const SendOrderCase& c, const std::string& order, int precursor_id,
                                         const std::vector<GateFrag>& frags, const ScanCommand& ctx)
  {
    Config cfg{gateConfig(c.mode, c.budget, c.fragment_charges, order)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedGateWinner(tracker, precursor_id, frags, /*ambiguous=*/false, c.sites);
    YieldGateConfig by_yield;
    by_yield.threshold = yieldOf(frags.back(), ctx) / 2.0;
    Ms3PlanGate gate = yieldGate(by_yield, ctx);
    gate.snapshot = c.snapshot;
    return tracker.planNextScans(precursor_id, &gate);
  }
}

START_TEST(ProteoformTracker_Trajectory, "$Id$")

/////////////////////////////////////////////////////////////
// Per-precursor trajectory: MS2 baseline row + additive (non-resetting) MS3 fold rows.
/////////////////////////////////////////////////////////////
START_SECTION(ms2_baseline_then_accumulating_ms3_folds)
{
  const int precursor_id = 7;

  // Log FOLDER embedded in the Config JSON so the engine opens its streams (Config.cpp reads
  // runtime.log_dir; IdaLogger joins the fixed basename pooled_identification.tsv onto it). Use a
  // RELATIVE folder in the test CWD (NOT an absolute temp path): an absolute Windows path contains
  // backslashes that would form invalid JSON escapes when embedded below. The streams open in
  // APPEND mode and IdaLogger never creates directories, so wipe-and-create up front.
  const std::string log_dir = "testlogs/pt_trajectory";
  File::removeDirRecursively(log_dir);
  File::makeDir(log_dir);
  const std::string pooled_path = log_dir + "/pooled_identification.tsv";

  std::string cfg_json = std::string(R"({
    "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 1, "max_charge": 50, "min_mass": 100, "max_mass": 50000, "tol": [10, 10, 10] },
    "precursor_selection": { "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 3 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
      "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25, "resolution": 120000 }
    },
    "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 } },
    "files": { "target_logs": [], "inclusion_list": "" },
    "identification": { "reference": { "sequence": "PEPTIDEK" }, "flashtnt": { "max_blind_mod_count": 2, "max_mod_mass": 700 }, "characterization": { "mode": "coverage", "max_targets": 10 } },
    "runtime": { "log_dir": ")") + log_dir + R"(" } })";

  Config cfg{cfg_json};
  IdaLogger logger(cfg);                 // opens the pooled stream (writes the header)
  ProteoformTracker tracker(cfg, logger);
  std::vector<ProteoformTracker::PooledUpdate> observed;   // the operator view's feed: one call per row written
  tracker.setPooledObserver([&](const ProteoformTracker::PooledUpdate& u) { observed.push_back(u); });

  // Real winner (PEPTIDEK) b/y ladder. A NON-winner MS2 scan is now re-matched against the winner
  // ladder (Change 1: winner-anchored pooling), so its raw peak mass MUST equal a real winner ion or it
  // is DROPPED. Compute the ladder with the SAME function the engine uses -> the fed peak matches at 0 ppm.
  std::map<char, std::vector<double>> winner_ladder;
  FragmentAnalysis::computePTMAdjustedFragmentMasses("PEPTIDEK", {}, {"b", "y"}, winner_ladder);
  const double y2_mass = winner_ladder['y'][1];   // y2 -> winner key ("y", 2): distinct from b6 / b3 / y3

  ScanCommand ms2_ctx = makeMs2Ctx();

  // ---------------------------------------------------------------------------------------------
  // MS2 phase (baseline): feed two MS2 scans, then finalize -> ONE "MS2"-tagged baseline row.
  //   scan 101 = WINNER (b6, mass 700): pooled verbatim via the winner path -> key ("b", 6).
  //   scan 102 = NON-winner (raw peak = real y2 mass): re-matched against the winner ladder (Change 1)
  //              -> key ("y", 2), a distinct position -> two distinct MappedFragments.
  // ---------------------------------------------------------------------------------------------
  {
    Ms2Params p;
    p.collision_energy = 29.0;
    p.activation_type = "HCD";
    p.reaction_time = 0.0;

    DeconvolvedSpectrum d101(101);
    d101.push_back(makeSyntheticPeakGroup(700.0 / 2.0 + 1.0, 700.0, 2));
    tracker.feedScan(precursor_id, 2, p, 101, d101, makeMatch("b", 6, 700.0), 1.0, ms2_ctx);

    DeconvolvedSpectrum d102(102);
    d102.push_back(makeSyntheticPeakGroup(y2_mass / 2.0 + 1.0, y2_mass, 2));
    // Non-winner: only its raw peak (y2_mass) is used -> re-matched to winner ("y", 2). The match arg's
    // fragments are ignored for a non-winner; score/sequence keep it a (losing) winner candidate.
    tracker.feedScan(precursor_id, 2, p, 102, d102, makeMatch("y", 2, y2_mass), 1.0, ms2_ctx);
  }

  tracker.finalizeMS2(precursor_id);

  // (A) after finalize: winner identified, two MS2 fragments mapped, update_index == 1.
  const ProteoformModel* mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL(mdl->proteoform_sequence, std::string(WINNER_SEQ))
  TEST_EQUAL(mdl->update_index, 1)
  TEST_EQUAL((int)mdl->fragments.size(), 2)   // b6, y6

  Ms2Params ms3_params;       // parent MS2 params carried with the MS3 observation
  ms3_params.collision_energy = 29.0;
  ms3_params.activation_type = "HCD";
  ms3_params.reaction_time = 0.0;

  // ---------------------------------------------------------------------------------------------
  // MS3 fold #1: stage one MS3 scan (equiv b3) then foldMs3 -> ADDITIVE, NO reset.
  // A DISTINCT equiv ion (b3) so the model gains a NEW fragment (size 2 -> 3). foldMs3 emits one row
  // tagged trigger = "b3", trigger_scan_id = "AAB", update_index = 2.
  // ---------------------------------------------------------------------------------------------
  {
    DeconvolvedSpectrum d201(201);
    d201.push_back(makeSyntheticPeakGroup(350.0 / 2.0 + 1.0, 350.0, 2));
    tracker.feedScan(precursor_id, 3, ms3_params, 201, d201, makeMatchMs3("b", 3, 350.0, 350.0), -1.0, ms2_ctx);
    tracker.foldMs3(precursor_id, "b3", "AAB");
  }

  // (A) after fold #1: fragment ADDED (NOT reset), update_index bumped.
  mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL(mdl->update_index, 2)
  TEST_EQUAL((int)mdl->fragments.size(), 3)     // b6, y6, + b3
  TEST_TRUE((int)mdl->fragments.size() > 2)      // a reset would DROP this to 1 -> this fails

  // ---------------------------------------------------------------------------------------------
  // MS3 fold #2: stage one MS3 scan (equiv y3) then foldMs3 -> ADDITIVE, NO reset.
  // trigger = "y3", trigger_scan_id = "AAC", update_index = 3, fragments 3 -> 4.
  // ---------------------------------------------------------------------------------------------
  {
    DeconvolvedSpectrum d202(202);
    d202.push_back(makeSyntheticPeakGroup(450.0 / 2.0 + 1.0, 450.0, 2));
    tracker.feedScan(precursor_id, 3, ms3_params, 202, d202, makeMatchMs3("y", 3, 450.0, 450.0), -1.0, ms2_ctx);
    tracker.foldMs3(precursor_id, "y3", "AAC");
  }

  // (A) after fold #2: another fragment ADDED, update_index bumped.
  mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL(mdl->update_index, 3)
  TEST_EQUAL((int)mdl->fragments.size(), 4)     // b6, y6, b3, + y3

  // ---------------------------------------------------------------------------------------------
  // (B) Pooled trajectory file: 3 data rows, one per emitRow_ call (finalize + 2x foldMs3).
  // Columns are resolved BY HEADER NAME (colOf, below), so this test is agnostic to the pooled column
  // order — a reorder in IdaLogger's pooled writer needs no change here. The pooled stream has 20 columns.
  // ---------------------------------------------------------------------------------------------
  std::ifstream f(pooled_path.c_str());
  TEST_TRUE(f.good())

  std::string header;
  TEST_TRUE(static_cast<bool>(std::getline(f, header)))
  const std::vector<std::string> header_cols = splitTab(header);
  TEST_EQUAL((int)header_cols.size(), 20)   // 14 base + grouped fragment-mass table (masses|ions|measured|theoretical|diff_da|diff_ppm) + protein (ADR-0061 push 4b)

  std::vector<std::vector<std::string>> rows;
  std::string line;
  while (std::getline(f, line))
  {
    if (line.empty()) continue;
    rows.push_back(splitTab(line));
  }
  TEST_EQUAL((int)rows.size(), 3)   // exactly one row per finalize + foldMs3 + foldMs3
  ABORT_IF((int)rows.size() != 3)

  // Resolve pooled columns BY HEADER NAME (order-agnostic to the pooled layout).
  auto colOf = [&](const std::string& name) {
    for (int j = 0; j < (int)header_cols.size(); ++j) if (header_cols[j] == name) return j;
    return -1;
  };
  const int CI_COV = colOf("coverage_pct"), CI_N_FRAG = colOf("n_fragments"),
            CI_MASSES = colOf("combined_ms2_frame_masses"), CI_IONS = colOf("combined_ms2_fragment_ions"),
            CI_UPDATE = colOf("update_index"), CI_PRECID = colOf("precursor_id"),
            CI_TRIGGER = colOf("trigger"), CI_TRIGGER_SCAN = colOf("trigger_scan_id");
  TEST_TRUE(CI_COV >= 0 && CI_N_FRAG >= 0 && CI_MASSES >= 0 && CI_IONS >= 0 &&
            CI_UPDATE >= 0 && CI_PRECID >= 0 && CI_TRIGGER >= 0 && CI_TRIGGER_SCAN >= 0)

  // Every data row must carry all 20 fields (the trailing protein is non-empty here).
  for (const std::vector<std::string>& r : rows) TEST_EQUAL((int)r.size(), 20)

  // combined_ms2_fragment_ions (col 10) aligns 1-for-1 with combined_ms2_frame_masses (col 9), and each
  // label looks like a fragment ion (ion-type letter + index).
  auto semiCount = [](const std::string& s) { if (s.empty()) return 0; int n = 1; for (char c : s) if (c == ';') ++n; return n; };
  for (const std::vector<std::string>& r : rows)
  {
    TEST_EQUAL(semiCount(r[CI_IONS]), semiCount(r[CI_MASSES]))
    TEST_TRUE(!r[CI_IONS].empty() && r[CI_IONS][0] >= 'a' && r[CI_IONS][0] <= 'z')
  }

  // row0 = MS2 baseline: trigger "MS2"; trigger_scan_id = encode(101) (winner = first eligible scan);
  //        update_index "1".
  TEST_EQUAL(rows[0][CI_TRIGGER], std::string("MS2"))
  TEST_EQUAL(rows[0][CI_TRIGGER_SCAN], ScanCommandQueue::encode(101))
  TEST_EQUAL(rows[0][CI_UPDATE], std::string("1"))
  TEST_EQUAL(rows[0][CI_PRECID], std::string("7"))

  // row1 = MS3 fold #1: trigger "b3"; trigger_scan_id "AAB"; update_index "2".
  TEST_EQUAL(rows[1][CI_TRIGGER], std::string("b3"))
  TEST_EQUAL(rows[1][CI_TRIGGER_SCAN], std::string("AAB"))
  TEST_EQUAL(rows[1][CI_UPDATE], std::string("2"))

  // row2 = MS3 fold #2: trigger "y3"; trigger_scan_id "AAC"; update_index "3".
  TEST_EQUAL(rows[2][CI_TRIGGER], std::string("y3"))
  TEST_EQUAL(rows[2][CI_TRIGGER_SCAN], std::string("AAC"))
  TEST_EQUAL(rows[2][CI_UPDATE], std::string("3"))

  // n_fragments parsed as int is non-decreasing across rows (here strictly 2 -> 3 -> 4). A reset
  // (old finalize) would log 1 on the MS3 rows, breaking the monotone growth.
  const int nf0 = std::stoi(rows[0][CI_N_FRAG]);
  const int nf1 = std::stoi(rows[1][CI_N_FRAG]);
  const int nf2 = std::stoi(rows[2][CI_N_FRAG]);
  TEST_EQUAL(nf0, 2)
  TEST_EQUAL(nf1, 3)
  TEST_EQUAL(nf2, 4)
  TEST_TRUE(nf1 >= nf0)
  TEST_TRUE(nf2 >= nf1)

  // coverage_pct (col 4) non-decreasing across the trajectory — the engine's LOGGED coverage, asserted
  // from the produced file rather than ProteoformModel::coveragePct() (the struct is not OPENMS_DLLAPI, so
  // its out-of-line method is not linkable from a test exe; the file column is the observable contract).
  const double cov0 = std::stod(rows[0][CI_COV]);
  const double cov1 = std::stod(rows[1][CI_COV]);
  const double cov2 = std::stod(rows[2][CI_COV]);
  TEST_TRUE(cov1 >= cov0)
  TEST_TRUE(cov2 >= cov1)

  // The pooled observer saw exactly the rows written, field for field (the board shows what the log carries).
  // FAILS UNDER: a call before the write, a missed emit site, or fields that differ from the row.
  const int CI_NOMINAL = colOf("nominal_mass");
  TEST_TRUE(CI_NOMINAL >= 0)
  TEST_EQUAL(observed.size(), rows.size())
  for (size_t k = 0; k < observed.size() && k < rows.size(); ++k)
  {
    TEST_EQUAL(std::to_string(observed[k].precursor_id), rows[k][CI_PRECID])
    TEST_EQUAL(std::to_string(observed[k].nominal_mass), rows[k][CI_NOMINAL])
    TEST_EQUAL(observed[k].trigger, rows[k][CI_TRIGGER])
    TEST_EQUAL(observed[k].trigger_scan_id, rows[k][CI_TRIGGER_SCAN])
    TEST_TRUE(std::fabs(observed[k].coverage - std::stod(rows[k][CI_COV])) < 1e-4)
  }
}
END_SECTION

START_SECTION(ms3_flipped_ion_takes_its_own_peak)
{
  // A flipped MS3 ion (a suffix sub-fragment of a b precursor, mapped to a prefix equivalent) has two masses:
  // the measured one, which is the peak in the MS3 scan, and the adjusted one, which is the equivalent ion's
  // mass in the proteoform frame and the mass the model pools. The fragment's own peak is found by the
  // MEASURED mass. Searched by the adjusted mass it is never in tolerance, and the closest-overall fallback
  // handed the fragment the m/z, charge and intensity of whatever unrelated peak lay nearest.
  const int precursor_id = 8;
  const std::string log_dir = "testlogs/pt_trajectory_flip";
  File::removeDirRecursively(log_dir);
  File::makeDir(log_dir);
  std::string cfg_json = std::string(R"({
    "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 1, "max_charge": 50, "min_mass": 100, "max_mass": 50000, "tol": [10, 10, 10] },
    "precursor_selection": { "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 3 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
      "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25, "resolution": 120000 }
    },
    "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 } },
    "files": { "target_logs": [], "inclusion_list": "" },
    "identification": { "reference": { "sequence": "PEPTIDEK" }, "flashtnt": { "max_blind_mod_count": 2, "max_mod_mass": 700 }, "characterization": { "mode": "coverage", "max_targets": 10 } },
    "runtime": { "log_dir": ")") + log_dir + R"(" } })";

  Config cfg{cfg_json};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  ScanCommand ms2_ctx = makeMs2Ctx();

  Ms2Params p;
  p.collision_energy = 29.0;
  p.activation_type = "HCD";
  p.reaction_time = 0.0;
  DeconvolvedSpectrum d101(101);
  d101.push_back(makeSyntheticPeakGroup(700.0 / 2.0 + 1.0, 700.0, 2));
  tracker.feedScan(precursor_id, 2, p, 101, d101, makeMatch("b", 6, 700.0), 1.0, ms2_ctx);
  tracker.finalizeMS2(precursor_id);

  // The MS3 scan: the fragment's own peak at its measured mass 620 (3+), and an unrelated peak at 352 (2+)
  // that happens to lie next to the fragment's ADJUSTED mass 350.
  DeconvolvedSpectrum d201(201);
  d201.push_back(makeEnvelope(620.0, {3}));
  d201.push_back(makeEnvelope(352.0, {2}));
  tracker.feedScan(precursor_id, 3, p, 201, d201, makeMatchMs3("b", 3, 350.0, 620.0), -1.0, ms2_ctx);
  tracker.foldMs3(precursor_id, "b6", "AAB");

  const ProteoformModel* mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  auto it = mdl->fragments.find(FragmentKey{"b", 3});
  TEST_TRUE(it != mdl->fragments.end())
  ABORT_IF(it == mdl->fragments.end())
  TEST_TRUE(it->second.best_ms3.has_value())
  ABORT_IF(!it->second.best_ms3.has_value())
  const FragmentObservation& o = *it->second.best_ms3;
  TEST_EQUAL(o.frag_charge, 3)                              // its own peak (was 2: the peak nearest 350)
  TEST_REAL_SIMILAR(o.frag_mz, (620.0 + 0.5 * ISOTOPE_STEP) / 3.0 + PROTON_MASS)   // the centre of its two 3+ isotopes
  TEST_REAL_SIMILAR(o.measured_mass, 620.0)
  TEST_REAL_SIMILAR(o.observed_mass, 350.0)                 // the model still pools the adjusted mass
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0066 G1. A null gate is today's plan, and a gate that admits everything changes nothing, under each of the
// three mapped objectives.
//
// FAILS UNDER: a gate path that re-orders, re-anchors or re-notches an admitted fragment (the anchor must stay best_ms2
// when best_ms2 is admitted), or one that drops a fragment it should admit.
/////////////////////////////////////////////////////////////
START_SECTION(gate_null_and_open_equal_today)
{
  const std::vector<GateFrag> two_amb{{"b", 6, 700.0, {2}}, {"y", 6, 900.0, {2}}};
  for (const std::string mode : {"ambiguity", "coverage", "fragments"})
  {
    Config cfg{gateConfig(mode, 10)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    if (mode == "ambiguity") feedGateWinner(tracker, 41, two_amb, /*ambiguous=*/true);
    else feedGateWinner(tracker, 41, FOUR_FRAGS);

    const std::vector<Ms3Target> today = tracker.planNextScans(41);
    const std::vector<Ms3Target> null_gate = tracker.planNextScans(41, nullptr);
    const Ms3PlanGate open = openGate();
    const std::vector<Ms3Target> open_gate = tracker.planNextScans(41, &open);
    TEST_EQUAL(today.empty(), false)   // non-vacuity: there is a plan to compare
    TEST_EQUAL(today.size(), (mode == "ambiguity") ? 2u : 4u)
    TEST_EQUAL(keysOf(null_gate) == keysOf(today), true)
    TEST_EQUAL(keysOf(open_gate) == keysOf(today), true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0066 G2. A refused fragment hands its slot to the next candidate in the same order -- under fragments, and
// through the ambiguity round-robin. The refusal is the 2x rule: the fragment the ungated plan takes was sent at
// intensity 1, and it is seen at 0 now.
//
// FAILS UNDER: a refusal applied after the budget is spent (the plan would be empty), or a refused fragment still
// taking the slot.
/////////////////////////////////////////////////////////////
START_SECTION(gate_refused_strongest_hands_its_slot_on)
{
  const std::vector<GateFrag> two_amb{{"b", 6, 700.0, {2}}, {"y", 6, 900.0, {2}}};
  for (const std::string mode : {"fragments", "ambiguity"})
  {
    Config cfg{gateConfig(mode, 1)};   // budget 1: only the first fragment fits
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    if (mode == "ambiguity") feedGateWinner(tracker, 42, two_amb, /*ambiguous=*/true);
    else feedGateWinner(tracker, 42, FOUR_FRAGS);

    const std::vector<Ms3Target> today = tracker.planNextScans(42);
    TEST_EQUAL(today.size(), 1)
    ABORT_IF(today.size() != 1)

    Ms3RegistrySnapshot snap;
    snap.found = true;
    snap.sent[Ms3Cleavage{today[0].ion_type, today[0].ion_index}][today[0].frag_charge] = 1.0;
    Ms3PlanGate gate = openGate();
    gate.snapshot = &snap;
    const std::vector<Ms3Target> gated = tracker.planNextScans(42, &gate);
    TEST_EQUAL(gated.size(), 1)   // the slot is still filled ...
    ABORT_IF(gated.size() != 1)
    TEST_EQUAL(gated[0].ion_type == today[0].ion_type && gated[0].ion_index == today[0].ion_index, false)   // ... by the next
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0066 G3. A charge whose window SNR is not above snr_min is refused. Every fixture charge here carries no
// signal (window SNR 0): with the dispatching spectrum and the default snr_min 1 no fragment enters a pool; with the
// same spectrum and snr_min -1 every fragment does, so the threshold is what refused them.
//
// FAILS UNDER: no SNR test in the gate (the first plan is not empty), or a test that ignores snr_min.
/////////////////////////////////////////////////////////////
START_SECTION(gate_window_snr_at_or_below_one_is_refused)
{
  Config cfg{gateConfig("fragments", 10)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 43, FOUR_FRAGS);

  // The dispatching MS2: one raw peak under every fragment's isolation window.
  MSSpectrum ms2;
  for (const GateFrag& g : FOUR_FRAGS)
  {
    Peak1D pk;
    pk.setMZ(g.mass / g.charges[0] + PROTON_MASS);
    pk.setIntensity(5000.0f);
    ms2.push_back(pk);
  }
  ms2.sortByPosition();

  Ms3PlanGate gate = openGate();
  gate.spectrum = &ms2;
  TEST_EQUAL(tracker.planNextScans(43, &gate).size(), 0)   // SNR 0 <= 1: all refused
  gate.snr_min = -1.0;
  TEST_EQUAL(tracker.planNextScans(43, &gate).size(), 4)   // the same spectrum, a threshold below 0: all admitted
  TEST_EQUAL(tracker.planNextScans(43).size(), 4)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0066 G4. Single mode, the strongest charge refused: the next admitted charge of the same fragment stands in
// (ISSUE 2). b6 is seen at 2+ (the representative) and 3+; a charge floor of 3 refuses 2+, a floor of 4 both.
//
// FAILS UNDER: an anchor that stays best_ms2 whatever the gate says (the gated target would be 2+), or a fallback
// that ignores the floor (the floor-4 plan would not be empty).
/////////////////////////////////////////////////////////////
START_SECTION(gate_single_mode_charge_fallback)
{
  Config cfg{gateConfig("fragments", 10)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 44, {{"b", 6, 700.0, {2, 3}}});

  const std::vector<Ms3Target> today = tracker.planNextScans(44);
  TEST_EQUAL(today.size(), 1)
  ABORT_IF(today.size() != 1)
  TEST_EQUAL(today[0].frag_charge, 2)

  Ms3PlanGate gate = openGate();
  gate.charge_floor = 3;
  const std::vector<Ms3Target> fallback = tracker.planNextScans(44, &gate);
  TEST_EQUAL(fallback.size(), 1)
  ABORT_IF(fallback.size() != 1)
  TEST_EQUAL(fallback[0].ion_type, std::string("b"))
  TEST_EQUAL(fallback[0].ion_index, 6)
  TEST_EQUAL(fallback[0].frag_charge, 3)
  // The 3+ window: the centre of the two 3+ isotopes, not the 2+ one.
  TEST_REAL_SIMILAR(fallback[0].frag_mz, (700.0 + 0.5 * ISOTOPE_STEP) / 3.0 + PROTON_MASS)
  TEST_EQUAL(fallback[0].notches.empty(), true)   // single mode never notches

  gate.charge_floor = 4;
  TEST_EQUAL(tracker.planNextScans(44, &gate).size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0066 G4b. Separate mode drops only the refused charges: b6 seen at 2+, 3+ and 4+ is three targets ungated; a
// charge floor of 3 refuses 2+ alone, so the governed plan is the 3+ and 4+ targets of the same fragment.
//
// FAILS UNDER: a gate that drops the whole fragment when its representative charge is refused, or one that still
// emits the refused charge from the fragment's full envelope.
/////////////////////////////////////////////////////////////
START_SECTION(gate_separate_drops_only_refused_charges)
{
  Config cfg{gateConfig("fragments", 10, "separate")};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 46, {{"b", 6, 700.0, {2, 3, 4}}});

  auto chargesOf = [](const std::vector<Ms3Target>& plan) {
    std::set<int> z;
    for (const Ms3Target& t : plan)
      if (t.ion_type == "b" && t.ion_index == 6) z.insert(t.frag_charge);
    return z;
  };
  const std::vector<Ms3Target> today = tracker.planNextScans(46);
  TEST_EQUAL(today.size(), 3)
  TEST_EQUAL(chargesOf(today) == std::set<int>({2, 3, 4}), true)   // non-vacuity: every charge is a target

  Ms3PlanGate gate = openGate();
  gate.charge_floor = 3;
  const std::vector<Ms3Target> gated = tracker.planNextScans(46, &gate);
  TEST_EQUAL(gated.size(), 2)
  TEST_EQUAL(chargesOf(gated) == std::set<int>({3, 4}), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0066 G4c. Multiplexed: the anchor and the notches come from the admitted charges only. Ungated, b6 is one target
// anchored at 2+ with 3+ and 4+ as notches; with 2+ refused the target is anchored at 3+ and notches 4+ alone.
//
// FAILS UNDER: notches drawn from the fragment's full envelope (2+ would come back as a notch), or an anchor that stays
// the refused representative charge.
/////////////////////////////////////////////////////////////
START_SECTION(gate_multiplexed_anchor_and_notches_are_admitted_charges)
{
  Config cfg{gateConfig("fragments", 10, "multiplexed")};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 47, {{"b", 6, 700.0, {2, 3, 4}}});

  auto notchCharges = [](const Ms3Target& t) {
    std::set<int> z;
    for (const NotchCandidate& n : t.notches) z.insert(n.charge);
    return z;
  };
  const std::vector<Ms3Target> today = tracker.planNextScans(47);
  TEST_EQUAL(today.size(), 1)
  ABORT_IF(today.size() != 1)
  TEST_EQUAL(today[0].frag_charge, 2)
  TEST_EQUAL(notchCharges(today[0]) == std::set<int>({3, 4}), true)   // non-vacuity: the envelope is notched

  Ms3PlanGate gate = openGate();
  gate.charge_floor = 3;
  const std::vector<Ms3Target> gated = tracker.planNextScans(47, &gate);
  TEST_EQUAL(gated.size(), 1)
  ABORT_IF(gated.size() != 1)
  TEST_EQUAL(gated[0].frag_charge, 3)
  TEST_EQUAL(notchCharges(gated[0]) == std::set<int>({4}), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0066 G5. The cleavage cap counts SENT cleavages, not re-acquisitions: with b6 sent and max_targets 2 the
// proteoform has room for one new cleavage, and b6's re-acquisition (sent at 0, seen at 0: admitted) takes none of it.
// With b6 and y6 sent and max_targets 2 there is no room: only the two re-acquisitions are planned.
//
// FAILS UNDER: a re-acquisition that consumes room (the first plan holds 1 target), a cap that counts planned targets
// instead of sent cleavages, or a capped fragment returning "budget full" (the callers would stop at the first one).
/////////////////////////////////////////////////////////////
START_SECTION(gate_cap_counts_sent_cleavages_not_reacquisitions)
{
  Config cfg{gateConfig("fragments", 10)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 45, FOUR_FRAGS);
  TEST_EQUAL(tracker.planNextScans(45).size(), 4)

  Ms3RegistrySnapshot snap;
  snap.found = true;
  snap.sent[Ms3Cleavage{"b", 6}][2] = 0.0;
  Ms3PlanGate gate = openGate();
  gate.snapshot = &snap;
  gate.max_targets = 2;
  const std::vector<Ms3Target> one_new = tracker.planNextScans(45, &gate);
  TEST_EQUAL(one_new.size(), 2)
  TEST_EQUAL(planHas(one_new, "b", 6), true)

  snap.sent[Ms3Cleavage{"y", 6}][2] = 0.0;
  const std::vector<Ms3Target> none_new = tracker.planNextScans(45, &gate);
  TEST_EQUAL(none_new.size(), 2)
  TEST_EQUAL(planHas(none_new, "b", 6), true)
  TEST_EQUAL(planHas(none_new, "y", 6), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0069: the yield gate in the mapped planner. Fixture PeakGroups carry no charge intensity (every I_frag clips to
// 372), so a fragment's predicted bond yield is set through its isotope cosine (GateFrag::iso_cos; the model's
// steepest input) and its ion index; every charge of one fragment then has the same yield. yieldCtx() is the
// dispatching MS2 every MS3 inherits stage 0 from.
/////////////////////////////////////////////////////////////

// T9. A yield gate whose threshold every target clears changes nothing, ungoverned (a gate with no snapshot) and
// governed (a snapshot with nothing sent), under each mapped objective.
// FAILS UNDER: a yield walk that re-orders or re-anchors an admitted fragment, or drops one it should admit.
START_SECTION(yield_gate_open_equals_no_gate)
{
  const std::vector<GateFrag> two_amb{{"b", 6, 700.0, {2}}, {"y", 6, 900.0, {2}}};
  const ScanCommand ctx = yieldCtx();
  YieldGateConfig open_yield;
  open_yield.threshold = 1e-12;
  for (const std::string mode : {"ambiguity", "coverage", "fragments"})
  {
    Config cfg{gateConfig(mode, 10)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    if (mode == "ambiguity") feedGateWinner(tracker, 141, two_amb, /*ambiguous=*/true);
    else feedGateWinner(tracker, 141, FOUR_FRAGS);

    const std::vector<Ms3Target> today = tracker.planNextScans(141);
    TEST_EQUAL(today.size(), (mode == "ambiguity") ? 2u : 4u)   // non-vacuity
    Ms3PlanGate ungoverned = yieldGate(open_yield, ctx);
    TEST_EQUAL(keysOf(tracker.planNextScans(141, &ungoverned)) == keysOf(today), true)
    Ms3RegistrySnapshot snap;
    snap.found = true;
    Ms3PlanGate governed = yieldGate(open_yield, ctx);
    governed.snapshot = &snap;
    Ms3PlanGate governed_only = openGate();
    governed_only.snapshot = &snap;
    TEST_EQUAL(keysOf(tracker.planNextScans(141, &governed)) == keysOf(tracker.planNextScans(141, &governed_only)), true)
  }
}
END_SECTION

// T10. A yield-refused fragment hands its place on: budget 2, the plan's first fragment set below the threshold, the
// plan is the next two, and one refused line names it.
// FAILS UNDER: a refusal applied after the budget is spent (the plan would hold one target), or a refused fragment
// still taking its place.
START_SECTION(yield_refused_strongest_hands_its_place_on)
{
  const ScanCommand ctx = yieldCtx();
  // The plan order, learnt with every fragment at the same isotope cosine (the order is the intensity order, which the
  // cosine does not touch).
  std::vector<std::pair<std::string, int>> order;
  {
    Config cfg{gateConfig("fragments", 10)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedGateWinner(tracker, 142, withIsoCos(FOUR_FRAGS, {}, 0.99));
    for (const Ms3Target& t : tracker.planNextScans(142)) order.emplace_back(t.ion_type, t.ion_index);
  }
  TEST_EQUAL(order.size(), 4)
  ABORT_IF(order.size() != 4)

  Config cfg{gateConfig("fragments", 2)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 142, withIsoCos(FOUR_FRAGS, {order[0]}, 0.85));
  const std::vector<Ms3Target> today = tracker.planNextScans(142);
  TEST_EQUAL(today.size(), 2)
  TEST_EQUAL(planHas(today, order[0].first, order[0].second), true)   // non-vacuity: ungated, it is planned

  YieldGateConfig yg;
  yg.rank_by = YieldRank::Intensity;   // the walk this section pins: strongest first, the next two take the place
  yg.threshold = thresholdBetween(FOUR_FRAGS, order[0], ctx);
  Ms3PlanGate gate = yieldGate(yg, ctx);
  std::vector<Ms3Target> gated;
  std::string log;
  {
    EventCapture cap;
    gated = tracker.planNextScans(142, &gate);
    log = cap.str();
  }
  TEST_EQUAL(gated.size(), 2)   // the place is still filled ...
  TEST_EQUAL(planHas(gated, order[0].first, order[0].second), false)   // ... not by the refused fragment ...
  TEST_EQUAL(planHas(gated, order[1].first, order[1].second), true)    // ... but by the next two
  TEST_EQUAL(planHas(gated, order[2].first, order[2].second), true)
  const std::string refused = "[MS3-YIELD] gate=characterization precursor_id=142 ion=" + order[0].first
                              + std::to_string(order[0].second) + " charge=2 verdict=refused";
  TEST_EQUAL(countOf(log, "[MS3-YIELD]"), 1)
  TEST_EQUAL(countOf(log, refused), 1)
}
END_SECTION

// T11. Governed, a yield-refused fragment takes no cleavage slot: one slot (max_targets 1), the first fragment refused,
// the slot goes to the second.
// FAILS UNDER: the yield test placed after room.admit (the slot would be spent on the refused fragment: empty plan).
START_SECTION(governed_yield_refusal_takes_no_cleavage_slot)
{
  const ScanCommand ctx = yieldCtx();
  std::vector<std::pair<std::string, int>> order;
  {
    Config cfg{gateConfig("fragments", 10)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedGateWinner(tracker, 143, withIsoCos(FOUR_FRAGS, {}, 0.99));
    for (const Ms3Target& t : tracker.planNextScans(143)) order.emplace_back(t.ion_type, t.ion_index);
  }
  ABORT_IF(order.size() != 4)

  Config cfg{gateConfig("fragments", 10)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 143, withIsoCos(FOUR_FRAGS, {order[0]}, 0.85));
  Ms3RegistrySnapshot snap;
  snap.found = true;
  Ms3PlanGate slots_only = openGate();
  slots_only.snapshot = &snap;
  slots_only.max_targets = 1;
  const std::vector<Ms3Target> one_slot = tracker.planNextScans(143, &slots_only);
  TEST_EQUAL(one_slot.size(), 1)
  ABORT_IF(one_slot.size() != 1)
  TEST_EQUAL(one_slot[0].ion_type == order[0].first && one_slot[0].ion_index == order[0].second, true)   // non-vacuity

  YieldGateConfig yg;
  yg.rank_by = YieldRank::Intensity;   // the walk this section pins: strongest first, the second takes the slot
  yg.threshold = thresholdBetween(FOUR_FRAGS, order[0], ctx);
  Ms3PlanGate gate = yieldGate(yg, ctx);
  gate.snapshot = &snap;
  gate.max_targets = 1;
  const std::vector<Ms3Target> gated = tracker.planNextScans(143, &gate);
  TEST_EQUAL(gated.size(), 1)
  ABORT_IF(gated.size() != 1)
  TEST_EQUAL(gated[0].ion_type == order[1].first && gated[0].ion_index == order[1].second, true)
}
END_SECTION

// T12. The yield walk runs on the charges the gate's other tests admitted, strongest first: b6 at 2+ and 3+, a charge
// floor of 3, the yield admits 3+ and it anchors. Every charge refused: the fragment is refused, one line per charge
// the walk tested (single mode tests until one passes, so here both).
// FAILS UNDER: a walk over the unfiltered envelope (2+ would anchor), or a refusal that writes no line per charge.
START_SECTION(yield_charge_walk_anchors_first_passing_charge)
{
  const ScanCommand ctx = yieldCtx();
  Config cfg{gateConfig("fragments", 10)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 144, withIsoCos({{"b", 6, 700.0, {2, 3}}}, {}, 0.99));
  YieldGateConfig open_yield;
  open_yield.threshold = 1e-12;
  Ms3PlanGate gate = yieldGate(open_yield, ctx);
  gate.charge_floor = 3;
  const std::vector<Ms3Target> walked = tracker.planNextScans(144, &gate);
  TEST_EQUAL(walked.size(), 1)
  ABORT_IF(walked.size() != 1)
  TEST_EQUAL(walked[0].frag_charge, 3)

  YieldGateConfig closed;
  closed.threshold = 1e9;
  Ms3PlanGate shut = yieldGate(closed, ctx);
  std::string log;
  {
    EventCapture cap;
    TEST_EQUAL(tracker.planNextScans(144, &shut).size(), 0)
    log = cap.str();
  }
  TEST_EQUAL(countOf(log, "[MS3-YIELD] gate=characterization precursor_id=144 ion=b6 charge=2 verdict=refused"), 1)
  TEST_EQUAL(countOf(log, "[MS3-YIELD] gate=characterization precursor_id=144 ion=b6 charge=3 verdict=refused"), 1)
}
END_SECTION

// T13. Separate tests every charge: b6 at 2+, 3+ and 4+ is three scans when the yield admits, none and three refused
// lines when it refuses. Single mode, admitted, writes no line at all.
// FAILS UNDER: a separate walk that stops at the first admitted charge (one scan), or a refused charge still acquired.
START_SECTION(separate_acquires_only_passing_charges)
{
  const ScanCommand ctx = yieldCtx();
  Config cfg{gateConfig("fragments", 10, "separate")};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 145, withIsoCos({{"b", 6, 700.0, {2, 3, 4}}}, {}, 0.99));
  YieldGateConfig open_yield;
  open_yield.threshold = 1e-12;
  Ms3PlanGate open = yieldGate(open_yield, ctx);
  std::string log;
  {
    EventCapture cap;
    TEST_EQUAL(tracker.planNextScans(145, &open).size(), 3)
    log = cap.str();
  }
  TEST_EQUAL(countOf(log, "[MS3-YIELD]"), 0)   // admitted MS3 write no line: their E is the row's
  YieldGateConfig closed;
  closed.threshold = 1e9;
  Ms3PlanGate shut = yieldGate(closed, ctx);
  {
    EventCapture cap;
    TEST_EQUAL(tracker.planNextScans(145, &shut).size(), 0)
    log = cap.str();
  }
  TEST_EQUAL(countOf(log, "verdict=refused"), 3)
}
END_SECTION

// T14. Multiplexed notches are not yield-tested: the walk stops at the anchor (2+), and 3+ and 4+ still ride as
// notches, exactly as without the yield gate.
// FAILS UNDER: notches drawn from the yield-admitted charges only (the target would carry none).
START_SECTION(multiplexed_notches_are_untested)
{
  const ScanCommand ctx = yieldCtx();
  Config cfg{gateConfig("fragments", 10, "multiplexed")};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 146, withIsoCos({{"b", 6, 700.0, {2, 3, 4}}}, {}, 0.99));
  YieldGateConfig open_yield;
  open_yield.threshold = 1e-12;
  Ms3PlanGate gate = yieldGate(open_yield, ctx);
  const std::vector<Ms3Target> plan = tracker.planNextScans(146, &gate);
  TEST_EQUAL(plan.size(), 1)
  ABORT_IF(plan.size() != 1)
  TEST_EQUAL(plan[0].frag_charge, 2)
  std::set<int> notches;
  for (const NotchCandidate& n : plan[0].notches) notches.insert(n.charge);
  TEST_EQUAL(notches == std::set<int>({3, 4}), true)
}
END_SECTION

// T15. The ungoverned yield gate (no snapshot, no spectrum, the dispatch's charge floor, the budget as its room) moves
// the charge floor into the walk: today a 2+ fragment takes the one place and dispatch then skips it below the floor;
// with the gate the 3+ fragment takes the place.
// FAILS UNDER: an ungoverned gate that ignores the floor (the 2+ fragment would still take the place).
START_SECTION(ungoverned_gate_moves_floor_into_the_walk)
{
  const ScanCommand ctx = yieldCtx();
  Config cfg{gateConfig("fragments", 1)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 147, withIsoCos({{"b", 5, 550.0, {2}}, {"b", 6, 700.0, {3}}}, {}, 0.99));
  const std::vector<Ms3Target> today = tracker.planNextScans(147);
  TEST_EQUAL(today.size(), 1)
  YieldGateConfig open_yield;
  open_yield.threshold = 1e-12;
  Ms3PlanGate gate = yieldGate(open_yield, ctx);
  gate.charge_floor = 3;
  gate.max_targets = 1;
  const std::vector<Ms3Target> gated = tracker.planNextScans(147, &gate);
  TEST_EQUAL(gated.size(), 1)
  ABORT_IF(gated.size() != 1)
  TEST_EQUAL(gated[0].frag_charge, 3)
  TEST_EQUAL(planHas(gated, "b", 6), true)
}
END_SECTION

// T16. rank_by yield (the default) walks a plan by predicted bond yield, highest first, and sends it in that order;
// rank_by intensity keeps the strongest-first walk. Four fragments whose yield order is their intensity order reversed,
// budget 2: yield plans the two highest yields, intensity the two strongest. With a threshold between the second and
// third yield and room for all four, the plan is exactly the two above it and the two below are refused.
// FAILS UNDER: a sort left on intensity (the yield plan would be the two strongest), or a ranking that lets a
// below-threshold fragment in.
START_SECTION(yield_rank_takes_the_highest_yield_first)
{
  const ScanCommand ctx = yieldCtx();
  // The intensity order, learnt with every fragment at the same isotope cosine.
  std::vector<std::pair<std::string, int>> order;
  {
    Config cfg{gateConfig("fragments", 10)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedGateWinner(tracker, 146, withIsoCos(FOUR_FRAGS, {}, 0.99));
    for (const Ms3Target& t : tracker.planNextScans(146)) order.emplace_back(t.ion_type, t.ion_index);
  }
  TEST_EQUAL(order.size(), 4)
  ABORT_IF(order.size() != 4)
  // Isotope cosines rising against the intensity order: 0.05 a place outweighs the fragments' lengths and masses, so
  // the yield order is the intensity order reversed. e[i] is the yield of order[i].
  std::vector<GateFrag> frags = FOUR_FRAGS;
  std::vector<double> e(4, 0.0);
  for (GateFrag& g : frags)
    for (size_t i = 0; i < 4; ++i)
      if (g.ion_type == order[i].first && g.ion_index == order[i].second)
      {
        g.iso_cos = 0.85 + 0.05 * static_cast<double>(i);
        e[i] = yieldOf(g, ctx);
      }
  const bool rising = e[0] > 0.0 && e[0] < e[1] && e[1] < e[2] && e[2] < e[3];
  TEST_EQUAL(rising, true)   // non-vacuity: the two orders are opposite
  ABORT_IF(!rising)
  auto is = [](const Ms3Target& t, const std::pair<std::string, int>& k) { return t.ion_type == k.first && t.ion_index == k.second; };

  for (const std::string mode : {"fragments", "coverage"})
    for (const bool governed : {false, true})
    {
      Config cfg{gateConfig(mode, 2)};
      IdaLogger logger(cfg);
      ProteoformTracker tracker(cfg, logger);
      feedGateWinner(tracker, 146, frags);
      Ms3RegistrySnapshot snap;
      snap.found = true;

      YieldGateConfig by_yield;
      by_yield.threshold = e[0] / 2.0;   // below every yield: the budget decides
      Ms3PlanGate gate = yieldGate(by_yield, ctx);
      if (governed) gate.snapshot = &snap;
      const std::vector<Ms3Target> ranked = tracker.planNextScans(146, &gate);
      TEST_EQUAL(ranked.size(), 2)
      ABORT_IF(ranked.size() != 2)
      TEST_EQUAL(is(ranked[0], order[3]), true)   // the highest yield, sent first
      TEST_EQUAL(is(ranked[1], order[2]), true)

      YieldGateConfig by_intensity = by_yield;
      by_intensity.rank_by = YieldRank::Intensity;
      Ms3PlanGate gate_i = yieldGate(by_intensity, ctx);
      if (governed) gate_i.snapshot = &snap;
      const std::vector<Ms3Target> strongest = tracker.planNextScans(146, &gate_i);
      TEST_EQUAL(strongest.size(), 2)
      ABORT_IF(strongest.size() != 2)
      TEST_EQUAL(is(strongest[0], order[0]), true)
      TEST_EQUAL(is(strongest[1], order[1]), true)
    }

  // A threshold between the second and third yield, room for all four.
  {
    Config cfg{gateConfig("fragments", 10)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedGateWinner(tracker, 147, frags);
    YieldGateConfig yg;
    yg.threshold = std::sqrt(e[1] * e[2]);
    Ms3PlanGate gate = yieldGate(yg, ctx);
    std::vector<Ms3Target> plan;
    std::string log;
    {
      EventCapture cap;
      plan = tracker.planNextScans(147, &gate);
      log = cap.str();
    }
    TEST_EQUAL(plan.size(), 2)
    ABORT_IF(plan.size() != 2)
    TEST_EQUAL(is(plan[0], order[3]), true)
    TEST_EQUAL(is(plan[1], order[2]), true)
    TEST_EQUAL(countOf(log, "verdict=refused"), 2)
  }
}
END_SECTION

// T17. Ambiguity keeps its round-robin and orders each modification's containers by yield: one ambiguous modification,
// two containers, the weaker one at the higher yield, budget 1.
// FAILS UNDER: the per-modification list left on intensity under rank_by yield.
START_SECTION(yield_rank_ambiguity_orders_each_modification_by_yield)
{
  const ScanCommand ctx = yieldCtx();
  const std::vector<GateFrag> two_amb{{"b", 6, 700.0, {2}}, {"y", 6, 900.0, {2}}};
  std::vector<std::pair<std::string, int>> order;
  {
    Config cfg{gateConfig("ambiguity", 10)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedGateWinner(tracker, 148, withIsoCos(two_amb, {}, 0.99), /*ambiguous=*/true);
    for (const Ms3Target& t : tracker.planNextScans(148)) order.emplace_back(t.ion_type, t.ion_index);
  }
  TEST_EQUAL(order.size(), 2)
  ABORT_IF(order.size() != 2)
  std::vector<GateFrag> frags = two_amb;
  double e_strong = 0.0, e_weak = 0.0;
  for (GateFrag& g : frags)
  {
    const bool strong = g.ion_type == order[0].first && g.ion_index == order[0].second;
    g.iso_cos = strong ? 0.85 : 1.0;
    (strong ? e_strong : e_weak) = yieldOf(g, ctx);
  }
  TEST_EQUAL(e_strong > 0.0 && e_strong < e_weak, true)   // non-vacuity: the weaker container has the higher yield
  ABORT_IF(!(e_strong > 0.0 && e_strong < e_weak))

  Config cfg{gateConfig("ambiguity", 1)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedGateWinner(tracker, 148, frags, /*ambiguous=*/true);
  YieldGateConfig by_yield;
  by_yield.threshold = e_strong / 2.0;
  Ms3PlanGate gate = yieldGate(by_yield, ctx);
  const std::vector<Ms3Target> ranked = tracker.planNextScans(148, &gate);
  TEST_EQUAL(ranked.size(), 1)
  ABORT_IF(ranked.size() != 1)
  TEST_EQUAL(ranked[0].ion_type == order[1].first && ranked[0].ion_index == order[1].second, true)

  YieldGateConfig by_intensity = by_yield;
  by_intensity.rank_by = YieldRank::Intensity;
  Ms3PlanGate gate_i = yieldGate(by_intensity, ctx);
  const std::vector<Ms3Target> strongest = tracker.planNextScans(148, &gate_i);
  TEST_EQUAL(strongest.size(), 1)
  ABORT_IF(strongest.size() != 1)
  TEST_EQUAL(strongest[0].ion_type == order[0].first && strongest[0].ion_index == order[0].second, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0073. The interleaved send order: the rule alone, ProteoformTracker::interleavedOrder. Targets are indexed in
// ranked order (0 is sent first under "ranked"); each turn goes to the group that has had the fewest fragments, and a
// tie to the higher-ranked target.
/////////////////////////////////////////////////////////////

// Two groups of three and two targets, nothing sent.
// FAILS UNDER: a ranked send (0 1 2 3 4).
START_SECTION(interleaved_order_alternates_two_groups)
{
  const int n_term = 0, c_term = 1;
  const std::vector<std::vector<int>> groups_of{{c_term}, {c_term}, {c_term}, {n_term}, {n_term}};
  TEST_STRING_EQUAL(orderText(ProteoformTracker::interleavedOrder(groups_of, {0, 0})), "0 3 1 4 2")
}
END_SECTION

// A group with one target: once it is sent, the other group sends the rest in ranked order.
// FAILS UNDER: a turn kept for a group that has nothing left.
START_SECTION(interleaved_order_runs_out_of_a_group)
{
  const int n_term = 0, c_term = 1;
  const std::vector<std::vector<int>> groups_of{{c_term}, {c_term}, {c_term}, {n_term}};
  TEST_STRING_EQUAL(orderText(ProteoformTracker::interleavedOrder(groups_of, {0, 0})), "0 3 1 2")
}
END_SECTION

// One group: the interleaved order is the ranked order.
// FAILS UNDER: any reordering inside a group.
START_SECTION(interleaved_order_of_one_group_is_ranked)
{
  const std::vector<std::vector<int>> groups_of{{0}, {0}, {0}};
  TEST_STRING_EQUAL(orderText(ProteoformTracker::interleavedOrder(groups_of, std::vector<int>(1, 0))), "0 1 2")
}
END_SECTION

// A target in two groups counts for both: target 2 serves A and B, so after it A has had two and B one, the next
// turn is B's (3), and the tie at two each goes to the higher-ranked target (1).
// FAILS UNDER: one owner per target (0 2 1 3 4).
START_SECTION(interleaved_order_counts_a_shared_target_for_each_group)
{
  const int mod_a = 0, mod_b = 1;
  const std::vector<std::vector<int>> groups_of{{mod_a}, {mod_a}, {mod_a, mod_b}, {mod_b}, {mod_b}};
  TEST_STRING_EQUAL(orderText(ProteoformTracker::interleavedOrder(groups_of, {0, 0})), "0 2 3 1 4")
}
END_SECTION

// The counts start from what was already sent: the C-terminal group has had one, so the N-terminal one goes first.
// FAILS UNDER: counts that start at zero (0 3 1 4 2).
START_SECTION(interleaved_order_starts_from_served)
{
  const int n_term = 0, c_term = 1;
  const std::vector<std::vector<int>> groups_of{{c_term}, {c_term}, {c_term}, {n_term}, {n_term}};
  std::vector<int> served(2, 0);
  served[c_term] = 1;
  TEST_STRING_EQUAL(orderText(ProteoformTracker::interleavedOrder(groups_of, served)), "3 0 4 1 2")
}
END_SECTION

// A target that is already among the sent takes its turn and serves nobody a second time: after target 0 the
// C-terminal group has still had one, so its next target (1) goes before the N-terminal group's second (4).
// FAILS UNDER: a re-acquisition counted twice (3 0 4 1 2).
START_SECTION(interleaved_order_does_not_count_a_counted_target)
{
  const int n_term = 0, c_term = 1;
  const std::vector<std::vector<int>> groups_of{{c_term}, {c_term}, {c_term}, {n_term}, {n_term}};
  std::vector<int> served(2, 0);
  served[c_term] = 1;
  const std::vector<bool> counted{true, false, false, false, false};
  TEST_STRING_EQUAL(orderText(ProteoformTracker::interleavedOrder(groups_of, served, counted)), "3 0 1 4 2")
}
END_SECTION

// A target in no group goes after every grouped one, however high it ranks.
// FAILS UNDER: a groupless target treated as the least served (0 1 2).
START_SECTION(interleaved_order_puts_a_groupless_target_last)
{
  const int n_term = 0, c_term = 1;
  const std::vector<std::vector<int>> groups_of{std::vector<int>(), {c_term}, {n_term}};
  TEST_STRING_EQUAL(orderText(ProteoformTracker::interleavedOrder(groups_of, {0, 0})), "1 2 0")
}
END_SECTION

// Three groups: every group has had one before any has two.
// FAILS UNDER: a two-group shortcut.
START_SECTION(interleaved_order_of_three_groups)
{
  const int mod_a = 0, mod_b = 1, mod_c = 2;
  const std::vector<std::vector<int>> groups_of{{mod_a}, {mod_a}, {mod_b}, {mod_c}};
  TEST_STRING_EQUAL(orderText(ProteoformTracker::interleavedOrder(groups_of, {0, 0, 0})), "0 2 3 1")
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0073. The interleaved send order in the planner, on PEPTIDEK (L = 8). Every section plans one fixture twice,
// characterization.order absent ("ranked") and "interleaved", behind a yield gate every fragment clears; rankedFrags
// makes the ranked order the list order, which each section asserts first. The two plans hold the same targets
// (decision 1); only the order differs.
/////////////////////////////////////////////////////////////

// Three C-terminal fragments ranked above two N-terminal ones: ranked sends a terminus out before the other starts,
// interleaved takes turns from the first fragment on. Under fragments and under coverage (each of the five holds an
// unwitnessed bond: 1, 5 and 7 are witnessed by no fragment).
// FAILS UNDER: planNextScans without the interleave call; groups that are not the two termini.
START_SECTION(interleaved_alternates_the_termini)
{
  const ScanCommand ctx = yieldCtx();
  const std::vector<GateFrag> frags = rankedFrags({{"y", 6, 0.0, {2}}, {"y", 5, 0.0, {2}}, {"y", 4, 0.0, {2}}, {"b", 3, 0.0, {2}}, {"b", 6, 0.0, {2}}});
  TEST_EQUAL(yieldsFall(frags, ctx), true)
  ABORT_IF(!yieldsFall(frags, ctx))
  for (const std::string mode : {"fragments", "coverage"})
  {
    SendOrderCase fixture;
    fixture.mode = mode;
    const std::vector<Ms3Target> ranked = planInSendOrder(fixture, "", 150, frags, ctx);
    const std::vector<Ms3Target> interleaved = planInSendOrder(fixture, "interleaved", 150, frags, ctx);
    TEST_STRING_EQUAL(ionsText(ranked), "y6 y5 y4 b3 b6")   // precondition: the imposed rank
    ABORT_IF(ionsText(ranked) != "y6 y5 y4 b3 b6")
    TEST_EQUAL(keysOf(ranked) == keysOf(interleaved), true)   // the same plan
    TEST_STRING_EQUAL(ionsText(interleaved), "y6 b3 y5 b6 y4")
  }
}
END_SECTION

// The order never changes what is planned: a budget of three takes the three highest-ranked fragments, all
// C-terminal, under both orders, and one group is sent in ranked order.
// FAILS UNDER: interleaving the choice instead of the send (b3 would take y4's place).
START_SECTION(interleaved_never_changes_the_plan)
{
  const ScanCommand ctx = yieldCtx();
  const std::vector<GateFrag> frags = rankedFrags({{"y", 6, 0.0, {2}}, {"y", 5, 0.0, {2}}, {"y", 4, 0.0, {2}}, {"b", 3, 0.0, {2}}, {"b", 6, 0.0, {2}}});
  TEST_EQUAL(yieldsFall(frags, ctx), true)
  ABORT_IF(!yieldsFall(frags, ctx))
  SendOrderCase fixture;
  fixture.mode = "fragments";
  fixture.budget = 3;
  const std::vector<Ms3Target> ranked = planInSendOrder(fixture, "", 151, frags, ctx);
  const std::vector<Ms3Target> interleaved = planInSendOrder(fixture, "interleaved", 151, frags, ctx);
  TEST_STRING_EQUAL(ionsText(ranked), "y6 y5 y4")   // precondition: the imposed rank, cut at the budget
  ABORT_IF(ionsText(ranked) != "y6 y5 y4")
  TEST_EQUAL(keysOf(ranked) == keysOf(interleaved), true)
  TEST_STRING_EQUAL(ionsText(interleaved), "y6 y5 y4")
}
END_SECTION

// Ambiguity: modification A on residues 2..3, B on 5..6. b3 and b4 contain A, y4 and y5 contain B, b7 contains both
// and counts for both; no fragment brackets either range, so nothing narrows before the plan. After b3 and b7, A has
// had two and B one: y4 is next, then the tie at two each goes to the higher-ranked b4.
// FAILS UNDER: one owner per container, or the order the round-robin chose the containers in (both: b3 b7 b4 y4 y5).
START_SECTION(interleaved_counts_a_container_of_two_modifications_for_both)
{
  const ScanCommand ctx = yieldCtx();
  const std::vector<GateFrag> frags = rankedFrags({{"b", 3, 0.0, {2}}, {"b", 4, 0.0, {2}}, {"b", 7, 0.0, {2}}, {"y", 4, 0.0, {2}}, {"y", 5, 0.0, {2}}});
  TEST_EQUAL(yieldsFall(frags, ctx), true)
  ABORT_IF(!yieldsFall(frags, ctx))
  SendOrderCase fixture;
  fixture.mode = "ambiguity";
  fixture.sites = {{2, 3}, {5, 6}};
  const std::vector<Ms3Target> ranked = planInSendOrder(fixture, "", 152, frags, ctx);
  const std::vector<Ms3Target> interleaved = planInSendOrder(fixture, "interleaved", 152, frags, ctx);
  TEST_STRING_EQUAL(ionsText(ranked), "b3 b4 b7 y4 y5")   // precondition: the imposed rank
  ABORT_IF(ionsText(ranked) != "b3 b4 b7 y4 y5")
  TEST_EQUAL(keysOf(ranked) == keysOf(interleaved), true)
  TEST_STRING_EQUAL(ionsText(interleaved), "b3 b7 y4 b4 y5")
}
END_SECTION

// A governed plan (per_proteoform) continues from the fragments the proteoform has sent. Nothing sent: the turns
// start level. y7 sent (not a fragment of this plan): the C-terminus has had one, so b3 goes first. y6 sent and
// planned again (a sent intensity of 0 admits the re-acquisition): b3 goes first, y6 takes the C-terminal turn and
// is not a second serving, so y5 follows it before b6.
// FAILS UNDER: counts that start at zero under a snapshot; a re-acquisition counted twice (b3 y6 b6 y5 y4).
START_SECTION(interleaved_continues_from_the_committed_fragments)
{
  const ScanCommand ctx = yieldCtx();
  const std::vector<GateFrag> frags = rankedFrags({{"y", 6, 0.0, {2}}, {"y", 5, 0.0, {2}}, {"y", 4, 0.0, {2}}, {"b", 3, 0.0, {2}}, {"b", 6, 0.0, {2}}});
  TEST_EQUAL(yieldsFall(frags, ctx), true)
  ABORT_IF(!yieldsFall(frags, ctx))
  struct Committed
  {
    std::vector<Ms3Cleavage> sent;
    std::string interleaved;
  };
  std::vector<Committed> cases(3);
  cases[0].interleaved = "y6 b3 y5 b6 y4";
  cases[1].sent.push_back(Ms3Cleavage{"y", 7});
  cases[1].interleaved = "b3 y6 b6 y5 y4";
  cases[2].sent.push_back(Ms3Cleavage{"y", 6});
  cases[2].interleaved = "b3 y6 y5 b6 y4";
  for (const Committed& committed : cases)
  {
    Ms3RegistrySnapshot snap;
    snap.found = true;
    for (const Ms3Cleavage& cleavage : committed.sent) snap.sent[cleavage][2] = 0.0;
    SendOrderCase fixture;
    fixture.mode = "fragments";
    fixture.snapshot = &snap;
    const std::vector<Ms3Target> ranked = planInSendOrder(fixture, "", 153, frags, ctx);
    const std::vector<Ms3Target> interleaved = planInSendOrder(fixture, "interleaved", 153, frags, ctx);
    TEST_STRING_EQUAL(ionsText(ranked), "y6 y5 y4 b3 b6")   // precondition: the imposed rank, whatever was sent
    ABORT_IF(ionsText(ranked) != "y6 y5 y4 b3 b6")
    TEST_EQUAL(keysOf(ranked) == keysOf(interleaved), true)
    TEST_STRING_EQUAL(ionsText(interleaved), committed.interleaved)
  }
}
END_SECTION

// Governed ambiguity: a sent fragment counts for the modifications its cover contains, read from its key (b7 and y6
// are not fragments of this identification). A on 2..3, B on 5..6; the plan is y5 (B), b3 (A), b4 (A), y4 (B).
// Nothing sent, b7 sent (one each for A and B) or a raw mass sent (no group): the turns start level and the order is
// the ranked one. y6 sent (it contains B only): A is behind, so b3 goes first.
// FAILS UNDER: a sent fragment's groups looked up in the model's fragments (y6 would count for nothing: y5 b3 b4 y4);
// a sent container counted for one modification, or a raw mass counted (both: y5 b3 y4 b4).
START_SECTION(interleaved_ambiguity_counts_committed_fragments_by_what_they_contain)
{
  const ScanCommand ctx = yieldCtx();
  const std::vector<GateFrag> frags = rankedFrags({{"y", 5, 0.0, {2}}, {"b", 3, 0.0, {2}}, {"b", 4, 0.0, {2}}, {"y", 4, 0.0, {2}}});
  TEST_EQUAL(yieldsFall(frags, ctx), true)
  ABORT_IF(!yieldsFall(frags, ctx))
  struct Committed
  {
    std::vector<Ms3Cleavage> sent;
    std::string interleaved;
  };
  std::vector<Committed> cases(4);
  cases[0].interleaved = "y5 b3 b4 y4";
  cases[1].sent.push_back(Ms3Cleavage{"b", 7});
  cases[1].interleaved = "y5 b3 b4 y4";
  cases[2].sent.push_back(Ms3Cleavage{"u", 1234});
  cases[2].interleaved = "y5 b3 b4 y4";
  cases[3].sent.push_back(Ms3Cleavage{"y", 6});
  cases[3].interleaved = "b3 y5 b4 y4";
  for (const Committed& committed : cases)
  {
    Ms3RegistrySnapshot snap;
    snap.found = true;
    for (const Ms3Cleavage& cleavage : committed.sent) snap.sent[cleavage][2] = 0.0;
    SendOrderCase fixture;
    fixture.mode = "ambiguity";
    fixture.sites = {{2, 3}, {5, 6}};
    fixture.snapshot = &snap;
    const std::vector<Ms3Target> ranked = planInSendOrder(fixture, "", 154, frags, ctx);
    const std::vector<Ms3Target> interleaved = planInSendOrder(fixture, "interleaved", 154, frags, ctx);
    TEST_STRING_EQUAL(ionsText(ranked), "y5 b3 b4 y4")   // precondition: the imposed rank, whatever was sent
    ABORT_IF(ionsText(ranked) != "y5 b3 b4 y4")
    TEST_EQUAL(keysOf(ranked) == keysOf(interleaved), true)
    TEST_STRING_EQUAL(ionsText(interleaved), committed.interleaved)
  }
}
END_SECTION

// fragment_charges separate: a fragment takes one turn and its charges stay together. y6 at 2+ and 3+, y5 at 2+, b3
// at 2+ and 3+: after the two scans of y6 the C-terminus has had one fragment, so both scans of b3 go before y5.
// FAILS UNDER: every scan taking a turn (y6 b3 y6 b3 y5), or a fragment counted once per charge.
START_SECTION(interleaved_keeps_a_fragments_charges_together)
{
  const ScanCommand ctx = yieldCtx();
  const std::vector<GateFrag> frags = rankedFrags({{"y", 6, 0.0, {2, 3}}, {"y", 5, 0.0, {2}}, {"b", 3, 0.0, {2, 3}}});
  TEST_EQUAL(yieldsFall(frags, ctx), true)
  ABORT_IF(!yieldsFall(frags, ctx))
  SendOrderCase fixture;
  fixture.mode = "fragments";
  fixture.fragment_charges = "separate";
  const std::vector<Ms3Target> ranked = planInSendOrder(fixture, "", 155, frags, ctx);
  const std::vector<Ms3Target> interleaved = planInSendOrder(fixture, "interleaved", 155, frags, ctx);
  TEST_STRING_EQUAL(ionsText(ranked), "y6 y6 y5 b3 b3")   // precondition: the imposed rank, a scan per charge
  ABORT_IF(ionsText(ranked) != "y6 y6 y5 b3 b3")
  TEST_EQUAL(keysOf(ranked) == keysOf(interleaved), true)
  TEST_STRING_EQUAL(ionsText(interleaved), "y6 y6 b3 b3 y5")
}
END_SECTION

END_TEST
