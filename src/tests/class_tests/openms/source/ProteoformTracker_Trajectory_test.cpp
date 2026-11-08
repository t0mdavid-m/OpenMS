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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
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

  // A single-fragment MS3 ProteoformMatch. mapScanOntoModel_'s ms_level == 3 branch
  // (ProteoformTracker.cpp:558-563) reads the MS3 fields equiv_type / equiv_index / adjusted_mass (NOT
  // ion_type / ion_index / observed_mass), so we set ONLY those. score = -1 and proteoform_sequence = ""
  // because foldMs3 never re-picks a winner: the model already carries the winner sequence, and
  // mapScanOntoModel_ uses the MODEL's sequence (m.proteoform_sequence), not the match's.
  FragmentAnalysis::ProteoformMatch makeMatchMs3(const std::string& equiv_type, int equiv_index, double adjusted_mass)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.score = -1.0;              // not a winner candidate (foldMs3 folds onto the existing baseline)
    m.region_start = -1;         // full-protein equiv index maps straight into the winner region
    m.region_end = -1;
    m.proteoform_sequence = "";  // model already holds WINNER_SEQ; mapScanOntoModel_ uses the model's
    FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
    fm.equiv_type = equiv_type;       // MS3: full-protein equivalent ion type ("b"/"y")
    fm.equiv_index = equiv_index;     // MS3: full-protein equivalent ion index
    fm.adjusted_mass = adjusted_mass; // MS3: offset-adjusted-to-full-protein mass
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

  // ---- ADR-0064: the per-proteoform plan gate (Ms3PlanGate) ---------------------------------------------------
  //
  // A hand-built PeakGroup has getChargeIntensity() == 0, so every observation below has intensity 0 and a window
  // SNR of 0 (FragmentAnalysis::windowSnr returns 0 for no signal). The sections that need a charge ADMITTED pass no
  // spectrum (no SNR test) or snr_min -1; the SNR section passes a spectrum and moves only snr_min. The 2x rule is
  // driven through the snapshot: a sent intensity of 0 admits a re-acquisition at 0, a sent intensity of 1 refuses it.

  /// The tracker config with characterization's objective and budget as parameters (no log folder: nothing is read).
  std::string gateConfig(const std::string& mode, int max_targets)
  {
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
      "characterization": { "mode": ")JSON") + mode + "\", \"max_targets\": " + std::to_string(max_targets) + " } } }";
  }

  const double PROTON_MASS = 1.00728;
  const double ISOTOPE_STEP = 1.00235;

  /// A fragment PeakGroup seen at every charge in @p charges with TWO isotopes each, so PeakRecord::by_charge holds a
  /// real (non-zero-width) record per charge and ms2_by_charge carries all of them. The representative charge is the
  /// lowest (getMaxIntensityAbsCharge with no intensity).
  PeakGroup makeEnvelope(double mono_mass, const std::vector<int>& charges)
  {
    const int zmin = *std::min_element(charges.begin(), charges.end());
    const int zmax = *std::max_element(charges.begin(), charges.end());
    PeakGroup pg(zmin, zmax, true);
    pg.setMonoisotopicMass(mono_mass);
    for (int z : charges)
    {
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
  };

  /// Feed ONE winner MS2 scan carrying @p frags (the winner match claims each, pooled verbatim) and finalize. With
  /// @p ambiguous the match carries one modification on residues 4..5, which b5, b6 and y6 contain and nothing brackets.
  void feedGateWinner(ProteoformTracker& tracker, int precursor_id, const std::vector<GateFrag>& frags, bool ambiguous = false)
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
      if (g.charges.size() == 1) d.push_back(makeSyntheticPeakGroup(g.mass / g.charges[0] + PROTON_MASS, g.mass, g.charges[0]));
      else d.push_back(makeEnvelope(g.mass, g.charges));
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
    tracker.feedScan(precursor_id, 3, ms3_params, 201, d201, makeMatchMs3("b", 3, 350.0), -1.0, ms2_ctx);
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
    tracker.feedScan(precursor_id, 3, ms3_params, 202, d202, makeMatchMs3("y", 3, 450.0), -1.0, ms2_ctx);
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
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0064 G1. A null gate is today's plan, and a gate that admits everything changes nothing, under each of the
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
// ADR-0064 G2. A refused fragment hands its slot to the next candidate in the same order -- under fragments, and
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
// ADR-0064 G3. A charge whose window SNR is not above snr_min is refused. Every fixture charge here carries no
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
// ADR-0064 G4. Single mode, the strongest charge refused: the next admitted charge of the same fragment stands in
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
// ADR-0064 G5. The cleavage cap counts SENT cleavages, not re-acquisitions: with b6 sent and max_targets 2 the
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

END_TEST
