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

#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
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
  // last column (trigger_scan_id) is always non-empty, so every data row yields all 19 fields; interior
  // empty fields (localized_mods/ambiguous_mods) sit between tabs and ARE preserved as "".
  std::vector<std::string> splitTab(const std::string& line)
  {
    std::vector<std::string> fields;
    std::istringstream ls(line);
    std::string fld;
    while (std::getline(ls, fld, '\t')) fields.push_back(fld);
    return fields;
  }
}

START_TEST(ProteoformTracker_Trajectory, "$Id$")

/////////////////////////////////////////////////////////////
// Per-precursor trajectory: MS2 baseline row + additive (non-resetting) MS3 fold rows.
/////////////////////////////////////////////////////////////
START_SECTION(ms2_baseline_then_accumulating_ms3_folds)
{
  const int precursor_id = 7;

  // Pooled-log path embedded in the Config JSON so the engine opens the pooled stream (Config.cpp:424
  // reads runtime.pooled_identification_log_path). Use a RELATIVE filename in the test CWD (NOT an
  // absolute temp path): an absolute Windows temp path contains backslashes that would form invalid JSON
  // escapes when embedded below. The pooled stream is opened in APPEND mode, so remove any stale file
  // first. (Mirrors the FLASHIda_LoggingFields_test relative-filename + std::remove pattern.)
  const std::string pooled_path = "pt_trajectory_pooled.tsv";
  std::remove(pooled_path.c_str());

  std::string cfg_json = std::string(R"({
    "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 1, "max_charge": 50, "min_mass": 100, "max_mass": 50000, "tol": [10, 10, 10] },
    "precursor_selection": { "RT_window": 180, "target_mode": 0, "AllCharges": false, "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1 },
    "flashtnt": { "min_length": 3, "max_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": [ { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 } ],
      "ms3": [ { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25, "resolution": 120000 } ]
    },
    "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 } },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "characterization": { "objective": "coverage", "protein_sequence": "PEPTIDEK" },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "intensity", "max_targets": 10 },
      "ms3": { "selection": "intensity", "max_targets": 10 }
    },
    "runtime": { "pooled_identification_log_path": ")") + std::string(pooled_path) + R"(" } })";

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
  // order — a reorder in IdaLogger's pooled writer needs no change here. The pooled stream has 19 columns.
  // ---------------------------------------------------------------------------------------------
  std::ifstream f(pooled_path.c_str());
  TEST_TRUE(f.good())

  std::string header;
  TEST_TRUE(static_cast<bool>(std::getline(f, header)))
  const std::vector<std::string> header_cols = splitTab(header);
  TEST_EQUAL((int)header_cols.size(), 19)   // 14 base + grouped fragment-mass table (masses|ions|measured|theoretical|diff_da|diff_ppm)

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

  // Every data row must carry all 19 fields (the trailing trigger_scan_id is non-empty here).
  for (const std::vector<std::string>& r : rows) TEST_EQUAL((int)r.size(), 19)

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

END_TEST
