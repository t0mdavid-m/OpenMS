// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Localization pin: an ambiguous PTM that a NON-bracketing MS2 fragment cannot localize IS narrowed
// by a bracketing MS3 fragment carrying the propagated localization verdict FragmentMatch::includes_ptm.
//
//   ProteoformTracker::finalizeMS2(pid)   seeds modifications from the winner match's ptm_sites, then
//                                      runs narrowModifications_ (Pass A: MS2 mass-test).
//   ProteoformTracker::foldMs3(pid,..) additively maps the staged MS3 scan, then re-runs
//                                      narrowModifications_ (Pass B: MS3 includes_ptm verdict).
//
// mapScanOntoModel_ (ProteoformTracker.cpp:641) copies the MS3 FragmentMatch::includes_ptm onto the
// FragmentObservation. narrowModifications_ Pass B (ProteoformTracker.cpp:788-807) consumes ONLY
// obs.includes_ptm + the fragment's winner-region cover_start/cover_end -- it does NOT re-derive the
// verdict from the folded adjusted_mass. So this test needs no real theoretical masses: adjusted_mass
// is a dummy; the narrowing is driven purely by includes_ptm and the bracketing geometry.
//
//   includes_ptm == true  => the PTM lies INSIDE this fragment's coverage -> tighten TOWARD it:
//                            prefix b_k -> tightenUpper(cover_end); suffix y_k -> tightenLower(cover_start).
//
// This is the observable contract of the propagate+consume change: BEFORE the fix the MS3 verdict was
// dropped / re-derived from mass, so the ambiguous range never narrowed from MS3 evidence and stayed
// [15,26]; AFTER the fix a bracketing prefix tightens it to [15,20] and a bracketing suffix to [20,26].
//
// CONSTRUCTION NOTES:
//  - region_start/region_end == -1 => winner region == full protein (ws == 0), so the full-protein mod
//    coordinates equal the winner-region frame: a mod seeded [15,26] stays [15,26] and narrows in place.
//  - The MS2 seed fragment (prefix b8, cover_end == 8) does NOT bracket [15,26] (prefix brackets iff
//    rs <= cover_end < re, i.e. 15 <= 8 -> false), so finalize leaves the mod ambiguous -> the MS3
//    fold is the ONLY thing that can narrow it (asserted).
//  - A hand-built synthetic PeakGroup has getChargeIntensity() == 0, so every observation's intensity is
//    0; mapScanOntoModel_ never drops a matched fragment for lack of a peak (closest-mass fallback), so
//    each fed fragment still lands in m.fragments. Pass B uses the (0) intensity only as tighten support,
//    and the first inward move on a fresh boundary (support 0) always applies.

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

#include <cstring>
#include <string>

using namespace OpenMS;

namespace
{
  // A 32-residue winner proteoform (L = 32), long enough to seed an ambiguous mod at [15,26] with room
  // for bracketing prefix (b20) and suffix (cover_start 20) cleavages on either side of the range.
  const char* WINNER_SEQ = "PEPTIDEKPEPTIDEKPEPTIDEKPEPTIDEK";  // 32 residues
  const int   SEQ_LEN = 32;

  // Minimal synthetic PeakGroup (getChargeIntensity == 0); the fed fragment still maps (closest-mass
  // fallback in mapScanOntoModel_). Mirrors ProteoformTracker_Trajectory_test::makeSyntheticPeakGroup.
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

  // A single-fragment MS2 ProteoformMatch for the winner proteoform, carrying one AMBIGUOUS PTM site
  // [start,end] (start < end) and one region-relative MS2 fragment (ion_type/ion_index/observed_mass).
  // score = 1.0 + non-empty sequence => finalize accepts it as the winner and seeds modifications from
  // ptm_sites. region_start/end = -1 => identity-frame mapping (ws == 0).
  FragmentAnalysis::ProteoformMatch makeMatchMs2WithMod(const std::string& ion_type, int ion_index, double observed_mass,
                                                        int ptm_start, int ptm_end, double mass_shift)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.score = 1.0;
    m.region_start = -1;
    m.region_end = -1;
    m.matched_protein = "synthetic";
    m.proteoform_sequence = WINNER_SEQ;

    FragmentAnalysis::PTMSite site;
    site.start_position = ptm_start;     // ambiguous when start < end
    site.end_position = ptm_end;
    site.position = (ptm_start + ptm_end) / 2;
    site.mass_shift = mass_shift;
    m.ptm_sites.push_back(site);

    FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
    fm.ion_type = ion_type;      // MS2: region-relative ion type/index
    fm.ion_index = ion_index;
    fm.observed_mass = observed_mass;
    m.fragments.push_back(fm);
    return m;
  }

  // A single-fragment MS3 ProteoformMatch. mapScanOntoModel_'s ms_level == 3 branch reads equiv_type /
  // equiv_index / adjusted_mass and copies includes_ptm onto the FragmentObservation. adjusted_mass is a
  // dummy: Pass B uses includes_ptm + cover_* only, never the mass. score = -1 (foldMs3 never re-picks).
  FragmentAnalysis::ProteoformMatch makeMatchMs3(const std::string& equiv_type, int equiv_index,
                                                 double adjusted_mass, bool includes_ptm,
                                                 bool is_complement_flip = false)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.score = -1.0;
    m.region_start = -1;
    m.region_end = -1;
    m.proteoform_sequence = "";  // model already holds WINNER_SEQ; mapScanOntoModel_ uses the model's

    FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
    fm.equiv_type = equiv_type;        // MS3: full-protein equivalent ion type ("b"/"y")
    fm.equiv_index = equiv_index;      // MS3: full-protein equivalent ion index
    fm.adjusted_mass = adjusted_mass;  // MS3: dummy (Pass B ignores the mass)
    fm.includes_ptm = includes_ptm;    // MS3: propagated localization verdict consumed by Pass B
    fm.is_complement_flip = is_complement_flip;  // MS3: complement-flip flag consumed by Pass B's mod-aware verdict
    m.fragments.push_back(fm);
    return m;
  }

  // A minimal 2-stage MS2 context command (only used as the feedScan ms2_ctx argument).
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

  // Config JSON with the winner sequence as the characterization protein_sequence. No pooled log path
  // (this test reads the model directly via tracker.getModel(pid), not the pooled file).
  std::string makeConfigJson()
  {
    return std::string(R"({
    "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 1, "max_charge": 50, "min_mass": 100, "max_mass": 50000, "tol": [10, 10, 10] },
    "precursor_selection": { "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 3 },
    "flashtnt": { "min_length": 3, "max_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
    "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
    "faims": { "cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
      "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25, "resolution": 120000 }
    },
    "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 } },
    "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
    "characterization": { "mode": "coverage", "max_targets": 10, "protein_sequence": ")") + std::string(WINNER_SEQ) + R"(" },
    "conditional_ms2": false,
    "runtime": {} })";
  }
}

START_TEST(ProteoformTracker_Localization, "$Id$")

/////////////////////////////////////////////////////////////
// Section 1: a bracketing PREFIX MS3 fragment (includes_ptm=true) tightens the UPPER boundary.
//            Ambiguous [15,26] -- non-bracketing MS2 b8 leaves it [15,26] -- MS3 b20 -> [15,20].
/////////////////////////////////////////////////////////////
START_SECTION(ms3_prefix_verdict_localizes_upper)
{
  const int precursor_id = 11;

  Config cfg{makeConfigJson()};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  ScanCommand ms2_ctx = makeMs2Ctx();

  Ms2Params p;
  p.collision_energy = 29.0;
  p.activation_type = "HCD";
  p.reaction_time = 0.0;

  // MS2 winner: ambiguous mod [15,26] (mass_shift 42.0) + one NON-bracketing prefix fragment b8
  // (cover_end == 8; prefix brackets iff 15 <= 8 < 26 -> false), so finalize leaves the mod [15,26].
  DeconvolvedSpectrum d101(101);
  d101.push_back(makeSyntheticPeakGroup(700.0 / 2.0 + 1.0, 700.0, 2));
  tracker.feedScan(precursor_id, 2, p, 101, d101, makeMatchMs2WithMod("b", 8, 700.0, 15, 26, 42.0), 1.0, ms2_ctx);

  tracker.finalizeMS2(precursor_id);

  // After finalize: winner identified, ONE mod, still ambiguous [15,26] (MS2 seed did not bracket it).
  const ProteoformModel* mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL(mdl->proteoform_sequence, std::string(WINNER_SEQ))
  TEST_EQUAL((int)mdl->modifications.size(), 1)
  ABORT_IF(mdl->modifications.size() != 1)
  TEST_EQUAL(mdl->modifications[0].candidate_start, 15)
  TEST_EQUAL(mdl->modifications[0].candidate_end, 26)   // unlocalized before MS3 evidence

  // MS3 fold: prefix b20 (equiv_index 20 -> winner_idx 20 -> cover_end 20; brackets [15,26] since
  // 15 <= 20 < 26) with includes_ptm=true -> Pass B tightenUpper(cover_end=20) -> [15,20].
  Ms2Params ms3_params;
  ms3_params.collision_energy = 29.0;
  ms3_params.activation_type = "HCD";
  ms3_params.reaction_time = 0.0;

  DeconvolvedSpectrum d201(201);
  d201.push_back(makeSyntheticPeakGroup(2000.0 / 2.0 + 1.0, 2000.0, 2));
  tracker.feedScan(precursor_id, 3, ms3_params, 201, d201, makeMatchMs3("b", 20, 2000.0, /*includes_ptm=*/true), -1.0, ms2_ctx);
  tracker.foldMs3(precursor_id, "b20", "AAB");

  // After fold: the prefix verdict narrowed the UPPER boundary 26 -> 20. RED before the fix (the MS3
  // verdict was dropped / re-derived from the folded mass, so the range stayed [15,26]); GREEN after.
  mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL((int)mdl->modifications.size(), 1)
  ABORT_IF(mdl->modifications.size() != 1)
  TEST_EQUAL(mdl->modifications[0].candidate_start, 15)
  TEST_EQUAL(mdl->modifications[0].candidate_end, 20)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Section 2: a bracketing SUFFIX MS3 fragment (includes_ptm=true) tightens the LOWER boundary.
//            Ambiguous [15,26] -- non-bracketing MS2 b8 leaves it [15,26] -- MS3 suffix -> [20,26].
//
// Suffix equiv_index derivation from mapScanOntoModel_ (ProteoformTracker.cpp:558-628, region -1 => ws=0):
//   MS3: winner_idx = equiv_index - (P - we) = equiv_index - (P - P) = equiv_index   (ws=0, we=P=L)
//   suffix cover_start = L - winner_idx + 1  = L - equiv_index + 1
//   Want cover_start == 20 with L == 32  ->  equiv_index = L - 20 + 1 = 32 - 20 + 1 = 13.
//   Check: winner_idx = 13; cover_start = 32 - 13 + 1 = 20; cover_end = L = 32.
//   Brackets [15,26]? suffix brackets iff rs < cover_start <= re -> 15 < 20 <= 26 -> true.
/////////////////////////////////////////////////////////////
START_SECTION(ms3_suffix_verdict_localizes_lower)
{
  const int precursor_id = 12;

  Config cfg{makeConfigJson()};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  ScanCommand ms2_ctx = makeMs2Ctx();

  Ms2Params p;
  p.collision_energy = 29.0;
  p.activation_type = "HCD";
  p.reaction_time = 0.0;

  // Same ambiguous mod [15,26] + same non-bracketing MS2 seed b8 -> finalize leaves it [15,26].
  DeconvolvedSpectrum d101(101);
  d101.push_back(makeSyntheticPeakGroup(700.0 / 2.0 + 1.0, 700.0, 2));
  tracker.feedScan(precursor_id, 2, p, 101, d101, makeMatchMs2WithMod("b", 8, 700.0, 15, 26, 42.0), 1.0, ms2_ctx);

  tracker.finalizeMS2(precursor_id);

  const ProteoformModel* mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL((int)mdl->modifications.size(), 1)
  ABORT_IF(mdl->modifications.size() != 1)
  TEST_EQUAL(mdl->modifications[0].candidate_start, 15)
  TEST_EQUAL(mdl->modifications[0].candidate_end, 26)   // unlocalized before MS3 evidence

  // MS3 fold: suffix y with equiv_index = 13 -> cover_start 20 (brackets [15,26] since 15 < 20 <= 26),
  // includes_ptm=true -> Pass B tightenLower(cover_start=20) -> [20,26].
  Ms2Params ms3_params;
  ms3_params.collision_energy = 29.0;
  ms3_params.activation_type = "HCD";
  ms3_params.reaction_time = 0.0;

  const int equiv_index = SEQ_LEN - 20 + 1;   // = 13; makes winner-region cover_start == 20

  DeconvolvedSpectrum d201(201);
  d201.push_back(makeSyntheticPeakGroup(2100.0 / 2.0 + 1.0, 2100.0, 2));
  tracker.feedScan(precursor_id, 3, ms3_params, 201, d201, makeMatchMs3("y", equiv_index, 2100.0, /*includes_ptm=*/true), -1.0, ms2_ctx);
  tracker.foldMs3(precursor_id, "y13", "AAC");

  // After fold: the suffix verdict narrowed the LOWER boundary 15 -> 20. [15,26] -> [20,26].
  mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL((int)mdl->modifications.size(), 1)
  ABORT_IF(mdl->modifications.size() != 1)
  TEST_EQUAL(mdl->modifications[0].candidate_start, 20)
  TEST_EQUAL(mdl->modifications[0].candidate_end, 26)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Section 3 (Leg 3): the pooled localization narrows ONE-DIRECTIONALLY across successive MS3 folds. Each
// fold's [candidate_start, candidate_end] is NESTED within the prior (never widens/flips outward), and once a
// mod localizes it stays localized. This is the cumulative-narrowing invariant the leaf<->pooled seam depends
// on (the leaf is per-scan and may be wider; the pooled is monotone across scans). Injected fragments ->
// fully deterministic, no deconvolution.
/////////////////////////////////////////////////////////////
START_SECTION(ms3_pooled_localization_is_monotonic)
{
  const int precursor_id = 13;

  Config cfg{makeConfigJson()};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  ScanCommand ms2_ctx = makeMs2Ctx();

  Ms2Params p;    p.collision_energy = 29.0; p.activation_type = "HCD"; p.reaction_time = 0.0;
  Ms2Params ms3p; ms3p.collision_energy = 29.0; ms3p.activation_type = "HCD"; ms3p.reaction_time = 0.0;

  // MS2 winner: ambiguous mod [15,26] + non-bracketing b8 -> finalize leaves it [15,26].
  DeconvolvedSpectrum d101(101);
  d101.push_back(makeSyntheticPeakGroup(700.0 / 2.0 + 1.0, 700.0, 2));
  tracker.feedScan(precursor_id, 2, p, 101, d101, makeMatchMs2WithMod("b", 8, 700.0, 15, 26, 42.0), 1.0, ms2_ctx);
  tracker.finalizeMS2(precursor_id);

  int lo = -1, hi = -1, plo = -1, phi = -1;
  auto readRange = [&]() {
    const ProteoformModel* m = tracker.getModel(precursor_id);
    TEST_TRUE(m != nullptr && m->modifications.size() == 1)
    if (m != nullptr && m->modifications.size() == 1)
    { lo = m->modifications[0].candidate_start; hi = m->modifications[0].candidate_end; }
  };
  readRange();
  TEST_EQUAL(lo, 15)
  TEST_EQUAL(hi, 26)   // wide before any MS3 evidence

  // Fold one injected MS3 fragment, then assert the new range is NESTED within the prior (monotone).
  int scan = 200;
  auto fold = [&](const std::string& eq_t, int eq_i, bool inc, bool flip, const char* tag) {
    plo = lo; phi = hi;
    DeconvolvedSpectrum d(++scan);
    d.push_back(makeSyntheticPeakGroup(2000.0 / 2.0 + 1.0, 2000.0, 2));
    tracker.feedScan(precursor_id, 3, ms3p, scan, d, makeMatchMs3(eq_t, eq_i, 2000.0, inc, flip), -1.0, ms2_ctx);
    tracker.foldMs3(precursor_id, tag, tag);
    readRange();
    TEST_TRUE(lo >= plo && hi <= phi)   // NESTED within the prior range: never widens
    TEST_TRUE(lo <= hi)                 // still a valid range
  };

  fold("b", 22, true, false, "b22");    // prefix cover_end 22 brackets [15,26] -> [15,22]
  TEST_EQUAL(lo, 15)
  TEST_EQUAL(hi, 22)
  fold("b", 18, true, false, "b18");    // -> [15,18]
  TEST_EQUAL(lo, 15)
  TEST_EQUAL(hi, 18)
  fold("b", 30, true, false, "b30");    // fully-covering (non-bracketing) -> must NOT widen -> stays [15,18]
  TEST_EQUAL(lo, 15)
  TEST_EQUAL(hi, 18)
  fold("y", SEQ_LEN - 16 + 1, true, false, "y17");   // suffix cover_start 16 -> tightenLower -> [16,18]
  TEST_EQUAL(lo, 16)
  TEST_EQUAL(hi, 18)
  fold("b", 16, true, false, "b16");    // -> localized [16,16]
  TEST_EQUAL(lo, 16)
  TEST_EQUAL(hi, 16)
  fold("b", 30, true, false, "b30b");   // localized-stays-localized (non-bracketing) -> [16,16]
  TEST_EQUAL(lo, 16)
  TEST_EQUAL(hi, 16)
  // (The complement-flip verdict that pooled Pass B applies is pinned deterministically against the production
  //  predicate coveredAmbiguousInEquivFrame in FragmentAnalysis_test's Leg-2 parity check.)
}
END_SECTION

END_TEST
