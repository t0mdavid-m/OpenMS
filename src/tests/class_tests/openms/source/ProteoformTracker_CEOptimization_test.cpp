// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Issue-2 pin: per-fragment / per-charge CE optimization (ADR-0003).
//
// The ADR-0003 mechanism carries each fragment's BEST-MS2 acquisition params (collision energy,
// activation, reaction time) forward into that fragment's MS3 isolation stage (stage[0]). The
// engine implements it end to end:
//   feedScan -> mapScanOntoModel_ stores obs.params = scan params on every FragmentObservation,
//               and keeps the STRICTLY-greater best_ms2 per fragment (ProteoformTracker.cpp:654-655);
//   planNextScans copies the winning observation's params into Ms3Target.stage0_params
//               (ProteoformTracker.cpp:475);
//   ScanCommandQueue::buildMS3 first copies the ctx stage[0] (ScanCommandQueue.cpp:336) and then,
//               iff stage0_params != nullptr, OVERRIDES stage[0].collision_energy with the
//               per-fragment value (ScanCommandQueue.cpp:341-347).
//
// Before this test nothing EXERCISED that override: every capture replays a static spectrum whose
// fragment intensities do not change with the commanded CE, so the first variant (CE = ce_min)
// always won and stage0 CE was uniformly ce_min. A regression that deleted the stage0_params
// override (always use the ctx CE) would pass every existing test. This test forces two fragments
// to win their best-MS2 at DIFFERENT collision energies and pins the per-fragment propagation.
//
// CONSTRUCTION NOTE (why per-scan fragment PRESENCE, not literal inverted intensity): a fragment
// observation's intensity is recovered from PeakGroup::getChargeIntensity(charge) (feedScan ->
// PeakRecord.intensity), which is 0 for hand-built synthetic PeakGroups (no scoring pass). With every
// observation intensity 0 the strictly-greater update (obs.intensity > best->intensity) keeps the FIRST
// observation of each fragment. We make each fragment win deterministically by having it appear in
// EXACTLY ONE scan: fragment A (b6) in the CE=20 scan (the WINNER, whose fragments are pooled verbatim),
// fragment B (y6) in the CE=35 scan (a NON-winner, re-matched against the winner ladder -- Change 1:
// winner-anchored pooling). Because B is re-matched, its raw peak MASS must equal the real winner y6 ion
// (both fragments CONTAIN the ambiguous mod, so the mod is applied); we compute it in-test with the same
// computePTMAdjustedFragmentMasses the engine uses. Deterministic outcome: A -> 20, B -> 35.

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms2Params.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include <cmath>
#include <cstring>
#include <map>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  // Ambiguity objective: planNextScans targets the best-MS2 fragments that fully CONTAIN a still-
  // ambiguous modification (round-robin, strongest container first). We seed ONE ambiguous mod that
  // BOTH fragments contain (and neither brackets, so it survives narrowModifications_), so both are
  // emitted as per-fragment MS3 targets. (The Coverage objective cannot re-select here: a coverage
  // "gap" is by definition uncovered, so no fragment's coverage can CONTAIN it -> empty plan.)
  // No runtime block / no log files -> IdaLogger writes nothing. MS3 config drives buildMS3 stage[1].
  // max_targets >= 2 so planNextScans can emit BOTH per-fragment MS3 targets (budget = level(2).max_targets).
  const char* tracker_config = R"({
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
    "characterization": { "objective": "ambiguity", "protein_sequence": "PEPTIDEK" },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "intensity", "max_targets": 10 },
      "ms3": { "selection": "intensity", "max_targets": 10 }
    }
  })";

  // 8-residue winner proteoform (L = 8). Two fragments both CONTAIN one ambiguous modification at
  // residues [4,5] (so the Ambiguity objective keeps both as MS3 targets) yet neither BRACKETS it (so
  // narrowModifications_ leaves it ambiguous):
  //   fragment A = prefix b6  -> covers region residues [1, 6]; contains [4,5] (1<=4 && 6>=5); does not
  //                              bracket (prefix brackets iff 4 <= 6 < 5 -> false).
  //   fragment B = suffix y6  -> covers region residues [3, 8] (cover_start = L - 6 + 1 = 3); contains
  //                              [4,5] (3<=4 && 8>=5); does not bracket (suffix brackets iff 4 < 3 <= 5 -> false).
  // With the full-sequence match region (region_start/end = -1) both ion indices map straight to the
  // winner-region frame (ProteoformTracker.cpp mapScanOntoModel_ identity reduction).
  const char* WINNER_SEQ = "PEPTIDEK";
  const int   AMB_MOD_START = 4;   // ambiguous modification region (full-protein 1-based) both fragments contain
  const int   AMB_MOD_END   = 5;

  // Build a minimal synthetic PeakGroup at a given (mz, mono_mass, charge). getChargeIntensity is 0
  // for such hand-built groups (no scoring pass), which is fine: best-MS2 differentiation here is
  // driven by per-scan fragment PRESENCE, not intensity (see file header).
  PeakGroup makeSyntheticPeakGroup(double mz, double mono_mass, int charge)
  {
    PeakGroup pg(charge, charge, true);
    pg.setMonoisotopicMass(mono_mass);
    FLASHHelperClasses::LogMzPeak lp;
    lp.mz = mz;
    lp.abs_charge = charge;
    lp.intensity = 1000.0f;  // recorded on the peak; getChargeIntensity stays 0 without a scoring pass
    pg.push_back(lp);
    return pg;
  }

  // A single-fragment ProteoformMatch for the winner proteoform: identifies the full sequence (score
  // >= 0 so finalize accepts it as the winner), carries exactly one matched fragment ion, and seeds
  // the ambiguous modification [AMB_MOD_START, AMB_MOD_END] via a PTMSite (finalize copies the WINNER
  // match's ptm_sites into the model's modifications; both fed matches carry it so either could win).
  FragmentAnalysis::ProteoformMatch makeMatch(const std::string& ion_type, int ion_index, double observed_mass)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.score = 1.0;                 // >= 0 + non-empty sequence => eligible winner (finalize)
    m.region_start = -1;           // full-sequence region (identity frame mapping)
    m.region_end = -1;
    m.matched_protein = "synthetic";
    m.proteoform_sequence = WINNER_SEQ;
    FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
    fm.ion_type = ion_type;        // MS2 fragment: region-relative ion type/index
    fm.ion_index = ion_index;
    fm.observed_mass = observed_mass;
    m.fragments.push_back(fm);
    FragmentAnalysis::PTMSite site;
    site.start_position = AMB_MOD_START;   // ambiguous region (start < end => candidate_start < candidate_end)
    site.end_position = AMB_MOD_END;
    site.position = (AMB_MOD_START + AMB_MOD_END) / 2;
    site.mass_shift = 42.0106;             // arbitrary shift; no bracketing fragment => never localized here
    m.ptm_sites.push_back(site);
    return m;
  }

  // A single-fragment winner match with NO modification (unlike makeMatch, which seeds an ambiguous
  // PTM). The Coverage objective is mod-agnostic (it reasons over backbone-cleavage sites, not mods),
  // and omitting the mod keeps every fragment mass THEORETICAL, so non-winner scans re-match verbatim
  // against the winner ladder computed with computePTMAdjustedFragmentMasses(seq, {}, ...).
  FragmentAnalysis::ProteoformMatch makeCoverageFrag(const std::string& ion_type, int ion_index, double observed_mass)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.score = 1.0;
    m.region_start = -1;   // full-sequence region (identity frame mapping)
    m.region_end = -1;
    m.matched_protein = "synthetic";
    m.proteoform_sequence = WINNER_SEQ;
    FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
    fm.ion_type = ion_type;
    fm.ion_index = ion_index;
    fm.observed_mass = observed_mass;
    m.fragments.push_back(fm);
    return m;   // no ptm_sites -> no modifications
  }

  // A 2-stage MS2 context command whose stage[0] carries the FIRST scan's CE (20). buildMS3 copies
  // this stage[0] as the fallback; the per-fragment stage0_params override (if present) replaces the CE.
  ScanCommand makeMs2Ctx(double ms2_ctx_ce)
  {
    ScanCommand ctx{};
    ctx.scan_id = 1;
    ctx.msn_level = 2;
    ctx.num_stages = 1;
    ctx.mono_mass = 900.0;
    ctx.stages[0].precursor_mz = 451.0;
    ctx.stages[0].isolation_width = 2.0;
    ctx.stages[0].charge_state = 2;
    ctx.stages[0].collision_energy = ms2_ctx_ce;
    std::strncpy(ctx.stages[0].activation_type, "HCD", sizeof(ctx.stages[0].activation_type) - 1);
    return ctx;
  }

  // Find the planned Ms3Target for a given ion (type + index); returns nullptr if absent.
  const Ms3Target* findTarget(const std::vector<Ms3Target>& targets, const std::string& ion_type, int ion_index)
  {
    for (const Ms3Target& t : targets)
      if (t.ion_type == ion_type && t.ion_index == ion_index) return &t;
    return nullptr;
  }
}

START_TEST(ProteoformTracker_CEOptimization, "$Id$")

/////////////////////////////////////////////////////////////
// Per-fragment CE optimization: best-MS2 CE is carried per fragment into the MS3 stage[0] command.
/////////////////////////////////////////////////////////////
START_SECTION(per_fragment_best_ms2_ce_propagates_to_ms3_stage0)
{
  const double CE_A = 20.0;   // fragment A wins its best MS2 here
  const double CE_B = 35.0;   // fragment B wins its best MS2 here
  const int precursor_id = 7;

  Config cfg{std::string(tracker_config)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  ScanCommandQueue queue(cfg);

  // Fragment masses = the REAL winner (PEPTIDEK) b6 / y6 with the ambiguous mod applied (both ions
  // CONTAIN [4,5], so the mod is inside their coverage). The NON-winner scan (B) is re-matched against
  // the winner ladder (Change 1: winner-anchored pooling), so its raw peak MUST equal a real winner ion;
  // compute it with the same computePTMAdjustedFragmentMasses the engine uses -> exact (0 ppm) match.
  std::map<char, std::vector<double>> winner_ladder;
  {
    FragmentAnalysis::PTMSite amb;
    amb.start_position = AMB_MOD_START;
    amb.end_position = AMB_MOD_END;
    amb.position = (AMB_MOD_START + AMB_MOD_END) / 2;
    amb.mass_shift = 42.0106;
    FragmentAnalysis::computePTMAdjustedFragmentMasses(WINNER_SEQ, {amb}, {"b", "y"}, winner_ladder);
  }
  const double MASS_A = winner_ladder['b'][5];   // b6 (prefix; contains [4,5] -> mod applied)
  const double MASS_B = winner_ladder['y'][5];   // y6 (suffix; contains [4,5] -> mod applied)

  // The MS2 context captured at the FIRST feedScan: its stage[0] CE = CE_A (= 20). This is the value
  // buildMS3 would use for stage[0] IF the per-fragment override were absent (adversarial baseline).
  ScanCommand ms2_ctx = makeMs2Ctx(CE_A);

  // --- Scan 1: CE = 20, contains ONLY fragment A (prefix b6) -> A's best MS2 is at CE 20. ---
  {
    Ms2Params p1;
    p1.collision_energy = CE_A;
    p1.activation_type = "HCD";
    p1.reaction_time = 0.0;

    DeconvolvedSpectrum d1(101);
    d1.push_back(makeSyntheticPeakGroup(MASS_A / 2.0 + 1.0, MASS_A, 2));   // matches fragment A's mass

    auto match1 = makeMatch("b", 6, MASS_A);
    tracker.feedScan(precursor_id, 2, p1, 101, d1, match1, 1.0, ms2_ctx);
  }

  // --- Scan 2: CE = 35, contains ONLY fragment B (suffix y6) -> B's best MS2 is at CE 35. ---
  {
    Ms2Params p2;
    p2.collision_energy = CE_B;
    p2.activation_type = "HCD";
    p2.reaction_time = 0.0;

    DeconvolvedSpectrum d2(102);
    d2.push_back(makeSyntheticPeakGroup(MASS_B / 2.0 + 1.0, MASS_B, 2));   // matches fragment B's mass

    auto match2 = makeMatch("y", 6, MASS_B);
    tracker.feedScan(precursor_id, 2, p2, 102, d2, match2, 1.0, ms2_ctx);
  }

  tracker.finalizeMS2(precursor_id);

  // The winner model identified the proteoform and mapped both fragments.
  const ProteoformModel* mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL(mdl->proteoform_sequence, std::string(WINNER_SEQ))

  // --- planNextScans: per-fragment MS3 targets, each carrying its OWN best-MS2 CE in stage0_params. ---
  std::vector<Ms3Target> targets = tracker.planNextScans(precursor_id);
  TEST_TRUE(targets.size() >= 2)   // both fragments selected (Ambiguity objective: both contain the mod; budget >= 2)

  const Ms3Target* tA = findTarget(targets, "b", 6);
  const Ms3Target* tB = findTarget(targets, "y", 6);
  TEST_TRUE(tA != nullptr)
  TEST_TRUE(tB != nullptr)
  ABORT_IF(tA == nullptr || tB == nullptr)

  // CORE ASSERTION 1: the per-fragment best-MS2 CE DIFFERS (A -> 20, B -> 35). A regression that
  // collapsed per-fragment params to a single value (or always used the ctx CE) breaks this.
  TEST_TRUE(std::abs(tA->stage0_params.collision_energy - CE_A) < 1e-6)
  TEST_TRUE(std::abs(tB->stage0_params.collision_energy - CE_B) < 1e-6)
  TEST_TRUE(std::abs(tA->stage0_params.collision_energy - tB->stage0_params.collision_energy) > 1.0)

  // --- buildMS3: the per-fragment stage0 CE must reach cmd.stages[0].collision_energy. ---
  // Pass each target's stage0_params into buildMS3 (the Exploration executor does this at runtime).
  FragmentAnalysis::FragmentScores fs;  // default-empty scores are fine for the stage[0] CE check.
  ABORT_IF(cfg.level(3).scans.empty())
  const ScanConfig ms3_config = cfg.level(3).scans[0];  // analyzer CID, fragmentation CE 25 (stage[1])

  ScanCommand cmdA = queue.buildMS3(ms2_ctx, ms3_config,
                                    tA->frag_mz, tA->frag_charge, tA->iso_width, ms2_ctx.scan_id,
                                    'b', tA->ion_index, 1, fs, &tA->stage0_params);
  ScanCommand cmdB = queue.buildMS3(ms2_ctx, ms3_config,
                                    tB->frag_mz, tB->frag_charge, tB->iso_width, ms2_ctx.scan_id,
                                    'y', tB->ion_index, 1, fs, &tB->stage0_params);

  // CORE ASSERTION 2: buildMS3 propagated each fragment's own MS2 CE into MS3 stage[0].
  TEST_TRUE(std::abs(cmdA.stages[0].collision_energy - CE_A) < 1e-6)
  TEST_TRUE(std::abs(cmdB.stages[0].collision_energy - CE_B) < 1e-6)

  // ADVERSARIAL PROPERTY (non-vacuity): fragment B's stage[0] CE (35) DIFFERS from the ctx stage[0]
  // CE (20). buildMS3 copies ctx.stages[0] first (ScanCommandQueue.cpp:336) and only the
  // stage0_params override (ScanCommandQueue.cpp:341) raises it to 35. If that override line were
  // DELETED, cmdB.stages[0].collision_energy would equal the ctx CE (20) and this assertion FAILS.
  TEST_TRUE(std::abs(cmdB.stages[0].collision_energy - ms2_ctx.stages[0].collision_energy) > 1.0)
  TEST_TRUE(std::abs(cmdB.stages[0].collision_energy - CE_B) < 1e-6)
  // (cmdA's stage[0] CE happens to equal the ctx CE because the ctx was captured from scan 1 (CE 20);
  //  the discriminating case is fragment B, whose best-MS2 CE differs from the ctx CE.)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Coverage objective (backbone-cleavage-site model): planNextScans targets observed fragments whose
// span-INTERIOR carries a still-unwitnessed bond (so an MS3 re-feed can witness it), strongest-MS2
// first, marginal-skip (a fragment whose interior is already fully witnessed contributes nothing).
//
// WINNER_SEQ = PEPTIDEK, L=8 -> backbone bonds 1..7. Each fed fragment witnesses ONE bond in MS2
// (prefix b_k -> bond k; suffix y_m -> bond L-m). Its MS3-reachable interior is prefix [1..k] -> bonds
// 1..k-1; suffix [L-m+1..L] -> bonds (L-m+1)..(L-1). Fragments feed one-per-scan (winner pooled
// verbatim, non-winners re-matched by mass against the theoretical PEPTIDEK ladder). Synthetic peak
// groups have getChargeIntensity==0, so these fixtures are built to be ORDER-INDEPENDENT (marginal-skip
// via already-witnessed interiors, not pick order); intensity-driven ordering among overlapping useful
// parents is exercised end-to-end by the C# ms3_coverage_cytc golden.
/////////////////////////////////////////////////////////////
START_SECTION(planNextScans_coverage_targets_gaps_and_skips_redundant)
{
  const int precursor_id = 11;

  // Coverage config = the ambiguity config with the objective flipped (single occurrence).
  std::string cov_json = std::string(tracker_config);
  const std::string amb = "\"ambiguity\"";
  const std::string::size_type pos = cov_json.find(amb);
  TEST_TRUE(pos != std::string::npos)
  ABORT_IF(pos == std::string::npos)
  cov_json.replace(pos, amb.size(), "\"coverage\"");
  Config cfg{cov_json};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);

  // Theoretical PEPTIDEK b/y ladder (no mods) -> exact (0 ppm) re-match for non-winner scans.
  std::map<char, std::vector<double>> ladder;
  FragmentAnalysis::computePTMAdjustedFragmentMasses(WINNER_SEQ, {}, {"b", "y"}, ladder);

  ScanCommand ms2_ctx = makeMs2Ctx(25.0);
  int sid = 200;
  auto feedFrag = [&](const std::string& it, int idx) {
    const double mass = (it == "b") ? ladder['b'][idx - 1] : ladder['y'][idx - 1];
    Ms2Params p;
    p.collision_energy = 25.0;
    p.activation_type = "HCD";
    p.reaction_time = 0.0;
    DeconvolvedSpectrum d(sid);
    d.push_back(makeSyntheticPeakGroup(mass / 2.0 + 1.0, mass, 2));
    auto match = makeCoverageFrag(it, idx, mass);
    tracker.feedScan(precursor_id, 2, p, sid, d, match, 1.0, ms2_ctx);
    ++sid;
  };

  // Witnessed sites: b1->1, b2->2, b3->3, b6->6, y3->(L-3)=5  =>  {1,2,3,5,6}; uncovered bonds {4,7}.
  //   b6 interior {1..5} covers bond 4  -> SELECTED (prefix gap-filler)
  //   y3 interior {6,7}  covers bond 7  -> SELECTED (suffix gap-filler)
  //   b2 interior {1}, b3 interior {1,2}: already witnessed -> marginal-skip (order-independent)
  //   b1 has no interior bond.
  feedFrag("b", 1);
  feedFrag("b", 2);
  feedFrag("b", 3);
  feedFrag("b", 6);
  feedFrag("y", 3);
  tracker.finalizeMS2(precursor_id);

  const ProteoformModel* mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)

  std::vector<Ms3Target> targets = tracker.planNextScans(precursor_id);
  TEST_TRUE(findTarget(targets, "b", 6) != nullptr)   // gap -> spanning prefix parent selected
  TEST_TRUE(findTarget(targets, "y", 3) != nullptr)   // suffix parent fills the C-terminal gap
  TEST_TRUE(findTarget(targets, "b", 2) == nullptr)   // marginal-skip: interior {1} already witnessed
  TEST_TRUE(findTarget(targets, "b", 3) == nullptr)   // marginal-skip: interior {1,2} already witnessed
  TEST_TRUE(findTarget(targets, "b", 1) == nullptr)   // no interior bonds
  TEST_EQUAL(static_cast<int>(targets.size()), 2)
}
END_SECTION

START_SECTION(planNextScans_coverage_fully_witnessed_emits_nothing)
{
  const int precursor_id = 12;

  std::string cov_json = std::string(tracker_config);
  const std::string amb = "\"ambiguity\"";
  const std::string::size_type pos = cov_json.find(amb);
  ABORT_IF(pos == std::string::npos)
  cov_json.replace(pos, amb.size(), "\"coverage\"");
  Config cfg{cov_json};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);

  std::map<char, std::vector<double>> ladder;
  FragmentAnalysis::computePTMAdjustedFragmentMasses(WINNER_SEQ, {}, {"b", "y"}, ladder);

  ScanCommand ms2_ctx = makeMs2Ctx(25.0);
  int sid = 300;
  auto feedFrag = [&](const std::string& it, int idx) {
    const double mass = (it == "b") ? ladder['b'][idx - 1] : ladder['y'][idx - 1];
    Ms2Params p;
    p.collision_energy = 25.0;
    p.activation_type = "HCD";
    p.reaction_time = 0.0;
    DeconvolvedSpectrum d(sid);
    d.push_back(makeSyntheticPeakGroup(mass / 2.0 + 1.0, mass, 2));
    auto match = makeCoverageFrag(it, idx, mass);
    tracker.feedScan(precursor_id, 2, p, sid, d, match, 1.0, ms2_ctx);
    ++sid;
  };

  // b1..b7 witness EVERY backbone bond 1..7 -> uncovered is empty -> Coverage correctly plans nothing
  // (nothing left to drill). This is the case the CI run hit for cytC under the OLD span model too, but
  // now it is reached by genuine full coverage, not by the span-union self-contradiction.
  for (int k = 1; k <= 7; ++k) feedFrag("b", k);
  tracker.finalizeMS2(precursor_id);

  const ProteoformModel* mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)

  std::vector<Ms3Target> targets = tracker.planNextScans(precursor_id);
  TEST_EQUAL(static_cast<int>(targets.size()), 0)   // fully witnessed -> no coverage target
}
END_SECTION

// The SOURCE contract behind the exploration MS3 render-seed fix (Exploration.cpp: the render seed
// proto_ctx / group.proteoform_ctx is now taken from tracker->buildWinnerProteoformContext instead of
// the exploration-metric winner's frag_result). Under a CE sweep the metric winner (e.g. mass_count)
// can carry a FUSED single blind-mod decomposition while the identification-best (highest match.score /
// flash_extender_score) winner is SPLIT into its true mods. finalizeMS2 picks by match.score
// (ProteoformTracker.cpp:285) and buildWinnerProteoformContext exposes exactly the finalized model's
// mods (ProteoformTracker.cpp:568-576), so neither the pooled model NOR the render seed may inherit a
// fused metric-winner. Pre-regression (winner-by-count, or a fused buildWinnerProteoformContext) would
// report ONE +526 site. The END-TO-END wiring -- that Exploration seeds the MS3 leaf/command from this
// context so identification.tsv == pooled_identification -- is pinned by the C# exploration_ms3
// leaf==pooled invariant in FLASHIdaLogGolden_test.cs.
START_SECTION(winner_context_reflects_id_best_split_not_metric_fused)
{
  const int precursor_id = 9;

  Config cfg{std::string(tracker_config)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);

  ScanCommand ms2_ctx = makeMs2Ctx(30.0);

  // A ProteoformMatch with a given match.score + explicit PTM sites (one fragment so it identifies).
  auto makeScored = [](double score, const std::vector<FragmentAnalysis::PTMSite>& sites)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.score = score;                 // finalizeMS2 winner = highest match.score
    m.region_start = -1;             // full-sequence region (identity frame mapping)
    m.region_end = -1;
    m.matched_protein = "synthetic";
    m.proteoform_sequence = WINNER_SEQ;   // PEPTIDEK (L = 8)
    FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
    fm.ion_type = "b"; fm.ion_index = 6; fm.observed_mass = 700.0;
    m.fragments.push_back(fm);
    m.ptm_sites = sites;
    return m;
  };

  // id-BEST SPLIT (higher score): -89.0302 localized to residue 1 (N-term loss) + 615.2512 ambiguous [4,6].
  FragmentAnalysis::PTMSite s_nterm; s_nterm.start_position = 1; s_nterm.end_position = 1; s_nterm.position = 1; s_nterm.mass_shift = -89.0302;
  FragmentAnalysis::PTMSite s_heme;  s_heme.start_position = 4;  s_heme.end_position = 6;  s_heme.position = 5;  s_heme.mass_shift = 615.2512;
  FragmentAnalysis::ProteoformMatch split = makeScored(3.0, {s_nterm, s_heme});

  // metric-winner FUSED (lower score): ONE blind mod = summed mass over the union range [1,6].
  FragmentAnalysis::PTMSite s_fused; s_fused.start_position = 1; s_fused.end_position = 6; s_fused.position = 3; s_fused.mass_shift = 526.2210;
  FragmentAnalysis::ProteoformMatch fused = makeScored(1.0, {s_fused});

  Ms2Params p; p.collision_energy = 30.0; p.activation_type = "HCD"; p.reaction_time = 0.0;
  DeconvolvedSpectrum d1(201); d1.push_back(makeSyntheticPeakGroup(351.0, 700.0, 2));
  DeconvolvedSpectrum d2(202); d2.push_back(makeSyntheticPeakGroup(351.0, 700.0, 2));
  tracker.feedScan(precursor_id, 2, p, 201, d1, split, 3.0, ms2_ctx);
  tracker.feedScan(precursor_id, 2, p, 202, d2, fused, 1.0, ms2_ctx);
  tracker.finalizeMS2(precursor_id);

  // The finalized model must be the SPLIT winner: two mods, NOT the one fused blind mod.
  const ProteoformModel* mdl = tracker.getModel(precursor_id);
  TEST_TRUE(mdl != nullptr)
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL(static_cast<int>(mdl->modifications.size()), 2)

  // The render seed the fix uses (buildWinnerProteoformContext) must expose that SPLIT decomposition:
  // one ~ -89.03 site + one ~ +615.25 site, and NO fused ~ +526.22 site.
  MS3FragmentMatcher::ProteoformContext ctx = tracker.buildWinnerProteoformContext(precursor_id);
  TEST_EQUAL(static_cast<int>(ctx.ptm_sites.size()), 2)
  bool has_nterm = false, has_heme = false, has_fused = false;
  for (const FragmentAnalysis::PTMSite& s : ctx.ptm_sites)
  {
    if (std::abs(s.mass_shift - (-89.0302)) < 0.01) has_nterm = true;
    if (std::abs(s.mass_shift - 615.2512) < 0.01) has_heme = true;
    if (std::abs(s.mass_shift - 526.2210) < 0.01) has_fused = true;
  }
  TEST_EQUAL(has_nterm, true)
  TEST_EQUAL(has_heme, true)
  TEST_EQUAL(has_fused, false)
}
END_SECTION

END_TEST
