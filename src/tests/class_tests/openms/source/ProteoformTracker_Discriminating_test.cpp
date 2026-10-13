// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// characterization.mode: discriminating (ADR-0049) -- the tracker side of Refine: the candidate record kept on the
// model, the discriminating planner, the pooled MS3 piece verdicts, the round decision, and narrowing that reads
// shift-carrying fragments only under characterization.method refine.
//
// Driven through the public sequence feedScan -> finalizeMS2 -> planNextScans -> recordDiscriminatingScan ->
// decideRefine. The MS2 identification is the REAL one: FragmentAnalysis::getTopFragmentMatches under method refine runs
// Detect + Assemble on the synthetic spectrum and fills ProteoformMatch::refine, exactly as the engine does.
//
// THE SYNTHETIC CASE (PEPTIDEKAGRK, b/y, every Refine setting at its default: fdr 0.2, 300 shuffles, both exclusions,
// cap 20, jump penalty 2). Two shift-carrying series that no single path explains:
//   shared   b1 b2 b3, y1 y2 y3 y4 (unshifted)                       7 masses
//   a-only   b4..b11 + X   (X = +79.97 on residue 4)                  8 masses
//   b-only   b6..b11 + Y   (Y = +42.01, anywhere on residues 4..6)   6 masses
// Detect finds two levels, +79.97 (8 prefix ions) and +42.01 (6). X and Y are 37.96 Da apart -- not a water/ammonia
// loss, not 1-3 neutrons -- so both survive the exclusions; both sit on the 0.01 Da grid so the level IS the shift.
// No level has a suffix ion, so the totals are 0, +X, +Y and Assemble ranks them
//   +X: 15 ions, 1 jump -> S 13      +Y: 13 ions, 1 jump -> S 11      0: 7 ions -> S 7
// The two best paths share the 7 unshifted masses and explain the other 14 one-sidedly. Intensities interleave the two
// sides (a 1000, b 950, a 900, b 850, ...), so the records alternate a-only / b-only and a budget of 3 dispatches
// a, b, a. After MS3 the recount is 7 + confirmed - 2 for each side: a confirmed piece of b with none of a flips.

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms2Params.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Refine.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  const char* SEQ = "PEPTIDEKAGRK";
  const double SHIFT_X = 79.97; // on residue 4
  const double SHIFT_Y = 42.01; // on residues 4..6
  const double PROTON = 1.007276;

  // The neighbouring tracker fixtures' config, with characterization's mode / method / budget as parameters. The raw
  // strings carry an explicit delimiter (an undelimited R"(...)" ends at the first `)"`).
  std::string trackerConfig(const std::string& mode, const std::string& method, int max_targets, const std::string& log_dir = "")
  {
    std::string s = R"JSON({
  "deconvolution": {
    "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 1, "max_charge": 50,
    "min_mass": 50, "max_mass": 50000, "tol": [10, 10, 10]
  },
  "flashtnt": {"min_length": 3, "max_length": 8},
  "faims": {"cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15},
  "scheduling": {
    "cycle_time": {"enabled": false, "value_ms": 60000},
    "scan_timeout": {"enabled": false, "value_ms": 30000}
  },
  "files": {"target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": ""},
  "conditional_ms2": false,
  "precursor_selection": {
    "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false,
    "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 3
  },
  "characterization": {
    "mode": ")JSON";
    s += mode + "\",\n    \"method\": \"" + method + "\",\n    \"protein_sequence\": \"" + SEQ + "\",\n    \"max_targets\": " +
         std::to_string(max_targets);
    s += R"JSON(
  },
  "ms_settings": {
    "ms1": {"analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246},
    "ms2": {"analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000},
    "ms3": {"analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25, "resolution": 120000}
  },
  "tagging": {},
  "quantification": {"enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4})JSON";
    if (!log_dir.empty()) s += ",\n  \"runtime\": {\"log_dir\": \"" + log_dir + "\"}";
    s += "\n}\n";
    return s;
  }

  // One deconvolved mass: a single charge-1 peak carrying the intensity, so PeakGroup::getIntensity() -- what
  // runRefineMatching_ hands Assemble and the tracker ranks the records by -- is that intensity.
  PeakGroup makePeak(double mass, double intensity)
  {
    PeakGroup pg(1, 1, true);
    FLASHHelperClasses::LogMzPeak lp;
    lp.mz = mass + PROTON;
    lp.intensity = static_cast<float>(intensity);
    lp.abs_charge = 1;
    lp.isotopeIndex = 0;
    pg.push_back(lp);
    pg.updateMonoMassAndIsotopeIntensities(0.0);
    pg.setMonoisotopicMass(mass);
    pg.setRepAbsCharge(1);
    pg.setQscore(0.9);
    return pg;
  }

  struct MassList
  {
    std::vector<double> masses, intensities;
    void add(double m, double i)
    {
      masses.push_back(m);
      intensities.push_back(i);
    }
  };

  // The two-totals case of the header. y_m = suffix[L - 1 - m] (Refine::Ladder: suffix[n-1] = y_{L-n}).
  MassList twoTotals()
  {
    const Refine::Ladder lad = Refine::ladder(SEQ, {"b", "y"});
    const int L = static_cast<int>(std::strlen(SEQ));
    MassList ml;
    for (int n = 1; n <= 3; ++n) ml.add(lad.prefix[n - 1], 50.0);
    for (int m = 1; m <= 4; ++m) ml.add(lad.suffix[L - 1 - m], 50.0);
    for (int k = 0; k < 8; ++k) ml.add(lad.prefix[3 + k] + SHIFT_X, 1000.0 - 100.0 * k); // b4..b11 + X
    for (int k = 0; k < 6; ++k) ml.add(lad.prefix[5 + k] + SHIFT_Y, 950.0 - 100.0 * k);  // b6..b11 + Y
    return ml;
  }

  // The best path's shift is NOT localized: b3 unshifted and b6..b11 + X (no b4/b5) put X anywhere on residues 4..6, so
  // the ambiguity branch has a window to plan containing fragments for. b7..b11 + Y keep the second path one ion behind:
  //   +X: 7 + 6 = 13 ions, 1 jump -> S 11      +Y: 7 + 5 = 12 ions, 1 jump -> S 10
  MassList ambiguousBest()
  {
    const Refine::Ladder lad = Refine::ladder(SEQ, {"b", "y"});
    const int L = static_cast<int>(std::strlen(SEQ));
    MassList ml;
    for (int n = 1; n <= 3; ++n) ml.add(lad.prefix[n - 1], 50.0);
    for (int m = 1; m <= 4; ++m) ml.add(lad.suffix[L - 1 - m], 50.0);
    for (int k = 0; k < 6; ++k) ml.add(lad.prefix[5 + k] + SHIFT_X, 1000.0 - 100.0 * k); // b6..b11 + X
    for (int k = 0; k < 5; ++k) ml.add(lad.prefix[6 + k] + SHIFT_Y, 950.0 - 100.0 * k);  // b7..b11 + Y
    return ml;
  }

  // trackerConfig with characterization.min_fragment_charge set: every synthetic fragment is charge 1, so a floor of 2
  // leaves no record isolatable (exhaustivePoolAdmits) while identification and the ambiguity branch are untouched.
  std::string withFragmentChargeFloor(std::string cfg, int floor)
  {
    const std::string key = "\"max_targets\": ";
    const size_t at = cfg.find(key);
    if (at == std::string::npos) return cfg;
    const size_t end = cfg.find_first_not_of("0123456789", at + key.size());
    cfg.insert(end, ",\n    \"min_fragment_charge\": " + std::to_string(floor));
    return cfg;
  }

  // Unshifted ions only: every mass is explained at dm = 0, the Signal is 0 everywhere, no level -> one total.
  MassList oneTotal()
  {
    const Refine::Ladder lad = Refine::ladder(SEQ, {"b", "y"});
    const int L = static_cast<int>(std::strlen(SEQ));
    MassList ml;
    for (int n = 1; n <= 9; ++n) ml.add(lad.prefix[n - 1], 100.0 + n);
    for (int m = 1; m <= 6; ++m) ml.add(lad.suffix[L - 1 - m], 200.0 + m);
    return ml;
  }

  // The order runRefineMatching_ and the tracker both use: ascending mass.
  void sortedCopy(const MassList& ml, std::vector<double>& masses, std::vector<double>& intensities)
  {
    std::vector<size_t> order(ml.masses.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return ml.masses[a] < ml.masses[b]; });
    masses.clear();
    intensities.clear();
    for (size_t i : order)
    {
      masses.push_back(ml.masses[i]);
      intensities.push_back(ml.intensities[i]);
    }
  }

  ScanCommand makeMs2Ctx()
  {
    ScanCommand ctx{};
    ctx.scan_id = 1;
    ctx.msn_level = 2;
    ctx.num_stages = 1;
    ctx.mono_mass = 1500.0;
    ctx.stages[0].precursor_mz = 751.0;
    ctx.stages[0].isolation_width = 2.0;
    ctx.stages[0].charge_state = 2;
    ctx.stages[0].collision_energy = 29.0;
    std::strncpy(ctx.stages[0].activation_type, "HCD", sizeof(ctx.stages[0].activation_type) - 1);
    return ctx;
  }

  // Identify @p ml with the engine's own MS2 path (method refine) and feed the result as the one MS2 scan of the
  // precursor, then finalize. One scan: it is the winner by construction.
  void feedAndFinalize(const Config& cfg, ProteoformTracker& tracker, int precursor_id, int scan_id, const MassList& ml,
                       const std::string& activation = "HCD")
  {
    DeconvolvedSpectrum ms2(scan_id);
    for (size_t i = 0; i < ml.masses.size(); ++i) ms2.push_back(makePeak(ml.masses[i], ml.intensities[i]));
    PeakGroup precursor;
    precursor.setMonoisotopicMass(Refine::ladder(SEQ, {"b", "y"}).sequence_mass + SHIFT_X);
    ms2.setPrecursorPeakGroup(precursor);

    FragmentAnalysis fa(cfg);
    double masses[64], qscores[64], wstarts[64], wends[64];
    int charges[64], frag_indices[64];
    char ion_types[64];
    FragmentAnalysis::ProteoformMatch match;
    fa.getTopFragmentMatches(SEQ, 64, masses, qscores, charges, wstarts, wends, ion_types, frag_indices, ms2, match, "HCD", 0.0, nullptr);

    Ms2Params p;
    p.collision_energy = 29.0;
    p.activation_type = activation; // the winner's stage-0 activation (winner_params); identification above stays b/y
    p.reaction_time = 0.0;
    tracker.feedScan(precursor_id, 2, p, scan_id, ms2, match, match.score, makeMs2Ctx());
    tracker.finalizeMS2(precursor_id);
  }

  // The shift-carrying piece pieceTest would predict for @p ion under @p path: the same subsequence, rebased sites and
  // MS3 ion types (Refine::pieceTest's predict), the heaviest piece whose mass differs from its bare backbone.
  double shiftCarryingPiece(const Refine::Path& path, const Refine::Path::Ion& ion)
  {
    const std::string seq(SEQ);
    const int L = static_cast<int>(seq.size());
    MS3FragmentMatcher::ProteoformContext ctx;
    ctx.region_start = 0;
    ctx.region_end = L;
    ctx.ptm_sites = Refine::ptmSites(path);
    const std::string sub = MS3FragmentMatcher::extractSubsequence(seq, ctx, ion.type, ion.index);
    const bool prefix = ion.type == 'a' || ion.type == 'b' || ion.type == 'c';
    const int start = prefix ? 1 : L - ion.index + 1;
    const auto rebased = MS3FragmentMatcher::rebasePTMSites(ctx.ptm_sites, start, ion.index);
    const auto types = MS3FragmentMatcher::getMS3IonTypes(ion.type);
    const auto theo = MS3FragmentMatcher::computeTheoreticalMasses(sub, types, rebased);
    const auto bare = MS3FragmentMatcher::computeTheoreticalMasses(sub, types, {});
    double best = 0.0;
    for (const auto& t : theo)
    {
      for (const auto& b : bare)
      {
        if (b.ion_type == t.ion_type && b.position == t.position && std::fabs(t.mass - b.mass) > t.mass * 10e-6 && t.mass > best)
          best = t.mass;
      }
    }
    return best;
  }

  DeconvolvedSpectrum ms3Of(const std::vector<double>& masses)
  {
    DeconvolvedSpectrum s(0);
    for (double m : masses)
    {
      PeakGroup pg;
      pg.setMonoisotopicMass(m);
      s.push_back(pg);
    }
    return s;
  }

  // An MS3 spectrum whose one piece fits the path that explains @p rec (a for an a-only record, b for a b-only one).
  DeconvolvedSpectrum fittingScan(const ProteoformModel& m, const DiscriminatingRecord& rec)
  {
    const auto& c = m.refine->candidates->assemble.candidates;
    if (rec.target.ion_a) return ms3Of({shiftCarryingPiece(c[0], *rec.target.ion_a)});
    return ms3Of({shiftCarryingPiece(c[1], *rec.target.ion_b)});
  }

  std::string label(const std::optional<Refine::Path::Ion>& ion)
  {
    return ion ? std::string(1, ion->type) + std::to_string(ion->index) : std::string("-");
  }

  int countDispatched(const ProteoformModel& m)
  {
    int n = 0;
    for (const auto& r : m.refine->records) n += r.dispatched;
    return n;
  }

  std::vector<std::string> splitTab(const std::string& line)
  {
    std::vector<std::string> out;
    std::string cell;
    std::stringstream ss(line);
    while (std::getline(ss, cell, '\t')) out.push_back(cell);
    if (!line.empty() && line.back() == '\t') out.push_back("");
    return out;
  }

  // planNextScans with its stdout captured: the [MS3-PLAN] line names the reason for a zero-target plan.
  std::vector<Ms3Target> planCapturing(ProteoformTracker& tracker, int precursor_id, std::string& out)
  {
    std::ostringstream captured;
    std::streambuf* const real_cout = std::cout.rdbuf(captured.rdbuf());
    std::vector<Ms3Target> plan;
    try
    {
      plan = tracker.planNextScans(precursor_id);
    }
    catch (...)
    {
      std::cout.rdbuf(real_cout);
      throw;
    }
    std::cout.rdbuf(real_cout);
    out = captured.str();
    std::cout << out; // keep the marker visible in the test log
    return plan;
  }

  const double LOOSE = MS3FragmentMatcher::LOOSE_TOLERANCE_PPM;
  const double TIGHT = 10.0;
} // namespace

START_TEST(ProteoformTracker_Discriminating, "$Id$")

/////////////////////////////////////////////////////////////
// 1. The records are Refine::discriminating(c0, c1) of the winner scan, intensity-desc, and the plan is their top
//    `budget`, labelled with the best path's ion; a mass only the second path explains goes out unassigned ('u'/0).
//
// FAILS UNDER: a tracker that does not keep ProteoformMatch::refine at finalize (no records); a mass list rebuilt in an
// order other than Refine's (the paths' mass_idx then point at the wrong masses); ranking by the representative-charge
// intensity (0 for these groups) instead of the group intensity (the order would become mass order); a planner that
// walks m.fragments (the ambiguity pool) instead of the records.
/////////////////////////////////////////////////////////////
START_SECTION(plan_targets_are_the_masses_the_two_best_paths_explain_differently)
{
  const int pid = 31;
  Config cfg{trackerConfig("discriminating", "refine", 3)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  const MassList ml = twoTotals();
  feedAndFinalize(cfg, tracker, pid, 601, ml);

  const ProteoformModel* m = tracker.getModel(pid);
  TEST_TRUE(m != nullptr)
  ABORT_IF(m == nullptr)
  TEST_TRUE(m->refine.has_value())
  ABORT_IF(!m->refine.has_value())
  const auto& c = m->refine->candidates->assemble.candidates;
  TEST_EQUAL(c.size(), 3) // totals 0, +X, +Y
  ABORT_IF(c.size() < 2)
  TEST_REAL_SIMILAR(c[0].total, SHIFT_X)
  TEST_REAL_SIMILAR(c[1].total, SHIFT_Y)
  TEST_EQUAL(c[0].score, 13)
  TEST_EQUAL(c[1].score, 11)

  std::vector<double> masses, intensities;
  sortedCopy(ml, masses, intensities);
  const std::vector<Refine::Discriminating> expected = Refine::discriminating(c[0], c[1], masses, intensities);
  const auto& recs = m->refine->records;
  TEST_EQUAL(recs.size(), 14)
  TEST_EQUAL(recs.size(), expected.size())
  ABORT_IF(recs.size() != expected.size())
  for (size_t i = 0; i < recs.size(); ++i)
  {
    TEST_REAL_SIMILAR(recs[i].target.mass, expected[i].mass)
    TEST_EQUAL(label(recs[i].target.ion_a), label(expected[i].ion_a))
    TEST_EQUAL(label(recs[i].target.ion_b), label(expected[i].ion_b))
    TEST_TRUE(recs[i].peak_index >= 0)
    if (recs[i].peak_index >= 0) TEST_REAL_SIMILAR(m->winner_peaks[recs[i].peak_index].mono_mass, recs[i].target.mass)
    if (i > 0) TEST_TRUE(recs[i].target.intensity <= recs[i - 1].target.intensity)
    TEST_EQUAL(recs[i].dispatched, false)
  }
  // interleaved: a 1000 (b4+X), b 950 (b6+Y), a 900 (b5+X)
  TEST_EQUAL(label(recs[0].target.ion_a), "b4")
  TEST_EQUAL(label(recs[0].target.ion_b), "-")
  TEST_EQUAL(label(recs[1].target.ion_a), "-")
  TEST_EQUAL(label(recs[1].target.ion_b), "b6")
  TEST_EQUAL(label(recs[2].target.ion_a), "b5")

  const std::vector<Ms3Target> plan = tracker.planNextScans(pid);
  TEST_EQUAL(plan.size(), 3)
  ABORT_IF(plan.size() != 3)
  m = tracker.getModel(pid);
  for (size_t i = 0; i < 3; ++i)
  {
    const auto& r = m->refine->records[i];
    TEST_REAL_SIMILAR(plan[i].frag_mass, r.target.mass)
    TEST_EQUAL(plan[i].ion_type, r.target.ion_a ? std::string(1, r.target.ion_a->type) : std::string("u"))
    TEST_EQUAL(plan[i].ion_index, r.target.ion_a ? r.target.ion_a->index : 0)
    TEST_TRUE(plan[i].frag_mz > 0.0)
    TEST_EQUAL(plan[i].frag_charge, 1)
    TEST_EQUAL(plan[i].stage0_params.activation_type, std::string("HCD"))
    TEST_EQUAL(r.dispatched, true)
    TEST_EQUAL(m->dispatched_nominal_masses.count(SpectralDeconvolution::getNominalMass(r.target.mass)), 1)
  }
  TEST_EQUAL(countDispatched(*m), 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 2. One total (Detect found no level): nothing tells two paths apart -> no records, no plan (reason single_candidate),
//    and the model is still identified.
//
// FAILS UNDER: building records from a single candidate (c[1] out of range); falling through to the ambiguity branch
// before a decision -- that branch also plans nothing here (no modification), so the REASON is what tells them apart:
// the discriminating planner's single_candidate, not ambiguity's all_mods_localized.
/////////////////////////////////////////////////////////////
START_SECTION(single_candidate_dispatches_nothing)
{
  const int pid = 32;
  Config cfg{trackerConfig("discriminating", "refine", 3)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedAndFinalize(cfg, tracker, pid, 602, oneTotal());

  const ProteoformModel* m = tracker.getModel(pid);
  TEST_TRUE(m != nullptr)
  ABORT_IF(m == nullptr)
  TEST_EQUAL(m->proteoform_sequence, std::string(SEQ))
  TEST_TRUE(m->refine.has_value())
  ABORT_IF(!m->refine.has_value())
  TEST_EQUAL(m->refine->candidates->assemble.candidates.size(), 1)
  TEST_EQUAL(m->refine->records.size(), 0)
  std::string log;
  TEST_EQUAL(planCapturing(tracker, pid, log).size(), 0)
  TEST_TRUE(log.find("objective=discriminating targets=0 reason=single_candidate") != std::string::npos)
  TEST_EQUAL(m->dispatched_nominal_masses.size(), 0)
  TEST_EQUAL(tracker.decideRefine(pid), false)
  TEST_EQUAL(m->refine->decided, false)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 3. A second plan before any MS3 returned takes the NEXT records, never one already dispatched.
//
// FAILS UNDER: a planner without the dispatched flag / nominal-mass memory (it would return records 0..2 again).
/////////////////////////////////////////////////////////////
START_SECTION(dispatched_masses_are_not_replanned)
{
  const int pid = 33;
  Config cfg{trackerConfig("discriminating", "refine", 3)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedAndFinalize(cfg, tracker, pid, 603, twoTotals());

  const std::vector<Ms3Target> first = tracker.planNextScans(pid);
  const std::vector<Ms3Target> second = tracker.planNextScans(pid);
  TEST_EQUAL(first.size(), 3)
  TEST_EQUAL(second.size(), 3)
  const ProteoformModel* m = tracker.getModel(pid);
  ABORT_IF(second.size() != 3 || m == nullptr)
  for (size_t i = 0; i < 3; ++i) TEST_REAL_SIMILAR(second[i].frag_mass, m->refine->records[3 + i].target.mass)
  std::set<int> nominal;
  for (const auto& t : first) nominal.insert(SpectralDeconvolution::getNominalMass(t.frag_mass));
  for (const auto& t : second) nominal.insert(SpectralDeconvolution::getNominalMass(t.frag_mass));
  TEST_EQUAL(nominal.size(), 6)
  TEST_EQUAL(countDispatched(*m), 6)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 3a. A winner acquired with an activation that cannot feed an MS3 (ETD) plans no discriminating MS3, though two
//     candidates and 14 records exist; the round can never complete, so it falls through (the reason names the
//     capability) and nothing is marked dispatched or stamped.
//
// FAILS UNDER: a discriminating planner without planExhaustive_'s capability guard (it would command 3 MS3 whose
// stage 0 replays the ETD winner).
/////////////////////////////////////////////////////////////
START_SECTION(etd_winner_plans_no_ms3)
{
  const int pid = 39;
  Config cfg{trackerConfig("discriminating", "refine", 3)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedAndFinalize(cfg, tracker, pid, 609, twoTotals(), "ETD");
  const ProteoformModel* m = tracker.getModel(pid);
  ABORT_IF(m == nullptr || !m->refine.has_value())
  TEST_EQUAL(m->winner_params.activation_type, std::string("ETD"))
  TEST_EQUAL(m->refine->candidates->assemble.candidates.size() >= 2, true)
  TEST_EQUAL(m->refine->records.size(), 14)

  std::string log;
  TEST_EQUAL(planCapturing(tracker, pid, log).size(), 0)
  TEST_TRUE(log.find("[DISCRIMINATING] precursor_id=39 decided=fallthrough reason=winner_scan_not_ms3_capable") != std::string::npos)
  TEST_EQUAL(m->refine->decided, true)
  TEST_EQUAL(countDispatched(*m), 0)
  TEST_EQUAL(m->dispatched_nominal_masses.size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 3b (Review Focus 3). The round completes only when EVERY dispatched mass has returned; a second scan of a mass already
//    returned does not stand in for a missing one.
//
// FAILS UNDER: completing on the first return, on a scan count, or on returned records only.
/////////////////////////////////////////////////////////////
START_SECTION(round_waits_for_every_dispatched_mass)
{
  const int pid = 34;
  Config cfg{trackerConfig("discriminating", "refine", 3)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedAndFinalize(cfg, tracker, pid, 604, twoTotals());
  TEST_EQUAL(tracker.planNextScans(pid).size(), 3)
  const ProteoformModel* m = tracker.getModel(pid);
  ABORT_IF(m == nullptr || !m->refine.has_value())
  const double m0 = m->refine->records[0].target.mass, m1 = m->refine->records[1].target.mass, m2 = m->refine->records[2].target.mass;
  const DeconvolvedSpectrum empty = ms3Of({});

  TEST_EQUAL(tracker.recordDiscriminatingScan(pid, m0, empty, LOOSE, TIGHT), false)
  TEST_EQUAL(tracker.recordDiscriminatingScan(pid, m1, empty, LOOSE, TIGHT), false)
  TEST_EQUAL(tracker.recordDiscriminatingScan(pid, m0, empty, LOOSE, TIGHT), false) // a repeat, not the missing third
  TEST_EQUAL(m->refine->records[0].scans, 2)
  TEST_EQUAL(tracker.recordDiscriminatingScan(pid, m2, empty, LOOSE, TIGHT), true)
  // a mass no record holds changes nothing
  TEST_EQUAL(tracker.recordDiscriminatingScan(pid, 5555.5, empty, LOOSE, TIGHT), false)
  // after the decision the round never reports complete again
  tracker.decideRefine(pid);
  TEST_EQUAL(tracker.recordDiscriminatingScan(pid, m2, empty, LOOSE, TIGHT), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 3c (decision 9). An MS3 that ARRIVED but deconvolved to nothing is received with no evidence: the engine's MS3 return
//    path hands recordDiscriminatingScan an empty DeconvolvedSpectrum(0), which marks the record returned with neither
//    path fitting. Two dispatched, both empty: the round completes on the second, and the recount over confirmed 0/0
//    keeps the best path.
//
// FAILS UNDER: an empty scan that does not mark its record returned (the round would never complete), one that throws,
// or one that counts as a fit (a false positive for either path).
/////////////////////////////////////////////////////////////
START_SECTION(empty_ms3_counts_as_returned_without_evidence)
{
  const int pid = 40;
  Config cfg{trackerConfig("discriminating", "refine", 2)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  // A first MS2 that identifies nothing, acquired at another energy: it fixes m.ms2_ctx (captured once, from the first
  // fed scan) but cannot win. The re-dispatch must replay the WINNER's command, so winner_cmd is the second scan's.
  {
    ScanCommand first_cmd = makeMs2Ctx();
    first_cmd.scan_id = 2;
    first_cmd.stages[0].collision_energy = 17.0;
    Ms2Params p17;
    p17.collision_energy = 17.0;
    p17.activation_type = "HCD";
    tracker.feedScan(pid, 2, p17, 609, ms3Of({}), FragmentAnalysis::ProteoformMatch{}, -1.0, first_cmd);
  }
  feedAndFinalize(cfg, tracker, pid, 610, twoTotals());
  TEST_EQUAL(tracker.planNextScans(pid).size(), 2)
  const ProteoformModel* m = tracker.getModel(pid);
  ABORT_IF(m == nullptr || !m->refine.has_value())
  ABORT_IF(m->winner_cmd == nullptr)
  TEST_EQUAL(m->winner_scan_id, 610)
  TEST_REAL_SIMILAR(m->winner_cmd->stages[0].collision_energy, 29.0) // the winner's (makeMs2Ctx), not the first scan's 17
  TEST_REAL_SIMILAR(m->ms2_ctx.stages[0].collision_energy, 17.0)     // ms2_ctx is the first fed scan's: why it cannot serve
  TEST_EQUAL(countDispatched(*m), 2)
  const auto& recs = m->refine->records;
  const DeconvolvedSpectrum no_evidence(0); // exactly what FLASHIda.cpp passes for a null/empty MS3

  TEST_EQUAL(tracker.recordDiscriminatingScan(pid, recs[0].target.mass, no_evidence, LOOSE, TIGHT), false)
  TEST_EQUAL(recs[0].returned, true)
  TEST_EQUAL(recs[0].fits_a, false)
  TEST_EQUAL(recs[0].fits_b, false)
  TEST_EQUAL(tracker.recordDiscriminatingScan(pid, recs[1].target.mass, no_evidence, LOOSE, TIGHT), true)
  TEST_EQUAL(recs[1].returned, true)
  TEST_EQUAL(recs[1].fits_a, false)
  TEST_EQUAL(recs[1].fits_b, false)

  TEST_EQUAL(tracker.decideRefine(pid), false) // confirmed 0/0: no flip
  TEST_EQUAL(m->refine->decided, true)
  TEST_EQUAL(m->refine->winner_candidate, 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 4. Pieces pool across the scans of one mass: a later scan that fits sets the verdict, a still later one that fits
//    nothing does not clear it.
//
// FAILS UNDER: assigning the latest verdict instead of OR-ing it in.
/////////////////////////////////////////////////////////////
START_SECTION(pieces_pool_across_scans_of_one_mass)
{
  const int pid = 35;
  Config cfg{trackerConfig("discriminating", "refine", 3)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedAndFinalize(cfg, tracker, pid, 605, twoTotals());
  TEST_EQUAL(tracker.planNextScans(pid).size(), 3)
  const ProteoformModel* m = tracker.getModel(pid);
  ABORT_IF(m == nullptr || !m->refine.has_value())
  const DiscriminatingRecord& r0 = m->refine->records[0]; // a-only (b4 + X)
  TEST_TRUE(r0.target.ion_a.has_value())

  tracker.recordDiscriminatingScan(pid, r0.target.mass, ms3Of({}), LOOSE, TIGHT);
  TEST_EQUAL(r0.fits_a, false)
  TEST_EQUAL(r0.returned, true)
  tracker.recordDiscriminatingScan(pid, r0.target.mass, fittingScan(*m, r0), LOOSE, TIGHT);
  TEST_EQUAL(r0.fits_a, true)
  TEST_EQUAL(r0.fits_b, false)
  tracker.recordDiscriminatingScan(pid, r0.target.mass, ms3Of({}), LOOSE, TIGHT);
  TEST_EQUAL(r0.fits_a, true)
  TEST_EQUAL(r0.scans, 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 5. The recount keeps the best path when its piece is confirmed: decided, no flip, the model untouched, and one pooled
//    row with trigger MS3-DECIDE.
//
// FAILS UNDER: a decision that does not mark the round decided or emits no row; a recount fed only the dispatched
// records (Refine::recount throws: the two shared counts disagree).
/////////////////////////////////////////////////////////////
START_SECTION(decide_keeps_the_best_when_confirmed)
{
  const int pid = 36;
  const std::string log_dir = "testlogs/pt_discriminating_keep";
  File::removeDirRecursively(log_dir);
  File::makeDir(log_dir);
  {
    Config cfg{trackerConfig("discriminating", "refine", 3, log_dir)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    feedAndFinalize(cfg, tracker, pid, 606, twoTotals());
    TEST_EQUAL(tracker.planNextScans(pid).size(), 3)
    const ProteoformModel* m = tracker.getModel(pid);
    ABORT_IF(m == nullptr || !m->refine.has_value())
    const auto mods_before = m->modifications;
    const size_t frags_before = m->fragments.size();
    const auto& recs = m->refine->records;
    tracker.recordDiscriminatingScan(pid, recs[0].target.mass, fittingScan(*m, recs[0]), LOOSE, TIGHT); // a confirmed
    tracker.recordDiscriminatingScan(pid, recs[1].target.mass, ms3Of({}), LOOSE, TIGHT);
    TEST_EQUAL(tracker.recordDiscriminatingScan(pid, recs[2].target.mass, ms3Of({}), LOOSE, TIGHT), true)

    TEST_EQUAL(tracker.decideRefine(pid), false)
    TEST_EQUAL(m->refine->decided, true)
    TEST_EQUAL(m->refine->winner_candidate, 0)
    TEST_EQUAL(m->update_index, 2)
    TEST_EQUAL(m->modifications.size(), mods_before.size())
    for (size_t i = 0; i < mods_before.size() && i < m->modifications.size(); ++i)
    {
      TEST_EQUAL(m->modifications[i].candidate_start, mods_before[i].candidate_start)
      TEST_EQUAL(m->modifications[i].candidate_end, mods_before[i].candidate_end)
    }
    TEST_EQUAL(m->fragments.size(), frags_before)
    TEST_EQUAL(tracker.decideRefine(pid), false) // fires once
    TEST_EQUAL(m->update_index, 2)
  }
  std::ifstream f((log_dir + "/pooled_identification.tsv").c_str());
  std::string header, line;
  TEST_TRUE(static_cast<bool>(std::getline(f, header)))
  const std::vector<std::string> cols = splitTab(header);
  int ci_trigger = -1, ci_update = -1;
  for (int j = 0; j < static_cast<int>(cols.size()); ++j)
  {
    if (cols[j] == "trigger") ci_trigger = j;
    if (cols[j] == "update_index") ci_update = j;
  }
  TEST_TRUE(ci_trigger >= 0 && ci_update >= 0)
  std::vector<std::vector<std::string>> rows;
  while (std::getline(f, line))
    if (!line.empty()) rows.push_back(splitTab(line));
  TEST_EQUAL(rows.size(), 2) // MS2 baseline + the decision
  ABORT_IF(rows.size() != 2 || ci_trigger < 0 || ci_update < 0)
  TEST_EQUAL(rows[0][ci_trigger], "MS2")
  TEST_EQUAL(rows[1][ci_trigger], "MS3-DECIDE")
  TEST_EQUAL(rows[1][ci_update], "2")
}
END_SECTION

/////////////////////////////////////////////////////////////
// 6. A confirmed piece of the second path, none of the first: the model flips -- modifications re-seeded from the second
//    path's windows, fragments rebuilt from the winner scan as that path explains it, update_index + 1.
//
// FAILS UNDER: a decision that never flips; a flip that keeps the first path's modifications or fragments.
/////////////////////////////////////////////////////////////
START_SECTION(decide_flips_and_reseeds_modifications)
{
  const int pid = 37;
  Config cfg{trackerConfig("discriminating", "refine", 3)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedAndFinalize(cfg, tracker, pid, 607, twoTotals());
  TEST_EQUAL(tracker.planNextScans(pid).size(), 3)
  const ProteoformModel* m = tracker.getModel(pid);
  ABORT_IF(m == nullptr || !m->refine.has_value())
  TEST_TRUE(m->fragments.count(FragmentKey{"b", 4}) == 1) // the first path's b4 + X
  const auto& recs = m->refine->records;
  TEST_TRUE(recs[1].target.ion_b.has_value())
  tracker.recordDiscriminatingScan(pid, recs[0].target.mass, ms3Of({}), LOOSE, TIGHT);
  tracker.recordDiscriminatingScan(pid, recs[1].target.mass, fittingScan(*m, recs[1]), LOOSE, TIGHT); // b confirmed
  TEST_EQUAL(recs[1].fits_b, true)
  TEST_EQUAL(tracker.recordDiscriminatingScan(pid, recs[2].target.mass, ms3Of({}), LOOSE, TIGHT), true)

  TEST_EQUAL(tracker.decideRefine(pid), true)
  TEST_EQUAL(m->refine->decided, true)
  TEST_EQUAL(m->refine->winner_candidate, 1)
  TEST_EQUAL(m->update_index, 2)
  const Refine::Path& b = m->refine->candidates->assemble.candidates[1];
  TEST_EQUAL(m->identification_score, b.score)
  const std::vector<FragmentAnalysis::PTMSite> sites = Refine::ptmSites(b);
  TEST_EQUAL(m->modifications.size(), sites.size())
  for (size_t i = 0; i < sites.size() && i < m->modifications.size(); ++i)
  {
    TEST_EQUAL(m->modifications[i].candidate_start, sites[i].start_position)
    TEST_EQUAL(m->modifications[i].candidate_end, sites[i].end_position)
    TEST_REAL_SIMILAR(m->modifications[i].mass_shift, sites[i].mass_shift)
  }
  // the fragments are path b's: b6 + Y in, b4 + X out, one key per explained ion
  TEST_TRUE(m->fragments.count(FragmentKey{"b", 6}) == 1)
  TEST_TRUE(m->fragments.count(FragmentKey{"b", 4}) == 0)
  TEST_EQUAL(m->fragments.size(), b.ions_explained.size())
}
END_SECTION

/////////////////////////////////////////////////////////////
// 7. After the decision the planner is the ambiguity branch: containing fragments of the decided model's windows. No
//    record is dispatched again, and a mass cannot buy two slots.
//
// FAILS UNDER: a planner that keeps planning discriminating masses after the decision, or stops planning at all.
/////////////////////////////////////////////////////////////
START_SECTION(after_decision_planner_falls_through_to_ambiguity)
{
  const int pid = 38;
  Config cfg{trackerConfig("discriminating", "refine", 3)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedAndFinalize(cfg, tracker, pid, 608, twoTotals());
  TEST_EQUAL(tracker.planNextScans(pid).size(), 3)
  const ProteoformModel* m = tracker.getModel(pid);
  ABORT_IF(m == nullptr || !m->refine.has_value())
  const auto& recs = m->refine->records;
  tracker.recordDiscriminatingScan(pid, recs[0].target.mass, ms3Of({}), LOOSE, TIGHT);
  tracker.recordDiscriminatingScan(pid, recs[1].target.mass, fittingScan(*m, recs[1]), LOOSE, TIGHT);
  tracker.recordDiscriminatingScan(pid, recs[2].target.mass, ms3Of({}), LOOSE, TIGHT);
  TEST_EQUAL(tracker.decideRefine(pid), true)
  ABORT_IF(m->modifications.empty())
  const ModificationState& mod = m->modifications[0];
  TEST_TRUE(mod.candidate_start < mod.candidate_end) // path b's jump has a window: something left to localize

  const std::vector<Ms3Target> plan = tracker.planNextScans(pid);
  TEST_EQUAL(plan.size(), 3)
  std::set<int> nominal;
  for (const Ms3Target& t : plan)
  {
    const auto f = m->fragments.find(FragmentKey{t.ion_type, t.ion_index});
    TEST_TRUE(f != m->fragments.end())
    if (f != m->fragments.end()) TEST_TRUE(f->second.cover_start <= mod.candidate_start && f->second.cover_end >= mod.candidate_end)
    nominal.insert(SpectralDeconvolution::getNominalMass(t.frag_mass));
  }
  TEST_EQUAL(nominal.size(), plan.size())
  TEST_EQUAL(countDispatched(*m), 3)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 7a (S1). A pool that is empty with nothing awaiting a return can never complete the round through MS3: the planner
//     decides (candidate 0 stands) and returns the AMBIGUITY plan in the same call -- the same targets an ambiguity-mode
//     tracker plans for the same scan. Control: a round whose pool is exhausted while a dispatched record has not yet
//     returned keeps waiting -- no plan, not decided -- and completes on that return.
//
// FAILS UNDER: the dead end (no plan and never decided: no MS3 of any kind for the precursor); a fall-through that does
// not plan in the same call; a fall-through that fires on a partial round.
/////////////////////////////////////////////////////////////
START_SECTION(exhausted_pool_falls_through_to_ambiguity)
{
  // Every record below the fragment-charge floor: nothing was ever dispatchable.
  const int pid = 41;
  Config cfg{withFragmentChargeFloor(trackerConfig("discriminating", "refine", 3), 2)};
  TEST_EQUAL(cfg.level(2).min_charge, 2)
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedAndFinalize(cfg, tracker, pid, 611, ambiguousBest());
  const ProteoformModel* m = tracker.getModel(pid);
  ABORT_IF(m == nullptr || !m->refine.has_value())
  TEST_EQUAL(m->refine->candidates->assemble.candidates.size() >= 2, true)
  TEST_EQUAL(m->refine->records.empty(), false)
  ABORT_IF(m->modifications.empty())
  TEST_TRUE(m->modifications[0].candidate_start < m->modifications[0].candidate_end) // the window ambiguity plans for

  Config amb_cfg{withFragmentChargeFloor(trackerConfig("ambiguity", "refine", 3), 2)};
  IdaLogger amb_logger(amb_cfg);
  ProteoformTracker amb_tracker(amb_cfg, amb_logger);
  feedAndFinalize(amb_cfg, amb_tracker, pid, 611, ambiguousBest());
  const std::vector<Ms3Target> expected = amb_tracker.planNextScans(pid);
  TEST_EQUAL(expected.empty(), false)

  std::string log;
  const std::vector<Ms3Target> plan = planCapturing(tracker, pid, log);
  TEST_TRUE(log.find("[DISCRIMINATING] precursor_id=41 decided=fallthrough reason=pool_exhausted") != std::string::npos)
  TEST_EQUAL(m->refine->decided, true)
  TEST_EQUAL(m->refine->winner_candidate, 0)
  TEST_EQUAL(countDispatched(*m), 0)
  TEST_EQUAL(plan.size(), expected.size())
  for (size_t i = 0; i < plan.size() && i < expected.size(); ++i)
  {
    TEST_EQUAL(plan[i].ion_type, expected[i].ion_type)
    TEST_EQUAL(plan[i].ion_index, expected[i].ion_index)
    TEST_REAL_SIMILAR(plan[i].frag_mass, expected[i].frag_mass)
  }

  // Control: every record dispatched, one not yet returned.
  const int pid2 = 42;
  Config cfg2{trackerConfig("discriminating", "refine", 20)};
  IdaLogger logger2(cfg2);
  ProteoformTracker tracker2(cfg2, logger2);
  feedAndFinalize(cfg2, tracker2, pid2, 612, twoTotals());
  const ProteoformModel* m2 = tracker2.getModel(pid2);
  ABORT_IF(m2 == nullptr || !m2->refine.has_value())
  const size_t n = m2->refine->records.size();
  TEST_EQUAL(tracker2.planNextScans(pid2).size(), n)
  TEST_EQUAL(countDispatched(*m2), static_cast<int>(n))
  ABORT_IF(n < 2 || countDispatched(*m2) != static_cast<int>(n))
  const DeconvolvedSpectrum empty = ms3Of({});
  for (size_t i = 0; i + 1 < n; ++i) tracker2.recordDiscriminatingScan(pid2, m2->refine->records[i].target.mass, empty, LOOSE, TIGHT);
  std::string log2;
  TEST_EQUAL(planCapturing(tracker2, pid2, log2).size(), 0)
  TEST_TRUE(log2.find("objective=discriminating targets=0 reason=pool_exhausted") != std::string::npos)
  TEST_TRUE(log2.find("decided=fallthrough") == std::string::npos)
  TEST_EQUAL(m2->refine->decided, false)
  TEST_EQUAL(tracker2.recordDiscriminatingScan(pid2, m2->refine->records[n - 1].target.mass, empty, LOOSE, TIGHT), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 7b (N1). A mass only the second path explains is dispatched UNASSIGNED ('u'/0, planExhaustive_'s marking for a mass
//     with no ion identity), so the MS3 return path does not match and fold path b's claim against the winner's context;
//     an a-only mass keeps its label. The returned scan still reaches its record by mass and is tested on path b's ion.
//
// FAILS UNDER: labelling a b-only target with path b's ion (the old `ion_a ? ion_a : ion_b`); unassigning every target;
// a record lookup that keys on the label instead of the mass.
/////////////////////////////////////////////////////////////
START_SECTION(b_only_target_is_dispatched_unassigned)
{
  const int pid = 43;
  Config cfg{trackerConfig("discriminating", "refine", 3)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  feedAndFinalize(cfg, tracker, pid, 613, twoTotals());
  const std::vector<Ms3Target> plan = tracker.planNextScans(pid);
  const ProteoformModel* m = tracker.getModel(pid);
  ABORT_IF(plan.size() != 3 || m == nullptr || !m->refine.has_value())
  const auto& recs = m->refine->records;
  // records 0 and 2 are a-only (b4 + X, b5 + X), record 1 b-only (b6 + Y)
  TEST_EQUAL(label(recs[1].target.ion_a), "-")
  TEST_EQUAL(label(recs[1].target.ion_b), "b6")
  TEST_EQUAL(plan[0].ion_type, std::string("b"))
  TEST_EQUAL(plan[0].ion_index, 4)
  TEST_EQUAL(plan[1].ion_type, std::string("u"))
  TEST_EQUAL(plan[1].ion_index, 0)
  TEST_REAL_SIMILAR(plan[1].frag_mass, recs[1].target.mass)
  TEST_EQUAL(plan[2].ion_type, std::string("b"))
  TEST_EQUAL(plan[2].ion_index, 5)

  tracker.recordDiscriminatingScan(pid, recs[1].target.mass, fittingScan(*m, recs[1]), LOOSE, TIGHT);
  TEST_EQUAL(recs[1].returned, true)
  TEST_EQUAL(recs[1].fits_b, true)
  TEST_EQUAL(recs[1].fits_a, false)
}
END_SECTION

/////////////////////////////////////////////////////////////
// 8. Under refine a fragment that matches WITHOUT the shift moves no boundary; a shift-carrying one still does. The same
//    model under flashtnt moves both (control).
//
// FAILS UNDER: an ungated Pass A (refine would narrow to 4..5 like flashtnt); a gate that also drops the shift-carrying
// branch (refine would stay at 2..6).
/////////////////////////////////////////////////////////////
START_SECTION(narrowing_ignores_bare_fragments_under_refine)
{
  std::map<char, std::vector<double>> lad;
  FragmentAnalysis::computePTMAdjustedFragmentMasses(SEQ, {}, {"b"}, lad);
  const double shift = 42.010565;
  // mod on residues 2..6; b3 matches WITHOUT the shift (brackets: 2 <= 3 < 6), b5 WITH it (2 <= 5 < 6)
  FragmentAnalysis::ProteoformMatch match;
  match.score = 1.0;
  match.region_start = -1;
  match.region_end = -1;
  match.matched_protein = "synthetic";
  match.proteoform_sequence = SEQ;
  match.ptm_sites.push_back(FragmentAnalysis::PTMSite{4, 2, 6, shift});
  FragmentAnalysis::ProteoformMatch::FragmentMatch b3, b5;
  b3.ion_type = "b";
  b3.ion_index = 3;
  b3.observed_mass = lad['b'][2];
  b5.ion_type = "b";
  b5.ion_index = 5;
  b5.observed_mass = lad['b'][4] + shift;
  match.fragments = {b3, b5};

  auto narrowed = [&](const std::string& method, bool with_b5) {
    Config cfg{trackerConfig("ambiguity", method, 3)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    FragmentAnalysis::ProteoformMatch mm = match;
    if (!with_b5) mm.fragments.pop_back();
    DeconvolvedSpectrum d(700);
    for (const auto& f : mm.fragments) d.push_back(makePeak(f.observed_mass, 100.0));
    Ms2Params p;
    p.activation_type = "HCD";
    tracker.feedScan(1, 2, p, 700, d, mm, 1.0, makeMs2Ctx());
    tracker.finalizeMS2(1);
    const ModificationState s = tracker.getModel(1)->modifications.at(0);
    return std::make_pair(s.candidate_start, s.candidate_end);
  };

  // the bare b3 alone
  TEST_EQUAL(narrowed("flashtnt", false).first, 4) // control: 2..6 -> 4..6
  TEST_EQUAL(narrowed("flashtnt", false).second, 6)
  TEST_EQUAL(narrowed("refine", false).first, 2)   // refine: unchanged
  TEST_EQUAL(narrowed("refine", false).second, 6)
  // with the shift-carrying b5
  TEST_EQUAL(narrowed("flashtnt", true).first, 4)
  TEST_EQUAL(narrowed("flashtnt", true).second, 5)
  TEST_EQUAL(narrowed("refine", true).first, 2)    // only b5 moved a boundary
  TEST_EQUAL(narrowed("refine", true).second, 5)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
