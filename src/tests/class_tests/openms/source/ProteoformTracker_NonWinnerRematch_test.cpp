// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Winner-anchored pooling (Change 1): a NON-winner MS2 scan no longer pools its own already-matched
// fragments; instead its RAW deconvolved masses are re-matched against the WINNER ladder. Only masses
// that land uniquely on a winner ion become fragments. This pins the five matching rules:
//
//   (a) no-match            -> the mass is DROPPED (never a fragment).
//   (b) new position        -> a real winner ion the winner didn't have ADDS a fragment (coverage up).
//   (c) shared position     -> a mass at a position the winner already has MERGES into the SAME key
//                              (n_ms2 == 2, one MappedFragment). [best_ms2 = max-intensity across scans
//                              is the SAME upsert path; it can't be driven here because hand-built
//                              PeakGroups have getChargeIntensity()==0, so with equal (0) intensities
//                              the tie keeps the winner's obs. The real intensity-driven selection is
//                              exercised by the C# golden AssertMs3Stage0CeNotCollapsed.]
//   (d) base AND base+shift -> a bracketing mass within tol of BOTH the mod-out and mod-in theoretical
//                              is DROPPED as ambiguous (constructed with a deliberately tiny shift).
//   (e) >=2 distinct ions   -> a mass within tol of two DISTINCT winner ions is assigned to the CLOSEST
//                              by ppm (constructed via a crafted ambiguous-mod shift so b4-with-shift
//                              lands ~5 ppm from b5-base; the peak sits exactly on b5 -> b5 wins).
//
// All ladders are computed in-test with the SAME FragmentAnalysis::computePTMAdjustedFragmentMasses the
// engine uses, so fed peak masses match the winner ladder at 0 ppm (or the crafted offset). Requires
// OPENMS_DATA_PATH (ResidueDB), as CI sets for ctest.

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

#include <cmath>
#include <cstring>
#include <map>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  const char* WINNER_SEQ = "PEPTIDEK";   // P = 8

  const char* tracker_config = R"({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": 0.9,
    "min_charge": 1,
    "max_charge": 50,
    "min_mass": 100,
    "max_mass": 50000,
    "tol": [
      10,
      10,
      10
    ]
  },
  "flashtnt": {
    "min_length": 3,
    "max_length": 8,
    "max_ptm_count": 3,
    "max_flanking_mass_diff": 50000
  },
  "faims": {
    "cv_values": [],
    "max_cv_skip": 0,
    "cv_precursor_threshold": 15
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false,
      "value_ms": 60000
    },
    "scan_timeout": {
      "enabled": false,
      "value_ms": 30000
    }
  },
  "files": {
    "target_logs": [],
    "fasta": "",
    "inclusion_list": "",
    "ptm_list": ""
  },
  "conditional_ms2": false,
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3
  },
  "characterization": {
    "mode": "coverage",
    "protein_sequence": "PEPTIDEK",
    "max_targets": 10
  },
  "ms_settings": {
    "ms1": {
      "analyzer": "Orbitrap",
      "first_mass": 100,
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
    },
    "ms3": {
      "analyzer": "Orbitrap",
      "activation": "CID",
      "collision_energy": 25,
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

  FragmentAnalysis::PTMSite ptmSite(int start, int end, double shift)
  {
    FragmentAnalysis::PTMSite s;
    s.start_position = start;
    s.end_position = end;
    s.position = (start + end) / 2;
    s.mass_shift = shift;
    return s;
  }

  ScanCommand makeMs2Ctx()
  {
    ScanCommand ctx{};
    ctx.scan_id = 1;
    ctx.msn_level = 2;
    ctx.num_stages = 1;
    std::strncpy(ctx.stages[0].activation_type, "HCD", sizeof(ctx.stages[0].activation_type) - 1);
    return ctx;
  }

  // Feed the WINNER MS2 scan (fed first, score 1.0). Its single fragment is pooled verbatim (winner path);
  // ptms seeds the model's modifications.
  void feedWinner(ProteoformTracker& t, int pid, int scan_id, const std::string& ion_type, int ion_index,
      double observed_mass, const std::vector<FragmentAnalysis::PTMSite>& ptms, const ScanCommand& ctx)
  {
    Ms2Params p; p.collision_energy = 20.0; p.activation_type = "HCD"; p.reaction_time = 0.0;
    DeconvolvedSpectrum d(scan_id);
    d.push_back(makeSyntheticPeakGroup(observed_mass / 2.0 + 1.0, observed_mass, 2));
    FragmentAnalysis::ProteoformMatch m;
    m.score = 1.0;
    m.region_start = -1; m.region_end = -1;
    m.matched_protein = "synthetic";
    m.proteoform_sequence = WINNER_SEQ;
    FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
    fm.ion_type = ion_type; fm.ion_index = ion_index; fm.observed_mass = observed_mass;
    m.fragments.push_back(fm);
    m.ptm_sites = ptms;
    t.feedScan(pid, 2, p, scan_id, d, m, 1.0, ctx);
  }

  // Feed a NON-winner MS2 scan (score 0.5 < winner). Only its raw peak masses are used (re-matched
  // against the winner ladder); the match's fragments are ignored.
  void feedNonWinner(ProteoformTracker& t, int pid, int scan_id, const std::vector<double>& peak_masses,
      const ScanCommand& ctx)
  {
    Ms2Params p; p.collision_energy = 35.0; p.activation_type = "HCD"; p.reaction_time = 0.0;
    DeconvolvedSpectrum d(scan_id);
    for (double mm : peak_masses) d.push_back(makeSyntheticPeakGroup(mm / 2.0 + 1.0, mm, 2));
    FragmentAnalysis::ProteoformMatch m;
    m.score = 0.5;                 // valid winner candidate but loses to the 1.0 winner
    m.region_start = -1; m.region_end = -1;
    m.proteoform_sequence = WINNER_SEQ;
    t.feedScan(pid, 2, p, scan_id, d, m, 0.5, ctx);
  }

  bool hasFrag(const ProteoformModel* mdl, const std::string& type, int idx)
  {
    return mdl->fragments.find(FragmentKey{type, idx}) != mdl->fragments.end();
  }
}

START_TEST(ProteoformTracker_NonWinnerRematch, "$Id$")

/////////////////////////////////////////////////////////////
// (a) A non-winner mass off the winner ladder is DROPPED.
/////////////////////////////////////////////////////////////
START_SECTION(off_ladder_mass_is_dropped)
{
  Config cfg{std::string(tracker_config)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  ScanCommand ctx = makeMs2Ctx();

  std::map<char, std::vector<double>> lad;
  FragmentAnalysis::computePTMAdjustedFragmentMasses(WINNER_SEQ, {}, {"b", "y"}, lad);

  feedWinner(tracker, 7, 101, "b", 6, lad['b'][5], {}, ctx);   // winner: b6
  feedNonWinner(tracker, 7, 102, {99999.0}, ctx);              // off-ladder mass
  tracker.finalizeMS2(7);

  const ProteoformModel* mdl = tracker.getModel(7);
  ABORT_IF(mdl == nullptr)
  TEST_TRUE(hasFrag(mdl, "b", 6))               // winner fragment kept
  TEST_EQUAL((int)mdl->fragments.size(), 1)     // the off-ladder mass added NOTHING
}
END_SECTION

/////////////////////////////////////////////////////////////
// (b) A non-winner mass at a NEW winner position adds a fragment.
/////////////////////////////////////////////////////////////
START_SECTION(new_position_adds_coverage)
{
  Config cfg{std::string(tracker_config)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  ScanCommand ctx = makeMs2Ctx();

  std::map<char, std::vector<double>> lad;
  FragmentAnalysis::computePTMAdjustedFragmentMasses(WINNER_SEQ, {}, {"b", "y"}, lad);

  feedWinner(tracker, 7, 101, "b", 6, lad['b'][5], {}, ctx);   // winner: b6
  feedNonWinner(tracker, 7, 102, {lad['y'][2]}, ctx);          // real y3 -> NEW position
  tracker.finalizeMS2(7);

  const ProteoformModel* mdl = tracker.getModel(7);
  ABORT_IF(mdl == nullptr)
  TEST_TRUE(hasFrag(mdl, "b", 6))
  TEST_TRUE(hasFrag(mdl, "y", 3))               // added by the non-winner re-match
  TEST_EQUAL((int)mdl->fragments.size(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
// (c) A non-winner mass at a SHARED winner position merges into ONE key.
/////////////////////////////////////////////////////////////
START_SECTION(shared_position_merges_into_one_key)
{
  Config cfg{std::string(tracker_config)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  ScanCommand ctx = makeMs2Ctx();

  std::map<char, std::vector<double>> lad;
  FragmentAnalysis::computePTMAdjustedFragmentMasses(WINNER_SEQ, {}, {"b", "y"}, lad);

  feedWinner(tracker, 7, 101, "b", 6, lad['b'][5], {}, ctx);   // winner: b6
  feedNonWinner(tracker, 7, 102, {lad['b'][5]}, ctx);          // real b6 again -> SHARED position
  tracker.finalizeMS2(7);

  const ProteoformModel* mdl = tracker.getModel(7);
  ABORT_IF(mdl == nullptr)
  TEST_EQUAL((int)mdl->fragments.size(), 1)     // ONE key, not a duplicate
  TEST_TRUE(hasFrag(mdl, "b", 6))
  const MappedFragment& f = mdl->fragments.at(FragmentKey{"b", 6});
  TEST_EQUAL(f.n_ms2, 2)                          // both observations merged onto the same key
  TEST_TRUE(f.best_ms2.has_value())
  ABORT_IF(!f.best_ms2.has_value())
  // Equal (0) intensities -> the strictly-greater tie keeps the winner's (first) observation.
  // The intensity-driven "higher wins" selection is the same path, covered by the C# golden.
  TEST_EQUAL(f.best_ms2->source_scan_id, 101)
}
END_SECTION

/////////////////////////////////////////////////////////////
// (d) A bracketing mass within tol of BOTH base and base+shift is DROPPED (ambiguous).
/////////////////////////////////////////////////////////////
START_SECTION(bracketing_base_and_shift_both_match_dropped)
{
  Config cfg{std::string(tracker_config)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  ScanCommand ctx = makeMs2Ctx();

  // Ambiguous mod over [3,5] with a DELIBERATELY TINY shift (0.005 Da) so that, for a bracketing ion,
  // base and base+shift both fall inside the 10 ppm tolerance of a mass placed between them.
  const double TINY = 0.005;
  std::vector<FragmentAnalysis::PTMSite> mods = { ptmSite(3, 5, TINY) };
  std::map<char, std::vector<double>> lad;
  FragmentAnalysis::computePTMAdjustedFragmentMasses(WINNER_SEQ, mods, {"b", "y"}, lad);

  // b7 covers [1,7] -> CONTAINS [3,5] but does NOT bracket it (so narrowModifications_ leaves the mod).
  feedWinner(tracker, 7, 101, "b", 7, lad['b'][6], mods, ctx);
  // b4 (cleave after 4) BRACKETS [3,5]; place the mass halfway between b4-base and b4-base+shift.
  feedNonWinner(tracker, 7, 102, {lad['b'][3] + TINY / 2.0}, ctx);
  tracker.finalizeMS2(7);

  const ProteoformModel* mdl = tracker.getModel(7);
  ABORT_IF(mdl == nullptr)
  TEST_TRUE(hasFrag(mdl, "b", 7))               // winner fragment kept
  TEST_TRUE(!hasFrag(mdl, "b", 4))              // ambiguous double-match DROPPED (dodgy)
}
END_SECTION

/////////////////////////////////////////////////////////////
// (e) A mass within tol of TWO distinct winner ions is assigned to the CLOSEST by ppm.
/////////////////////////////////////////////////////////////
START_SECTION(two_distinct_ions_closest_ppm_wins)
{
  Config cfg{std::string(tracker_config)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  ScanCommand ctx = makeMs2Ctx();

  // Bare ladder (no mod) to size the crafted shift.
  std::map<char, std::vector<double>> bare;
  FragmentAnalysis::computePTMAdjustedFragmentMasses(WINNER_SEQ, {}, {"b", "y"}, bare);
  const double b4 = bare['b'][3];
  const double b5 = bare['b'][4];
  const double offset = b5 * 5.0e-6;            // ~5 ppm (< 10 ppm tol)
  const double DELTA = (b5 - b4) - offset;      // so b4-base+DELTA = b5 - offset (~5 ppm below b5)

  // Ambiguous mod over [3,6]: both b4 and b5 bracket it, so both carry a base / base+shift candidate.
  std::vector<FragmentAnalysis::PTMSite> mods = { ptmSite(3, 6, DELTA) };

  // Winner y7 covers [2,8] -> contains [3,6], does NOT bracket it.
  std::map<char, std::vector<double>> lad;
  FragmentAnalysis::computePTMAdjustedFragmentMasses(WINNER_SEQ, mods, {"b", "y"}, lad);
  feedWinner(tracker, 7, 101, "y", 7, lad['y'][6], mods, ctx);

  // Peak sits EXACTLY on b5-base (0 ppm) and ~5 ppm from b4-base+shift -> closest-ppm must pick b5.
  feedNonWinner(tracker, 7, 102, {b5}, ctx);
  tracker.finalizeMS2(7);

  const ProteoformModel* mdl = tracker.getModel(7);
  ABORT_IF(mdl == nullptr)
  TEST_TRUE(hasFrag(mdl, "y", 7))               // winner fragment kept
  TEST_TRUE(hasFrag(mdl, "b", 5))               // closest ion (0 ppm) wins
  TEST_TRUE(!hasFrag(mdl, "b", 4))              // farther ion (~5 ppm) loses
}
END_SECTION

END_TEST
