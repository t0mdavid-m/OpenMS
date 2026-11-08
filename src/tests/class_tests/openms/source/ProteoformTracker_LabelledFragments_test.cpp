// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// ADR-0064: ProteoformTracker::labelledFragments, the pool a quantification MS3 is chosen from. A labelled fragment
// is a mapped fragment, of ANY ion type (D-C), with a best-MS2 observation, shorter than the region (never the
// full-length ion) and spanning at least one cysteine. The window-SNR filter, the budget and the outcome live in
// Exploration::dispatchQuant and are pinned end to end in FLASHIda_ProcessScan_test on the real iodoTMT fixture.
// Also ScanCommandQueue::buildMS3's role marker, which makes a quantification MS3 'Q'.
//
// Construction follows ProteoformTracker_CEOptimization_test: one synthetic winner scan whose fragments are pooled
// verbatim. Hand-built PeakGroups carry getChargeIntensity 0, so every candidate ties on intensity here and the order
// is the documented tie-break (ion type, then index); the intensity ranking itself is pinned on real data.

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms2Params.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanRole.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <cstring>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;

namespace
{
  // No runtime block / no log files -> IdaLogger writes nothing.
  const char* kConfig = R"({
  "deconvolution": { "min_charge": 1, "max_charge": 50, "min_mass": 100, "max_mass": 50000, "tol": [10, 10, 10] },
  "precursor_selection": { "rank_by": "qscore", "max_precursors": 3 },
  "identification": { "reference": { "sequence": "PEPTCDEK" } },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000 },
    "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 30 }
  }
})";

  PeakGroup peakAt(double mono_mass, int charge)
  {
    PeakGroup pg(charge, charge, true);
    pg.setMonoisotopicMass(mono_mass);
    FLASHHelperClasses::LogMzPeak lp;
    lp.mz = mono_mass / charge + 1.007276;
    lp.abs_charge = charge;
    lp.intensity = 1000.0f;
    pg.push_back(lp);
    return pg;
  }

  ScanCommand ms2Ctx()
  {
    ScanCommand ctx{};
    ctx.scan_id = 1;
    ctx.msn_level = 2;
    ctx.num_stages = 1;
    ctx.mono_mass = 900.0;
    ctx.stages[0].precursor_mz = 451.0;
    ctx.stages[0].isolation_width = 2.0;
    ctx.stages[0].charge_state = 2;
    ctx.stages[0].collision_energy = 30;
    std::strncpy(ctx.stages[0].activation_type, "HCD", sizeof(ctx.stages[0].activation_type) - 1);
    return ctx;
  }

  // Feed one winner scan of `seq` carrying the matched ions `ions` (type, index) at their theoretical masses, finalize,
  // and return the tracker's labelled fragments.
  std::vector<Ms3Target> labelledOf(const std::string& seq, const std::vector<std::pair<std::string, int>>& ions,
                                    bool* has_cys)
  {
    Config cfg{std::string(kConfig)};
    IdaLogger logger(cfg);
    ProteoformTracker tracker(cfg, logger);
    std::map<char, std::vector<double>> ladder;
    FragmentAnalysis::computePTMAdjustedFragmentMasses(seq, {}, {"b", "y", "c"}, ladder);

    FragmentAnalysis::ProteoformMatch m;
    m.score = 1.0;
    m.identification_score = 1.0;
    m.region_start = -1;
    m.region_end = -1;
    m.matched_protein = "synthetic";
    m.proteoform_sequence = seq;
    DeconvolvedSpectrum d(101);
    for (const auto& [type, index] : ions)
    {
      const double mass = ladder[type[0]][index - 1];
      FragmentAnalysis::ProteoformMatch::FragmentMatch fm;
      fm.ion_type = type;
      fm.ion_index = index;
      fm.observed_mass = mass;
      m.fragments.push_back(fm);
      d.push_back(peakAt(mass, 1));
    }
    d.sort();
    Ms2Params p;
    p.collision_energy = 30;
    p.activation_type = "HCD";
    tracker.feedScan(7, 2, p, 101, d, m, 1.0, ms2Ctx());
    tracker.finalizeMS2(7);
    return tracker.labelledFragments(7, 'C', has_cys);
  }

  bool has(const std::vector<Ms3Target>& ts, const std::string& type, int index)
  {
    for (const auto& t : ts) if (t.ion_type == type && t.ion_index == index) return true;
    return false;
  }
}

START_TEST(ProteoformTracker_LabelledFragments, "$Id$")

START_SECTION(only_fragments_spanning_a_cysteine_any_ion_type_never_full_length)
{
  // PEPTCDEK, L = 8, C at residue 5.  b4 PEPT (no C), b5 PEPTC, c6 PEPTCD, y4 CDEK, y3 DEK (no C), b8 the whole region.
  bool has_cys = false;
  const auto ts = labelledOf("PEPTCDEK", {{"b", 4}, {"b", 5}, {"c", 6}, {"y", 4}, {"y", 3}, {"b", 8}}, &has_cys);
  TEST_EQUAL(has_cys, true)
  TEST_EQUAL(ts.size(), 3)
  TEST_EQUAL(has(ts, "b", 5), true)
  TEST_EQUAL(has(ts, "c", 6), true)    // D-C: a c ion is a labelled fragment too
  TEST_EQUAL(has(ts, "y", 4), true)
  TEST_EQUAL(has(ts, "b", 4), false)   // spans no cysteine
  TEST_EQUAL(has(ts, "y", 3), false)
  TEST_EQUAL(has(ts, "b", 8), false)   // the full-length ion is the precursor (ADR-0052)
  // Equal intensities: ion type, then index.
  if (ts.size() == 3)
  {
    TEST_STRING_EQUAL(ts[0].ion_type, "b")
    TEST_STRING_EQUAL(ts[1].ion_type, "c")
    TEST_STRING_EQUAL(ts[2].ion_type, "y")
  }
  // Each target carries its best-MS2 observation: a positive isolation m/z and charge.
  for (const auto& t : ts)
  {
    TEST_EQUAL(t.frag_mz > 0.0, true)
    TEST_EQUAL(t.frag_charge, 1)
  }
}
END_SECTION

START_SECTION(no_cysteine_in_the_region_and_no_fragment_reaching_it)
{
  bool has_cys = true;
  const auto none = labelledOf("PEPTIDEK", {{"b", 4}, {"y", 4}}, &has_cys);
  TEST_EQUAL(has_cys, false)   // dispatchQuant logs no_cysteine
  TEST_EQUAL(none.empty(), true)

  // PEPTIDEKCAR: C at residue 9, outside b4 (PEPT) and y2 (AR).
  const auto unreached = labelledOf("PEPTIDEKCAR", {{"b", 4}, {"y", 2}}, &has_cys);
  TEST_EQUAL(has_cys, true)    // dispatchQuant logs no_labelled_fragment
  TEST_EQUAL(unreached.empty(), true)

  // No model at all.
  Config cfg{std::string(kConfig)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);
  TEST_EQUAL(tracker.labelledFragments(99, 'C', &has_cys).empty(), true)
  TEST_EQUAL(has_cys, false)
}
END_SECTION

// ScanCommandQueue::buildMS3's role marker: 'R' unless told otherwise (every characterization call site, unchanged),
// 'Q' for a quantification MS3 -- the one character ScanRole decodes.
START_SECTION(buildms3_marker_defaults_to_R_and_honours_Q)
{
  Config cfg{std::string(kConfig)};
  ScanCommandQueue queue(cfg);
  const ScanCommand ctx = ms2Ctx();
  ScanConfig ms3;
  ms3.activation = "HCD";
  ms3.collision_energy = 80;
  const ScanCommand r = queue.buildMS3(ctx, ms3, 600.0, 2, 2.0, ctx.scan_id, 'b', 5);
  const ScanCommand q = queue.buildMS3(ctx, ms3, 600.0, 2, 2.0, ctx.scan_id, 'b', 5, 1, {}, nullptr, {}, nullptr, 'Q');
  TEST_EQUAL(r.scan_description[3], 'R')
  TEST_EQUAL(q.scan_description[3], 'Q')
  TEST_EQUAL(roleOf(q) == ScanRole::Quantification, true)
  TEST_EQUAL(roleOf(r) == ScanRole::Identification, true)
  // Everything but the marker (and the minted id) is the same command.
  TEST_REAL_SIMILAR(q.stages[1].precursor_mz, r.stages[1].precursor_mz)
  TEST_REAL_SIMILAR(q.stages[1].collision_energy, 80.0)
  TEST_STRING_EQUAL(std::string(q.scan_description + 4), std::string(r.scan_description + 4))
}
END_SECTION

END_TEST
