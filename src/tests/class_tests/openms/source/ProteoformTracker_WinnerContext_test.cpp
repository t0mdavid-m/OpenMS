// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// buildWinnerProteoformContext: renders the LIVE winner model as an MS3FragmentMatcher::ProteoformContext
// so MS3 is scored against the WINNER proteoform (not the triggering scan). Pins:
//   (1) no finalized/non-empty winner  -> EMPTY context (region -1/-1, no PTMs) => MS3 matches nothing;
//   (2) full-region winner             -> region 0..P, ALL mods (localized AND ambiguous) as
//                                         region-1-based PTMSites (localized point PTMs are NOT skipped);
//   (3) truncated-region winner        -> resolved region bounds + PTM start/end rebased by ws.
//
// CONSTRUCTION NOTE: each winner match carries NO fragments, so finalize maps nothing and
// narrowModifications_ returns early (m.fragments empty) -> the modifications stay EXACTLY as seeded,
// isolating the region/PTM conversion under test from any narrowing.

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms2Params.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  const char* WINNER_SEQ = "PEPTIDEK";   // P = 8

  // Minimal characterization config; no runtime log paths (IdaLogger writes nothing).
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
    "characterization": { "objective": "coverage", "protein_sequence": "PEPTIDEK" },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "intensity", "max_targets": 10 },
      "ms3": { "selection": "intensity", "max_targets": 10 }
    }
  })";

  FragmentAnalysis::PTMSite ptm(int start, int end, double shift)
  {
    FragmentAnalysis::PTMSite s;
    s.start_position = start;
    s.end_position = end;
    s.position = (start + end) / 2;
    s.mass_shift = shift;
    return s;
  }

  // A winner ProteoformMatch: score >= 0 + non-empty sequence (finalize accepts it as winner), the given
  // region (0-based; -1 => full), the given PTM sites (FULL-PROTEIN 1-based, as m.modifications expects),
  // and NO fragments (so narrowModifications_ no-ops and the seeded mods are inspected verbatim).
  FragmentAnalysis::ProteoformMatch makeWinner(int region_start, int region_end,
      const std::vector<FragmentAnalysis::PTMSite>& ptms)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.score = 1.0;
    m.region_start = region_start;
    m.region_end = region_end;
    m.matched_protein = "synthetic";
    m.proteoform_sequence = WINNER_SEQ;
    m.ptm_sites = ptms;
    return m;
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
}

START_TEST(ProteoformTracker_WinnerContext, "$Id$")

/////////////////////////////////////////////////////////////
// (1) No finalized / non-empty winner -> empty context.
/////////////////////////////////////////////////////////////
START_SECTION(no_winner_yields_empty_context)
{
  Config cfg{std::string(tracker_config)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);

  // Unknown precursor: no model -> empty context.
  MS3FragmentMatcher::ProteoformContext c0 = tracker.buildWinnerProteoformContext(999);
  TEST_EQUAL(c0.region_start, -1)
  TEST_EQUAL(c0.region_end, -1)
  TEST_EQUAL((int)c0.ptm_sites.size(), 0)

  // Finalized but NO winner (fed match scores < 0): model stays empty (no proteoform) -> empty context.
  Ms2Params p; p.activation_type = "HCD";
  ScanCommand ms2_ctx = makeMs2Ctx();
  DeconvolvedSpectrum d(101);
  FragmentAnalysis::ProteoformMatch noid;
  noid.score = -1.0;               // not a winner candidate
  noid.proteoform_sequence = "";
  tracker.feedScan(5, 2, p, 101, d, noid, -1.0, ms2_ctx);
  tracker.finalizeMS2(5);
  MS3FragmentMatcher::ProteoformContext c1 = tracker.buildWinnerProteoformContext(5);
  TEST_EQUAL(c1.region_start, -1)
  TEST_EQUAL(c1.region_end, -1)
  TEST_EQUAL((int)c1.ptm_sites.size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// (2) Full-region winner: region 0..P, BOTH localized and ambiguous mods present.
/////////////////////////////////////////////////////////////
START_SECTION(full_region_includes_localized_and_ambiguous_mods)
{
  Config cfg{std::string(tracker_config)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);

  Ms2Params p; p.activation_type = "HCD";
  ScanCommand ms2_ctx = makeMs2Ctx();
  DeconvolvedSpectrum d(101);
  // full region (-1,-1); localized point PTM @2 + ambiguous PTM @[4,6].
  std::vector<FragmentAnalysis::PTMSite> mods = { ptm(2, 2, 79.9663), ptm(4, 6, 42.0106) };
  tracker.feedScan(7, 2, p, 101, d, makeWinner(-1, -1, mods), 1.0, ms2_ctx);
  tracker.finalizeMS2(7);

  MS3FragmentMatcher::ProteoformContext c = tracker.buildWinnerProteoformContext(7);
  TEST_EQUAL(c.region_start, 0)              // full region resolved to 0
  TEST_EQUAL(c.region_end, 8)                // .. and P (8)
  TEST_EQUAL((int)c.ptm_sites.size(), 2)     // BOTH mods present -- localized NOT skipped
  ABORT_IF((int)c.ptm_sites.size() != 2)
  // region-1-based == full since ws = 0
  TEST_EQUAL(c.ptm_sites[0].start_position, 2)
  TEST_EQUAL(c.ptm_sites[0].end_position, 2)     // localized point PTM preserved (start == end)
  TEST_TRUE(std::abs(c.ptm_sites[0].mass_shift - 79.9663) < 1e-6)
  TEST_EQUAL(c.ptm_sites[1].start_position, 4)
  TEST_EQUAL(c.ptm_sites[1].end_position, 6)     // ambiguous range preserved
  TEST_TRUE(std::abs(c.ptm_sites[1].mass_shift - 42.0106) < 1e-6)
}
END_SECTION

/////////////////////////////////////////////////////////////
// (3) Truncated-region winner: bounds resolved, PTM start/end rebased by ws.
/////////////////////////////////////////////////////////////
START_SECTION(truncated_region_rebases_ptm_by_ws)
{
  Config cfg{std::string(tracker_config)};
  IdaLogger logger(cfg);
  ProteoformTracker tracker(cfg, logger);

  Ms2Params p; p.activation_type = "HCD";
  ScanCommand ms2_ctx = makeMs2Ctx();
  DeconvolvedSpectrum d(101);
  // truncated region [2,7) (ws = 2); localized mod at FULL-PROTEIN residue 4.
  std::vector<FragmentAnalysis::PTMSite> mods = { ptm(4, 4, 42.0106) };
  tracker.feedScan(9, 2, p, 101, d, makeWinner(2, 7, mods), 1.0, ms2_ctx);
  tracker.finalizeMS2(9);

  MS3FragmentMatcher::ProteoformContext c = tracker.buildWinnerProteoformContext(9);
  TEST_EQUAL(c.region_start, 2)
  TEST_EQUAL(c.region_end, 7)
  TEST_EQUAL((int)c.ptm_sites.size(), 1)
  ABORT_IF((int)c.ptm_sites.size() != 1)
  TEST_EQUAL(c.ptm_sites[0].start_position, 2)   // full-protein 4 - ws 2 = region-1-based 2
  TEST_EQUAL(c.ptm_sites[0].end_position, 2)
}
END_SECTION

END_TEST
