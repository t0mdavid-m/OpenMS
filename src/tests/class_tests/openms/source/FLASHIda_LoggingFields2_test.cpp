// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// The second half of FLASHIda_LoggingFields_test's sections, split off on 2026-11-10 so that ctest -j runs the two
// halves at once (the whole suite was the cpp-tests job's floor). Same division of labour and coverage notes as
// FLASHIda_LoggingFields_test.cpp; helpers in FLASHIda_LoggingFields_TestHelpers.h. A section that another file
// names as FLASHIda_LoggingFields_test::<section> stayed there; the sections here were moved byte for byte.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include "FLASHIda_TestHelpers.h"
#include "FLASHIda_LoggingFields_TestHelpers.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;

START_TEST(FLASHIda_LoggingFields2, "$Id$")

// The faims_enabled column earns its recapture cost only if it says something faims_cv cannot.
// It does, and the [0] case below is the proof: ADR-0012 made cv_values: [0] mean "FAIMS on at
// compensation voltage 0", which logs faims_cv = 0 -- byte-identical to a run with no FAIMS at
// all. Before this column the two were indistinguishable in every log stream.
//
// Drives the same acquisition three times, changing ONLY cv_values, so the column is the only
// thing that can differ.
START_SECTION(commands_faims_enabled_column)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  // Each invocation builds its OWN engine, so each takes its own log folder (tag).
  auto driveWithCvs = [&](const std::string& cv_list, const std::string& tag) {
    const std::string dir = freshLogDir(tag);
    std::string json = buildJsonWithLogDir(dir, true);
    const std::string needle = "\"cv_values\": []";
    const size_t p = json.find(needle);
    if (p != std::string::npos) json.replace(p, needle.size(), "\"cv_values\": " + cv_list);
    FLASHIda ida(const_cast<char*>(json.c_str()));
    runFullCycle(&ida, ms1, ms2);
    return TSVFile::parse(dir + "/scan_commands.tsv");
  };

  // Empty cv_values -> FAIMS off.
  {
    auto t = driveWithCvs("[]", "lf_faims_off");
    ABORT_IF(t.rows.empty())
    bool enabled_ok = true, cv_ok = true;
    for (const auto& row : t.rows)
    {
      enabled_ok = enabled_ok && (cell(t, row, "faims_enabled") == "0");
      cv_ok      = cv_ok      && (std::abs(toD(cell(t, row, "faims_cv"))) < 1e-9);
    }
    TEST_TRUE(enabled_ok)
    TEST_TRUE(cv_ok)
  }

  // One CV -> FAIMS on at a fixed voltage. Under the pre-ADR-0012 `size() > 1` rule this run
  // reported no FAIMS and logged faims_cv = 0; both columns now say otherwise.
  {
    auto t = driveWithCvs("[-45]", "lf_faims_on");
    ABORT_IF(t.rows.empty())
    bool enabled_ok = true, cv_ok = true;
    for (const auto& row : t.rows)
    {
      enabled_ok = enabled_ok && (cell(t, row, "faims_enabled") == "1");
      cv_ok      = cv_ok      && (std::abs(toD(cell(t, row, "faims_cv")) + 45.0) < 1e-9);
    }
    TEST_TRUE(enabled_ok)
    TEST_TRUE(cv_ok)
  }

  // THE case the column exists for: FAIMS on, at CV 0. faims_cv is 0 here, exactly as in the
  // FAIMS-off run above -- so faims_enabled is the only column that tells them apart.
  {
    auto t = driveWithCvs("[0]", "lf_faims_cv0");
    ABORT_IF(t.rows.empty())
    bool enabled_ok = true, cv_is_zero = true;
    for (const auto& row : t.rows)
    {
      enabled_ok = enabled_ok && (cell(t, row, "faims_enabled") == "1");
      cv_is_zero = cv_is_zero && (std::abs(toD(cell(t, row, "faims_cv"))) < 1e-9);
    }
    TEST_TRUE(enabled_ok)   // ISSUE(pre-ADR-0012): reported 0 -- indistinguishable from FAIMS off
    TEST_TRUE(cv_is_zero)   // and faims_cv alone genuinely cannot disambiguate it
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// §R4 -- scan_results.tsv : MS3-normal 2-stage CE/act/rt + proteoform/protein/fragment/tic
/////////////////////////////////////////////////////////////
START_SECTION(results_ms3_normal_columns)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_r4");
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto cycle = runFullCycle(&ida, ms1, ms2);
  TEST_TRUE(cycle.ms3_cmds.size() > 0)

  auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
  auto res = TSVFile::parse(dir + "/scan_results.tsv");
  auto level = commandLevels(cmds);

  bool found_ms3 = false;
  for (const auto& row : res.rows)
  {
    std::string tid = cell(res, row, "tracking_id");
    auto it = level.find(tid);
    if (it == level.end() || it->second != 3) continue;  // MS3 rows only
    found_ms3 = true;

    // terminal scan: no children
    TEST_EQUAL(cell(res, row, "commands_pushed"), std::string("0"))
    TEST_TRUE(cell(res, row, "child_ids").empty())
    TEST_TRUE(isTrackingId(cell(res, row, "parent_tracking_id")))  // == its MS2

    // 2-stage CE/act/rt (MS2 isolation ; MS3 fragmentation). FAIL-CLOSED: the per-stage split MUST
    // yield exactly two tokens on every MS3 result row (a regression that logs a 1-token CE/act/rt
    // -- e.g. the Issue-1 MS3-exploration single-stage bug, or num_stages<2 -> "0" -- now FAILS here
    // instead of silently skipping the CE/act/rt checks).
    auto ce  = splitTokens(cell(res, row, "collision_energy"), ';');
    auto act = splitTokens(cell(res, row, "activation_type"), ';');
    auto rt  = splitTokens(cell(res, row, "reaction_time"), ';');
    TEST_EQUAL(ce.size(), 2)
    TEST_EQUAL(act.size(), 2)
    TEST_EQUAL(rt.size(), 2)
    TEST_TRUE(posFinite(toD(ce[0])) && posFinite(toD(ce[1])))
    TEST_TRUE(std::abs(toD(ce[0]) - 29.0) < 1.0)  // MS2 isolation CE (buildJsonWithLogDir fixed MS2 CE)
    TEST_TRUE(std::abs(toD(ce[1]) - 35.0) < 1.0)  // MS3 fragmentation CE
    TEST_TRUE(inActivationSet(act[0]) && inActivationSet(act[1]))
    TEST_TRUE(nonNegFinite(toD(rt[0])) && nonNegFinite(toD(rt[1])))
    // deconv_qscores reaches MS3 too -- it is one value per deconvolved mass of THIS scan, so unlike the
    // 2-stage scalars above it never gains the ';'-stage axis; it stays index-aligned with deconv_masses.
    // Unguarded on purpose: an MS3 that deconvolved nothing logs both cells empty, and 0 == 0 still holds.
    TEST_EQUAL(splitTokens(cell(res, row, "deconv_qscores"), ';').size(),
               splitTokens(cell(res, row, "deconv_masses"), ';').size())
    // NOTE: the MS3 identification payload (matched_protein/proteoform_sequence/tic_coverage/fragment_count/
    // tag_count) was removed from scan_results in the slim-down. The MS3 fragment proteoform now lives on
    // scan_commands.ms3_proteoform (locked in §T9/§F3) and tic on identification.tsv; this section retains the
    // match-independent scan_results checks (2-stage CE/act/rt + terminal structure) on every MS3 row.
  }
  TEST_TRUE(found_ms3)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §R5 -- scan_results.tsv : deconv_qscores logs THE sort key, so an MS1 row's list is non-increasing.
//        PrecursorSelection::filterPeakGroupsUsingMassExclusion_ orders the survey spectrum with
//        DeconvolvedSpectrum::sortByQscore(), which compares PeakGroup::getQscore() verbatim, and the
//        writer walks that same spectrum in place -- so a monotone list is the observable consequence of
//        the column holding that exact score. §R1 pins the column's SHAPE (one value per mass, in [0,1]);
//        this pins WHICH score, and it is the only assertion that fails if the column is ever re-pointed
//        at getBestQScore() or getQscore2D() -- both are per-PeakGroup too, so they would keep the shape
//        and the mass order intact while silently no longer explaining the ranking.
/////////////////////////////////////////////////////////////
START_SECTION(results_ms1_deconv_qscore_order)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  ABORT_IF(ms1.empty())

  const std::string dir = freshLogDir("lf_r5");
  // FIXTURE, and it must stay this one: buildJsonWithLogDir(dir) takes the DEFAULT enable_ms3 = false,
  // which gives targeting "none" and an EMPTY inclusion list. Passing `true` selects "inclusion" with
  // ../../FlashIDA/test-data/configs/inclusion_cytc.txt, whose priority-1 entry for mass 12351.3 makes
  // target_priority_map_ non-empty; the std::stable_sort at PrecursorSelection.cpp:261-272 then re-orders
  // any pair whose qscores differ by less than tie_threshold (0.1) by target priority, LEGITIMATELY
  // breaking monotonicity. So do not "helpfully" switch this section to the MS3 recipe -- under it the
  // assertion below would be wrong about the engine, not the other way round. The config's
  // rank_by "qscore" matters for the same reason: "intensity" would sort the survey by intensity instead.
  std::string json = buildJsonWithLogDir(dir);
  FLASHIda ida(const_cast<char*>(json.c_str()));

  // MS1-only drive, exactly as §L1: the canonical interleaved driver feeds each survey back under the
  // engine's OWN emitted tracking id (the always-on MS1 gate rejects fabricated ones). ms2 selection is
  // "none" here, so MS2 commands are recorded but never fed -- no MS2 fixture is needed.
  AcqResult acq = runInterleaved(&ida, ms1, std::vector<ScanData>{});
  TEST_TRUE(acq.ms2_cmds.size() > 0)

  auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
  auto res = TSVFile::parse(dir + "/scan_results.tsv");
  auto level = commandLevels(cmds);

  bool saw_multi = false, non_increasing = true;
  for (const auto& row : res.rows)
  {
    auto it = level.find(cell(res, row, "tracking_id"));
    if (it != level.end() && it->second != 1) continue;  // MS1 rows only (classified as in §R1)
    auto qs = splitTokens(cell(res, row, "deconv_qscores"), ';');
    if (qs.size() < 2) continue;                         // a single-mass list cannot be out of order
    saw_multi = true;
    for (size_t i = 0; i + 1 < qs.size(); ++i)
      non_increasing = non_increasing && (toD(qs[i]) >= toD(qs[i + 1]));
  }
  TEST_TRUE(saw_multi)        // ANTI-VACUITY: a survey with <=1 deconvolved mass per scan proves nothing
  TEST_TRUE(non_increasing)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §R3 -- scan_results.tsv : MS2-exploration columns + child_ids (B7)
/////////////////////////////////////////////////////////////
START_SECTION(results_ms2_exploration_columns)
{
  auto ms1 = loadTsvScans(FI_MS1_STD);
  auto ms2 = loadTsvScans(FI_MS2_HCD);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_r3");
  std::string json = injectRuntime(explorationConfig(), dir);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  TEST_TRUE(r.group_commands > 0)

  auto res = TSVFile::parse(dir + "/scan_results.tsv");

  bool found_expl = false;
  for (const auto& row : res.rows)
  {
    int gid = std::atoi(cell(res, row, "exploration_group_id").c_str());
    if (gid < 0) continue;  // exploration rows only (>=1; non-expl pins -1)
    found_expl = true;
    TEST_TRUE(gid >= 1)
    int metric = std::atoi(cell(res, row, "exploration_metric").c_str());
    TEST_TRUE(metric >= 1 && metric <= 4)   // 4 = kinetic_model (ADR-0047)
    int vi = std::atoi(cell(res, row, "variant_index").c_str());
    int tv = std::atoi(cell(res, row, "total_variants").c_str());
    if (vi == -1) continue;  // defensive: since ADR-0048 no sweep acquires a reference row, so none carries -1
    TEST_TRUE(vi >= 0 && tv > 0 && vi < tv)
    TEST_TRUE(posFinite(toD(cell(res, row, "collision_energy"))))   // variant CE (single stage at MS2)
    TEST_TRUE(inActivationSet(cell(res, row, "activation_type")))
    TEST_TRUE(finiteVal(toD(cell(res, row, "exploration_score"))))
    // B7: child_ids count matches commands_pushed and each is a well-formed tracking id.
    // (Children are space-split; they may not be dequeued into the commands TSV yet, so we
    // check id FORMAT here -- the full child->command join is covered by §R1 / join_integrity.)
    auto kids = splitTokens(cell(res, row, "child_ids"), ' ');
    int pushed = std::atoi(cell(res, row, "commands_pushed").c_str());
    TEST_EQUAL((int)kids.size(), pushed)
    bool kids_ok = true;
    for (const auto& k : kids) kids_ok = kids_ok && isTrackingId(k);
    TEST_TRUE(kids_ok)
  }
  TEST_TRUE(found_expl)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F5 -- scan_results.tsv : every exploration variant logs ITS OWN sweep point (variant_index spans
//        0..total_variants-1 exactly once per group) and exactly ONE row per completed group names the
//        winner via the trailing winner_tracking_id column. Catches F5 — the winner-overwrite at
//        Exploration.cpp re-stamped the LAST-resolved variant's row with the WINNER's index/CE/score, so
//        a variant_index was duplicated and another vanished, and the winner was otherwise unidentifiable.
//        Set-membership / cross-row resolution only (no captured CE/score) -> drift-stable.
/////////////////////////////////////////////////////////////
START_SECTION(results_exploration_per_variant_and_winner)
{
  auto ms1 = loadTsvScans(FI_MS1_STD);
  auto ms2 = loadTsvScans(FI_MS2_HCD);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_f5");
  std::string json = injectRuntime(explorationConfig(), dir);            // mass_count -> 5 variants, no baseline
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);                  // feeds all variants -> group completes
  ABORT_IF(r.group_commands == 0)

  auto res = TSVFile::parse(dir + "/scan_results.tsv");

  std::map<int, std::set<int>> seen_vi;     // group_id -> set of variant_index logged
  std::map<int, int> tv_of, winner_count;   // group_id -> total_variants / # rows naming a winner
  std::map<int, std::string> winner_of;     // group_id -> the winner's tracking_id
  std::set<std::string> all_row_ids;        // every exploration row's own tracking_id
  bool found_expl = false;
  for (const auto& row : res.rows)
  {
    int gid = std::atoi(cell(res, row, "exploration_group_id").c_str());
    if (gid < 0) continue;                  // exploration rows only (non-expl pin -1)
    found_expl = true;
    int vi = std::atoi(cell(res, row, "variant_index").c_str());
    int tv = std::atoi(cell(res, row, "total_variants").c_str());
    if (vi == -1) continue;  // CE-0 baseline row: excluded from total_variants; not a sweep point
    TEST_TRUE(vi >= 0 && tv > 0 && vi < tv)
    seen_vi[gid].insert(vi);
    tv_of[gid] = tv;
    all_row_ids.insert(cell(res, row, "tracking_id"));
    std::string w = cell(res, row, "winner_tracking_id");
    if (!w.empty()) { winner_of[gid] = w; winner_count[gid]++; }
  }
  TEST_TRUE(found_expl)

  for (const auto& kv : seen_vi)
  {
    int gid = kv.first;
    const std::set<int>& vis = kv.second;
    // ISSUE(F5): pre-fix variant_index did NOT self-describe (e.g. {0,0,1,2,3} — winner index duplicated,
    // the completing variant's own index missing) because the completing row was overwritten with the winner.
    for (int k = 0; k < tv_of[gid]; ++k) TEST_EQUAL(vis.count(k), 1)
    TEST_EQUAL((int)vis.size(), tv_of[gid])
    // ISSUE(F5): pre-fix there was NO winner column, so once per-variant overwriting stopped the winner
    // would be unidentifiable. Exactly one row per completed group names a real, in-set winner.
    TEST_EQUAL(winner_count[gid], 1)
    TEST_TRUE(isTrackingId(winner_of[gid]) && all_row_ids.count(winner_of[gid]) == 1)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// §I4a -- scan_results.tsv : exploration FOLLOW-UP. A non-tolerance override makes the MS2-exploration
//         winner re-acquire as a production MS2 (then MS2->MS3) on inclusion-pinned cytC. Mirrors C++
//         ms2_exploration_production_winner_then_ms3 + checks the LOGGED exploration columns (matching data).
/////////////////////////////////////////////////////////////
START_SECTION(results_exploration_followup_columns)
{
  auto ms1 = loadTsvScans(FI_MS1_CYTC2);
  auto ms2 = loadTsvScans(FI_MS2_CYTC2);
  ABORT_IF(ms1.empty() || ms2.empty())

  // Inclusion-pin the cytC precursor + M-start proteoform; add a non-tolerance override ("analyzer")
  // so config_.level(2).overrides stays non-empty -> the winner is a production-MS2 re-acquisition.
  std::string cfg = inclusionPinCytc2(ms3SelectionOnlyConfig());
  {
    auto p = cfg.find("\"metric\": \"mass_count\"");
    ABORT_IF(p == std::string::npos)
    cfg.insert(p, "\"overrides\": { \"analyzer\": \"Orbitrap\" }, ");
  }
  const std::string dir = freshLogDir("lf_i4followup");
  cfg = injectRuntime(cfg, dir);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  ABORT_IF(r.group_commands == 0)

  // Override non-empty -> winner re-acquires as production MS2; feeding it drives the MS2->MS3 path.
  TEST_TRUE(r.found_production_ms2)
  TEST_TRUE(r.found_ms3)

  // Logged exploration result rows (inclusion-pinned cytC => they MATCH, unlike §R3's E. coli data).
  auto res = TSVFile::parse(dir + "/scan_results.tsv");
  bool found_expl = false;
  for (const auto& row : res.rows)
  {
    int gid = std::atoi(cell(res, row, "exploration_group_id").c_str());
    if (gid < 0) continue;
    found_expl = true;
    TEST_TRUE(gid >= 1)
    int vi = std::atoi(cell(res, row, "variant_index").c_str());
    int tv = std::atoi(cell(res, row, "total_variants").c_str());
    if (vi == -1) continue;  // CE-0 baseline row (variant_index -1): a reference, not a sweep point
    TEST_TRUE(vi >= 0 && tv > 0 && vi < tv)
    TEST_TRUE(finiteVal(toD(cell(res, row, "exploration_score"))))
    auto kids = splitTokens(cell(res, row, "child_ids"), ' ');
    int pushed = std::atoi(cell(res, row, "commands_pushed").c_str());
    TEST_EQUAL((int)kids.size(), pushed)
  }
  TEST_TRUE(found_expl)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §I4b -- scan_results.tsv : MS3-level EXPLORATION. A two-level MS2->MS3 CE-sweep cascade on inclusion-
//         pinned cytC. Mirrors C++ ms2_then_ms3_exploration_acquires_ms3 + checks the LOGGED exploration
//         columns. (MS3 result rows themselves are golden-locked in the C# Golden_Exploration_MS3_CytC.)
/////////////////////////////////////////////////////////////
START_SECTION(results_ms3_exploration_columns)
{
  auto ms1 = loadTsvScans(FI_MS1_CYTC2);
  auto ms2 = loadTsvScans(FI_MS2_CYTC2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_i4ms3");
  std::string cfg = injectRuntime(inclusionPinCytc2(ms3ExplorationConfig()), dir);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  ABORT_IF(r.group_commands == 0)

  // The MS2-exploration winner cascades to an MS3 (two-stage) command.
  TEST_TRUE(r.found_ms3)
  TEST_EQUAL(r.ms3_num_stages, 2)

  auto res = TSVFile::parse(dir + "/scan_results.tsv");
  bool found_expl = false;
  bool found_ms3_expl = false;  // §I4b(Issue-1): >=1 MS3-exploration result row whose CE/act is 2-stage
  for (const auto& row : res.rows)
  {
    int gid = std::atoi(cell(res, row, "exploration_group_id").c_str());
    if (gid < 0) continue;
    found_expl = true;
    TEST_TRUE(gid >= 1)
    int metric = std::atoi(cell(res, row, "exploration_metric").c_str());
    TEST_TRUE(metric >= 1 && metric <= 4)   // 4 = kinetic_model (ADR-0047)
    int vi = std::atoi(cell(res, row, "variant_index").c_str());
    int tv = std::atoi(cell(res, row, "total_variants").c_str());
    if (vi == -1) continue;  // defensive: since ADR-0048 no sweep acquires a reference row, so none carries -1
    TEST_TRUE(vi >= 0 && tv > 0 && vi < tv)

    // §I4b ISSUE(Issue-1): the MS3-exploration *result* rows used to log SINGLE-stage CE/activation
    // (only the MS3 fragmentation half), even though the matching scan_commands row logged the correct
    // 2-stage "ms2CE;ms3CE" / "HCD;CID". The Issue-1 engine fix carries the MS2 isolation stage through
    // FeedResultInfo, so the MS3 result row is now 2-stage too. Gate on the scan_results ms_level column
    // (col 1; this section does not join to scan_commands). ms3ExplorationConfig() SWEEPS BOTH
    // stages -- MS2 CE in {20,25,30,35,40} (ce_min 20, ce_max 40, step 5) and MS3 CE in {15,20,25,30,35}
    // (ce_min 15, ce_max 35, step 5) -- so DO NOT pin single values: assert exactly 2 ';'-tokens with
    // ce[0] (MS2) in [20,40] and ce[1] (MS3) in [15,35], activation == "HCD;CID".
    // PRE-FIX: collision_energy/activation_type split to size()==1 -> these TEST_EQUALs FAIL.
    if (cell(res, row, "ms_level") == "3")
    {
      auto ce  = splitTokens(cell(res, row, "collision_energy"), ';');
      auto act = splitTokens(cell(res, row, "activation_type"), ';');
      TEST_EQUAL(ce.size(), 2)
      TEST_EQUAL(act.size(), 2)
      TEST_TRUE(toD(ce[0]) >= 20.0 - 1.0 && toD(ce[0]) <= 40.0 + 1.0)  // MS2 isolation CE in swept range
      TEST_TRUE(toD(ce[1]) >= 15.0 - 1.0 && toD(ce[1]) <= 35.0 + 1.0)  // MS3 fragmentation CE in swept range
      TEST_EQUAL(act[0], std::string("HCD"))   // MS2 isolation activation
      TEST_EQUAL(act[1], std::string("CID"))   // MS3 fragmentation activation (ms3ExplorationConfig)
      found_ms3_expl = true;
    }
  }
  TEST_TRUE(found_expl)
  TEST_TRUE(found_ms3_expl)   // >=1 MS3-exploration result row was checked (else the 2-stage pin proves nothing)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F1 -- scan_commands.tsv : every MS3 exploration command's parent resolves to an MS2 (not the MS1
//        survey). Catches F1 — the feedResult push-loop used to overwrite buildMS3's correct MS2 parent
//        with the group's MS1 originating id. Cross-stream/level relationship only -> drift-stable.
/////////////////////////////////////////////////////////////
START_SECTION(commands_exploration_ms3_parent_is_ms2)
{
  auto ms1 = loadTsvScans(FI_MS1_CYTC2);
  auto ms2 = loadTsvScans(FI_MS2_CYTC2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_f1");
  std::string cfg = injectRuntime(inclusionPinCytc2(ms3ExplorationConfig()), dir);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  ABORT_IF(r.group_commands == 0)
  ABORT_IF(!r.found_ms3)   // need >=1 MS3 command emitted to check its parent

  auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
  auto level = commandLevels(cmds);   // tracking_id -> ms_level (from scan_commands)
  int checked = 0;
  for (const auto& row : cmds.rows)
  {
    std::string tid = cell(cmds, row, "tracking_id");
    auto it = level.find(tid);
    if (it == level.end() || it->second != 3) continue;   // MS3 commands only
    std::string p = cell(cmds, row, "parent_tracking_id");
    // ISSUE(F1): pre-fix p resolved to an MS1 (the group's grandparent survey) because the
    // exploration push-loop overwrote buildMS3's correct MS2 parent with info.parent_scan_id.
    auto pit = level.find(p);
    TEST_TRUE(pit != level.end() && pit->second == 2)
    checked++;
  }
  TEST_TRUE(checked >= 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F2 -- scan_commands.tsv : exploration MS3 commands carry REAL stage-1 (fragment) scores. Catches F2 —
//        the exploration buildMS3 used to get a default-empty FragmentScores, so the ';'-joined stage-1
//        token of the scoring columns was 0 (and identification ms3_window_snr/ms3_charge_intensity = 0).
//        The single-group driver can't feed MS3 spectra (no MS3 identification rows), but it DOES drain+log
//        the MS3 commands, whose stage-1 scalars are the upstream value the bug zeroed. Range only -> drift-stable.
/////////////////////////////////////////////////////////////
START_SECTION(commands_exploration_ms3_stage1_scores_nonzero)
{
  auto ms1 = loadTsvScans(FI_MS1_CYTC2);
  auto ms2 = loadTsvScans(FI_MS2_CYTC2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_f2");
  std::string cfg = injectRuntime(inclusionPinCytc2(ms3ExplorationConfig()), dir);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  ABORT_IF(r.group_commands == 0)
  ABORT_IF(!r.found_ms3)   // need >=1 MS3 command to inspect its stage-1 scores

  auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
  int checked = 0;
  for (const auto& row : cmds.rows)
  {
    if (std::atoi(cell(cmds, row, "ms_level").c_str()) != 3) continue;   // MS3 commands only
    // MS3 scoring columns are 2-stage ';'-joined "MS2;fragment"; token[1] is the stage-1 (fragment) value.
    auto pin = splitTokens(cell(cmds, row, "precursor_intensity"), ';');
    auto csn = splitTokens(cell(cmds, row, "charge_snr"), ';');
    TEST_EQUAL(pin.size(), 2)
    TEST_EQUAL(csn.size(), 2)
    // ISSUE(F2): pre-fix the stage-1 tokens were 0 because exploration buildMS3 got the default-empty
    // FragmentScores; now they carry the selected fragment's real deconvolution scores.
    TEST_TRUE(toD(pin[1]) > 0.0)
    TEST_TRUE(finiteVal(toD(csn[1])) && toD(csn[1]) >= 0.0)
    checked++;
  }
  TEST_TRUE(checked >= 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F3 -- MS3 exploration ('E') identification rows produce a proteoform, and the wide MS3 fragment lives on
//        scan_commands.ms3_proteoform (the identification 'E' proteoform is NARROWED ⊆ the scan_commands wide
//        fragment; T-N2 exploration sink). Needs the MS3 variants FED (the batch must complete), so this drives
//        the inclusion-pinned MS3-exploration config with the real per-ion MS3 manifest (skips cleanly when
//        no ms3_cytc fixtures exist; the exact values are golden-locked in C# Golden_Exploration_MS3).
/////////////////////////////////////////////////////////////
START_SECTION(identification_ms3_exploration_proteoform)
{
  auto manifest = buildMs3IonManifest("../../FlashIDA/test-data/spectra");
  if (manifest.empty())
  {
    TEST_TRUE(true)   // no real MS3 fragment fixtures -> skip cleanly (never fabricate)
  }
  else
  {
    auto ms1 = loadTsvScans(CYTC_MS1);
    auto ms2 = loadTsvScans(CYTC_MS2);
    ABORT_IF(ms1.empty() || ms2.empty())

    const std::string dir = freshLogDir("lf_f3");
    std::string cfg = injectRuntime(inclusionPinCytc(ms3ExplorationConfig()), dir);
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    runFullCycle(&ida, ms1, ms2, &manifest);   // feed the MS3 variants so MS3-'E' rows are produced

    auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
    auto idf = TSVFile::parse(dir + "/identification.tsv");

    // F3-id: a matched MS3-'E' identification row has a non-empty ProForma proteoform.
    int matched = 0;
    for (const auto& row : idf.rows)
    {
      if (std::atoi(cell(idf, row, "ms_level").c_str()) != 3) continue;
      if (cell(idf, row, "scan_mode") != "E") continue;
      if (cell(idf, row, "ms3_fragments").empty()) continue;   // only rows that matched fragments
      // ISSUE(F3-id): pre-fix proteoform was blank (use_ctx_proteoform was 'R'-only).
      TEST_TRUE(isProForma(cell(idf, row, "proteoform")))
      // C2: MS3 per-ion coverage is a fraction in [0,1] on a matched MS3 row.
      { double cov = toD(cell(idf, row, "ms3_fragment_coverage")); TEST_TRUE(cov >= 0.0 && cov <= 1.0) }
      matched++;
    }
    // F3-cmd: the wide MS3 fragment proteoform now lives on scan_commands.ms3_proteoform (moved off
    // scan_results). Every MS3 command row that carries it must be a CLIPPED FRAGMENT (bare-residue count
    // < full cytC 105); the tracking_id -> proteoform map feeds the T-N2 subset check below.
    std::map<std::string, std::string> cmd_ms3pf;
    for (const auto& row : cmds.rows)
    {
      if (std::atoi(cell(cmds, row, "ms_level").c_str()) != 3) continue;
      const std::string pf = cell(cmds, row, "ms3_proteoform");
      if (pf.empty()) continue;   // MS3 command without a fragment context -> tolerated
      cmd_ms3pf[cell(cmds, row, "tracking_id")] = pf;
      TEST_TRUE(isProForma(pf))
      int bare = 0; for (char ch : pf) if (ch >= 'A' && ch <= 'Z') ++bare;
      TEST_TRUE(bare > 0 && bare < 105)   // ISSUE(C1): the wide clipped fragment, not the 105-mer parent
    }
    TEST_TRUE(matched >= 1)   // the inclusion-pinned MS3-exploration cascade must yield >=1 matched MS3-'E' row

    // T-N2 (exploration sink): under Change L the identification 'E' MS3 proteoform is a per-scan EQUIV-frame
    // bracket that MAY exceed the wide scan_commands base AND may MERGE co-observed mods (fewer ranges). The old
    // id-subset-of-cmd (and the exact size-equality) are retired; what holds is: the leaf never SPLITS a cmd mod
    // (leaf mod count <= cmd mod count -- merges only reduce) and every leaf range is well-formed.
    int e_cmp_rows = 0;
    for (const auto& row : idf.rows)
    {
      if (std::atoi(cell(idf, row, "ms_level").c_str()) != 3) continue;
      if (cell(idf, row, "scan_mode") != "E") continue;
      if (cell(idf, row, "ms3_fragments").empty()) continue;
      auto cit = cmd_ms3pf.find(cell(idf, row, "tracking_id"));
      if (cit == cmd_ms3pf.end()) continue;                 // need the paired scan_commands proteoform
      auto id_ranges  = parseModRanges(cell(idf, row, "proteoform"));
      auto cmd_ranges = parseModRanges(cit->second);        // wide fragment from scan_commands
      TEST_TRUE(id_ranges.size() <= cmd_ranges.size())      // Change L merges co-observed mods; never splits
      for (size_t k = 0; k < id_ranges.size(); ++k)
        TEST_TRUE(id_ranges[k].first >= 1 && id_ranges[k].first <= id_ranges[k].second)   // well-formed range
      ++e_cmp_rows;
    }
    TEST_TRUE(e_cmp_rows >= 1)   // >=1 MS3-'E' row joined id<->scan_commands -> exploration sink covered
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F4 -- MS2 exploration ('E') identification rows source ms1_precursor_* from the GROUP's own deconvolved
//        precursor, not the originating MS1 survey command. Catches F4 — buildMS2ContextForVariant filled
//        ms1_precursor_* from group.originating_cmd ONLY when num_stages>0, but for an MS2 exploration group
//        the originating_cmd is the MS1 survey (num_stages==0), so the guard failed and the fields stayed 0.
//        Drift-stable cross-stream ppm join: an MS2-'E' id row's ms1_precursor_mass must equal the SAME
//        tracking_id's scan_commands.mono_mass (both = the group's deconvolved precursor) within ppm.
/////////////////////////////////////////////////////////////
START_SECTION(identification_ms2_exploration_ms1_precursor_from_group)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_f4");
  // single-group drive feeds the MS2 variants back -> feedResult writes the MS2-'E' identification rows.
  std::string cfg = injectRuntime(inclusionPinCytc(ms3ExplorationConfig()), dir);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  auto r = driveOneExplorationGroup(&ida, ms1, ms2[0]);
  ABORT_IF(r.group_commands == 0)

  auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
  auto idf  = TSVFile::parse(dir + "/identification.tsv");

  // tracking_id -> mono_mass from scan_commands (an MS2 command's mono_mass is a single value == the
  // precursor mono mass; only MS3 commands ';'-join two stages, which these MS2-'E' ids never resolve to).
  std::map<std::string, double> cmd_mono;
  for (const auto& row : cmds.rows)
    cmd_mono[cell(cmds, row, "tracking_id")] = toD(cell(cmds, row, "mono_mass"));

  int checked = 0;
  for (const auto& row : idf.rows)
  {
    if (std::atoi(cell(idf, row, "ms_level").c_str()) != 2) continue;   // MS2 identification rows only
    if (cell(idf, row, "scan_mode") != "E") continue;                   // exploration variants only
    auto it = cmd_mono.find(cell(idf, row, "tracking_id"));
    if (it == cmd_mono.end() || it->second <= 0.0) continue;            // need the matching command row
    double id_mass = toD(cell(idf, row, "ms1_precursor_mass"));
    // ISSUE(F4): pre-fix ms1_precursor_mass came from group.originating_cmd, which for an MS2 exploration
    // group is the MS1 survey (num_stages==0) -> the >0 guard failed -> the field stayed 0.0.
    TEST_TRUE(id_mass > 0.0)
    double ppm = std::abs(id_mass - it->second) / it->second * 1e6;
    TEST_TRUE(ppm < 50.0)
    checked++;
  }
  TEST_TRUE(checked >= 1)   // >=1 MS2-'E' id row joined to its command within ppm
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F8-incl -- inclusion (cytC) selects the LISTED precursor within the target window. Catches the
//   F8-inclusion bug: the target was the cytC AVERAGE mass (~12358.31), ~570-700 ppm off the engine's
//   MONOISOTOPIC deconvolution (~12351.3), so the ~0.247 Da match window never matched -> silent DDA
//   fall-through (DDA picks cytC anyway, masking the miss). inclusion_cytc.txt is now the monoisotopic
//   12351.3. NON-STRICT: >=1 MS2 within the window (a floor — DDA would also pick cytC). STRICT: selection
//   is RESTRICTED to the list, so EVERY MS2 is within the window AND >=1 exists (the real lever — on the
//   wrong target strict matches nothing -> 0 MS2). Window membership only -> drift-stable.
/////////////////////////////////////////////////////////////
START_SECTION(inclusion_ms2_within_target_window)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const double T = 12351.3, HW = 2.0 * 10.0 * T * 1e-6;   // 10 ppm tol, doubled at the match site (~0.247 Da)
  auto onList = [&](double m){ return std::abs(m - T) < HW; };

  // Two engines are built below (non-strict, then strict), so each gets its OWN log folder.
  auto driveCountMs2 = [&](const std::string& json, const std::string& dir, int& ms2_total, int& ms2_on_list) {
    FLASHIda ida(const_cast<char*>(json.c_str()));
    runFullCycle(&ida, ms1, ms2);
    auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
    ms2_total = 0; ms2_on_list = 0;
    for (const auto& row : cmds.rows)
    {
      if (std::atoi(cell(cmds, row, "ms_level").c_str()) != 2) continue;   // MS2 commands only
      ms2_total++;
      if (onList(toD(cell(cmds, row, "mono_mass")))) ms2_on_list++;        // mono_mass is single-valued at MS2
    }
  };

  // NON-STRICT: the cytC precursor IS selected (>=1 MS2 within the window); DDA may add off-list MS2 too.
  int ns_total = 0, ns_on = 0;
  const std::string ns_dir = freshLogDir("lf_f8incl_1");
  driveCountMs2(buildJsonWithLogDir(ns_dir, true), ns_dir, ns_total, ns_on);
  // ISSUE(F8-incl): pre-fix ns_on==0 — the target was the cytC AVERAGE mass, ~570-700 ppm off the
  // monoisotopic ~12351.3, so the window never matched (silent DDA fall-through).
  TEST_TRUE(ns_on >= 1)

  // STRICT: selection RESTRICTED to the list -> EVERY MS2 is within the window, and >=1 exists.
  const std::string s_dir = freshLogDir("lf_f8incl_2");
  std::string strict_json = buildJsonWithLogDir(s_dir, true);
  {
    const std::string from = "\"strict_inclusion\": false";
    auto p = strict_json.find(from);
    ABORT_IF(p == std::string::npos)
    strict_json.replace(p, from.size(), "\"strict_inclusion\": true");
  }
  int s_total = 0, s_on = 0;
  driveCountMs2(strict_json, s_dir, s_total, s_on);
  // ISSUE(F8-incl-strict): pre-fix s_total==0 (nothing matched the wrong target); any leaked MS2 would be an
  // off-list DDA pick -> s_on != s_total. With the corrected monoisotopic target every strict MS2 is on-list.
  TEST_TRUE(s_total >= 1 && s_on == s_total)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F8b -- target-log reader : a >= 10 target count must not be read as "0 targets"
/////////////////////////////////////////////////////////////
START_SECTION(target_log_double_digit_target_count_is_not_skipped)
{
  // PrecursorSelection's target_log_files loader is the SECOND reader of the ida.log grammar, and it
  // is the in-scope one: feeding an engine-written ida.log back in as a target log is the round-trip
  // this covers. It carried the same unanchored `line.find("0 targets")` as parseFLASHIdaLog, so a
  // "- 10 targets" or "- 20 targets" header was dropped -- and here the consequence is not a mis-keyed
  // map but a WRONG RETENTION TIME: `rt` keeps the previous header's value, or 0.0 at the head of the
  // file, and mode 2's suppression is gated on |survey_rt - prt| < rt_window
  // (PrecursorSelection.cpp:332).
  //
  // Its Mass= path had no test at all before this: the only existing target-log test
  // (exclusion_mode2_tqscore_suppresses_target_mass, above) uses the older space-delimited
  // "Mass <v>" grammar and headers with no "targets" text, so it cannot see this defect.
  const std::string FI_MS1_ECOLI = "../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt";
  auto ms1 = loadTsvScans(FI_MS1_ECOLI);
  auto ms2 = loadTsvScans(FI_MS2_HCD);
  ABORT_IF(ms1.empty() || ms2.empty())

  // rt_window 60, NOT the 180 the section above uses, and that is the whole reason this test can
  // fail. The ecoli surveys sit at rt 151-160; the log header below is written at RT 2.6 -> prt 156.
  //   fixed : |151..160 - 156| <= 5  <  60  -> in window  -> suppression applies
  //   broken: header skipped -> prt 0 -> |151 - 0| = 151 >= 60 -> out of window -> NO suppression
  // At 180 both cases are in window and the assertion would hold under the bug -- i.e. vacuous.
  // 60 still spans the whole 9-unit survey range, so mass exclusion is unaffected.
  auto cfg = [](const char* targeting, int max_precursors, const std::string& target_log_json) {
    std::ostringstream o;
    o << R"({
      "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10] },
      "precursor_selection": { "rt_window": 60, "targeting": ")" << targeting << R"(", "consider_all_charges": false, "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": )" << max_precursors << R"( },
      "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
      "faims": { "cv_values": [], "max_cv_skip": 0 },
      "ms_settings": { "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 }, "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 } },
      "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 }, "agc_interval_seconds": 999999 },
      "files": { "target_logs": [)" << target_log_json << R"(], "inclusion_list": "" }
    })";
    return o.str();
  };

  auto selectedMs2Masses = [&](const std::string& json, const std::string& tag) {
    const std::string dir = freshLogDir(tag);
    std::string full = injectRuntime(json, dir);
    FLASHIda ida(const_cast<char*>(full.c_str()));
    runFullCycle(&ida, ms1, ms2);
    auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
    std::vector<double> masses;
    for (const auto& row : cmds.rows)
      if (std::atoi(cell(cmds, row, "ms_level").c_str()) == 2)
        masses.push_back(toD(cell(cmds, row, "mono_mass")));
    return masses;
  };
  auto hasNear = [](const std::vector<double>& v, double t) {
    for (double m : v) if (std::abs(m - t) < 0.5) return true;   // suppression is per-NOMINAL mass
    return false;
  };

  // 1) DDA baseline names the mass to suppress. >= 2 entries = the survey is productive.
  auto dda = selectedMs2Masses(cfg("none", 1, ""), "lf_f8b_dda");
  ABORT_IF(dda.size() < 2)
  const double target = dda[0];
  TEST_TRUE(hasNear(dda, target))

  // 2) The same log content under two headers that differ ONLY in the target count. Two observations
  //    at qscore 0.7 give tqscore = 1 - (1-0.7)^2 = 0.91 > the 0.9 threshold, so the mass is
  //    de-prioritized out of the single MS2 slot -- provided the header's RT was read.
  auto writeLog = [&](const std::string& path, const char* count_token) {
    std::ofstream lg(path);
    for (int i = 0; i < 2; i++)
    {
      lg << "MS1 Scan# " << (7 + i) << " RT 2.6000 (Access ID !!\") - " << count_token << " targets\n";
      lg << "Mass=" << std::to_string(target) << "\tZ=4\tScore=0.70000"
         << "\tWindow=[500.0000-501.0000]\tPrecursorIntensity=1.00000"
         << "\tPrecursorMassIntensity=2.00000\tFeatures=[0.9,1.0,0.9,1.0,0.9,1.0]"
         << "\tChargeRange=[4-6]\tHCD=0\n";
    }
  };

  // CONTROL: a single-digit count was never affected by the substring bug, so this suppresses both
  // before and after the fix. It proves the suppression mechanism is live here at rt_window 60 --
  // without it, the assertion below could pass merely because nothing is ever suppressed.
  const std::string log_single = "lf_f8b_single.log";
  writeLog(log_single, "4");
  auto suppressed_single = selectedMs2Masses(cfg("in_depth", 1, "\"" + log_single + "\""), "lf_f8b_single");
  TEST_TRUE(!hasNear(suppressed_single, target))

  // THE REGRESSING CASE: identical but for the count. Before the fix the header is swallowed by
  // `find("0 targets")`, rt stays 0.0, the observations are filed 151 s away from every survey, the
  // suppression never applies and `target` keeps the slot.
  const std::string log_double = "lf_f8b_double.log";
  writeLog(log_double, "10");
  auto suppressed_double = selectedMs2Masses(cfg("in_depth", 1, "\"" + log_double + "\""), "lf_f8b_double");
  TEST_TRUE(!hasNear(suppressed_double, target))

  std::remove(log_single.c_str());
  std::remove(log_double.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// §F8c -- target-log reader : a superseded decision is not an observation (survey re-planning)
/////////////////////////////////////////////////////////////
START_SECTION(target_log_skips_superseded_masses)
{
  // A survey's re-plan can drop a waiting MS2 an older survey decided. That decision stays in the older
  // entry's Mass= line (ida.log is a decision-time record), and the newer entry names it on a Superseded=
  // line. The target_logs loader must not count it: in_depth's 1 - prod(1 - qscore) would otherwise read a
  // command that was never acquired as an extra observation, and push the mass over tqscore_threshold.
  //
  // Same drive as the section above: rich ecoli surveys, rt_window 60, in_depth, ONE MS2 slot. Two
  // observations at 0.7 give 0.91 > 0.9 and suppress `target`; one gives 0.7 and does not.
  const std::string FI_MS1_ECOLI = "../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt";
  auto ms1 = loadTsvScans(FI_MS1_ECOLI);
  auto ms2 = loadTsvScans(FI_MS2_HCD);
  ABORT_IF(ms1.empty() || ms2.empty())

  auto cfg = [](const char* targeting, const std::string& target_log_json) {
    std::ostringstream o;
    o << R"({
      "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10] },
      "precursor_selection": { "rt_window": 60, "targeting": ")" << targeting << R"(", "consider_all_charges": false, "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 1 },
      "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
      "faims": { "cv_values": [], "max_cv_skip": 0 },
      "ms_settings": { "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 }, "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 } },
      "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 }, "agc_interval_seconds": 999999 },
      "files": { "target_logs": [)" << target_log_json << R"(], "inclusion_list": "" }
    })";
    return o.str();
  };
  auto selectedMs2Masses = [&](const std::string& json, const std::string& tag) {
    const std::string dir = freshLogDir(tag);
    std::string full = injectRuntime(json, dir);
    FLASHIda ida(const_cast<char*>(full.c_str()));
    runFullCycle(&ida, ms1, ms2);
    auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
    std::vector<double> masses;
    for (const auto& row : cmds.rows)
      if (std::atoi(cell(cmds, row, "ms_level").c_str()) == 2)
        masses.push_back(toD(cell(cmds, row, "mono_mass")));
    return masses;
  };
  auto hasNear = [](const std::vector<double>& v, double t) {
    for (double m : v) if (std::abs(m - t) < 0.5) return true;   // suppression is per-NOMINAL mass
    return false;
  };

  auto dda = selectedMs2Masses(cfg("none", ""), "lf_f8c_dda");
  ABORT_IF(dda.size() < 2)
  const double target = dda[0];

  // Two entries, each deciding `target` under its own command id (!!# then !!%). With @p supersede the
  // second entry names the first decision on its Superseded= line, right after its header -- where the
  // writer puts it.
  auto writeLog = [&](const std::string& path, bool supersede) {
    std::ofstream lg(path);
    const char* ids[2] = {"!!#", "!!%"};
    for (int i = 0; i < 2; i++)
    {
      lg << "MS1 Scan# " << (7 + i) << " RT 2.6000 (Access ID !!" << (i == 0 ? "\"" : "$") << ") - 1 targets\n";
      if (supersede && i == 1) lg << "Superseded=!!#\n";
      lg << "Mass=" << std::to_string(target) << "\tZ=4\tScore=0.70000"
         << "\tWindow=[500.0000-501.0000]\tPrecursorIntensity=1.00000"
         << "\tPrecursorMassIntensity=2.00000\tFeatures=[0.9,1.0,0.9,1.0,0.9,1.0]"
         << "\tChargeRange=[4-6]\tHCD=0\tID=" << ids[i] << "\n";
      lg << "AllMass=" << std::to_string(target) << "\n";
    }
  };

  // CONTROL: both decisions stand -> two observations -> suppressed. Proves the ID= suffix leaves the
  // Mass= line readable and that suppression is live here; without it the assertion below is vacuous.
  const std::string log_kept = "lf_f8c_kept.log";
  writeLog(log_kept, false);
  TEST_TRUE(!hasNear(selectedMs2Masses(cfg("in_depth", "\"" + log_kept + "\""), "lf_f8c_kept"), target))

  // The first decision superseded -> one observation at 0.7 -> below threshold -> `target` keeps the slot.
  const std::string log_superseded = "lf_f8c_superseded.log";
  writeLog(log_superseded, true);
  TEST_TRUE(hasNear(selectedMs2Masses(cfg("in_depth", "\"" + log_superseded + "\""), "lf_f8c_superseded"), target))

  std::remove(log_kept.c_str());
  std::remove(log_superseded.c_str());
}
END_SECTION

/////////////////////////////////////////////////////////////
// [REMOVED §F9] results_ms3_fragment_and_tag_count_are_sentinel — asserted scan_results MS3
//   fragment_count==-1 && tag_count==-1; both columns were removed in the scan_results slim-down
//   (fragment_count derivable from identification.tsv ms3_fragments; tag_count dropped entirely).
//   MS3 fragment-match coverage is still asserted in §T9 (ms3_id_with_frags>=1).
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
// §X1 -- scan_results.tsv : duration semantics (ordering + decomposition)
/////////////////////////////////////////////////////////////
START_SECTION(results_duration_semantics)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_x1");
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto res = TSVFile::parse(dir + "/scan_results.tsv");
  TEST_TRUE(res.rows.size() > 0)
  for (const auto& row : res.rows)
  {
    double dur = toD(cell(res, row, "duration_ms"));
    TEST_TRUE(inDuration(dur))
    double q = toD(cell(res, row, "queue_duration_ms"));
    double inst = toD(cell(res, row, "instrument_duration_ms"));
    double proc = toD(cell(res, row, "processing_duration_ms"));
    TEST_TRUE(nonNegFinite(q) && nonNegFinite(inst) && nonNegFinite(proc))
    double resolve = toD(cell(res, row, "resolve_ts"));
    double recv = toD(cell(res, row, "received_ts"));
    double deq = toD(cell(res, row, "dequeue_ts"));
    // ordering when timestamps are present (>0)
    if (deq > 0 && recv > 0) TEST_TRUE(recv >= deq)
    if (recv > 0 && resolve > 0) TEST_TRUE(resolve >= recv)
    // decomposition when fully tracked
    if (q > 0 && inst > 0 && proc > 0) TEST_TRUE(std::abs(dur - (q + inst + proc)) < 2.0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// §X2 -- scan_results.tsv : parent lineage (MS2->MS1, MS3->MS2)
/////////////////////////////////////////////////////////////
START_SECTION(results_parent_lineage)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_x2");
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
  auto res = TSVFile::parse(dir + "/scan_results.tsv");
  auto level = commandLevels(cmds);

  bool checked_ms2 = false, checked_ms3 = false;
  for (const auto& row : res.rows)
  {
    std::string tid = cell(res, row, "tracking_id");
    std::string parent = cell(res, row, "parent_tracking_id");
    auto it = level.find(tid);
    if (it == level.end() || it->second == 1) { TEST_TRUE(parent.empty()); continue; }  // MS1 survey/input
                                                          // (now in scan_commands.tsv with ms_level=1)
    // MS2/MS3 result rows must carry a well-formed parent that is itself a known id
    TEST_TRUE(isTrackingId(parent))
    if (it->second == 2) { checked_ms2 = true; }
    if (it->second == 3) { checked_ms3 = true; TEST_TRUE(level.count(parent) && level[parent] == 2); }
  }
  TEST_TRUE(checked_ms2)
  TEST_TRUE(checked_ms3)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §X4 -- results vs identification agree for the same tracking_id (#2-4)
/////////////////////////////////////////////////////////////
START_SECTION(results_identification_agreement)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_x4");
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto idf = TSVFile::parse(dir + "/identification.tsv");
  ABORT_IF(idf.rows.size() == 0)   // identification must fire, else this gate proves nothing

  // scan_results no longer carries a proteoform (moved to scan_commands.ms3_proteoform for MS3, and the
  // MS2 proteoform lives only on identification.tsv). Cross-stream scan_commands<->identification agreement
  // is locked in §T9/§F3; here we assert every non-empty identification proteoform is a well-formed ProForma.
  int checked = 0;
  for (const auto& row : idf.rows)
  {
    const std::string p = cell(idf, row, "proteoform");
    if (p.empty()) continue;
    TEST_TRUE(isProForma(p))
    ++checked;
  }
  TEST_TRUE(checked > 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §I1 -- identification.tsv : universal gates on all (R-mode) rows
/////////////////////////////////////////////////////////////
START_SECTION(identification_R_rows_universal)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_i1");
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  auto idf = TSVFile::parse(dir + "/identification.tsv");
  TEST_TRUE(idf.rows.size() > 0)  // cytC (non-exploration) => all rows scan_mode 'R'
  for (const auto& row : idf.rows)
  {
    int lvl = std::atoi(cell(idf, row, "ms_level").c_str());
    TEST_TRUE(lvl == 2 || lvl == 3)
    std::string mode = cell(idf, row, "scan_mode");
    TEST_TRUE(mode == "R" || mode == "E")
    TEST_TRUE(isTrackingId(cell(idf, row, "tracking_id")))
    TEST_TRUE(isProForma(cell(idf, row, "proteoform")))
    TEST_TRUE(std::atoi(cell(idf, row, "start_pos").c_str()) >= 0)
    TEST_TRUE(std::atoi(cell(idf, row, "end_pos").c_str()) >= std::atoi(cell(idf, row, "start_pos").c_str()))
    TEST_TRUE(finiteVal(toD(cell(idf, row, "ppm_offset"))))
    double corr = toD(cell(idf, row, "correction_factor"));
    TEST_TRUE(corr > 0.0 && std::abs(corr - 1.0) < 0.05)
    TEST_TRUE(inMass(toD(cell(idf, row, "ms1_precursor_mass"))))
    TEST_TRUE(posFinite(toD(cell(idf, row, "ms1_precursor_mz"))))
    TEST_TRUE(inChargeD(toD(cell(idf, row, "ms1_precursor_charge"))))
    // ms2_fragments / ms2_fragment_masses token-aligned
    auto f = splitTokens(cell(idf, row, "ms2_fragments"), ';');
    auto m = splitTokens(cell(idf, row, "ms2_fragment_masses"), ';');
    TEST_EQUAL(f.size(), m.size())

    // I2 (isolation-window cols): every matched row reports the MS2 precursor's commanded window width,
    // the SNR over that actual window, and the selected-charge intensity. Plausibility ranges only;
    // the exact values are golden-locked in the C# suite. window_snr is bounded (~1000) for a pure window.
    double w2 = toD(cell(idf, row, "ms2_isolation_width"));
    TEST_TRUE(w2 > 0.0 && w2 < 50.0)
    double s2 = toD(cell(idf, row, "ms2_window_snr"));
    TEST_TRUE(finiteVal(s2) && s2 >= 0.0 && s2 <= 1001.0)
    TEST_TRUE(posFinite(toD(cell(idf, row, "ms2_charge_intensity"))))

    // §I2 folded in: MS3-R rows additionally carry the MS2-fragment precursor (#1: proteoform/
    // region sourced from ctx) plus a non-empty, token-aligned MS3 fragment list. Asserted WHEN
    // an MS3-R row is present -- those require the MS3 fragment match to fire (which depends on
    // real MS3 fragment data, golden-locked exactly in the C# suite, not this cycle's shortcut).
    if (lvl == 3 && mode == "R")
    {
      std::string pion = cell(idf, row, "ms2_precursor_ion");
      TEST_TRUE(pion.size() >= 2 && ionTypeOk(std::string(1, pion[0])))
      TEST_TRUE(posFinite(toD(cell(idf, row, "ms2_precursor_mass"))))
      // the MS2 fragment chosen as MS3 precursor: a fragment charge, bounded by the MS1 precursor charge
      TEST_TRUE(inFragCharge(toD(cell(idf, row, "ms2_precursor_charge")), toD(cell(idf, row, "ms1_precursor_charge"))))
      auto f3 = splitTokens(cell(idf, row, "ms3_fragments"), ';');
      auto m3 = splitTokens(cell(idf, row, "ms3_fragment_masses"), ';');
      TEST_TRUE(!f3.empty())
      TEST_EQUAL(f3.size(), m3.size())

      // I1: start_pos/end_pos must be the MS3 FRAGMENT sub-range, not the parent proteoform's full range.
      // A b/y precursor ion of index k spans exactly k residues, so end_pos - start_pos == k (the
      // ms2_precursor_ion index). Proves the sub-range was stored, with no hard-coded coordinates.
      int frag_idx = std::atoi(pion.c_str() + 1);  // pion = e.g. "b3" -> 3
      int sp = std::atoi(cell(idf, row, "start_pos").c_str());
      int ep = std::atoi(cell(idf, row, "end_pos").c_str());
      TEST_TRUE(frag_idx > 0 && (ep - sp) == frag_idx)

      // I2 (isolation-window cols): MS3 rows additionally report the MS3 fragment precursor's window.
      double w3 = toD(cell(idf, row, "ms3_isolation_width"));
      TEST_TRUE(w3 >= 2.0 && w3 < 50.0)   // engine floors the MS3 isolation window at 2.0 Da
      double s3 = toD(cell(idf, row, "ms3_window_snr"));
      TEST_TRUE(finiteVal(s3) && s3 >= 0.0 && s3 <= 1001.0)
      TEST_TRUE(posFinite(toD(cell(idf, row, "ms3_charge_intensity"))))
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// §I1b -- identification.tsv : ADR-0049 refine columns
/////////////////////////////////////////////////////////////
// Under characterization.method flashtnt (the default every committed config runs) no match carries a
// RefineCandidates record, so every identification row -- MS2 and MS3 alike -- carries the seven sentinels:
// annotation_score -1, total_shift empty, second_total_shift empty, second_score -1, cut -1, n_levels -1,
// precursor_mass_difference empty. Asserted by header name on the RAW line split at every tab, because a
// tab-splitting parser can drop an empty field and the empty cells are exactly what is being pinned. The five
// ADR-0061 columns after them carry per-row values since push 4c; new_columns_carry_their_values pins those.
START_SECTION(identification_refine_columns_are_sentinels_under_flashtnt)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_i1b");
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2);

  std::ifstream in(dir + "/identification.tsv");
  std::string line;
  // Every field, empty ones included (a stream-based split would drop a trailing empty field).
  auto split_all = [](const std::string& l) {
    std::vector<std::string> out;
    size_t b = 0;
    for (size_t e; (e = l.find('\t', b)) != std::string::npos; b = e + 1) out.push_back(l.substr(b, e - b));
    out.push_back(l.substr(b));
    return out;
  };
  std::getline(in, line);  // header
  if (!line.empty() && line.back() == '\r') line.pop_back();
  const std::vector<std::string> header = split_all(line);
  const std::vector<std::pair<std::string, std::string>> sentinels = {
    {"annotation_score", "-1"}, {"total_shift", ""}, {"second_total_shift", ""}, {"second_score", "-1"},
    {"cut", "-1"}, {"n_levels", "-1"}, {"precursor_mass_difference", ""}};
  std::vector<size_t> sentinel_idx;
  for (const auto& kv : sentinels)
  {
    const auto it = std::find(header.begin(), header.end(), kv.first);
    TEST_TRUE(it != header.end())
    if (it != header.end()) sentinel_idx.push_back(static_cast<size_t>(it - header.begin()));
  }
  const auto lvl = std::find(header.begin(), header.end(), std::string("ms_level"));
  ABORT_IF(sentinel_idx.size() != sentinels.size() || lvl == header.end())
  const size_t lvl_idx = static_cast<size_t>(lvl - header.begin());
  int n_rows = 0, n_ms2 = 0;
  while (std::getline(in, line))
  {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    ++n_rows;
    // 52 columns = 51 tabs, whatever the row type (ADR-0064: + quant_ms3; ADR-0074: + bond_coverage, explained_signal; ADR-0075: + chance_p, chance_e, confirmation_p)
    TEST_EQUAL(std::count(line.begin(), line.end(), '\t'), 51)
    const std::vector<std::string> f = split_all(line);
    ABORT_IF(f.size() != header.size())
    for (size_t k = 0; k < sentinels.size(); ++k) TEST_EQUAL(f[sentinel_idx[k]], sentinels[k].second)
    if (f[lvl_idx] == "2") ++n_ms2;
  }
  TEST_TRUE(n_rows > 0)
  TEST_TRUE(n_ms2 > 0)  // the MS2 identification row is the site Refine fills under method refine
}
END_SECTION

// The writer's formatting of the seven columns, pinned on a descriptor built by hand from Refine's decision on
// FI_REFINE_CYTC at the engine defaults (Refine_test: best S 61 at lambda 526.14, second S 46 at 655.19, cut 5,
// 5 levels kept). Refine itself is not run here -- a mass list is not a spectrum, so processScan cannot reach it;
// what is under test is that the values land in the right columns at the stream's float convention
// (std::fixed, 4 decimals, set explicitly so no neighbouring precision is inherited) and that a single
// candidate leaves second_total_shift empty and second_score -1.
START_SECTION(identification_refine_columns_carry_the_candidates)
{
  const std::string dir = freshLogDir("lf_i1c");
  Config cfg{buildJsonWithLogDir(dir, true)};
  {
    IdaLogger logger(cfg);

    IdaLogger::IdRowDescriptor row;
    row.tracking_id = "!!#";
    row.ms_level = 2;
    row.scan_mode = 'R';
    row.match.proteoform_sequence = "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";
    FragmentAnalysis::ProteoformMatch::FragmentMatch f;
    f.ion_type = "b";
    f.ion_index = 10;
    f.observed_mass = 1000.5;
    row.match.fragments.push_back(f);
    row.tic_coverage = 0.5f;
    row.refine = true;
    row.annotation_score = 61;
    row.total_shift = 526.14;
    row.has_second = true;
    row.second_total_shift = 655.19;
    row.second_score = 46;
    row.cut = 5;
    row.n_levels = 5;
    row.precursor_mass_difference = -0.98;
    // distinct values, so the two swapped columns are pinned by value, not only by header (2026-11-12)
    row.flash_extender_score = 61;
    row.match.identification_score = 429.0;
    logger.writeIdentificationRow(row);

    // the same decision with a single candidate
    row.has_second = false;
    row.second_total_shift = 655.19;  // ignored without a second candidate
    row.second_score = 46;            // ignored without a second candidate
    logger.writeIdentificationRow(row);
  }

  std::ifstream in(dir + "/identification.tsv");
  std::string header, two, one;
  std::getline(in, header);
  std::getline(in, two);
  std::getline(in, one);
  auto endsWith = [](const std::string& s, const std::string& tail) {
    return s.size() >= tail.size() && s.compare(s.size() - tail.size(), tail.size(), tail) == 0;
  };
  // tic_coverage, then the seven
  // (ADR-0061 push 4b: then its five; a hand-built row chose no protein, so protein is empty; ADR-0064: then quant_ms3;
  // ADR-0074: then bond_coverage and explained_signal, -1 on a hand-built match; ADR-0075: then chance_p, chance_e, confirmation_p, -1 likewise)
  TEST_TRUE(endsWith(header, "\ttic_coverage\tannotation_score\ttotal_shift\tsecond_total_shift\tsecond_score\tcut\tn_levels\tprecursor_mass_difference"
                             "\tprotein\traw_score\tflash_extender_score\tclosed_out\tinferred_precursor_mass\tquant_ms3\tbond_coverage\texplained_signal\tchance_p\tchance_e\tconfirmation_p"))
  TEST_TRUE(header.find("\tproteoform\tidentification_score\ttag_count\t") != std::string::npos)
  TEST_EQUAL(std::count(two.begin(), two.end(), '\t'), 51)
  TEST_EQUAL(std::count(one.begin(), one.end(), '\t'), 51)
  TEST_TRUE(endsWith(two, "\t0.5000\t61\t526.1400\t655.1900\t46\t5\t5\t-0.9800\t\t-1.0000\t61.0000\t0\t-1.0000\tnone\t-1.0000\t-1.0000\t-1\t-1\t-1"))
  TEST_TRUE(endsWith(one, "\t0.5000\t61\t526.1400\t\t-1\t5\t5\t-0.9800\t\t-1.0000\t61.0000\t0\t-1.0000\tnone\t-1.0000\t-1.0000\t-1\t-1\t-1"))
  auto field = [](const std::string& l, int k) {   // the k-th tab-separated field (0-based)
    size_t b = 0;
    for (int n = 0; n < k && b != std::string::npos; ++n) b = l.find('\t', b) == std::string::npos ? std::string::npos : l.find('\t', b) + 1;
    if (b == std::string::npos) return std::string();
    return l.substr(b, l.find('\t', b) - b);
  };
  TEST_EQUAL(field(two, 7), "429.0000")   // identification_score, next to proteoform
  TEST_EQUAL(field(one, 7), "429.0000")
}
END_SECTION

/////////////////////////////////////////////////////////////
// [REMOVED §R2t] results_tag_targeting (E4) — asserted a tag-targeted MS2 scan_results row logs
//   tag_count>0. tag_count was dropped entirely from all logs in the scan_results slim-down, so there
//   is no re-point target (identification.tsv has no tag_count column). Tag-based targeting itself
//   was later retired (ADR-0061).
/////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////
// §E6rt -- scan_commands.tsv scan_description column == the drained ScanCommand buffer
//
// X5 only inspects the DRAINED cmd.scan_description; this section closes the E6 loop by
// proving the LOGGED scan_description column (index 28) is byte-identical to the descriptor
// the engine actually drained for the same command (joined on the 3-char tracking-id prefix).
/////////////////////////////////////////////////////////////
START_SECTION(commands_scan_description_roundtrip)
{
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_e6");
  std::string json = buildJsonWithLogDir(dir, true);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  auto cycle = runFullCycle(&ida, ms1, ms2);
  TEST_TRUE(cycle.ms2_cmds.size() > 0)

  // Map each drained command's 3-char id prefix -> its raw scan_description buffer.
  std::map<std::string, std::string> drained_desc;
  for (const auto& c : cycle.ms2_cmds) { std::string d(c.scan_description); if (d.size() >= 3) drained_desc[d.substr(0, 3)] = d; }
  for (const auto& c : cycle.ms3_cmds) { std::string d(c.scan_description); if (d.size() >= 3) drained_desc[d.substr(0, 3)] = d; }

  auto t = TSVFile::parse(dir + "/scan_commands.tsv");
  ABORT_IF(t.colIndex("scan_description") != 28)  // E6: scan_description at col 28 (P5 appended precursor_id after it)

  bool matched_any = false, equal_ok = true;
  for (const auto& row : t.rows)
  {
    std::string tid = cell(t, row, "tracking_id");
    auto it = drained_desc.find(tid);
    if (it == drained_desc.end()) continue;  // MS1/AGC rows we did not capture as drained MS2/MS3
    matched_any = true;
    equal_ok = equal_ok && (cell(t, row, "scan_description") == it->second);
  }
  TEST_TRUE(matched_any)  // at least the MS2 commands must round-trip
  TEST_TRUE(equal_ok)     // logged col[28] == drained descriptor buffer, byte-for-byte
}
END_SECTION

/////////////////////////////////////////////////////////////
// §S4B -- the ADR-0061 push 4b columns carry their values (push 4c: the proteome search's)
//
// scan_results: search / search_wait_ms / search_ms. Every regular MS2 of an identifying run is searched on
// the pool (a backlog drop cannot happen here: the drive waits for the pool after every feed), so each MS2
// row says "searched" with two non-negative integer timings; MS1 and MS3 never reach the pool ("none", -1).
// identification: protein / raw_score / identification_score / closed_out / inferred_precursor_mass. The
// only protein is the inline reference, whose accession is "sequence": an MS2 self-match carries it on its
// match, an MS3 row through its MS2 context. An MS2 self-match is FLASHTnT's PrSM: its raw score is a
// -log10 probability (>= 0), and its identification score = the ExtensionScore (flash_extender_score) + 8 x
// max(0, ...), so never below it. skip_precursor_inference is on, so the inferred mass is -1 or a real mass.
// The drive authors no identification.threshold, so nothing is closed out (closed_out 0, never closed_out_skip;
// FLASHIda_ProcessScan_test's close_out_* sections drive a threshold). pooled: protein, the accession the model's
// winner chose.
//
// The §T9 drive (identifying cytC recipe + the real per-ion MS3 fixtures), because it is the one that yields
// MS2 rows, MS3 identification rows and pooled rows in one run. Each kind is required, so no check is vacuous.
/////////////////////////////////////////////////////////////
START_SECTION(new_columns_carry_their_values)
{
  auto manifest = buildMs3IonManifest("../../FlashIDA/test-data/spectra");
  TEST_TRUE(!manifest.empty())   // the MS3 rows below need the committed fixtures
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(manifest.empty() || ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_s4b");
  std::string json = buildJsonWithLogDir(dir, true, "HCD", 50);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2, &manifest);

  const std::string accession = "sequence";   // IdentificationConfig::kSequenceAccession, spelled out
  auto isNonNegInt = [](const std::string& s) {
    return !s.empty() && s.find_first_not_of("0123456789") == std::string::npos;
  };

  // --- scan_results ---------------------------------------------------------------------------------
  auto res = TSVFile::parse(dir + "/scan_results.tsv");
  int n_ms1 = 0, n_ms2 = 0, n_ms3 = 0;
  for (const auto& row : res.rows)
  {
    TEST_EQUAL(row.size(), res.headers.size())
    const int lvl = std::atoi(cell(res, row, "ms_level").c_str());
    if (lvl == 2)
    {
      ++n_ms2;
      TEST_EQUAL(cell(res, row, "search"), std::string("searched"))
      TEST_TRUE(isNonNegInt(cell(res, row, "search_wait_ms")))
      TEST_TRUE(isNonNegInt(cell(res, row, "search_ms")))
    }
    else
    {
      n_ms1 += (lvl == 1);
      n_ms3 += (lvl == 3);
      TEST_EQUAL(cell(res, row, "search"), std::string("none"))
      TEST_EQUAL(cell(res, row, "search_wait_ms"), std::string("-1"))
      TEST_EQUAL(cell(res, row, "search_ms"), std::string("-1"))
    }
  }
  TEST_TRUE(n_ms1 > 0)
  TEST_TRUE(n_ms2 > 0)
  TEST_TRUE(n_ms3 > 0)

  // --- identification -------------------------------------------------------------------------------
  auto idf = TSVFile::parse(dir + "/identification.tsv");
  int n_ms2_self = 0, n_ms3_id = 0;
  for (const auto& row : idf.rows)
  {
    TEST_EQUAL(row.size(), idf.headers.size())
    const int lvl = std::atoi(cell(idf, row, "ms_level").c_str());
    if (lvl == 2 && toD(cell(idf, row, "flash_extender_score")) >= 0.0)   // the scan's own match (-1 = re-matched)
    {
      ++n_ms2_self;
      TEST_EQUAL(cell(idf, row, "protein"), accession)
      const double extension = toD(cell(idf, row, "flash_extender_score"));
      const double identification = toD(cell(idf, row, "identification_score"));
      const double inferred = toD(cell(idf, row, "inferred_precursor_mass"));
      TEST_TRUE(toD(cell(idf, row, "raw_score")) >= 0.0)
      TEST_TRUE(identification > -1.0)
      TEST_TRUE(identification >= extension)
      TEST_EQUAL(cell(idf, row, "closed_out"), std::string("0"))
      TEST_TRUE(inferred == -1.0 || inferred > 0.0)
    }
    else if (lvl == 3)
    {
      ++n_ms3_id;
      TEST_EQUAL(cell(idf, row, "protein"), accession)
    }
  }
  TEST_TRUE(n_ms2_self > 0)
  TEST_TRUE(n_ms3_id > 0)

  // --- pooled_identification: a row is written only for an identified model ----------------------------
  auto pooled = TSVFile::parse(dir + "/pooled_identification.tsv");
  TEST_TRUE(pooled.rows.size() > 0)
  for (const auto& row : pooled.rows)
  {
    TEST_EQUAL(row.size(), pooled.headers.size())
    TEST_EQUAL(cell(pooled, row, "protein"), accession)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// §BC -- identification.tsv : bond_coverage and explained_signal (ADR-0074)
//
// The search's FLASHTnT annotates every PrSM on the MS2's raw spectrum and puts RawCoverage (a percentage) on the
// hit; the engine logs it as bond_coverage, a FRACTION, on the scan's own MS2 row. explained_signal is computed only
// for a sweep under that metric, so it is -1 on every row of this regular run. An MS3 row is a fragment match, never
// a search PrSM: both are -1. The §S4B drive.
START_SECTION(bond_coverage_is_logged_on_a_searched_ms2)
{
  auto manifest = buildMs3IonManifest("../../FlashIDA/test-data/spectra");
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(manifest.empty() || ms1.empty() || ms2.empty())

  const std::string dir = freshLogDir("lf_bc");
  std::string json = buildJsonWithLogDir(dir, true, "HCD", 50);
  FLASHIda ida(const_cast<char*>(json.c_str()));
  runFullCycle(&ida, ms1, ms2, &manifest);

  auto idf = TSVFile::parse(dir + "/identification.tsv");
  int n_ms2_self = 0, n_ms3_id = 0;
  for (const auto& row : idf.rows)
  {
    TEST_EQUAL(row.size(), idf.headers.size())
    TEST_EQUAL(cell(idf, row, "explained_signal"), std::string("-1.0000"))
    const int lvl = std::atoi(cell(idf, row, "ms_level").c_str());
    if (lvl == 2 && toD(cell(idf, row, "flash_extender_score")) >= 0.0)   // the scan's own match (-1 = re-matched)
    {
      ++n_ms2_self;
      const double coverage = toD(cell(idf, row, "bond_coverage"));
      TEST_TRUE(coverage > 0.0)
      TEST_TRUE(coverage <= 1.0)
    }
    else if (lvl == 3)
    {
      ++n_ms3_id;
      TEST_EQUAL(cell(idf, row, "bond_coverage"), std::string("-1.0000"))
    }
  }
  TEST_TRUE(n_ms2_self > 0)
  TEST_TRUE(n_ms3_id > 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// §FASTA -- a FASTA reference identifies and characterizes its protein (ADR-0061 push 4c)
//
// The §S4B drive twice: with buildJsonWithLogDir's cytC sequence inline, and with a two-entry FASTA holding that
// sequence (as sp|P00004|CYC_HORSE) and one unrelated E. coli K-12 protein (APHA_ECOLI, copied from OpenMS's K-12
// test FASTA). The search chooses cytC in both, so every identification row names its accession in the FASTA run
// and "sequence" in the other, the Identification scores are the same row for row, and the two runs plan the same
// MS3 targets. With a FASTA there is no configured sequence for any reader to fall back on, so this checks the
// chosen protein travelling from the match through the model and the MS3 context to every reader (push 4b) end to
// end. Non-vacuous: at least one identification row and one MS3 command.
/////////////////////////////////////////////////////////////
START_SECTION(a_fasta_reference_identifies_and_characterizes_its_protein)
{
  auto manifest = buildMs3IonManifest("../../FlashIDA/test-data/spectra");
  auto ms1 = loadTsvScans(CYTC_MS1);
  auto ms2 = loadTsvScans(CYTC_MS2);
  ABORT_IF(manifest.empty() || ms1.empty() || ms2.empty())

  // buildJsonWithLogDir's inline reference, spelled out: the FASTA run replaces exactly this
  const std::string cytc = "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";
  const std::string cytc_accession = "sp|P00004|CYC_HORSE";
  const std::string aphA = "MRKITQAISAVCLLFALNSSAVALASSPSPLNPGTNVARLAEQAPIHWVSVAQIENSLAGRPPMAVGFDIDDTVLFSSPGFWRGKKTFSPESEDYLKNPVFWEKMNNGWDEF"
                           "SIPKEVARQLIDMHVRRGDAIFFVTGRSPTKTETVSKTLADNFHIPATNMNPVIFAGDKPGQNTKSQWLQDKNIRIFYGDSDNDITAARDVGARGIRILRASNSTYKPLPQAG"
                           "AFGEEVIVNSEY";
  const std::string inline_ref = R"("reference": { "sequence": ")" + cytc + R"(" })";

  struct Logs
  {
    TSVFile idf;
    std::vector<std::string> ms3_mz;   // scan_commands MS3 rows' precursor_mz, in dequeue order
  };
  auto drive = [&](const std::string& tag, bool fasta) {
    const std::string dir = freshLogDir(tag);
    std::string json = buildJsonWithLogDir(dir, true, "HCD", 50);
    const auto at = json.find(inline_ref);
    TEST_TRUE(at != std::string::npos)
    if (fasta && at != std::string::npos)
    {
      const std::string path = dir + "/reference.fasta";
      std::ofstream out(path);
      out << ">" << cytc_accession << " Cytochrome c OS=Equus caballus\n" << cytc << "\n"
          << ">sp|P0AE22|APHA_ECOLI Class B acid phosphatase OS=Escherichia coli (strain K12) OX=83333 GN=aphA PE=1 SV=1\n"
          << aphA << "\n";
      out.close();
      json.replace(at, inline_ref.size(), R"("reference": { "fasta": ")" + path + R"(" })");
    }
    {
      FLASHIda ida(const_cast<char*>(json.c_str()));
      runFullCycle(&ida, ms1, ms2, &manifest);
    }   // the engine flushes and closes its streams
    Logs logs;
    logs.idf = TSVFile::parse(dir + "/identification.tsv");
    auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
    for (const auto& row : cmds.rows)
      if (std::atoi(cell(cmds, row, "ms_level").c_str()) == 3) logs.ms3_mz.push_back(cell(cmds, row, "precursor_mz"));
    return logs;
  };
  const Logs by_sequence = drive("lf_ref_seq", false);
  const Logs by_fasta = drive("lf_ref_fasta", true);

  TEST_TRUE(by_sequence.idf.rows.size() > 0)
  TEST_EQUAL(by_fasta.idf.rows.size(), by_sequence.idf.rows.size())
  for (const auto& row : by_sequence.idf.rows)
    TEST_EQUAL(cell(by_sequence.idf, row, "protein"), std::string("sequence"))
  for (const auto& row : by_fasta.idf.rows)
    TEST_EQUAL(cell(by_fasta.idf, row, "protein"), cytc_accession)
  for (Size i = 0; i < by_sequence.idf.rows.size() && i < by_fasta.idf.rows.size(); ++i)
  {
    TEST_EQUAL(cell(by_fasta.idf, by_fasta.idf.rows[i], "identification_score"),
               cell(by_sequence.idf, by_sequence.idf.rows[i], "identification_score"))
  }

  TEST_TRUE(by_sequence.ms3_mz.size() > 0)
  TEST_EQUAL(by_fasta.ms3_mz.size(), by_sequence.ms3_mz.size())
  TEST_EQUAL(by_fasta.ms3_mz == by_sequence.ms3_mz, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ADR-0047 d.9 / ADR-0048 d.4 -- scan_results.tsv: precursor_share on every MSn row, the Kinetic
// Model's fit on the one row that completes the group, sentinels everywhere else
/////////////////////////////////////////////////////////////
START_SECTION(results_kinetic_columns_only_on_the_completing_row)
{
  // Drive a kinetic_model MS2 sweep through the ENGINE on the cytC 17+ fixtures: inclusion-pinned
  // survey, six trap pre-scans at 1/11/21/31/41/50, one Orbitrap follow-up at the modeled energy.
  // Every CE 0..50 is mapped so the follow-up resolves whatever the fit decides (no fallback).
  auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc_z17.txt");
  ABORT_IF(ms1.empty())
  std::map<int, ScanData> ce_map;
  for (int ce = 0; ce <= 50; ++ce)
  {
    auto s = loadTsvScans(std::string("../../FlashIDA/test-data/spectra/") + "ms2_cytc_z17_ce" + std::to_string(ce) + ".txt");
    ABORT_IF(s.empty())
    ce_map[ce] = s[0];
  }
  const std::string dir = freshLogDir("lf_kin");
  std::string cfg = injectRuntime(kineticModelConfig(), dir);
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult a = runInterleaved(&ida, ms1, std::vector<ScanData>{ce_map.at(31)}, nullptr, 600, true, 0.0, &ce_map);
  TEST_TRUE(a.first_group_commands >= 6)

  auto res = TSVFile::parse(dir + "/scan_results.tsv");
  int completing = 0, with_source = 0, msn_rows = 0, ms1_rows = 0;
  std::string modeled_ce_cell;
  for (const auto& row : res.rows)
  {
    const std::string lvl = cell(res, row, "ms_level");
    if (lvl == "1")
    {
      ++ms1_rows;
      TEST_EQUAL(cell(res, row, "precursor_share"), std::string("-1"))
    }
    else
    {
      ++msn_rows;
      const double s = toD(cell(res, row, "precursor_share"));      // every MSn row: pre-scan AND follow-up
      TEST_TRUE(s >= 0.0 && s <= 1.0)
    }
    const std::string src = cell(res, row, "energy_source");
    if (!src.empty()) { ++with_source; TEST_TRUE(src.rfind("modeled", 0) == 0) }
    if (std::atoi(cell(res, row, "fit_points").c_str()) > 0)
    {
      ++completing;
      TEST_EQUAL(cell(res, row, "exploration_metric"), std::string("4"))
      TEST_TRUE(toD(cell(res, row, "fit_p0")) > 0.0)
      TEST_TRUE(toD(cell(res, row, "fit_width")) >= 0.25 && toD(cell(res, row, "fit_width")) <= 25.0)
      TEST_TRUE(toD(cell(res, row, "fit_rss")) >= 0.0)
      TEST_TRUE(toD(cell(res, row, "k2")) > 0.0)
      TEST_TRUE(toD(cell(res, row, "modeled_b")) > 0.0)
      TEST_TRUE(std::atoi(cell(res, row, "modeled_ce").c_str()) >= 1)
      TEST_EQUAL(cell(res, row, "winner_tracking_id"), std::string(""))   // no winner under kinetic_model
      modeled_ce_cell = cell(res, row, "modeled_ce");
    }
    else
    {
      TEST_EQUAL(cell(res, row, "modeled_ce"), std::string("0"))
      TEST_EQUAL(cell(res, row, "k2"), std::string("-1"))
      TEST_EQUAL(cell(res, row, "fit_p0"), std::string("-1"))
      TEST_EQUAL(cell(res, row, "energy_source"), std::string(""))
    }
  }
  TEST_TRUE(ms1_rows >= 1)
  TEST_TRUE(msn_rows >= 7)              // 6 pre-scans + the follow-up, at least
  TEST_EQUAL(completing, 1)
  TEST_EQUAL(with_source, 1)
  // The follow-up command carries the modeled energy: exactly one MS2 'recording' row at level 2 whose
  // collision energy equals the completing row's modeled_ce.
  auto cmds = TSVFile::parse(dir + "/scan_commands.tsv");
  int r_rows = 0;
  bool r_at_modeled = false;
  for (const auto& row : cmds.rows)
  {
    if (cell(cmds, row, "scan_type") != "recording" || cell(cmds, row, "ms_level") != "2") continue;
    ++r_rows;
    r_at_modeled = std::to_string((int)std::lround(toD(cell(cmds, row, "hcd_energy")))) == modeled_ce_cell;
  }
  TEST_EQUAL(r_rows, 1)
  TEST_TRUE(r_at_modeled)
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
