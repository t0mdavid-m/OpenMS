// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license.
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// ADR-0065 push 2: the operator view. PrecursorBoard (the model) and ConsoleView's pure renderers are checked
// directly; one cytC drive checks that the board shows what identification.tsv logs, and that the engine prints its
// run summary when it is destroyed. The terminal driver (VT, code page, scroll region) is not reachable here: ctest's
// stdout is a pipe, so ConsoleView::start returns nullptr.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ConsoleView.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorBoard.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>

#include "FLASHIda_EventCapture.h"
#include "FLASHIda_TestAccess.h"
#include "FLASHIda_TestHelpers.h"

#include <nlohmann/json.hpp>  // the status section sets identification.search_threads on the parsed recipe

#include <cmath>
#include <map>
#include <optional>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  ScanCommand ms2At(double mass, int charge)
  {
    ScanCommand c{};
    c.msn_level = 2;
    c.mono_mass = mass;
    c.num_stages = 1;
    c.stages[0].charge_state = charge;
    return c;
  }

  ScanCommand ms3()
  {
    ScanCommand c{};
    c.msn_level = 3;
    return c;
  }

  FragmentAnalysis::ProteoformMatch identifiedMatch(double score, const std::string& protein, const std::string& seq)
  {
    FragmentAnalysis::ProteoformMatch m;
    m.proteoform_sequence = seq;
    m.fragments.resize(5);
    m.search_identification_score = score;
    m.matched_protein = protein;
    return m;
  }

  const BoardRow* rowOf(const BoardSnapshot& s, double mass)
  {
    const int n = SpectralDeconvolution::getNominalMass(mass);
    for (const auto& r : s.live)
      if (r.nominal == n) return &r;
    return nullptr;
  }

  BoardRow syntheticRow(int i, PrecursorState state, size_t proforma_len)
  {
    BoardRow r;
    r.nominal = 10000 + 100 * i;
    r.mass = r.nominal + 5.0;
    r.charge = 10;
    r.state = state;
    r.score = state == PrecursorState::ClosedOut ? 500.0 + i : -1.0;
    r.ms2 = 1;
    r.last_rt_s = 1000.0 - i;
    r.accession = "P" + std::to_string(i);
    r.proforma = std::string(proforma_len, 'A');
    return r;
  }
}

START_TEST(FLASHIda_ConsoleView, "$Id$")

/////////////////////////////////////////////////////////////

// One Precursor through acquiring -> searching -> closed out; a second searched to no ID; a third returned to nothing
// that searches it. FAILS UNDER: a state that skips or regresses, or a counter that does not follow the rows.
START_SECTION(board_moves_a_precursor_through_its_states)
{
  PrecursorBoard b;
  b.setRt(60.0);
  b.onSent(ms2At(12358.2, 14), 7);
  TEST_EQUAL(rowOf(b.snapshot(180, 60, 0, 0), 12358.2)->state == PrecursorState::Acquiring, true)
  b.onReturned(12358.2, true);
  TEST_EQUAL(rowOf(b.snapshot(180, 60, 0, 0), 12358.2)->state == PrecursorState::Searching, true)
  b.onResult(12358.2, 14, 7, identifiedMatch(821.0, "CYC_HORSE", "GDVEKGKK"), true, std::optional<double>(500.0));

  b.onSent(ms2At(8441.0, 11), 8);
  b.onReturned(8441.0, true);
  b.onResult(8441.0, 11, 8, FragmentAnalysis::ProteoformMatch(), false, std::optional<double>(500.0));

  b.onSent(ms2At(6230.0, 9), 9);
  b.onReturned(6230.0, false);

  const BoardSnapshot s = b.snapshot(180, 60, 0, 0);
  const BoardRow* a = rowOf(s, 12358.2);
  const BoardRow* n = rowOf(s, 8441.0);
  const BoardRow* q = rowOf(s, 6230.0);
  ABORT_IF(a == nullptr || n == nullptr || q == nullptr)
  TEST_EQUAL(a->state == PrecursorState::ClosedOut, true)
  TEST_REAL_SIMILAR(a->score, 821.0)
  TEST_EQUAL(a->accession, "CYC_HORSE")
  TEST_EQUAL(a->charge, 14)
  TEST_EQUAL(a->ms2, 1)
  TEST_EQUAL(n->state == PrecursorState::NoId, true)
  TEST_EQUAL(q->state == PrecursorState::Acquired, true)
  TEST_EQUAL(s.precursors, 3)
  TEST_EQUAL(s.closed_out, 1)
  TEST_EQUAL(s.identified, 0)   // exclusive: the one identified row is closed out, and counted there
  TEST_EQUAL(s.no_id, 1)
  TEST_EQUAL(s.ms2, 3)
  TEST_EQUAL(s.id_rts_s.size(), 1)
}
END_SECTION

// A re-detection of one nominal mass (a second selection, a new precursor id) is one row; an MS3 counts on the row of
// its own Precursor. FAILS UNDER: a re-detection splitting into two rows, or an MS3 counted nowhere.
START_SECTION(board_merges_a_redetection_and_counts_ms3_on_its_precursor)
{
  PrecursorBoard b;
  b.setRt(30.0);
  b.onSent(ms2At(12358.20, 14), 7);
  b.setRt(90.0);
  b.onSent(ms2At(12358.31, 15), 9);   // same nominal mass
  b.onSent(ms3(), 9);
  const BoardSnapshot s = b.snapshot(180, 60, 0, 0);
  TEST_EQUAL(s.live.size(), 1)
  ABORT_IF(s.live.empty())
  TEST_EQUAL(s.live[0].ms2, 2)
  TEST_EQUAL(s.live[0].ms3, 1)
  TEST_REAL_SIMILAR(s.live[0].last_rt_s, 90.0)
  TEST_EQUAL(s.ms3, 1)
}
END_SECTION

// Every value of an identified row describes its best identification (CONTEXT.md "Precursor board"): two detections of
// one nominal mass, the first scoring 180 at 34 % coverage, the second 961 at 5 %. The row shows the second's mass,
// charge, score, coverage and ProForma; a later, lower-scoring result and a later MS2 send change none of them, and an
// MS3 that raises the first detection's coverage leaves the row at 5 %.
// FAILS UNDER: mass or charge left at the first MS2's, coverage taken from another detection, or a lower result or a
// send overwriting the row.
START_SECTION(board_row_follows_its_best_identification)
{
  PrecursorBoard b;
  b.setRt(60.0);
  b.onSent(ms2At(7304.70, 7), 5);
  b.onResult(7304.70, 7, 5, identifiedMatch(180.0, "CSPC_ECOLI", "MSGKMTGIVK"), false, std::nullopt);
  b.onPooled(5, 7301, "MS2", "!aa", 0.34);
  b.setRt(120.0);
  b.onSent(ms2At(7304.90, 8), 9);   // a re-detection 0.2 Da away: the same nominal mass, 7301
  b.onResult(7304.90, 8, 9, identifiedMatch(961.0, "CSPC_ECOLI", "MSGKMTGIVKWF"), false, std::nullopt);
  b.onPooled(9, 7301, "MS2", "!bb", 0.05);
  b.onResult(7304.70, 7, 5, identifiedMatch(200.0, "OTHER", "PEPTIDE"), false, std::nullopt);   // lower than 961
  b.onSent(ms2At(7304.70, 12), 5);                                                                // a later send
  b.onPooled(5, 7301, "b39", "!cc", 0.41);                                                         // detection 5 rises

  const BoardSnapshot s = b.snapshot(180, 60, 0, 0);
  TEST_EQUAL(s.live.size(), 1)
  ABORT_IF(s.live.size() != 1)
  const BoardRow& r = s.live[0];
  TEST_REAL_SIMILAR(r.mass, 7304.90)
  TEST_EQUAL(r.charge, 8)
  TEST_REAL_SIMILAR(r.score, 961.0)
  TEST_EQUAL(r.precursor_id, 9)
  TEST_REAL_SIMILAR(r.coverage, 0.05)
  TEST_EQUAL(r.accession, "CSPC_ECOLI")
  TEST_EQUAL(r.ms2, 3)

  // a row with no identification keeps today's behaviour: its first mass, its latest charge, no coverage
  b.onSent(ms2At(9000.6, 9), 11);
  b.onSent(ms2At(9000.8, 10), 12);   // both nominal 8996
  const BoardRow* u = rowOf(b.snapshot(180, 60, 0, 0), 9000.6);
  ABORT_IF(u == nullptr)
  TEST_REAL_SIMILAR(u->mass, 9000.6)
  TEST_EQUAL(u->charge, 10)
  TEST_REAL_SIMILAR(u->coverage, -1.0)

  // the panel draws coverage in whole percent after the score and its E-value ('-' here: no E), and '-' without a model
  const auto lines = ConsoleView::renderPanel(b.snapshot(180, 60, 0, 0), 120, 8, false, false);
  int with_cov = 0, without = 0;
  for (const auto& l : lines)
  {
    if (l.find("   961       -    5%") != std::string::npos) ++with_cov;
    if (l.find("     -     -") != std::string::npos) ++without;
  }
  TEST_EQUAL(with_cov, 1)
  TEST_EQUAL(without, 1)
}
END_SECTION

// ADR-0075: the board shows each row's E-value in its own column after the score, right-aligned under the header's E,
// one significant digit; '-' for a row without one. The board takes it from the identification it takes the score from.
// FAILS UNDER: no column, a column out of line with its header, an E kept from a lower-scoring identification.
START_SECTION(board_shows_e_after_score)
{
  PrecursorBoard b;
  b.setRt(60.0);
  auto high = identifiedMatch(212.0, "P99999", "MSGKMTGIVK");
  high.chance_log10_e = std::log10(3e-12);
  auto low = identifiedMatch(100.0, "P99999", "MSGKMTGIVK");
  low.chance_log10_e = std::log10(4e-3);
  b.onSent(ms2At(7304.70, 7), 5);
  b.onResult(7304.70, 7, 5, high, false, std::nullopt);
  b.onResult(7304.70, 7, 5, low, false, std::nullopt);   // lower: the row keeps 212 and its E
  b.onSent(ms2At(9000.6, 9), 11);                        // no identification
  const BoardSnapshot s = b.snapshot(180, 60, 0, 0);
  const BoardRow* r = rowOf(s, 7304.70);
  ABORT_IF(r == nullptr)
  TEST_EQUAL(r->log10_e.has_value(), true)
  TEST_REAL_SIMILAR(*r->log10_e, std::log10(3e-12))

  // in line with its header as the score is with "score": the same offset between a label's last letter and its values
  const auto lines = ConsoleView::renderPanel(s, 120, 8, false, false);
  const size_t score_label_end = lines[2].find("score") + 4;
  const size_t e_label_end = lines[2].find("  E  ") + 2;
  const std::string id_nominal = std::to_string(r->nominal) + " ";
  const std::string bare_nominal = std::to_string(SpectralDeconvolution::getNominalMass(9000.6)) + " ";
  int with_e = 0, without_e = 0;
  for (size_t k = 3; k < lines.size(); ++k)
  {
    const std::string& l = lines[k];
    if (l.find(id_nominal) != std::string::npos && l.find("3e-12") != std::string::npos && l.find("212") != std::string::npos)
    {
      const size_t score_end = l.find("212") + 2, e_end = l.find("3e-12") + 4;
      if (score_label_end - score_end == e_label_end - e_end) ++with_e;
    }
    if (l.find(bare_nominal) != std::string::npos)
    {
      const size_t offset = score_label_end - (l.find(" -") + 1);   // the score's '-'
      if (e_label_end - offset < l.size() && l[e_label_end - offset] == '-') ++without_e;
    }
  }
  TEST_EQUAL(with_e, 1)
  TEST_EQUAL(without_e, 1)
}
END_SECTION

// ADR-0075: the run summary's best-score line carries that row's E-value, and leaves it out when there is none.
// FAILS UNDER: no E, an E from another row, or ", E -" for a row without one.
START_SECTION(summary_best_score_carries_e)
{
  BoardSnapshot s;
  s.proteins = 1;
  s.best_score = 212.0;
  s.best_accession = "P99999";
  s.best_log10_e = std::log10(3e-12);
  std::string proteins;
  for (const auto& l : ConsoleView::renderSummary(s, false))
    if (l.rfind("proteins", 0) == 0) proteins = l;
  TEST_EQUAL(proteins.find("best score 212 (P99999, E 3e-12)") != std::string::npos, true)

  s.best_log10_e.reset();
  proteins.clear();
  for (const auto& l : ConsoleView::renderSummary(s, false))
    if (l.rfind("proteins", 0) == 0) proteins = l;
  TEST_EQUAL(proteins.find("best score 212 (P99999)") != std::string::npos, true)
}
END_SECTION

// An MS3 fold that raises a detection's Pooled coverage returns one operator line naming the MS3, its nominal mass,
// the coverage before and after, and that detection's best score; nothing else returns one.
// FAILS UNDER: a line for an MS2 baseline, an MS3-DECIDE re-pick, an unchanged or lower coverage, or a bracket
// carrying another detection's score.
START_SECTION(pooled_rise_prints_one_line)
{
  PrecursorBoard b;
  b.onResult(7304.70, 7, 5, identifiedMatch(180.0, "CSPC_ECOLI", "MSGKMTGIVK"), false, std::nullopt);
  b.onResult(7304.90, 8, 9, identifiedMatch(961.0, "CSPC_ECOLI", "MSGKMTGIVKWF"), false, std::nullopt);
  TEST_EQUAL(b.onPooled(5, 7301, "MS2", "!aa", 0.342).has_value(), false)   // the baseline
  const auto rise = b.onPooled(5, 7301, "b39", "!Xy", 0.41);
  TEST_EQUAL(rise.has_value(), true)
  ABORT_IF(!rise)
  TEST_EQUAL(*rise, "#!Xy  7301  MS3 b39  cov 34.2% -> 41.0%  (score 180)")
  TEST_EQUAL(b.onPooled(5, 7301, "y12", "!Xz", 0.41).has_value(), false)    // unchanged
  TEST_EQUAL(b.onPooled(5, 7301, "y13", "!X0", 0.40).has_value(), false)    // lower
  TEST_EQUAL(b.onPooled(5, 7301, "MS3-DECIDE", "!aa", 0.50).has_value(), false)
  TEST_EQUAL(b.onPooled(9, 7301, "b20", "!X1", 0.30).has_value(), false)    // no baseline seen for 9
}
END_SECTION

// The live set is the retention-time window, most recently active first; the run counters keep every row.
// FAILS UNDER: an expired row left on the board, or the wrong order.
START_SECTION(board_live_set_is_the_retention_time_window)
{
  PrecursorBoard b;
  b.setRt(10.0);
  b.onSent(ms2At(5000.0, 5), 1);
  b.setRt(300.0);
  b.onSent(ms2At(7000.0, 7), 2);
  b.setRt(400.0);
  b.onSent(ms2At(9000.0, 9), 3);
  const BoardSnapshot s = b.snapshot(180, 60, 0, 0);   // cut-off 220 s
  TEST_EQUAL(s.live.size(), 2)
  ABORT_IF(s.live.size() != 2)
  TEST_EQUAL(s.live[0].nominal, SpectralDeconvolution::getNominalMass(9000.0))
  TEST_EQUAL(s.live[1].nominal, SpectralDeconvolution::getNominalMass(7000.0))
  TEST_EQUAL(s.precursors, 3)
}
END_SECTION

// The panel is exactly the window: H lines, none wider than W; rows newest first; a ProForma wraps under its row,
// indented 9; the rows that do not fit are counted in the footer by state.
// FAILS UNDER: a frame wider or taller than the window, a lost row, or a footer that miscounts.
START_SECTION(render_panel_fits_the_window)
{
  BoardSnapshot s;
  s.rt_s = 1000.0;
  s.duration_min = 60.0;
  for (int i = 0; i < 12; ++i)
    s.live.push_back(syntheticRow(i, i % 2 == 0 ? PrecursorState::ClosedOut : PrecursorState::Acquiring, i % 2 == 0 ? 140 : 0));
  const int W = 100, H = 20;
  const auto lines = ConsoleView::renderPanel(s, W, H, false, false);
  TEST_EQUAL(lines.size(), static_cast<size_t>(H))
  int too_wide = 0;
  for (const auto& l : lines)
    if (ConsoleView::columns(l) > W) ++too_wide;
  TEST_EQUAL(too_wide, 0)

  // rows in order: the first row line is lines[3]; a 140-residue ProForma at wrap 91 takes 2 lines
  TEST_EQUAL(lines[3].find(std::to_string(s.live[0].nominal)) != std::string::npos, true)
  TEST_EQUAL(lines[4].substr(0, 9), std::string(9, ' '))
  TEST_EQUAL(lines[4].substr(9, 5), "AAAAA")
  TEST_EQUAL(lines[6].find(std::to_string(s.live[1].nominal)) != std::string::npos, true)

  // 16 lines for rows (20 - status - timeline - header - footer): rows of 3,1,3,1,3,1,3,1 lines = 16 -> 8 shown
  int shown = 0;
  for (const auto& r : s.live)
    for (size_t k = 3; k + 1 < lines.size(); ++k)
      if (lines[k].rfind(" ", 0) == 0 && lines[k].find(" " + std::to_string(r.nominal) + " ") != std::string::npos) { ++shown; break; }
  TEST_EQUAL(shown, 8)
  TEST_EQUAL(lines.back(), " + 4 more live (2 acquiring | 2 closed out)")
  TEST_EQUAL(lines[1].back(), '|')   // the timeline ends at its frame: the status line carries the counts
}
END_SECTION

// ASCII glyphs mean ASCII bytes, for a console whose code page could not be switched to UTF-8.
// FAILS UNDER: a UTF-8 glyph leaking into the fallback.
START_SECTION(render_ascii_has_no_non_ascii_bytes)
{
  BoardSnapshot s;
  s.rt_s = 1800.0;
  s.duration_min = 60.0;
  s.identified = 2;
  s.id_rts_s = {600.0, 620.0};
  s.search_wait_ms = {3, 5};
  s.live.push_back(syntheticRow(0, PrecursorState::ClosedOut, 30));
  int high = 0;
  for (const auto& l : ConsoleView::renderPanel(s, 120, 20, false, false))
    for (unsigned char c : l) high += c >= 0x80;
  for (const auto& l : ConsoleView::renderSummary(s, false))
    for (unsigned char c : l) high += c >= 0x80;
  TEST_EQUAL(high, 0)
}
END_SECTION

// Each cell spans duration / cells minutes and is shaded relative to the busiest cell; cells after now are blank.
// FAILS UNDER: identifications in the wrong cell, a busiest cell not drawn full, or the future drawn.
START_SECTION(timeline_shades_by_identifications_and_stops_at_now)
{
  BoardSnapshot s;
  s.duration_min = 60.0;
  s.rt_s = 31.0 * 60.0;
  s.id_rts_s = {600.0, 600.0, 600.0, 660.0, 1800.0};   // minutes 10, 10, 10, 11, 30
  const std::string t = ConsoleView::timelineCells(s, 40, false);   // 1.5 min per cell
  TEST_EQUAL(t.size(), 40)
  ABORT_IF(t.size() != 40)
  TEST_EQUAL(t[6], '@')    // minute 10: the busiest cell, full
  TEST_EQUAL(t[7] != ' ' && t[7] != '.', true)   // minute 11
  TEST_EQUAL(t[20] != ' ' && t[20] != '.', true)   // minute 30
  TEST_EQUAL(t[0], ' ')    // no identification, in the past
  TEST_EQUAL(t.substr(21), std::string(19, '.'))   // after minute 31
}
END_SECTION

// The summary's search line: nearest-rank median and 95th percentile of the search times, then of the waits; the
// precursors line counts identified and closed out apart; the timeline carries no count.
// FAILS UNDER: another percentile definition, the waits shown as the times, a missing line, or a closed-out row
// counted as identified too.
START_SECTION(summary_reports_search_wait_percentiles)
{
  BoardSnapshot s;
  for (int i = 100; i >= 1; --i) s.search_wait_ms.push_back(i);
  for (int i = 100; i >= 1; --i) s.search_ms.push_back(10 * i);
  s.precursors = 10;
  s.identified = 3;
  s.closed_out = 2;
  s.no_id = 5;
  std::string search, precursors, timeline;
  for (const auto& l : ConsoleView::renderSummary(s, false))
  {
    if (l.rfind("search", 0) == 0) search = l;
    if (l.rfind("precursors", 0) == 0) precursors = l;
    if (l.rfind("timeline", 0) == 0) timeline = l;
  }
  TEST_EQUAL(search, "search      median 500 ms | p95 950 ms | wait median 50 ms | p95 95 ms | n 100")
  TEST_EQUAL(precursors, "precursors  10 | identified 3 | closed out 2 | no-id 5")
  TEST_EQUAL(timeline.empty(), false)
  TEST_EQUAL(timeline.back(), '|')
}
END_SECTION

// The status line's `searching` is the searches WAITING against the backlog's capacity, read from the pool: one search
// parked at the gate, two more fed, so two wait while the capacity is 4 (one worker). No pool, no item.
// FAILS UNDER: the snapshot carrying the capacity as the waiting count (the bug this replaced: `backlog` was
// SearchPool::backlog(), a constant), or the item drawn without a pool.
START_SECTION(status_counts_waiting_searches_not_capacity)
{
  BoardSnapshot plain;
  plain.search_waiting = 3;
  plain.search_capacity = 8;
  TEST_EQUAL(ConsoleView::renderPanel(plain, 200, 6, false, false)[0].find("searching 3/8") != std::string::npos, true)
  plain.search_capacity = 0;
  TEST_EQUAL(ConsoleView::renderPanel(plain, 200, 6, false, false)[0].find("searching") == std::string::npos, true)

  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  nlohmann::json j = nlohmann::json::parse(buildJsonWithLogDir(freshLogDir("console_view_waiting"), true));
  j["identification"]["search_threads"] = 1;
  const std::string cfg = j.dump();
  FLASHIda ida(const_cast<char*>(cfg.c_str()));
  AcqResult acq = runInterleaved(&ida, ms1_scans, {});
  ABORT_IF(acq.ms2_cmds.size() < 3)

  SearchGate gate;
  FLASHIdaTestAccess::onBeforeSearch(ida, gate.hook());
  const auto& ms2 = ms2_scans[0];
  auto feed = [&](const ScanCommand& c) {
    ida.processScan(ms2.mzs.data(), ms2.ints.data(), (int)ms2.mzs.size(), ms2.rt, 2, c.scan_description, 0.0,
                    instrumentScanNumberOf(ms2));
  };
  feed(acq.ms2_cmds[0]);   // taken by the worker, parked before its search
  const bool parked = gate.reached();
  if (parked)
  {
    feed(acq.ms2_cmds[1]);
    feed(acq.ms2_cmds[2]);
  }
  const BoardSnapshot held = FLASHIdaTestAccess::board(ida);
  gate.release();
  ida.waitForSearchIdle();
  ABORT_IF(!parked)
  TEST_EQUAL(held.search_waiting, 2)
  TEST_EQUAL(held.search_capacity, 4)
  TEST_EQUAL(FLASHIdaTestAccess::board(ida).search_waiting, 0)
}
END_SECTION

// Logged value = decision value: a cytC drive's board shows an identified row whose ProForma, score, mass and charge
// are one MS2 row of identification.tsv for the same nominal mass.
// FAILS UNDER: the board rendering its own ProForma, or keeping a different score, mass or charge than the log.
START_SECTION(engine_board_matches_the_identification_log)
{
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  const std::string dir = freshLogDir("console_view_board");
  BoardSnapshot snap;
  {
    std::string json = buildJsonWithLogDir(dir, true);
    FLASHIda ida(const_cast<char*>(json.c_str()));
    runFullCycle(&ida, ms1_scans, ms2_scans);
    snap = FLASHIdaTestAccess::board(ida);
  }
  const BoardRow* id_row = nullptr;
  for (const auto& r : snap.live)
    if ((r.state == PrecursorState::Identified || r.state == PrecursorState::ClosedOut) && !r.proforma.empty()) { id_row = &r; break; }
  TEST_EQUAL(id_row != nullptr, true)
  ABORT_IF(id_row == nullptr)

  auto tsv = TSVFile::parse(dir + "/identification.tsv");
  int matching = 0;
  for (const auto& row : tsv.rows)
  {
    if (cell(tsv, row, "ms_level") != "2") continue;
    if (SpectralDeconvolution::getNominalMass(std::stod(cell(tsv, row, "ms1_precursor_mass"))) != id_row->nominal) continue;
    if (cell(tsv, row, "proteoform") != id_row->proforma) continue;
    if (std::fabs(std::stod(cell(tsv, row, "identification_score")) - id_row->score) >= 1e-3) continue;
    if (std::fabs(std::stod(cell(tsv, row, "ms1_precursor_mass")) - id_row->mass) >= 1e-3) continue;
    if (std::stoi(cell(tsv, row, "ms1_precursor_charge")) != id_row->charge) continue;
    ++matching;
  }
  TEST_EQUAL(matching > 0, true)
  TEST_EQUAL(snap.ms2 > 0, true)
}
END_SECTION

// Logged value = decision value: on a cytC MS2+MS3 drive, every operator coverage line is one pooled row whose MS3
// fold raised its precursor id's coverage_pct, and a live row's coverage is the last coverage_pct of its precursor id.
// The drive produces such rises (13 on this fixture), so the section cannot pass on none.
// FAILS UNDER: the observer not wired, a line for a fold that did not raise coverage, a missing line, or the board's
// coverage taken from another precursor id than the one whose identification it shows.
START_SECTION(ms3_coverage_rise_reaches_board_and_console)
{
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  auto ms2_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
  ABORT_IF(ms1_scans.empty() || ms2_scans.empty())
  const std::string dir = freshLogDir("console_view_coverage");
  EventCapture cap;
  BoardSnapshot snap;
  {
    std::string json = buildJsonWithLogDir(dir, true);
    FLASHIda ida(const_cast<char*>(json.c_str()));
    runFullCycle(&ida, ms1_scans, ms2_scans);
    snap = FLASHIdaTestAccess::board(ida);
  }
  const std::string events = cap.str();

  auto pooled = TSVFile::parse(dir + "/pooled_identification.tsv");
  std::map<int, double> last;
  int rises = 0, found = 0;
  for (const auto& row : pooled.rows)
  {
    const int pid = std::stoi(cell(pooled, row, "precursor_id"));
    const double cov = std::stod(cell(pooled, row, "coverage_pct"));
    const std::string trigger = cell(pooled, row, "trigger");
    const auto before = last.find(pid);
    if (trigger != "MS2" && trigger != "MS3-DECIDE" && before != last.end() && cov > before->second + 1e-9)
    {
      ++rises;
      const std::string head = "#" + cell(pooled, row, "trigger_scan_id") + "  " + cell(pooled, row, "nominal_mass")   // parse unquotes ids
                               + "  MS3 " + trigger + "  cov ";
      if (events.find(head) != std::string::npos) ++found;
    }
    last[pid] = cov;
  }
  int lines = 0;
  for (size_t p = events.find("[OP-MS3] #"); p != std::string::npos; p = events.find("[OP-MS3] #", p + 1))
  {
    const size_t eol = events.find('\n', p);
    if (events.substr(p, eol - p).find("  cov ") != std::string::npos) ++lines;
  }
  TEST_EQUAL(rises > 0, true)
  TEST_EQUAL(found, rises)
  TEST_EQUAL(lines, rises)

  int with_coverage = 0, agree = 0;
  for (const auto& r : snap.live)
  {
    if (r.coverage < 0.0) continue;
    ++with_coverage;
    const auto l = last.find(r.precursor_id);
    if (l != last.end() && std::fabs(l->second - r.coverage) < 1e-4) ++agree;
  }
  TEST_EQUAL(with_coverage > 0, true)
  TEST_EQUAL(agree, with_coverage)
}
END_SECTION

// The engine prints its run summary when it is destroyed, through the event log, in any mode.
// FAILS UNDER: no summary, or a summary that bypasses the event log.
START_SECTION(engine_destructor_prints_the_run_summary)
{
  auto ms1_scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  ABORT_IF(ms1_scans.empty())
  EventCapture cap;
  {
    std::string json = buildJsonWithLogDir(freshLogDir("console_view_summary"));
    FLASHIda ida(const_cast<char*>(json.c_str()));
    runFullCycle(&ida, ms1_scans, {});
  }
  const std::string events = cap.str();
  TEST_EQUAL(events.find("[OP-SUMMARY] scans") != std::string::npos, true)
  TEST_EQUAL(events.find("[OP-SUMMARY] precursors") != std::string::npos, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
