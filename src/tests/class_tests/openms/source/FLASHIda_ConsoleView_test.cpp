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

#include "FLASHIda_TestAccess.h"
#include "FLASHIda_TestHelpers.h"

#include <cmath>
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
  TEST_EQUAL(rowOf(b.snapshot(180, 60, 0), 12358.2)->state == PrecursorState::Acquiring, true)
  b.onReturned(12358.2, true);
  TEST_EQUAL(rowOf(b.snapshot(180, 60, 0), 12358.2)->state == PrecursorState::Searching, true)
  b.onResult(12358.2, identifiedMatch(821.0, "CYC_HORSE", "GDVEKGKK"), true, std::optional<double>(500.0));

  b.onSent(ms2At(8441.0, 11), 8);
  b.onReturned(8441.0, true);
  b.onResult(8441.0, FragmentAnalysis::ProteoformMatch(), false, std::optional<double>(500.0));

  b.onSent(ms2At(6230.0, 9), 9);
  b.onReturned(6230.0, false);

  const BoardSnapshot s = b.snapshot(180, 60, 0);
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
  TEST_EQUAL(s.identified, 1)
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
  const BoardSnapshot s = b.snapshot(180, 60, 0);
  TEST_EQUAL(s.live.size(), 1)
  ABORT_IF(s.live.empty())
  TEST_EQUAL(s.live[0].ms2, 2)
  TEST_EQUAL(s.live[0].ms3, 1)
  TEST_REAL_SIMILAR(s.live[0].last_rt_s, 90.0)
  TEST_EQUAL(s.ms3, 1)
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
  const BoardSnapshot s = b.snapshot(180, 60, 0);   // cut-off 220 s
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

// The summary's search line: nearest-rank median and 95th percentile of the waits.
// FAILS UNDER: another percentile definition, or the line missing.
START_SECTION(summary_reports_search_wait_percentiles)
{
  BoardSnapshot s;
  for (int i = 100; i >= 1; --i) s.search_wait_ms.push_back(i);
  std::string search;
  for (const auto& l : ConsoleView::renderSummary(s, false))
    if (l.rfind("search", 0) == 0) search = l;
  TEST_EQUAL(search.find("median 50 ms") != std::string::npos, true)
  TEST_EQUAL(search.find("p95 95 ms") != std::string::npos, true)
}
END_SECTION

// Logged value = decision value: a cytC drive's board shows an identified row whose ProForma and score are an MS2 row
// of identification.tsv for the same nominal mass.
// FAILS UNDER: the board rendering its own ProForma, or keeping a different score than the log.
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
    if (std::fabs(std::stod(cell(tsv, row, "identification_score")) - id_row->score) < 1e-3) ++matching;
  }
  TEST_EQUAL(matching > 0, true)
  TEST_EQUAL(snap.ms2 > 0, true)
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
