// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

// KineticModel (ADR-0047): the pure functions behind the kinetic_model exploration metric --
// precursor share, the bounded Gompertz fit, the fragment-pool recursion and the modeled maximum.
// Layer A of the validation: the pure sections below pin the math on synthetic curves; the table
// sections (Task 4) replay the study's 108 sweeps against the vendored pipeline's expected values.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/KineticModel.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::KineticModel;

namespace
{
  struct Row { std::vector<std::string> cells; };

  // Reads a '#'-headed TSV: the first '#' line is the provenance, the second the column list; data rows follow.
  std::vector<Row> readKineticTsv(const std::string& path, std::string& provenance)
  {
    std::ifstream in(path);
    std::vector<Row> rows;
    std::string line;
    int hashes = 0;
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (!line.empty() && line[0] == '#')
      {
        if (hashes++ == 0) provenance = line;
        continue;
      }
      if (line.empty()) continue;
      Row r;
      std::stringstream ss(line);
      std::string c;
      while (std::getline(ss, c, '\t')) r.cells.push_back(c);
      rows.push_back(r);
    }
    return rows;
  }

  double toD(const std::string& s) { return s.empty() ? std::numeric_limits<double>::quiet_NaN() : std::stod(s); }

  // The unique sweep key is (dataset, charge, sweep): two runs of one protein repeat charge + sweep index.
  std::string sweepKey(const Row& r) { return r.cells[1] + "|" + r.cells[2] + "|" + r.cells[3]; }

  const std::string kShares   = "../../FlashIDA/test-data/kinetic/sweep_shares.tsv";
  const std::string kExpected = "../../FlashIDA/test-data/kinetic/kinetic_expected.tsv";
} // namespace

START_TEST(KineticModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(precursorShare_is_window_over_whole_scan)
{
  const double mz[] = {790.0, 800.0, 810.0, 900.0};
  const double in[] = {100.0, 500.0, 200.0, 200.0};
  TEST_REAL_SIMILAR(precursorShare(mz, in, 4, 799.6, 800.4), 0.5)     // 500 / 1000
  TEST_REAL_SIMILAR(precursorShare(mz, in, 4, 700.0, 1000.0), 1.0)    // everything in the window
  TEST_REAL_SIMILAR(precursorShare(mz, in, 4, 850.0, 860.0), 0.0)     // nothing in the window, total > 0
}
END_SECTION

START_SECTION(precursorShare_returns_minus_one_on_an_empty_scan)
{
  // Review Focus 2: an empty trap read must be "no measurement" (-1), never a share of 0.
  const double mz[] = {800.0};
  const double zero[] = {0.0};
  TEST_REAL_SIMILAR(precursorShare(mz, zero, 1, 799.0, 801.0), -1.0)
  TEST_REAL_SIMILAR(precursorShare(nullptr, nullptr, 0, 799.0, 801.0), -1.0)
}
END_SECTION

START_SECTION(fitGompertz_recovers_an_exact_curve_from_six_points)
{
  // P0 0.9, inflection 31, width 3.3 at the study's offset-0 grid 1, 11, 21, 31, 41, 50.
  std::vector<SharePoint> pts;
  for (double ce : {1.0, 11.0, 21.0, 31.0, 41.0, 50.0})
    pts.push_back({ce, 0.9 * std::exp(-std::exp((ce - 31.0) / 3.3))});
  GompertzFit f = fitGompertz(pts);
  TEST_EQUAL(f.ok, true)
  TEST_EQUAL(f.n_points, 6)
  // A lattice fit: the inflection lattice is 0.5 NCE coarse then 0.02 fine, the width lattice ~12 % coarse then
  // ~1.3 % fine, so the exact curve is recovered to that resolution, not to machine precision.
  TEST_EQUAL(std::abs(f.p0 - 0.9) <= 0.005, true)
  TEST_EQUAL(std::abs(f.ce_max_s - 31.0) <= 0.05, true)
  TEST_EQUAL(std::abs(f.width - 3.3) <= 0.05, true)
  TEST_EQUAL(f.rss < 1e-5, true)
  // The curve's own readings at its fitted parameters.
  TEST_REAL_SIMILAR(f.share(f.ce_max_s), f.p0 / std::exp(1.0))         // the inflection leaves P0/e
  TEST_REAL_SIMILAR(f.rate(f.ce_max_s), (f.p0 / std::exp(1.0)) / f.width)   // s = P/w * e^0 at the inflection
}
END_SECTION

START_SECTION(fitGompertz_needs_three_finite_points)
{
  std::vector<SharePoint> two = {{1.0, 0.9}, {50.0, 0.0}};
  TEST_EQUAL(fitGompertz(two).ok, false)
  std::vector<SharePoint> three_with_junk = {{1.0, 0.9}, {11.0, -1.0}, {21.0, 0.6}, {31.0, 0.3},
                                             {41.0, std::numeric_limits<double>::quiet_NaN()}};
  GompertzFit f = fitGompertz(three_with_junk);
  TEST_EQUAL(f.ok, true)
  TEST_EQUAL(f.n_points, 3)                    // -1 and NaN are excluded, never counted
}
END_SECTION

START_SECTION(fitGompertz_respects_the_bounds)
{
  // Shares read on the falling flank only: a free fit would inflate P0 above 1. The cap is 1.1 x 0.5.
  std::vector<SharePoint> flank = {{30.0, 0.5}, {35.0, 0.3}, {40.0, 0.12}, {45.0, 0.03}, {50.0, 0.005}};
  GompertzFit f = fitGompertz(flank);
  TEST_EQUAL(f.ok, true)
  TEST_EQUAL(f.p0 <= 1.1 * 0.5 + 1e-9, true)
  TEST_EQUAL(f.width >= kWidthMin && f.width <= kWidthMax, true)
}
END_SECTION

START_SECTION(fragmentPool_matches_the_recursion)
{
  std::vector<double> s = {0.0, 0.1, 0.2, 0.1, 0.0};
  const double k2 = 0.0388, a = std::exp(-k2), g = (1.0 - std::exp(-k2)) / k2;
  auto B = fragmentPool(s, k2);
  TEST_EQUAL((int)B.size(), 5)
  double b = 0.0;
  for (size_t k = 0; k < s.size(); ++k)
  {
    b = a * b + s[k] * g;
    TEST_REAL_SIMILAR(B[k], b)
  }
}
END_SECTION

START_SECTION(modeledMaximum_is_an_integer_argmax_with_edge_flags)
{
  std::vector<SharePoint> pts;
  for (double ce : {1.0, 11.0, 21.0, 31.0, 41.0, 50.0})
    pts.push_back({ce, 0.9 * std::exp(-std::exp((ce - 31.0) / 3.3))});
  ModeledMaximum m = modeledMaximum(fitGompertz(pts), 0.0388, 50);
  TEST_EQUAL(m.ce >= 32 && m.ce <= 38, true)   // past the inflection, before the pool decays away
  TEST_EQUAL(m.at_lower_edge, false)
  TEST_EQUAL(m.at_upper_edge, false)
  // A curve whose inflection sits far beyond the range peaks at the top edge (Review Focus 1).
  std::vector<SharePoint> late;
  for (double ce : {1.0, 11.0, 21.0, 31.0, 41.0, 50.0})
    late.push_back({ce, 0.9 * std::exp(-std::exp((ce - 90.0) / 3.3))});
  ModeledMaximum e = modeledMaximum(fitGompertz(late), 0.0388, 50);
  TEST_EQUAL(e.ce, 50)
  TEST_EQUAL(e.at_upper_edge, true)
  GompertzFit none;
  TEST_EQUAL(modeledMaximum(none, 0.0388, 50).ce, 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Layer A against the study's own data (ADR-0047 decision 10): the committed tables under
// FlashIDA/test-data/kinetic/ were generated by the vendored pipeline (test-scripts/kinetic/).
// Helpers live above START_TEST (which opens main()).
/////////////////////////////////////////////////////////////

START_SECTION(the_committed_tables_carry_their_provenance)
{
  std::string prov;
  auto shares = readKineticTsv(kShares, prov);
  TEST_EQUAL(prov.find("make_kinetic_tables.py") != std::string::npos, true)
  TEST_EQUAL(prov.find("corpus eed88d7") != std::string::npos, true)
  // 108 sweeps; 107 run to NCE 50 and alphaLactalbumin 10+ sweep 2 ends at 41.
  std::map<std::string, int> per_sweep;
  for (const auto& r : shares) ++per_sweep[sweepKey(r)];
  TEST_EQUAL((int)per_sweep.size(), 108)
  TEST_EQUAL((int)shares.size(), 107 * 50 + 41)
  auto expected = readKineticTsv(kExpected, prov);
  TEST_EQUAL(prov.find("make_kinetic_tables.py") != std::string::npos, true)
  TEST_EQUAL((int)expected.size(), 108 * 10)
}
END_SECTION

START_SECTION(the_engine_fit_reproduces_the_studys_modeled_maxima)
{
  // Summary: 1080 fits (108 sweeps x 10 grid offsets at spacing 10, top edge always acquired) against the
  // vendored pipeline's grid-search fit. Purpose: layer A of ADR-0047 d.10 -- the engine's fit (the lattice port, after the d.3 amendment) is "equivalent" under
  // the tolerances the owner set (Q15): CE* within +-1 on EVERY fit, identical on >= 95 %, mean regret within
  // 0.005 of the pipeline's, parameters within the grid's own resolution on >= 95 %.
  std::string prov;
  auto shares   = readKineticTsv(kShares, prov);
  auto expected = readKineticTsv(kExpected, prov);
  std::map<std::string, std::map<int, std::pair<double, double>>> curves;   // key -> ce -> (share, reading)
  for (const auto& r : shares)
    curves[sweepKey(r)][std::stoi(r.cells[4])] = {toD(r.cells[5]), toD(r.cells[6])};

  int n = 0, both_fit = 0, identical = 0, within1 = 0, p0_ok = 0, c_ok = 0, w_ok = 0;
  double regret_engine = 0.0, regret_pipeline = 0.0;
  int worst_gap = 0;
  std::string worst_key;
  for (const auto& e : expected)
  {
    const int off = std::stoi(e.cells[4]);
    const auto& curve = curves.at(sweepKey(e));
    const int ce_max = curve.rbegin()->first;
    std::vector<SharePoint> pts;
    for (const auto& [ce, v] : curve)                                              // pipeline.acquired
      if ((((ce - 1 - off) % 10 == 0) && ce >= 1 + off) || ce == ce_max) pts.push_back({(double)ce, v.first});
    GompertzFit f = fitGompertz(pts);
    const double k2 = toD(e.cells[10]);
    ModeledMaximum m = modeledMaximum(f, k2, ce_max);
    const int ce_star_pipeline = std::stoi(e.cells[11]);
    ++n;
    if (ce_star_pipeline == 0)
    {
      TEST_EQUAL(m.ce, 0)                                                            // the pipeline had no fit either
      continue;
    }
    ++both_fit;
    TEST_EQUAL(f.ok, true)
    const int gap = std::abs(m.ce - ce_star_pipeline);
    if (gap == 0) ++identical;
    if (gap <= 1) ++within1;
    if (gap > worst_gap) { worst_gap = gap; worst_key = sweepKey(e) + " off " + e.cells[4]; }
    if (gap > 1)
      STATUS("gap " << gap << " " << e.cells[0] << " z" << e.cells[2] << " sw" << e.cells[3] << " off" << off
             << " engine p0/c/w/rss/ce " << f.p0 << " " << f.ce_max_s << " " << f.width << " " << f.rss << " " << m.ce
             << " pipeline " << e.cells[6] << " " << e.cells[7] << " " << e.cells[8] << " " << e.cells[9] << " " << ce_star_pipeline
             << " points " << pts.size())
    // regret of the ENGINE's pick, read off the committed reading column exactly as pipeline.score does
    double best = -1.0;
    for (const auto& [ce, v] : curve) if (std::isfinite(v.second)) best = std::max(best, v.second);
    const auto it = curve.find(m.ce);
    const double read_at = (it != curve.end()) ? it->second.second : std::numeric_limits<double>::quiet_NaN();
    regret_engine   += std::isfinite(read_at) ? 1.0 - read_at / best : 1.0;
    regret_pipeline += toD(e.cells[12]);
    if (std::abs(f.p0 - toD(e.cells[6])) <= 0.02) ++p0_ok;
    if (std::abs(f.ce_max_s - toD(e.cells[7])) <= 0.5) ++c_ok;
    if (std::abs(f.width / toD(e.cells[8]) - 1.0) <= 0.15) ++w_ok;
  }
  STATUS("fits " << n << ", both fit " << both_fit << ", identical " << identical << ", within 1 " << within1
         << ", worst gap " << worst_gap << " at " << worst_key
         << ", mean regret engine " << regret_engine / both_fit << " pipeline " << regret_pipeline / both_fit
         << ", p0/c/w within grid resolution " << p0_ok << "/" << c_ok << "/" << w_ok)
  TEST_EQUAL(n, 1080)
  TEST_EQUAL(within1, both_fit)                                                      // every fit within +-1 NCE
  TEST_EQUAL(identical * 100 >= both_fit * 95, true)                                 // >= 95 % identical
  TEST_EQUAL(std::abs(regret_engine / both_fit - regret_pipeline / both_fit) <= 0.005, true)
  TEST_EQUAL(std::abs(regret_pipeline / both_fit - 0.1485) <= 0.003, true)           // the report's six-scan arm
  TEST_EQUAL(p0_ok * 100 >= both_fit * 95, true)
  TEST_EQUAL(c_ok * 100 >= both_fit * 95, true)
  TEST_EQUAL(w_ok * 100 >= both_fit * 95, true)
}
END_SECTION

START_SECTION(ubiquitin_9_plus_at_offset_zero_is_the_reports_example)
{
  std::string prov;
  auto shares = readKineticTsv(kShares, prov);
  std::vector<SharePoint> pts;
  for (const auto& r : shares)
    if (r.cells[0] == "Ubiquitin" && r.cells[2] == "9" && r.cells[3] == "1")      // one ubiquitin run, so protein is unique here
    {
      const int ce = std::stoi(r.cells[4]);
      if (((ce - 1) % 10 == 0) || ce == 50) pts.push_back({(double)ce, toD(r.cells[5])});
    }
  TEST_EQUAL((int)pts.size(), 6)                        // 1, 11, 21, 31, 41, 50
  GompertzFit f = fitGompertz(pts);
  TEST_EQUAL(f.ok, true)
  TEST_EQUAL(std::abs(f.p0 - 0.906) <= 0.02, true)
  TEST_EQUAL(std::abs(f.ce_max_s - 31.08) <= 0.5, true)
  TEST_EQUAL(std::abs(f.width - 3.32) <= 0.5, true)
  ModeledMaximum m = modeledMaximum(f, 0.039046, 50);   // ubiquitin's fold k2, every protein once
  TEST_EQUAL(m.ce == 35 || m.ce == 36, true)            // report: 35 from six scans, 36 from fifty, measured 37
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
