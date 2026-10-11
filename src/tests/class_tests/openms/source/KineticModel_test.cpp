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
  TEST_REAL_SIMILAR(f.p0, 0.9)
  TEST_REAL_SIMILAR(f.ce_max_s, 31.0)
  TEST_REAL_SIMILAR(f.width, 3.3)
  TEST_EQUAL(f.rss < 1e-8, true)
  TEST_REAL_SIMILAR(f.share(31.0), 0.9 / std::exp(1.0))                // the inflection leaves P0/e
  TEST_REAL_SIMILAR(f.rate(31.0), (0.9 / std::exp(1.0)) / 3.3)         // s = P/w * e^0 at the inflection
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
/////////////////////////////////////////////////////////////
END_TEST
