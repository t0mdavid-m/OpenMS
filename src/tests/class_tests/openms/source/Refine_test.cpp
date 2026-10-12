// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

// Refine (ADR-0049): the annotation method beside FLASHTnT. Detect -- the shift scan (the curve),
// Expected, the shuffled-sequence FDR and the levels -- is pinned here on synthetic inputs; the
// rules are the owner's prototype (trunc_core.py / detect_fdr.py), one to one.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Refine.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  struct MassList
  {
    std::vector<double> masses, intensities;
    std::vector<int> charges;
  };

  // Reads a deconvolved mass list `mass\tintensity\t#charges`; '#' comment lines and a non-numeric
  // header line are skipped.
  [[maybe_unused]] MassList loadMassList(const std::string& path)
  {
    MassList out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty() || line[0] == '#') continue;
      std::stringstream ss(line);
      std::string m, i, z;
      std::getline(ss, m, '\t');
      std::getline(ss, i, '\t');
      std::getline(ss, z, '\t');
      char* end = nullptr;
      const double mass = std::strtod(m.c_str(), &end);
      if (end == m.c_str()) continue; // the header row
      out.masses.push_back(mass);
      out.intensities.push_back(i.empty() ? 0.0 : std::stod(i));
      out.charges.push_back(z.empty() ? 0 : std::stoi(z));
    }
    return out;
  }

  std::vector<double> concat(const std::vector<double>& a, const std::vector<double>& b)
  {
    std::vector<double> out(a);
    out.insert(out.end(), b.begin(), b.end());
    return out;
  }
} // namespace

START_TEST(Refine, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(curve_counts_each_ion_once_and_respects_tolerance)
{
  // grid -1..+1 Da in steps of 0.001 Da: index 1000 is dm = 0, index 1030 is dm = +0.03
  const std::vector<double> theo = {1000.0, 2000.0};
  const std::vector<double> masses = {1000.005, 1000.009, 2000.0};
  const std::vector<int> f = Refine::curve(theo, masses, 10.0, -1.0, 0.001, 2001);
  TEST_EQUAL(f.size(), 2001)
  // ion 1000 matches two masses but counts once; ion 2000 matches -> 2 (3 if the ion counted per mass)
  TEST_EQUAL(f[1000], 2)
  // +0.03: 1000.03 is 0.021 Da from 1000.009 (tol 0.010), 2000.03 is 0.03 Da from 2000.0 (tol 0.020)
  TEST_EQUAL(f[1030], 0)
  // +0.02: 2000.02 sits exactly tol(2000.0) = 0.020 Da away -- the edge counts; 1000.02 is 0.011 Da out
  TEST_EQUAL(f[1020], 1)
  TEST_EQUAL(*std::max_element(f.begin(), f.end()), 2)
}
END_SECTION

START_SECTION(region_counts_one_pass_equals_naive)
{
  const std::vector<int> f = {0, 1, 3, 3, 2, 0, 4, 1, 0};
  const std::vector<int> n = Refine::regionCounts(f, 5);
  TEST_EQUAL(n.size(), 6) // index h = 1..hmax
  TEST_EQUAL(n[1], 2)
  TEST_EQUAL(n[2], 2)
  TEST_EQUAL(n[3], 2)
  TEST_EQUAL(n[4], 1)
  TEST_EQUAL(n[5], 0)
}
END_SECTION

START_SECTION(expected_removes_internal_coincidences)
{
  const Refine::Ladder lad = Refine::ladder("PEPTIDEK", {"b", "y"});
  const double b3 = lad.prefix[2], b5 = lad.prefix[4];
  const double dm = b5 - b3; // the mass of residues 4..5 (T + I)
  std::vector<double> masses = {b3, b5, dm};
  std::sort(masses.begin(), masses.end());
  const std::vector<double> T = concat(lad.prefix, lad.suffix);
  const std::vector<bool> ex = Refine::explained(T, masses, 10.0);
  std::vector<double> m0;
  for (size_t i = 0; i < masses.size(); ++i)
    if (ex[i]) m0.push_back(masses[i]);
  TEST_EQUAL(m0.size(), 2) // b3 and b5; the internal mass is not a ladder ion
  // grid dm-1 .. dm+1 in steps of 0.001: index 1000 is dm
  const std::vector<int> O = Refine::curve(T, masses, 10.0, dm - 1.0, 0.001, 2001);
  const std::vector<int> E = Refine::curve(T, m0, 10.0, dm - 1.0, 0.001, 2001);
  TEST_EQUAL(O[1000] >= 1, true) // b3 + dm == b5 is observed ...
  TEST_EQUAL(O[1000] - E[1000], 0) // ... and expected, so the Signal is 0 there
}
END_SECTION

START_SECTION(neutral_loss_levels_are_flagged_not_dropped)
{
  std::vector<Refine::Level> levels(2);
  levels[0].shift = 508.13; levels[0].signal = 5;  levels[0].lo = 508.07; levels[0].hi = 508.16;
  levels[1].shift = 526.14; levels[1].signal = 23; levels[1].lo = 526.05; levels[1].hi = 526.20;
  Refine::flagNeutralLosses(levels);
  TEST_EQUAL(levels.size(), 2)            // flagged, not dropped
  TEST_EQUAL(levels[0].neutral_loss, true) // 526.14 - 18.0106 (water) lies in 508.07..508.16
  TEST_EQUAL(levels[1].neutral_loss, false)

  // no stronger level than this one, but the reference 0 is always there: ammonia below 0
  std::vector<Refine::Level> alone(1);
  alone[0].shift = -17.03; alone[0].signal = 6; alone[0].lo = -17.07; alone[0].hi = -16.97;
  Refine::flagNeutralLosses(alone);
  TEST_EQUAL(alone[0].neutral_loss, true)
}
END_SECTION

START_SECTION(no_level_gives_cut_zero)
{
  // the masses are the ladder itself: everything is explained, so Observed == Expected and S == 0
  const Refine::Ladder lad = Refine::ladder("PEPTIDEK", {"b", "y"});
  std::vector<double> masses = concat(lad.prefix, lad.suffix);
  std::sort(masses.begin(), masses.end());
  const Refine::Detect d = Refine::detect(masses, lad, RefineSettings());
  TEST_EQUAL(d.observed_at_zero, 14)
  TEST_EQUAL(d.cut, 0)
  TEST_EQUAL(d.levels.size(), 0)
  TEST_EQUAL(d.fdr_by_height.size(), 61)
  TEST_EQUAL(std::count(d.fdr_by_height.begin() + 1, d.fdr_by_height.end(), -1.0), 60)
}
END_SECTION

START_SECTION(etd_scan_detects_on_c_and_z)
{
  const Refine::Ladder lad = Refine::ladder("PEPTIDEK", {"c", "z"});
  TEST_EQUAL(lad.prefix.size(), 7)
  TEST_EQUAL(lad.suffix.size(), 7)
  TEST_EQUAL(lad.residues.size(), 8)
  TOLERANCE_RELATIVE(1.0 + 1e-9) // absolute 1e-4 Da governs, not the default 1e-5 relative (8 mDa here)
  TOLERANCE_ABSOLUTE(1e-4)
  TEST_REAL_SIMILAR(lad.prefix[0], 97.052764 + 17.026549)   // c1 = P + NH3
  TEST_REAL_SIMILAR(lad.suffix[0], 812.391599 + 0.984016)   // z7 = EPTLDEK + O - N - H (OpenMS z = y - NH3)
  TEST_REAL_SIMILAR(lad.sequence_mass, 927.454928)          // sum of residues + water

  std::vector<double> masses;
  for (double c : lad.prefix) masses.push_back(c + 42.0106);
  const Refine::Detect d = Refine::detect(masses, lad, RefineSettings());
  TEST_EQUAL(d.observed_at_zero, 0)
  TEST_EQUAL(d.cut >= 3, true)
  const auto it = std::find_if(d.levels.begin(), d.levels.end(), [](const Refine::Level& l) { return std::fabs(l.shift - 42.01) < 0.006; });
  TEST_EQUAL(it != d.levels.end(), true)
  if (it != d.levels.end())
  {
    TEST_EQUAL(it->signal, 7)
    TEST_EQUAL(it->n_prefix, 7)
    TEST_EQUAL(it->n_suffix, 0)
    TEST_EQUAL(it->neutral_loss, false)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
