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

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;

// |a - b| <= tol, absolute only (TEST_REAL_SIMILAR passes on relative OR absolute agreement); prints both on failure.
#define TEST_REAL_SIMILAR_ABS(a, b, tol)                                                                          \
  {                                                                                                               \
    const double va_ = (a), vb_ = (b);                                                                            \
    if (!(std::fabs(va_ - vb_) <= (tol))) std::cout << "    |" << va_ << " - " << vb_ << "| > " << (tol) << "\n"; \
    TEST_EQUAL(std::fabs(va_ - vb_) <= (tol), true)                                                               \
  }

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

  // The report's cytochrome c scan (150 masses) and its supplied sequence: Met 1 plus the 104 residues of horse cytochrome c.
  const std::string FI_REFINE_CYTC = "../../FlashIDA/test-data/refine/cytc_grp1_prec1.tsv";
  const std::string FI_REFINE_JULY = "../../FlashIDA/test-data/refine/july_ms3_pieces.tsv";
  const String CYTC105 = "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";

  // july_ms3_pieces.tsv: `target_ion\tmass\tintensity`, '#' lines and the header skipped; one spectrum per target ion.
  std::map<std::string, DeconvolvedSpectrum> loadJulyPieces(const std::string& path)
  {
    std::map<std::string, DeconvolvedSpectrum> out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty() || line[0] == '#') continue;
      std::stringstream ss(line);
      std::string ion, m, i;
      std::getline(ss, ion, '\t');
      std::getline(ss, m, '\t');
      std::getline(ss, i, '\t');
      char* end = nullptr;
      const double mass = std::strtod(m.c_str(), &end);
      if (end == m.c_str()) continue; // the header row
      auto it = out.emplace(ion, DeconvolvedSpectrum(0)).first;
      PeakGroup pg; // a mono mass is all the piece test reads (makeSyntheticPeakGroup's shape, FLASHIda_exploration_test)
      pg.setMonoisotopicMass(mass);
      it->second.push_back(pg);
    }
    return out;
  }

  std::string label(const std::optional<Refine::Path::Ion>& ion)
  {
    return ion ? std::string(1, ion->type) + std::to_string(ion->index) : std::string("-");
  }

  bool isPrefixIon(char t)
  {
    return t == 'a' || t == 'b' || t == 'c';
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

// ---- the report's cytochrome c scan: Detect, Assemble (column_dp_v85.py), the targets (targets_v88.py), the recount ----

const MassList fx = loadMassList(FI_REFINE_CYTC);
RefineSettings cytc_settings;
cytc_settings.tolerance_ppm = 10; // defaults otherwise: fdr 0.5, K 300, grid 0.01, +-700, min 3
const Refine::Ladder cytc_ladder = Refine::ladder(CYTC105, {"b", "y"});
const Refine::Detect cytc_detect = Refine::detect(fx.masses, cytc_ladder, cytc_settings);
const Refine::Assemble cytc_assemble = Refine::assemble(fx.masses, fx.intensities, cytc_ladder, cytc_detect, cytc_settings);

START_SECTION(cytc_detect_reproduces_the_report)
{
  TEST_EQUAL(fx.masses.size(), 150)
  TEST_EQUAL(cytc_detect.observed_at_zero, 36)
  TEST_EQUAL(cytc_detect.cut, 4) // 73.7 % at 3, 44.4 % at 4 (K = 10,000); K = 300 gives the same cut
  TEST_REAL_SIMILAR_ABS(cytc_detect.fdr_by_height[3], 0.737, 0.05)
  TEST_REAL_SIMILAR_ABS(cytc_detect.fdr_by_height[4], 0.444, 0.05)
  TEST_EQUAL(cytc_detect.levels.size(), 26)
  int kept = 0, n_b = 0, n_y = 0;
  for (const auto& l : cytc_detect.levels)
  {
    if (l.neutral_loss) continue;
    ++kept;
    n_b += l.n_prefix > 0;
    n_y += l.n_suffix > 0;
  }
  TEST_EQUAL(kept, 22) // -18.02, -17.02, +508.12, +637.18 flagged
  TEST_EQUAL(n_b, 11)
  TEST_EQUAL(n_y, 15)
  auto level = [&](double shift) {
    return *std::min_element(cytc_detect.levels.begin(), cytc_detect.levels.end(),
                             [&](const Refine::Level& a, const Refine::Level& b) { return std::fabs(a.shift - shift) < std::fabs(b.shift - shift); });
  };
  const Refine::Level nt = level(-89.03), heme = level(526.14), nh3 = level(-17.02);
  TEST_REAL_SIMILAR_ABS(nt.shift, -89.03, 0.006)
  TEST_EQUAL(nt.signal, 7)
  TEST_EQUAL(nt.n_prefix, 7)
  TEST_EQUAL(nt.n_suffix, 0)
  TEST_REAL_SIMILAR_ABS(heme.shift, 526.14, 0.006)
  TEST_EQUAL(heme.signal, 23)
  TEST_EQUAL(heme.n_prefix, 23)
  TEST_EQUAL(heme.n_suffix, 0)
  TEST_REAL_SIMILAR_ABS(nh3.shift, -17.02, 0.006)
  TEST_EQUAL(nh3.neutral_loss, true)
  TEST_EQUAL(nh3.n_suffix, 6)
  int flagged_lt = 0; // the four flagged levels, by shift
  for (double x : {-18.02, -17.02, 508.12, 637.18}) flagged_lt += level(x).neutral_loss && std::fabs(level(x).shift - x) < 0.006;
  TEST_EQUAL(flagged_lt, 4)
}
END_SECTION

START_SECTION(cytc_assemble_reproduces_the_report)
{
  const auto& c = cytc_assemble.candidates;
  TEST_EQUAL(c.size(), 182) // one best path per lambda_k
  const auto& best = c[0];
  TEST_REAL_SIMILAR_ABS(best.total, 526.14, 0.01)
  TEST_EQUAL(best.score, 63)
  TEST_EQUAL(best.jumps, 2)
  TEST_EQUAL(best.ions, 65)
  TEST_EQUAL(best.levels.size(), 3)
  TEST_REAL_SIMILAR_ABS(best.levels[0], 0.0, 1e-9)
  TEST_REAL_SIMILAR_ABS(best.levels[1], -89.03, 0.01)
  TEST_REAL_SIMILAR_ABS(best.levels[2], 526.14, 0.01)
  TEST_EQUAL(best.windows.size(), 2)
  TEST_EQUAL(best.windows[0].first, 1)
  TEST_EQUAL(best.windows[0].second, 8)
  TEST_EQUAL(best.windows[1].first, 14)
  TEST_EQUAL(best.windows[1].second, 22)
  int nb = 0, ny = 0;
  for (const auto& ion : best.ions_explained) (isPrefixIon(ion.type) ? nb : ny)++;
  TEST_EQUAL(nb, 29)
  TEST_EQUAL(ny, 36)
  TEST_EQUAL(best.ions_explained.size(), 65)

  const auto& second = c[1];
  TEST_REAL_SIMILAR_ABS(second.total, 655.19, 0.01)
  TEST_EQUAL(second.score, 50)
  TEST_EQUAL(second.jumps, 5)
  TEST_EQUAL(second.ions, 55) // 55 ions from 52 masses (targets_v88.py)
  TEST_EQUAL(second.levels.size(), 6)
  const double second_levels[] = {0.0, -89.03, 526.14, 799.56, 171.18, 655.19};
  for (size_t j = 0; j < 6; ++j) TEST_REAL_SIMILAR_ABS(second.levels[j], second_levels[j], 0.01)
  const std::pair<int, int> second_windows[] = {{1, 8}, {14, 22}, {35, 35}, {38, 42}, {44, 45}};
  TEST_EQUAL(second.windows.size(), 5)
  for (size_t j = 0; j < 5 && j < second.windows.size(); ++j)
  {
    TEST_EQUAL(second.windows[j].first, second_windows[j].first)
    TEST_EQUAL(second.windows[j].second, second_windows[j].second)
  }

  // third..fifth: +412.04 / +525.12 / +655.22, all S 49, N 52, 3 mass shifts; their order is intensity_sum, descending
  std::multiset<long> totals;
  for (size_t i = 2; i < 5; ++i)
  {
    totals.insert(std::lround(c[i].total * 100));
    TEST_EQUAL(c[i].score, 49)
    TEST_EQUAL(c[i].ions, 52)
    TEST_EQUAL(c[i].jumps, 3)
  }
  TEST_EQUAL(totals == std::multiset<long>({41204, 52512, 65522}), true)
  TEST_EQUAL(c[2].intensity_sum >= c[3].intensity_sum && c[3].intensity_sum >= c[4].intensity_sum, true)
  TEST_EQUAL(c[5].score < 49 || (c[5].score == 49 && c[5].jumps > 3), true) // nothing else ties the three
}
END_SECTION

const std::vector<Refine::Discriminating> cytc_targets =
  Refine::discriminating(cytc_assemble.candidates[0], cytc_assemble.candidates[1], fx.masses, fx.intensities);

START_SECTION(cytc_discriminating_masses_reproduce_the_report)
{
  const auto& t = cytc_targets;
  TEST_EQUAL(t.size(), 21)
  int only_a = 0, only_b = 0, both = 0;
  for (const auto& x : t)
  {
    only_a += x.ion_a && !x.ion_b;
    only_b += !x.ion_a && x.ion_b;
    both += x.ion_a && x.ion_b;
  }
  TEST_EQUAL(only_a, 15)
  TEST_EQUAL(only_b, 2)
  TEST_EQUAL(both, 4)
  TEST_REAL_SIMILAR_ABS(t[0].mass, 8111.25, 0.01)
  TEST_EQUAL(label(t[0].ion_a), "y71")
  TEST_EQUAL(label(t[0].ion_b), "-")
  // the whole list of targets_v88.txt, most intense first: mass, reported path's ion, +655.19 path's ion
  const struct
  {
    double mass;
    const char* a;
    const char* b;
  } report[] = {{8111.25, "y71", "-"},  {6036.92, "b51", "-"},  {7649.77, "b65", "y68"}, {7435.65, "b63", "b62"}, {5255.60, "b44", "-"},
                {7780.82, "b66", "-"},  {7177.56, "b61", "-"},  {5921.90, "b50", "-"},   {7909.86, "b67", "b66"}, {4297.06, "b35", "-"},
                {5657.79, "b48", "-"},  {8072.92, "b68", "-"},  {9441.70, "b80", "-"},   {5556.74, "b47", "-"},   {8054.22, "y70", "-"},
                {8315.04, "b70", "b69"}, {9960.94, "b85", "-"}, {6107.95, "b52", "-"},   {6222.01, "b53", "-"},   {5538.72, "-", "b46"},
                {7650.78, "-", "y62"}};
  for (size_t i = 0; i < 21 && i < t.size(); ++i)
  {
    TEST_REAL_SIMILAR_ABS(t[i].mass, report[i].mass, 0.01)
    TEST_EQUAL(label(t[i].ion_a), report[i].a)
    TEST_EQUAL(label(t[i].ion_b), report[i].b)
  }
  // masses both paths explain as the same ion are shared, not targets: y39 on 4570.49 and y40 on 4701.52 Da
  for (const auto& x : t)
  {
    TEST_EQUAL(std::fabs(x.mass - 4570.49) > 0.01, true)
    TEST_EQUAL(std::fabs(x.mass - 4701.52) > 0.01, true)
  }
}
END_SECTION

START_SECTION(cytc_recount_after_july_ms3_reproduces_the_report)
{
  // july_ms3_pieces.tsv: for each of b47, b48, b51 the pooled deconvolved masses of its 14 July MS3 scans
  const std::map<std::string, DeconvolvedSpectrum> july = loadJulyPieces(FI_REFINE_JULY);
  TEST_EQUAL(july.size(), 3)
  const auto& a = cytc_assemble.candidates[0];
  const auto& b = cytc_assemble.candidates[1];
  std::vector<std::pair<bool, bool>> confirmed(cytc_targets.size(), {false, false});
  int tested = 0;
  for (size_t i = 0; i < cytc_targets.size(); ++i)
  {
    const auto it = july.find(label(cytc_targets[i].ion_a));
    if (it == july.end()) continue; // every other target: not isolated in the July run, untested
    ++tested;
    const Refine::PieceVerdict v = Refine::pieceTest(it->second, CYTC105, a, b, cytc_targets[i], 500.0, 10.0);
    TEST_EQUAL(v.fits_a, true)
    TEST_EQUAL(v.fits_b, false)
    confirmed[i] = {v.fits_a, v.fits_b};
  }
  TEST_EQUAL(tested, 3)
  const Refine::Recount r = Refine::recount(a, b, cytc_targets, confirmed);
  // shared 46 = 65−19 = 52−6; the 49/45 of the V79 script was against +637.18, excluded as a water loss since V85
  TEST_EQUAL(r.score_a, 47) // 46 + 3 - 2
  TEST_EQUAL(r.score_b, 41) // 46 + 0 - 5
  TEST_EQUAL(r.confirmed_a, 3)
  TEST_EQUAL(r.confirmed_b, 0)
  TEST_EQUAL(r.a_wins, true)
}
END_SECTION

START_SECTION(tie_on_score_and_shifts_resolves_by_intensity)
{
  // PEPTIDEK, Lambda_b = {+10}, Lambda_y = {+20}: the totals are 0, +20, +10, +30 in construction order. Total +10 explains
  // b1, b2, b6+10, b7+10, y2, y1; total +20 explains b1, b2, y7+20, y6+20, y2, y1: S 5 with 1 mass shift each. The +10 masses are
  // the intense ones, so +10 ranks first although +20 was built first.
  const Refine::Ladder lad = Refine::ladder("PEPTIDEK", {"b", "y"});
  std::vector<std::pair<double, double>> mi = {{lad.prefix[0], 1.0},        {lad.prefix[1], 1.0},        {lad.prefix[5] + 10.0, 100.0},
                                               {lad.prefix[6] + 10.0, 100.0}, {lad.suffix[0] + 20.0, 1.0}, {lad.suffix[1] + 20.0, 1.0},
                                               {lad.suffix[5], 1.0},         {lad.suffix[6], 1.0}};
  std::sort(mi.begin(), mi.end());
  std::vector<double> masses, intensities;
  for (const auto& p : mi)
  {
    masses.push_back(p.first);
    intensities.push_back(p.second);
  }
  Refine::Detect d;
  d.cut = 3;
  d.levels.resize(2);
  d.levels[0].shift = 10.0;
  d.levels[0].signal = 3;
  d.levels[0].n_prefix = 2;
  d.levels[1].shift = 20.0;
  d.levels[1].signal = 3;
  d.levels[1].n_suffix = 2;
  const Refine::Assemble a = Refine::assemble(masses, intensities, lad, d, RefineSettings());
  TEST_EQUAL(a.candidates.size(), 4)
  TEST_REAL_SIMILAR_ABS(a.candidates[0].total, 10.0, 1e-9)
  TEST_REAL_SIMILAR_ABS(a.candidates[1].total, 20.0, 1e-9)
  TEST_EQUAL(a.candidates[0].score, 5)
  TEST_EQUAL(a.candidates[1].score, 5)
  TEST_EQUAL(a.candidates[0].jumps, 1)
  TEST_EQUAL(a.candidates[1].jumps, 1)
  TEST_REAL_SIMILAR_ABS(a.candidates[0].intensity_sum, 204.0, 1e-9)
  TEST_REAL_SIMILAR_ABS(a.candidates[1].intensity_sum, 6.0, 1e-9)
}
END_SECTION

START_SECTION(piece_shared_when_masses_coincide)
{
  // Two paths put a mass shift on residue 1 of PEPTIDEK: +42.0106 (path a) and +42.0106 + CO + 1 ppm (path b). Path a's
  // piece b3 and path b's piece a3 (a = b - CO) then sit 1 ppm apart: one piece, labelled two ways. Seen alone it is
  // shared -- it calibrates, but it fits neither path.
  const String seq = "PEPTIDEK";
  const Refine::Ladder lad = Refine::ladder(seq, {"b", "y"});
  const double delta = 42.010565, co = 27.994915;
  const double b3 = lad.prefix[2] + delta;
  Refine::Path pa, pb;
  pa.levels = {0.0, delta};
  pa.windows = {{1, 1}};
  pa.jumps = 1;
  pb.levels = {0.0, delta + co + b3 * 1e-6};
  pb.windows = {{1, 1}};
  pb.jumps = 1;
  Refine::Discriminating target;
  target.mass_idx = 0;
  target.mass = lad.prefix[5] + delta;
  target.intensity = 1.0;
  target.ion_a = Refine::Path::Ion{'b', 6, 0, delta, lad.prefix[5] + delta};
  target.ion_b = Refine::Path::Ion{'b', 6, 0, delta + co, lad.prefix[5] + delta + co};
  DeconvolvedSpectrum ms3(0);
  PeakGroup pg;
  pg.setMonoisotopicMass(b3);
  ms3.push_back(pg);
  const Refine::PieceVerdict v = Refine::pieceTest(ms3, seq, pa, pb, target, 500.0, 10.0);
  TEST_EQUAL(v.shared_seen, 1)
  TEST_REAL_SIMILAR_ABS(v.correction, 1.0, 1e-12)
  TEST_EQUAL(v.fits_a, false)
  TEST_EQUAL(v.fits_b, false)
}
END_SECTION

START_SECTION(no_level_gives_one_candidate)
{
  // Detect found nothing (cut 0, no level): one matrix, lambda_k = 0, the path without a mass shift
  const Refine::Ladder lad = Refine::ladder("PEPTIDEK", {"b", "y"});
  std::vector<double> masses = {lad.prefix[1], lad.prefix[3], lad.suffix[2], lad.suffix[4] + 5.0};
  std::sort(masses.begin(), masses.end());
  const std::vector<double> intensities(masses.size(), 1.0);
  const Refine::Detect d; // cut 0, no levels
  const Refine::Assemble a = Refine::assemble(masses, intensities, lad, d, RefineSettings());
  TEST_EQUAL(a.candidates.size(), 1)
  TEST_REAL_SIMILAR_ABS(a.candidates[0].total, 0.0, 1e-12)
  TEST_EQUAL(a.candidates[0].jumps, 0)
  TEST_EQUAL(a.candidates[0].ions, 3) // b2, b4, y5 (site 3); the +5 Da mass is unexplained
  TEST_EQUAL(a.candidates[0].score, a.candidates[0].ions)
  TEST_EQUAL(a.candidates[0].levels.size(), 1)
  TEST_EQUAL(a.candidates[0].windows.size(), 0)
  TEST_REAL_SIMILAR_ABS(a.candidates[0].intensity_sum, 3.0, 1e-12)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
