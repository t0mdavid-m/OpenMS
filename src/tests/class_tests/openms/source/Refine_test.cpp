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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/CONCEPT/Constants.h>

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

  // refine_expected_cytc.tsv (make_refine_tables.py --cytc, at the engine's defaults): `quantity\tvalue`, '#' lines and the
  // header skipped. An absent file gives an empty map, so every lookup below fails loudly rather than passing.
  const std::string FI_REFINE_EXPECTED_CYTC = "../../FlashIDA/test-data/refine/refine_expected_cytc.tsv";
  std::map<std::string, std::string> loadQuantities(const std::string& path)
  {
    std::map<std::string, std::string> out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty() || line[0] == '#') continue;
      const size_t tab = line.find('\t');
      if (tab == std::string::npos) continue;
      const std::string key = line.substr(0, tab);
      if (key == "quantity") continue; // the header row
      out[key] = line.substr(tab + 1);
    }
    return out;
  }

  std::vector<std::string> words(const std::string& x)
  {
    std::vector<std::string> out;
    std::stringstream ss(x);
    std::string w;
    while (ss >> w) out.push_back(w);
    return out;
  }

  std::string exclusionName(Refine::Level::Exclusion e)
  {
    switch (e)
    {
      case Refine::Level::Exclusion::None: return "none";
      case Refine::Level::Exclusion::NeutralLoss: return "neutral_loss";
      case Refine::Level::Exclusion::NeutronSatellite: return "neutron_satellite";
      case Refine::Level::Exclusion::Capped: return "capped";
    }
    return "?";
  }

  Refine::Level makeLevel(double shift, int signal, double half_width = 0.03)
  {
    Refine::Level l;
    l.shift = shift;
    l.signal = signal;
    l.lo = shift - half_width;
    l.hi = shift + half_width;
    return l;
  }

  std::string label(const std::optional<Refine::Path::Ion>& ion)
  {
    return ion ? std::string(1, ion->type) + std::to_string(ion->index) : std::string("-");
  }

  bool isPrefixIon(char t)
  {
    return t == 'a' || t == 'b' || t == 'c';
  }

  // The corpus (criterion a): one deconvolved-mass fixture per proteoform under corpus/, and one expectation row per key in
  // refine_expected_corpus.tsv, both written by FlashIDA/test-scripts/refine/make_refine_tables.py --corpus.
  const std::string FI_REFINE_CORPUS = "../../FlashIDA/test-data/refine/corpus/";
  const std::string FI_REFINE_EXPECTED_CORPUS = "../../FlashIDA/test-data/refine/refine_expected_corpus.tsv";

  struct ExpectedShift
  {
    double mass = 0;
    int lo = 0, hi = 0; ///< residues (1-based) the jump's window must hold; 0 = no site
  };
  struct CorpusRow
  {
    std::string backbone, activation;
    double tolerance_ppm = 0, expected_total = 0, precursor_mass_observed = 0, proto_total = 0;
    int proto_score = 0, proto_jumps = 0; ///< the prototype's S and the number of its best path's jumps (proto_shifts)
    std::vector<ExpectedShift> shifts;
  };

  // refine_expected_corpus.tsv: '#' lines skipped, columns by header name; expected_shifts is `mass@lo-hi;...`, `mass@-` for a
  // shift without a site, `-` for none. An absent file gives an empty map (the sections then fail on .at()).
  std::map<std::string, CorpusRow> loadExpectations(const std::string& path)
  {
    std::map<std::string, CorpusRow> out;
    std::ifstream in(path);
    std::string line;
    std::map<std::string, size_t> col;
    auto split = [](const std::string& x, char sep) {
      std::vector<std::string> f;
      std::stringstream ss(x);
      std::string c;
      while (std::getline(ss, c, sep)) f.push_back(c);
      return f;
    };
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty() || line[0] == '#') continue;
      const std::vector<std::string> f = split(line, '\t');
      if (col.empty())
      {
        for (size_t i = 0; i < f.size(); ++i) col[f[i]] = i;
        continue;
      }
      auto at = [&](const char* name) { return f.at(col.at(name)); };
      CorpusRow r;
      r.backbone = at("backbone");
      r.activation = at("activation");
      r.tolerance_ppm = std::stod(at("tolerance_ppm"));
      r.expected_total = std::stod(at("expected_total"));
      r.precursor_mass_observed = std::stod(at("precursor_mass_observed"));
      r.proto_total = std::stod(at("proto_total"));
      r.proto_score = std::stoi(at("proto_score"));
      r.proto_jumps = at("proto_shifts") == "-" ? 0 : (int)split(at("proto_shifts"), ';').size();
      if (at("expected_shifts") != "-")
      {
        for (const std::string& s : split(at("expected_shifts"), ';'))
        {
          const size_t a = s.find('@');
          ExpectedShift e;
          e.mass = std::stod(s.substr(0, a));
          const std::string site = s.substr(a + 1);
          if (site != "-")
          {
            const size_t d = site.find('-');
            e.lo = std::stoi(site.substr(0, d));
            e.hi = std::stoi(site.substr(d + 1));
          }
          r.shifts.push_back(e);
        }
      }
      out[at("key")] = r;
    }
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
  TEST_EQUAL(exclusionName(levels[0].excluded), std::string("neutral_loss")) // 526.14 - 18.0106 (water) lies in 508.07..508.16
  TEST_EQUAL(exclusionName(levels[1].excluded), std::string("none"))

  // no stronger level than this one, but the reference 0 is always there: ammonia below 0
  std::vector<Refine::Level> alone(1);
  alone[0].shift = -17.03; alone[0].signal = 6; alone[0].lo = -17.07; alone[0].hi = -16.97;
  Refine::flagNeutralLosses(alone);
  TEST_EQUAL(exclusionName(alone[0].excluded), std::string("neutral_loss"))
}
END_SECTION

START_SECTION(neutron_satellites_are_excluded_not_merged)
{
  // +655.14 and +654.13 sit 2 and 3 neutrons (1.003355 Da) below +657.17; +42.01 is near nothing
  std::vector<Refine::Level> levels = {makeLevel(42.01, 7), makeLevel(654.13, 11), makeLevel(655.14, 15), makeLevel(657.17, 16)};
  Refine::flagNeutralLosses(levels);
  Refine::flagNeutronSatellites(levels);
  TEST_EQUAL(levels.size(), 4) // marked, not merged into +657.17 and not dropped
  TEST_EQUAL(exclusionName(levels[0].excluded), std::string("none"))
  TEST_EQUAL(exclusionName(levels[1].excluded), std::string("neutron_satellite"))
  TEST_EQUAL(exclusionName(levels[2].excluded), std::string("neutron_satellite"))
  TEST_EQUAL(exclusionName(levels[3].excluded), std::string("none"))
  TEST_REAL_SIMILAR_ABS(levels[2].shift, 655.14, 1e-12) // the shift is not moved onto the reference

  // the references are the levels of higher signal, excluded or not: +507.13 is one neutron below +508.13, itself a water
  // loss of +526.14, and nothing else explains it
  std::vector<Refine::Level> chain = {makeLevel(507.13, 4), makeLevel(508.13, 10), makeLevel(526.14, 23)};
  Refine::flagNeutralLosses(chain);
  Refine::flagNeutronSatellites(chain);
  TEST_EQUAL(exclusionName(chain[0].excluded), std::string("neutron_satellite"))
  TEST_EQUAL(exclusionName(chain[1].excluded), std::string("neutral_loss"))
  TEST_EQUAL(exclusionName(chain[2].excluded), std::string("none"))

  // a satellite is only ever beside a STRONGER level (or 0): the weaker of two neutron-spaced levels is the satellite
  std::vector<Refine::Level> pair = {makeLevel(300.00, 9), makeLevel(301.00, 3)};
  Refine::flagNeutronSatellites(pair);
  TEST_EQUAL(exclusionName(pair[0].excluded), std::string("none"))
  TEST_EQUAL(exclusionName(pair[1].excluded), std::string("neutron_satellite"))
  // and 0 is always a reference: +2.01 is two neutrons above the unshifted ladder
  std::vector<Refine::Level> near0 = {makeLevel(2.01, 30)};
  Refine::flagNeutronSatellites(near0);
  TEST_EQUAL(exclusionName(near0[0].excluded), std::string("neutron_satellite"))
}
END_SECTION

START_SECTION(a_neutral_loss_that_is_also_a_satellite_is_a_neutral_loss)
{
  // the two rules apart: +525.14 is one neutron below +526.14 (a satellite only), +508.13 a water loss of it (a loss only)
  std::vector<Refine::Level> apart = {makeLevel(508.13, 5), makeLevel(525.14, 6), makeLevel(526.14, 23)};
  Refine::flagNeutralLosses(apart);
  Refine::flagNeutronSatellites(apart);
  TEST_EQUAL(exclusionName(apart[0].excluded), std::string("neutral_loss"))
  TEST_EQUAL(exclusionName(apart[1].excluded), std::string("neutron_satellite"))
  TEST_EQUAL(exclusionName(apart[2].excluded), std::string("none"))

  // both at once: the region 508.10..508.20 holds 526.14 - 18.0106 = 508.1294 (water) AND 509.20 - 1.003355 = 508.1966
  // (one neutron below a stronger level) -- the first rule that matches decides, so it is a neutral loss
  std::vector<Refine::Level> both = {makeLevel(508.15, 5, 0.05), makeLevel(509.20, 10), makeLevel(526.14, 23)};
  Refine::flagNeutralLosses(both);
  Refine::flagNeutronSatellites(both);
  TEST_EQUAL(exclusionName(both[0].excluded), std::string("neutral_loss"))
  TEST_EQUAL(exclusionName(both[1].excluded), std::string("none"))
  TEST_EQUAL(exclusionName(both[2].excluded), std::string("none"))
  // ... and on its own the neutron rule would have taken it (the overlap is real, not a miss of the satellite rule)
  std::vector<Refine::Level> sat_only = {makeLevel(508.15, 5, 0.05), makeLevel(509.20, 10)};
  Refine::flagNeutronSatellites(sat_only);
  TEST_EQUAL(exclusionName(sat_only[0].excluded), std::string("neutron_satellite"))
}
END_SECTION

START_SECTION(cap_keeps_the_strongest_levels_and_marks_the_rest)
{
  // 25 levels with distinct signals 1..25 plus one neutral loss stronger than all: the cap keeps the 20 strongest of the
  // 25, marks the 5 weakest Capped, and neither counts nor touches the neutral loss
  std::vector<Refine::Level> levels;
  for (int i = 0; i < 25; ++i) levels.push_back(makeLevel(100.0 + 10.0 * i, i + 1));
  levels.push_back(makeLevel(-18.01, 40));
  levels.back().excluded = Refine::Level::Exclusion::NeutralLoss;
  const int kept = Refine::capLevels(levels, std::vector<double>(levels.size(), 0.0), 20);
  TEST_EQUAL(kept, 20)
  int none = 0, capped = 0;
  for (size_t i = 0; i < 25; ++i)
  {
    none += levels[i].excluded == Refine::Level::Exclusion::None;
    capped += levels[i].excluded == Refine::Level::Exclusion::Capped;
    TEST_EQUAL(exclusionName(levels[i].excluded), std::string(levels[i].signal <= 5 ? "capped" : "none"))
  }
  TEST_EQUAL(none, 20)
  TEST_EQUAL(capped, 5)
  TEST_EQUAL(levels.size(), 26) // marked, not dropped
  TEST_EQUAL(exclusionName(levels[25].excluded), std::string("neutral_loss"))

  // fewer levels than the cap: unchanged
  std::vector<Refine::Level> few;
  for (int i = 0; i < 15; ++i) few.push_back(makeLevel(100.0 + 10.0 * i, i + 1));
  TEST_EQUAL(Refine::capLevels(few, std::vector<double>(few.size(), 0.0), 20), 15)
  for (const auto& l : few) TEST_EQUAL(exclusionName(l.excluded), std::string("none"))

  // a tie on signal goes to the level whose matched masses sum higher; a tie on both to the smaller |shift|
  std::vector<Refine::Level> tie = {makeLevel(200.0, 5), makeLevel(300.0, 5), makeLevel(400.0, 9)};
  TEST_EQUAL(Refine::capLevels(tie, {100.0, 200.0, 1.0}, 2), 2)
  TEST_EQUAL(exclusionName(tie[0].excluded), std::string("capped"))
  TEST_EQUAL(exclusionName(tie[1].excluded), std::string("none"))
  TEST_EQUAL(exclusionName(tie[2].excluded), std::string("none"))
  std::vector<Refine::Level> tie2 = {makeLevel(-300.0, 5), makeLevel(200.0, 5)};
  TEST_EQUAL(Refine::capLevels(tie2, {7.0, 7.0}, 1), 1)
  TEST_EQUAL(exclusionName(tie2[0].excluded), std::string("capped"))
  TEST_EQUAL(exclusionName(tie2[1].excluded), std::string("none"))

  // end to end through detect, where the tie-break is the summed intensity of the masses each level's ions match: +42.01 on
  // b1, b3, b5 (intensity 1 each) and +80.00 on b2, b4, b6 (intensity 100 each), both of signal 3 (both shifts on the grid)
  const Refine::Ladder lad = Refine::ladder("GASPVTLNDKEF", {"b", "y"});
  std::vector<std::pair<double, double>> mi = {{lad.prefix[0] + 42.01, 1.0},   {lad.prefix[2] + 42.01, 1.0},
                                               {lad.prefix[4] + 42.01, 1.0},   {lad.prefix[1] + 80.00, 100.0},
                                               {lad.prefix[3] + 80.00, 100.0}, {lad.prefix[5] + 80.00, 100.0}};
  std::sort(mi.begin(), mi.end());
  std::vector<double> masses, intensities;
  for (const auto& p : mi)
  {
    masses.push_back(p.first);
    intensities.push_back(p.second);
  }
  RefineSettings s;
  const Refine::Detect all = Refine::detect(masses, intensities, lad, s);
  TEST_EQUAL(all.cut, 3)
  TEST_EQUAL(all.levels.size(), 2)
  TEST_EQUAL(all.n_kept, 2)
  s.max_levels = 1;
  const Refine::Detect one = Refine::detect(masses, intensities, lad, s);
  TEST_EQUAL(one.levels.size(), 2)
  TEST_EQUAL(one.n_kept, 1)
  for (const auto& l : one.levels)
  {
    TEST_EQUAL(l.signal, 3)
    TEST_EQUAL(exclusionName(l.excluded), std::string(std::fabs(l.shift - 80.00) < 0.006 ? "none" : "capped"))
  }
}
END_SECTION

START_SECTION(no_level_gives_cut_zero)
{
  // the masses are the ladder itself: everything is explained, so Observed == Expected and S == 0
  const Refine::Ladder lad = Refine::ladder("PEPTIDEK", {"b", "y"});
  std::vector<double> masses = concat(lad.prefix, lad.suffix);
  std::sort(masses.begin(), masses.end());
  const Refine::Detect d = Refine::detect(masses, std::vector<double>(masses.size(), 1.0), lad, RefineSettings());
  TEST_EQUAL(d.observed_at_zero, 14)
  TEST_EQUAL(d.cut, 0)
  TEST_EQUAL(d.levels.size(), 0)
  TEST_EQUAL(d.n_kept, 0)
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
  const Refine::Detect d = Refine::detect(masses, std::vector<double>(masses.size(), 1.0), lad, RefineSettings());
  TEST_EQUAL(d.observed_at_zero, 0)
  TEST_EQUAL(d.cut >= 3, true)
  const auto it = std::find_if(d.levels.begin(), d.levels.end(), [](const Refine::Level& l) { return std::fabs(l.shift - 42.01) < 0.006; });
  TEST_EQUAL(it != d.levels.end(), true)
  if (it != d.levels.end())
  {
    TEST_EQUAL(it->signal, 7)
    TEST_EQUAL(it->n_prefix, 7)
    TEST_EQUAL(it->n_suffix, 0)
    TEST_EQUAL(exclusionName(it->excluded), std::string("none"))
  }
}
END_SECTION

// ---- the report's cytochrome c scan at the engine's DEFAULTS: Detect, Assemble, the targets, the recount ----
// Every expected value comes from refine_expected_cytc.tsv, written by make_refine_tables.py --cytc: the vendored prototype with
// the 2026-10-12 amendments (alpha 0.2, neutral-loss + neutron-satellite exclusions, cap 20, penalty 2) layered over it. These
// are not the report's settings (alpha 0.5, penalty 1), so the numbers are not the report's either.

const MassList fx = loadMassList(FI_REFINE_CYTC);
const std::map<std::string, std::string> cytc_expected = loadQuantities(FI_REFINE_EXPECTED_CYTC);
RefineSettings cytc_settings;
cytc_settings.tolerance_ppm = 10; // and nothing else: every other setting is the engine default
const Refine::Ladder cytc_ladder = Refine::ladder(CYTC105, {"b", "y"});
const Refine::Detect cytc_detect = Refine::detect(fx.masses, fx.intensities, cytc_ladder, cytc_settings);
const Refine::Assemble cytc_assemble = Refine::assemble(fx.masses, fx.intensities, cytc_ladder, cytc_detect, cytc_settings);

// the table's value for a key ("" when absent, which every caller turns into a failure), as words / an int / a double
auto expected = [&](const std::string& key) {
  const auto it = cytc_expected.find(key);
  if (it == cytc_expected.end()) std::cout << "    refine_expected_cytc.tsv has no '" << key << "'\n";
  return it == cytc_expected.end() ? std::string() : it->second;
};
auto expectedInt = [&](const std::string& key) {
  const std::string v = expected(key);
  return v.empty() ? -999999 : std::stoi(v);
};
// "total 526.14 S 61 shifts 2 N 65" -> {total, S, shifts, N}
struct BestRow
{
  double total = 0;
  int score = -1, jumps = -1, ions = -1;
};
auto expectedBest = [&](const std::string& key) {
  const std::vector<std::string> w = words(expected(key));
  BestRow r;
  if (w.size() == 8)
  {
    r.total = std::stod(w[1]);
    r.score = std::stoi(w[3]);
    r.jumps = std::stoi(w[5]);
    r.ions = std::stoi(w[7]);
  }
  return r;
};

START_SECTION(cytc_detect_matches_the_defaults_table)
{
  TEST_EQUAL(cytc_expected.empty(), false)
  TEST_EQUAL(fx.masses.size(), 150)
  TEST_EQUAL(cytc_detect.observed_at_zero, expectedInt("explained_unshifted"))
  TEST_EQUAL(cytc_detect.cut, expectedInt("cut")) // 5 at alpha 0.2 (FDR 0.44 at 4, 0.18 at 5)
  int fdr_rows = 0;
  for (int h = 1; h <= 60; ++h)
  {
    const auto it = cytc_expected.find("fdr_h" + std::to_string(h));
    if (it == cytc_expected.end()) continue;
    ++fdr_rows;
    // the decoys are another shuffle stream (std::mt19937_64 against numpy's), so the FDR agrees only within the noise of K 300
    TEST_REAL_SIMILAR_ABS(cytc_detect.fdr_by_height[h], std::stod(it->second), 0.05)
  }
  TEST_EQUAL(fdr_rows >= 2, true)
  const int n_levels = expectedInt("levels");
  TEST_EQUAL(static_cast<int>(cytc_detect.levels.size()), n_levels)
  TEST_EQUAL(cytc_detect.n_kept, expectedInt("levels_kept"))
  int none = 0, nl = 0, sat = 0, capped = 0, n_b = 0, n_y = 0;
  for (const auto& l : cytc_detect.levels)
  {
    none += l.excluded == Refine::Level::Exclusion::None;
    nl += l.excluded == Refine::Level::Exclusion::NeutralLoss;
    sat += l.excluded == Refine::Level::Exclusion::NeutronSatellite;
    capped += l.excluded == Refine::Level::Exclusion::Capped;
    if (l.excluded != Refine::Level::Exclusion::None) continue;
    n_b += l.n_prefix > 0;
    n_y += l.n_suffix > 0;
  }
  TEST_EQUAL(none, expectedInt("levels_kept"))
  TEST_EQUAL(nl, expectedInt("excluded_neutral_loss"))
  TEST_EQUAL(sat, expectedInt("excluded_neutron_satellite"))
  TEST_EQUAL(capped, expectedInt("excluded_capped"))
  TEST_EQUAL(n_b, expectedInt("lambda_b"))
  TEST_EQUAL(n_y, expectedInt("lambda_y"))
  // every level, ascending shift: shift, signal, typed counts, exclusion
  for (int i = 0; i < n_levels && i < static_cast<int>(cytc_detect.levels.size()); ++i)
  {
    const std::vector<std::string> w = words(expected("level" + std::to_string(i + 1)));
    TEST_EQUAL(w.size(), 5)
    if (w.size() != 5) continue;
    const Refine::Level& l = cytc_detect.levels[i];
    TEST_REAL_SIMILAR_ABS(l.shift, std::stod(w[0]), 0.006)
    TEST_EQUAL(l.signal, std::stoi(w[1]))
    TEST_EQUAL(l.n_prefix, std::stoi(w[2]))
    TEST_EQUAL(l.n_suffix, std::stoi(w[3]))
    TEST_EQUAL(exclusionName(l.excluded), w[4])
  }
}
END_SECTION

START_SECTION(cytc_assemble_matches_the_defaults_table)
{
  const auto& c = cytc_assemble.candidates;
  TEST_EQUAL(static_cast<int>(c.size()), expectedInt("totals")) // one best path per lambda_k
  TEST_EQUAL(c.size() >= 6, true)
  if (c.size() >= 6)
  {
    // the best two, exactly; S = ions - jump_penalty x jumps
    for (int i = 0; i < 2; ++i)
    {
      const std::string k = "best" + std::to_string(i + 1);
      const BestRow e = expectedBest(k);
      TEST_REAL_SIMILAR_ABS(c[i].total, e.total, 0.01)
      TEST_EQUAL(c[i].score, e.score)
      TEST_EQUAL(c[i].jumps, e.jumps)
      TEST_EQUAL(c[i].ions, e.ions)
      TEST_EQUAL(c[i].score, c[i].ions - cytc_settings.jump_penalty * c[i].jumps)
      const std::vector<std::string> lv = words(expected(k + "_levels"));
      TEST_EQUAL(c[i].levels.size(), lv.size())
      for (size_t j = 0; j < lv.size() && j < c[i].levels.size(); ++j) TEST_REAL_SIMILAR_ABS(c[i].levels[j], std::stod(lv[j]), 0.01)
      const std::vector<std::string> wd = words(expected(k + "_windows"));
      TEST_EQUAL(c[i].windows.size(), wd.size())
      for (size_t j = 0; j < wd.size() && j < c[i].windows.size(); ++j)
      {
        const size_t d = wd[j].find('-');
        TEST_EQUAL(c[i].windows[j].first, std::stoi(wd[j].substr(0, d)))
        TEST_EQUAL(c[i].windows[j].second, std::stoi(wd[j].substr(d + 1)))
      }
      const std::vector<std::string> by = words(expected(k + "_ions_b_y"));
      int nb = 0, ny = 0;
      for (const auto& ion : c[i].ions_explained) (isPrefixIon(ion.type) ? nb : ny)++;
      TEST_EQUAL(by.size(), 2)
      if (by.size() == 2)
      {
        TEST_EQUAL(nb, std::stoi(by[0]))
        TEST_EQUAL(ny, std::stoi(by[1]))
      }
      TEST_EQUAL(static_cast<int>(c[i].ions_explained.size()), c[i].ions)
    }
    // third..fifth tie on S and jumps (the table says how many totals share that pair): the same set, ordered by intensity_sum
    const BestRow b3 = expectedBest("best3");
    const int ties = expectedInt("ties_with_best3");
    TEST_EQUAL(ties, 3) // the table's tie block is exactly best3..best5; a wider tie would need another assertion shape
    std::multiset<long> want, got;
    for (int i = 2; i < 5; ++i)
    {
      const BestRow e = expectedBest("best" + std::to_string(i + 1));
      TEST_EQUAL(e.score, b3.score)
      TEST_EQUAL(e.jumps, b3.jumps)
      want.insert(std::lround(e.total * 100));
      got.insert(std::lround(c[i].total * 100));
      TEST_EQUAL(c[i].score, e.score)
      TEST_EQUAL(c[i].jumps, e.jumps)
      TEST_EQUAL(c[i].ions, e.ions)
    }
    TEST_EQUAL(got == want, true)
    TEST_EQUAL(c[2].intensity_sum >= c[3].intensity_sum && c[3].intensity_sum >= c[4].intensity_sum, true)
    TEST_EQUAL(c[5].score < b3.score || (c[5].score == b3.score && c[5].jumps > b3.jumps), true) // nothing else ties the three
  }
}
END_SECTION

const std::vector<Refine::Discriminating> cytc_targets =
  Refine::discriminating(cytc_assemble.candidates[0], cytc_assemble.candidates[1], fx.masses, fx.intensities);

START_SECTION(cytc_discriminating_masses_match_the_defaults_table)
{
  const auto& t = cytc_targets;
  const int n = expectedInt("targets");
  TEST_EQUAL(static_cast<int>(t.size()), n)
  int only_a = 0, only_b = 0, both = 0;
  for (const auto& x : t)
  {
    only_a += x.ion_a && !x.ion_b;
    only_b += !x.ion_a && x.ion_b;
    both += x.ion_a && x.ion_b;
  }
  TEST_EQUAL(only_a, expectedInt("targets_only_a"))
  TEST_EQUAL(only_b, expectedInt("targets_only_b"))
  TEST_EQUAL(both, expectedInt("targets_both"))
  // the whole list, most intense first: mass, the best path's ion, the second's
  for (int i = 0; i < n && i < static_cast<int>(t.size()); ++i)
  {
    const std::vector<std::string> w = words(expected("target" + std::to_string(i + 1)));
    TEST_EQUAL(w.size(), 3)
    if (w.size() != 3) continue;
    TEST_REAL_SIMILAR_ABS(t[i].mass, std::stod(w[0]), 0.01)
    TEST_EQUAL(label(t[i].ion_a), w[1])
    TEST_EQUAL(label(t[i].ion_b), w[2])
  }
}
END_SECTION

START_SECTION(cytc_recount_after_july_ms3_matches_the_defaults_table)
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
  const Refine::Recount r = Refine::recount(a, b, cytc_targets, confirmed, cytc_settings.jump_penalty);
  // S after MS3 = shared + confirmed - 2 x jumps, the best path (2 jumps) against the second (3 jumps)
  TEST_EQUAL(r.score_a, expectedInt("recount_a"))
  TEST_EQUAL(r.score_b, expectedInt("recount_b"))
  TEST_EQUAL(r.confirmed_a, expectedInt("recount_confirmed_a"))
  TEST_EQUAL(r.confirmed_b, 0)
  TEST_EQUAL(r.score_a - r.confirmed_a + cytc_settings.jump_penalty * a.jumps, expectedInt("recount_shared"))
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
  RefineSettings s;
  s.jump_penalty = 1; // the case is built at penalty 1: at 2, total 0 (4 ions, no jump) would tie both on S and win on jumps
  const Refine::Assemble a = Refine::assemble(masses, intensities, lad, d, s);
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

START_SECTION(piece_test_calibrates_on_shared_pieces_only)
{
  // Both paths explain the target as b6 of PEPTIDEK with +42.0106 Da, path a on residue 1, path b on residue 3. So b1+d and
  // b2+d are path a's discriminating, shift-carrying pieces, while b3+d and b4+d are shared. The MS3 masses read 30 ppm high.
  const String seq = "PEPTIDEK";
  const Refine::Ladder lad = Refine::ladder(seq, {"b", "y"});
  const double delta = 42.010565, off = 1.0 + 30e-6;
  Refine::Path pa, pb;
  pa.levels = {0.0, delta};
  pa.windows = {{1, 1}};
  pa.jumps = 1;
  pb.levels = {0.0, delta};
  pb.windows = {{3, 3}};
  pb.jumps = 1;
  Refine::Discriminating target;
  target.mass_idx = 0;
  target.mass = lad.prefix[5] + delta;
  target.intensity = 1.0;
  target.ion_a = Refine::Path::Ion{'b', 6, 0, delta, lad.prefix[5] + delta};
  target.ion_b = target.ion_a;
  auto spectrum = [&](const std::vector<double>& theo) {
    DeconvolvedSpectrum s(0);
    for (double t : theo)
    {
      PeakGroup pg;
      pg.setMonoisotopicMass(t * off);
      s.push_back(pg);
    }
    return s;
  };
  // no shared piece seen: nothing calibrates, b2+d sits 30 ppm out of the 10 ppm window
  const Refine::PieceVerdict alone = Refine::pieceTest(spectrum({lad.prefix[1] + delta}), seq, pa, pb, target, 500.0, 10.0);
  TEST_EQUAL(alone.shared_seen, 0)
  TEST_REAL_SIMILAR_ABS(alone.correction, 1.0, 1e-12)
  TEST_EQUAL(alone.fits_a, false)
  TEST_EQUAL(alone.fits_b, false)
  // two shared pieces seen: their +30 ppm calibrates the copy and pulls b2+d into the window
  const Refine::PieceVerdict cal =
    Refine::pieceTest(spectrum({lad.prefix[1] + delta, lad.prefix[2] + delta, lad.prefix[3] + delta}), seq, pa, pb, target, 500.0, 10.0);
  TEST_EQUAL(cal.shared_seen, 2) // one entry per coincident piece, not one per path
  TEST_REAL_SIMILAR_ABS(cal.correction, 1.0 / off, 1e-9)
  TEST_EQUAL(cal.fits_a, true)
  TEST_EQUAL(cal.fits_b, false)
}
END_SECTION

START_SECTION(piece_test_ignores_an_unannotatable_side)
{
  // piece_test_calibrates_on_shared_pieces_only with ONE change: the target is c6, not b6 (an ETD winner's ion). A c
  // precursor's MS3 pieces cannot be annotated (ADR-0054), so neither side predicts anything: no shared piece is seen,
  // nothing calibrates, nothing fits -- on the very spectrum that gives b6 two shared pieces and a fit for path a.
  //
  // FAILS UNDER: no annotatability guard in pieceTest's predict -- getMS3IonTypes('c') yields the same a/b pieces of
  // the same subsequence as 'b', so this returns shared_seen 2 and fits_a true, exactly the b6 section's verdict.
  const String seq = "PEPTIDEK";
  const Refine::Ladder lad = Refine::ladder(seq, {"b", "y"});
  const double delta = 42.010565, off = 1.0 + 30e-6;
  Refine::Path pa, pb;
  pa.levels = {0.0, delta};
  pa.windows = {{1, 1}};
  pa.jumps = 1;
  pb.levels = {0.0, delta};
  pb.windows = {{3, 3}};
  pb.jumps = 1;
  Refine::Discriminating target;
  target.mass_idx = 0;
  target.mass = lad.prefix[5] + delta;
  target.intensity = 1.0;
  target.ion_a = Refine::Path::Ion{'c', 6, 0, delta, lad.prefix[5] + delta};
  target.ion_b = target.ion_a;
  DeconvolvedSpectrum s(0);
  for (double t : {lad.prefix[1] + delta, lad.prefix[2] + delta, lad.prefix[3] + delta})
  {
    PeakGroup pg;
    pg.setMonoisotopicMass(t * off);
    s.push_back(pg);
  }
  const Refine::PieceVerdict v = Refine::pieceTest(s, seq, pa, pb, target, 500.0, 10.0);
  TEST_EQUAL(v.shared_seen, 0)
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

START_SECTION(jump_penalty_two_rejects_a_two_ion_excursion)
{
  // A 10-residue sequence, b ions only. Unshifted: b1, b2 and b8. The excursion 0 -> +42.0106 (sites 3, 4) -> +100.0 (sites 5
  // on): b3+d, b4+d, b5+l, b6+l. Total 0 stays at 0: b1, b2, b8 = 3 ions, 0 jumps. Total +100 through +42.01: b1, b2, b3+d,
  // b4+d, b5+l, b6+l = 6 ions, 2 jumps -- the excursion buys exactly 3 ions (it loses b8). The alternatives: total +42.01 has
  // 4 ions with 1 jump, total +100 straight from 0 has 4 ions with 1 jump.
  const Refine::Ladder lad = Refine::ladder("GASPVTLNDK", {"b", "y"});
  const double d = 42.0106, l = 100.0;
  std::vector<double> masses = {lad.prefix[0],     lad.prefix[1],     lad.prefix[2] + d, lad.prefix[3] + d,
                                lad.prefix[4] + l, lad.prefix[5] + l, lad.prefix[7]};
  std::sort(masses.begin(), masses.end());
  const std::vector<double> intensities(masses.size(), 1.0);
  Refine::Detect det;
  det.cut = 2;
  det.levels = {makeLevel(d, 2), makeLevel(l, 2)};
  det.levels[0].n_prefix = 2;
  det.levels[1].n_prefix = 2;
  det.n_kept = 2;
  auto at = [](const Refine::Assemble& a, double total) {
    return *std::find_if(a.candidates.begin(), a.candidates.end(), [total](const Refine::Path& p) { return std::fabs(p.total - total) < 1e-6; });
  };

  RefineSettings s1;
  s1.jump_penalty = 1;
  const Refine::Assemble a1 = Refine::assemble(masses, intensities, lad, det, s1);
  TEST_EQUAL(a1.candidates.size(), 3) // totals 0, +42.01, +100
  TEST_REAL_SIMILAR_ABS(a1.candidates[0].total, l, 1e-6)
  TEST_EQUAL(a1.candidates[0].jumps, 2)
  TEST_EQUAL(a1.candidates[0].ions, 6)
  TEST_EQUAL(a1.candidates[0].score, 4) // 6 - 1 x 2, against 3 (total 0) and 4 - 1 = 3 (total +42.01)
  TEST_EQUAL(a1.candidates[0].levels.size(), 3)
  TEST_EQUAL(at(a1, 0.0).score, 3)
  TEST_EQUAL(at(a1, d).score, 3)

  RefineSettings s2; // the default, 2
  TEST_EQUAL(s2.jump_penalty, 2)
  const Refine::Assemble a2 = Refine::assemble(masses, intensities, lad, det, s2);
  TEST_EQUAL(a2.candidates.size(), 3)
  TEST_REAL_SIMILAR_ABS(a2.candidates[0].total, 0.0, 1e-12)
  TEST_EQUAL(a2.candidates[0].jumps, 0)
  TEST_EQUAL(a2.candidates[0].ions, 3)
  TEST_EQUAL(a2.candidates[0].score, 3)
  // total +100 now prefers the straight path (4 ions, 1 jump: 4 - 2 = 2) over the excursion (6 - 4 = 2) on the tie -- fewer jumps
  const Refine::Path p100 = at(a2, l);
  TEST_EQUAL(p100.jumps, 1)
  TEST_EQUAL(p100.ions, 4)
  TEST_EQUAL(p100.score, 2)
  for (const auto& p : a2.candidates) TEST_EQUAL(p.score, p.ions - 2 * p.jumps)
}
END_SECTION

START_SECTION(recount_charges_the_penalty)
{
  // path a explains masses 0..10, path b masses 0..9 and 11: shared 10; the two targets are mass 10 (a only) and mass 11
  // (b only); MS3 confirmed a's. a has 2 jumps, b 3.
  Refine::Path a, b;
  a.jumps = 2;
  b.jumps = 3;
  for (int m = 0; m <= 10; ++m) a.ions_explained.push_back(Refine::Path::Ion{'b', m + 1, m, 0.0, 1000.0 + m});
  for (int m = 0; m <= 9; ++m) b.ions_explained.push_back(Refine::Path::Ion{'b', m + 1, m, 0.0, 1000.0 + m});
  b.ions_explained.push_back(Refine::Path::Ion{'y', 3, 11, 0.0, 1011.0});
  std::vector<Refine::Discriminating> targets(2);
  targets[0].mass_idx = 10;
  targets[0].ion_a = a.ions_explained.back();
  targets[1].mass_idx = 11;
  targets[1].ion_b = b.ions_explained.back();
  const std::vector<std::pair<bool, bool>> confirmed = {{true, false}, {false, false}};
  const Refine::Recount r2 = Refine::recount(a, b, targets, confirmed, 2);
  TEST_EQUAL(r2.confirmed_a, 1)
  TEST_EQUAL(r2.confirmed_b, 0)
  TEST_EQUAL(r2.score_a, 7) // 10 + 1 - 2 x 2
  TEST_EQUAL(r2.score_b, 4) // 10 + 0 - 2 x 3
  TEST_EQUAL(r2.a_wins, true)
  const Refine::Recount r1 = Refine::recount(a, b, targets, confirmed, 1);
  TEST_EQUAL(r1.score_a, 9) // 10 + 1 - 2
  TEST_EQUAL(r1.score_b, 7) // 10 + 0 - 3
}
END_SECTION

// ---- the corpus: seven characterized proteoforms, criterion a (the best total and every sited shift; the prototype's S and jumps).
// The other seven rows of the corpus carry intact disulfides, outside the step-function model (ADR-0049 decision 11). ----

const std::map<std::string, CorpusRow> corpus_rows = loadExpectations(FI_REFINE_EXPECTED_CORPUS);

// One proteoform: Detect + Assemble on its fixture at the run's MS2 tolerance (every other setting the default), the ion
// alphabet of its activation. The best path's total lies within 2 tol of the expected total, every expected shift with a site
// is a jump within tol whose window holds the site (a shift without a site counts through the total), and S and the number of
// jumps equal the prototype's exactly. tol = tolerance_ppm at the observed precursor mass.
auto corpusCriterionA = [&](const std::string& key) {
  const CorpusRow& e = corpus_rows.at(key);
  const MassList m = loadMassList(FI_REFINE_CORPUS + key + ".tsv");
  TEST_EQUAL(m.masses.empty(), false)
  TEST_EQUAL(std::is_sorted(m.masses.begin(), m.masses.end()), true)
  RefineSettings s;
  s.tolerance_ppm = e.tolerance_ppm;
  const Refine::Ladder lad = Refine::ladder(String(e.backbone), FragmentAnalysis::getIonTypesForFragmentationMethod(String(e.activation)));
  TEST_EQUAL(lad.residues.size(), e.backbone.size())
  const Refine::Detect d = Refine::detect(m.masses, m.intensities, lad, s);
  const Refine::Assemble a = Refine::assemble(m.masses, m.intensities, lad, d, s);
  TEST_EQUAL(a.candidates.empty(), false)
  if (a.candidates.empty()) return;
  const Refine::Path& best = a.candidates[0];
  const double tol = e.tolerance_ppm * 1e-6 * e.precursor_mass_observed;
  TEST_REAL_SIMILAR_ABS(best.total, e.expected_total, 2 * tol)
  for (const ExpectedShift& x : e.shifts)
  {
    if (x.lo == 0) continue; // no site: the total carries it
    bool found = false;
    for (size_t j = 0; j + 1 < best.levels.size() && j < best.windows.size(); ++j)
    {
      found = found || (std::fabs(best.levels[j + 1] - best.levels[j] - x.mass) <= tol && best.windows[j].first <= x.lo &&
                        best.windows[j].second >= x.hi);
    }
    if (!found) std::cout << "    " << key << ": no jump " << x.mass << " with a window holding " << x.lo << "-" << x.hi << "\n";
    TEST_EQUAL(found, true)
  }
  TEST_EQUAL(best.score, e.proto_score)
  TEST_EQUAL(best.jumps, e.proto_jumps)
};

START_SECTION(corpus_ca_meets_criterion_a)
{
  corpusCriterionA("ca");
}
END_SECTION

START_SECTION(corpus_cytc_meets_criterion_a)
{
  corpusCriterionA("cytc");
}
END_SECTION

START_SECTION(corpus_myoglobin_meets_criterion_a)
{
  corpusCriterionA("myoglobin");
}
END_SECTION

START_SECTION(corpus_tia_red_meets_criterion_a)
{
  corpusCriterionA("tia_red");
}
END_SECTION

START_SECTION(corpus_ubiquitin_meets_criterion_a)
{
  corpusCriterionA("ubiquitin");
}
END_SECTION

START_SECTION(corpus_lactalbumin_red_meets_criterion_a)
{
  corpusCriterionA("lactalbumin_red");
}
END_SECTION

START_SECTION(corpus_betacasein_meets_criterion_a)
{
  corpusCriterionA("betacasein");
}
END_SECTION

// ---- characterization.method switches FragmentAnalysis::getTopFragmentMatches (ADR-0049) ----
// Driven through FragmentAnalysis itself: the refine section reproduces the defaults table's best path through the contract,
// the flashtnt section shows the same call under the default method never takes the Refine path.

// A loadable config: mode off, the cytochrome c sequence, the named method; @p refine_block is spliced after the sequence.
auto methodConfigJson = [](const std::string& method, const std::string& refine_block) {
  return std::string(R"({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": { "rank_by": "qscore", "max_precursors": 5 },
    "characterization": { "mode": "off", "method": ")") + method + R"(", "protein_sequence": ")" + std::string(CYTC105) + "\"" +
         refine_block + R"( },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 }
    }
  })";
};

// One PeakGroup per fixture mass: its mono mass and intensity, and a representative charge of mass / 1000 (at least 1 -- the
// fixture's third column counts charges, it names none). One LogMzPeak at that charge carries the intensity and the m/z
// window the out-arrays read; qscore 0.9 so the fragments' qscore column is seen to come from the PeakGroup.
auto spectrumFromMassList = [](const MassList& m, int scan, double precursor_mono) {
  DeconvolvedSpectrum ds(scan);
  for (size_t i = 0; i < m.masses.size(); ++i)
  {
    const int z = std::max(1, static_cast<int>(std::lround(m.masses[i] / 1000.0)));
    PeakGroup pg(z, z, true);
    FLASHHelperClasses::LogMzPeak lp;
    lp.mz = m.masses[i] / z + Constants::PROTON_MASS_U;
    lp.intensity = static_cast<float>(m.intensities[i]);
    lp.abs_charge = z;
    lp.isotopeIndex = 0;
    pg.push_back(lp);
    pg.updateMonoMassAndIsotopeIntensities(0.0); // the group's intensity = the peak's
    pg.setMonoisotopicMass(m.masses[i]);         // and the mass exactly as the fixture writes it
    pg.setRepAbsCharge(z);
    pg.setQscore(0.9);
    ds.push_back(pg);
  }
  PeakGroup precursor;
  precursor.setMonoisotopicMass(precursor_mono);
  ds.setPrecursorPeakGroup(precursor);
  return ds;
};

START_SECTION(getTopFragmentMatches_under_refine_fills_the_contract)
{
  const Config cfg(methodConfigJson("refine", R"(, "refine": {})"));
  TEST_EQUAL((int)cfg.characterization().method, (int)AnnotationMethod::Refine)
  FragmentAnalysis fa(cfg);
  DeconvolvedSpectrum ms2 = spectrumFromMassList(fx, 4242, 12351.30);
  TEST_EQUAL(ms2.size(), 150)
  double masses[200], qscores[200], wstarts[200], wends[200];
  int charges[200], frag_indices[200];
  char ion_types[200];
  FragmentAnalysis::ProteoformMatch r;
  r.tag_count = 7; // anything but -1: the path must write it
  const int found = fa.getTopFragmentMatches(CYTC105, 200, masses, qscores, charges, wstarts, wends, ion_types, frag_indices, ms2, r, "HCD",
                                             0.0, nullptr);

  const BestRow b1 = expectedBest("best1"); // total 526.14 S 61 shifts 2 N 65
  TEST_EQUAL(b1.ions > 0, true)
  TEST_EQUAL(r.score, b1.score)
  TEST_EQUAL(r.tag_count, -1)
  TEST_EQUAL(r.region_start, 0)
  TEST_EQUAL(r.region_end, static_cast<int>(CYTC105.size()))
  TEST_EQUAL(r.proteoform_sequence, std::string(CYTC105))
  TEST_EQUAL(r.matched_protein, "input_protein")
  TEST_EQUAL(static_cast<int>(r.fragments.size()), b1.ions)
  TEST_EQUAL(r.total_match_count, b1.ions)
  TEST_EQUAL(found, b1.ions)

  // ptm_sites = the best path's windows, each with its level step: [1-8, -89.03], [14-22, 526.14 - (-89.03)]
  const std::vector<std::string> lv = words(expected("best1_levels"));
  const std::vector<std::string> wd = words(expected("best1_windows"));
  TEST_EQUAL(r.ptm_sites.size(), wd.size())
  TEST_EQUAL(lv.size(), wd.size() + 1)
  for (size_t j = 0; j < wd.size() && j < r.ptm_sites.size() && j + 1 < lv.size(); ++j)
  {
    const size_t d = wd[j].find('-');
    TEST_EQUAL(r.ptm_sites[j].start_position, std::stoi(wd[j].substr(0, d)))
    TEST_EQUAL(r.ptm_sites[j].end_position, std::stoi(wd[j].substr(d + 1)))
    TEST_REAL_SIMILAR_ABS(r.ptm_sites[j].mass_shift, std::stod(lv[j + 1]) - std::stod(lv[j]), 0.02)
  }

  // every fragment: the explaining mass within tolerance of its shifted ladder ion, the qscore read off its PeakGroup
  int nb = 0, ny = 0;
  for (const auto& f : r.fragments)
  {
    TEST_EQUAL(std::fabs(f.diff_ppm) <= 10.0 + 1e-6, true)
    TEST_REAL_SIMILAR_ABS(f.diff_da, f.observed_mass - f.theoretical_mass, 1e-9)
    TEST_REAL_SIMILAR_ABS(f.qscore, 0.9, 1e-9)
    TEST_EQUAL(f.ion_index >= 1 && f.ion_index < static_cast<int>(CYTC105.size()), true)
    (f.ion_type == "b" ? nb : ny)++;
  }
  const std::vector<std::string> by = words(expected("best1_ions_b_y"));
  TEST_EQUAL(by.size(), 2)
  if (by.size() == 2)
  {
    TEST_EQUAL(nb, std::stoi(by[0]))
    TEST_EQUAL(ny, std::stoi(by[1]))
  }
  // the out-arrays come from the same ions: every reported mass is an explained fragment's, at the group's charge
  for (int i = 0; i < found; ++i)
  {
    const bool known = std::any_of(r.fragments.begin(), r.fragments.end(), [&](const auto& f) {
      return f.observed_mass == masses[i] && f.ion_index == frag_indices[i] && f.ion_type == std::string(1, ion_types[i]);
    });
    TEST_EQUAL(known, true)
    TEST_EQUAL(charges[i], std::max(1, static_cast<int>(std::lround(masses[i] / 1000.0))))
    TEST_EQUAL(wstarts[i] < wends[i], true)
  }

  // the candidate record: every lambda_k, the cut, the kept levels, the precursor mass difference
  TEST_EQUAL(r.refine != nullptr, true)
  if (r.refine)
  {
    TEST_EQUAL(static_cast<int>(r.refine->assemble.candidates.size()), expectedInt("totals"))
    TEST_EQUAL(r.refine->cut, expectedInt("cut"))
    TEST_EQUAL(r.refine->n_levels, expectedInt("levels_kept"))
    TEST_EQUAL(static_cast<int>(r.refine->detect.levels.size()), expectedInt("levels"))
    TEST_REAL_SIMILAR_ABS(r.refine->precursor_mass_difference, 12351.30 - (cytc_ladder.sequence_mass + b1.total), 0.05)
  }
}
END_SECTION

START_SECTION(getTopFragmentMatches_under_flashtnt_is_unchanged)
{
  // FI_MS2_CYTC ("Spec scan=N\tRT" then "mz\tintensity"), deconvolved by the engine: the existing path's input.
  std::vector<double> mzs, ints;
  double rt = 0;
  {
    std::ifstream in("../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt");
    std::string line;
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      const size_t tab = line.find('\t');
      if (tab == std::string::npos) continue;
      if (line.rfind("Spec", 0) == 0)
      {
        rt = std::stod(line.substr(tab + 1));
        continue;
      }
      mzs.push_back(std::stod(line.substr(0, tab)));
      ints.push_back(std::stod(line.substr(tab + 1)));
    }
  }
  TEST_EQUAL(mzs.empty(), false)

  double masses[50], qscores[50], wstarts[50], wends[50];
  int charges[50], frag_indices[50];
  char ion_types[50];
  for (const std::string& method : {std::string("flashtnt"), std::string("refine")})
  {
    const Config cfg(methodConfigJson(method, method == "refine" ? R"(, "refine": {})" : ""));
    Deconvolution deconv(cfg, {10.0, 10.0, 10.0});
    TEST_EQUAL(deconv.deconvolveMSn(mzs.data(), ints.data(), static_cast<int>(mzs.size()), rt, 0.0, 0) > 0, true)
    FragmentAnalysis fa(cfg);
    FragmentAnalysis::ProteoformMatch r;
    r.tag_count = -7; // a sentinel neither path writes
    fa.getTopFragmentMatches(CYTC105, 50, masses, qscores, charges, wstarts, wends, ion_types, frag_indices, deconv.storedMS2(), r);
    if (method == "flashtnt")
    {
      TEST_EQUAL(r.refine == nullptr, true) // the tag path never builds a candidate record ...
      TEST_EQUAL(r.tag_count >= 0, true)    // ... and records the tagger's tag count
    }
    else
    {
      TEST_EQUAL(r.refine != nullptr, true) // the same call under refine takes the other branch
      TEST_EQUAL(r.tag_count, -1)
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
