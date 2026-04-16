// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <cmath>
#include <string>
#include <vector>

using namespace OpenMS;

START_TEST(MS3FragmentMatcher, "$Id$")

START_SECTION(getMS3IonTypes)
{
  // b-precursor: a, b (same direction) + yb, ya (cross-direction)
  auto b_types = MS3FragmentMatcher::getMS3IonTypes('b');
  TEST_EQUAL(b_types.size(), 4)
  TEST_EQUAL(b_types[0], "a")
  TEST_EQUAL(b_types[1], "b")
  TEST_EQUAL(b_types[2], "yb")
  TEST_EQUAL(b_types[3], "ya")

  // a-precursor: same as b-precursor
  auto a_types = MS3FragmentMatcher::getMS3IonTypes('a');
  TEST_EQUAL(a_types.size(), 4)

  // y-precursor: a, b, y (standard)
  auto y_types = MS3FragmentMatcher::getMS3IonTypes('y');
  TEST_EQUAL(y_types.size(), 3)
  TEST_EQUAL(y_types[0], "a")
  TEST_EQUAL(y_types[1], "b")
  TEST_EQUAL(y_types[2], "y")

  // z-precursor: same as y-precursor
  auto z_types = MS3FragmentMatcher::getMS3IonTypes('z');
  TEST_EQUAL(z_types.size(), 3)

  // Unknown: defaults to y-precursor behavior
  auto unk_types = MS3FragmentMatcher::getMS3IonTypes('?');
  TEST_EQUAL(unk_types.size(), 3)
}
END_SECTION

START_SECTION(isPrefixIonType)
{
  TEST_TRUE(MS3FragmentMatcher::isPrefixIonType("a"))
  TEST_TRUE(MS3FragmentMatcher::isPrefixIonType("b"))
  TEST_TRUE(! MS3FragmentMatcher::isPrefixIonType("y"))
  TEST_TRUE(! MS3FragmentMatcher::isPrefixIonType("yb"))
  TEST_TRUE(! MS3FragmentMatcher::isPrefixIonType("ya"))
}
END_SECTION

START_SECTION(getIonShift)
{
  double co = 27.994915;
  double water = 18.010565;
  TEST_REAL_SIMILAR(MS3FragmentMatcher::getIonShift("a"), -co)
  TEST_REAL_SIMILAR(MS3FragmentMatcher::getIonShift("b"), 0.0)
  TEST_REAL_SIMILAR(MS3FragmentMatcher::getIonShift("y"), water)
  TEST_REAL_SIMILAR(MS3FragmentMatcher::getIonShift("yb"), 0.0)
  TEST_REAL_SIMILAR(MS3FragmentMatcher::getIonShift("ya"), -co)
}
END_SECTION

START_SECTION(computeTheoreticalMasses_no_ptm)
{
  // Use "ACDEF" — 5 residues, expect 4 positions per ion type
  std::string seq = "ACDEF";
  int n = 5;

  // Get residue masses from ResidueDB for verification
  std::vector<double> res(n);
  for (int i = 0; i < n; ++i)
    res[i] = ResidueDB::getInstance()->getResidue(seq[i])->getMonoWeight(Residue::Internal);

  double co = 27.994915;
  double water = 18.010565;

  // Test with just b-ions first
  auto masses = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"b"});
  TEST_EQUAL(masses.size(), 4) // n-1 = 4 positions

  // b-ions: cumulative prefix sums, shift = 0
  double cumul = 0.0;
  for (int i = 0; i < 4; ++i)
  {
    cumul += res[i];
    TEST_REAL_SIMILAR(masses[i].mass, cumul)
    TEST_EQUAL(masses[i].position, i + 1)
    TEST_EQUAL(masses[i].ion_type, "b")
    TEST_TRUE(! masses[i].includes_ptm)
  }

  // a-ions: b - CO
  auto a_masses = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"a"});
  TEST_EQUAL(a_masses.size(), 4)
  for (int i = 0; i < 4; ++i)
    TEST_REAL_SIMILAR(a_masses[i].mass, masses[i].mass - co)

  // y-ions: cumulative suffix sums + water
  auto y_masses = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"y"});
  TEST_EQUAL(y_masses.size(), 4)
  double cumul_suffix = 0.0;
  for (int i = 0; i < 4; ++i)
  {
    cumul_suffix += res[n - 1 - i];
    TEST_REAL_SIMILAR(y_masses[i].mass, cumul_suffix + water)
    TEST_EQUAL(y_masses[i].position, i + 1)
    TEST_EQUAL(y_masses[i].ion_type, "y")
  }

  // yb-ions: same cumulative suffix but no water (shift = 0)
  auto yb_masses = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"yb"});
  TEST_EQUAL(yb_masses.size(), 4)
  for (int i = 0; i < 4; ++i)
    TEST_REAL_SIMILAR(yb_masses[i].mass, y_masses[i].mass - water)

  // ya-ions: cumulative suffix - CO (shift = -CO)
  auto ya_masses = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"ya"});
  TEST_EQUAL(ya_masses.size(), 4)
  for (int i = 0; i < 4; ++i)
    TEST_REAL_SIMILAR(ya_masses[i].mass, y_masses[i].mass - water - co)

  // Multiple ion types: should get 4 * num_types entries
  auto all_masses = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"a", "b", "y"});
  TEST_EQUAL(all_masses.size(), 12)

  // Verify monotonically increasing within each type
  for (int i = 1; i < 4; ++i)
    TEST_TRUE(masses[i].mass > masses[i - 1].mass)
  for (int i = 1; i < 4; ++i)
    TEST_TRUE(y_masses[i].mass > y_masses[i - 1].mass)
}
END_SECTION

START_SECTION(computeTheoreticalMasses_fixed_ptm)
{
  std::string seq = "ACDEF";
  double ptm_shift = 42.0106; // acetylation

  // Fixed PTM at position 3 (1-based) = residue D
  FragmentAnalysis::PTMSite fixed_ptm;
  fixed_ptm.start_position = 3;
  fixed_ptm.end_position = 3;
  fixed_ptm.position = 3;
  fixed_ptm.mass_shift = ptm_shift;

  auto masses_no_ptm = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"b"});
  auto masses_ptm = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"b"}, {fixed_ptm});

  // b1, b2: before PTM position 3 — should be identical
  TEST_REAL_SIMILAR(masses_ptm[0].mass, masses_no_ptm[0].mass)
  TEST_REAL_SIMILAR(masses_ptm[1].mass, masses_no_ptm[1].mass)

  // b3, b4: at or past PTM position — should include PTM mass
  TEST_REAL_SIMILAR(masses_ptm[2].mass, masses_no_ptm[2].mass + ptm_shift)
  TEST_REAL_SIMILAR(masses_ptm[3].mass, masses_no_ptm[3].mass + ptm_shift)

  // Test suffix ions: y-ions
  auto y_no_ptm = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"y"});
  auto y_ptm = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"y"}, {fixed_ptm});

  // y1 (covers pos 5), y2 (covers 4-5): PTM at pos 3 is outside — identical
  TEST_REAL_SIMILAR(y_ptm[0].mass, y_no_ptm[0].mass)
  TEST_REAL_SIMILAR(y_ptm[1].mass, y_no_ptm[1].mass)

  // y3 (covers 3-5), y4 (covers 2-5): PTM at pos 3 is inside — includes shift
  TEST_REAL_SIMILAR(y_ptm[2].mass, y_no_ptm[2].mass + ptm_shift)
  TEST_REAL_SIMILAR(y_ptm[3].mass, y_no_ptm[3].mass + ptm_shift)
}
END_SECTION

START_SECTION(computeTheoreticalMasses_ambiguous_ptm)
{
  std::string seq = "ACDEF";
  double ptm_shift = 79.966; // phosphorylation

  // Ambiguous PTM spanning positions 2-4 (1-based) = residues C, D, E
  FragmentAnalysis::PTMSite amb_ptm;
  amb_ptm.start_position = 2;
  amb_ptm.end_position = 4;
  amb_ptm.position = 3;
  amb_ptm.mass_shift = ptm_shift;

  auto masses = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"b"}, {amb_ptm});

  // b1 (covers pos 1): PTM range [2,4] outside — single entry, no PTM
  TEST_EQUAL(masses[0].ion_type, "b")
  TEST_TRUE(! masses[0].includes_ptm)

  // b2 (covers pos 1-2): PTM range [2,4] partially overlaps — dual entries
  // Find the two b2 entries
  std::vector<MS3FragmentMatcher::TheoreticalMass> b2_entries;
  for (const auto& m : masses)
    if (m.position == 2) b2_entries.push_back(m);
  TEST_EQUAL(b2_entries.size(), 2)
  // One with PTM, one without
  bool found_with = false, found_without = false;
  for (const auto& e : b2_entries)
  {
    if (e.includes_ptm) found_with = true;
    else found_without = true;
  }
  TEST_TRUE(found_with)
  TEST_TRUE(found_without)
  // Mass difference should be ptm_shift
  double diff = std::abs(b2_entries[0].mass - b2_entries[1].mass);
  TEST_REAL_SIMILAR(diff, ptm_shift)

  // b4 (covers pos 1-4): PTM range [2,4] fully covered — single entry with PTM included
  std::vector<MS3FragmentMatcher::TheoreticalMass> b4_entries;
  for (const auto& m : masses)
    if (m.position == 4) b4_entries.push_back(m);
  TEST_EQUAL(b4_entries.size(), 1)
}
END_SECTION

START_SECTION(matchSpectrum)
{
  // Build theoretical masses from a known sequence
  std::string seq = "ACDEF";
  auto theoretical = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"b", "y"});
  // 4 b-ions + 4 y-ions = 8 theoretical masses

  // Create a synthetic spectrum with 3 exact matches (b1, b2, y1) + 1 non-match
  std::vector<double> obs_masses;
  obs_masses.push_back(theoretical[0].mass);  // b1 — exact match
  obs_masses.push_back(theoretical[1].mass);  // b2 — exact match
  obs_masses.push_back(theoretical[4].mass);  // y1 — exact match
  obs_masses.push_back(99999.0);              // garbage — no match

  DeconvolvedSpectrum spec(0);
  for (double m : obs_masses)
  {
    PeakGroup pg(1, 1, true);
    pg.setMonoisotopicMass(m);
    spec.push_back(pg);
  }

  int count = MS3FragmentMatcher::matchSpectrum(spec, theoretical, 10.0);
  TEST_EQUAL(count, 3)

  // Test with ppm errors
  std::vector<double> ppm_errors;
  int count2 = MS3FragmentMatcher::matchSpectrum(spec, theoretical, 10.0, &ppm_errors);
  TEST_EQUAL(count2, 3)
  TEST_EQUAL(ppm_errors.size(), 3)
  // Exact matches should have ~0 ppm error
  for (double e : ppm_errors)
    TEST_TRUE(std::abs(e) < 0.01)

  // Test with a systematic shift (50 ppm)
  double shift_factor = 1.0 + 50.0e-6;
  DeconvolvedSpectrum shifted_spec(0);
  for (double m : obs_masses)
  {
    PeakGroup pg(1, 1, true);
    pg.setMonoisotopicMass(m * shift_factor);
    shifted_spec.push_back(pg);
  }

  // At 10 ppm tolerance, shifted masses should NOT match
  int count_tight = MS3FragmentMatcher::matchSpectrum(shifted_spec, theoretical, 10.0);
  TEST_EQUAL(count_tight, 0)

  // At 100 ppm tolerance, shifted masses should match
  std::vector<double> shift_errors;
  int count_loose = MS3FragmentMatcher::matchSpectrum(shifted_spec, theoretical, 100.0, &shift_errors);
  TEST_EQUAL(count_loose, 3)
  // ppm errors should be ~50
  for (double e : shift_errors)
    TEST_REAL_SIMILAR(e, 50.0)

  // Empty inputs
  DeconvolvedSpectrum empty_spec(0);
  TEST_EQUAL(MS3FragmentMatcher::matchSpectrum(empty_spec, theoretical, 10.0), 0)
  TEST_EQUAL(MS3FragmentMatcher::matchSpectrum(spec, {}, 10.0), 0)
}
END_SECTION

START_SECTION(extractSubsequence)
{
  std::string protein = "ABCDEFGHIJ"; // 10 residues

  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 2;  // 0-based: proteoform is residues 2-7 = "CDEFGH"
  ctx.region_end = 8;

  // b3: first 3 residues of proteoform = "CDE"
  std::string b3 = MS3FragmentMatcher::extractSubsequence(protein, ctx, 'b', 3);
  TEST_EQUAL(b3, "CDE")

  // y3: last 3 residues of proteoform = "FGH"
  std::string y3 = MS3FragmentMatcher::extractSubsequence(protein, ctx, 'y', 3);
  TEST_EQUAL(y3, "FGH")

  // Full proteoform as b6
  std::string b6 = MS3FragmentMatcher::extractSubsequence(protein, ctx, 'b', 6);
  TEST_EQUAL(b6, "CDEFGH")

  // Invalid: index too large
  std::string bad = MS3FragmentMatcher::extractSubsequence(protein, ctx, 'b', 7);
  TEST_EQUAL(bad, "")

  // a-precursor treated same as b-precursor
  std::string a3 = MS3FragmentMatcher::extractSubsequence(protein, ctx, 'a', 3);
  TEST_EQUAL(a3, "CDE")
}
END_SECTION

START_SECTION(rebasePTMSites)
{
  // Proteoform has 10 residues, PTM at positions 3-5 (1-based, ambiguous)
  FragmentAnalysis::PTMSite ptm;
  ptm.start_position = 3;
  ptm.end_position = 5;
  ptm.position = 4;
  ptm.mass_shift = 80.0;

  // Fixed PTM at position 7
  FragmentAnalysis::PTMSite fixed;
  fixed.start_position = 7;
  fixed.end_position = 7;
  fixed.position = 7;
  fixed.mass_shift = 42.0;

  std::vector<FragmentAnalysis::PTMSite> ptms = {ptm, fixed};

  // b5 subsequence: positions 1-5 in proteoform
  auto rebased = MS3FragmentMatcher::rebasePTMSites(ptms, 1, 5);
  TEST_EQUAL(rebased.size(), 1) // only the ambiguous PTM overlaps; fixed at 7 is outside
  TEST_EQUAL(rebased[0].start_position, 3) // 3 - 1 + 1 = 3
  TEST_EQUAL(rebased[0].end_position, 5)   // 5 - 1 + 1 = 5
  TEST_REAL_SIMILAR(rebased[0].mass_shift, 80.0)

  // y5 subsequence: positions 6-10 in proteoform
  auto rebased_y = MS3FragmentMatcher::rebasePTMSites(ptms, 6, 5);
  TEST_EQUAL(rebased_y.size(), 1) // only fixed PTM at 7 overlaps
  TEST_EQUAL(rebased_y[0].start_position, 2) // 7 - 6 + 1 = 2
  TEST_EQUAL(rebased_y[0].end_position, 2)
  TEST_REAL_SIMILAR(rebased_y[0].mass_shift, 42.0)

  // Partial overlap: subsequence positions 4-8
  auto rebased_mid = MS3FragmentMatcher::rebasePTMSites(ptms, 4, 5);
  TEST_EQUAL(rebased_mid.size(), 2) // both PTMs overlap
  // Ambiguous PTM [3,5] clipped to [4,5], rebased to [1,2]
  TEST_EQUAL(rebased_mid[0].start_position, 1)
  TEST_EQUAL(rebased_mid[0].end_position, 2)
  // Fixed PTM [7,7] rebased to [4,4]
  TEST_EQUAL(rebased_mid[1].start_position, 4)
  TEST_EQUAL(rebased_mid[1].end_position, 4)
}
END_SECTION

START_SECTION(calibrateAndScore)
{
  // Build a known sequence and theoretical masses
  std::string protein = "ACDEFGHIKLMNPQRSTVWY"; // 20 residues
  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 0;
  ctx.region_end = 20;

  // Pretend we're fragmenting b10: first 10 residues
  char frag_type = 'b';
  int frag_index = 10;

  std::string subseq = MS3FragmentMatcher::extractSubsequence(protein, ctx, frag_type, frag_index);
  TEST_EQUAL(subseq, "ACDEFGHIKL")

  auto ion_types = MS3FragmentMatcher::getMS3IonTypes(frag_type);
  auto theoretical = MS3FragmentMatcher::computeTheoreticalMasses(subseq, ion_types);

  // Create two variant spectra:
  // Variant 0: 5 matching masses with +50 ppm shift (simulating systematic instrument drift)
  // Variant 1: 3 matching masses with the same shift
  double shift_factor = 1.0 + 50.0e-6;

  DeconvolvedSpectrum var0(0);
  for (int i = 0; i < 5 && i < static_cast<int>(theoretical.size()); ++i)
  {
    PeakGroup pg(1, 1, true);
    pg.setMonoisotopicMass(theoretical[i].mass * shift_factor);
    var0.push_back(pg);
  }

  DeconvolvedSpectrum var1(0);
  for (int i = 0; i < 3 && i < static_cast<int>(theoretical.size()); ++i)
  {
    PeakGroup pg(1, 1, true);
    pg.setMonoisotopicMass(theoretical[i].mass * shift_factor);
    var1.push_back(pg);
  }

  std::vector<const DeconvolvedSpectrum*> variants = {&var0, &var1};

  auto scores = MS3FragmentMatcher::calibrateAndScore(
    variants, protein, ctx, frag_type, frag_index,
    MS3FragmentMatcher::LOOSE_TOLERANCE_PPM, 10.0);

  TEST_EQUAL(scores.size(), 2)
  // After calibration (corrects the 50 ppm shift), matching at 10 ppm should succeed
  TEST_TRUE(scores[0] >= 4.0)  // variant 0 had 5 masses, should match most
  TEST_TRUE(scores[1] >= 2.0)  // variant 1 had 3 masses
  TEST_TRUE(scores[0] > scores[1]) // variant 0 should score higher

  // Without calibration (tight tolerance only), same masses would NOT match
  // Verify by testing matchSpectrum directly at 10 ppm with unshifted theoretical
  int direct_count = MS3FragmentMatcher::matchSpectrum(var0, theoretical, 10.0);
  TEST_EQUAL(direct_count, 0) // 50 ppm shift exceeds 10 ppm tolerance
}
END_SECTION

START_SECTION(matchSpectrum_match_details)
{
  // Same setup as matchSpectrum test: "ACDEF" with b+y ions
  std::string seq = "ACDEF";
  auto theoretical = MS3FragmentMatcher::computeTheoreticalMasses(seq, {"b", "y"});

  // Synthetic spectrum: 3 exact matches (b1, b2, y1) + 1 garbage
  DeconvolvedSpectrum spec(0);
  PeakGroup pg1(1, 1, true); pg1.setMonoisotopicMass(theoretical[0].mass); spec.push_back(pg1);
  PeakGroup pg2(1, 1, true); pg2.setMonoisotopicMass(theoretical[1].mass); spec.push_back(pg2);
  PeakGroup pg3(1, 1, true); pg3.setMonoisotopicMass(theoretical[4].mass); spec.push_back(pg3);
  PeakGroup pg4(1, 1, true); pg4.setMonoisotopicMass(99999.0); spec.push_back(pg4);

  std::vector<MS3FragmentMatcher::MatchDetail> details;
  int count = MS3FragmentMatcher::matchSpectrum(spec, theoretical, 10.0, nullptr, &details);

  TEST_EQUAL(count, 3)
  TEST_EQUAL(details.size(), 3)

  // Verify each detail has correct fields
  for (const auto& md : details)
  {
    TEST_TRUE(! md.ion_type.empty())
    TEST_TRUE(md.position >= 1)
    TEST_TRUE(std::abs(md.ppm_error) < 0.01) // exact matches
    TEST_REAL_SIMILAR(md.observed_mass, md.theoretical_mass) // exact matches
    TEST_TRUE(! md.includes_ptm) // no PTMs in this sequence
  }

  // Check that we got the expected ion types (b1, b2, y1 — order depends on spectrum iteration)
  int b_count = 0, y_count = 0;
  for (const auto& md : details)
  {
    if (md.ion_type == "b") ++b_count;
    if (md.ion_type == "y") ++y_count;
  }
  TEST_EQUAL(b_count, 2)
  TEST_EQUAL(y_count, 1)

  // Test with 50 ppm shift at 100 ppm tolerance
  double shift_factor = 1.0 + 50.0e-6;
  DeconvolvedSpectrum shifted_spec(0);
  PeakGroup sp1(1, 1, true); sp1.setMonoisotopicMass(theoretical[0].mass * shift_factor); shifted_spec.push_back(sp1);
  PeakGroup sp2(1, 1, true); sp2.setMonoisotopicMass(theoretical[1].mass * shift_factor); shifted_spec.push_back(sp2);

  std::vector<MS3FragmentMatcher::MatchDetail> shift_details;
  int shift_count = MS3FragmentMatcher::matchSpectrum(shifted_spec, theoretical, 100.0, nullptr, &shift_details);
  TEST_EQUAL(shift_count, 2)
  TEST_EQUAL(shift_details.size(), 2)
  for (const auto& md : shift_details)
    TEST_REAL_SIMILAR(md.ppm_error, 50.0)
}
END_SECTION

END_TEST
