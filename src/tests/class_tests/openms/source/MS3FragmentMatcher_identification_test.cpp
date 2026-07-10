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

// Helper to get internal residue mass
static double resMass(char aa)
{
  const Residue* r = ResidueDB::getInstance()->getResidue(aa);
  return r ? r->getMonoWeight(Residue::Internal) : 0.0;
}

START_TEST(MS3FragmentMatcher_identification, "$Id$")

/////////////////////////////////////////////////////////////
// T1: computeProteinPrefixMasses
/////////////////////////////////////////////////////////////

START_SECTION(computeProteinPrefixMasses)
{
  // Verify prefix sums for "ACDE" (4 residues)
  std::string seq = "ACDE";
  auto prefix = MS3FragmentMatcher::computeProteinPrefixMasses(seq);

  // prefix has size n+1 = 5
  TEST_EQUAL(prefix.size(), 5)

  // prefix[0] = 0
  TEST_REAL_SIMILAR(prefix[0], 0.0)

  // prefix[i] = cumulative sum of internal residue masses for [0, i)
  double mA = resMass('A');
  double mC = resMass('C');
  double mD = resMass('D');
  double mE = resMass('E');

  TEST_REAL_SIMILAR(prefix[1], mA)
  TEST_REAL_SIMILAR(prefix[2], mA + mC)
  TEST_REAL_SIMILAR(prefix[3], mA + mC + mD)
  TEST_REAL_SIMILAR(prefix[4], mA + mC + mD + mE)
}
END_SECTION

/////////////////////////////////////////////////////////////
// T2: computeEquivalentIon - same direction b from b-precursor
//     b2 from b4 of proteoform at start=0 -> equiv b2, offset=0
/////////////////////////////////////////////////////////////

START_SECTION(computeEquivalentIon_same_direction_b_from_b)
{
  std::string protein = "ACDE";
  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 0;
  ctx.region_end = 4;

  auto prefix = MS3FragmentMatcher::computeProteinPrefixMasses(protein);

  // b2 theoretical mass in the subsequence "ACDE" = mass(A) + mass(C) + shift(b)=0
  double ms3_theo = resMass('A') + resMass('C');

  std::string equiv_type;
  int equiv_index = 0;
  double mass_offset = 0.0;

  MS3FragmentMatcher::computeEquivalentIon(
    protein, ctx,
    'b', 4,       // precursor: b4
    "b", 2,       // MS3 fragment: b2
    ms3_theo,
    prefix,
    equiv_type, equiv_index, mass_offset);

  // Same direction: equiv type = "b", index = start + ms3_ion_index = 0 + 2 = 2
  TEST_EQUAL(equiv_type, "b")
  TEST_EQUAL(equiv_index, 2)
  // Since start=0, the MS3 b2 is the same as full-protein b2 -> offset = 0
  TEST_REAL_SIMILAR(mass_offset, 0.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// T3: computeEquivalentIon - opposite direction yb from b-precursor
//     yb2 from b4 of "ACDE" -> equiv b2, offset = (mA+mC) - (mD+mE)
/////////////////////////////////////////////////////////////

START_SECTION(computeEquivalentIon_opposite_direction_yb_from_b)
{
  std::string protein = "ACDE";
  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 0;
  ctx.region_end = 4;

  auto prefix = MS3FragmentMatcher::computeProteinPrefixMasses(protein);

  // yb2 of subsequence "ACDE": last 2 residues = D, E; shift(yb) = 0
  double ms3_theo = resMass('D') + resMass('E');

  std::string equiv_type;
  int equiv_index = 0;
  double mass_offset = 0.0;

  MS3FragmentMatcher::computeEquivalentIon(
    protein, ctx,
    'b', 4,       // precursor: b4
    "yb", 2,      // MS3 fragment: yb2
    ms3_theo,
    prefix,
    equiv_type, equiv_index, mass_offset);

  // Opposite direction: equiv_type = "b", equiv_index = start + N - ms3_ion_index = 0 + 4 - 2 = 2
  TEST_EQUAL(equiv_type, "b")
  TEST_EQUAL(equiv_index, 2)

  // theo_equiv = prefix[2] + shift(b)=0 = mA + mC
  // offset = (mA + mC) - (mD + mE)
  double expected_offset = (resMass('A') + resMass('C')) - (resMass('D') + resMass('E'));
  TEST_REAL_SIMILAR(mass_offset, expected_offset)
}
END_SECTION

/////////////////////////////////////////////////////////////
// T4: computeEquivalentIon - same direction y from y-precursor
//     y2 from y3 of protein "ACDEFG" (P=6, start=0, end=6) -> equiv y2, offset=0
/////////////////////////////////////////////////////////////

START_SECTION(computeEquivalentIon_same_direction_y_from_y)
{
  std::string protein = "ACDEFG";
  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 0;
  ctx.region_end = 6;

  auto prefix = MS3FragmentMatcher::computeProteinPrefixMasses(protein);

  // y3 precursor -> subsequence = last 3 residues of proteoform = "EFG"
  // y2 of "EFG" = last 2 residues = F, G; shift(y) = water = 18.010565
  double water = 18.010565;
  double ms3_theo = resMass('F') + resMass('G') + water;

  std::string equiv_type;
  int equiv_index = 0;
  double mass_offset = 0.0;

  MS3FragmentMatcher::computeEquivalentIon(
    protein, ctx,
    'y', 3,       // precursor: y3
    "y", 2,       // MS3 fragment: y2
    ms3_theo,
    prefix,
    equiv_type, equiv_index, mass_offset);

  // Same direction y from y-precursor: equiv_type = "y", equiv_index = P - end + ms3_ion_index = 6 - 6 + 2 = 2
  TEST_EQUAL(equiv_type, "y")
  TEST_EQUAL(equiv_index, 2)
  // Full-protein y2 = last 2 residues (F, G) + water = same as MS3 y2 from y3
  TEST_REAL_SIMILAR(mass_offset, 0.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// T5: computeEquivalentIon - opposite direction b from y-precursor
//     b2 from y4 of "ACDEFG" -> equiv y2, offset = theo_y2 - theo_b2
/////////////////////////////////////////////////////////////

START_SECTION(computeEquivalentIon_opposite_direction_b_from_y)
{
  std::string protein = "ACDEFG";
  int P = 6;
  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 0;
  ctx.region_end = 6;

  auto prefix = MS3FragmentMatcher::computeProteinPrefixMasses(protein);

  // y4 precursor -> subsequence = last 4 residues = "DEFG"
  // "b" is a prefix ion type in y-precursor context -> treated as opposite direction
  // b2 of "DEFG" = D + E; shift(b) = 0
  double ms3_theo = resMass('D') + resMass('E');

  std::string equiv_type;
  int equiv_index = 0;
  double mass_offset = 0.0;

  MS3FragmentMatcher::computeEquivalentIon(
    protein, ctx,
    'y', 4,       // precursor: y4
    "b", 2,       // MS3 fragment: b2
    ms3_theo,
    prefix,
    equiv_type, equiv_index, mass_offset);

  // b is a prefix ion type -> goes to the else branch (opposite for y-precursor)
  // equiv_type = "y", equiv_index = P - end + N - ms3_ion_index = 6 - 6 + 4 - 2 = 2
  TEST_EQUAL(equiv_type, "y")
  TEST_EQUAL(equiv_index, 2)

  // theo_equiv = (prefix[P] - prefix[P - equiv_index]) + water
  //            = (prefix[6] - prefix[4]) + water
  //            = mass(E) + mass(F) + water
  double water = 18.010565;
  double theo_y2 = (prefix[P] - prefix[P - 2]) + water;
  double expected_offset = theo_y2 - ms3_theo;
  TEST_REAL_SIMILAR(mass_offset, expected_offset)
  // Verify offset is nonzero (different ion masses)
  TEST_TRUE(std::abs(mass_offset) > 0.01)
}
END_SECTION

/////////////////////////////////////////////////////////////
// T6: computeEquivalentIon - truncated proteoform
//     Protein "ACDEFGHIJ" (P=9), proteoform [2,7)="DEFGH",
//     b3 precursor, b2 fragment -> equiv b4 (start=2+2), offset > 0
/////////////////////////////////////////////////////////////

START_SECTION(computeEquivalentIon_truncated_proteoform)
{
  std::string protein = "ACDEFGHIJ";
  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 2;  // proteoform starts at residue index 2 = "DEFGH" (indices 2,3,4,5,6)
  ctx.region_end = 7;    // exclusive

  auto prefix = MS3FragmentMatcher::computeProteinPrefixMasses(protein);

  // b3 precursor of proteoform "DEFGH" -> subsequence = "DEF" (first 3 residues)
  // b2 of "DEF" = D + E; shift(b) = 0
  double ms3_theo = resMass('D') + resMass('E');

  std::string equiv_type;
  int equiv_index = 0;
  double mass_offset = 0.0;

  MS3FragmentMatcher::computeEquivalentIon(
    protein, ctx,
    'b', 3,       // precursor: b3
    "b", 2,       // MS3 fragment: b2
    ms3_theo,
    prefix,
    equiv_type, equiv_index, mass_offset);

  // Same direction b from b-precursor: equiv_type = "b", equiv_index = start + ms3_ion_index = 2 + 2 = 4
  TEST_EQUAL(equiv_type, "b")
  TEST_EQUAL(equiv_index, 4)

  // theo_equiv = prefix[4] + shift(b)=0 = mA + mC + mD + mE
  // ms3_theo = mD + mE
  // offset = (mA + mC + mD + mE) - (mD + mE) = mA + mC
  double expected_offset = resMass('A') + resMass('C');
  TEST_REAL_SIMILAR(mass_offset, expected_offset)
  TEST_TRUE(mass_offset > 0.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// T7: calibrateAndScore with detailed results
//     Synthetic spectrum with one peak at theoretical b2 mass
//     for a b4 subsequence, verify MatchResult is populated
//     with correct FragmentMatch data
/////////////////////////////////////////////////////////////

START_SECTION(calibrateAndScore_detailed_results)
{
  std::string protein = "ACDEFGHI";
  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 0;
  ctx.region_end = 8;

  // b4 precursor -> subsequence = "ACDE"
  char frag_type = 'b';
  int frag_index = 4;

  std::string subseq = MS3FragmentMatcher::extractSubsequence(protein, ctx, frag_type, frag_index);
  TEST_EQUAL(subseq, "ACDE")

  // Compute theoretical masses for b-type ion types: a, b, yb, ya
  auto ion_types = MS3FragmentMatcher::getMS3IonTypes(frag_type);
  auto theoretical = MS3FragmentMatcher::computeTheoreticalMasses(subseq, ion_types);

  // Find the b2 theoretical mass: ion_type == "b", position == 2
  double theo_b2 = 0.0;
  for (const auto& tm : theoretical)
  {
    if (tm.ion_type == "b" && tm.position == 2)
    {
      theo_b2 = tm.mass;
      break;
    }
  }
  TEST_TRUE(theo_b2 > 0.0)

  // Create a spectrum with one peak at exactly the b2 mass
  DeconvolvedSpectrum ds(0);
  PeakGroup pg(1, 1, true);
  pg.setMonoisotopicMass(theo_b2);
  ds.push_back(pg);

  std::vector<const DeconvolvedSpectrum*> variants = {&ds};
  std::vector<FragmentAnalysis::ProteoformMatch> detailed;

  auto scores = MS3FragmentMatcher::calibrateAndScore(
    variants, protein, ctx, frag_type, frag_index,
    MS3FragmentMatcher::LOOSE_TOLERANCE_PPM, 10.0,
    &detailed);

  TEST_EQUAL(scores.size(), 1)
  TEST_TRUE(scores[0] >= 1.0)  // at least 1 match

  TEST_EQUAL(detailed.size(), 1)
  TEST_TRUE(detailed[0].fragments.size() >= 1)

  // Find the b2 match in the detailed results
  bool found_b2 = false;
  for (const auto& fm : detailed[0].fragments)
  {
    if (fm.ion_type == "b" && fm.ion_index == 2)
    {
      found_b2 = true;
      // Observed mass should be close to theo_b2 (after calibration correction)
      TEST_REAL_SIMILAR(fm.observed_mass, theo_b2)
      // Equivalent: since start=0, b2 in MS3 maps to b2 in full protein
      TEST_EQUAL(fm.equiv_type, "b")
      TEST_EQUAL(fm.equiv_index, 2)
      // Adjusted mass should be close to observed + offset (offset ~0 since start=0)
      TEST_TRUE(fm.adjusted_mass > 0.0)
      break;
    }
  }
  TEST_TRUE(found_b2)
}
END_SECTION

/////////////////////////////////////////////////////////////
// T8: PPM-HONEST frame projection -- same-direction b-from-b, cross-frame.
//     The MS3->MS2 projection is multiplicative (scale the equivalent theoretical by the
//     sub-fragment's measured fractional error), NOT Da-additive. A sub-fragment measured at a
//     known +10 ppm must report +10 ppm in the MS2 frame -- the OLD Da-additive form deflated it
//     to 10*theo_sub/theo_equiv (~4.9 ppm here), so these assertions fail under the old code.
//     Calibration trap: a LONE injected peak is absorbed by the Pass-1 median correction
//     (correction_factor cancels it) -> diff_ppm ~ 0 under BOTH old and new. Neutralize it with an
//     ODD, majority-0-ppm set so median = 0, correction_factor = 1, and the +10 ppm survives.
/////////////////////////////////////////////////////////////

START_SECTION(calibrateAndScore_ppm_frame_projection_same_dir_b)
{
  // protein indices: A0 C1 D2 E3 F4 G5 H6 I7 K8 L9
  std::string protein = "ACDEFGHIKL";
  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 3;                 // proteoform "EFGHIKL" -> offset > 0 (cross-frame)
  ctx.region_end = 10;

  char frag_type = 'b';
  int frag_index = 4;                   // b4 precursor -> subseq = protein[3:7] = "EFGH"
  std::string subseq = MS3FragmentMatcher::extractSubsequence(protein, ctx, frag_type, frag_index);
  TEST_EQUAL(subseq, "EFGH")

  // Sub-frame b-ion masses (prefix sums of the subseq; computeTheoreticalMasses reproduces these exactly).
  double theo_b1 = resMass('E');
  double theo_b2 = resMass('E') + resMass('F');
  double theo_b3 = resMass('E') + resMass('F') + resMass('G');

  const double kPpm = 10.0;
  DeconvolvedSpectrum ds(0);
  auto addPeak = [&](double m) { PeakGroup pg(1, 1, true); pg.setMonoisotopicMass(m); ds.push_back(pg); };
  addPeak(theo_b1);                          // 0 ppm
  addPeak(theo_b2 * (1.0 + kPpm * 1e-6));    // +10 ppm (target)
  addPeak(theo_b3);                          // 0 ppm

  std::vector<const DeconvolvedSpectrum*> variants = {&ds};
  std::vector<FragmentAnalysis::ProteoformMatch> detailed;
  MS3FragmentMatcher::calibrateAndScore(
    variants, protein, ctx, frag_type, frag_index,
    MS3FragmentMatcher::LOOSE_TOLERANCE_PPM, 50.0, &detailed);

  TEST_EQUAL(detailed.size(), 1)
  // Calibration neutralized: median of {0, 0, +10} = 0 -> correction_factor = 1 (the +10 ppm survives).
  TEST_REAL_SIMILAR(detailed[0].correction_factor, 1.0)

  bool found = false;
  for (const auto& fm : detailed[0].fragments)
  {
    if (fm.ion_type == "b" && fm.ion_index == 2)
    {
      found = true;
      // Same-direction b: equiv b(start + 2) = b5, offset > 0.
      TEST_EQUAL(fm.equiv_type, "b")
      TEST_EQUAL(fm.equiv_index, 5)
      // PPM-honest: the MS2-frame residual EQUALS the sub-frame measured +10 ppm (fails under Da-additive).
      TEST_REAL_SIMILAR(fm.diff_ppm, kPpm)
      // Multiplicative: adjusted = equivalent theoretical * (1 + 10 ppm).
      TEST_REAL_SIMILAR(fm.adjusted_mass, fm.theoretical_mass * (1.0 + kPpm * 1e-6))
      // diff_da is the equiv-frame projected Da residual (theoretical * (ratio - 1)), larger than the raw sub-frame Da.
      TEST_REAL_SIMILAR(fm.diff_da, fm.theoretical_mass * (kPpm * 1e-6))
      // fm.theoretical_mass is the full equivalent (b5) theoretical = prefix[5].
      double theo_b5 = resMass('A') + resMass('C') + resMass('D') + resMass('E') + resMass('F');
      TEST_REAL_SIMILAR(fm.theoretical_mass, theo_b5)
      break;
    }
  }
  TEST_TRUE(found)
}
END_SECTION

/////////////////////////////////////////////////////////////
// T9: PPM-HONEST frame projection -- complement-flip yb-from-b, with an ambiguous PTM on the
//     COMPLEMENT of the sub-fragment. coveredAmbiguousInEquivFrame folds the phospho into the
//     equivalent ion even though the yb2 sub-fragment does not cover it, so equiv_ambiguous (+phospho)
//     != md.ambiguous_included (0). The invariant must still hold: diff_ppm == the sub-frame +10 ppm
//     (offset AND equiv_ambiguous cancel), and the WHOLE mod-inclusive equivalent theoretical scales.
/////////////////////////////////////////////////////////////

START_SECTION(calibrateAndScore_ppm_frame_projection_complement_flip_yb)
{
  // protein indices: A0 C1 D2 E3 F4 G5 H6 I7 K8 L9 M10 N11 P12 Q13 R14 S15
  std::string protein = "ACDEFGHIKLMNPQRS";
  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 3;                  // proteoform "EFGHIKLM" (region positions 1..8)
  ctx.region_end = 11;
  // Ambiguous PTM on region positions 3-4 (= "GH", absolute residues 5-6): on the COMPLEMENT of the
  // yb2 sub-fragment, but INSIDE the equivalent b7 ion.
  const double kPtm = 79.96633;          // phospho
  FragmentAnalysis::PTMSite site;
  site.start_position = 3; site.end_position = 4; site.position = 3; site.mass_shift = kPtm;
  ctx.ptm_sites.push_back(site);

  char frag_type = 'b';
  int frag_index = 6;                    // b6 precursor -> subseq = protein[3:9] = "EFGHIK"
  std::string subseq = MS3FragmentMatcher::extractSubsequence(protein, ctx, frag_type, frag_index);
  TEST_EQUAL(subseq, "EFGHIK")

  // Clean calibration ions: prefix b1,b2 cover subseq [1..2] = "EF" (no PTM overlap).
  double theo_b1 = resMass('E');
  double theo_b2 = resMass('E') + resMass('F');
  // Target: suffix yb2 covers subseq [5..6] = "IK" (absolute 7-8), does NOT cover the PTM at [3,4].
  double theo_yb2 = resMass('I') + resMass('K');

  const double kPpm = 10.0;
  DeconvolvedSpectrum ds(0);
  auto addPeak = [&](double m) { PeakGroup pg(1, 1, true); pg.setMonoisotopicMass(m); ds.push_back(pg); };
  addPeak(theo_b1);                          // 0 ppm
  addPeak(theo_b2);                          // 0 ppm
  addPeak(theo_yb2 * (1.0 + kPpm * 1e-6));   // +10 ppm (target)

  std::vector<const DeconvolvedSpectrum*> variants = {&ds};
  std::vector<FragmentAnalysis::ProteoformMatch> detailed;
  MS3FragmentMatcher::calibrateAndScore(
    variants, protein, ctx, frag_type, frag_index,
    MS3FragmentMatcher::LOOSE_TOLERANCE_PPM, 50.0, &detailed);

  TEST_EQUAL(detailed.size(), 1)
  TEST_REAL_SIMILAR(detailed[0].correction_factor, 1.0)

  bool found = false;
  for (const auto& fm : detailed[0].fragments)
  {
    if (fm.ion_type == "yb" && fm.ion_index == 2)
    {
      found = true;
      // Complement flip: yb2 (suffix) maps to the prefix equivalent b7 (covers absolute [0,6]).
      TEST_EQUAL(fm.equiv_type, "b")
      TEST_EQUAL(fm.equiv_index, 7)
      // Invariant survives equiv_ambiguous (+phospho) != md.ambiguous_included (0): diff_ppm == sub-frame +10 ppm.
      TEST_REAL_SIMILAR(fm.diff_ppm, kPpm)
      TEST_REAL_SIMILAR(fm.adjusted_mass, fm.theoretical_mass * (1.0 + kPpm * 1e-6))
      // Q3: the WHOLE mod-inclusive equivalent theoretical = b7 backbone (prefix[7]) + the ambiguous phospho.
      double prefix7 = resMass('A') + resMass('C') + resMass('D') + resMass('E')
                     + resMass('F') + resMass('G') + resMass('H');
      TEST_REAL_SIMILAR(fm.theoretical_mass, prefix7 + kPtm)
      break;
    }
  }
  TEST_TRUE(found)
}
END_SECTION

/////////////////////////////////////////////////////////////
// T10: PPM-HONEST frame projection -- cross map, b-from-y-precursor, truncated proteoform
//      (region_end < P so the P-end y-frame offset is exercised). A prefix b ion of a y-precursor
//      maps to a y equivalent; diff_ppm must equal the injected sub-frame +10 ppm.
/////////////////////////////////////////////////////////////

START_SECTION(calibrateAndScore_ppm_frame_projection_cross_map_y)
{
  // protein indices: A0 C1 D2 E3 F4 G5 H6 I7 K8 L9  (P = 10)
  std::string protein = "ACDEFGHIKL";
  MS3FragmentMatcher::ProteoformContext ctx;
  ctx.region_start = 1;                  // proteoform "CDEFGHI" (region_end 8 < P 10)
  ctx.region_end = 8;

  char frag_type = 'y';
  int frag_index = 4;                    // y4 precursor -> subseq = protein[4:8] = "FGHI"
  std::string subseq = MS3FragmentMatcher::extractSubsequence(protein, ctx, frag_type, frag_index);
  TEST_EQUAL(subseq, "FGHI")

  // Prefix b-ions of the subseq (a y-precursor emits a,b,y; a prefix b maps CROSS to a y equivalent).
  double theo_b1 = resMass('F');
  double theo_b2 = resMass('F') + resMass('G');
  double theo_b3 = resMass('F') + resMass('G') + resMass('H');

  const double kPpm = 10.0;
  DeconvolvedSpectrum ds(0);
  auto addPeak = [&](double m) { PeakGroup pg(1, 1, true); pg.setMonoisotopicMass(m); ds.push_back(pg); };
  addPeak(theo_b1);                          // 0 ppm
  addPeak(theo_b2 * (1.0 + kPpm * 1e-6));    // +10 ppm (target)
  addPeak(theo_b3);                          // 0 ppm

  std::vector<const DeconvolvedSpectrum*> variants = {&ds};
  std::vector<FragmentAnalysis::ProteoformMatch> detailed;
  MS3FragmentMatcher::calibrateAndScore(
    variants, protein, ctx, frag_type, frag_index,
    MS3FragmentMatcher::LOOSE_TOLERANCE_PPM, 50.0, &detailed);

  TEST_EQUAL(detailed.size(), 1)
  TEST_REAL_SIMILAR(detailed[0].correction_factor, 1.0)

  bool found = false;
  for (const auto& fm : detailed[0].fragments)
  {
    if (fm.ion_type == "b" && fm.ion_index == 2)
    {
      found = true;
      // Cross map: prefix b from a y-precursor -> y equivalent, index = P - end + N - idx = 10 - 8 + 4 - 2 = 4.
      TEST_EQUAL(fm.equiv_type, "y")
      TEST_EQUAL(fm.equiv_index, 4)
      TEST_REAL_SIMILAR(fm.diff_ppm, kPpm)
      TEST_REAL_SIMILAR(fm.adjusted_mass, fm.theoretical_mass * (1.0 + kPpm * 1e-6))
      break;
    }
  }
  TEST_TRUE(found)
}
END_SECTION

END_TEST
