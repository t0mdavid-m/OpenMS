// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <cmath>
#include <cstdint>
#include <functional>
#include <set>
#include <tuple>
#include <vector>
///////////////////////////

using namespace OpenMS;

// The label-blind rules of FLASHTnT that can be checked without a search: the entry order that breaks candidate-cut ties, the protein
// evidence rule (option protein_evidence) and the blind modification count it reads, the grouping of small proteoform classes (option
// merge_small_classes), the window masses of DeconvMassInfo (option precursor_reselect), the rounding of the hot loops
// (FLASHHelperClasses::roundHalfAway), and the defaults of the options
START_TEST(FLASHTnTAlgorithm, "$Id$")

START_SECTION((static bool proteinEvidenceSupported(Size proteoforms, bool strong)))
{
  // two distinct proteoforms, or one strong proteoform
  TEST_EQUAL(FLASHTnTAlgorithm::proteinEvidenceSupported(1, false), false)
  TEST_EQUAL(FLASHTnTAlgorithm::proteinEvidenceSupported(2, false), true)
  TEST_EQUAL(FLASHTnTAlgorithm::proteinEvidenceSupported(1, true), true)
  TEST_EQUAL(FLASHTnTAlgorithm::proteinEvidenceSupported(0, false), false)
}
END_SECTION

START_SECTION((static Size blindModCount(const ProteinHit& hit)))
{
  ProteinHit h;
  h.setMetaValue("Modifications", DoubleList {15.99, 42.01});
  TEST_EQUAL(FLASHTnTAlgorithm::blindModCount(h), 2) // without the kinds, every modification counts as blind
  h.setMetaValue("ModificationVariable", IntList {1, 0}); // one variable, one blind
  TEST_EQUAL(FLASHTnTAlgorithm::blindModCount(h), 1)
  h.setMetaValue("ModificationVariable", IntList {1, 2}); // one variable, one fixed by the database entry
  TEST_EQUAL(FLASHTnTAlgorithm::blindModCount(h), 0)
  h.setMetaValue("BlindModCount", 1); // the search's own count wins
  TEST_EQUAL(FLASHTnTAlgorithm::blindModCount(h), 1)
}
END_SECTION

START_SECTION((static std::uint64_t entryOrder(std::uint64_t index)))
{
  // deterministic, collision-free on the first entries, and not the database order: about half the consecutive pairs ascend
  std::set<std::uint64_t> seen;
  Size ascending = 0;
  for (std::uint64_t i = 0; i < 1000; i++)
  {
    TEST_EQUAL(FLASHTnTAlgorithm::entryOrder(i), FLASHTnTAlgorithm::entryOrder(i))
    seen.insert(FLASHTnTAlgorithm::entryOrder(i));
    if (i > 0 && FLASHTnTAlgorithm::entryOrder(i) > FLASHTnTAlgorithm::entryOrder(i - 1)) ascending++;
  }
  TEST_EQUAL(seen.size(), 1000)
  TEST_EQUAL(ascending > 400 && ascending < 600, true)
  TEST_EQUAL(FLASHTnTAlgorithm::entryOrder(0) == 0, false)
}
END_SECTION

START_SECTION((static std::vector<std::vector<int>> groupClasses(const std::vector<int>& class_order, const std::vector<Size>& class_size, Size min_size, bool by_kind, const std::function<int(int)>& kind_of)))
{
  // classes as in FLASHTnT: 0 determined, 1 both termini, 2 open terminus, 3 precursor ladder, 4 terminal ladder
  const std::vector<int> class_order = {0, 1, 2, 3, 4};
  const std::vector<int> kind_order = {0, 3, 4, 1, 2};
  const auto kind_of = [](int c) { return c == 0 ? 0 : 1; };
  // E. coli-like run: a small undetermined group (93 + 5 < 100) before a large precursor-ladder class
  std::vector<Size> size = {760, 93, 5, 156, 577};
  auto g = FLASHTnTAlgorithm::groupClasses(class_order, size, 100, false, kind_of);
  TEST_EQUAL(g.size(), 3)
  TEST_EQUAL(g[1].size(), 3) // true: the undetermined group absorbs the precursor-ladder class
  TEST_EQUAL(g[1][2], 3)
  g = FLASHTnTAlgorithm::groupClasses(kind_order, size, 100, true, kind_of);
  TEST_EQUAL(g.size(), 3) // kind: determined | precursor ladder | terminal ladder + the small undetermined group (last, joins the group before it)
  TEST_EQUAL(g[1].size(), 1)
  TEST_EQUAL(g[1][0], 3)
  TEST_EQUAL(g[2].size(), 3)
  TEST_EQUAL(g[2][0], 4)
  // an undetermined group of at least 100 stays its own group, in both walks
  size = {760, 93, 10, 156, 577};
  g = FLASHTnTAlgorithm::groupClasses(class_order, size, 100, false, kind_of);
  TEST_EQUAL(g.size(), 4)
  g = FLASHTnTAlgorithm::groupClasses(kind_order, size, 100, true, kind_of);
  TEST_EQUAL(g.size(), 4) // determined | precursor ladder | terminal ladder | undetermined
  TEST_EQUAL(g[3].size(), 2)
  TEST_EQUAL(g[3][0], 1)
  // empty classes take no part
  size = {760, 93, 10, 0, 577};
  g = FLASHTnTAlgorithm::groupClasses(kind_order, size, 100, true, kind_of);
  TEST_EQUAL(g.size(), 3)
  TEST_EQUAL(g[1].size(), 1)
  TEST_EQUAL(g[1][0], 4)
  // a tiny run: the determined class alone is small and joins the other kind's first group
  size = {30, 20, 5, 10, 50};
  g = FLASHTnTAlgorithm::groupClasses(kind_order, size, 100, true, kind_of);
  TEST_EQUAL(g.size(), 1)
  TEST_EQUAL(g[0].size(), 5)
}
END_SECTION

START_SECTION((static String DeconvolvedSpectrum::formatWindowMasses(const std::vector<std::tuple<double, double>>& map)))
{
  // the isolation-window masses FLASHDeconv writes into DeconvMassInfo ("windowmasses=") and FLASHTnT reads back
  std::vector<std::tuple<double, double>> map = {{0.0, 1e6}, {8565.4321, 1234.5}, {8579.4477, 99.0}};
  const String s = DeconvolvedSpectrum::formatWindowMasses(map);
  TEST_EQUAL(s.hasSubstring("8565.432"), true)
  TEST_EQUAL(s.hasSubstring("1e+06"), false) // an entry without a mass (<= 0) is left out
  TEST_EQUAL(s.hasSubstring(";"), false)     // no field separator of DeconvMassInfo
  const auto back = DeconvolvedSpectrum::parseWindowMasses(s);
  TEST_EQUAL(back.size(), 2)
  TEST_REAL_SIMILAR(std::get<0>(back[0]), 8565.4321)
  TEST_REAL_SIMILAR(std::get<1>(back[1]), 99.0)
  TEST_EQUAL(DeconvolvedSpectrum::parseWindowMasses("").size(), 0)
  TEST_EQUAL(DeconvolvedSpectrum::parseWindowMasses("garbage,1:2").size(), 1)
}
END_SECTION

START_SECTION((static double FLASHHelperClasses::roundHalfAway(double x)))
{
  // std::round without the libm call: the same value for halfway cases, negative values, the doubles next to a half, values at and
  // beyond 2^52 (integers), beyond 2^63 (converted to Size at some call sites) and a sweep of quarter values
  std::vector<double> xs = {0.0, -0.0, 0.5, -0.5, 1.5, -1.5, 2.5, -2.5, 0.49999999999999994, -0.49999999999999994, std::nextafter(0.5, 1.0),
                            std::nextafter(2.5, 0.0), 1e-300, 4503599627370495.5, -4503599627370495.5, 4503599627370496.0, 9007199254740993.0,
                            9.3e18, 1.7e19, 1e300, -1e300};
  for (int i = -2000; i <= 2000; i++) xs.push_back(i * 0.25 + i * 1e-9);
  int differ = 0;
  for (double x : xs) differ += ! (FLASHHelperClasses::roundHalfAway(x) == std::round(x));
  TEST_EQUAL(differ, 0)
  TEST_EQUAL((Size)FLASHHelperClasses::roundHalfAway(1.7e19), (Size)std::round(1.7e19))
  TEST_EQUAL((int)FLASHHelperClasses::roundHalfAway(-2.5), -3)
  TEST_EQUAL((int)FLASHHelperClasses::roundHalfAway(2.4999999999999996), 2)
}
END_SECTION

START_SECTION((FLASHTnTAlgorithm()))
{
  // the defaults of the FDR and mass options, and the values each of them takes
  const Param p = FLASHTnTAlgorithm().getDefaults();
  TEST_STRING_EQUAL(p.getValue("terminal_ladder_class").toString(), "true")
  TEST_EQUAL(p.getValidStrings("terminal_ladder_class").size(), 2)
  TEST_EQUAL(p.hasTag("terminal_ladder_class", "advanced"), true)
  TEST_STRING_EQUAL(p.getValue("protein_evidence").toString(), "proteoforms")
  TEST_EQUAL(p.getValidStrings("protein_evidence").size(), 2)
  TEST_REAL_SIMILAR((double)p.getValue("protein_evidence_min_mass"), 3000.0)
  TEST_STRING_EQUAL(p.getValue("merge_small_classes").toString(), "kind")
  TEST_EQUAL(p.getValidStrings("merge_small_classes").size(), 2)
  TEST_STRING_EQUAL(p.getValue("precursor_reselect").toString(), "true")
  TEST_STRING_EQUAL(p.getValue("one_proteoform_per_feature").toString(), "fitted")
  TEST_EQUAL(p.getValidStrings("one_proteoform_per_feature").size(), 3)
  TEST_EQUAL(p.hasTag("one_proteoform_per_feature", "advanced"), true)
  // FLASHTnT's extender parameters differ from the extender's own defaults
  TEST_STRING_EQUAL(p.getValue("ex:skip_precursor_inference").toString(), "corroborated")
  TEST_EQUAL(p.getValidStrings("ex:skip_precursor_inference").size(), 3)
  TEST_STRING_EQUAL(p.getValue("ex:ladder_end_penalty").toString(), "false")
  TEST_EQUAL(p.getValidStrings("ex:ladder_end_penalty").size(), 2)
  TEST_EQUAL(p.hasTag("ex:ladder_end_penalty", "advanced"), true)
}
END_SECTION

END_TEST
