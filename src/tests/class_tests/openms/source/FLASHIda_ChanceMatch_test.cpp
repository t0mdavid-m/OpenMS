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
// ADR-0075: the chance-match p-value's pure functions (FLASHIda/ChanceMatch.h), and the window count the E-value shares
// with FLASHTnT's raw-score penalty (FLASHTnTAlgorithm::shiftedCandidateCount).

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ChanceMatch.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <cmath>
#include <optional>
#include <string>
#include <vector>

using namespace OpenMS;

START_TEST(FLASHIda_ChanceMatch, "$Id$")

/////////////////////////////////////////////////////////////

// log10 P(X >= k) against the direct sum of Poisson terms, in both branches: k <= lambda (one minus the head) and
// k > lambda (the upward series).
// FAILS UNDER: an off-by-one in the tail's first term, or a head that is not subtracted from 1.
START_SECTION(log10_poisson_tail_matches_the_direct_sum)
{
  TEST_REAL_SIMILAR(ChanceMatch::log10PoissonTail(0, 0.5), 0.0)
  TEST_REAL_SIMILAR(ChanceMatch::log10PoissonTail(-3, 0.5), 0.0)

  // k <= lambda: 1 - e^-5 (1 + 5 + 12.5)
  TEST_REAL_SIMILAR(ChanceMatch::log10PoissonTail(3, 5.0), std::log10(1.0 - std::exp(-5.0) * (1.0 + 5.0 + 12.5)))

  // k > lambda: the lesson's example, 20 matches at lambda 0.1236 (about 2.7e-37)
  const double lambda = 0.1236;
  long double direct = 0.0L;
  for (int j = 20; j <= 60; ++j) direct += std::exp(static_cast<long double>(-lambda + j * std::log(lambda) - std::lgamma(j + 1.0)));
  TEST_REAL_SIMILAR(ChanceMatch::log10PoissonTail(20, lambda), static_cast<double>(std::log10(direct)))
  TEST_EQUAL(ChanceMatch::log10PoissonTail(20, lambda) < -36.5 && ChanceMatch::log10PoissonTail(20, lambda) > -36.7, true)
}
END_SECTION

// A long ladder with many matches puts p below the smallest double; the log10 stays finite.
// FAILS UNDER: the tail computed as a double and logged afterwards (log10(0) = -inf).
START_SECTION(log10_poisson_tail_does_not_underflow)
{
  const double lp = ChanceMatch::log10PoissonTail(150, 0.12);
  TEST_EQUAL(std::isfinite(lp), true)
  TEST_EQUAL(lp < -300.0, true)
  TEST_EQUAL(std::isinf(ChanceMatch::log10PoissonTail(1, 0.0)), true)   // no masses: a match cannot be chance
}
END_SECTION

// lambda = N x the share of the chain's mass range the tolerance windows cover; non-positive masses carry no window.
// FAILS UNDER: a window of tol instead of 2 x tol, or a missing N.
START_SECTION(chance_lambda_is_n_times_window_share)
{
  TEST_REAL_SIMILAR(ChanceMatch::chanceLambda(10, {1000.0, 2000.0, -5.0}, 10.0, 5000.0), 10.0 * 2.0 * 10.0e-6 * 3000.0 / 5000.0)
  TEST_REAL_SIMILAR(ChanceMatch::chanceLambda(10, {1000.0}, 10.0, 0.0), 0.0)
}
END_SECTION

// The chain's mass: residues (I as L), shifts and one water.
// FAILS UNDER: a missing water or a dropped shift.
START_SECTION(chain_mass_is_residues_shifts_and_water)
{
  const double water = EmpiricalFormula("H2O").getMonoWeight();
  const double gly = ResidueDB::getInstance()->getResidue('G')->getMonoWeight(Residue::Internal);
  const double leu = ResidueDB::getInstance()->getResidue('L')->getMonoWeight(Residue::Internal);
  TEST_REAL_SIMILAR(ChanceMatch::chainMass("GGI", {79.966331}, FixedModifications::none()), 2 * gly + leu + 79.966331 + water)
}
END_SECTION

// identification.tsv's form: three significant digits, written from the log10, with the mantissa's carry; -1 absent.
// FAILS UNDER: a double formatted directly (4.10e-401 prints as 0), or 9.996e-5 written as "10.00e-05".
START_SECTION(format_log10_scientific)
{
  TEST_EQUAL(ChanceMatch::formatLog10Scientific(std::log10(2.71e-37)), "2.71e-37")
  TEST_EQUAL(ChanceMatch::formatLog10Scientific(-400.38722), "4.10e-401")
  TEST_EQUAL(ChanceMatch::formatLog10Scientific(0.0), "1.00e+00")
  TEST_EQUAL(ChanceMatch::formatLog10Scientific(std::log10(9.996e-5)), "1.00e-04")
  TEST_EQUAL(ChanceMatch::formatLog10Scientific(std::log10(250.0)), "2.50e+02")
  TEST_EQUAL(ChanceMatch::formatLog10Scientific(std::nullopt), "-1")
  TEST_EQUAL(ChanceMatch::formatLog10Short(std::log10(3e-12)), "3e-12")
  TEST_EQUAL(ChanceMatch::formatLog10Short(std::nullopt), "-")
}
END_SECTION

// The E-value's window count is FLASHTnT's own: every window with no mass, and the windows within the shift bound of it.
// FAILS UNDER: a copy of the count that drifts from the one FLASHTnT's penalty uses.
START_SECTION(shifted_candidate_count_is_flashtnt_s)
{
  TEST_REAL_SIMILAR(FLASHTnTAlgorithm::shiftedCandidateCount("ACDEF", -1.0, 100.0, FixedModifications::none()), 15.0)
  const double water = EmpiricalFormula("H2O").getMonoWeight();
  const double gly = ResidueDB::getInstance()->getResidue('G')->getMonoWeight(Residue::Internal);
  TEST_REAL_SIMILAR(FLASHTnTAlgorithm::shiftedCandidateCount("GGGG", 2 * gly + water, 0.5, FixedModifications::none()), 3.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
