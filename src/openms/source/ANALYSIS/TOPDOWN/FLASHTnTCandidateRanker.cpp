// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTCandidateRanker.h>
#include "FLASHTnTCandidateModel.h"
#include <boost/math/special_functions/gamma.hpp>
#include <algorithm>
#include <cmath>

namespace OpenMS
{
namespace
{
// the values the model was trained with (TD_survey/scripts/features_pf.py)
constexpr double match_ppm = 10.0, h2o = 18.010565, isotope = 1.00235, error_cap = 9.99, window = 0.05;
// a cut before the start plus one of these is a residue boundary: none, acetyl, formyl, methyl, dimethyl, trimethyl, myristoyl,
// palmitoyl, pyro-Glu from Gln and from Glu; after the end minus one of these: none, amidation
constexpr std::array<double, 10> n_mods {0.0, 42.010565, 27.994915, 14.01565, 28.0313, 42.04695, 210.198366, 238.229666, -17.026549, -18.010565};
constexpr std::array<double, 2> c_mods {0.0, 0.984016};

// the offset of a side: the median of the shifts in the first window of 0.05 Da (each way) around a shift holding the most of them
double sideOffset(const double* shifts, Size count)
{
  std::vector<double> sorted(shifts, shifts + count);
  std::sort(sorted.begin(), sorted.end());
  Size best = 0, lo_best = 0, hi_best = 0;
  for (Size k = 0; k < sorted.size(); k++)
  {
    const Size lo = std::lower_bound(sorted.begin(), sorted.end(), sorted[k] - window) - sorted.begin();
    const Size hi = std::upper_bound(sorted.begin(), sorted.end(), sorted[k] + window) - sorted.begin();
    if (hi - lo > best)
    {
      best = hi - lo;
      lo_best = lo;
      hi_best = hi;
    }
  }
  const Size n = hi_best - lo_best, mid = lo_best + n / 2;
  return n % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
}

// the smallest distance of a residue sum from the cut plus (sign 1) or minus (sign -1) a modification, with +-1 isotope
template<Size N> double boundaryError(const std::vector<double>& sums, double cut, const std::array<double, N>& mods, double sign)
{
  double best = error_cap;
  for (double mod : mods)
  {
    for (double iso : {0.0, isotope, -isotope})
    {
      const double x = cut + sign * (mod + iso);
      const Size j = std::lower_bound(sums.begin(), sums.end(), x) - sums.begin();
      const double lo = sums[j == 0 ? 0 : j - 1], hi = sums[std::min(j, sums.size() - 1)];
      best = std::min(best, std::min(std::abs(x - lo), std::abs(hi - x)));
    }
  }
  return best;
}

// -log10 P(Poisson >= the masses within match_ppm of the side's ladder): the ladder from the residue sums at the cut (prefix ions above
// 50 Da; suffix ions, H2O included, above H2O + 50 Da), lambda from its density
double ladderChance(const std::vector<double>& masses, const std::vector<double>& sums, double cut, double offset, bool suffix)
{
  std::vector<double> ions;
  for (Size k = suffix ? 0 : 1; k < sums.size(); k++)
  {
    const double ion = sums[k] - cut + offset + (suffix ? h2o : 0.0);
    if (ion > (suffix ? h2o + 50 : 50)) ions.push_back(ion);
  }
  if (ions.empty()) return 0;
  int k = 0;
  double lambda = 0;
  const double top = ions.back(), density = (double)ions.size() / std::max(top, 1.0);
  for (double m : masses)
  {
    const Size j = std::lower_bound(ions.begin(), ions.end(), m) - ions.begin();
    const double lo = ions[j == 0 ? 0 : j - 1], hi = ions[std::min(j, ions.size() - 1)];
    k += std::min(std::abs(m - lo), std::abs(hi - m)) <= match_ppm * 1e-6 * m;
    if (m <= top) lambda += std::min(1.0, 2 * match_ppm * 1e-6 * m * density);
  }
  if (k == 0) return 0;
  return -std::log10(std::max(boost::math::gamma_p((double)k, std::max(lambda, 1e-12)), 1e-300));
}
} // namespace

std::array<double, FLASHTnTCandidateRanker::feature_count> FLASHTnTCandidateRanker::features(const Spectrum& spectrum, const Candidate& candidate)
{
  const auto& pre = *candidate.prefix;
  const auto& suf = *candidate.suffix;
  double best_chance = 0, n_bound = error_cap, c_bound = error_cap, c_intact = 0;
  if (candidate.n_count > 0)
  {
    const double cut = sideOffset(candidate.n_shifts, candidate.n_count) + spectrum.prefix_offset; // residue mass before the start
    n_bound = boundaryError(pre, cut, n_mods, 1.0);
    best_chance = std::max(best_chance, ladderChance(spectrum.masses, pre, cut, spectrum.prefix_offset, false));
  }
  if (candidate.c_count > 0)
  {
    const double cut = sideOffset(candidate.c_shifts, candidate.c_count) + h2o + spectrum.suffix_offset; // residue mass after the end
    c_bound = boundaryError(suf, cut, c_mods, -1.0);
    c_intact = std::abs(cut) < 1.1;
    best_chance = std::max(best_chance, ladderChance(spectrum.masses, suf, cut, spectrum.suffix_offset, true));
  }
  return {candidate.score, best_chance, n_bound, c_bound, c_intact};
}

double FLASHTnTCandidateRanker::score(const std::array<double, feature_count>& features)
{
  using namespace FLASHTnTCandidateModel;
  static_assert(FLASHTnTCandidateModel::feature_count == FLASHTnTCandidateRanker::feature_count, "the model reads the features of features()");
  double s = 0;
  for (int t = 0; t < tree_count; t++)
  {
    int n = tree_root[t];
    while (n >= 0) n = features[nodes[n].feature] <= nodes[n].threshold ? nodes[n].left : nodes[n].right;
    s += leaf_values[~n];
  }
  return s;
}
} // namespace OpenMS
