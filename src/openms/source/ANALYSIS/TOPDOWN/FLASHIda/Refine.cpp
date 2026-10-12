// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Refine.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <random>

namespace OpenMS
{
  namespace
  {
    constexpr double kWater = 18.010565;
    constexpr double kAmmonia = 17.026549;
    constexpr double kLossSlack = 0.005; ///< Da; the neutral-loss test's slack on either side of a region
    constexpr int kHMax = 60;            ///< the highest peak height FDR is evaluated at
    constexpr double kEps = 1e-9;        ///< grid units; a match exactly on a grid point counts

    double tolerance(double mass, double tol_ppm)
    {
      return mass * tol_ppm * 1e-6;
    }

    /// The prototype's ladder_from_residues for one ion pair: prefix = cumulative residues + prefix_shift,
    /// suffix (site-indexed) = the remaining residues + suffix_shift.
    void ladderFromResidues(const std::vector<double>& res, double prefix_shift, double suffix_shift, std::vector<double>& prefix,
                            std::vector<double>& suffix)
    {
      const size_t L = res.size();
      prefix.clear();
      suffix.clear();
      if (L < 2) return;
      prefix.resize(L - 1);
      suffix.resize(L - 1);
      double cum = 0;
      for (size_t i = 0; i + 1 < L; ++i)
      {
        cum += res[i];
        prefix[i] = cum + prefix_shift;
      }
      cum = 0;
      for (size_t i = L - 1; i >= 1; --i) // the suffix of cleavage site i holds residues i..L-1
      {
        cum += res[i];
        suffix[i - 1] = cum + suffix_shift;
      }
    }

    /// The prototype's Scan.matches: the number of ions t with a mass m, among the two neighbours of t + shift,
    /// such that |m - (t + shift)| <= tol(m).
    int countMatching(const std::vector<double>& theo, const std::vector<double>& masses_sorted, double tol_ppm, double shift)
    {
      int count = 0;
      for (double t : theo)
      {
        const double x = t + shift;
        const auto it = std::lower_bound(masses_sorted.begin(), masses_sorted.end(), x);
        bool hit = false;
        if (it != masses_sorted.end()) hit = std::fabs(*it - x) <= tolerance(*it, tol_ppm);
        if (!hit && it != masses_sorted.begin()) hit = std::fabs(*(it - 1) - x) <= tolerance(*(it - 1), tol_ppm);
        if (hit) ++count;
      }
      return count;
    }

    std::vector<double> concat(const std::vector<double>& a, const std::vector<double>& b)
    {
      std::vector<double> out;
      out.reserve(a.size() + b.size());
      out.insert(out.end(), a.begin(), a.end());
      out.insert(out.end(), b.begin(), b.end());
      return out;
    }
  } // namespace

  Refine::Ladder Refine::ladder(const String& sequence, const std::vector<std::string>& ion_types)
  {
    std::string prefix_type, suffix_type;
    for (const auto& ion : ion_types)
    {
      if (prefix_type.empty() && (ion == "a" || ion == "b" || ion == "c")) prefix_type = ion;
      if (suffix_type.empty() && (ion == "x" || ion == "y" || ion == "z")) suffix_type = ion;
    }
    if (prefix_type.empty()) prefix_type = "b";
    if (suffix_type.empty()) suffix_type = "y";

    Ladder out;
    const size_t L = sequence.size();
    if (L < 2)
    {
      // No cleavage site, so no ladder to difference: the residue mass straight from the engine's table.
      for (char aa : sequence)
      {
        const Residue* r = ResidueDB::getInstance()->getResidue(aa == 'I' ? 'L' : aa);
        out.residues.push_back(r != nullptr && aa != 'X' ? r->getMonoWeight(Residue::Internal) : 0.0);
      }
      for (double r : out.residues) out.sequence_mass += r;
      out.sequence_mass += kWater;
      return out;
    }

    std::map<char, std::vector<double>> by_type;
    FragmentAnalysis::computePTMAdjustedFragmentMasses(sequence, {}, {prefix_type, suffix_type}, by_type);
    const std::vector<double>& p = by_type[prefix_type[0]]; // p[i] = the prefix ion with i+1 residues, i = 0..L-1
    const std::vector<double>& q = by_type[suffix_type[0]]; // q[i] = the suffix ion with i+1 residues

    out.prefix.assign(p.begin(), p.begin() + (L - 1));
    out.suffix.resize(L - 1);
    for (size_t n = 1; n < L; ++n) out.suffix[n - 1] = q[L - n - 1];

    // Residues by differencing the engine's series; residue 0 = (all L) - (the last L-1) on the suffix side.
    out.residues.resize(L);
    out.residues[0] = q[L - 1] - q[L - 2];
    for (size_t i = 1; i < L; ++i) out.residues[i] = p[i] - p[i - 1];
    for (double r : out.residues) out.sequence_mass += r;
    out.sequence_mass += kWater;
    return out;
  }

  std::vector<int> Refine::curve(const std::vector<double>& theo, const std::vector<double>& masses_sorted, double tol_ppm, double lo,
                                 double step, int n)
  {
    std::vector<int> f(std::max(n, 0), 0);
    if (n <= 0 || theo.empty() || masses_sorted.empty()) return f;
    std::vector<long long> diff(n + 1, 0);
    for (double t : theo)
    {
      long long prev_e = 0;
      bool first = true;
      for (double m : masses_sorted)
      {
        const double tol = tolerance(m, tol_ppm);
        long long s = static_cast<long long>(std::ceil((m - tol - t - lo) / step - kEps));
        long long e = static_cast<long long>(std::floor((m + tol - t - lo) / step + kEps));
        s = std::clamp<long long>(s, 0, n); // n = past the end
        e = std::clamp<long long>(e, -1, n - 1);
        // one ion counts once: along one ion the ends ascend, so an interval starts behind the end of the one before
        if (!first) s = std::max(s, prev_e + 1);
        first = false;
        prev_e = e;
        if (s <= e)
        {
          ++diff[s];
          --diff[e + 1];
        }
      }
    }
    long long run = 0;
    for (int i = 0; i < n; ++i)
    {
      run += diff[i];
      f[i] = static_cast<int>(run);
    }
    return f;
  }

  std::vector<bool> Refine::explained(const std::vector<double>& theo, const std::vector<double>& masses_sorted, double tol_ppm)
  {
    std::vector<double> sorted(theo);
    std::sort(sorted.begin(), sorted.end());
    std::vector<bool> out(masses_sorted.size(), false);
    for (size_t k = 0; k < masses_sorted.size(); ++k)
    {
      const double m = masses_sorted[k];
      const double tol = tolerance(m, tol_ppm);
      const auto it = std::lower_bound(sorted.begin(), sorted.end(), m);
      bool hit = it != sorted.end() && std::fabs(m - *it) <= tol;
      if (!hit && it != sorted.begin()) hit = std::fabs(m - *(it - 1)) <= tol;
      out[k] = hit;
    }
    return out;
  }

  std::vector<int> Refine::regionCounts(const std::vector<int>& curve, int hmax)
  {
    std::vector<int> counts(std::max(hmax, 0) + 1, 0);
    if (hmax <= 0) return counts;
    // A run of f >= h starts at i exactly when f[i-1] < h <= f[i] (f[-1] = 0): one difference array over h.
    std::vector<int> diff(hmax + 2, 0);
    int prev = 0;
    for (int v : curve)
    {
      if (v > prev)
      {
        const int h_lo = std::max(prev + 1, 1);
        const int h_hi = std::min(v, hmax);
        if (h_lo <= h_hi)
        {
          ++diff[h_lo];
          --diff[h_hi + 1];
        }
      }
      prev = v;
    }
    int run = 0;
    for (int h = 1; h <= hmax; ++h)
    {
      run += diff[h];
      counts[h] = run;
    }
    return counts;
  }

  void Refine::flagNeutralLosses(std::vector<Level>& levels)
  {
    const double losses[] = {kWater, kAmmonia};
    for (Level& level : levels)
    {
      std::vector<double> refs = {0.0};
      for (const Level& other : levels)
      {
        if (other.signal > level.signal) refs.push_back(other.shift);
      }
      for (double r : refs)
      {
        for (double v : losses)
        {
          if (level.lo - kLossSlack <= r - v && r - v <= level.hi + kLossSlack) level.neutral_loss = true;
        }
      }
    }
  }

  Refine::Detect Refine::detect(const std::vector<double>& masses_sorted, const Ladder& ladder, const RefineSettings& s)
  {
    Detect out;
    out.fdr_by_height.assign(kHMax + 1, -1.0);

    const double lo = -s.max_shift_da;
    const double step = s.grid_da;
    const int n = static_cast<int>(std::lround(2.0 * s.max_shift_da / s.grid_da)) + 1;

    const std::vector<double> T = concat(ladder.prefix, ladder.suffix);
    const std::vector<bool> ex = explained(T, masses_sorted, s.tolerance_ppm);
    std::vector<double> m0, unexplained;
    for (size_t k = 0; k < masses_sorted.size(); ++k) (ex[k] ? m0 : unexplained).push_back(masses_sorted[k]);
    out.observed_at_zero = static_cast<int>(m0.size());

    // The Signal: Observed - Expected.
    std::vector<int> S = curve(T, masses_sorted, s.tolerance_ppm, lo, step, n);
    const std::vector<int> E = curve(T, m0, s.tolerance_ppm, lo, step, n);
    for (int i = 0; i < n; ++i) S[i] -= E[i];
    const std::vector<int> nS = regionCounts(S, kHMax);

    // The null: K shuffled sequences against the unexplained masses, one stream.
    std::vector<double> decoy_total(kHMax + 1, 0.0);
    const size_t L = ladder.residues.size();
    if (s.shuffles > 0 && !unexplained.empty() && L >= 2)
    {
      // The ion offsets of this alphabet, read back from the ladder: prefix_1 holds residue 0, suffix_{L-1} residue L-1.
      const double prefix_shift = ladder.prefix[0] - ladder.residues[0];
      const double suffix_shift = ladder.suffix[L - 2] - ladder.residues[L - 1];
      std::mt19937_64 rng(RefineSettings::kSeed);
      std::vector<double> residues(ladder.residues);
      std::vector<double> prefix, suffix;
      for (int k = 0; k < s.shuffles; ++k)
      {
        std::shuffle(residues.begin(), residues.end(), rng);
        ladderFromResidues(residues, prefix_shift, suffix_shift, prefix, suffix);
        const std::vector<int> nD = regionCounts(curve(concat(prefix, suffix), unexplained, s.tolerance_ppm, lo, step, n), kHMax);
        for (int h = 1; h <= kHMax; ++h) decoy_total[h] += nD[h];
      }
    }
    const double K = s.shuffles > 0 ? static_cast<double>(s.shuffles) : 1.0;
    for (int h = 1; h <= kHMax; ++h)
    {
      if (nS[h] > 0) out.fdr_by_height[h] = (decoy_total[h] / K) / nS[h];
    }

    for (int h = std::max(s.min_ions, 1); h <= kHMax; ++h)
    {
      if (nS[h] > 0 && out.fdr_by_height[h] <= s.fdr)
      {
        out.cut = h;
        break;
      }
    }
    if (out.cut == 0) return out;

    // The levels: maximal runs of S >= cut.
    for (int i = 0; i < n;)
    {
      if (S[i] < out.cut)
      {
        ++i;
        continue;
      }
      int j = i;
      int top = S[i];
      while (j + 1 < n && S[j + 1] >= out.cut)
      {
        ++j;
        top = std::max(top, S[j]);
      }
      std::vector<int> best; // the grid points at the maximum; the middle one is the shift
      for (int k = i; k <= j; ++k)
      {
        if (S[k] == top) best.push_back(k);
      }
      Level level;
      level.shift = lo + step * best[best.size() / 2];
      level.signal = top;
      level.lo = lo + step * i;
      level.hi = lo + step * j;
      level.n_prefix = countMatching(ladder.prefix, masses_sorted, s.tolerance_ppm, level.shift);
      level.n_suffix = countMatching(ladder.suffix, masses_sorted, s.tolerance_ppm, level.shift);
      out.levels.push_back(level);
      i = j + 1;
    }
    flagNeutralLosses(out.levels);
    return out;
  }
} // namespace OpenMS
