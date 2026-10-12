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

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <set>

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

  // ======================================================================================================================
  // ASSEMBLE -- column_dp_v85.py on lenient_core.py's matrices, one to one
  // ======================================================================================================================

  namespace
  {
    constexpr long long kNeg = -1000000000LL; ///< the prototype's NEG = -10 ** 9
    constexpr long long kNegHalf = kNeg / 2;  ///< NEG // 2 (exact: NEG is even)
    constexpr long long kIon = 1000;          ///< a' = 1000 a: on a tie the path with fewer mass shifts ...
    constexpr long long kJump = 1001;         ///< ... a jump costs 1001 (column_dp_v85.py, "on a tie, fewer mass shifts")

    /// Python's floor division and modulo for a positive divisor.
    long long floorDiv(long long a, long long b)
    {
      long long q = a / b;
      if ((a % b != 0) && ((a < 0) != (b < 0))) --q;
      return q;
    }
    long long floorMod(long long a, long long b)
    {
      return a - floorDiv(a, b) * b;
    }

    bool same(double x, double y)
    {
      return std::fabs(x - y) <= RefineSettings::kSame;
    }

    bool isin(double x, const std::vector<double>& S)
    {
      return std::any_of(S.begin(), S.end(), [x](double s) { return same(x, s); });
    }

    /// lenient_core.idx_of: the index of the observed mass within tolerance of x, or -1 (the higher of two neighbours when
    /// both match).
    int idxOf(double x, const std::vector<double>& M, double tol_ppm)
    {
      const long long i = std::lower_bound(M.begin(), M.end(), x) - M.begin();
      int best = -1;
      for (long long j : {i - 1, i})
      {
        if (j >= 0 && j < static_cast<long long>(M.size()) && std::fabs(M[j] - x) <= tolerance(M[j], tol_ppm)) best = static_cast<int>(j);
      }
      return best;
    }

    /// The ion types of a ladder, read back from its offsets (prefix_1 holds residue 0, suffix_{L-1} residue L-1).
    char prefixType(const Refine::Ladder& lad)
    {
      const double off = lad.prefix.front() - lad.residues.front();
      if (std::fabs(off - 17.026549) < 0.01) return 'c';
      if (std::fabs(off + 27.994915) < 0.01) return 'a';
      return 'b';
    }
    char suffixType(const Refine::Ladder& lad)
    {
      const double off = lad.suffix.back() - lad.residues.back();
      if (std::fabs(off - kWater) < 0.01) return 'y';
      if (std::fabs(off - 43.989829) < 0.01) return 'x';
      return 'z';
    }

    bool isPrefix(char t)
    {
      return t == 'a' || t == 'b' || t == 'c';
    }

    /// path_dp of column_dp_v85.py over one order of levels: Fp(n, j) with site n, Gp(n, j) without it, N unscaled.
    /// Aj[j][n] = a(n, level j) for n = 1..NS (Aj[j][0] = 0).
    void pathDp(const std::vector<std::vector<long long>>& Aj, int NS, std::vector<std::vector<long long>>& Fp,
                std::vector<std::vector<long long>>& Gp)
    {
      const int k = static_cast<int>(Aj.size()) - 1;
      Fp.assign(NS + 1, std::vector<long long>(k + 1, kNeg));
      Fp[0][0] = 0;
      for (int n = 1; n <= NS; ++n)
      {
        for (int j = 0; j <= k; ++j)
        {
          long long v = Fp[n - 1][j];
          if (j >= 1) v = std::max(v, Fp[n - 1][j - 1]);
          if (v > kNegHalf) Fp[n][j] = v + Aj[j][n];
        }
      }
      Gp.assign(NS + 1, std::vector<long long>(k + 1, kNeg));
      Gp[NS][k] = 0;
      for (int n = NS - 1; n >= 0; --n)
      {
        for (int j = 0; j <= k; ++j)
        {
          long long v = Gp[n + 1][j] > kNegHalf ? Gp[n + 1][j] + Aj[j][n + 1] : kNeg;
          if (j < k && Gp[n + 1][j + 1] > kNegHalf) v = std::max(v, Gp[n + 1][j + 1] + Aj[j + 1][n + 1]);
          Gp[n][j] = v;
        }
      }
    }

    /// The first placement r_1 < ... < r_k of itertools.product(*windows) (targets_v88.py's pl0).
    bool firstPlacement(const std::vector<std::vector<int>>& w, size_t j, int after, std::vector<int>& pl)
    {
      if (j == w.size()) return true;
      for (int r : w[j])
      {
        if (r <= after) continue;
        pl.push_back(r);
        if (firstPlacement(w, j + 1, r, pl)) return true;
        pl.pop_back();
      }
      return false;
    }
  } // namespace

  Refine::Assemble Refine::assemble(const std::vector<double>& masses_sorted, const std::vector<double>& intensities, const Ladder& ladder,
                                    const Detect& d, const RefineSettings& s)
  {
    Assemble out;
    const int NS = static_cast<int>(ladder.prefix.size()); // cleavage sites 1..L-1
    if (NS < 1 || ladder.suffix.size() != ladder.prefix.size()) return out;
    const int L = NS + 1;
    const char ptype = prefixType(ladder), stype = suffixType(ladder);

    // Lambda_b, Lambda_y: the levels that are not neutral losses with a prefix / suffix ion, ascending (column_dp_v85.py LB2/LY2)
    std::vector<double> LB, LY;
    for (const Level& l : d.levels)
    {
      if (l.neutral_loss) continue;
      if (l.n_prefix > 0) LB.push_back(l.shift);
      if (l.n_suffix > 0) LY.push_back(l.shift);
    }
    std::sort(LB.begin(), LB.end());
    std::sort(LY.begin(), LY.end());

    // lenient_core.matrices: the totals lambda_k = u + w, u in {0} u Lambda_b, w in {0} u Lambda_y, two within SAME one
    std::vector<double> UB = {0.0}, UY = {0.0};
    UB.insert(UB.end(), LB.begin(), LB.end());
    UY.insert(UY.end(), LY.begin(), LY.end());
    std::vector<double> totals;
    for (double u : UB)
    {
      for (double w : UY)
      {
        if (!isin(u + w, totals)) totals.push_back(u + w);
      }
    }

    for (double lk : totals)
    {
      // the columns D = {0, lambda_k} u Lambda_b u (lambda_k - Lambda_y), in the prototype's order
      std::vector<double> D = {0.0};
      std::vector<double> cand = {lk};
      cand.insert(cand.end(), LB.begin(), LB.end());
      for (double y : LY) cand.push_back(lk - y);
      for (double x : cand)
      {
        if (!isin(x, D)) D.push_back(x);
      }
      const int nd = static_cast<int>(D.size());
      const int i0 = 0;
      int ik = 0;
      for (int i = 0; i < nd; ++i)
      {
        if (same(D[i], lk))
        {
          ik = i;
          break;
        }
      }

      // lenient_core.cells: bi[i][n-1], yi[i][n-1] = the masses explaining the prefix / suffix ion of site n on column i
      std::vector<std::vector<int>> bi(nd, std::vector<int>(NS)), yi(nd, std::vector<int>(NS));
      std::vector<std::vector<long long>> Am(nd, std::vector<long long>(NS + 1, 0)); // amat(scale=1000): Am[i][n], n = 1..NS
      for (int i = 0; i < nd; ++i)
      {
        for (int n = 0; n < NS; ++n)
        {
          bi[i][n] = idxOf(ladder.prefix[n] + D[i], masses_sorted, s.tolerance_ppm);
          yi[i][n] = idxOf(ladder.suffix[n] + lk - D[i], masses_sorted, s.tolerance_ppm);
          Am[i][n + 1] = kIon * ((bi[i][n] >= 0 ? 1 : 0) + (yi[i][n] >= 0 ? 1 : 0));
        }
      }

      // fwd_bwd: F (prefix scores, site n included), G (suffix scores, site n excluded)
      std::vector<std::vector<long long>> F(NS + 1, std::vector<long long>(nd, kNeg)), G(NS + 1, std::vector<long long>(nd, kNeg));
      F[0][i0] = 0;
      for (int n = 1; n <= NS; ++n)
      {
        std::vector<long long> best = F[n - 1];
        if (ik != i0)
        {
          std::vector<long long> col = F[n - 1];
          col[ik] = kNeg; // no jump out of lambda_k
          for (int l = 0; l < nd; ++l)
          {
            long long other = kNeg; // the best lambda' != lambda
            for (int m = 0; m < nd; ++m)
            {
              if (m != l) other = std::max(other, col[m]);
            }
            long long jump = other > kNegHalf ? other - kJump : kNeg;
            if (l == i0) jump = kNeg; // no jump into 0
            best[l] = std::max(best[l], jump);
          }
        }
        for (int l = 0; l < nd; ++l) F[n][l] = best[l] > kNegHalf ? best[l] + Am[l][n] : kNeg;
      }
      G[NS][ik] = 0;
      for (int n = NS - 1; n >= 0; --n)
      {
        std::vector<long long> nxt(nd);
        for (int l = 0; l < nd; ++l) nxt[l] = G[n + 1][l] > kNegHalf ? G[n + 1][l] + Am[l][n + 1] : kNeg;
        std::vector<long long> best = nxt;
        if (ik != i0)
        {
          std::vector<long long> col = nxt;
          col[i0] = kNeg;
          for (int l = 0; l < nd; ++l)
          {
            long long other = kNeg;
            for (int m = 0; m < nd; ++m)
            {
              if (m != l) other = std::max(other, col[m]);
            }
            long long jump = other > kNegHalf ? other - kJump : kNeg;
            if (l == ik) jump = kNeg;
            best[l] = std::max(best[l], jump);
          }
        }
        G[n] = best;
      }
      const long long S = F[NS][ik];

      // tied_orders: every order of columns among the paths at the best score; orders[0] = the lexicographically smallest
      auto opt = [&](int n, int l) { return F[n][l] > kNegHalf && G[n][l] > kNegHalf && F[n][l] + G[n][l] == S; };
      std::vector<std::set<std::vector<int>>> next(nd), cur(nd);
      next[ik].insert({ik});
      for (int n = NS - 1; n >= 0; --n)
      {
        for (int l = 0; l < nd; ++l)
        {
          cur[l].clear();
          if (!opt(n, l)) continue;
          for (int mu = 0; mu < nd; ++mu)
          {
            if (!opt(n + 1, mu)) continue;
            if (mu == l)
            {
              if (F[n][l] + Am[mu][n + 1] + G[n + 1][mu] == S) cur[l].insert(next[mu].begin(), next[mu].end());
            }
            else if (l != ik && mu != i0 && F[n][l] - kJump + Am[mu][n + 1] + G[n + 1][mu] == S)
            {
              for (const auto& sq : next[mu])
              {
                std::vector<int> t = {l};
                t.insert(t.end(), sq.begin(), sq.end());
                cur[l].insert(std::move(t));
              }
            }
          }
        }
        std::swap(cur, next);
      }
      if (next[i0].empty()) continue; // unreachable lambda_k (never for lambda_k != 0: a jump reaches it)
      const std::vector<int> order = *next[i0].begin();

      Path p;
      p.total = lk;
      p.score = static_cast<int>(-floorDiv(-S, 1000)); // rows: (-((-S) // 1000), (-S) % 1000, ...)
      p.jumps = static_cast<int>(floorMod(-S, 1000));
      for (int c : order) p.levels.push_back(D[c]);

      // path_dp and windows over the chosen order, unscaled
      const int k = static_cast<int>(order.size()) - 1;
      std::vector<std::vector<long long>> Aj(k + 1, std::vector<long long>(NS + 1, 0)), Fp, Gp;
      for (int j = 0; j <= k; ++j)
      {
        for (int n = 1; n <= NS; ++n) Aj[j][n] = Am[order[j]][n] / kIon;
      }
      pathDp(Aj, NS, Fp, Gp);
      const long long N = Fp[NS][k];
      p.ions = static_cast<int>(N);
      std::vector<std::vector<int>> w(k);
      for (int j = 1; j <= k; ++j)
      {
        for (int r = 1; r <= NS; ++r)
        {
          if (Fp[r - 1][j - 1] > kNegHalf && Gp[r][j] > kNegHalf && Fp[r - 1][j - 1] + Aj[j][r] + Gp[r][j] == N) w[j - 1].push_back(r);
        }
        p.windows.emplace_back(w[j - 1].empty() ? 0 : w[j - 1].front(), w[j - 1].empty() ? 0 : w[j - 1].back());
      }

      // the explained ions at the first placement (targets_v88.py path(): pl0 and its labels)
      std::vector<int> pl;
      if (firstPlacement(w, 0, 0, pl))
      {
        std::vector<bool> seen(masses_sorted.size(), false);
        for (int n = 1; n <= NS; ++n)
        {
          const int i = static_cast<int>(std::count_if(pl.begin(), pl.end(), [n](int r) { return r <= n; }));
          const int c = order[i];
          if (bi[c][n - 1] >= 0)
          {
            p.ions_explained.push_back(Path::Ion{ptype, n, bi[c][n - 1], D[c], ladder.prefix[n - 1] + D[c]});
          }
          if (yi[c][n - 1] >= 0)
          {
            p.ions_explained.push_back(Path::Ion{stype, L - n, yi[c][n - 1], lk - D[c], ladder.suffix[n - 1] + lk - D[c]});
          }
        }
        for (const auto& ion : p.ions_explained)
        {
          if (seen[ion.mass_idx]) continue;
          seen[ion.mass_idx] = true;
          if (ion.mass_idx < static_cast<int>(intensities.size())) p.intensity_sum += intensities[ion.mass_idx];
        }
      }
      out.candidates.push_back(std::move(p));
    }

    // rows.sort(key=(-S, shifts)), stable over the totals' order; the summed intensity breaks a tie on both
    std::stable_sort(out.candidates.begin(), out.candidates.end(), [](const Path& x, const Path& y) {
      if (x.score != y.score) return x.score > y.score;
      if (x.jumps != y.jumps) return x.jumps < y.jumps;
      return x.intensity_sum > y.intensity_sum;
    });
    return out;
  }

  // ======================================================================================================================
  // The targets -- targets_v88.py
  // ======================================================================================================================

  std::vector<Refine::Discriminating> Refine::discriminating(const Path& a, const Path& b, const std::vector<double>& masses,
                                                             const std::vector<double>& intensities)
  {
    using Ion = Path::Ion;
    auto byMass = [](const Path& p) {
      std::map<int, std::vector<Ion>> m;
      for (const Ion& ion : p.ions_explained) m[ion.mass_idx].push_back(ion);
      // sorted(...) of the prototype's ("b35", 35) tuples: by the label string
      for (auto& kv : m)
      {
        std::stable_sort(kv.second.begin(), kv.second.end(), [](const Ion& x, const Ion& y) {
          return std::string(1, x.type) + std::to_string(x.index) < std::string(1, y.type) + std::to_string(y.index);
        });
      }
      return m;
    };
    auto sameIon = [](const Ion& x, const Ion& y) { return x.type == y.type && x.index == y.index; };
    // one ion per mass and path: the one in the series of the other path's ion, else the prefix ion
    auto one = [](const std::vector<Ion>& own, const std::vector<Ion>& other) -> std::optional<Ion> {
      if (own.empty()) return std::nullopt;
      if (own.size() < 2) return own.front();
      std::vector<Ion> k;
      for (const Ion& x : own)
      {
        if (std::any_of(other.begin(), other.end(), [&x](const Ion& y) { return y.type == x.type; })) k.push_back(x);
      }
      if (k.size() == 1) return k.front();
      for (const Ion& x : own)
      {
        if (isPrefix(x.type)) return x;
      }
      return own.front();
    };

    const auto A = byMass(a), B = byMass(b);
    std::set<int> keys;
    for (const auto& kv : A) keys.insert(kv.first);
    for (const auto& kv : B) keys.insert(kv.first);
    const std::vector<Ion> none;
    std::vector<Discriminating> out;
    for (int m : keys)
    {
      const auto ia_it = A.find(m), ib_it = B.find(m);
      const std::vector<Ion>& ia = ia_it == A.end() ? none : ia_it->second;
      const std::vector<Ion>& ib = ib_it == B.end() ? none : ib_it->second;
      const bool shared = std::any_of(ia.begin(), ia.end(), [&](const Ion& x) {
        return std::any_of(ib.begin(), ib.end(), [&](const Ion& y) { return sameIon(x, y); });
      });
      if (shared) continue;
      Discriminating t;
      t.mass_idx = m;
      t.mass = m < static_cast<int>(masses.size()) ? masses[m] : 0.0;
      t.intensity = m < static_cast<int>(intensities.size()) ? intensities[m] : 0.0;
      t.ion_a = one(ia, ib);
      t.ion_b = one(ib, ia);
      out.push_back(t);
    }
    std::stable_sort(out.begin(), out.end(), [](const Discriminating& x, const Discriminating& y) { return x.intensity > y.intensity; });
    return out;
  }

  std::vector<FragmentAnalysis::PTMSite> Refine::ptmSites(const Path& p)
  {
    std::vector<FragmentAnalysis::PTMSite> out;
    for (size_t j = 1; j < p.levels.size() && j - 1 < p.windows.size(); ++j)
    {
      const auto& w = p.windows[j - 1];
      out.push_back(FragmentAnalysis::PTMSite{(w.first + w.second) / 2, w.first, w.second, p.levels[j] - p.levels[j - 1]});
    }
    return out;
  }

  // ======================================================================================================================
  // The MS3 piece test (ADR-0049 decision 8): the pieces MS3FragmentMatcher predicts for a target under each path
  // ======================================================================================================================

  Refine::PieceVerdict Refine::pieceTest(const DeconvolvedSpectrum& ms3, const String& sequence, const Path& a, const Path& b,
                                         const Discriminating& target, double loose_ppm, double tight_ppm)
  {
    using TM = MS3FragmentMatcher::TheoreticalMass;
    PieceVerdict out;
    const int L = static_cast<int>(sequence.size());
    const std::string seq(sequence);

    struct Side
    {
      bool has = false;
      MS3FragmentMatcher::ProteoformContext ctx;
      char type = 0;
      int index = 0;
      std::vector<TM> theo, bare;
    };
    // the pieces exactly as calibrateAndScore computes them: subsequence, rebased PTM sites, MS3 ion types
    auto predict = [&](const Path& p, const std::optional<Path::Ion>& ion) {
      Side sd;
      if (!ion) return sd;
      sd.ctx.region_start = 0;
      sd.ctx.region_end = L;
      sd.ctx.ptm_sites = ptmSites(p);
      sd.type = ion->type;
      sd.index = ion->index;
      const std::string sub = MS3FragmentMatcher::extractSubsequence(seq, sd.ctx, sd.type, sd.index);
      if (sub.empty()) return sd;
      const int start = isPrefix(sd.type) ? 1 : L - sd.index + 1;
      const auto rebased = MS3FragmentMatcher::rebasePTMSites(sd.ctx.ptm_sites, start, sd.index);
      const auto types = MS3FragmentMatcher::getMS3IonTypes(sd.type);
      sd.theo = MS3FragmentMatcher::computeTheoreticalMasses(sub, types, rebased);
      sd.bare = MS3FragmentMatcher::computeTheoreticalMasses(sub, types, {}); // the bare backbone of every piece
      sd.has = !sd.theo.empty();
      return sd;
    };
    const Side sa = predict(a, target.ion_a), sb = predict(b, target.ion_b);
    auto within = [tight_ppm](double x, double y) { return std::fabs(x - y) <= std::max(x, y) * tight_ppm * 1e-6; };
    auto nearAny = [&](double x, const std::vector<TM>& set) {
      return std::any_of(set.begin(), set.end(), [&](const TM& t) { return within(x, t.mass); });
    };

    // the shared pieces: predicted masses (both paths) within tight_ppm of each other, one entry per coincident piece
    // (both paths' copies of one piece would let one observed mass count twice, in shared_seen and in the calibration)
    std::vector<TM> both;
    for (const TM& t : sa.theo)
    {
      if (nearAny(t.mass, sb.theo)) both.push_back(t);
    }
    for (const TM& t : sb.theo)
    {
      if (nearAny(t.mass, sa.theo)) both.push_back(t);
    }
    std::stable_sort(both.begin(), both.end(), [](const TM& x, const TM& y) { return x.mass < y.mass; });
    std::vector<TM> shared;
    for (const TM& t : both)
    {
      if (shared.empty() || !within(t.mass, shared.back().mass)) shared.push_back(t);
    }

    // calibration: a loose pass against the shared pieces only, the median signed ppm
    std::vector<double> ppm;
    MS3FragmentMatcher::matchSpectrum(ms3, shared, loose_ppm, &ppm);
    if (!ppm.empty())
    {
      std::sort(ppm.begin(), ppm.end());
      const size_t mid = ppm.size() / 2;
      const double median = ppm.size() % 2 == 0 ? (ppm[mid - 1] + ppm[mid]) / 2.0 : ppm[mid];
      out.correction = 1.0 / (1.0 + median * 1e-6);
    }
    DeconvolvedSpectrum corrected(0);
    for (Size i = 0; i < ms3.size(); ++i)
    {
      PeakGroup pg = ms3[i];
      pg.setMonoisotopicMass(pg.getMonoMass() * out.correction);
      corrected.push_back(pg);
    }
    out.shared_seen = MS3FragmentMatcher::matchSpectrum(corrected, shared, tight_ppm);

    // per path: the matched pieces that are discriminating (no piece of the other path within tight_ppm) and carry a shift
    auto fits = [&](const Side& own, const Side& other) {
      if (!own.has) return false;
      std::vector<FragmentAnalysis::ProteoformMatch> detail;
      // loose = 0: calibrateAndScore's own pass 1 then matches nothing and leaves its correction at 1.0, so the only
      // calibration is the shared-piece scaling above -- a path is never calibrated on the pieces whose fit is tested
      MS3FragmentMatcher::calibrateAndScore({&corrected}, seq, own.ctx, own.type, own.index, 0.0, tight_ppm, &detail);
      if (detail.empty()) return false;
      int count = 0;
      for (const auto& fm : detail[0].fragments)
      {
        const auto t = std::find_if(own.theo.begin(), own.theo.end(), [&fm](const TM& x) {
          return x.ion_type == fm.ion_type && x.position == fm.ion_index && x.includes_ptm == fm.includes_ptm;
        });
        const auto bare = std::find_if(own.bare.begin(), own.bare.end(), [&fm](const TM& x) {
          return x.ion_type == fm.ion_type && x.position == fm.ion_index;
        });
        if (t == own.theo.end() || bare == own.bare.end()) continue;
        const bool discriminating = !nearAny(t->mass, other.theo);
        const bool shift_carrying = std::fabs(t->mass - bare->mass) > t->mass * tight_ppm * 1e-6;
        count += discriminating && shift_carrying;
      }
      return count >= 1;
    };
    out.fits_a = fits(sa, sb);
    out.fits_b = fits(sb, sa);
    return out;
  }

  // ======================================================================================================================
  // The recount -- ms3_recount_v79.py
  // ======================================================================================================================

  Refine::Recount Refine::recount(const Path& a, const Path& b, const std::vector<Discriminating>& targets,
                                  const std::vector<std::pair<bool, bool>>& confirmed)
  {
    // shared = the masses a path explains minus the targets it explains (the prototype's N - len(only[k]), one ion per mass)
    auto masses = [](const Path& p) {
      std::set<int> m;
      for (const auto& ion : p.ions_explained) m.insert(ion.mass_idx);
      return static_cast<int>(m.size());
    };
    int targets_a = 0, targets_b = 0;
    Recount r;
    for (size_t i = 0; i < targets.size(); ++i)
    {
      targets_a += targets[i].ion_a.has_value();
      targets_b += targets[i].ion_b.has_value();
      if (i < confirmed.size())
      {
        r.confirmed_a += confirmed[i].first && targets[i].ion_a.has_value();
        r.confirmed_b += confirmed[i].second && targets[i].ion_b.has_value();
      }
    }
    const int shared_a = masses(a) - targets_a, shared_b = masses(b) - targets_b;
    if (shared_a != shared_b)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Refine::recount: the two paths share " + std::to_string(shared_a) + " and " + std::to_string(shared_b) +
                                          " masses; the targets are not theirs");
    }
    r.score_a = shared_a + r.confirmed_a - a.jumps;
    r.score_b = shared_b + r.confirmed_b - b.jumps;
    r.a_wins = r.score_a > r.score_b || (r.score_a == r.score_b && a.jumps <= b.jumps);
    return r;
  }
} // namespace OpenMS
