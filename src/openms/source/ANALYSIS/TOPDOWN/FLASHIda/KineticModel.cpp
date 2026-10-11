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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/KineticModel.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace OpenMS::KineticModel
{
  namespace
  {
    /// curve_families.U clips the exponent to +-60 so exp(exp(u)) cannot overflow on a wild grid point.
    inline double clampU(double u) { return std::min(std::max(u, -60.0), 60.0); }

    /// numpy.linspace(a, b, n): a + i * (b - a) / (n - 1), the last element set to b exactly.
    std::vector<double> linspace(double a, double b, int n)
    {
      std::vector<double> v(static_cast<size_t>(n));
      const double step = (b - a) / (n - 1);
      for (int i = 0; i < n; ++i) v[static_cast<size_t>(i)] = a + i * step;
      v.back() = b;
      return v;
    }

    /// numpy.geomspace(a, b, n): 10 ** linspace(log10 a, log10 b, n), the endpoints set exactly.
    std::vector<double> geomspace(double a, double b, int n)
    {
      std::vector<double> v = linspace(std::log10(a), std::log10(b), n);
      for (double& x : v) x = std::pow(10.0, x);
      v.front() = a;
      v.back() = b;
      return v;
    }

    struct GridBest
    {
      double A = -1.0, c = -1.0, w = -1.0, sse = std::numeric_limits<double>::infinity();
    };

    /// curve_families.fit's best_on for the Gompertz family: for every (c, w) of the lattice the plateau A is
    /// the closed-form least-squares scale, clipped to [0, cap]; the lattice point with the smallest residual
    /// wins, the FIRST such point on a tie (np.argmin over the c-major mesh). Verbatim behaviour, so the
    /// engine's fit is the study's fit, tie-breaks included (ADR-0047 decision 3, as amended at
    /// implementation: see KineticModel.h).
    GridBest bestOn(const std::vector<SharePoint>& pts, double cap, const std::vector<double>& cs, const std::vector<double>& ws)
    {
      GridBest best;
      std::vector<double> g(pts.size());
      for (double c : cs)
      {
        for (double w : ws)
        {
          double gp = 0.0, gg = 0.0;
          for (size_t i = 0; i < pts.size(); ++i)
          {
            g[i] = std::exp(-std::exp(clampU((pts[i].ce - c) / w)));
            gp += g[i] * pts[i].share;
            gg += g[i] * g[i];
          }
          double A = gp / std::max(gg, 1e-30);
          A = std::min(std::max(A, 0.0), cap);
          double sse = 0.0;
          for (size_t i = 0; i < pts.size(); ++i)
          {
            const double r = pts[i].share - A * g[i];
            sse += r * r;
          }
          if (sse < best.sse)
          {
            best.A = A;
            best.c = c;
            best.w = w;
            best.sse = sse;
          }
        }
      }
      return best;
    }
  } // namespace

  double GompertzFit::share(double ce) const
  {
    if (!ok) return 0.0;
    return p0 * std::exp(-std::exp(clampU((ce - ce_max_s) / width)));
  }

  double GompertzFit::rate(double ce) const
  {
    if (!ok) return 0.0;
    const double e = std::exp(clampU((ce - ce_max_s) / width));
    return std::max(0.0, p0 * e * std::exp(-e) / width);   // curve_families "Gompertz": -dg/dx
  }

  double precursorShare(const double* mzs, const double* ints, int n, double mz_lo, double mz_hi)
  {
    if (mzs == nullptr || ints == nullptr || n <= 0) return -1.0;
    double total = 0.0, window = 0.0;
    for (int i = 0; i < n; ++i)
    {
      total += ints[i];
      if (mzs[i] >= mz_lo && mzs[i] <= mz_hi) window += ints[i];
    }
    return total > 0.0 ? window / total : -1.0;
  }

  GompertzFit fitGompertz(const std::vector<SharePoint>& points)
  {
    GompertzFit out;
    std::vector<SharePoint> pts;
    for (const auto& p : points)
      if (std::isfinite(p.share) && p.share >= 0.0 && std::isfinite(p.ce)) pts.push_back(p);
    std::sort(pts.begin(), pts.end(), [](const SharePoint& a, const SharePoint& b) { return a.ce < b.ce; });
    out.n_points = static_cast<int>(pts.size());
    if (out.n_points < kMinPoints) return out;

    double hi = 0.0;
    for (const auto& p : pts) hi = std::max(hi, p.share);
    if (hi <= 0.0) return out;
    const double cap = kPlateauFactor * hi;   // curve_families.plateau_bound

    // The study's two passes (curve_families.fit, family "Gompertz"): the coarse lattice -- inflection
    // 0..60 in steps of 0.5, 41 log-spaced widths over [kWidthMin, kWidthMax] -- then one refinement
    // around the winner: 51 inflections within +-kRefineC, 41 widths within a factor kRefineW either way,
    // the width floor at half the coarse floor exactly as the study has it.
    std::vector<double> cs;
    for (int i = 0; i <= 120; ++i) cs.push_back(0.5 * i);
    const GridBest coarse = bestOn(pts, cap, cs, geomspace(kWidthMin, kWidthMax, 41));
    if (!std::isfinite(coarse.sse)) return out;
    const GridBest fine = bestOn(pts, cap,
                                 linspace(std::max(coarse.c - kRefineC, cs.front()), coarse.c + kRefineC, 51),
                                 geomspace(std::max(coarse.w / kRefineW, kWidthMin / 2.0), coarse.w * kRefineW, 41));
    if (!std::isfinite(fine.sse)) return out;

    out.ok = true;
    out.p0 = fine.A;
    out.ce_max_s = fine.c;
    out.width = fine.w;
    out.rss = fine.sse;
    return out;
  }

  std::vector<double> fragmentPool(const std::vector<double>& rate, double k2)
  {
    std::vector<double> B(rate.size(), 0.0);
    const double a = std::exp(-k2), g = -std::expm1(-k2) / k2;   // pipeline.direct, verbatim
    double b = 0.0;
    for (size_t k = 0; k < rate.size(); ++k)
    {
      b = a * b + rate[k] * g;
      B[k] = b;
    }
    return B;
  }

  ModeledMaximum modeledMaximum(const GompertzFit& fit, double k2, int ce_max)
  {
    ModeledMaximum m;
    if (!fit.ok || ce_max < 1 || !(k2 > 0.0)) return m;
    std::vector<double> s(static_cast<size_t>(ce_max));
    for (int ce = 1; ce <= ce_max; ++ce) s[static_cast<size_t>(ce - 1)] = fit.rate(ce);
    const std::vector<double> B = fragmentPool(s, k2);
    int best = 0;
    for (int k = 1; k < ce_max; ++k)
      if (B[static_cast<size_t>(k)] > B[static_cast<size_t>(best)]) best = k;   // strict: the lowest CE on a tie
    m.ce = best + 1;
    m.b = B[static_cast<size_t>(best)];
    m.at_lower_edge = (m.ce == 1);
    m.at_upper_edge = (m.ce == ce_max);
    return m;
  }
} // namespace OpenMS::KineticModel
