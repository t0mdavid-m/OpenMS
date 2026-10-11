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

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

#include <algorithm>
#include <cmath>

namespace OpenMS::KineticModel
{
  namespace
  {
    inline double sigmoid(double t) { return 1.0 / (1.0 + std::exp(-t)); }
    inline double logit(double p)
    {
      p = std::min(std::max(p, 1e-6), 1.0 - 1e-6);
      return std::log(p / (1.0 - p));
    }
    /// curve_families.U clips the exponent to +-60 so exp(exp(u)) cannot overflow on a wild step.
    inline double clampU(double u) { return std::min(std::max(u, -60.0), 60.0); }

    /// The coordinates LM walks: x(0) -> p0 through cap * sigmoid, x(1) = c free, x(2) -> w through the
    /// log box. Each map is a smooth bijection of the real line onto the open box, so every point the
    /// study's grid could return is reachable and nothing outside it is; the objective's value at a point
    /// is unchanged, so its minimum over the box and over (x0, x1, x2) coincide (ADR-0047 decision 3).
    struct Decoded { double p0, c, w; };
    inline Decoded decode(const Eigen::VectorXd& x, double cap)
    {
      return {cap * sigmoid(x(0)), x(1), kWidthMin * std::pow(kWidthMax / kWidthMin, sigmoid(x(2)))};
    }
    inline Eigen::VectorXd encode(double p0, double c, double w, double cap)
    {
      Eigen::VectorXd x(3);
      x(0) = logit(p0 / cap);
      x(1) = c;
      x(2) = logit(std::log(w / kWidthMin) / std::log(kWidthMax / kWidthMin));
      return x;
    }

    // The functor shape Eigen's NumericalDiff / LevenbergMarquardt expect (same idiom as
    // GumbelMaxLikelihoodFitter.cpp).
    template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
    struct Functor
    {
      typedef _Scalar Scalar;
      enum { InputsAtCompileTime = NX, ValuesAtCompileTime = NY };
      typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
      typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
      typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;
      int m_inputs, m_values;
      Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}
      int inputs() const { return m_inputs; }
      int values() const { return m_values; }
    };

    struct GompertzFunctor : Functor<double>
    {
      const std::vector<SharePoint>& pts;
      double cap;
      GompertzFunctor(const std::vector<SharePoint>& p, double c) : Functor<double>(3, static_cast<int>(p.size())), pts(p), cap(c) {}
      int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& f) const
      {
        const Decoded d = decode(x, cap);
        for (int i = 0; i < values(); ++i)
          f(i) = pts[i].share - d.p0 * std::exp(-std::exp(clampU((pts[i].ce - d.c) / d.w)));
        return 0;
      }
    };

    /// The landmark read (the report's "read without a fit"): p0 = highest share; c = the first crossing
    /// of p0/e after the maximum, linearly interpolated; w from the 0.85 p0 crossing,
    /// w = (ce_85 - c) / ln(-ln 0.85) (ln(-ln 0.85) = -1.817, and ce_85 < c, so w > 0). Fallbacks:
    /// c = the last energy, w = kSeedWidth.
    Decoded landmarkSeed(const std::vector<SharePoint>& pts)
    {
      int imax = 0;
      for (int i = 1; i < static_cast<int>(pts.size()); ++i)
        if (pts[i].share > pts[imax].share) imax = i;
      const double p0 = pts[imax].share;
      auto crossing = [&](double level, bool& found) {
        for (int i = imax; i + 1 < static_cast<int>(pts.size()); ++i)
        {
          if (pts[i].share >= level && pts[i + 1].share < level)
          {
            found = true;
            const double t = (pts[i].share - level) / (pts[i].share - pts[i + 1].share);
            return pts[i].ce + t * (pts[i + 1].ce - pts[i].ce);
          }
        }
        found = false;
        return pts.back().ce;
      };
      bool fc = false, f85 = false;
      const double c = crossing(p0 / std::exp(1.0), fc);
      const double c85 = crossing(0.85 * p0, f85);
      double w = kSeedWidth;
      if (fc && f85) w = (c85 - c) / std::log(-std::log(0.85));
      // Keep the seed strictly inside the box: logit() is undefined at its walls.
      w = std::min(std::max(w, kWidthMin * 1.2), kWidthMax / 1.2);
      return {std::max(p0, 1e-3), c, w};
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
    const double cap = kPlateauFactor * hi;

    const Decoded seed = landmarkSeed(pts);
    Eigen::VectorXd x = encode(std::min(seed.p0, 0.99 * cap), seed.c, seed.w, cap);

    GompertzFunctor functor(pts, cap);
    Eigen::NumericalDiff<GompertzFunctor> nd(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<GompertzFunctor>, double> lm(nd);
    lm.parameters.maxfev = 2000;
    const Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(x);
    // GaussFitter's rule: every status but these two is a good termination state.
    if (status == Eigen::LevenbergMarquardtSpace::ImproperInputParameters
        || status == Eigen::LevenbergMarquardtSpace::TooManyFunctionEvaluation)
      return out;

    const Decoded d = decode(x, cap);
    Eigen::VectorXd f(pts.size());
    functor(x, f);
    const double rss = f.squaredNorm();
    if (!std::isfinite(d.p0) || !std::isfinite(d.c) || !std::isfinite(d.w) || !std::isfinite(rss)) return out;

    out.ok = true;
    out.p0 = d.p0;
    out.ce_max_s = d.c;
    out.width = d.w;
    out.rss = rss;
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
