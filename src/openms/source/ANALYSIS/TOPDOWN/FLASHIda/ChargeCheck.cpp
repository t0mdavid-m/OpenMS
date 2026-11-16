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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ChargeCheck.h>

#include <unsupported/Eigen/FFT>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <vector>

namespace OpenMS
{
  namespace ChargeCheck
  {
    namespace
    {
      // The fixed parameters of the spec (4.1). They are the prototype's (flashida-paper w08 / w14d), whose results the replay
      // fixtures reproduce; changing one changes what FLASHIda_ChargeCheck_test's replay sections pin.
      constexpr double GRID_STEP = 0.0002;      // Th
      constexpr double REBUILD_MARGIN = 0.05;   // Th: centroids this far outside the window still add their tails
      constexpr double SUPPORT_FWHM = 2.0;      // each Gaussian spans +- 2 FWHM
      constexpr double FWHM_PER_SIGMA = 2.3548;
      constexpr double MAX_LAG_SPAN = 0.6;      // Th
      constexpr double LAG_SPAN_FRACTION = 0.9; // of the window width
      constexpr double HALF_CHARGE_RATIO = 0.8;
      constexpr int REFINE_POINTS = 201;
      constexpr size_t MIN_GRID_POINTS = 10;

      double fwhmAt(double mz, double r200)
      {
        return mz / (r200 * std::sqrt(200.0 / mz));
      }

      // numpy.searchsorted(grid, v) on the uniform grid lo + i * step, i < n: the first index whose point is >= v.
      long firstAtOrAbove(double v, double lo, size_t n)
      {
        const double p = std::ceil((v - lo) / GRID_STEP);
        if (p <= 0.0) return 0;
        if (p >= static_cast<double>(n)) return static_cast<long>(n);
        return static_cast<long>(p);
      }

      // numpy.median: the mean of the two middle values for an even count.
      double median(std::vector<double> v)
      {
        const size_t n = v.size();
        const size_t mid = n / 2;
        std::nth_element(v.begin(), v.begin() + mid, v.end());
        const double upper = v[mid];
        if (n % 2 == 1) return upper;
        const double lower = *std::max_element(v.begin(), v.begin() + mid);
        return 0.5 * (lower + upper);
      }

      // The comb sum at spacing s: the mean autocorrelation at k * s, k = 1..floor(span / s), linearly interpolated (numpy.interp,
      // clamped at the last lag). One span for every spacing, so twice the charge cannot score on every second peak (spec 4.4).
      double combAt(const std::vector<double>& ac, double s, double span)
      {
        const int K = static_cast<int>(span / s);
        if (K < 1) return 0.0;
        const size_t n = ac.size();
        double sum = 0.0;
        for (int k = 1; k <= K; ++k)
        {
          const double p = k * s / GRID_STEP;
          const size_t i = static_cast<size_t>(p);
          if (i + 1 >= n)
          {
            sum += ac[n - 1];
            continue;
          }
          const double f = p - static_cast<double>(i);
          sum += ac[i] * (1.0 - f) + ac[i + 1] * f;
        }
        return sum / K;
      }
    } // namespace

    double resolvability(double mz, int z, double r200)
    {
      if (z <= 0 || mz <= 0.0 || r200 <= 0.0) return 0.0;
      return (ISO_SPACING / z) / fwhmAt(mz, r200);
    }

    Measurement comb(const MSSpectrum& survey, double lo, double hi, double r200, int z_min, int z_max)
    {
      Measurement m;
      m.status = Status::Empty;
      if (!(hi > lo) || r200 <= 0.0) return m;
      z_min = std::max(z_min, 1);
      if (z_max < z_min) return m;

      // 1. The window's profile, rebuilt from its centroids on the grid lo + i * step, i < n (numpy.arange(lo, hi, step)).
      const size_t n = static_cast<size_t>(std::ceil((hi - lo) / GRID_STEP));
      if (n < MIN_GRID_POINTS) return m;
      std::vector<double> y(n, 0.0);
      bool any = false;
      for (auto it = survey.MZBegin(lo - REBUILD_MARGIN); it != survey.end() && it->getMZ() <= hi + REBUILD_MARGIN; ++it)
      {
        const double mz = it->getMZ();
        const double a = it->getIntensity();
        if (a <= 0.0) continue;
        const double fwhm = fwhmAt(mz, r200);
        const double sigma = fwhm / FWHM_PER_SIGMA;
        const long a0 = firstAtOrAbove(mz - SUPPORT_FWHM * fwhm, lo, n);
        const long b0 = firstAtOrAbove(mz + SUPPORT_FWHM * fwhm, lo, n);
        for (long i = a0; i < b0; ++i)
        {
          const double d = (lo + static_cast<double>(i) * GRID_STEP - mz) / sigma;
          y[i] += a * std::exp(-0.5 * d * d);
          any = true;
        }
      }
      if (!any) return m;

      // 2. Baseline: subtract the median, clamp at 0.
      const double base = median(y);
      for (double& v : y) v = std::max(0.0, v - base);

      // 3. Autocorrelation by FFT, zero-padded to a power of two >= 2n: linear, not circular, at every lag below n.
      size_t nfft = 1;
      while (nfft < 2 * n) nfft <<= 1;
      std::vector<std::complex<double>> x(nfft, std::complex<double>(0.0, 0.0));
      for (size_t i = 0; i < n; ++i) x[i] = std::complex<double>(y[i], 0.0);
      Eigen::FFT<double> fft;
      std::vector<std::complex<double>> spectrum;
      fft.fwd(spectrum, x);
      for (auto& c : spectrum) c = std::complex<double>(std::norm(c), 0.0);
      std::vector<std::complex<double>> r;
      fft.inv(r, spectrum);
      const double ac0 = r[0].real();
      if (!(ac0 > 0.0)) return m;
      std::vector<double> ac(n);
      for (size_t i = 0; i < n; ++i) ac[i] = r[i].real() / ac0;

      // 4. The comb over the candidate charges. Of those within HALF_CHARGE_RATIO of the best, the highest charge: half the
      //    charge (twice the spacing) also matches every second peak and scores high.
      const double span = std::min((hi - lo) * LAG_SPAN_FRACTION, MAX_LAG_SPAN);
      std::vector<double> S;
      S.reserve(static_cast<size_t>(z_max - z_min + 1));
      double best = -1.0;
      for (int z = z_min; z <= z_max; ++z)
      {
        S.push_back(combAt(ac, ISO_SPACING / z, span));
        best = std::max(best, S.back());
      }
      int z_meas = z_min;
      for (int z = z_max; z >= z_min; --z)
      {
        if (S[static_cast<size_t>(z - z_min)] >= HALF_CHARGE_RATIO * best)
        {
          z_meas = z;
          break;
        }
      }

      // 5. Refine the spacing within +- 0.5 charge; the first maximum, as numpy.argmax.
      const double s_lo = ISO_SPACING / (z_meas + 0.5);
      const double s_hi = z_meas > 1 ? ISO_SPACING / (z_meas - 0.5) : ISO_SPACING * 1.5;
      double s_best = s_lo;
      double c_best = -1.0;
      for (int k = 0; k < REFINE_POINTS; ++k)
      {
        const double s = s_lo + (s_hi - s_lo) * k / (REFINE_POINTS - 1);
        const double c = combAt(ac, s, span);
        if (c > c_best)
        {
          c_best = c;
          s_best = s;
        }
      }

      m.status = Status::NotChecked;
      m.z_meas = z_meas;
      m.z_cont = ISO_SPACING / s_best;
      m.score = S[static_cast<size_t>(z_meas - z_min)];
      return m;
    }

    Measurement measure(const MSSpectrum& survey, double lo, double hi, double mz, int z_cmd, double r200, int z_min, int z_max,
                        const Settings& settings)
    {
      // The comb runs whatever the gates decide, so the log carries the measured charge of every checked window.
      Measurement m = comb(survey, lo, hi, r200, z_min, z_max);
      m.resolvability = resolvability(mz, z_cmd, r200);
      if (m.resolvability < settings.min_resolvability)
      {
        m.status = Status::UnverifiedResolution;
        return m;
      }
      if (m.status == Status::Empty || m.z_meas == 0)
      {
        m.status = Status::Empty;
        return m;
      }
      if (m.score < settings.min_score)
      {
        m.status = Status::UnverifiedLowScore;
        return m;
      }
      const int dz = std::abs(m.z_meas - z_cmd);
      if (dz == 0) m.status = Status::Confirmed;
      else if (dz <= settings.max_dz) m.status = Status::Changed;
      else m.status = Status::Inconsistent;
      return m;
    }

    const char* toString(Status status)
    {
      switch (status)
      {
        case Status::NotChecked: return "not_checked";
        case Status::Confirmed: return "confirmed";
        case Status::Changed: return "changed";
        case Status::UnverifiedResolution: return "unverified_resolution";
        case Status::UnverifiedLowScore: return "unverified_low_score";
        case Status::UnverifiedNoSpecies: return "unverified_no_species";
        case Status::Inconsistent: return "inconsistent";
        case Status::Empty: return "empty";
        case Status::Swapped: return "swapped";
        case Status::Remeasured: return "remeasured";
      }
      return "not_checked";
    }
  } // namespace ChargeCheck
} // namespace OpenMS
