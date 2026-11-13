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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ChanceMatch.h>

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace OpenMS
{
  namespace ChanceMatch
  {
    double log10PoissonTail(int k, double lambda)
    {
      if (k <= 0) return 0.0;
      if (!(lambda > 0.0)) return -std::numeric_limits<double>::infinity();
      const double ln10 = std::log(10.0);
      const double ln_lambda = std::log(lambda);
      auto lnTerm = [&](int j) { return -lambda + j * ln_lambda - std::lgamma(static_cast<double>(j) + 1.0); };

      if (static_cast<double>(k) <= lambda)
      {
        // the tail holds most of the mass: 1 minus the head (k terms), each term from its own log
        double head = 0.0;
        for (int j = 0; j < k; ++j) head += std::exp(lnTerm(j));
        const double tail = 1.0 - head;
        if (tail > 0.0) return std::log10(tail);
        // the head rounded to 1: fall through to the series, which is exact for any k
      }

      // the upward series from term k, factored out in log space: terms fall by lambda / j once j > lambda
      double sum = 1.0, ratio = 1.0;
      for (int j = k + 1; j < k + 100000; ++j)
      {
        ratio *= lambda / static_cast<double>(j);
        sum += ratio;
        if (j > lambda && ratio < 1e-17 * sum) break;
      }
      return (lnTerm(k) + std::log(sum)) / ln10;
    }

    double chanceLambda(std::size_t n_masses, const std::vector<double>& theoretical, double tol_ppm, double chain_mass)
    {
      if (!(chain_mass > 0.0)) return 0.0;
      double share = 0.0;
      for (double t : theoretical)
      {
        if (t > 0.0) share += 2.0 * tol_ppm * 1e-6 * t;
      }
      return static_cast<double>(n_masses) * share / chain_mass;
    }

    double chainMass(const std::string& sequence, const std::vector<double>& shifts, const FixedModifications& fixed_mods)
    {
      static const double water = EmpiricalFormula("H2O").getMonoWeight();
      double mass = water;
      for (char c : sequence)
      {
        const char aa = (c == 'I') ? 'L' : c;   // as computePTMAdjustedFragmentMasses counts it
        if (aa == 'X' || ResidueDB::getInstance()->getResidue(aa) == nullptr) continue;
        mass += fixed_mods.mass(aa);
      }
      for (double s : shifts) mass += s;
      return mass;
    }

    namespace
    {
      // mantissa and exponent of 10^x, the mantissa rounded to @p decimals places with its carry into the exponent
      std::string scientificFromLog10(double x, int decimals)
      {
        int exponent = static_cast<int>(std::floor(x));
        double mantissa = std::pow(10.0, x - exponent);
        char m[32];
        std::snprintf(m, sizeof(m), "%.*f", decimals, mantissa);
        if (std::atof(m) >= 10.0)
        {
          ++exponent;
          std::snprintf(m, sizeof(m), "%.*f", decimals, 1.0);
        }
        char out[64];
        std::snprintf(out, sizeof(out), "%se%+03d", m, exponent);
        return out;
      }
    } // namespace

    std::string formatLog10Scientific(const std::optional<double>& log10_value)
    {
      if (!log10_value.has_value()) return "-1";
      if (std::isinf(*log10_value) && *log10_value < 0) return "0";
      return scientificFromLog10(*log10_value, 2);
    }

    std::string formatLog10Short(const std::optional<double>& log10_value)
    {
      if (!log10_value.has_value()) return "-";
      if (std::isinf(*log10_value) && *log10_value < 0) return "0";
      return scientificFromLog10(*log10_value, 0);
    }
  } // namespace ChanceMatch
} // namespace OpenMS
