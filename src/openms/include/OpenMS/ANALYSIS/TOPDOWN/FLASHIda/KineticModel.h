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

#pragma once

#include <OpenMS/config.h>

#include <vector>

/**
  @brief The Kinetic Model of fragmentation along collision energy (ADR-0047).

  Pure functions, no engine state. The progress report of 2026-09-30 ("Precursor Tuned Fragmentation")
  models fragments as a pool fed at the FRAGMENTATION RATE s(CE) = -dP/dCE of the PRECURSOR SHARE
  P(CE) and drained at the OVERFRAGMENTATION RATE k2:

      B(CE) = beta * integral_0^CE s(u) exp(-k2 (CE - u)) du,     CE* = argmax B(CE)

  beta only scales B and cancels from the argmax. P is read from a few pre-scans and carried by a
  GOMPERTZ CURVE P(CE) = P0 exp(-exp((CE - c)/w)) -- the one assumption being that the relative
  fragmentation rate s/P rises exponentially along CE -- so six pre-scans suffice. The study
  (108 HCD sweeps of 11 proteins) reaches a mean regret of 0.149 against 0.357 for the best fixed
  energy with that six-point fit; this file reproduces its fit, its pool recursion and its argmax.

  The reference implementation is the study's own pipeline, vendored at
  FlashIDA/test-scripts/kinetic/{pipeline.py,curve_families.py}; KineticModel_test replays all 108
  sweeps x 10 grid offsets against the tables it generated.
*/
namespace OpenMS::KineticModel
{
  /// The study's bounded-Gompertz fit box and rules. Constants of the validated algorithm, not configuration.
  constexpr double kWidthMin      = 0.25;  ///< NCE; the study's width grid floor
  constexpr double kWidthMax      = 25.0;  ///< NCE; the study's width grid ceiling
  constexpr double kPlateauFactor = 1.1;   ///< P0 <= kPlateauFactor x highest acquired share (the "bounded" arm)
  constexpr int    kMinPoints     = 3;     ///< fewer finite shares than this -> no fit
  constexpr double kSeedWidth     = 4.2;   ///< NCE; median width of the 108 study sweeps, the seed when no 0.85 crossing exists

  /// One pre-scan's reading: its collision energy and its precursor share. A share < 0 or non-finite
  /// means "no measurement" and is dropped from any fit without counting toward kMinPoints.
  struct SharePoint
  {
    double ce = 0.0;
    double share = -1.0;
  };

  /// A fitted Gompertz curve P(ce) = p0 * exp(-exp((ce - ce_max_s) / width)). `ok == false` means no fit:
  /// too few finite points, or the solver did not converge. Every field keeps its -1 sentinel then.
  struct OPENMS_DLLAPI GompertzFit
  {
    bool   ok = false;
    double p0 = -1.0;        ///< initial precursor share, the plateau
    double ce_max_s = -1.0;  ///< the inflection point: the CE of the highest fragmentation rate, where P = p0/e
    double width = -1.0;     ///< the CE step over which the relative fragmentation rate s/P grows by e
    double rss = -1.0;       ///< residual sum of squares in share units over the fitted points
    int    n_points = 0;     ///< finite points the fit was given (also set when ok == false)

    double share(double ce) const;  ///< P(ce); 0 when !ok
    double rate(double ce) const;   ///< s(ce) = -dP/dce >= 0; 0 when !ok
  };

  /// The energy the model chooses. `ce == 0` means none (no fit, no k2, no range).
  struct ModeledMaximum
  {
    int    ce = 0;                 ///< integer collision energy in [1, ce_max]; 0 = none
    double b = -1.0;               ///< the pool value at ce (up to the factor beta)
    bool   at_lower_edge = false;  ///< ce == 1: the grid started too high
    bool   at_upper_edge = false;  ///< ce == ce_max: the curve's fall lies beyond the range
  };

  /// Precursor share of one scan: sum of @p ints whose mz lies in [mz_lo, mz_hi] over the sum of all
  /// @p ints. Both sums carry the scan's ion supply, so the supply cancels and no reference scan is
  /// needed. Returns -1.0 when n <= 0 or the total is <= 0 (no measurement, not a share of 0).
  OPENMS_DLLAPI double precursorShare(const double* mzs, const double* ints, int n, double mz_lo, double mz_hi);

  /// Bounded Gompertz fit by Levenberg-Marquardt (ADR-0047 decision 3): least squares in share units on
  /// the finite points, with the study's box imposed exactly by reparametrization -- width through a log
  /// box onto (kWidthMin, kWidthMax), p0 through a sigmoid onto (0, kPlateauFactor x highest share), the
  /// inflection free -- seeded from the landmark read of the points (p0 = highest share, inflection at
  /// the interpolated p0/e crossing, width from the 0.85 p0 crossing). Fewer than kMinPoints finite
  /// points, or non-convergence, returns a fit with ok == false.
  OPENMS_DLLAPI GompertzFit fitGompertz(const std::vector<SharePoint>& points);

  /// The study's pool recursion (pipeline.direct): B(k) = e^{-k2} B(k-1) + rate[k] (1 - e^{-k2}) / k2,
  /// B before the first energy = 0. Returns one value per entry of @p rate.
  OPENMS_DLLAPI std::vector<double> fragmentPool(const std::vector<double>& rate, double k2);

  /// Rate on the integer energies 1..ce_max from @p fit, the pool with @p k2, and the integer argmax
  /// (the lowest CE on a tie). The maximum may lie below the first acquired pre-scan; it is never clamped
  /// to the acquired range, only flagged when it sits at either end of [1, ce_max].
  OPENMS_DLLAPI ModeledMaximum modeledMaximum(const GompertzFit& fit, double k2, int ce_max);
}
