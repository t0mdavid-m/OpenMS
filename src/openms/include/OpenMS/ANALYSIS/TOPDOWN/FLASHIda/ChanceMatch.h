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

#include <OpenMS/ANALYSIS/TOPDOWN/FixedModifications.h>
#include <OpenMS/config.h>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
   * @brief The chance-match p-value of a fragment match (ADR-0075; CONTEXT.md "Chance-match p-value").
   *
   * Were a spectrum's N masses unrelated to a chain of mass M, each would fall within tolerance of one of its
   * theoretical fragment masses t with probability q = sum 2 * tol * t / M, so the chance matches are Poisson with
   * lambda = N * q. The p-value of k matches is P(X >= k). Every value here is a log10: a long ladder with many matches
   * reaches p far below the smallest double. Pure functions; reported only, never decided on.
   */
  namespace ChanceMatch
  {
    /// log10 P(X >= k) for X ~ Poisson(@p lambda). 0 for k <= 0; -infinity for lambda <= 0 with k > 0.
    OPENMS_DLLAPI double log10PoissonTail(int k, double lambda);

    /// lambda = @p n_masses * sum over @p theoretical (masses > 0) of 2 * @p tol_ppm * 1e-6 * t / @p chain_mass; 0 when the
    /// chain mass is not positive.
    OPENMS_DLLAPI double chanceLambda(std::size_t n_masses, const std::vector<double>& theoretical, double tol_ppm, double chain_mass);

    /// The neutral mass of the chain the fragments are cut from: its residues (fixed modifications included, I as L, as the
    /// fragment ladders count them), every shift in @p shifts, and one water.
    OPENMS_DLLAPI double chainMass(const std::string& sequence, const std::vector<double>& shifts, const FixedModifications& fixed_mods);

    /// identification.tsv's form of a log10 value: three significant digits in scientific notation ("2.71e-37",
    /// "4.10e-401", "1.00e+00"), written from the log10 so nothing underflows. "-1" when absent, "0" for -infinity.
    OPENMS_DLLAPI std::string formatLog10Scientific(const std::optional<double>& log10_value);

    /// The console's form: one significant digit ("3e-12", "4e-401"); "-" when absent, "0" for -infinity.
    OPENMS_DLLAPI std::string formatLog10Short(const std::optional<double>& log10_value);
  } // namespace ChanceMatch
} // namespace OpenMS
