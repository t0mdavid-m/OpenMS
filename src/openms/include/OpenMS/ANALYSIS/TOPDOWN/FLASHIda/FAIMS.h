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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>

#include <vector>

namespace OpenMS
{

  /**
   * @brief FAIMS CV cycling state machine for FLASHIda.
   *
   * Manages the per-CV adaptive skip state (skip amounts, skip counts,
   * current CV index) and the round-robin CV cycling logic.
   * Configuration values (cv_values, max_cv_skip, precursor_threshold)
   * are read from FAIMSConfig at construction time.
   */
  class OPENMS_DLLAPI FAIMS
  {
  public:
    /// Construct from Config; initialises all per-CV state vectors.
    explicit FAIMS(const Config& config);

    /**
     * @brief Returns true if FAIMS is in use at all, i.e. any CV is configured.
     *
     * Distinct from isCycling(): a single configured CV is a perfectly ordinary FAIMS run at a
     * fixed compensation voltage. This used to be `cv_values.size() > 1`, which conflated "FAIMS
     * is off" with "there is nothing to cycle between", so a single-CV method silently acquired
     * with FAIMS at whatever the instrument method carried. See ADR-0012.
     */
    bool isEnabled() const;

    /**
     * @brief Returns true if there is more than one CV to rotate between.
     *
     * Guards the CV-transition MS1 push and the adaptive skip policy. Both are meaningless with a
     * single CV, and pushing a CV-transition MS1 after every MS1 when the "next" CV is always the
     * current one would double the survey rate.
     */
    bool isCycling() const;

    /// Returns the CV value at the current cycling index.
    double currentCV() const;

    /**
     * @brief Update per-CV adaptive skip policy after MS1 deconvolution.
     * @param cv       CV value of the scan that was just processed.
     * @param precursor_count  Number of MS2 commands pushed for this scan.
     */
    void updateSkip(double cv, int precursor_count);

    /**
     * @brief Advance to the next non-skipped CV and return its value.
     *
     * Increments the cycling index first, wrapping at the end of the list.
     * CVs whose skip counter has not yet been exhausted are bypassed.
     * Falls back to the current CV if all entries are being skipped.
     */
    double advanceToNextCV();

    /// Returns the current skip amount (spacing) for CV at \p index.
    int cvSkipAmount(size_t index) const;

  private:
    std::vector<double> cv_values_;
    std::vector<int>    cv_skip_amount_;   ///< per-CV: current skip spacing
    std::vector<int>    cv_skip_count_;    ///< per-CV: current skip counter
    int current_cv_index_ = 0;
    int max_cv_skip_ = 0;
    int precursor_threshold_ = 15;
    bool enabled_ = false;
  };

} // namespace OpenMS
