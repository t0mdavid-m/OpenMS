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
#include <OpenMS/config.h>

namespace OpenMS
{

  /**
   * @brief Isobaric quantification component for FLASHIda.
   *
   * Owns the TMT reporter-ion differential abundance test used to decide
   * whether a precursor warrants a quantification follow-up MS2 scan.
   */
  class OPENMS_DLLAPI Quantification
  {
  public:
    /// Constructor — holds a reference to the shared Config object
    explicit Quantification(const Config& config);

    /**
     * @brief Determine whether a spectrum shows differential abundance across TMT channels.
     *
     * Extracts TMT6plex reporter-ion channels and compares the mean intensities of
     * the first three channels (sample 1) against the last three (sample 2).
     *
     * @param mzs        m/z values of the input spectrum
     * @param ints       intensity values of the input spectrum
     * @param length     number of entries in mzs / ints
     * @param rt         retention time in seconds
     * @param ms_level   MS level of the spectrum
     * @param name       spectrum name (passed to MSSpectrum)
     * @param reporter_mz_tol      reporter ion m/z tolerance
     * @param fold_change_threshold fold-change threshold above which a feature is
     *                              considered differentially abundant
     * @param only_one_condition   when true, treat a completely missing condition
     *                              as differentially abundant
     * @return true if the spectrum is differentially abundant, false otherwise
     */
    bool isDifferentiallyAbundant(const double* mzs,
                                  const double* ints,
                                  int length,
                                  double rt,
                                  int ms_level,
                                  const char* name,
                                  double reporter_mz_tol,
                                  double fold_change_threshold,
                                  bool only_one_condition);

  private:
    const Config& config_;
  };

} // namespace OpenMS
