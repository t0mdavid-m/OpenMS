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

#include <array>
#include <string>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Isobaric quantification component for FLASHIda (ADR-0038).
   *
   * Measures the reporter-ion channels of a returning QUANTIFICATION scan -- the `'Q'`-labelled
   * `ms_settings.ms2_quant`, rostered once per selected precursor -- and reports both the numbers
   * and the verdict. A differential verdict is what buys the identification scan.
   *
   * This is the only place reporter ions are measured, and it is deliberately NOT the scan that
   * gets bought. The reverse arrangement shipped until ADR-0038 and could not work: it screened a
   * scan whose activation cannot release the reporter, then acquired the one that can and never
   * read it.
   */
  class OPENMS_DLLAPI Quantification
  {
  public:
    /// Constructor — holds a reference to the shared Config object
    explicit Quantification(const Config& config);

    /// Outcome of measuring one quantification scan. Every field is logged (scan_results.tsv), so
    /// the value reported is by construction the value the engine decided on.
    struct Result
    {
      enum class Verdict
      {
        Differential,       ///< conditions differ beyond fold_change_threshold, OR one is wholly absent
        NotDifferential,    ///< measured cleanly, ratio inside the threshold
        IncompleteChannels, ///< an ASSIGNED channel was empty; no honest ratio exists
        ExtractionFailed    ///< the extractor produced no usable channel set for this spectrum
      };

      Verdict verdict = Verdict::ExtractionFailed;
      /// All N corrected channel intensities, in getChannelInformation() order. Every channel is
      /// reported, including any assigned to no condition -- gating ignores those, the log does not.
      std::vector<double> channels;
      /// Mean of each condition, in `quantification.conditions` order (== the ratio direction).
      std::array<double, 2> condition_means {{0.0, 0.0}};
      /// mean(conditions[0]) / mean(conditions[1]); -1.0 when a condition is wholly absent, so the
      /// ratio has no finite value. Callers distinguish that from "not measured" by the presence of
      /// condition_means, never by this field alone.
      double fold_change = -1.0;
    };

    /// Measure the reporter channels of one quantification scan.
    ///
    /// Reads every parameter but the spectrum from `config_.quantification()` -- the labelling
    /// scheme, the tolerance, the two conditions and the correction matrix. That reference was
    /// held and never read before ADR-0038.
    ///
    /// @param mzs       m/z values of the returning spectrum
    /// @param ints      intensity values, parallel to @p mzs
    /// @param length    number of entries in @p mzs / @p ints
    /// @param rt        retention time
    /// @param ms_level  MS level of the spectrum
    /// @param name      spectrum name (passed through to MSSpectrum)
    Result measure(const double* mzs, const double* ints, int length,
                   double rt, int ms_level, const char* name) const;

    /// Render a Verdict as the string logged in scan_results.tsv's `quant_verdict`.
    static const char* verdictName(Result::Verdict v);

    /// The accepted `quantification.labelling` values -- TopDownIsobaricQuantification's own valid
    /// strings MINUS "none", because `quantification.enabled` is the switch (ADR-0038).
    static const std::vector<std::string>& labellingNames();

    /// The channel NAMES of one labelling scheme, in getChannelInformation() order, or empty for an
    /// unknown scheme. Config calls this at load to resolve authored condition channel names to
    /// ordinals, so an unknown channel fails loudly at load rather than silently reading the wrong
    /// intensity at acquisition time.
    static std::vector<std::string> channelNames(const std::string& labelling);

  private:
    const Config& config_;
  };

} // namespace OpenMS
