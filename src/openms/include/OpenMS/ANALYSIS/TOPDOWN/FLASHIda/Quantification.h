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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/config.h>

#include <array>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Isobaric quantification component for FLASHIda (ADR-0038).
   *
   * Measures the reporter-ion channels of a returning QUANTIFICATION scan -- the `'Q'`-labelled
   * `ms_settings.ms2_quant`, rostered once per selected precursor -- and reports both the numbers
   * and the verdict. A verdict satisfying `quantification.identify` is what buys the identification
   * scan -- `differential` by default, which is ADR-0038's original rule (ADR-0039).
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

    // ---- the QUANTIFICATION GROUP, and its vote (ADR-0040) -------------------------------------
    //
    // Under precursor_charges "separate" one species gets one 'Q' per charge state. The reporter ion
    // is charge-INDEPENDENT, so a clean species must give the same ratio at every charge -- which is
    // what makes cross-charge disagreement evidence of chimericity, and what makes pooling REPORTER
    // measurements sound where pooling FRAGMENT evidence across charges is not (CONTEXT.md).
    //
    // The group is what buys the identification scan, never a member on its own.

    /// One quantification scan of a group: the command that produced it, and what it measured.
    struct GroupMember
    {
      int    tracking_id  = 0;
      int    precursor_id = 0;   ///< the WINNER's is what the bought 'R' inherits
      int    charge       = 0;
      double window_snr   = 0.0; ///< ranks the ID charge and breaks a tied vote (ADR-0040 d.4/d.5)
      double intensity    = 0.0; ///< this charge's own; weights the aggregate
      ScanCommand cmd {};        ///< handed to buildFollowUp as ctx, so the 'R' inherits its geometry
      bool   received     = false;
      Result result {};
    };

    /// What a completed group decided. Carries the winner outright so the caller needs no second
    /// lookup and the group can be closed immediately.
    struct Consensus
    {
      /// The synthesized measurement. Deliberately a plain Result, so quantBuysIdentification --
      /// and therefore `identify` and `enriched_in` -- work on it UNCHANGED (ADR-0039).
      Result result {};

      bool   has_winner = false;      ///< false = too few ballots; nothing is bought
      int    winner_charge = 0;
      int    winner_precursor_id = 0;
      ScanCommand winner_cmd {};

      int    agreeing = 0;            ///< ballots in the winning bloc
      int    total_ballots = 0;       ///< ballots cast; abstentions are NOT counted
      /// Every balloting charge, and whether it agreed. Index-parallel; drives consensus_charges.
      std::vector<int>  ballot_charges;
      std::vector<bool> ballot_agreed;
    };

    /// Register one 'Q' command as a member of @p group_id. Called at MS1, once per emitted 'Q'.
    void addMember(int group_id, const ScanCommand& cmd, int precursor_id);

    /// Record a returning scan's measurement. Returns its group_id, or 0 if the scan is not a
    /// member of any group (which is every scan in a non-quant run).
    int deposit(int tracking_id, const Result& result);

    /// Has every member of @p group_id returned? Mirrors Exploration's all_of(variants, received)
    /// -- and, like it, has NO timeout: a scan the instrument never returns leaves its group open.
    bool isComplete(int group_id) const;

    /// Run the vote. See Quantification.cpp for the policy; the short version is a DIRECTIONAL
    /// majority (verdict + which condition is enriched), abstentions excluded, ties broken by
    /// window_snr, and an intensity-weighted aggregate over the AGREEING members only.
    Consensus decide(int group_id) const;

    /// Drop a group and its tracking-id index entries. Call once the consensus has been used.
    void closeGroup(int group_id);

    /// Test seam: how many groups are currently open.
    size_t openGroupCount() const { return groups_.size(); }

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

    /// group_id -> its members, in emit order. std::map so iteration is deterministic; a run holds
    /// only the groups whose scans are still in flight, which is a handful.
    ///
    /// NO LOCK, and that is not an oversight: every access is from processScan, which the host
    /// pipeline serialises against itself (MaxDegreeOfParallelism = 1). The drain never touches
    /// this. That is exactly why it does NOT live beside FLASHIda::precursor_id_by_tracking_, whose
    /// leaf mutex exists solely because the drain reads that map -- putting it there would advertise
    /// a sharing that does not exist.
    std::map<int, std::vector<GroupMember>> groups_;

    /// tracking_id -> group_id, so a returning scan finds its group in one hop.
    std::unordered_map<int, int> group_by_tracking_;
  };

} // namespace OpenMS
