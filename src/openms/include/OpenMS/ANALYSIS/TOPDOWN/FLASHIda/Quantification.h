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
      /// For a GROUP result (ADR-0056) these are the aggregated channels the group was decided on.
      std::vector<double> channels;
      /// GEOMETRIC mean of each condition, 2^mean(log2), in `quantification.conditions` order (== the
      /// ratio direction) -- ADR-0056. 0 for any condition with an empty assigned channel, where a
      /// geometric mean is undefined.
      std::array<double, 2> condition_means {{0.0, 0.0}};
      /// condition_means[0] / condition_means[1] (the geometric ratio, ADR-0056); -1.0 when no
      /// finite ratio exists. Callers tell "measured, one condition wholly absent" from "not
      /// measured" by condition_means, never by this field alone.
      double fold_change = -1.0;
      /// Two-sided Welch p on the log2 channels (ADR-0056), computed under EVERY test whenever both
      /// conditions are complete with at least two channels each; -1.0 otherwise.
      double p_value = -1.0;
      /// Every assigned channel was measured. Only a complete result is a USABLE group member.
      bool complete = false;
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

    /// The channel -> verdict step (ADR-0056), shared by measure() on one scan and decide() on a
    /// group's aggregate, so a group is decided by exactly the function a scan is. A channel below
    /// @p empty_floor (or not above 0) is empty: kEmptyChannel for raw intensities, 0 for an
    /// aggregate of shares, which can legitimately sit below 1e-3.
    static Result evaluate(const std::vector<double>& channels, const QuantConfig& cfg, double empty_floor);

    /// Render a Verdict as the string logged in scan_results.tsv's `quant_verdict`.
    static const char* verdictName(Result::Verdict v);

    // ---- the QUANTIFICATION GROUP (ADR-0040, decided by ADR-0056) -------------------------------
    //
    // Under precursor_charges "separate" one species gets one 'Q' per charge state. The reporter ion
    // is charge-INDEPENDENT, so the members measure the same population through different windows,
    // and pooling REPORTER measurements is sound where pooling FRAGMENT evidence is not (CONTEXT.md).
    //
    // The group is what buys the identification scan, never a member on its own. Since ADR-0056 it
    // is decided by ONE test on the aggregate of its usable members -- there is no vote.

    /// One quantification scan of a group: the command that produced it, and what it measured.
    struct GroupMember
    {
      int    tracking_id  = 0;
      int    precursor_id = 0;   ///< the WINNER's is what the bought 'R' inherits
      int    charge       = 0;
      double window_snr   = 0.0; ///< ranks the identification member (the purest window wins)
      double intensity    = 0.0; ///< tiebreak only, after window_snr
      ScanCommand cmd {};        ///< handed to buildFollowUp as ctx, so the 'R' inherits its geometry
      /// Survey re-planning: every charge this Q isolates (ScanCommandQueue::isolatedCharges: its anchor, then
      /// its notches), and the survey that decided it. The bought R may move only among the charges its group's
      /// Q scans covered (owner 2026-10-23).
      std::vector<int> isolated_charges;
      int    survey_id    = 0;
      bool   received     = false;
      Result result {};
    };

    /// What a completed group decided. Carries the winner outright so the caller needs no second
    /// lookup and the group can be closed immediately.
    struct GroupDecision
    {
      /// The group's result. Deliberately a plain Result, so quantBuysIdentification -- and
      /// therefore `identify` and `enriched_in` -- work on it UNCHANGED (ADR-0039).
      Result result {};

      /// True whenever any member returned (ADR-0056 decision 9): even an undecided group names the
      /// window an identification scan would be built from, so `identify: all` can buy it.
      bool   has_winner = false;
      int    winner_charge = 0;
      int    winner_precursor_id = 0;
      ScanCommand winner_cmd {};
      /// The union of isolated_charges over every member whose Q returned, sorted and unique: the charges the
      /// bought R may move to.
      std::vector<int> q_covered_charges;
      int    winner_survey_id = 0;           ///< the winner's GroupMember::survey_id

      int    usable = 0;              ///< members whose every assigned channel was measured
      int    members = 0;             ///< members registered for the group
      std::vector<int> usable_charges; ///< drives group_charges
    };

    /// Register one 'Q' command as a member of @p group_id. Called at MS1, once per emitted 'Q'. The charges it
    /// isolates are read off @p cmd; @p survey_id is the survey that decided it (GroupMember).
    void addMember(int group_id, const ScanCommand& cmd, int precursor_id, int survey_id = 0);

    /// Record a returning scan's measurement. Returns its group_id, or 0 if the scan is not a
    /// member of any group (which is every scan in a non-quant run).
    int deposit(int tracking_id, const Result& result);

    /// Has every member of @p group_id returned? Mirrors Exploration's all_of(variants, received)
    /// -- and, like it, has NO timeout: a scan the instrument never returns leaves its group open.
    bool isComplete(int group_id) const;

    /// Decide the group (ADR-0056). Usable members are scaled to reporter shares, combined per
    /// channel, and evaluate()d once; with no usable member, wholly-absent-against-full members
    /// decide through quantification.empty_condition; both paths need min_usable_scans members.
    GroupDecision decide(int group_id) const;

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
