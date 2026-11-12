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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

namespace OpenMS
{
  /// Where a Precursor stands in its cascade, as the operator view shows it (ADR-0065, CONTEXT.md "Precursor state").
  enum class PrecursorState
  {
    Acquiring,   ///< an MS2 handed to the instrument, none returned
    Searching,   ///< an MS2 returned, its search not finished
    NoId,        ///< searched, no proteoform (or its search was dropped from the backlog)
    Identified,  ///< a proteoform below the close-out threshold, or no threshold authored
    ClosedOut,   ///< an Identification score at or above the threshold
    Acquired     ///< returned, and nothing searches it (no identification configured, a Q scan, a trap variant)
  };

  /// One Precursor of the board: one row per committed nominal mass.
  struct OPENMS_DLLAPI BoardRow
  {
    int nominal = 0;
    double mass = 0.0;              ///< the commanded mass of its best identification; until then, of its first MS2
    int charge = 0;                 ///< the anchor charge of its best identification; until then, of its latest MS2
    PrecursorState state = PrecursorState::Acquiring;
    double score = -1.0;            ///< the Identification score the close-out reads; -1 = none
    int ms2 = 0;                    ///< MS2 scans sent for it
    int ms3 = 0;                    ///< MS3 scans sent for it
    double last_rt_s = 0.0;         ///< retention time (s) of its last activity
    std::string accession;
    std::string proforma;           ///< as identification.tsv writes it for the MS2 row
    bool below_threshold = false;   ///< identified, a threshold is authored and the score is below it
    double threshold = 0.0;         ///< that threshold, for "(below 70)"
    int precursor_id = 0;           ///< the precursor id of its best identification; 0 = none
    double coverage = -1.0;         ///< the Pooled coverage (fraction) of that identification's model; -1 = no model
  };

  /// Everything the view draws, copied under the board's lock.
  struct OPENMS_DLLAPI BoardSnapshot
  {
    double rt_s = 0.0;              ///< the latest retention time seen (s)
    double duration_min = 0.0;      ///< global.duration (min); 0 = scale to rt_s
    std::vector<BoardRow> live;     ///< the live set, most recently active first
    int ms1 = 0, ms2 = 0, ms3 = 0;  ///< scans sent, run totals
    int planned = 0;                ///< the last survey's plan units
    /// run totals, by row and by state: identified counts the Identified state only, so a closed-out row counts once
    int precursors = 0, identified = 0, closed_out = 0, no_id = 0;
    std::vector<double> id_rts_s;   ///< retention time (s) of every identification (a row's first result with a proteoform)
    std::vector<int> search_wait_ms;
    std::vector<int> search_ms;     ///< duration of every search, as scan_results.tsv's search_ms
    size_t search_waiting = 0;      ///< searches waiting on the pool now (SearchPool::waiting)
    size_t search_capacity = 0;     ///< how many may wait before the oldest is dropped (SearchPool::backlog); 0 = no pool
    std::string log_dir;
    std::string best_accession;     ///< the best-scoring row's protein, for the summary
    double best_score = -1.0;
    int proteins = 0;               ///< distinct accessions identified
  };

  /**
   * @brief The operator view's model (ADR-0065 push 2): one row per committed Precursor, keyed by nominal mass.
   *
   * Fed at the engine's own decision sites -- the drain's send, an MS2's return, the search's apply, the survey -- with
   * the values those sites log, and read by the view's render thread through snapshot(). One mutex; it takes no other
   * lock, so it may be called under analysis_mutex_ (processScan, the applies) and without it (the drain).
   */
  class OPENMS_DLLAPI PrecursorBoard
  {
  public:
    /// The latest retention time seen (s): processScan's entry.
    void setRt(double rt_s);

    /// The run folder, for the summary.
    void setLogDir(const std::string& log_dir);

    /// A command handed to the instrument (the drain). MS1 and AGC count; an MS2 creates or updates its Precursor's
    /// row (Acquiring unless it has a result) and records precursor_id -> nominal; an MS3 counts on that row.
    void onSent(const ScanCommand& cmd, int precursor_id);

    /// The units of the latest survey's plan.
    void onPlanned(int units);

    /// An MS2 of @p commanded_mass returned: Searching when it is submitted to the search, else Acquired (unless the
    /// row already has a result).
    void onReturned(double commanded_mass, bool searched);

    /// A searched MS2's result, from the apply: the match, the ProForma the identification row writes, whether it
    /// closed the species out, and the authored threshold. A result never lowers a row: ClosedOut is sticky, a
    /// proteoform beats none, the score is the best seen.
    /// Every value an identified row shows comes from its best identification (CONTEXT.md "Precursor board"): when this
    /// result becomes it, the row takes its @p commanded_mass, anchor @p charge and @p precursor_id with its score.
    void onResult(double commanded_mass, int charge, int precursor_id, const FragmentAnalysis::ProteoformMatch& match,
                  bool closed_out, const std::optional<double>& threshold);

    /// A pooled_identification.tsv row (ProteoformTracker's observer): records that precursor id's Pooled coverage, and
    /// returns the operator line when an MS3 fold (@p trigger an ion, not "MS2" or "MS3-DECIDE") raised it:
    /// "#<trigger_scan_id>  <nominal>  MS3 <ion>  cov 34.2% -> 41.0%  (score 180)", the score that id's best.
    std::optional<std::string> onPooled(int precursor_id, int nominal, const std::string& trigger,
                                        const std::string& trigger_scan_id, double coverage);

    /// A search dropped from the backlog: the row ends NoId unless it already has a proteoform.
    void onDropped(double commanded_mass);

    /// The wait (ms) of one search, from submit to start.
    void onSearchWait(int ms);

    /// The duration (ms) of one search: scan_results.tsv's search_ms.
    void onSearchTime(int ms);

    /// The live set (last activity within @p rt_window_s of the latest retention time), most recent first, and the
    /// run's counters.
    BoardSnapshot snapshot(double rt_window_s, double duration_min, size_t search_waiting, size_t search_capacity) const;

    /// The ProForma the board shows for @p match (the MS2 identification row's rendering).
    static std::string proformaOf(const FragmentAnalysis::ProteoformMatch& match);

  private:
    BoardRow& rowFor_(double mass);   // under mutex_

    mutable std::mutex mutex_;
    std::map<int, BoardRow> rows_;
    std::map<int, int> nominal_by_precursor_;
    std::map<int, double> coverage_by_precursor_;     ///< latest Pooled coverage per precursor id
    std::map<int, double> best_score_by_precursor_;   ///< best Identification score per precursor id
    std::vector<double> id_rts_s_;
    std::vector<int> search_wait_ms_;
    std::vector<int> search_ms_;
    double rt_s_ = 0.0;
    int ms1_ = 0, ms2_ = 0, ms3_ = 0, planned_ = 0;
    std::string log_dir_;
  };
} // namespace OpenMS
