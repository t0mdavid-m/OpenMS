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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ConsoleView.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Exploration.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FAIMS.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorBoard.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteomeSearch.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Quantification.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/SearchPool.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <map>
#include <set>
#include <memory>
#include <mutex>

namespace OpenMS
{

  /**
   * @brief FLASHIda class for real time deconvolution
   * This class contains functions to perform deconvolution (by SpectralDeconvolution) for the spectrum received from Thermo iAPI.
   * Also precursor selection is done in this class.
   * The functions in this class are invoked in C# Thermo iAPI side through the functions in FLASHIdaBridgeFunctions class
   * @see FLASHIdaBridgeFunction, https://stackoverflow.com/questions/31417688/passing-a-vector-array-from-unmanaged-c-to-c-sharp
   */
  class OPENMS_DLLAPI FLASHIda
  {
  public:
    typedef FLASHHelperClasses::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHHelperClasses::LogMzPeak LogMzPeak;

    using PTMSite = FragmentAnalysis::PTMSite;

    /// constructor that takes string input argument
    explicit FLASHIda(char *arg);

    /// destructor
    ~FLASHIda();

    /// copy constructor (deleted: ofstreams are non-copyable)
    FLASHIda(const FLASHIda&) = delete;

    /// move constructor
    FLASHIda(FLASHIda&& other) = default;

    /// assignment operator (deleted: ofstreams are non-copyable)
    FLASHIda& operator=(const FLASHIda&) = delete;

    /// Process an incoming scan and enqueue resulting commands.
    /// @param rt retention time in SECONDS -- the engine's unit (rt_window, inclusion rows, target
    ///        logs). The C export converts the host's minutes; a direct caller passes seconds. The two
    ///        log streams write it back in minutes (ADR-0051).
    /// @param instrument_scan_number the number the INSTRUMENT assigned this scan (Header["Scan"] on
    ///        the C# side). <= 0 means "not supplied" and the ida.log falls back to the tracking id
    ///        (ADR-0035). faims_cv deliberately LOSES its default here: a defaulted trailing
    ///        parameter would let every existing call site keep compiling while silently passing 0,
    ///        which collapses every ida.log entry onto one map key with nothing to notice it.
    int processScan(const double* mzs, const double* ints, int length,
                    double rt, int ms_level, const char* scan_description,
                    double faims_cv, int instrument_scan_number);

    /// Dequeue the next scan command by priority. Returns 1 if command filled, 0 if error.
    int getNextScanCommand(ScanCommand& out);

    /// Get the next monotonically increasing tracking ID (thread-safe)
    int getNextTrackingId();

    /// Block until the search pool is idle -- nothing waiting, nothing running -- and return the commands
    /// the search workers' applies pushed since the previous call (ADR-0062 d.7). For HARNESSES only: a
    /// drive that waits after every feed sees exactly the order processScan alone produced before the pool
    /// existed. Production never calls it. Must not be called while analysis_mutex_ is held (a worker's
    /// apply would wait for the caller). 0 at once when the engine has no pool (no identification section).
    int waitForSearchIdle();

    /// Test-only access (push / queue / queueSize / explorationActive) lives in the test header
    /// FLASHIda_TestAccess.h via this friend, so test scaffolding stays out of the production API.
    friend struct FLASHIdaTestAccess;

  private:
    /// Configuration object (owns all parsed config values)
    Config config_;

    /// Owns all output logging: the 4 streams, TSV headers, row writers, and the ida_log parser
    IdaLogger logger_;

    /// The operator view's model (ADR-0065 push 2): one row per committed Precursor. Its own mutex; fed under
    /// analysis_mutex_ (processScan, the applies) and without it (the drain).
    PrecursorBoard board_;

    /// Scan command queue (owns queue, tracking IDs, pending map, command builders)
    ScanCommandQueue queue_;

    /// Deconvolution engine (owns SpectralDeconvolution, MS1 result, MS2 result)
    Deconvolution deconv_;

    /// The MONITOR scan's OWN deconvolution engine (ADR-0042). Null unless a level enables
    /// monitor_ms1, so the feature costs nothing -- not even a second averagine table -- when off.
    ///
    /// Private, and that is the point. deconv_ holds ONE SpectralDeconvolution shared between MS1
    /// and MSn, and deconvolveMS1 restores the global charge window on it (see the LOAD-BEARING
    /// comment in Deconvolution.cpp). Routing a monitor scan through deconv_ mid-sweep would reset
    /// that window under the sweep. It self-corrects, because deconvolveMSn re-narrows on every
    /// call -- but "self-corrects" is an inference across two call sites, and a private instance
    /// makes it impossible instead. Built with the same toleranceList() as deconv_, so a monitor
    /// reading is numerically comparable with a survey's.
    ///
    /// Precedent: Exploration::exploration_deconv_ does the same thing, incidentally. Here it is
    /// deliberate. Do NOT "simplify" this away in favour of deconv_.
    std::unique_ptr<Deconvolution> monitor_deconv_;

    /// Fragment analysis (tag-based matching, MS2 mass queries, PTM ambiguity)
    FragmentAnalysis fragments_;

    /// Precursor selection, targeting, mass exclusion (owns all selection state)
    PrecursorSelection selection_;

    /// Isobaric quantification (TMT reporter-ion differential abundance test)
    Quantification quant_;

    /// Mutex protecting analysis state: deconv_, selection_, exploration_, faims_, quant_, fragments_,
    /// tracker_, ms2_context_cache_ and next_precursor_id_.
    ///
    /// TWO acquisition sites (ADR-0062): processScan (function-scoped, so it is held across the whole
    /// deconvolution), and a search worker's APPLY step (applyRegularSearch_ / applyVariantSearch_ run
    /// under it; the search itself holds no lock). getNextScanCommand does NOT acquire it and must not
    /// start -- the instrument event thread calls the drain, and making it wait on a deconvolution to
    /// write one TSV row and read one map entry is what this arrangement exists to prevent.
    ///
    /// It no longer serialises the logger. IdaLogger owns a mutex per stream, which matters rather
    /// than being tidier: one shared logger lock would put the drain's writeScanCommandRow behind
    /// processScan's writeScanResultRow / writeIdentificationRow / writeIDALogEntry -- each ending in
    /// a synchronous flush -- and reintroduce exactly the coupling removed here, through a different
    /// mutex, in a way the drain-blocking test cannot see because that test is scoped by mutex name.
    ///
    /// Pinned by FLASHIda_ProcessScan_test's drain_completes_while_analysis_mutex_held (the drain
    /// must not wait on it) and process_scan_still_blocks_while_analysis_mutex_held (processScan must
    /// still take it -- deleting THIS lock is the wrong way to satisfy the first test).
    mutable std::mutex analysis_mutex_;

    /// Atomic flag: true when any exploration group is active (set by processScan, read by getNextScanCommand)
    std::atomic<bool> exploration_active_{false};

    /// Per-level active-group counts, [0] = level 2, [1] = level 3 (ADR-0042). Same producer, same
    /// consumer and the same release/acquire discipline as exploration_active_ above -- written by
    /// processScan at the four sites that write that flag, read lock-free by getNextScanCommand.
    ///
    /// A SEPARATE atomic per level, not one packed value: the drain reads them independently and
    /// they are only ever compared to zero, so no cross-level consistency is needed and no reader
    /// can observe a torn pair. They deliberately do NOT replace exploration_active_ -- that flag
    /// still gates the cycle-time MS1, which is a whole-run question, not a per-level one.
    std::atomic<int> exploration_active_lvl_[2]{};

    /// Atomic FAIMS CV: current CV value (set by processScan after advanceToNextCV, read by getNextScanCommand)
    std::atomic<double> current_faims_cv_{0.0};

    FAIMS faims_;   ///< FAIMS CV cycling state machine

    /// Exploration CE sweep engine (owns groups, variants, scoring)
    Exploration exploration_;

    /// Per-precursor proteoform tracking model (accumulates fragment observations across MS2/MS3)
    ProteoformTracker tracker_;

    /// MS2 context cache keyed by MS3 tracking ID (for non-exploration identification)
    std::unordered_map<int, Exploration::MS2Context> ms2_context_cache_;

    /// P5: monotonic, deterministic-per-run, base-10 precursor identity counter. A NEW id is allocated
    /// for each precursor selected at MS1 (allocPrecursorId_); carried by all of that precursor's
    /// MS2 / exploration-variant / MS3 scans via precursor_id_by_tracking_. Logged as a plain decimal
    /// (NOT base-94 like tracking_id). Starts at 1 (0 = "no precursor / untracked").
    int next_precursor_id_ = 1;

    /// P5: tracking_id -> precursor_id, the propagation map. Stamped when an MS2 command is issued at MS1
    /// (allocPrecursorId_ per precursor) and when child commands inherit their parent's precursor_id.
    /// Read on scan return to key the ProteoformTracker model and source the precursor_id log column.
    ///
    /// Guarded by its OWN leaf mutex, NOT by analysis_mutex_. This is the only mutable state the
    /// instrument-event-thread drain and the deconvolution thread genuinely share, and giving it a
    /// dedicated lock is what stops the drain waiting out a whole deconvolution to log one row.
    ///
    /// The lock is for the CONTAINER, not for the value. Value visibility is already guaranteed
    /// without it: every write is immediately followed by queue_.push(), which RELEASES queue_mutex_,
    /// and the drain's dequeue() ACQUIRES the same mutex -- a happens-before edge, so a dequeued
    /// command's entry is final by the time the drain reads it. What genuinely needs protecting is
    /// the concurrent INSERT: an insert that rehashes while find() is walking a bucket is undefined
    /// behaviour whichever value would have won.
    ///
    /// Deliberately unreachable except through the accessors below -- a raw subscript no longer
    /// compiles, which is a stronger guarantee than any test we can write for that rehash race on
    /// this toolchain (MSVC ships no thread sanitizer).
    std::unordered_map<int, int> precursor_id_by_tracking_;

    /// ADR-0061 d.7 / ADR-0062 d.6: one closed-out species -- the RT (seconds) of its latest close-out and the
    /// precursor ids whose MS2 closed it out. Analysis state (read and written under analysis_mutex_).
    struct ClosedOut
    {
      double rt = 0.0;
      std::set<int> closing_precursor_ids;
    };
    /// nominal mass -> its close-out. A later scan of the species within rt_window is not searched (scan_results
    /// search = closed_out_skip) -- except the closing Precursors' own scans (owner 2026-11-08: their sweep variants,
    /// follow-up and bought R are still searched, so characterization still gets its winner, ADR-0061 d.9).
    /// Expired entries are erased as they are read.
    std::unordered_map<int, ClosedOut> closed_out_;

    /// ADR-0068: one rescue scan per DETECTION -- one species of one survey's selection. Under `separate` a detection's
    /// charges are separate precursor ids, so the cap keys on (selecting survey's tracking id, nominal mass), recorded
    /// per precursor id when a survey plans it (only with ms_settings.ms2_rescue set). A rescue carries its trigger's
    /// precursor id, so its own return finds the detection taken and buys nothing. Analysis state (under
    /// analysis_mutex_); never erased.
    std::unordered_map<int, std::pair<int, int>> detection_of_precursor_;
    std::set<std::pair<int, int>> rescued_detections_;

    /// Leaf lock for precursor_id_by_tracking_.
    ///
    /// Lock hierarchy: analysis_mutex_ -> { queue_mutex_, precursor_map_mutex_, the IdaLogger stream
    /// mutexes, the SearchPool's mutex }, all leaves. This one is taken ONLY inside the two accessors below, and neither of
    /// them calls anything -- that is what makes it a leaf and the hierarchy acyclic. Keep it that
    /// way: a queue callback that consulted this map would invert the edge.
    ///
    /// std::mutex is NOT recursive, so the accessors must never call each other. The obvious-looking
    /// "read the existing id first" version of the setter self-deadlocks on the first MS2 of the run.
    mutable std::mutex precursor_map_mutex_;

    /// P5: allocate the next precursor_id (monotonic). Call once per MS1-selected precursor.
    /// processScan-only, and processScan is serialised against itself by the host pipeline, so this
    /// needs no lock -- it is not part of the shared surface.
    int allocPrecursorId_() { return next_precursor_id_++; }

    /// ADR-0040: monotonic quantification-GROUP id. One per species per survey, shared by every 'Q'
    /// that species emits -- which under precursor_charges "separate" is one per charge state.
    ///
    /// Deliberately NOT a second precursor_id. precursor_id means "one Precursor, one proteoform
    /// model"; this means "these scans measure the same reporter population". Collapsing them would
    /// pool FRAGMENT evidence across charges along with the reporters, which is the defect
    /// CONTEXT.md's Precursor entry warns about. The group->members map itself lives on
    /// Quantification, not here: the drain never reads it, so it needs none of the locking that
    /// precursor_id_by_tracking_ above does.
    int next_quant_group_id_ = 1;
    int allocQuantGroupId_() { return next_quant_group_id_++; }

    /// P5: look up the precursor_id for a tracking_id, or 0 if untracked.
    /// Self-locking so that no call site can forget it -- including the drain's, which is the whole
    /// reason the lock exists.
    int precursorIdForTracking_(int tracking_id) const
    {
      std::lock_guard<std::mutex> lk(precursor_map_mutex_);
      auto it = precursor_id_by_tracking_.find(tracking_id);
      return (it != precursor_id_by_tracking_.end()) ? it->second : 0;
    }

    /// P5: stamp a tracking_id with its precursor_id. Self-locking; must NOT call the reader above.
    void setPrecursorForTracking_(int tracking_id, int precursor_id)
    {
      std::lock_guard<std::mutex> lk(precursor_map_mutex_);
      precursor_id_by_tracking_[tracking_id] = precursor_id;
    }

    /// Survey re-planning: forget the precursor ids of commands a re-plan took out of the queue -- dropped,
    /// or reconciled away before they were ever inserted. Self-locking, and a leaf like the two above.
    void forgetTracking_(const std::vector<int>& tracking_ids)
    {
      std::lock_guard<std::mutex> lk(precursor_map_mutex_);
      for (int id : tracking_ids) precursor_id_by_tracking_.erase(id);
    }

    /// Stamp a command with its precursor_id and enqueue it -- in that order, always.
    ///
    /// The ORDER is the entire point of this helper. The drain reads the map for a command it has
    /// just dequeued, and dequeue() acquires the same queue_mutex_ that push() releases, so stamping
    /// BEFORE the push puts the write on the happens-before side of that edge and the drain is
    /// guaranteed a final value. Stamping after would be a genuine defect -- one that no compiler and
    /// no CPU can introduce, and only an editor can. Six call sites each remembering the rule is six
    /// chances to get it wrong; one helper is none.
    void stampAndPush_(ScanCommand& cmd, int precursor_id)
    {
      setPrecursorForTracking_(cmd.scan_id, precursor_id);
      queue_.push(cmd);
    }

    /// Publish exploration state for the drain's lock-free reads: the whole-run flag AND both
    /// per-level counts, together (ADR-0042).
    ///
    /// One helper for the same reason stampAndPush_ is one: there are four sites that must publish
    /// this, and asking each of them to remember three atomics is twelve chances to get it wrong.
    /// A per-level count that lags its flag would strand a monitor scan's cadence -- armed against a
    /// sweep the drain no longer believes is running, or the reverse.
    void publishExplorationState_()
    {
      exploration_active_.store(exploration_.activeGroupCount() > 0, std::memory_order_release);
      exploration_active_lvl_[0].store(exploration_.activeGroupCount(2), std::memory_order_release);
      exploration_active_lvl_[1].store(exploration_.activeGroupCount(3), std::memory_order_release);
    }

    /**
     * @brief The MS1 SURVEY arm of processScan: deconvolve -> select -> push MS2 -> FAIMS cycle.
     *
     * Extracted verbatim from the `ms_level == 1` branch (ADR-0042) so that branch can fork into a
     * deciding arm and an observing one with NOTHING AFTER THE FORK. Every acquisition decision an
     * MS1 makes lives in here; a side effect added to this function is therefore unreachable from a
     * scan that took the other arm, by construction rather than by a guard.
     *
     * @return commands pushed; becomes processScan's return_code.
     */
    int runSurveyMS1_(const double* mzs, const double* ints, int length, double rt,
                      double faims_cv, int parent_tracking_id, const ScanCommand& parent_ctx,
                      int instrument_scan_number,
                      IdaLogger::ScanRowDescriptor& results_row,
                      std::vector<std::string>& child_ids);

    /**
     * @brief The MS1 OBSERVING arm: deconvolve, record the masses, decide nothing (ADR-0042).
     *
     * `static` IS LOAD-BEARING. It is the entire neutrality guarantee, not a style choice.
     *
     * A static member function has no `this`, so inside it `selection_`, `queue_`, `faims_`,
     * `tracker_`, `quant_` and `exploration_` are not merely discouraged -- naming any of them is
     * `error: invalid use of non-static data member`. A monitor scan therefore cannot age an
     * exclusion map, mint a precursor id, push a command or advance the FAIMS wheel, because there
     * is no object in scope through which to do it. Every parameter is passed explicitly, so this
     * signature IS the function's complete list of capabilities, and widening it is a header diff a
     * reviewer sees.
     *
     * DO NOT remove `static` to "simplify". DO NOT add a FLASHIda& parameter. Either change silently
     * converts a compile-time guarantee back into a promise, and the failure it stops -- a monitor
     * scan quietly corrupting acquisition memory -- is invisible in every log.
     *
     * @param d MUST be monitor_deconv_, never deconv_ (see that member's comment).
     */
    static void observeMonitoringMS1_(Deconvolution& d, const Config& cfg,
                                      const double* mzs, const double* ints, int length,
                                      double rt, double faims_cv,
                                      IdaLogger::ScanRowDescriptor& results_row);

    /**
     * @brief Close a identification.method refine discriminating round (ADR-0049): decide, then re-dispatch MS3 from
     *        the decided model through Exploration::dispatchFromModel, replaying the WINNER MS2's own command as stage 0.
     *
     * Pushes the new MS3 (each stamped with @p precursor_id and appended to @p child_ids), caches their MS2 contexts so
     * they are identified on return, and writes commands_pushed / child_ids into @p results_row. Without a winner
     * spectrum or command nothing is re-dispatched and a warning line is printed.
     *
     * @return commands pushed; becomes processScan's return_code.
     */
    int refineRoundClose_(int precursor_id, const ScanCommand& parent_ctx, std::vector<std::string>& child_ids,
                          IdaLogger::ScanRowDescriptor& results_row);

    // ===== Per-proteoform MS3 re-planning (ADR-0066, characterization.per_proteoform) ===========================

    /// Swap in a governed MS3 plan (Exploration::GovernedMs3) for @p precursor_id, the one site that calls
    /// queue_.replaceMs3Plan: stamp every new id's precursor id first (the drain may dequeue it at once), run the
    /// test seam, swap in one queue_mutex_ hold, then -- outside it -- forget the withdrawn and reconciled ids, drop the
    /// withdrawn ids' MS2 contexts, cache the inserted ids' and undo the withdrawn and reconciled items. Appends the
    /// inserted ids to @p child_ids and prints `[MS3-REPLAN]`. Caller holds analysis_mutex_. Returns the commands
    /// inserted.
    int pushGovernedMs3_(Exploration::GovernedMs3& gov, int precursor_id, std::vector<std::string>& child_ids);

    /// Give back what unsent governed MS3 items took at plan time (ADR-0066): per owner, tracker_.withdrawMs3 (a
    /// discriminating round it completes is decided, not re-dispatched); every sweep group is dropped (an
    /// `[MS3-REPLAN-INVARIANT]` line if one cannot be), then the exploration state is republished.
    /// @p committed_owner_nominals: (precursor id, nominal) of items that went out or stay queued -- their siblings
    /// keep the stamp. Caller holds analysis_mutex_.
    void undoMs3Items_(const std::vector<Ms3PlanItem>& items, const std::set<std::pair<int, int>>& committed_owner_nominals);

    /// Take the queue's expired governed items (cleanupExpired's inbox) and undo them with their entries' committed
    /// nominals (a committed sibling keeps the stamp): their ids are forgotten, their MS2 contexts dropped. A no-op with
    /// characterization.per_proteoform off. Caller holds analysis_mutex_.
    void undoExpiredMs3_();

    /// Test seam between a governed MS3 plan and its swap (replaceMs3Plan), so a test can commit an MS3 in exactly
    /// that gap. Empty in production; only FLASHIda_TestAccess.h sets it.
    std::function<void()> before_ms3_replan_;
    void beforeMs3Replan_()
    {
      if (before_ms3_replan_) before_ms3_replan_();
    }

    /// Steady-clock reference for timestamps
    std::chrono::steady_clock::time_point engine_start_time_;

    /// Survey re-planning: the last survey's survey-level stamps -- writes of candidates that cleared
    /// exclusion and were then refused, so no plan unit owns them. Folded at the start of the next survey,
    /// beside the stamps the drain committed since. Nothing reads acquisition memory in between, so holding
    /// them until then is identical to applying them at once.
    std::vector<MemoryStamp> survey_level_stamps_;

    /// Test seam between a survey's plan and its swap (replacePlan), so a test can commit work in exactly
    /// that gap. Empty in production; only FLASHIda_TestAccess.h sets it.
    std::function<void()> before_replan_;
    void beforeReplan_()
    {
      if (before_replan_) before_replan_();
    }

    // ===== The search pool (ADR-0062) ===========================================================
    // processScan copies what a search needs under analysis_mutex_ and returns; a worker searches with
    // no lock held, then applies under analysis_mutex_ (the second acquisition site). The two job kinds
    // are exactly the matcher's two callers on an MS2: the regular MS2 (Exploration::searchNextLevel)
    // and the MS2 exploration variant (Exploration::searchVariant). Nothing else is searched here.

    /// When a job was submitted and when its search() ran (ADR-0061 push 4b): set by processScan and by the
    /// worker, read by the apply on that same worker, which writes them as scan_results' search_wait_ms /
    /// search_ms. Steady clock, so wall-clock steps cannot make either negative.
    struct SearchTiming
    {
      std::chrono::steady_clock::time_point submitted{};
      std::chrono::steady_clock::time_point started{};
      std::chrono::steady_clock::time_point finished{};
    };

    /// A regular MS2's search: the copies its apply needs once deconv_ has moved on to later scans.
    struct RegularSearch
    {
      int precursor_id = 0;
      std::string parent_id_str;
      ScanCommand parent_ctx{};
      double rt = 0.0;
      int instrument_scan_number = 0;
      bool has_spectrum = false;                ///< deconv_.hasStoredMS2() when the job was made
      DeconvolvedSpectrum spectrum;             ///< deconv_.storedMS2(), copied
      IdaLogger::ScanRowDescriptor results_row{};   ///< the half-filled row; deconv_spectrum re-pointed at `spectrum`
      std::vector<std::string> child_ids;
      int commands_pushed = 0;
      Exploration::NextLevelSearch result;      ///< filled by the worker
      SearchTiming timing;
    };

    /// An MS2 exploration variant's search.
    struct VariantSearch
    {
      int precursor_id = 0;
      std::string parent_id_str;
      int instrument_scan_number = 0;           ///< the scan the proteome search names its spectrum by (ADR-0061 push 4c)
      Exploration::VariantInput input;          ///< prepareVariant's copy (spectrum, share, slot); input.prsm by the worker
      IdaLogger::ScanRowDescriptor results_row{};
      FragmentAnalysis::ProteoformMatch match;  ///< filled by the worker
      SearchTiming timing;
    };

    /// The regular-MS2 tail from MS3 targeting to the row fields (processScan inline, or a worker's apply).
    /// @p searched: the worker's search; nullptr runs initiateNextLevel inline, as before the pool.
    /// @p dropped: the backlog dropped the search -- nothing is matched, planned or identified, and the
    /// identification columns keep their not-run sentinels (ADR-0062 d.5). Returns commands pushed.
    /// @p match_spec: the spectrum matched against; @p row_spec: the one the row reports (nullptr = none stored).
    int finishRegularMS2_(const Exploration::NextLevelSearch* searched, bool dropped, const DeconvolvedSpectrum& match_spec,
                          const DeconvolvedSpectrum* row_spec, const ScanCommand& parent_ctx, int precursor_id,
                          const std::string& parent_id_str, double rt, int instrument_scan_number, bool is_quant_scan,
                          int commands_pushed, std::vector<std::string>& child_ids,
                          IdaLogger::ScanRowDescriptor& results_row, std::vector<IdaLogger::IdRowDescriptor>& id_rows);

    /// The exploration-variant tail: push the group's commands, seed the MS3 contexts, fill the row and the
    /// identification rows. @p mass_count / @p spectrum: the variant's own deconvolution. Returns commands pushed.
    int finishVariantMS2_(Exploration::FeedResultInfo& expl_result, int mass_count, const DeconvolvedSpectrum* spectrum,
                          int precursor_id, const std::string& parent_id_str, bool dropped,
                          IdaLogger::ScanRowDescriptor& results_row, std::vector<IdaLogger::IdRowDescriptor>& id_rows);

    /// ADR-0061 d.7-d.10 / ADR-0062 d.3: when @p match's FLASHTnT Identification score reaches
    /// identification.threshold, close its species out -- arm the mass and m/z bars keyed on the commanded precursor
    /// (@p commanded_mass, @p isolation_mz) stamped at @p rt (seconds), plus an inferred mass that is the same species,
    /// withdraw the species' wholly unsent plan units at every CV, and record it for the skip. A species whose
    /// inclusion row names charges is never closed out. Caller holds analysis_mutex_. True when it closed out.
    bool closeOutIfIdentified_(const FragmentAnalysis::ProteoformMatch& match, double commanded_mass,
                               double isolation_mz, double rt, int precursor_id);

    /// Whether @p match reaches identification.threshold for a commanded precursor of @p commanded_mass: the score half
    /// of closeOutIfIdentified_, with no side effect (ADR-0064 buys its quantification MS3 on the same call).
    bool reachesThreshold_(const FragmentAnalysis::ProteoformMatch& match, double commanded_mass) const;

    /// ADR-0064: whether the species of @p mass is closed out within rt_window of @p rt by ANY Precursor, the closing
    /// one included -- so a species is quantified at MS3 once per close-out, not once per closing scan.
    bool speciesClosedOut_(double mass, double rt) const;

    /// ADR-0068: record that the survey @p survey_tracking_id planned @p precursor_id for the species of @p mono_mass
    /// (a no-op without ms_settings.ms2_rescue). Called at both allocPrecursorId_ sites.
    void noteDetection_(int precursor_id, int survey_tracking_id, double mono_mass);

    /// ADR-0068: claim the one rescue of @p precursor_id's detection; false if that detection already bought one.
    bool rescueDetectionFree_(int precursor_id, double mono_mass);

    /// Whether a scan of @p mass belonging to @p precursor_id is skipped (ADR-0062 d.6): its species was closed out
    /// within rt_window of @p rt by ANOTHER Precursor. Analysis state.
    bool isClosedOut_(double mass, double rt, int precursor_id);

    /// A worker's apply step: take analysis_mutex_, finish the scan, write its rows. Returns commands pushed.
    int applyRegularSearch_(RegularSearch& job);
    int applyVariantSearch_(VariantSearch& job);

    /// Hand a job to the pool; if the backlog drops the oldest waiting job, run its drop() here, under the
    /// caller's analysis_mutex_. Caller holds analysis_mutex_.
    void submitSearch_(std::unique_ptr<SearchPool::Job> job);

    /// Test seam run on the worker before every search (FLASHIda_TestAccess.h: block a worker, observe that
    /// the search runs while a test holds analysis_mutex_). Set before the first feed; empty in production.
    std::function<void()> before_search_;

    /// The Protein reference, loaded and prepared once per ion-type set (ADR-0061 push 4c): identification.reference.fasta
    /// read as authored, or the inline sequence as one entry with the accession "sequence". FLASHTnT builds its mass
    /// ladders from ion_type when it prepares, so each set of ion types an MS2 activation of this method has
    /// (FragmentAnalysis::getIonTypesForFragmentationMethod; b/y always) gets its own prepared search, keyed by that set.
    /// Empty without an identification section. Immutable after the constructor, so every search worker reads it
    /// unlocked; each worker thread owns one ProteomeSearch::Worker per entry (FLASHIda.cpp, thread_local). Declared
    /// BEFORE search_pool_, so the pool -- and with it every worker thread and its Workers -- is gone before this is.
    std::map<std::vector<std::string>, std::unique_ptr<ProteomeSearch>> proteomes_;

    /// The pool. Null without an identification section (nothing is ever searched then). Stopped and joined
    /// in ~FLASHIda's body, before any member a running apply touches is destroyed.
    std::unique_ptr<SearchPool> search_pool_;

    /// The operator view (ADR-0065 push 2): null unless stdout is a VT-capable console. Stopped in ~FLASHIda's
    /// body before the pool, since its render thread reads the pool's backlog.
    std::unique_ptr<ConsoleView> view_;

    /// The board as the view draws it: the live set over precursor_selection.rt_window, global.duration, the
    /// search backlog. Takes no analysis_mutex_ (the render thread calls it).
    BoardSnapshot boardSnapshot_() const;
  };
}
