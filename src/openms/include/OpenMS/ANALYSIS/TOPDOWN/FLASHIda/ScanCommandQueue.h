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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms2Params.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/NotchSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PlanUnit.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanRole.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <chrono>
#include <deque>
#include <map>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OpenMS
{

  /// A bought R rewritten in place by a newer survey (it is adjustable, never droppable).
  struct AdjustRequest
  {
    int scan_id = 0;
    ScanCommand rewritten{};
    int decision_survey_id = 0;
    double isolation_purity = -1.0;
  };

  /// What a survey decided a command on, written when the command is planned and taken by the drain.
  struct DecisionInfo
  {
    int decision_survey_id = 0;
    double isolation_purity = -1.0;
  };

  /// What one replacePlan call did. The caller closes / drops the groups AFTER the queue lock is released.
  struct ReplanResult
  {
    std::vector<int> dropped_scan_ids;     ///< for ida.log Superseded= and forgetTracking_
    std::vector<int> dropped_quant_groups; ///< caller: quant_.closeGroup
    std::vector<int> dropped_sweep_groups; ///< caller: exploration_.dropGroup
    std::vector<int> reconciled_unit_ids;  ///< new units refused: species committed meanwhile, or a remainder with nothing to replace
    std::vector<int> reconciled_scan_ids;  ///< every command reconcile kept out of the queue (forgetTracking_)
    std::vector<int> reconciled_quant_groups; ///< groups of the refused units: caller closes them
    std::vector<int> reconciled_sweep_groups; ///< groups of the refused units: caller drops them
    int kept = 0;     ///< committed units at the CV whose queued commands stay
    int added = 0;    ///< units inserted
    int adjusted = 0; ///< adjustable commands rewritten
  };

  /**
   * @brief Manages scan command building, priority queuing, and tracking ID encoding for FLASHIda.
   *
   * Fully thread-safe: every public method that touches mutable state acquires queue_mutex_ internally.
   * Callers never need to hold an external lock. Building methods (buildMS2, buildMS3, etc.) produce
   * ScanCommand values and return them by value; they do NOT enqueue and do NOT register in the
   * pending scan map. Registration happens in dequeue() and in explicit registerPending() calls,
   * which FLASHIda.cpp — the sole enqueuer — makes itself on the drain path.
   *
   * Methods that do NOT acquire queue_mutex_ (safe without locking):
   * - makeMS1(), makeAGC(): const, only read config_ (immutable after construction)
   * - encode(): static pure function
   * - decode(): reads static const tracking_alphabet_ only
   */
  class OPENMS_DLLAPI ScanCommandQueue
  {
  public:
    /// Construct with a reference to the shared Config
    explicit ScanCommandQueue(const Config& config);

    // --- Building ---

    /// Build MS2 ScanCommand from a PeakGroup + ScanConfig (unified factory).
    /// parent_scan_id is MANDATORY (0 = no parent / root); stamped onto the command at creation.
    /// @p acquisition_set, when non-null and non-empty, is the ACQUISITION CHARGE SET admission
    /// recorded (PrecursorSelection::triggerAcquisitionSets): under Multiplexed exactly its members
    /// become the notches, in its order. Without one the notches are selected here from config with
    /// the same rule and cap (ADR-0059).
    /// @p marker is the scan_description type character (index 3): 'R' for an ordinary
    /// identification MS2, 'Q' for the quantification scan (ADR-0038). Defaulted to 'R' so every
    /// non-quant caller is unchanged. It is a parameter rather than something derived from
    /// @p scan_config because a scan config carries instrument parameters, never an acquisition
    /// role (ADR-0009) -- two roster entries can be byte-identical configs in different roles.
    ScanCommand buildMS2(const PeakGroup& pg, int charge, const ScanConfig& scan_config, int priority, int parent_scan_id,
                         const std::vector<int>* acquisition_set = nullptr, char marker = 'R');

    /// Re-isolate @p r on a newer survey (survey re-planning: the bought R a newer survey adjusts). A copy of
    /// @p r -- same scan id, parent, scan parameters and FAIMS CV -- whose stage-0 anchor is `picks[0]` and whose
    /// stage-0 notches are `picks[1..]`, each measured on the PeakGroup that resolved it with the geometry
    /// buildMS2 uses (@p margin is optimal_window_margin_). The precursor fields and the scan_description are
    /// rewritten from the anchor's PeakGroup, keeping @p r's 3-char id: a stale `@z` would otherwise reach the
    /// instrument and every log. Static and lock-free -- it mints nothing. Empty @p picks returns @p r unchanged.
    static ScanCommand reisolate(const ScanCommand& r, const std::vector<std::pair<int, const PeakGroup*>>& picks, double margin,
                                 char marker);

    /// Build MS3 ScanCommand from MS2 context + fragment target + MS3 config for CE/activation.
    /// parent_scan_id is MANDATORY (the immediate parent's tracking id); stamped at creation.
    /// stage0_params: if non-null, overrides the stage[0] CE/activation/reaction_time copied from
    /// ms2_ctx with per-ion optimised MS2 parameters (caller wiring is a later task).
    /// @p ms3_proteoform: the wide clipped b/y fragment ProForma of this MS3 target, rendered by the
    /// caller (which holds the ProteoformContext). Stashed in a scan_id-keyed side-map and surfaced on the
    /// scan_commands.tsv row via takeMS3Proteoform at drain time. "" = not applicable (leaves the log cell empty).
    ScanCommand buildMS3(const ScanCommand& ms2_ctx, const ScanConfig& ms3_config,
                         double frag_mz, int frag_charge, double iso_width, int parent_scan_id,
                         char ion_type = '\0', int frag_index = 0, int priority = 1,
                         const FragmentAnalysis::FragmentScores& frag_scores = {},
                         const Ms2Params* stage0_params = nullptr,
                         const MS3FragmentMatcher::ProteoformContext& proto_ctx = {},
                         /// EXTRA co-isolated charge states of the fragment, for
                         /// characterization.fragment_charges == Multiplexed. Selected by the tracker
                         /// (Ms3Target::notches) because this builder receives only scalars and cannot
                         /// see the model's per-charge observations. May be truncated: stage-0's
                         /// inherited notches are written first and share the same slot budget.
                         const std::vector<NotchCandidate>* stage1_notches = nullptr);

    /// Create an MS1 survey scan command from current config
    ScanCommand makeMS1() const;

    /// Create an AGC calibration scan command
    ScanCommand makeAGC() const;

    /// Build follow-up MS2 at priority 0 using the given scan config and description suffix.
    /// Returns the command; caller pushes.
    ScanCommand buildFollowUp(const ScanCommand& ctx, const ScanConfig& follow_up_config, char suffix, int priority = 0);

    // --- Queue ops (all thread-safe) ---

    /// Push command into appropriate priority queue
    void push(ScanCommand cmd);

    /// Dequeue next command by priority (0 = highest). Returns nullopt if all queues empty.
    std::optional<ScanCommand> dequeue();

    /// Push @p cmd and take the next command in ONE lock acquisition.
    ///
    /// push() and dequeue() are each thread-safe; the SEQUENCE is not. Between them the lock is
    /// released, so a concurrent drain's dequeue() can take the command this thread just pushed.
    /// Harmless at the Step 2 cycle-time push, which simply takes whatever else is queued -- but at
    /// the Step 5 idle survey the caller falls back to logging its own copy, so the command is
    /// emitted TWICE while the thief, having consumed a command it never built, never allocates a
    /// tracking id at all. That is a duplicate in the tracking-id channel, which is the join key
    /// from scan_commands.tsv to the other four streams (ADR-0008).
    /// Measured on gcc/Linux: 19-42 duplicate ids per 1000 concurrent drains, at EVERY core count
    /// including --cpus 1. Never observed on MSVC in ~45 CI runs since 2026-08-14.
    std::optional<ScanCommand> pushAndDequeue(ScanCommand cmd);

    /// Register a bypass command (AGC, etc.) in pending_scan_map_ without queuing it
    void registerPending(const ScanCommand& cmd);

    /// Look up and remove a pending scan by tracking ID. Returns nullopt if not found.
    std::optional<ScanCommand> resolvePending(int id);

    /// Look up a pending command by tracking ID without removing it
    std::optional<ScanCommand> peekPending(int id) const
    {
      std::lock_guard<std::mutex> lock(queue_mutex_);
      auto it = pending_scan_map_.find(id);
      if (it == pending_scan_map_.end()) return std::nullopt;
      return it->second;
    }

    /// Remove commands that have waited in the queued deques longer than timeout_ms. A plan unit's expired
    /// member leaves the registry; an uncommitted unit that empties is discarded with its stamps.
    void cleanupExpired();

    /// Take (find + erase + return) the wide MS3 fragment ProForma stashed for @p scan_id at buildMS3 time;
    /// returns "" if absent (non-MS3 command, or already taken). Thread-safe (acquires queue_mutex_).
    std::string takeMS3Proteoform(int scan_id);

    /// Cancel commands by scan_id: remove any matching entries from the priority
    /// queues AND the pending (in-flight) map. Returns the scan_ids actually removed.
    /// No production caller: the baseline abort it was written for went with the baseline (ADR-0048).
    std::vector<int> cancelByScanIds(const std::vector<int>& scan_ids);

    // --- Survey re-planning: the plan registry (all thread-safe) ---
    //
    // A survey's MS2 decisions stay provisional while they wait here. A survey hands its whole plan for one CV
    // to replacePlan, which replaces that CV's still-UNCOMMITTED units; the drain's dequeue commits a unit and
    // moves its stamps to an inbox the analysis thread folds (takeCommitted). A committed unit is never dropped.
    // Commands no survey planned (follow-ups, MS3, MS1-family) never enter the registry and are never touched.

    /// Atomically: drop every UNCOMMITTED unit at @p cv; insert @p units (minus reconciled ones) at the back of their
    /// priority deques; apply @p adjust in place. Stamps enqueue_timestamp_ms on everything inserted or adjusted.
    /// A new uncommitted unit is reconciled away when the stamps committed since the last fold would have barred
    /// its anchor (wouldBar). A unit arriving already committed is a re-planned remainder: it inherits the
    /// sent_charges the unit it replaces (same nominal mass) holds AT THE SWAP, a command of it whose charge is among
    /// them is not inserted, and the whole unit is reconciled away when that unit is gone or nothing of it is left.
    /// Every reconciled unit and command is reported in the ReplanResult.
    /// @p replace_remainder_of: committed units whose QUEUED commands are dropped because @p units carries their
    /// re-planned remainder (non-quant `separate` R only).
    /// An empty unit is skipped. A unit's unit_id is kept when it is > 0 and not in use, otherwise one is minted.
    ReplanResult replacePlan(double cv, std::vector<PlanUnit> units, const std::vector<AdjustRequest>& adjust,
                             const std::vector<int>& replace_remainder_of);

    /// Register a queued bought R as adjustable (never droppable). @p resolved_charges = charges resolved on the
    /// Q's survey, @p q_survey_id = that survey, which stays its decision survey until an adjust says otherwise.
    void registerAdjustable(int scan_id, std::vector<int> resolved_charges, int q_survey_id);

    /// Commands of committed units at a CV that are still queued, grouped by unit, and the adjustable bought Rs
    /// still queued with their resolved charges. Snapshot copies.
    struct Waiting
    {
      std::vector<PlanUnit> committed_units;
      std::vector<std::pair<ScanCommand, std::vector<int>>> adjustable;
    };
    Waiting waitingAt(double cv) const;

    /// Move out every stamp committed since the last call, in commit order (the analysis thread folds them).
    std::vector<MemoryStamp> takeCommitted();

    /// Take (find + erase) the decision info written for @p scan_id; nullopt for commands no survey decided.
    std::optional<DecisionInfo> takeDecisionInfo(int scan_id);

    /// Check if an AGC scan is needed based on agc_interval_ms
    bool needsAGC() const;

    // --- Tracking IDs (thread-safe) ---

    /// Get next monotonically increasing tracking ID (thread-safe)
    int nextTrackingId();

    /// Encode an integer as a 3-character base-94 string (all printable ASCII 0x21-0x7E)
    static std::string encode(int value);

    /// Operator-facing scan_description mass token: the kDa mass with the most significant digits that
    /// fit the 15-char scan_description budget, plus a 'k' (kDa) unit suffix. Reserves room for the
    /// trailing ion (so it is never truncated) given the actual charge + ion fields of the command.
    static std::string formatMassToken(double mass_kda, int charge, char ion_type, int ion_index);

    /**
     * @brief Write a bare MS1-family scan_description: `{3-char tracking id}{role marker}`.
     *
     * Requires @p cmd.scan_id to be assigned already. Funnels the four sites that used to spell the
     * marker as a literal ("%sS" / "%sA"), so the marker is sourced from ScanRole.h's table rather
     * than from four independent string literals (ADR-0042).
     *
     * DELIBERATELY NOT used for the COMPOUND descriptions -- buildMS2 / buildMS3 / buildFollowUp and
     * Exploration's two variant writers append a mass token, charge and ion, and scan_description is
     * a golden column, so funnelling those would risk moving 28 files for no gain.
     *
     * Also not used for makeMS1()/makeAGC()'s bare "S"/"A" placeholders: those run before a scan_id
     * exists and are always overwritten by one of this function's call sites.
     */
    static void writeMs1Description(ScanCommand& cmd, ScanRole role);

    /// Decode a 3-character base-94 string back to integer
    int decode(const std::string& s) const;

    // --- Timing ---

    /// Record that an MS1 scan was just performed
    void recordMS1Time();

    /// Record that an AGC scan was just performed
    void recordAGCTime();

    /// Milliseconds since last MS1 scan (for cycle time enforcement)
    uint64_t msSinceLastMS1() const;

    /**
     * @brief Is a MONITOR scan due for one exploration level right now? (ADR-0042)
     *
     * A CONSUME, not a peek: it tests and updates the level's clock under ONE queue_mutex_
     * acquisition, so two concurrent drains can never both be told the same monitor scan is due.
     * Call it exactly once per level per drain and act on the answer.
     *
     * MUST BE CALLED ON EVERY DRAIN, including drains where nothing is sweeping -- the
     * `sweeping == false` path is what RE-ARMS the level for its next episode. Skip it and the arm
     * stays set from the previous sweep, so the next sweep silently loses its anchor scan.
     *
     * Cadence: one anchor the moment a level starts sweeping (`armed` false), then one every
     * @p interval_ms while it keeps sweeping.
     *
     * @param idx 0 for level 2, 1 for level 3
     * @param sweeping is that level currently sweeping AND configured to monitor?
     * @param interval_ms authored spacing; Config::validate guarantees > 0 whenever enabled
     */
    bool takeMonitorDue(int idx, bool sweeping, double interval_ms);

    /// Test-only: back-date a level's monitor clock so a wall-clock cadence is deterministic.
    /// There is no production caller and must not be one.
    void backdateMonitorClock(int idx, double ms);

    // --- Introspection ---

    /// Number of entries in pending_scan_map_
    size_t pendingScanMapSize() const;

    /// Queue size for a given priority level (0-3)
    size_t queueSize(int priority) const;

  private:
    const Config& config_;

    /// Priority queues: index 0 = highest priority, 3 = lowest
    std::deque<ScanCommand> queues_[4];

    /// Mutex protecting queues_, tracking_id_counter_, pending_scan_map_
    mutable std::mutex queue_mutex_;

    /// Monotonically increasing tracking ID counter.
    ///
    /// **Starts at 1, and 0 is never issued.** 0 is the "no parent / root" sentinel that buildMS2
    /// tests for (`if (parent_scan_id > 0)`) and that FLASHIda.cpp passes literally for a root MS2,
    /// so a command actually HOLDING id 0 is indistinguishable from one with no parent — its
    /// children ship with an empty parent_tracking_id and the lineage graph loses a generation.
    ///
    /// This was latent for the whole life of the code: the drained-queue AGC prescan was the first
    /// command minted on every fresh engine, so it absorbed id 0 and no survey ever held it.
    /// Deleting that prescan (ADR-0031) handed 0 to the first survey and broke MS2 parentage on the
    /// instrument as well as in tests. ADR-0008 independently reserves 0 on the sibling identity
    /// channel, for the same "it looks like an absence" reason.
    ///
    /// Pinned by FLASHIda_ProcessScan_test::tracking_id_zero_is_never_issued, which covers the wrap
    /// as well as the fresh engine.
    int tracking_id_counter_ = 1;

    /// Map of tracking ID -> ScanCommand for pending (in-flight) scans
    std::unordered_map<int, ScanCommand> pending_scan_map_;

    /// scan_id -> wide MS3 fragment ProForma, stashed at buildMS3 time, drained by takeMS3Proteoform for
    /// the scan_commands.tsv ms3_proteoform column. Only MS3 commands ever have an entry. Guarded by queue_mutex_.
    std::unordered_map<int, std::string> ms3_cmd_proteoform_;

    // --- Survey re-planning: the plan registry. All guarded by queue_mutex_. ---

    /// A bought R a newer survey may rewrite: the charges its Q's survey resolved, and that survey.
    struct AdjustableR
    {
      std::vector<int> resolved_charges;
      int q_survey_id = 0;
    };

    /// unit_id -> unit. A unit leaves when its last queued command is dequeued, dropped or expired.
    std::map<int, PlanUnit> units_;

    /// scan_id -> unit_id, for registered commands still QUEUED (erased on dequeue, drop and expiry).
    std::unordered_map<int, int> unit_of_;

    /// scan_id -> adjustable bought R, while it is queued.
    std::unordered_map<int, AdjustableR> adjustable_;

    /// Stamps committed by dequeue_ since the last takeCommitted, in COMMIT order -- wouldBar takes the last
    /// MassQscore by position, so neither may ever be sorted or merged out of order.
    std::vector<MemoryStamp> inbox_;
    std::vector<MemoryStamp> committed_since_fold_;

    /// scan_id -> what the survey decided it on, written by replacePlan / registerAdjustable, taken by the drain.
    std::unordered_map<int, DecisionInfo> decision_info_;

    /// Last unit id minted by replacePlan.
    int unit_id_counter_ = 0;

    /// Timestamp of last MS1 scan (for cycle time logic)
    std::chrono::steady_clock::time_point last_ms1_time_ = std::chrono::steady_clock::now();

    /// Monitor-scan cadence state, per exploration level ([0] = level 2, [1] = level 3). ADR-0042.
    /// `armed` false means "this level's next drain owes an EPISODE ANCHOR"; see takeMonitorDue.
    std::chrono::steady_clock::time_point last_monitor_time_[2];
    bool monitor_armed_[2] = {false, false};

    /// Timestamp of last AGC scan (for AGC interval logic)
    mutable std::chrono::steady_clock::time_point last_agc_time_;

    /// All 94 printable ASCII characters (0x21-0x7E) used as tracking ID alphabet
    static const std::string tracking_alphabet_;

    /// Get next tracking ID (not thread-safe; caller must hold queue_mutex_)
    int nextTrackingIdInt_();

    /// push()/dequeue() bodies WITHOUT the lock; caller must hold queue_mutex_.
    /// Trailing underscore = unlocked, the same convention nextTrackingIdInt_() uses.
    /// They exist so pushAndDequeue() can run both under a single acquisition -- queue_mutex_ is
    /// non-recursive, so it could not simply call the public wrappers.
    void push_(ScanCommand cmd);
    std::optional<ScanCommand> dequeue_();

    /// Remove the registered, QUEUED command @p scan_id from the registry (not from the deques): its unit loses a
    /// command and is erased when none is left. Caller must hold queue_mutex_.
    void forgetQueued_(int scan_id);

    /// Every charge @p cmd isolates at its precursor stage: the anchor first, then stage 0's notches (multiplexed).
    static std::vector<int> isolatedCharges_(const ScanCommand& cmd);

    /// buildMS2's stage-0 anchor geometry: the centre of @p pg's measured m/z range at @p charge, that range
    /// widened by @p margin on each side, and the charge. Shared with reisolate.
    static void applyIsolation_(ScanCommand& cmd, const PeakGroup& pg, int charge, double margin);

    /// The notches for an explicit charge list: each member of @p charges but @p anchor_charge, in @p charges
    /// order, taken from @p candidates; a charge absent from them is skipped.
    static std::vector<NotchCandidate> notchesFor_(const std::vector<NotchCandidate>& candidates, int anchor_charge,
                                                   const std::vector<int>& charges);

    /// buildMS2's precursor scoring fields, from @p pg at @p charge.
    static void applyPrecursorScores_(ScanCommand& cmd, const PeakGroup& pg, int charge);

    /// buildMS2's scan_description: `{id}{marker}{mass token}@{charge}`.
    static void writeMs2Description_(ScanCommand& cmd, const std::string& id_str, char marker, double mono_mass, int charge);

    /// Test scaffolding reads the plan registry through this friend rather than through a public getter.
    friend struct FLASHIdaTestAccess;
  };

} // namespace OpenMS
