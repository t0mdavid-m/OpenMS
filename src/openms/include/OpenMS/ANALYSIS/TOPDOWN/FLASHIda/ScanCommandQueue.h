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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanRole.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <chrono>
#include <deque>
#include <mutex>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{

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
    /// @p allowed_charges, when non-null and non-empty, is the ACQUISITION CHARGE SET the caller
    /// already resolved (PrecursorSelection::triggerAuthoredCharges). Co-isolation is confined to
    /// it. Defaulted, so callers that never restrict compile and behave exactly as before.
    /// @p marker is the scan_description type character (index 3): 'R' for an ordinary
    /// identification MS2, 'Q' for the quantification scan (ADR-0038). Defaulted to 'R' so every
    /// non-quant caller is unchanged. It is a parameter rather than something derived from
    /// @p scan_config because a scan config carries instrument parameters, never an acquisition
    /// role (ADR-0009) -- two roster entries can be byte-identical configs in different roles.
    ScanCommand buildMS2(const PeakGroup& pg, int charge, const ScanConfig& scan_config, int priority, int parent_scan_id,
                         const std::vector<int>* allowed_charges = nullptr, char marker = 'R');

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

    /// Remove expired commands from pending_scan_map_ using timeout_ms
    void cleanupExpired();

    /// Take (find + erase + return) the wide MS3 fragment ProForma stashed for @p scan_id at buildMS3 time;
    /// returns "" if absent (non-MS3 command, or already taken). Thread-safe (acquires queue_mutex_).
    std::string takeMS3Proteoform(int scan_id);

    /// Cancel commands by scan_id: remove any matching entries from the priority
    /// queues AND the pending (in-flight) map. Returns the scan_ids actually removed.
    /// Used to abort an exploration group's child scans when its baseline fails.
    std::vector<int> cancelByScanIds(const std::vector<int>& scan_ids);

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

  };

} // namespace OpenMS
