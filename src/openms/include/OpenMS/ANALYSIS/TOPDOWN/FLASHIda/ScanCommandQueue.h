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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
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
   * ScanCommand values, register them in the pending scan map, and return by value.
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
    ScanCommand buildMS2(const PeakGroup& pg, int charge, const ScanConfig& scan_config, int priority, int parent_scan_id);

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
                         const MS3FragmentMatcher::ProteoformContext& proto_ctx = {});

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

    /// Decode a 3-character base-94 string back to integer
    int decode(const std::string& s) const;

    // --- Timing ---

    /// Record that an MS1 scan was just performed
    void recordMS1Time();

    /// Record that an AGC scan was just performed
    void recordAGCTime();

    /// Milliseconds since last MS1 scan (for cycle time enforcement)
    uint64_t msSinceLastMS1() const;

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

    /// Monotonically increasing tracking ID counter
    int tracking_id_counter_ = 0;

    /// Map of tracking ID -> ScanCommand for pending (in-flight) scans
    std::unordered_map<int, ScanCommand> pending_scan_map_;

    /// scan_id -> wide MS3 fragment ProForma, stashed at buildMS3 time, drained by takeMS3Proteoform for
    /// the scan_commands.tsv ms3_proteoform column. Only MS3 commands ever have an entry. Guarded by queue_mutex_.
    std::unordered_map<int, std::string> ms3_cmd_proteoform_;

    /// Timestamp of last MS1 scan (for cycle time logic)
    std::chrono::steady_clock::time_point last_ms1_time_ = std::chrono::steady_clock::now();

    /// Timestamp of last AGC scan (for AGC interval logic)
    mutable std::chrono::steady_clock::time_point last_agc_time_;

    /// All 94 printable ASCII characters (0x21-0x7E) used as tracking ID alphabet
    static const std::string tracking_alphabet_;

    /// Get next tracking ID (not thread-safe; caller must hold queue_mutex_)
    int nextTrackingIdInt_();

  };

} // namespace OpenMS
