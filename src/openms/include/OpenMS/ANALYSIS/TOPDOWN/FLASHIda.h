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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Exploration.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FAIMS.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Quantification.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iomanip>
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
    /// @param instrument_scan_number the number the INSTRUMENT assigned this scan (Header["Scan"] on
    ///        the C# side). <= 0 means "not supplied" and the ida.log falls back to the tracking id
    ///        (ADR-0035). faims_cv deliberately LOSES its default here: a defaulted trailing
    ///        parameter would let every existing call site keep compiling while silently passing 0,
    ///        which collapses every ida.log entry onto one map key with nothing to notice it.
    int processScan(const double* mzs, const double* ints, int length,
                    double rt_min, int ms_level, const char* scan_description,
                    double faims_cv, int instrument_scan_number);

    /// Dequeue the next scan command by priority. Returns 1 if command filled, 0 if error.
    int getNextScanCommand(ScanCommand& out);

    /// Get the next monotonically increasing tracking ID (thread-safe)
    int getNextTrackingId();

    /// Test-only access (push / queue / queueSize / explorationActive) lives in the test header
    /// FLASHIda_TestAccess.h via this friend, so test scaffolding stays out of the production API.
    friend struct FLASHIdaTestAccess;

  private:
    /// Configuration object (owns all parsed config values)
    Config config_;

    /// Owns all output logging: the 4 streams, TSV headers, row writers, and the ida_log parser
    IdaLogger logger_;

    /// Scan command queue (owns queue, tracking IDs, pending map, command builders)
    ScanCommandQueue queue_;

    /// Deconvolution engine (owns SpectralDeconvolution, MS1 result, MS2 result)
    Deconvolution deconv_;

    /// Fragment analysis (tag-based matching, MS2 mass queries, PTM ambiguity)
    FragmentAnalysis fragments_;

    /// Precursor selection, targeting, mass exclusion (owns all selection state)
    PrecursorSelection selection_;

    /// Isobaric quantification (TMT reporter-ion differential abundance test)
    Quantification quant_;

    /// Mutex protecting analysis state: deconv_, selection_, exploration_, faims_, quant_, fragments_,
    /// tracker_, ms2_context_cache_ and next_precursor_id_.
    ///
    /// ONE acquisition site: processScan (function-scoped, so it is held across the whole
    /// deconvolution). getNextScanCommand does NOT acquire it and must not start -- the instrument
    /// event thread calls the drain, and making it wait on a deconvolution to write one TSV row and
    /// read one map entry is what this arrangement exists to prevent.
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

    /// Leaf lock for precursor_id_by_tracking_.
    ///
    /// Lock hierarchy: analysis_mutex_ -> { queue_mutex_, precursor_map_mutex_, the IdaLogger stream
    /// mutexes }, all leaves. This one is taken ONLY inside the two accessors below, and neither of
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

    /// Steady-clock reference for timestamps
    std::chrono::steady_clock::time_point engine_start_time_;

  };
}
