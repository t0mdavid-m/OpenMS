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
// $Authors: Kyowon Jeong, Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Exploration.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cstdint>
#include <fstream>
#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Owns all FLASHIda output logging: the five append-only streams
   * (ida.log, scan_commands.tsv, scan_results.tsv, identification.tsv,
   * pooled_identification.tsv), their TSV headers, the row writers, and the
   * ida_log parser.
   *
   * Constructed from Config, which supplies ONE value — runtime().log_dir, the
   * already-resolved and already-created run folder. Each stream is
   * <log_dir>/<fixed basename> (the five k*Name constants below). An empty
   * log_dir opens nothing; see Config::RuntimeConfig and ADR-0015 for why that
   * meaning is the opposite of the authored method.json layer's.
   *
   * This class never creates a directory: the host does that before constructing
   * the engine. A log_dir that does not exist therefore leaves every stream
   * closed and every writer an early-return no-op.
   *
   * The writers operate purely on their arguments (ScanCommand /
   * DeconvolvedSpectrum / the two row descriptors) and the owned streams — they
   * hold no engine state.
   *
   * SELF-SYNCHRONISING, one mutex PER STREAM. Locking used to be the caller's
   * job, borrowed from FLASHIda::analysis_mutex_ — which meant the instrument
   * event thread waited out a whole deconvolution to write one row. It does not
   * any more, so the streams protect themselves.
   *
   * Per-stream rather than one shared lock, and that is not tidiness. The two
   * threads write DISJOINT streams: writeScanCommandRow is the sole writer of
   * commands_tsv_stream_ and is called only from getNextScanCommand (the event
   * thread); the other four are written only from processScan (the pool thread).
   * A single logger mutex would therefore couple them back together — the drain's
   * row would queue behind writeScanResultRow / writeIdentificationRow, each
   * ending in a synchronous flush, and behind writeIDALogEntry, which walks every
   * peak group of an MS1. That reintroduces exactly the stall the split removed,
   * through a different mutex, and the drain-blocking test cannot see it because
   * that test is scoped by mutex identity.
   *
   * The locks are uncontended in production precisely BECAUSE the streams are
   * disjoint. They are here so that disjointness is a fact rather than an
   * assumption: nothing anywhere guarantees the Thermo iAPI raises MsScanArrived
   * on a single thread. Pinned by FLASHIda_Logging_test's
   * concurrent_drain_writes_one_wellformed_row_per_call, which without them
   * produced 17 intact rows out of 1000.
   *
   * Headers are written in the CONSTRUCTOR, before any thread exists, so they
   * need no lock and cannot race a writer.
   */
  class OPENMS_DLLAPI IdaLogger
  {
  public:
    /// Stream basenames, joined onto runtime().log_dir. These are the wire contract with the
    /// C# golden comparer (LogGoldenComparer.cs) — changing one requires a golden recapture.
    static constexpr const char* kIdaLogName = "ida.log";
    static constexpr const char* kScanCommandsName = "scan_commands.tsv";
    static constexpr const char* kScanResultsName = "scan_results.tsv";
    static constexpr const char* kIdentificationName = "identification.tsv";
    static constexpr const char* kPooledIdentificationName = "pooled_identification.tsv";

    /// Construct from Config: open the five streams (append) under runtime().log_dir and write
    /// the four TSV headers. Opens nothing when log_dir is empty.
    explicit IdaLogger(const Config& config);

    /// copy constructor (deleted: ofstreams are non-copyable)
    IdaLogger(const IdaLogger&) = delete;

    /// move constructor (needed so the owning FLASHIda stays movable)
    IdaLogger(IdaLogger&& other) = default;

    /// assignment operator (deleted: ofstreams are non-copyable)
    IdaLogger& operator=(const IdaLogger&) = delete;

    /// One scan_results.tsv row (a pure acquisition-event log after the slim-down — no identification
    /// payload), filled by a processScan branch and written once at the bottom. Field set == the
    /// scan_results columns; defaults == the sentinels the non-applicable paths log, so a branch only
    /// assigns the fields it actually owns.
    struct ScanRowDescriptor
    {
      std::string tracking_id;
      int ms_level = 0;
      double rt = 0.0;
      int mass_count = 0;
      int commands_pushed = 0;
      std::vector<std::string> child_ids;
      uint64_t enqueue_ts = 0, dequeue_ts = 0, received_ts = 0;
      // INVARIANT: raw pointer into ENGINE-OWNED storage (selection_.deconvolvedMS1() /
      // deconv_.storedMS2() / exploration_.exploration_deconv_->storedMS2()) — all FLASHIda members
      // that outlive processScan. Exactly one ms-level branch runs per call, so the spectrum is never
      // re-deconvolved between fill and the bottom write. Do NOT re-deconvolve after a fill or it dangles.
      const DeconvolvedSpectrum* deconv_spectrum = nullptr;
      std::string parent_tracking_id;
      int exploration_group_id = -1;
      int exploration_metric = 0;
      int variant_index = -1;
      int total_variants = 0;
      std::string collision_energy = "0";
      double exploration_score = -1.0;
      double remaining_ratio = -1.0;
      std::string activation_type;
      std::string reaction_time = "0";
      // Per-scan identification YIELD. Sentinels are load-bearing and the defaults encode them:
      //   -1  no tagger ran on this spectrum   (every MS1 row, and every MS3 row)
      //    0  it ran and read nothing          (for a real protein, a meaningful negative result)
      //   >0  real count
      // Collapsing the first two to a plain 0 recreates precisely the ambiguity ADR-0012 had to add a
      // whole column (faims_enabled) to undo, once faims_cv = 0.0 turned out to mean two things.
      // tag_count is the identification count from ProteoformMatch, taken before any protein is
      // consulted -- NOT PrecursorSelection's FASTA-gated targeting return, which reports 0 both when
      // no tags existed and when tags existed but matched nothing, and so is a gate, not a measurement.
      int tag_count = -1;
      int fragment_count = -1;
      /// Matched-fragment TIC coverage. -1.0 when identification did not run, so that "ran and matched
      /// nothing" (0.0) stays distinguishable. NOTE this column also exists on identification.tsv, where
      /// it is per-ID-row; the two are written by different writers from the same source value.
      float tic_coverage = -1.0f;
      // F5: encoded id of the winning variant; "" on every non-completing / non-exploration row.
      std::string winner_tracking_id;
      // Isobaric quantification (ADR-0038). Measured on the 'Q' quantification scan and on nothing
      // else, so every other row logs the sentinels these defaults encode.
      //
      // quant_fold_change carries -1 for TWO states -- "not measured", and "measured, but one
      // condition was wholly empty so the ratio has no finite value" -- and quant_condition_means is
      // what separates them: "" means not measured, two numbers mean measured. Keeping both states
      // representable is the point: a species present in one condition and absent in the other is the
      // strongest result the experiment can produce, and collapsing it into a rejection is what the
      // (unreachable) only_one_condition flag existed to avoid.
      /// All N corrected channel intensities in getChannelInformation() order; EMPTY = not measured.
      /// Held as values and joined by the writer, the same division of labour as child_ids — the
      /// producer should not have to know the delimiter.
      std::vector<double> quant_channels;
      /// The two condition means, in quantification.conditions order (== the ratio direction);
      /// EMPTY = not measured. This is also what disambiguates quant_fold_change's -1: empty means
      /// "not measured", two values with -1 means "measured, but a condition was wholly absent so
      /// the ratio has no finite value".
      std::vector<double> quant_condition_means;
      double quant_fold_change = -1.0;    ///< mean(conditions[0]) / mean(conditions[1]); -1 = no finite ratio
      std::string quant_verdict;          ///< differential | not_differential | incomplete_channels | extraction_failed
    };

    /// One identification.tsv row. All members are held BY VALUE because the sources (info.*, a local
    /// ms2_ctx, or a reference into ms2_context_cache_) go out of scope / are erased before the bottom
    /// write. 0..N per scan; empty for MS1 and MS3 non-exploration.
    struct IdRowDescriptor
    {
      std::string tracking_id;
      int ms_level = 0;
      char scan_mode = '\0';
      Exploration::MS2Context ctx;
      FragmentAnalysis::ProteoformMatch match;
      /// Per-MS1-selection precursor identity (plain decimal); 0 when unknown (MS1, untracked id).
      int precursor_id = 0;
      /// Per-scan TIC / matched-fragment coverage (moved here from scan_results); appended LAST.
      float tic_coverage = 0.0f;
      /// C: FLASHExtender/FLASHTnT identification score of THIS scan's own match; -1.0 = the scan did NOT
      /// self-identify — it is a winner-re-matched non-winner row (fragments matched the WINNER ladder).
      /// Score presence is the sole distinguisher between self-ID rows and re-matched rows. Appended LAST.
      double flash_extender_score = -1.0;
    };

    /// One pooled_identification.tsv row: the current state of a ProteoformModel after each update.
    /// All members held by value (model state is discarded between updates).
    struct PooledModelDescriptor
    {
      /// Per-MS1-selection precursor identity (the model key; plain decimal).
      int precursor_id = 0;
      int nominal_mass = 0;
      double mono_mass = 0.0;
      std::string proforma;
      double score = -1.0;
      double coverage_pct = 0.0;
      int n_fragments = 0;
      std::vector<std::string> localized_mods;
      std::vector<std::string> ambiguous_mods;
      std::vector<int> contributing_scan_ids;
      std::vector<double> combined_masses;
      std::vector<std::string> combined_ms2_fragment_ions;  ///< ion label (ion_type+index, e.g. "b22") per entry, aligned with combined_masses
      std::vector<double> combined_measured;      ///< raw measured, aligned index-for-index with combined_masses
      std::vector<double> combined_theoretical;   ///< proteoform theoretical (MS2 & MS3), aligned
      std::vector<double> combined_diff_da;       ///< adjusted - theoretical, aligned
      std::vector<double> combined_diff_ppm;      ///< 1e6*diff_da/theoretical, aligned
      int update_index = 0;
      std::string trigger;            ///< "MS2" for the baseline row; folded fragment ion (e.g. "y7") for an MS3 update row
      std::string trigger_scan_id;    ///< Encoded tracking-id of the scan that drove this row (MS2 winner / MS3 production-or-variant)
    };

    /// Write IDA log entry for MS1 deconvolution results.
    ///
    /// The two identifiers are deliberately separate and neither substitutes for the other
    /// (ADR-0035, CONTEXT.md "Language - instrument control"):
    ///  - @p tracking_id is engine-minted. It is encoded base-94 into the "Access ID" token and is
    ///    the join key to the other four streams.
    ///  - @p instrument_scan_number is what the INSTRUMENT assigned. It becomes "MS1 Scan#", which
    ///    is what FLASHDeconv matches against the mzML native id. <= 0 = not supplied, in which case
    ///    the tracking id is written instead so the field never carries a fabricated scan.
    ///
    /// The old single `scan_number` parameter carried the tracking id under a name that claimed
    /// otherwise; that is the confusion this signature exists to remove.
    /// @p ms2_sources is INDEX-PARALLEL to @p ms2_commands: entry i is the PeakGroup command i was
    /// built from, and is what the line's `ChargeRange` reports (ADR-0035 decision 6). Passing it
    /// rather than looking the mass up in @p all_peak_groups is deliberate -- several PeakGroups
    /// routinely share one mass within a survey (ADR-0036), each with a different charge subset, so
    /// a lookup would pick one of several answers with nothing to notice. A short vector or a null
    /// entry falls back to the trigger charge, i.e. the pre-ADR-0035 output.
    void writeIDALogEntry(double rt, int tracking_id, int instrument_scan_number,
                          const std::vector<ScanCommand>& ms2_commands,
                          const std::vector<const PeakGroup*>& ms2_sources,
                          const DeconvolvedSpectrum& all_peak_groups);

    /// Write one TSV row for a dequeued scan command. @p precursor_id is the per-MS1-selection
    /// identity for this command (sourced from the engine-side tracking_id->precursor_id map);
    /// 0 for commands with no precursor (MS1 / AGC) or an untracked id. No ScanCommand ABI change.
    void writeScanCommandRow(const ScanCommand& cmd, int precursor_id = 0, const std::string& ms3_proteoform = "");

    /// Write one TSV row for a processScan result (ms_level is logged at scan_results col 1)
    void writeScanResultRow(const ScanRowDescriptor& row);

    /// Write one identification.tsv row for an MS2 or MS3 scan with matched fragments
    void writeIdentificationRow(const IdRowDescriptor& row);

    /// Write one pooled_identification.tsv row for a ProteoformModel update (trajectory log)
    void writePooledModelRow(const PooledModelDescriptor& r);

    /**
           @brief parse FLASHIda log file
           @param in_log_file input log file
           @return parsed information : scan number - percursor information
    **/
    static std::map<int, std::vector<std::vector<float>>> parseFLASHIdaLog(const String& in_log_file);

  private:
    /// Configuration object (owns the runtime output paths)
    const Config& config_;

    // --- Logging file streams (append-only, crash-safe), each with its own lock ---
    //
    // One mutex per stream, declared beside the stream it guards so the pairing cannot drift. Every
    // writer takes its own and nothing else, so these are LEAVES of the lock hierarchy
    // (analysis_mutex_ -> {queue_mutex_, precursor_map_mutex_, these}) and no writer may grow a call
    // into another writer without inverting an edge.
    std::ofstream ida_log_stream_;
    mutable std::mutex ida_log_mutex_;
    std::ofstream commands_tsv_stream_;
    mutable std::mutex commands_tsv_mutex_;
    std::ofstream results_tsv_stream_;
    mutable std::mutex results_tsv_mutex_;
    std::ofstream identification_tsv_stream_;
    mutable std::mutex identification_tsv_mutex_;
    std::ofstream pooled_stream_;
    mutable std::mutex pooled_mutex_;

    /// One-shot latch for the "no instrument scan number was supplied" warning. Guarded by
    /// ida_log_mutex_, which writeIDALogEntry already holds. ONE warning per run, deliberately: it
    /// is emitted with std::endl, i.e. a flush, under that lock -- per-scan it would be thousands of
    /// flushes racing the host's own console appender, which is exactly the pathology ADR-0025 and
    /// ADR-0033 were about.
    bool warned_missing_scan_number_ = false;

    /// Derive scan_type string from scan_description
    static std::string scanTypeFromDescription_(const ScanCommand& cmd);
  };
}
