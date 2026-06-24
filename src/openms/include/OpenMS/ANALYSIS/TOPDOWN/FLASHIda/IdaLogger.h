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
#include <string>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Owns all FLASHIda output logging: the four append-only streams
   * (ida_log, scan_commands.tsv, scan_results.tsv, identification.tsv), their
   * TSV headers, the row writers, and the ida_log parser.
   *
   * Constructed from Config (reads the four runtime() paths and opens + headers
   * the streams). The writers operate purely on their arguments (ScanCommand /
   * DeconvolvedSpectrum / the two row descriptors) and the owned streams — they
   * hold no engine state. Locking is the CALLER's responsibility: FLASHIda holds
   * analysis_mutex_ across writeScanResultRow/writeIdentificationRow/writeIDALogEntry
   * (called from processScan), and wraps writeScanCommandRow in the same lock at
   * its getNextScanCommand call sites. IdaLogger itself is lock-agnostic.
   */
  class OPENMS_DLLAPI IdaLogger
  {
  public:
    /// Construct from Config: open the four streams (append) and write the three TSV headers.
    explicit IdaLogger(const Config& config);

    /// copy constructor (deleted: ofstreams are non-copyable)
    IdaLogger(const IdaLogger&) = delete;

    /// move constructor (needed so the owning FLASHIda stays movable)
    IdaLogger(IdaLogger&& other) = default;

    /// assignment operator (deleted: ofstreams are non-copyable)
    IdaLogger& operator=(const IdaLogger&) = delete;

    /// One scan_results.tsv row, filled by a processScan branch and written once at the bottom.
    /// Field set == the scan_results columns; defaults == the sentinels the non-applicable paths log,
    /// so a branch only assigns the fields it actually owns. NOTE: tag_count/fragment_count default to
    /// the MS3 sentinel (-1) — MS1 must set both to 0 and MS2 to its real counts.
    struct ScanRowDescriptor
    {
      std::string tracking_id;
      int ms_level = 0;
      double rt = 0.0;
      int mass_count = 0;
      int commands_pushed = 0;
      std::vector<std::string> child_ids;
      int tag_count = -1;
      std::string matched_protein;
      std::string proteoform_sequence;
      uint64_t enqueue_ts = 0, dequeue_ts = 0, received_ts = 0;
      // INVARIANT: raw pointer into ENGINE-OWNED storage (selection_.deconvolvedMS1() /
      // deconv_.storedMS2() / exploration_.exploration_deconv_->storedMS2()) — all FLASHIda members
      // that outlive processScan. Exactly one ms-level branch runs per call, so the spectrum is never
      // re-deconvolved between fill and the bottom write. Do NOT re-deconvolve after a fill or it dangles.
      const DeconvolvedSpectrum* deconv_spectrum = nullptr;
      std::string parent_tracking_id;
      float tic_coverage = 0.0f;
      int fragment_count = -1;
      int exploration_group_id = -1;
      int exploration_metric = 0;
      int variant_index = -1;
      int total_variants = 0;
      std::string collision_energy = "0";
      double exploration_score = -1.0;
      double remaining_ratio = -1.0;
      std::string activation_type;
      std::string reaction_time = "0";
      // F5: encoded id of the winning variant; "" on every non-completing / non-exploration row.
      std::string winner_tracking_id;
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
    };

    /// One pooled_identification.tsv row: the current state of a ProteoformModel after each update.
    /// All members held by value (model state is discarded between updates).
    struct PooledModelDescriptor
    {
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
      int update_index = 0;
    };

    /// Write IDA log entry for MS1 deconvolution results
    void writeIDALogEntry(double rt, int scan_number, const std::string& tracking_id,
                          const std::vector<ScanCommand>& ms2_commands,
                          const DeconvolvedSpectrum& all_peak_groups);

    /// Write one TSV row for a dequeued scan command
    void writeScanCommandRow(const ScanCommand& cmd);

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

    // --- Logging file streams (append-only, crash-safe) ---
    std::ofstream ida_log_stream_;
    std::ofstream commands_tsv_stream_;
    std::ofstream results_tsv_stream_;
    std::ofstream identification_tsv_stream_;
    std::ofstream pooled_stream_;

    /// Derive scan_type string from scan_description
    static std::string scanTypeFromDescription_(const ScanCommand& cmd);
  };
}
