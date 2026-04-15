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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Quantification.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
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
    using TagMatch = PrecursorSelection::TagMatch;

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

    // Fragment analysis methods delegate to fragments_ component.
    // C-pointer overloads are kept as non-inline methods for binary compatibility.

    int getBestMS2Masses(int n, double* masses, double* qscores, int* charges,
                         double* window_starts, double* window_ends);

    int getTopFragmentMatches(const String& protein_sequence, int n,
                              double* masses, double* qscores, int* charges,
                              double* window_starts, double* window_ends,
                              char* ion_types, int* fragment_indices,
                              const String& fragmentation_method = "HCD");

    int getTerminalFragmentIons(const String& protein_sequence, int n,
                                double* masses, double* qscores, int* charges,
                                double* window_starts, double* window_ends,
                                char* ion_types, int* fragment_indices,
                                const String& fragmentation_method = "HCD");

    int getAmbiguityEnclosingIons(const String& protein_sequence, int n,
                                  double* masses, double* qscores, int* charges,
                                  double* window_starts, double* window_ends,
                                  char* ion_types, int* fragment_indices,
                                  const String& fragmentation_method = "HCD");

    // ──────────────────────────────────────────────────
    // Python API forwards -- delegates to components.
    // Update FLASHIda.pxd if signatures change.
    // ──────────────────────────────────────────────────
    int getBestMS2MassesPy(int n, std::vector<double>& masses, std::vector<double>& qscores,
                           std::vector<int>& charges, std::vector<double>& window_starts,
                           std::vector<double>& window_ends)
    {
      return fragments_.getBestMS2MassesPy(n, masses, qscores, charges, window_starts, window_ends,
                                            deconv_.storedMS2());
    }

    int getTopFragmentMatchesPy(const String& protein_sequence, int n,
                                std::vector<double>& masses, std::vector<double>& qscores,
                                std::vector<int>& charges, std::vector<double>& window_starts,
                                std::vector<double>& window_ends,
                                std::vector<int>& is_b_ions, std::vector<int>& fragment_indices)
    {
      return fragments_.getTopFragmentMatchesPy(protein_sequence, n, masses, qscores, charges,
                                                 window_starts, window_ends, is_b_ions,
                                                 fragment_indices, deconv_.storedMS2());
    }

    int getTerminalFragmentIonsPy(const String& protein_sequence, int n,
                                  std::vector<double>& masses, std::vector<double>& qscores,
                                  std::vector<int>& charges, std::vector<double>& window_starts,
                                  std::vector<double>& window_ends,
                                  std::vector<int>& is_b_ions, std::vector<int>& fragment_indices)
    {
      return fragments_.getTerminalFragmentIonsPy(protein_sequence, n, masses, qscores, charges,
                                                   window_starts, window_ends, is_b_ions,
                                                   fragment_indices, deconv_.storedMS2());
    }

    int getAmbiguityEnclosingIonsPy(const String& protein_sequence, int n,
                                    std::vector<double>& masses, std::vector<double>& qscores,
                                    std::vector<int>& charges, std::vector<double>& window_starts,
                                    std::vector<double>& window_ends,
                                    std::vector<int>& is_b_ions, std::vector<int>& fragment_indices)
    {
      return fragments_.getAmbiguityEnclosingIonsPy(protein_sequence, n, masses, qscores, charges,
                                                     window_starts, window_ends, is_b_ions,
                                                     fragment_indices, deconv_.storedMS2());
    }
    // ──────────────────────────────────────────────────

    /// Retrieve an integer config value by key (for bridge functions)
    int getConfigInt(const std::string& key) const;

    /// Retrieve a double config value by key (for bridge functions)
    double getConfigDouble(const std::string& key) const;

    /// Process an incoming scan and enqueue resulting commands
    int processScan(const double* mzs, const double* ints, int length,
                    double rt_min, int ms_level, const char* scan_description,
                    double faims_cv = 0.0);

    /// Dequeue the next scan command by priority. Returns 1 if command filled, 0 if error.
    int getNextScanCommand(ScanCommand& out);

    /// Get the next monotonically increasing tracking ID (thread-safe)
    int getNextTrackingId();

    /// Test-only: push a command into the priority queue (delegates to queue_)
    void pushCommandForTest(ScanCommand cmd)
    {
      queue_.push(cmd);
    }

    /// Test-only accessor: access the ScanCommandQueue directly
    ScanCommandQueue& getQueueForTest() { return queue_; }

    /**
           @brief parse FLASHIda log file
           @param in_log_file input log file
           @return parsed information : scan number - percursor information
    **/
    static std::map<int, std::vector<std::vector<float>>> parseFLASHIdaLog(const String& in_log_file);

    /**
     * @brief Process stored MS2 deconvolution for protein family detection and inclusion list expansion
     *
     * Delegates to PrecursorSelection::processMS2ForTagBasedTargeting().
     *
     * @param precursor_mass monoisotopic mass of the precursor (from iAPI)
     * @return true if target protein detected and targets expanded, false otherwise
     */
    bool processMS2ForTagBasedTargeting(double precursor_mass)
    {
      return selection_.processMS2ForTagBasedTargeting(precursor_mass);
    }

  private:
    /// Configuration object (owns all parsed config values)
    Config config_;

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

    /// Mutex protecting analysis state: deconv_, selection_, exploration_, faims_, quant_, fragments_.
    /// Also protects logging streams: writeScanResultRow_/writeIDALogEntry_ are called with this lock
    /// already held (from processScan). writeScanCommandRow_ acquires it internally (called from
    /// getNextScanCommand which does NOT hold this lock).
    mutable std::mutex analysis_mutex_;

    /// Atomic flag: true when any exploration group is active (set by processScan, read by getNextScanCommand)
    std::atomic<bool> exploration_active_{false};

    /// Atomic FAIMS CV: current CV value (set by processScan after advanceToNextCV, read by getNextScanCommand)
    std::atomic<double> current_faims_cv_{0.0};

    FAIMS faims_;   ///< FAIMS CV cycling state machine

    /// Exploration CE sweep engine (owns groups, variants, scoring)
    Exploration exploration_;

    // --- Logging file streams (append-only, crash-safe) ---
    std::ofstream ida_log_stream_;
    std::ofstream commands_tsv_stream_;
    std::ofstream results_tsv_stream_;

    /// Steady-clock reference for timestamps
    std::chrono::steady_clock::time_point engine_start_time_;

    /// Write IDA log entry for MS1 deconvolution results
    void writeIDALogEntry_(double rt, const std::string& tracking_id,
                           const std::vector<ScanCommand>& ms2_commands,
                           const DeconvolvedSpectrum& all_peak_groups);

    /// Write one TSV row for a dequeued scan command
    void writeScanCommandRow_(const ScanCommand& cmd);

    /// Write one TSV row for a processScan result
    void writeScanResultRow_(const std::string& tracking_id, double rt,
                             int mass_count, int commands_pushed,
                             const std::vector<std::string>& child_ids,
                             int tag_count, const std::string& matched_protein,
                             const std::string& proteoform_sequence,
                             uint64_t enqueue_ts, uint64_t dequeue_ts, uint64_t received_ts,
                             const DeconvolvedSpectrum* deconv_spectrum,
                             const std::string& parent_tracking_id,
                             float tic_coverage = 0.0f, int fragment_count = 0,
                             int exploration_group_id = -1, int exploration_metric = 0,
                             int variant_index = -1, int total_variants = 0,
                             double collision_energy = 0.0, double exploration_score = -1.0,
                             double remaining_ratio = -1.0);

    /// Derive scan_type string from scan_description
    static std::string scanTypeFromDescription_(const ScanCommand& cmd);

  public:
    /// Test-only: get level config for a given MSn level
    const MSLevelConfig& getLevelConfigForTest(int level) const { return config_.level(level); }

    /// Test-only: access the Config object directly
    const Config& getConfigForTest() const { return config_; }

    /// Test-only: get queue size for a given priority (delegates to queue_)
    size_t getQueueSizeForTest(int priority) const
    {
      return queue_.queueSize(priority);
    }

  };
}
