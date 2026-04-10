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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Exploration.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FAIMS.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <deque>
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

    /// Structure representing a PTM site detected by FLASHExtender
    struct PTMSite {
      int position;        ///< Position in protein sequence (1-based, midpoint)
      int start_position;  ///< Start of the region where PTM could be localized (1-based)
      int end_position;    ///< End of the region where PTM could be localized (1-based)
      double mass_shift;   ///< Observed mass shift (modification mass)
    };

    /// constructor that takes string input argument
    explicit FLASHIda(char *arg);

    /// destructor
    ~FLASHIda() = default;

    /// copy constructor
    FLASHIda(const FLASHIda& ) = default;

    /// move constructor
    FLASHIda(FLASHIda&& other) = default;

    /// assignment operator
    FLASHIda& operator=(const FLASHIda& fd) = default;

    /**
           @brief get peak groups (deconvolved masses) from input spectrum, specified by mzs and intensities (due to C# interface it is necessary)
           @param mzs mz values of the input spectrum
           @param intensities intensities of the input spectrum
           @param length length of mzs and ints
           @param rt Retention time in seconds
           @param ms_level ms level
           @param name spectrum name
           @param reporter_mz_tol reporter ion mz tolerance for isobaric quantification
           @param fold_change_threshold the threshold value for when something is considered a fold change 
           @param only_one_condition whether or not a missing condition should be considered as differentially abundant
           @return number of acquired peak groups
      */
    bool isDifferentiallyAbundant(const double* mzs,
                            const double* ints,
                            const int length,
                            const double rt,
                            const int ms_level,
                            const char* name,
                            double reporter_mz_tol,
                            double fold_change_threshold,
                            bool only_one_condition);

    /**
           @brief get peak groups (deconvolved masses) from input spectrum, specified by mzs and intensities (due to C# interface it is necessary)
           @param mzs mz values of the input spectrum
           @param intensities intensities of the input spectrum
           @param length length of mzs and ints
           @param rt Retention time in seconds
           @param ms_level ms level
           @param name spectrum name
           @param cv CV values when FAIMS is used
           @return number of acquired peak groups
      */
    int getPeakGroups(const double *mzs,
                      const double *intensities,
                      int length,
                      double rt,
                      int ms_level,
                      const char *name,
                      const char *cv)
    {
      return selection_.filterAndRank(mzs, intensities, length, rt, ms_level, cv);
    }

    /**
           @brief get isolation windows using FLASHDeconv algorithm. Many parameters are in primitive types so they can be passed to C# FLASHIda side.
           All parameters are for isolation windows.
           @param window_start window start mzs
           @param window_end window end mzs
           @param qscores QScores of windows
           @param charges charges of windows
           @param min_charges minimum charges
           @param max_charges maximum charges
           @param mono_masses monoisotopic masses
           @param charge_cos charge cosine scores
           @param charge_snrs charge SNRs or precursor SNRs
           @param iso_cos mass cosine scores
           @param snrs mass SNRs
           @param charge_scores charge distribution scores
           @param ppm_errors average PPM errors
           @param precursor_intensities precursor peak intensities
           @param peakgroup_intensities precursor mass intensities
           @param ids precursor IDs
      */
    void getIsolationWindows(double *window_start,
                             double *window_end,
                             double *qscores,
                             int *charges,
                             int *min_charges,
                             int *max_charges,
                             double *mono_masses,
                             double *charge_cos,
                             double *charge_snrs,
                             double *iso_cos,
                             double *snrs, double *charge_scores,
                             double *ppm_errors,
                             double *precursor_intensities,
                             double *peakgroup_intensities,
                             int* hcds,
                             int *ids)
    {
      selection_.getIsolationWindows(window_start, window_end, qscores, charges,
                                     min_charges, max_charges, mono_masses, charge_cos,
                                     charge_snrs, iso_cos, snrs, charge_scores,
                                     ppm_errors, precursor_intensities, peakgroup_intensities,
                                     hcds, ids);
    }
    /**
           @brief Remove a given precursor from the exclusion list by id (needed for FAIMS)
           @param id id of precursor
      */
    void removeFromExlusionList(int id)
    {
      selection_.removeFromExclusionList(id);
    }

    double getRepresentativeMass();

    void getAllMonoisotopicMasses(double *masses, int length);

    int GetAllPeakGroupSize();

    /**
     * @brief Get the top N MS2 masses with isolation window info
     *
     * Requires prior call to deconvolveMS2(). Returns masses sorted by qscore.
     *
     * @param n maximum number of masses to return
     * @param masses output: monoisotopic masses
     * @param qscores output: quality scores
     * @param charges output: representative charges
     * @param window_starts output: isolation window start m/z
     * @param window_ends output: isolation window end m/z
     * @return actual number of masses returned
     */
    int getBestMS2Masses(int n,
                         double* masses,
                         double* qscores,
                         int* charges,
                         double* window_starts,
                         double* window_ends);

    /**
     * @brief Python-friendly overload of getBestMS2Masses
     */
    int getBestMS2MassesPy(int n,
                           std::vector<double>& masses,
                           std::vector<double>& qscores,
                           std::vector<int>& charges,
                           std::vector<double>& window_starts,
                           std::vector<double>& window_ends);

    /**
     * @brief Get top fragment ion matches against a protein sequence, sorted by qscore
     *
     * Performs direct fragment matching using stored MS2 deconvolution results.
     * Uses MS2 tolerance from config (config_.level(2).tolerance_ppm).
     * Only returns peaks that match theoretical b/y fragments of the sequence.
     * Requires deconvolveMS2() to be called first.
     *
     * @param protein_sequence the protein sequence to match against
     * @param n maximum number of matches to return
     * @param masses output: observed monoisotopic masses
     * @param qscores output: qscores of matched peaks
     * @param charges output: charges of matched peaks
     * @param window_starts output: isolation window start m/z values
     * @param window_ends output: isolation window end m/z values
     * @return number of matches found (may be less than n)
     */
    int getTopFragmentMatches(const String& protein_sequence,
                              int n,
                              double* masses,
                              double* qscores,
                              int* charges,
                              double* window_starts,
                              double* window_ends,
                              char* ion_types,
                              int* fragment_indices,
                              const String& fragmentation_method = "HCD");

    /**
     * @brief Python-friendly overload of getTopFragmentMatches
     */
    int getTopFragmentMatchesPy(const String& protein_sequence,
                                int n,
                                std::vector<double>& masses,
                                std::vector<double>& qscores,
                                std::vector<int>& charges,
                                std::vector<double>& window_starts,
                                std::vector<double>& window_ends,
                                std::vector<int>& is_b_ions,
                                std::vector<int>& fragment_indices);

    /**
     * @brief Get unique fragment ions that enclose PTM ambiguity regions
     *
     * Identifies PTM sites and returns the best fragment ions that bracket
     * ambiguity regions. Returns deduplicated ions sorted by qscore.
     * Useful for MS3 targeting to resolve PTM localization.
     * Requires deconvolveMS2() to be called first.
     *
     * @param protein_sequence the protein sequence to analyze
     * @param n maximum number of ions to return
     * @param masses output: monoisotopic masses
     * @param qscores output: quality scores
     * @param charges output: representative charges
     * @param window_starts output: isolation window start m/z
     * @param window_ends output: isolation window end m/z
     * @return number of enclosing ions found
     */
    int getAmbiguityEnclosingIons(const String& protein_sequence,
                                  int n,
                                  double* masses,
                                  double* qscores,
                                  int* charges,
                                  double* window_starts,
                                  double* window_ends,
                                  char* ion_types,
                                  int* fragment_indices,
                                  const String& fragmentation_method = "HCD");

    /**
     * @brief Python-friendly overload of getAmbiguityEnclosingIons
     */
    int getAmbiguityEnclosingIonsPy(const String& protein_sequence,
                                    int n,
                                    std::vector<double>& masses,
                                    std::vector<double>& qscores,
                                    std::vector<int>& charges,
                                    std::vector<double>& window_starts,
                                    std::vector<double>& window_ends,
                                    std::vector<int>& is_b_ions,
                                    std::vector<int>& fragment_indices);

    /**
     * @brief Get terminal (innermost) fragment ions sorted by sequence position
     *
     * Returns fragment ions that extend furthest toward the opposite terminus:
     * - b-ions: sorted by highest fragment_index (rightmost, closest to C-terminus)
     * - y-ions: sorted by highest fragment_index (leftmost, closest to N-terminus)
     *
     * Output is interleaved: [top_b, top_y, 2nd_b, 2nd_y, ...]
     * Requires DeconvolveMS2() to be called first.
     *
     * @param protein_sequence The protein sequence to match against
     * @param n Maximum number of results to return
     * @param masses Output: observed monoisotopic masses
     * @param qscores Output: quality scores
     * @param charges Output: charge states
     * @param window_starts Output: isolation window start m/z
     * @param window_ends Output: isolation window end m/z
     * @param is_b_ions Output: true if b-ion, false if y-ion
     * @return Number of results returned (may be < n)
     */
    int getTerminalFragmentIons(const String& protein_sequence,
                                int n,
                                double* masses,
                                double* qscores,
                                int* charges,
                                double* window_starts,
                                double* window_ends,
                                char* ion_types,
                                int* fragment_indices,
                                const String& fragmentation_method = "HCD");

    /**
     * @brief Python-friendly overload of getTerminalFragmentIons
     *
     * @param is_b_ions Output: 1 if b-ion, 0 if y-ion (uses int for Python compatibility)
     * @param fragment_indices Output: 1-based fragment index (e.g., b3=3, y5=5)
     */
    int getTerminalFragmentIonsPy(const String& protein_sequence,
                                  int n,
                                  std::vector<double>& masses,
                                  std::vector<double>& qscores,
                                  std::vector<int>& charges,
                                  std::vector<double>& window_starts,
                                  std::vector<double>& window_ends,
                                  std::vector<int>& is_b_ions,
                                  std::vector<int>& fragment_indices);

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

    // SelectionMetric, ExplorationMetric, MSLevelConfig are now in FLASHIda/Config.h

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

    // TagMatch, InclusionTarget, TargetPTM structs moved to PrecursorSelection

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

    /// Internal struct for fragment match results from tag-based matching
    struct TagBasedFragmentMatch {
      int peak_index;           ///< Index in ms2_deconvolved_spectrum_
      double observed_mass;     ///< Monoisotopic mass
      double theoretical_mass;  ///< Theoretical mass
      double qscore;            ///< Quality score from PeakGroup
      int charge;               ///< Charge state
      int fragment_index;       ///< 1-based position in protein sequence
      char ion_type;            ///< Ion type: 'a', 'b', 'c', 'x', 'y', or 'z'
      double ppm_error;         ///< ppm error from match
    };

    /**
     * @brief Run FLASHTagger+FLASHExtender workflow to get matched fragments
     *
     * Uses the stored MS2 deconvolution (ms2_deconvolved_spectrum_) to:
     * 1. Generate sequence tags via FLASHTagger
     * 2. Validate tags match the protein sequence
     * 3. Run FLASHExtender for path-based matching
     * 4. Extract matched fragments with qscores, charges, and positions
     *
     * @param protein_sequence the protein sequence to match against
     * @param matches output: vector of matched fragments sorted by qscore descending
     * @param ptm_sites optional output: PTM sites detected by FLASHExtender
     * @return number of matches found
     */
    int runTagBasedFragmentMatching_(const String& protein_sequence,
                                      std::vector<TagBasedFragmentMatch>& matches,
                                      std::vector<PTMSite>* ptm_sites = nullptr,
                                      const String& fragmentation_method = "HCD");

    // Mass exclusion maps, targeting maps, filter/parse methods moved to PrecursorSelection

    /// Deconvolution engine (owns SpectralDeconvolution, MS1 result, MS2 result)
    Deconvolution deconv_;

    /// Precursor selection, targeting, mass exclusion (owns all selection state)
    PrecursorSelection selection_;

    std::map<double, std::vector<double>> cv_to_mass_ = {
      {-80.0, {2400.0, 2900.0}},
      {-70.0, {3500.0, 4000.0}},
      {-60.0, {4500.0, 5000.0}},
      {-50.0, {5100.0, 6500.0}},
      {-40.0, {7500.0, 10000.0}},
      {-30.0, {11000.0, 14000.0}},
      {-20.0, {12000.0, 15000.0}},
      {-10.0, {13000.0, 15500.0}},
      {-0.0, {14000.0, 16000.0}},
      {10.0, {15000.0, 16500.0}},
    };

    // parseJSONConfig_ removed: parsing logic moved to Config class

    // Phase 3 scan command queue infrastructure moved to ScanCommandQueue

    // Phase 4 command building moved to ScanCommandQueue (except selectMS3Targets_)

    /// MS3Target typedef from ScanCommandQueue (used by selectMS3Targets_)
    using MS3Target = ScanCommandQueue::MS3Target;

    /// Select MS3 fragment targets from last MS2 deconvolution
    std::vector<MS3Target> selectMS3Targets_();

    /// Process MS2 scan: tracking resolution, deconv, routing
    int processMS2Path_(const double* mzs, const double* ints, int length, double rt_min, const char* scan_desc);

    // Queue data members moved to ScanCommandQueue

    /// Mutex protecting processMS2Path_ and getNextScanCommand (replaces old queue_mutex_)
    mutable std::mutex mutex_;

    // Phase 4 config values moved to Config; AGC runtime state moved to ScanCommandQueue

    // level_configs_, exploration_enabled_, quant_enabled_, reporter_mz_tol_,
    // fold_change_threshold_, faims_cv_values_, max_cv_skip_, cv_precursor_threshold_,
    // faims_enabled_ moved to Config; FAIMS runtime state moved to FAIMS class

    FAIMS faims_;   ///< FAIMS CV cycling state machine

    /// Exploration CE sweep engine (owns groups, variants, scoring)
    Exploration exploration_;

    // Exploration structs, state, and methods moved to Exploration class

  public:
    // --- Phase 7: Test-only helpers ---

    /// Test-only: get number of active exploration groups (delegates to exploration_)
    int getActiveExplorationGroupCountForTest() const
    {
      std::lock_guard<std::mutex> lock(mutex_);
      return exploration_.activeGroupCount();
    }

    /// Test-only: get exploration group by ID (delegates to exploration_)
    Exploration::ExplorationGroup getExplorationGroupForTest(int group_id) const
    {
      std::lock_guard<std::mutex> lock(mutex_);
      return exploration_.getGroup(group_id);
    }

    /// Test-only: get level config for a given MSn level
    const MSLevelConfig& getLevelConfigForTest(int level) const { return config_.level(level); }

    /// Test-only: access the Config object directly
    const Config& getConfigForTest() const { return config_; }

    /// Test-only: get queue size for a given priority (delegates to queue_)
    size_t getQueueSizeForTest(int priority) const
    {
      return queue_.queueSize(priority);
    }

    /// Test-only: directly call exploration_.initiate() and push results
    void initiateExplorationForTest(int msn_level, double mz, double mass, int charge, double cv)
    {
      std::lock_guard<std::mutex> lock(mutex_);
      auto cmds = exploration_.initiate(msn_level, mz, mass, charge, cv, queue_);
      for (auto& c : cmds) queue_.push(c);
    }

    /// Test-only: directly call exploration_.feedResult() and push results
    void feedExplorationResultForTest(int group_id, int variant_index,
                                      const DeconvolvedSpectrum& ds, double rt)
    {
      // feedResult takes a tracking_id; for backward compat we look up
      // the variant's tracking_id from the group.
      std::lock_guard<std::mutex> lock(mutex_);
      auto group = exploration_.getGroup(group_id);
      if (variant_index < 0 || variant_index >= static_cast<int>(group.variants.size())) return;
      int tracking_id = queue_.decode(group.variants[variant_index].tracking_id);
      auto cmds = exploration_.feedResult(tracking_id, ds, rt, queue_);
      for (auto& c : cmds) queue_.push(c);
    }

    /// Test-only: access the Exploration object directly
    Exploration& getExplorationForTest() { return exploration_; }

  };
}
