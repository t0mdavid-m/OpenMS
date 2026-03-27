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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <set>

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
                      const char *cv);

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
                             int *ids);
    /**
           @brief Remove a given precursor from the exclusion list by id (needed for FAIMS)
           @param id id of precursor
      */
    void removeFromExlusionList(int id);

    double getRepresentativeMass();

    void getAllMonoisotopicMasses(double *masses, int length);

    int GetAllPeakGroupSize();

    /**
     * @brief Deconvolve an MS2 spectrum and store the result for subsequent operations
     *
     * After calling this, getSequenceTagsAndMatches() and identifyProteoform()
     * will automatically use the stored deconvolution instead of re-deconvolving.
     *
     * @param mzs m/z values of the input spectrum
     * @param ints intensities of the input spectrum
     * @param length number of peaks
     * @param rt retention time in seconds
     * @param precursor_mass precursor monoisotopic mass (from MS1 deconvolution), <= 0 means no precursor
     * @param precursor_charge precursor charge state (positive for positive mode), 0 means no precursor
     * @return number of peak groups found
     */
    int deconvolveMS2(const double* mzs,
                      const double* ints,
                      int length,
                      double rt,
                      double precursor_mass,
                      int precursor_charge);

    /**
     * @brief Python-friendly overload of deconvolveMS2
     */
    int deconvolveMS2Py(const std::vector<double>& mzs,
                        const std::vector<double>& ints,
                        double rt,
                        double precursor_mass,
                        int precursor_charge);

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
     * @brief Check if MS2 deconvolution data is available
     */
    bool hasMS2Deconvolution() const;

    /**
     * @brief Get number of peak groups in stored MS2 deconvolution
     */
    int getMS2PeakGroupCount() const;

    /**
     * @brief Clear stored MS2 deconvolution
     */
    void clearMS2Deconvolution();

    /**
     * @brief Get top fragment ion matches against a protein sequence, sorted by qscore
     *
     * Performs direct fragment matching using stored MS2 deconvolution results.
     * Uses MS2 tolerance from constructor (tol_[1]).
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

    /**
           @brief parse FLASHIda log file
           @param in_log_file input log file
           @return parsed information : scan number - percursor information
    **/
    static std::map<int, std::vector<std::vector<float>>> parseFLASHIdaLog(const String& in_log_file);

    /**
     * @brief Structure to hold a sequence tag match result
     */
    struct TagMatch
    {
      String tag_sequence;        ///< The sequence tag string
      double n_term_mass;         ///< N-terminal flanking mass
      double c_term_mass;         ///< C-terminal flanking mass
      double tag_score;           ///< Score of the tag
      int protein_index;          ///< Index in the FASTA file
      String protein_accession;   ///< Protein accession
      int match_position;         ///< Position in the protein sequence where tag matches
      double flanking_mass_diff;  ///< Difference between tag flanking mass and protein flanking mass
    };

    /**
     * @brief Structure for TSV-based inclusion target
     *
     * Represents a single target from a TSV inclusion list with mass, charge,
     * RT range, and priority for tie-breaking during precursor selection.
     */
    struct InclusionTarget
    {
      double mass;        ///< Target monoisotopic mass
      int charge;         ///< Charge state (-1 = any charge)
      double rt_start;    ///< RT window start (seconds)
      double rt_end;      ///< RT window end (seconds)
      int priority;       ///< Tie-breaking priority (higher = preferred)

      /// Check if current RT is within the target's active window
      bool isActiveAt(double rt) const { return rt >= rt_start && rt <= rt_end; }

      /// Check if charge matches (true if target charge is -1 or matches)
      bool matchesCharge(int c) const { return charge == -1 || charge == c; }
    };

    /**
     * @brief Structure for PTM modifications with combinatorial limits for target expansion
     *
     * Used by tag-based targeting to generate PTM mass combinations for dynamic inclusion lists.
     */
    struct TargetPTM
    {
      String name;      ///< Modification name (e.g., "acetylation")
      double mass;      ///< Delta mass in Da
      int max_count;    ///< Maximum occurrences on single proteoform
    };

    /**
     * @brief Process stored MS2 deconvolution for protein family detection and inclusion list expansion
     *
     * This DLL bridge function performs real-time tag-based targeting:
     * 1. Uses stored MS2 deconvolution results (requires deconvolveMS2() first)
     * 2. Extracts sequence tags (minimum length 3)
     * 3. Matches tags against the target protein database
     * 4. If match found: Expands target masses using PTM combinations from precursor mass
     * 5. Adds expanded masses to dynamic inclusion list
     *
     * @param precursor_mass monoisotopic mass of the precursor (from iAPI)
     * @return true if target protein detected and targets expanded, false otherwise
     */
    bool processMS2ForTagBasedTargeting(double precursor_mass);

  private:
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

    /// PeakGroup comparator for soring by QScore
    /*struct
    {
      bool operator()(const PeakGroup& a, const PeakGroup& b) const
      {
        return a.getQScore() > b.getQScore();
      }
    } QscoreComparator_;
*/
    /// Maps that are necessary for mass exclusion
    std::unordered_map<int, double> tqscore_exceeding_mz_rt_map_; /// integer mz value vs. retention time with tqscore exceeding total qscore threshold
    std::unordered_map<int, double> tqscore_exceeding_mass_rt_map_; /// integer mass value vs. retention time with tqscore exceeding total qscore threshold
    std::unordered_map<int, double> all_mass_rt_map_; /// mz value vs. retention time for all acquired precursors
    std::unordered_map<int, double> mass_qscore_map_; /// mass value vs. total qscore for all acquired precursors



    /// Maps that are neccessary for selectively disabling mass exclusion (needed for FAIMS support)
    std::unordered_map<int, int> id_mass_map_;
    std::unordered_map<int, int> id_mz_map_;
    std::unordered_map<int, double> id_qscore_map_;


    /**
         @brief discard peak groups using mass exclusion
         @param ms_level MS level
         @param rt Retention time
    */
    void filterPeakGroupsUsingMassExclusion_(int ms_level, double rt);

    /**
         @brief parse TSV inclusion list file
         @param filename path to TSV file with columns: mass, charge, rt_start, rt_end, priority
    */
    void parseInclusionListTSV_(const String& filename);

    /**
         @brief parse TSV file containing PTM modifications for target expansion
         @param filename path to TSV file with columns: name, mass, max_count
    */
    void parseTargetPTMsTSV_(const String& filename);

    /**
         @brief Generate all PTM mass combinations for a base protein mass
         @param base_mass base protein mass
         @param ptms vector of PTMs with max counts
         @return vector of all possible modified masses
    */
    std::vector<double> generatePTMCombinations_(double base_mass,
                                                 const std::vector<TargetPTM>& ptms) const;

    /**
         @brief Add expanded target masses to the dynamic inclusion list
         @param masses vector of masses to add
         @param rt current retention time (targets active for rt_window_ seconds)
         @param priority priority for the new targets
    */
    void addDynamicTargets_(const std::vector<double>& masses,
                            double rt,
                            int priority);

    /**
         @brief generate MSSpectrum class using mzs and intensities. mzs and intensities and other information are
         provided by Thermo iAPI
         @param mzs m/z values
         @param ints intensities
         @param length number of peaks
         @param rt Retention time
         @param ms_level MS level
         @param name spectrum name
    */
    static MSSpectrum makeMSSpectrum_(const double *mzs,
                                      const double *ints,
                                      int length,
                                      double rt,
                                      int ms_level,
                                      const char *name);

    /// deconvolved spectrum that contains the peak group
    DeconvolvedSpectrum deconvolved_spectrum_;
    /// selected peak groups out of deconvolved_spectrum_
    DeconvolvedSpectrum selected_peak_groups_;

    /// MS2 deconvolved spectrum storage (mirrors deconvolved_spectrum_ for MS1)
    DeconvolvedSpectrum ms2_deconvolved_spectrum_;
    /// Retention time of stored MS2 deconvolution (for validation)
    double ms2_deconv_rt_ = -1.0;
    /// Flag indicating if MS2 deconvolution is valid
    bool ms2_deconv_valid_ = false;

    /// peakGroup charges to be triggered
    std::vector<int> trigger_charges;
    std::vector<int> trigger_hcds;
    std::vector<float> trigger_scores;
    /// peakGroup isolation window ranges
    std::vector<double> trigger_left_isolation_mzs_;
    std::vector<double> trigger_right_isolation_mzs_;
    std::vector<int> trigger_ids_;

    /// SpectralDeconvolution class for deconvolution
    SpectralDeconvolution fd_;

    /// total QScore threshold
    //double tqscore_threshold = .99;
    double tqscore_threshold;

    /// q score threshold - determined from C# side
    double qscore_threshold_;
    /// retention time window - determined from C# side
    double rt_window_;
    /// how many masses will be selected per ms level? - determined from C# side
    IntList mass_count_;

    int targeting_mode_ = 0; /// 0 no targeting 1 inclusive 2 exclusive 3 deep

    /// maps for global inclusion targeting
    std::map<double, std::vector<double>> target_mass_rt_map_;
    std::map<double, std::vector<double>> target_mass_qscore_map_;
    std::map<double, std::vector<int>> target_mass_charge_map_;
    std::vector<double> target_masses_; /// current target masses. updated per spectrum

    // For the possibility of removal each window is given an id, starting at zero (needed for FAIMS support)
    int window_id_ = 0;

    /// maps for global exclusion
    std::map<double, std::vector<double>> exclusion_rt_masses_map_; /// if rt == 0, its mapped masses are always excluded.
    std::vector<double> excluded_masses_; /// current target masses. updated per spectrum

    /// TSV-based inclusion targets
    std::vector<InclusionTarget> inclusion_targets_;  ///< All targets loaded from TSV file
    std::vector<const InclusionTarget*> active_targets_;  ///< Current active targets (filtered by RT)
    std::map<int, int> target_priority_map_;  ///< nominal_mass → priority for tie-breaking
    double tie_threshold_ = 0.1;  ///< qscore difference threshold for priority tie-breaking
    bool strict_inclusion_ = true;  ///< If true, only acquire targets in inclusion mode; if false, non-targets can fill remaining slots

    /// Tag-based protein family targeting
    std::vector<FASTAFile::FASTAEntry> target_protein_database_;  ///< Target protein family entries
    std::vector<TargetPTM> target_ptms_;                          ///< PTM modifications for mass expansion
    bool tag_based_targeting_enabled_ = false;                    ///< Flag indicating tag-based mode active
    int min_tag_length_for_targeting_ = 3;                        ///< Minimum tag length for matching
    int max_tag_length_for_targeting_ = 8;                        ///< Maximum tag length for matching
    double tag_matching_tolerance_ppm_ = 10.0;                    ///< Mass tolerance for tag matching
    double max_flanking_mass_diff_ = 50000.0;                       ///< Max flanking mass diff for tag matches
    std::set<double> expanded_target_masses_;                     ///< Track already-expanded masses (avoid duplicates)
    int max_total_ptm_count_ = 3;                                 ///< Maximum total PTM modifications per proteoform

    /// precursor SNR threshold
    double snr_threshold_ = 1;
    bool use_idscore_ = false;
    bool consider_all_Charge_states_ = false;
    bool ms3_all_charges_ = false;  ///< If true, output all charge states for MS3 fragment ions
    int hcd_energy_ = -1;

    /// mass tolerance
    DoubleList tol_;

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

    /// Parse JSON configuration string (Phase 1: replaces legacy space-delimited format)
    void parseJSONConfig_(const std::string& json_str);

    // --- Phase 1 JSON-only member variables (stored for future phases) ---

    // ms_settings (stored for Phase 3+ scan command construction)
    std::string ms1_analyzer_;
    double ms1_first_mass_ = 0, ms1_last_mass_ = 0, ms1_max_it_ = 0;
    int ms1_resolution_ = 0, ms1_agc_target_ = 0;

    struct MS2ConfigJson
    {
      std::string analyzer, activation;
      int collision_energy = 0, resolution = 0;
    };
    std::vector<MS2ConfigJson> ms2_configs_;

    // scheduling (stored for Phase 3)
    bool cycle_time_enabled_ = false, timeout_enabled_ = false;
    double cycle_time_ms_ = 60000.0, timeout_ms_ = 30000.0;

    // exploration (stored for Phase 7)
    bool exploration_enabled_ = false;
    int exploration_max_depth_ = 1, exploration_max_variants_ = 5;

    // quantification (stored — currently passed per-call via IsDifferentiallyAbundant)
    bool quant_enabled_ = false;
    double reporter_mz_tol_ = 0.002, fold_change_threshold_ = 1.4;

    // FAIMS
    std::vector<double> faims_cv_values_;
    int max_cv_skip_ = 0;

  };
}
