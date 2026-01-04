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
     * @brief Deconvolute a spectrum and find sequence tags with database matches
     *
     * This method performs spectrum deconvolution using the FLASHDeconv algorithm,
     * generates sequence tags using FLASHTagger, and matches them against a protein database.
     *
     * @param mzs m/z values of the input spectrum
     * @param ints intensities of the input spectrum
     * @param length number of peaks in the spectrum
     * @param rt retention time in seconds
     * @param ms_level MS level of the spectrum
     * @param fasta_entries protein database entries to match against
     * @param tagger_param parameters for the FLASHTagger algorithm (optional, uses defaults if empty)
     * @param tags output vector of detected sequence tags
     * @param matches output vector of tag matches to database entries
     * @param ppm_tolerance mass tolerance in ppm for tag matching (default 10.0)
     * @param max_flanking_mass_diff maximum allowed flanking mass difference for a match (default 500.0 Da)
     * @return number of tags found
     */
    int getSequenceTagsAndMatches(const double* mzs,
                                  const double* ints,
                                  int length,
                                  double rt,
                                  int ms_level,
                                  const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                  const Param& tagger_param,
                                  std::vector<FLASHHelperClasses::Tag>& tags,
                                  std::vector<TagMatch>& matches,
                                  double ppm_tolerance = 10.0,
                                  double max_flanking_mass_diff = 500.0);

    /**
     * @brief Identify proteoform from MS2 spectrum against a single protein sequence
     *
     * Implements the core FLASHTnT identification workflow:
     * 1. Deconvolves the MS2 spectrum to get monoisotopic masses
     * 2. Calculates theoretical fragment masses from protein sequence
     * 3. Matches observed masses against theoretical fragments
     * 4. Identifies PTM positions based on mass differences
     *
     * @param mzs m/z values of the input MS2 spectrum
     * @param ints intensities of the input MS2 spectrum
     * @param length number of peaks in the spectrum
     * @param rt retention time in seconds
     * @param protein_sequence the protein sequence to match against
     * @param ppm_tolerance mass tolerance in ppm (default 10.0)
     * @param ion_types ion types to consider (default {"b", "y"})
     * @param ptm_mass_threshold minimum mass shift to consider as PTM (default 5.0 Da)
     * @param matched_fragment_indices output: indices of matched fragment ions (1-based positions)
     * @param ptm_start_positions output: start positions of PTM localization ranges (1-based)
     * @param ptm_end_positions output: end positions of PTM localization ranges (1-based)
     * @param ptm_masses output: mass shifts at each PTM position
     * @return number of matched fragment ions
     */
    int identifyProteoform(const double* mzs,
                           const double* ints,
                           int length,
                           double rt,
                           const String& protein_sequence,
                           double ppm_tolerance,
                           const std::vector<String>& ion_types,
                           double ptm_mass_threshold,
                           std::vector<int>& matched_fragment_indices,
                           std::vector<int>& ptm_start_positions,
                           std::vector<int>& ptm_end_positions,
                           std::vector<double>& ptm_masses);

    /**
     * @brief Extended proteoform identification with detailed output
     *
     * Extended version providing more detailed match information including
     * peak indices, theoretical masses, ion types, and PTM localization ranges.
     *
     * @param mzs m/z values of the input MS2 spectrum
     * @param ints intensities of the input MS2 spectrum
     * @param length number of peaks in the spectrum
     * @param rt retention time in seconds
     * @param protein_sequence the protein sequence to match against
     * @param ppm_tolerance mass tolerance in ppm
     * @param ion_types ion types to consider
     * @param max_ptm_count maximum number of PTMs to consider
     * @param max_ptm_mass maximum mass shift for a single PTM
     * @param matched_peak_indices output: indices of matched peaks in deconvolved spectrum
     * @param matched_theoretical_masses output: theoretical masses that were matched
     * @param matched_ion_types output: ion types (true=N-terminal/prefix, false=C-terminal/suffix)
     * @param ptm_start_positions output: start of region for each PTM
     * @param ptm_end_positions output: end of region for each PTM
     * @param ptm_masses output: mass shift for each PTM
     * @param coverage output: sequence coverage (0.0-1.0)
     * @param total_score output: total identification score
     * @return number of matched ions
     */
    int identifyProteoformExtended(const double* mzs,
                                   const double* ints,
                                   int length,
                                   double rt,
                                   const String& protein_sequence,
                                   double ppm_tolerance,
                                   const std::vector<String>& ion_types,
                                   int max_ptm_count,
                                   double max_ptm_mass,
                                   std::vector<int>& matched_peak_indices,
                                   std::vector<double>& matched_theoretical_masses,
                                   std::vector<bool>& matched_ion_types,
                                   std::vector<int>& ptm_start_positions,
                                   std::vector<int>& ptm_end_positions,
                                   std::vector<double>& ptm_masses,
                                   double& coverage,
                                   double& total_score);

    // ============ Python-friendly overloads using vectors ============

    /**
     * @brief Python-friendly overload of getSequenceTagsAndMatches using vectors
     */
    int getSequenceTagsAndMatchesPy(const std::vector<double>& mzs,
                                    const std::vector<double>& ints,
                                    double rt,
                                    int ms_level,
                                    const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                    const Param& tagger_param,
                                    std::vector<FLASHHelperClasses::Tag>& tags,
                                    std::vector<TagMatch>& matches,
                                    double ppm_tolerance = 10.0,
                                    double max_flanking_mass_diff = 500.0);

    /**
     * @brief Python-friendly overload of identifyProteoform using vectors
     */
    int identifyProteoformPy(const std::vector<double>& mzs,
                             const std::vector<double>& ints,
                             double rt,
                             const String& protein_sequence,
                             double ppm_tolerance,
                             const std::vector<String>& ion_types,
                             double ptm_mass_threshold,
                             std::vector<int>& matched_fragment_indices,
                             std::vector<int>& ptm_start_positions,
                             std::vector<int>& ptm_end_positions,
                             std::vector<double>& ptm_masses);

    /**
     * @brief Python-friendly overload of identifyProteoformExtended using vectors
     */
    int identifyProteoformExtendedPy(const std::vector<double>& mzs,
                                     const std::vector<double>& ints,
                                     double rt,
                                     const String& protein_sequence,
                                     double ppm_tolerance,
                                     const std::vector<String>& ion_types,
                                     int max_ptm_count,
                                     double max_ptm_mass,
                                     std::vector<int>& matched_peak_indices,
                                     std::vector<double>& matched_theoretical_masses,
                                     std::vector<bool>& matched_ion_types,
                                     std::vector<int>& ptm_start_positions,
                                     std::vector<int>& ptm_end_positions,
                                     std::vector<double>& ptm_masses,
                                     double& coverage,
                                     double& total_score);

    /**
     * @brief Identify proteoform from pre-deconvolved masses (core identification function)
     *
     * This function performs fragment ion matching without deconvolution,
     * allowing direct testing with known masses. It takes deconvolved monoisotopic
     * masses directly and matches them against theoretical fragment masses.
     *
     * @param observed_masses deconvolved monoisotopic masses
     * @param mass_scores quality scores for each mass (use 1.0 if unknown)
     * @param protein_sequence the protein sequence to match against
     * @param ppm_tolerance mass tolerance in ppm
     * @param ion_types ion types to consider (e.g., {"b", "y"})
     * @param ptm_mass_threshold minimum mass shift to consider as PTM
     * @param matched_fragment_indices output: matched fragment indices (1-based positions)
     * @param matched_ion_types output: true=prefix/b-ion, false=suffix/y-ion
     * @param matched_observed_masses output: observed mass for each match
     * @param matched_theoretical_masses output: theoretical mass for each match
     * @param matched_ppm_errors output: ppm error for each match
     * @param ptm_start_positions output: start positions of PTM localization ranges (1-based)
     * @param ptm_end_positions output: end positions of PTM localization ranges (1-based)
     * @param ptm_masses output: mass shifts at PTM positions
     * @return number of matched fragment ions
     */
    int identifyProteoformFromMasses(const std::vector<double>& observed_masses,
                                     const std::vector<double>& mass_scores,
                                     const String& protein_sequence,
                                     double ppm_tolerance,
                                     const std::vector<String>& ion_types,
                                     double ptm_mass_threshold,
                                     std::vector<int>& matched_fragment_indices,
                                     std::vector<bool>& matched_ion_types,
                                     std::vector<double>& matched_observed_masses,
                                     std::vector<double>& matched_theoretical_masses,
                                     std::vector<double>& matched_ppm_errors,
                                     std::vector<int>& ptm_start_positions,
                                     std::vector<int>& ptm_end_positions,
                                     std::vector<double>& ptm_masses);

    /**
     * @brief Python-friendly overload of identifyProteoformFromMasses
     *
     * Same as identifyProteoformFromMasses but with all vector parameters for Python binding.
     */
    int identifyProteoformFromMassesPy(std::vector<double>& observed_masses,
                                       std::vector<double>& mass_scores,
                                       const String& protein_sequence,
                                       double ppm_tolerance,
                                       std::vector<String>& ion_types,
                                       double ptm_mass_threshold,
                                       std::vector<int>& matched_fragment_indices,
                                       std::vector<bool>& matched_ion_types,
                                       std::vector<double>& matched_observed_masses,
                                       std::vector<double>& matched_theoretical_masses,
                                       std::vector<double>& matched_ppm_errors,
                                       std::vector<int>& ptm_start_positions,
                                       std::vector<int>& ptm_end_positions,
                                       std::vector<double>& ptm_masses);

    /**
     * @brief Calculate theoretical fragment masses for a protein sequence (static, for Python)
     *
     * @param sequence the protein sequence
     * @param ion_types ion types to consider (e.g., {"b", "y"})
     * @param prefix_masses output: cumulative masses from N-terminus
     * @param suffix_masses output: cumulative masses from C-terminus
     */
    static void calculateTheoreticalFragmentMassesPy(const String& sequence,
                                                     const std::vector<String>& ion_types,
                                                     std::vector<double>& prefix_masses,
                                                     std::vector<double>& suffix_masses);

    /**
     * @brief Process MS2 spectrum for protein family detection and inclusion list expansion
     *
     * This DLL bridge function performs real-time tag-based targeting:
     * 1. Deconvolves the MS2 spectrum
     * 2. Extracts sequence tags (minimum length 3)
     * 3. Matches tags against the target protein database
     * 4. If match found: Expands target masses using PTM combinations from precursor mass
     * 5. Adds expanded masses to dynamic inclusion list
     *
     * @param mzs m/z values of the input spectrum
     * @param ints intensities of the input spectrum
     * @param length number of peaks in the spectrum
     * @param rt retention time in seconds
     * @param ms_level MS level of the spectrum (must be 2)
     * @param name spectrum name
     * @param cv CV value for FAIMS (can be nullptr)
     * @param precursor_mass monoisotopic mass of the precursor (from iAPI)
     * @return true if target protein detected and targets expanded, false otherwise
     */
    bool processMS2ForTagBasedTargeting(const double* mzs,
                                        const double* ints,
                                        int length,
                                        double rt,
                                        int ms_level,
                                        const char* name,
                                        const char* cv,
                                        double precursor_mass);

  private:
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
    double tie_threshold_ = 0.01;  ///< qscore difference threshold for priority tie-breaking
    bool strict_inclusion_ = true;  ///< If true, only acquire targets in inclusion mode; if false, non-targets can fill remaining slots

    /// Tag-based protein family targeting
    std::vector<FASTAFile::FASTAEntry> target_protein_database_;  ///< Target protein family entries
    std::vector<TargetPTM> target_ptms_;                          ///< PTM modifications for mass expansion
    bool tag_based_targeting_enabled_ = false;                    ///< Flag indicating tag-based mode active
    int min_tag_length_for_targeting_ = 3;                        ///< Minimum tag length for matching
    double tag_matching_tolerance_ppm_ = 10.0;                    ///< Mass tolerance for tag matching
    double max_flanking_mass_diff_ = 500.0;                       ///< Max flanking mass diff for tag matches
    std::set<double> expanded_target_masses_;                     ///< Track already-expanded masses (avoid duplicates)
    int max_total_ptm_count_ = 3;                                 ///< Maximum total PTM modifications per proteoform

    /// precursor SNR threshold
    double snr_threshold_ = 1;
    bool use_idscore_ = false;
    bool consider_all_Charge_states_ = false;
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

  };
}
