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

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <map>
#include <set>
#include <unordered_map>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Precursor selection, targeting, and mass exclusion for FLASHIda.
   *
   * Responsible for:
   * - Filtering and ranking deconvolved MS1 peak groups
   * - Mass exclusion (mz-based and mass-based)
   * - Inclusion/exclusion targeting (TSV-based, log-file-based, tag-based)
   * - Dynamic target expansion via PTM combinations
   * - Loading targeting data (FASTA, TSV, log files) from config paths
   */
  class OPENMS_DLLAPI PrecursorSelection
  {
  public:
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

    /// Constructor: loads targeting data (FASTA, TSV, log files) from config
    PrecursorSelection(const Config& config, Deconvolution& deconv);

    /**
     * @brief Filter and rank deconvolved MS1 spectrum. Returns number of selected precursors.
     *
     * Populates trigger_charges, trigger_hcds, trigger_scores internally.
     * Requires deconv_.deconvolveMS1() to have been called first.
     *
     * @param mzs m/z values of the input spectrum
     * @param ints intensities of the input spectrum
     * @param length number of peaks
     * @param rt retention time in seconds
     * @param ms_level MS level
     * @param faims_cv FAIMS CV value (0.0 if FAIMS not enabled)
     * @return number of selected peak groups
     */
    int filterAndRank(const double* mzs, const double* ints, int length,
                      double rt, int ms_level, double faims_cv);

    /// Access selected peak groups after filterAndRank (const)
    const DeconvolvedSpectrum& selectedPeakGroups() const { return selected_peak_groups_; }

    /// Access selected peak groups after filterAndRank (mutable)
    DeconvolvedSpectrum& selectedPeakGroups() { return selected_peak_groups_; }

    /// Access full deconvolved MS1 spectrum (all peak groups, not just selected)
    const DeconvolvedSpectrum& deconvolvedMS1() const { return deconv_.deconvolvedMS1(); }

    /// Access trigger charges populated by filterAndRank
    const std::vector<int>& triggerCharges() const { return trigger_charges_; }

    /// Access trigger HCD energies populated by filterAndRank
    const std::vector<int>& triggerHcds() const { return trigger_hcds_; }

    /// Access trigger scores populated by filterAndRank
    const std::vector<float>& triggerScores() const { return trigger_scores_; }

    /// Access trigger isolation window left m/z values
    const std::vector<double>& triggerLeftIsolationMzs() const { return trigger_left_isolation_mzs_; }

    /// Access trigger isolation window right m/z values
    const std::vector<double>& triggerRightIsolationMzs() const { return trigger_right_isolation_mzs_; }

    /// Access trigger IDs
    const std::vector<int>& triggerIds() const { return trigger_ids_; }

    /**
     * @brief Process stored MS2 deconvolution for protein family detection and inclusion list expansion
     *
     * Performs real-time tag-based targeting:
     * 1. Uses stored MS2 deconvolution results (requires deconvolveMSn() first)
     * 2. Extracts sequence tags (minimum length 3)
     * 3. Matches tags against the target protein database
     * 4. If match found: Expands target masses using PTM combinations from precursor mass
     * 5. Adds expanded masses to dynamic inclusion list
     *
     * @param precursor_mass monoisotopic mass of the precursor (from iAPI)
     * @return true if target protein detected and targets expanded, false otherwise
     */
    bool processMS2ForTagBasedTargeting(double precursor_mass, const std::string& activation_type);

    /**
     * @brief Remove a given precursor from the exclusion list by id (needed for FAIMS)
     * @param id id of precursor
     */
    void removeFromExclusionList(int id);

    /**
     * @brief Fill legacy isolation window arrays from the last filterAndRank call
     */
    void getIsolationWindows(double* window_start, double* window_end,
                             double* qscores, int* charges,
                             int* min_charges, int* max_charges,
                             double* mono_masses, double* charge_cos,
                             double* charge_snrs, double* iso_cos,
                             double* snrs, double* charge_scores,
                             double* ppm_errors, double* precursor_intensities,
                             double* peakgroup_intensities, int* hcds, int* ids);

    /// Check if tag-based targeting database is loaded (non-empty)
    bool hasTargetProteinDatabase() const { return !target_protein_database_.empty(); }

  private:
    const Config& config_;
    Deconvolution& deconv_;

    /// Selected peak groups out of MS1 deconvolution (filtered subset for triggering)
    DeconvolvedSpectrum selected_peak_groups_{0};

    /// Trigger arrays populated by filterAndRank
    std::vector<int> trigger_charges_;
    std::vector<int> trigger_hcds_;
    std::vector<float> trigger_scores_;
    std::vector<double> trigger_left_isolation_mzs_;
    std::vector<double> trigger_right_isolation_mzs_;
    std::vector<int> trigger_ids_;

    /// Mass exclusion maps
    std::unordered_map<int, double> tqscore_exceeding_mz_rt_map_;
    std::unordered_map<int, double> tqscore_exceeding_mass_rt_map_;
    std::unordered_map<int, double> all_mass_rt_map_;
    std::unordered_map<int, double> mass_qscore_map_;

    /// Per-(nominal_mass, charge) cross-scan exclusion set (charge_based_exclusion flag).
    std::set<std::pair<int, int>> tqscore_exceeding_mass_charge_set_;

    /// Per-(nominal_mass, charge) qscore accumulator, parallel to mass_qscore_map_
    /// but one level deeper. Only touched when charge_based_exclusion is on.
    std::map<std::pair<int, int>, double> mass_charge_qscore_map_;

    /// Maps for selectively disabling mass exclusion (needed for FAIMS support)
    std::unordered_map<int, int> id_mass_map_;
    std::unordered_map<int, int> id_mz_map_;
    std::unordered_map<int, double> id_qscore_map_;
    std::unordered_map<int, int> id_charge_map_;

    /// Window ID counter for FAIMS exclusion list removal
    int window_id_ = 0;

    /// Maps for global inclusion targeting
    std::map<double, std::vector<double>> target_mass_rt_map_;
    std::map<double, std::vector<double>> target_mass_qscore_map_;
    std::map<double, std::vector<int>> target_mass_charge_map_;
    std::vector<double> target_masses_;

    /// Maps for global exclusion
    std::map<double, std::vector<double>> exclusion_rt_masses_map_;
    std::vector<double> excluded_masses_;

    /// TSV-based inclusion targets
    std::vector<InclusionTarget> inclusion_targets_;
    std::vector<const InclusionTarget*> active_targets_;
    std::map<int, int> target_priority_map_;

    /// Tag-based protein family targeting
    std::vector<FASTAFile::FASTAEntry> target_protein_database_;
    std::vector<TargetPTM> target_ptms_;
    std::set<double> expanded_target_masses_;

    /// Filter peak groups using mass exclusion logic
    void filterPeakGroupsUsingMassExclusion_(int ms_level, double rt);

    /// Parse TSV inclusion list file
    void parseInclusionListTSV_(const String& filename);

    /// Parse TSV file containing PTM modifications for target expansion
    void parseTargetPTMsTSV_(const String& filename);

    /// Generate all PTM mass combinations for a base protein mass
    std::vector<double> generatePTMCombinations_(double base_mass,
                                                 const std::vector<TargetPTM>& ptms) const;

    /// Add expanded target masses to the dynamic inclusion list
    void addDynamicTargets_(const std::vector<double>& masses, double rt, int priority);
  };

} // namespace OpenMS
