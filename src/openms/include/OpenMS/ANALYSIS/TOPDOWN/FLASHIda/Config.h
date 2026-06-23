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

#include <OpenMS/config.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{

  /// Selection metric: how targets are ranked for MSn+1
  enum class SelectionMetric
  {
    None = 0,              ///< No selection at this level -- don't select targets for MSn+1
    Intensity,             ///< Rank by raw intensity
    QScore,                ///< Rank by deconvolution quality score
    TerminalFragments,     ///< Innermost b/y ions, interleaved
    AmbiguityResolution    ///< PTM-site bracketing ions
  };

  /// Exploration metric: what to optimize during CE sweep (MS2+ only)
  enum class ExplorationMetric
  {
    None = 0,            ///< No exploration at this level (default)
    MassCount,           ///< Optimize for most deconvolved masses
    RemainingPrecursor,  ///< Optimize for least remaining precursor intensity
    FragmentCount        ///< Optimize for most fragment ions
  };

  /// Configuration for a single scan type (MS1 survey or one MS2 config slot)
  struct OPENMS_DLLAPI ScanConfig
  {
    std::string analyzer = "Orbitrap";
    std::string activation;       ///< empty for MS1
    int collision_energy = 0;
    int resolution = 0;
    int agc_target = 0;
    double first_mass = 0;
    double last_mass = 0;
    double max_it = 0;
    int microscans = 0;
    double rf_lens = 0;
    double source_cid = 0;
    double source_cid_scaling = 0;
    std::string data_type;
    std::string scan_rate;
    double reaction_time = 0;          ///< Ion/ion reaction time (ms), 0 = not used (ETD)
    double reagent_max_it = 0;         ///< Reagent max injection time (ms), 0 = not used (ETD)
    int reagent_agc_target = 0;        ///< Reagent AGC target, 0 = not used (ETD)

    /// Apply string-keyed overrides to matching fields (exploration overrides)
    void applyOverrides(const std::unordered_map<std::string, std::string>& overrides);
  };

  /// Per-MSn-level configuration: selection + optional exploration
  struct OPENMS_DLLAPI MSLevelConfig
  {
    std::vector<ScanConfig> scans;  ///< [0]=primary, [1]=conditional follow-up
    SelectionMetric selection = SelectionMetric::Intensity;
    int max_targets = 10;
    int min_charge = 0;  ///< Minimum charge for target selection (0 = no filter)

    ExplorationMetric exploration = ExplorationMetric::None;
    double ce_min = 20.0;
    double ce_max = 40.0;
    double ce_step = 5.0;
    double rt_min = 0;                 ///< Reaction time sweep min (ms), 0 = no RT sweep
    double rt_max = 0;                 ///< Reaction time sweep max (ms)
    double rt_step = 1.0;              ///< Reaction time sweep step (ms)
    std::vector<std::string> activations;  ///< Activation types to sweep (e.g. {"HCD","ETD"})
    std::unordered_map<std::string, std::string> overrides;
    double tolerance_ppm = 10.0;
    double exploration_tolerance_ppm = 10.0;  ///< Resolved exploration tolerance (from overrides or base tol)
    double remaining_precursor_target = 0.1;  ///< Target remaining precursor ratio (0.1 = 10%)
  };

  /// Deconvolution engine parameters
  struct OPENMS_DLLAPI DeconvolutionConfig
  {
    int min_charge = 4;
    int max_charge = 50;
    double min_mass = 500.0;
    double max_mass = 50000.0;
  };

  /// FAIMS configuration values (runtime state lives on FLASHIda)
  struct OPENMS_DLLAPI FAIMSConfig
  {
    std::vector<double> cv_values;
    int max_cv_skip = 0;
    int precursor_threshold = 15;
    bool enabled = false;
  };

  /// Scheduling and timing configuration
  struct OPENMS_DLLAPI SchedulingConfig
  {
    bool cycle_time_enabled = false;
    bool timeout_enabled = false;
    double cycle_time_ms = 60000.0;
    double timeout_ms = 30000.0;
    uint64_t agc_interval_ms = 30000;
  };

  /// Targeting and precursor selection configuration
  struct OPENMS_DLLAPI TargetingConfig
  {
    int mode = 0;                      ///< 0=none, 1=inclusive, 2=exclusive, 3=deep
    bool strict_inclusion = false;
    double tie_threshold = 0.1;
    double rt_window = 180.0;
    bool consider_all_charges = false;
    bool charge_based_exclusion = false;  ///< Treat each (mass, charge) as an independent exclusion target (developer flag).
    int hcd_energy = -1;
    double qscore_threshold = 0.0;
    double tqscore_threshold = 0.9;
    double snr_threshold = 1.0;
    std::string protein_sequence;
    bool conditional_ms2_enabled = false; ///< Explicit conditional_ms2 flag from JSON
    std::vector<std::string> target_log_files;
    std::string fasta_file;
    std::string inclusion_list_file;
    std::string ptm_list_file;
    bool tag_based_enabled = false;
    int min_tag_length = 3;
    int max_tag_length = 8;
    double tag_matching_tolerance_ppm = 10.0;
    double max_flanking_mass_diff = 50000.0;
    int max_total_ptm_count = 3;
    ScanConfig tagging_follow_up_scan;  ///< Follow-up scan config for conditional MS2
  };

  /// Isobaric quantification configuration
  struct OPENMS_DLLAPI QuantConfig
  {
    bool enabled = false;
    double reporter_mz_tol = 0.002;
    double fold_change_threshold = 1.4;
    ScanConfig follow_up_scan;  ///< Follow-up scan config for quant follow-up MS2
  };

  /// Runtime file paths (set by C# or user override in JSON)
  struct OPENMS_DLLAPI RuntimeConfig
  {
    std::string ida_log_path;
    std::string scan_commands_path;
    std::string scan_results_path;
    std::string identification_path;
  };

  /**
   * @brief Owns all configuration parsing and storage for FLASHIda.
   *
   * Constructed from a JSON string (the same format produced by C# Parameter.ToJSON()).
   * Throws std::invalid_argument if the input is not JSON.
   * File-loading code (FASTA, log file parsing, TSV parsing) stays in FLASHIda.
   */
  class OPENMS_DLLAPI Config
  {
  public:
    /// Construct from JSON configuration string. Throws std::invalid_argument for non-JSON input.
    explicit Config(const std::string& json_str);

    /// Access per-level config; returns default if level is not configured
    const MSLevelConfig& level(int msn_level) const;

    /// Returns true if exploration is enabled at the given MSn level
    bool hasExploration(int msn_level) const;

    /// Returns true if any level has exploration enabled
    bool explorationEnabled() const { return exploration_enabled_; }

    /// Tolerance-ppm list across configured MS levels (ascending level order), for Deconvolution construction
    DoubleList toleranceList() const;

    /// Exploration-tolerance-ppm list across configured MS levels (ascending level order)
    DoubleList explorationToleranceList() const;

    /// Accessors for config sub-sections
    const DeconvolutionConfig& deconvolution() const { return deconv_; }
    const FAIMSConfig& faims() const { return faims_; }
    const SchedulingConfig& scheduling() const { return scheduling_; }
    const TargetingConfig& targeting() const { return targeting_; }
    const QuantConfig& quantification() const { return quant_; }
    const RuntimeConfig& runtime() const { return runtime_; }

    /// Validate configuration consistency; throws std::invalid_argument on conflict
    void validate() const;

    /// Access the full levels map (for iteration)
    const std::map<int, MSLevelConfig>& levels() const { return levels_; }

  private:
    std::map<int, MSLevelConfig> levels_;
    DeconvolutionConfig deconv_;
    FAIMSConfig faims_;
    SchedulingConfig scheduling_;
    TargetingConfig targeting_;
    QuantConfig quant_;
    RuntimeConfig runtime_;
    bool exploration_enabled_ = false;

    static const MSLevelConfig default_level_;
  };

} // namespace OpenMS
