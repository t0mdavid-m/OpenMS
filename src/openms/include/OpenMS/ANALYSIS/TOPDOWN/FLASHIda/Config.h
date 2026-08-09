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
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{

  /// Objective driving characterization scan planning for a proteoform
  enum class CharacterizationObjective
  {
    Ambiguity, ///< Resolve PTM site ambiguity
    Coverage   ///< Extend sequence coverage
  };

  /// Selection metric: how targets are ranked for MSn+1
  enum class SelectionMetric
  {
    None = 0,              ///< No selection at this level -- don't select targets for MSn+1
    Intensity,             ///< Rank by raw intensity
    QScore,                ///< Rank by deconvolution quality score
    TerminalFragments,     ///< Innermost b/y ions, interleaved
    AmbiguityResolution    ///< PTM-site bracketing ions
  };

  /// How many charge states of one species a scan acquires, and how (ADR-0016).
  ///
  /// The unit of acquisition is a Precursor's *acquisition charge set*, not a single charge.
  /// Membership is SNR-gated: a charge joins only if its own envelope rises above noise, because a
  /// charge contributing no signal still consumes part of the scan's ion budget.
  ///
  /// This is an acquisition-geometry question, deliberately separate from
  /// `charge_based_exclusion`, which is an exclusion-KEYING question: that flag decides whether a
  /// mass already fragmented at one charge stays eligible at another on a LATER survey (ADR-0018).
  enum class ChargeAcquisitionMode
  {
    Single = 0,   ///< One charge per detection -- the representative / best-qscore charge. Default.
    Separate,     ///< One scan PER charge state; each is its own Precursor with its own model.
    Multiplexed   ///< ONE scan co-isolating the whole set as notches; one Precursor, one model.
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
    /// The DISPATCH ROSTER for this level, in dispatch order -- NOT every scan config defined for
    /// it. Config materialises this from ms_settings.msN plus the names listed in
    /// precursor_selection.additional_scans, so a block that exists only to back a follow-up
    /// reference is absent here and never fires unconditionally.
    /// (The old comment claimed "[1]=conditional follow-up". That was never true: FLASHIda.cpp:211
    /// range-fors the whole vector, so entry [1] was a second UNCONDITIONAL MS2.)
    std::vector<ScanConfig> scans;

    /// Defaults to None, not Intensity. This is load-bearing: ms_settings.msN materialises
    /// levels_[N] BEFORE any selection is parsed (Config.cpp), so an in-class default of Intensity
    /// meant that merely DEFINING a scan config switched that level on. Failing closed here is what
    /// makes characterization.mode the single MS3 switch. Every level is assigned explicitly by
    /// applyCharacterizationMode_() regardless, so this default is a safety net, not the mechanism.
    SelectionMetric selection = SelectionMetric::None;
    int max_targets = 10;
    int min_charge = 0;  ///< Minimum charge for target selection (0 = no filter)

    ExplorationMetric exploration = ExplorationMetric::None;
    double ce_min = 20.0;
    double ce_max = 40.0;
    double ce_step = 5.0;
    /// Ion-ion REACTION time (ms), not retention time. Wire keys are reaction_time_{min,max,step};
    /// "rt" elsewhere in this codebase means retention time (TargetingConfig::rt_window, seconds).
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
    /// 0=none, 1=inclusion, 2=in_depth, 3=exclusion_masses -- from precursor_selection.targeting.
    /// This comment said "2=exclusive, 3=deep" for a long time, and it was BACKWARDS. Mode 2 loads
    /// `Mass` lines into target_mass_rt_map_ and de-prioritizes them via the tqscore product map;
    /// mode 3 loads `AllMass` lines and hard-skips matches. Two committed configs were named for the
    /// mode they did not set because of this. Trust the parse in Config.cpp and the branches in
    /// PrecursorSelection.cpp, not prose.
    int mode = 0;
    bool strict_inclusion = false;
    double tie_threshold = 0.1;
    double rt_window = 180.0;
    bool consider_all_charges = false;
    bool charge_based_exclusion = false;  ///< Treat each (mass, charge) as an independent exclusion target (developer flag).
    /// How many charge states of a selected precursor one MS2 acquires -- from
    /// precursor_selection.precursor_charges. Orthogonal to charge_based_exclusion above: that flag
    /// keys EXCLUSION per (mass, charge) and makes a later survey fall back to the next unexcluded
    /// charge; this decides the isolation GEOMETRY of a single scan (ADR-0016/0018).
    ChargeAcquisitionMode precursor_charges = ChargeAcquisitionMode::Single;
    // NOTE: there is no hcd_energy here any more. The precursor_selection.HCDEnergy key was deleted
    // (ADR-0014) because its only export, PrecursorSelection::getIsolationWindows(), had zero callers
    // repo-wide. Do not confuse it with ScanCommand::hcd_energy, which is alive and unrelated: that
    // one mirrors stages[0].collision_energy into the log schema and is part of the 2048-byte ABI.
    double qscore_threshold = 0.0;
    double tqscore_threshold = 0.9;
    double snr_threshold = 1.0;
    bool conditional_ms2_enabled = false; ///< Explicit conditional_ms2 flag from JSON
    std::vector<std::string> target_log_files;
    std::string fasta_file;
    std::string inclusion_list_file;
    std::string ptm_list_file;
    int min_tag_length = 3;
    int max_tag_length = 8;
    double tag_matching_tolerance_ppm = 10.0;
    double max_flanking_mass_diff = 50000.0;
    int max_total_ptm_count = 3;
    // --- FLASHTnT (FLASHTagger/FLASHExtender) tuning, from the `flashtnt` config block ---
    bool allow_gap = false;                 ///< FLASHTagger allow_gap (passed as "true"/"false")
    int max_aa_in_gap = 2;                   ///< FLASHTagger max_aa_in_gap
    std::vector<std::string> fixed_mod;      ///< Fixed mods, passed VERBATIM to tagger+extender; empty = empty Param list, NOT the declared {""}
    int max_blind_mod_count = 2;             ///< FLASHExtender max_blind_mod_count
    double max_mod_mass = 700.0;             ///< FLASHExtender max_mod_mass. 700 preserves current MS2 behavior; NOT the extender's own 500 default.
    ScanConfig tagging_follow_up_scan;  ///< Follow-up scan config for conditional MS2, RESOLVED from a name
    /// Whether tagging.follow_up_scan named anything. Replaces the old
    /// `tagging_follow_up_scan.activation.empty()` sentinel: with the block now defined elsewhere
    /// and merely referenced, "is one configured" is a property of the REFERENCE, not of the
    /// resolved block's contents.
    bool has_tagging_follow_up = false;
    std::string tagging_follow_up_name;  ///< the referenced name, for diagnostics only
  };

  /// Characterization configuration: objective + protein sequence
  /// The single MS3 switch. `Off` means no MS3 is emitted at all; the two on-values ARE the
  /// objectives, so there is no separate enable flag and no way to express "on but with no
  /// objective". Supersedes ADR-0004, which decided there would be no enable flag at all.
  enum class CharacterizationMode
  {
    Off,
    Ambiguity,
    Coverage
  };

  struct OPENMS_DLLAPI CharacterizationConfig
  {
    CharacterizationMode mode = CharacterizationMode::Off;

    /// Kept as a separate field so every existing objective read site (ProteoformTracker.cpp:353
    /// and friends) is untouched. It is derived from `mode` at parse time and is meaningless when
    /// mode == Off, because nothing downstream runs in that case.
    CharacterizationObjective objective = CharacterizationObjective::Ambiguity;

    std::string protein_sequence;

    /// The MS3 budget, and the fragment-charge floor for MS3 targeting. Moved here from
    /// selection_strategy.ms2.{max_targets,min_charge}, which named them one level off their
    /// effect. Projected into levels_[2] so the existing reads (ProteoformTracker.cpp:354,
    /// Exploration.cpp:733/:800) keep working unchanged.
    int max_targets = 3;
    int min_fragment_charge = 0;

    /// How many charge states of a target FRAGMENT one MS3 acquires -- from
    /// characterization.fragment_charges. Replaces the bool `ms3_all_charges`, whose two states are
    /// now Single (the fragment's best-MS2 charge) and Separate (one MS3 per observed charge);
    /// Multiplexed co-isolates them into one MS3, which is the free direction because synchronous
    /// precursor selection is one simultaneous waveform rather than N sequential fills.
    ChargeAcquisitionMode fragment_charges = ChargeAcquisitionMode::Single;
  };

  /// Isobaric quantification configuration
  struct OPENMS_DLLAPI QuantConfig
  {
    bool enabled = false;
    double reporter_mz_tol = 0.002;
    double fold_change_threshold = 1.4;
    ScanConfig follow_up_scan;  ///< Follow-up scan config for quant follow-up MS2, RESOLVED from a name
    bool has_follow_up = false;      ///< whether quantification.follow_up_scan named anything
    std::string follow_up_name;      ///< the referenced name, for diagnostics only
  };

  /// Where IdaLogger writes its five streams.
  ///
  /// This is the FULLY RESOLVED run folder, not a base directory: the host composes it (base
  /// directory + per-run timestamp) and creates it before constructing the engine. IdaLogger
  /// joins its own fixed basenames onto it and never creates a directory itself.
  ///
  /// EMPTY MEANS OPEN NOTHING, and that is load-bearing -- a config with no "runtime" section
  /// or with "runtime": {} opens no streams at all, which several tests rely on. Note the
  /// authored method.json layer gives "" the OPPOSITE meaning ("." = the working directory);
  /// C# resolves that before emitting, so an empty value never reaches here while logging is
  /// on. See ADR-0015 before collapsing the two meanings into one.
  struct OPENMS_DLLAPI RuntimeConfig
  {
    std::string log_dir;
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
    const CharacterizationConfig& characterization() const { return characterization_; }

    /**
     * @brief The complete set of keys any scan object may carry (ms1, ms2[], ms3[], follow_up_scan).
     *
     * This is schema, not an implementation detail: it is what rejectUnknownKeys validates against
     * and what parseScanConfig reads, and the loader's error message already points users at the
     * generated FlashIDA/test-data/config_schema_reference.json. Exposed so the parity test can
     * assert that the C#-generated reference carries every one of them -- the C# emit set and this
     * list drifting apart is precisely how rf_lens/source_cid/source_cid_scaling/scan_rate became
     * unreachable from method.json while C++ went on parsing them (ADR-0011).
     */
    static const std::set<std::string>& scanKeys();

    /// Validate configuration consistency; throws std::invalid_argument on conflict
    void validate() const;

    /// Access the full levels map (for iteration)
    const std::map<int, MSLevelConfig>& levels() const { return levels_; }

  private:
    /**
     * @brief Derive every level's selection state from characterization.mode, after the whole
     *        document has been parsed.
     *
     * This is the mechanism that makes `mode` the single MS3 switch. It replaces the old
     * selection_strategy block, which let each level's on/off state be set independently and so
     * allowed states the engine cannot actually express (MS3 on with MS2 off, for instance).
     *
     * MUST run after the characterization and precursor_selection parses and BEFORE the
     * deconvolution.tol loop, because that loop sizes itself from levels_.
     *
     * MUST assign all THREE levels. Level 1 is the easy one to forget and the expensive one: with
     * MSLevelConfig::selection now defaulting to None, omitting level 1 leaves it None, and
     * FLASHIda.cpp:168 then short-circuits every MS1 before selection -- the instrument acquires
     * nothing at all, silently.
     */
    void applyCharacterizationMode_(const std::string& rank_by, int max_precursors,
                                    int min_precursor_charge);

    // Scan-name resolution (resolving a `<section>.follow_up_scan` reference, validating a
    // user-authored name as an identifier, listing what is defined for an error message) lives in
    // Config.cpp's anonymous namespace, NOT here. Those helpers need nlohmann::json in their
    // signatures, and this header deliberately does not know that type exists -- the public
    // constructor takes a `const std::string&` of raw JSON precisely so the JSON library stays an
    // implementation detail. None of them touches member state, so none needs to be a member.

    std::map<int, MSLevelConfig> levels_;
    DeconvolutionConfig deconv_;
    FAIMSConfig faims_;
    SchedulingConfig scheduling_;
    TargetingConfig targeting_;
    CharacterizationConfig characterization_;
    QuantConfig quant_;
    RuntimeConfig runtime_;
    bool exploration_enabled_ = false;

    /// Names listed in precursor_selection.additional_scans, kept only so the unreferenced-block
    /// warning can tell a dispatched definition from an orphaned one.
    std::set<std::string> additional_scan_names_;

    static const MSLevelConfig default_level_;
  };

} // namespace OpenMS
