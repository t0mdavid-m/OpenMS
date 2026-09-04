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
    Coverage,  ///< Extend sequence coverage
    Exhaustive ///< Fragment every deconvolved mass of the winner MS2 scan, mapped or not (ADR-0023)
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
  /// This is the SOLE source of acquisition geometry (ADR-0021). It used to share the question with
  /// `charge_based_exclusion`, an exclusion-KEYING developer flag: the list `Separate` walked was
  /// built multi-charge only when that flag was on, so `Separate` silently equalled `Single` in every
  /// config that left it at its default. The flag is gone; geometry has one owner.
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

  /**
    @brief True for a metric that weighs bulk signal without matching fragments (ADR-0020).

    A **measuring** metric (MassCount, RemainingPrecursor) scores a variant from peak or
    intensity counts alone, so it never produces a ProteoformMatch and its pre-scans leave no
    identification behind. FragmentCount is the sole **reading** metric: it matches every
    variant against the proteoform in order to score it, so identification is a side effect of
    scoring.

    At MS3 the distinction decides whether an exploration group produces any evidence at all,
    which is why a measuring MS3 sweep must always be closed by a production scan
    (`Exploration.cpp`, the post-winner dispatch). At MS2 it does not apply — every variant is
    matched there regardless of metric.

    `None` is not a sweep, so it is not measuring.
  */
  inline bool isMeasuringMetric(ExplorationMetric metric)
  {
    return metric == ExplorationMetric::MassCount || metric == ExplorationMetric::RemainingPrecursor;
  }

  /**
    @brief Activation-coupled parameter predicates: does @p act give this parameter meaning?

    Definitions live in Config.cpp; declared here because Exploration asks the same question when
    it decides which axis of a sweep to vary and which one its baseline zeroes. It used to
    re-inline the activation literals instead, so the answer had two definitions that could
    disagree silently.

    Mirrored in C# by @c ScanFactory.NeedsReactionTime, which decides whether the ReactionTime key
    is emitted at all. Both sides are pinned by a drift-guard pair — @c Config_SchemaProjection_test
    (C++) and @c ScanFactoryTests (C#) — so change one and you must change the other.
  */
  OPENMS_DLLAPI bool needsCollisionEnergy(const std::string& act);
  OPENMS_DLLAPI bool needsReactionTime(const std::string& act);

  /**
    @brief The smallest ion-ion reaction time the instrument will accept (ms).

    THE TWO COUPLED AXES DO NOT SHARE A "FRAGMENTATION OFF" VALUE, and that asymmetry is the
    instrument's, not ours. A collision energy of 0 is commandable — it simply does not fragment —
    but a reaction time of 0 is **rejected outright**, so 0.03 is what "no reaction" has to mean on
    that axis. 0.03 ms is ~300x shorter than a working ETD time, so it is still an un-fragmented
    reference.

    Two consumers, both of which would otherwise command a value that cannot be acquired:
    Exploration's synthesized ETD/EThcD baseline (ADR-0029 decision 2), and Config::validate's
    rejection of an authored @c reaction_time_min below it.

    Owner-confirmed against the hardware. It is not derivable from anything in this repository —
    the iAPI headers are not distributed and no test can exercise the instrument — so treat it as
    an external fact rather than something to re-derive.
  */
  constexpr double MIN_REACTION_TIME_MS = 0.03;

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

    /// MONITOR SCAN (ADR-0042): a periodic MS1 acquired WHILE THIS LEVEL IS SWEEPING, so the
    /// operator can watch the source. Deconvolved and logged; it selects nothing, excludes nothing,
    /// and no later acquisition decision reads anything it wrote.
    ///
    /// Per level because the two exploration blocks are authored separately and their sweeps differ
    /// in length by an order of magnitude: an MS3 sweep is targets x (points + 2) scans and is the
    /// one that leaves the operator blind, while an MS2 sweep is often over in seconds. Enabling it
    /// on characterization.exploration alone is the intended production shape.
    ///
    /// Off by default, so no committed config emits a monitor scan and no golden moves.
    bool monitor_ms1_enabled = false;
    /// Wall-clock spacing. Config::validate THROWS on <= 0 when enabled: at 0 every drain would mint
    /// another priority-0 MS1 ahead of the sweep and the sweep would never progress.
    double monitor_ms1_interval_ms = 30000.0;
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
    /// Interval between scheduled AGC prescans. This is the ONLY thing that emits one (ADR-0031):
    /// the drained-queue path used to fabricate a prescan as filler AND reset this timer, so the
    /// authored value only ever bound a run whose queue stayed busy for a whole interval. 1 s is
    /// close to the cadence that filler was actually delivering; it is a domain choice, not a
    /// measured one, and no test run lasts long enough for CI to observe it.
    uint64_t agc_interval_ms = 1000;
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
    /// How many charge states of a selected precursor one MS2 acquires -- from
    /// precursor_selection.precursor_charges. Decides the isolation GEOMETRY of a single scan
    /// (ADR-0016), and since ADR-0021 it is the only thing that does. Exclusion is mass-keyed.
    ChargeAcquisitionMode precursor_charges = ChargeAcquisitionMode::Single;
    /// precursor_selection.max_charge_states -- an upper bound on the ACQUISITION CHARGE SET
    /// (anchor + notches), 0 = no bound (ADR-0040). It can only LOWER the cap; the hard ceiling
    /// stays MAX_NOTCHES_PER_STAGE + 1, so the 2048-byte ABI is never at risk.
    ///
    /// It exists because `separate` had no bound of its own: it never writes a notch
    /// (ScanCommandQueue.cpp:316 gates that on Multiplexed) yet was clamped by the notch array's
    /// capacity, so a SCAN COUNT was limited by a geometry constant. Measured on the
    /// separate_charges golden: 10 of 18 species sat at exactly 10 charges, 100 of its 135 MS2.
    int max_charge_states = 0;
    /// precursor_selection.min_charge_states -- the minimum ACQUISITION CHARGE SET size for a
    /// species to be selected at all (ADR-0040). 1 = no floor, and 1 is the default, so this is
    /// inert unless authored. Inert under `single` too, where the set is always size 1.
    ///
    /// Refused in admitCandidate, so a species below the floor costs no max_precursors slot.
    /// A SIBLING is exempt: its set is REDUCED by what the species already isolated (ADR-0036),
    /// so the last charge of a split envelope presents size 1 and the floor would refuse exactly
    /// the charge that ADR needed to recover.
    int min_charge_states = 1;
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
  /// The single MS3 switch. `Off` means no MS3 is emitted at all; the on-values ARE the
  /// objectives, so there is no separate enable flag and no way to express "on but with no
  /// objective". Supersedes ADR-0004, which decided there would be no enable flag at all.
  enum class CharacterizationMode
  {
    Off,
    Ambiguity,
    Coverage,
    /// Target every deconvolved mass of the winner MS2 scan, not only the ones that matched a
    /// theoretical fragment (ADR-0023). An unmatched mass is acquired identically and logged
    /// rather than matched.
    Exhaustive
  };

  /// Which quantification verdicts buy the identification scan (ADR-0039).
  ///
  /// A CUT POINT on the verdict quality ladder, not a set:
  ///   Differential < NotDifferential < IncompleteChannels < ExtractionFailed
  /// IncompleteChannels and ExtractionFailed are deliberately not separate cut points -- both mean
  /// "no usable quant number" and differ only in why, which quant_verdict already reports.
  ///
  /// This replaces the hardcoded `verdict == Differential` that was the whole of the decision until
  /// ADR-0039. Note what it transitively governs: tagging and MS3 targeting are suppressed on a 'Q'
  /// scan and ride the bought 'R', so this enum also decides which species are CHARACTERIZED.
  enum class QuantIdentify
  {
    /// Differential only, subject to `enriched_in`. THE DEFAULT, and exactly ADR-0038's behaviour.
    Differential,
    /// Differential | NotDifferential -- any cleanly measured species.
    Quantified,
    /// Every screened precursor, verdict irrelevant: quantification becomes pure annotation.
    All,
    /// Nothing is ever bought. The run measures and never identifies, which halves the per-precursor
    /// cost and is what a labelled survey wants when identification comes from elsewhere.
    None
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

    /// Pool floor for characterization.mode == Exhaustive: a deconvolved mass below this (Da) is
    /// not a target. 0 = off, and off is the default deliberately -- the mode does exactly what
    /// its name says until told otherwise.
    ///
    /// NOT inheritable from deconvolution.min_mass: that floor is not applied to MSn output. The
    /// reference config sets min_mass 500 / min_charge 4 and its MS2 spectra still contain 248 Da
    /// and charge-1 species, so this is a genuinely new floor rather than a duplicate of one that
    /// already reaches here.
    double min_target_mass = 0.0;

    /// How many charge states of a target FRAGMENT one MS3 acquires -- from
    /// characterization.fragment_charges. Replaces the bool `ms3_all_charges`, whose two states are
    /// now Single (the fragment's best-MS2 charge) and Separate (one MS3 per observed charge);
    /// Multiplexed co-isolates them into one MS3, which is the free direction because synchronous
    /// precursor selection is one simultaneous waveform rather than N sequential fills.
    ChargeAcquisitionMode fragment_charges = ChargeAcquisitionMode::Single;
  };

  /// Isobaric quantification configuration (ADR-0038).
  ///
  /// The QUANTIFICATION SCAN is `ms_settings.ms2_quant`: rostered once per selected precursor,
  /// labelled `'Q'`, and the only scan whose reporter ions are ever measured. A differential
  /// verdict then BUYS the IDENTIFICATION scan, `ms_settings.ms2`, labelled `'R'` -- which is why
  /// `ms2` leaves the unconditional roster whenever this is enabled.
  ///
  /// The reverse arrangement -- measure the base MS2, buy a "quant follow-up" -- is what shipped
  /// until ADR-0038 and could never work: it screened for reporter ions on a scan whose activation
  /// (ETD, in the only config that enabled the feature) cannot release them, then acquired the HCD
  /// scan that could and never looked at it.
  struct OPENMS_DLLAPI QuantConfig
  {
    /// One of the two experimental conditions the fold change is taken over. Channels are held as
    /// ORDINALS into the selected method's getChannelInformation(), resolved from the authored
    /// names once at load, so no consumer downstream needs the name->m/z table -- and so an
    /// unknown channel name fails at load rather than silently reading the wrong intensity.
    struct Condition
    {
      std::string name;
      std::vector<size_t> channels;
    };

    bool enabled = false;
    /// One of the seven isobaric methods OpenMS ships, spelled as
    /// TopDownIsobaricQuantification's own valid strings. `"none"` is deliberately NOT among them:
    /// `enabled` is the switch, so `enabled: true, labelling: "none"` is unwritable.
    std::string labelling = "tmt6plex";
    double reporter_mz_tol = 0.002;
    double fold_change_threshold = 1.4;
    /// ADR-0039. Which verdicts buy `identification_scan`. The default reproduces ADR-0038 exactly.
    QuantIdentify identify = QuantIdentify::Differential;
    /// ADR-0039. Ordinal into `conditions` that a species must be ENRICHED IN for a Differential
    /// verdict to buy the identification scan. `-1` = either direction (the default, and the
    /// symmetric test ADR-0038 shipped). Resolved from the authored CONDITION NAME at load, so no
    /// consumer downstream needs the name table -- the same treatment channel names get.
    ///
    /// Named by condition, NEVER as up/down. fold_change = mean(conditions[0]) / mean(conditions[1]),
    /// so "up" would mean *enriched in conditions[0]* -- on method_quant.json, enriched in CONTROL,
    /// the opposite of what anyone writing "up" intends -- and would silently invert the experiment
    /// the moment someone reordered the array.
    ///
    /// ⚠ Consumers must read `Result::condition_means`, NEVER `fold_change`: a wholly-absent
    /// condition is Differential with fold_change == -1.0, a SENTINEL rather than a ratio, so a
    /// `fold_change > threshold` direction test silently drops every species present in one
    /// condition and absent in the other -- the strongest result the experiment can produce.
    int enriched_in = -1;
    /// EXACTLY TWO, and the ARRAY ORDER IS THE RATIO DIRECTION:
    ///   fold_change = mean(conditions[0]) / mean(conditions[1])
    /// Authored as a JSON array rather than an object because nlohmann's object_t is a std::map,
    /// so an object would sort the two conditions alphabetically and silently decide which is the
    /// numerator -- the same trap Config.cpp already documents for additional_ms2 dispatch order.
    std::vector<Condition> conditions;
    /// Isotope-impurity correction matrix, one entry per channel, Thermo data-sheet spelling.
    /// EMPTY = use the selected method's stock matrix (IsobaricQuantifier's own default, which
    /// corrects impurity and does nothing else because `normalization` defaults false). An
    /// all-zero matrix yields the identity, i.e. correction off.
    std::vector<std::string> correction_matrix;
    /// `ms_settings.ms2`: the scan a differential verdict buys. In a quant config it is
    /// deliberately NOT on levels_[2].scans -- the quantification scan holds that slot, and that
    /// swap is the whole of ADR-0038.
    ///
    /// Both scan configs below are copied here UNCONDITIONALLY, whatever `enabled` says. Only the
    /// ROSTER decision is conditional. That split is what lets the generated schema reference --
    /// which has every key populated but quantification off -- still assert that ms2_quant's keys
    /// round-trip; tying the copy to `enabled` would have made the only new scan slot in this
    /// schema unreachable from the one fixture whose job is to prove keys reachable.
    ScanConfig identification_scan;
    /// `ms_settings.ms2_quant`: the quantification scan, rostered and measured when enabled.
    ScanConfig quantification_scan;
    /// Whether `ms_settings.ms2_quant` was authored. Distinct from "levels_[2].scans is non-empty",
    /// which stays true when ms2_quant is absent but ms2 is present -- the exact state validate()
    /// has to reject, since it is plain DDA that silently quantifies nothing.
    bool has_quant_scan = false;
    /// Whether `ms_settings.ms2` was authored (ADR-0039). Same reason `has_quant_scan` exists:
    /// `identification_scan` is copied unconditionally and is DEFAULT-CONSTRUCTED when absent, so
    /// its contents cannot answer "was one authored". validate() has to reject that state -- a
    /// quant config with no ms2 builds the bought scan from defaults, which is latent while only a
    /// Differential verdict reaches it and fires on every precursor under `identify: "all"`.
    bool has_identification_scan = false;
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
