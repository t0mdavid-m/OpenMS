// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3WitnessModel.h>

#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace OpenMS
{
  /**
    @brief The W gate of one MS3 plan (ADR-0076): the inputs a candidate MS3 is scored on, the walk's running rank, and
    the records the [MS3-W] line is written from.

    Every input is the value the engine's own logs show for it, read back at the precision the log prints it at: W was
    fitted on values read from scan_commands.tsv, scan_results.tsv and identification.tsv, and its split thresholds sit
    at values as those files print them (the ms3_esm package's ENGINE.md and inputs_from_log.py, 8.2). The parent-level
    values are filled once per plan where the parent MS2 is finished (FLASHIda::finishRegularMS2_); the per-candidate
    values where the planner walks a target's charges.
  */

  /// @p v as a default-formatted ostream prints it (6 significant digits), read back.
  OPENMS_DLLAPI double asLoggedSig6(double v);
  /// @p v as std::fixed / setprecision(4) prints it, read back.
  OPENMS_DLLAPI double asLoggedFixed4(double v);

  /// The residue letters of a ProForma string: the leading <...> global modification removed, every [...] removed, and
  /// every other non-letter (parentheses of a ranged shift) dropped.
  OPENMS_DLLAPI std::string proformaLetters(const std::string& proforma);
  /// 1 if a fragment's ProForma string, without its leading <...>, writes a mass shift ('['), else 0.
  OPENMS_DLLAPI bool proformaHasShift(const std::string& proforma);
  /// Whether a fragment of residues @p fragment_letters lies on the parent's identified stretch @p parent_letters at
  /// its terminus: the first L residues for an N-terminal (a, b, c) fragment, the last L for a C-terminal one.
  OPENMS_DLLAPI bool witnessPlaced(const std::string& parent_letters, const std::string& fragment_letters, bool n_terminal);

  /// The proteoform an MS2 identification row names: (protein, start_pos, end_pos), compared as identification.tsv
  /// writes them. The key of n_ms2_prior and n_ms3_prior.
  struct OPENMS_DLLAPI WitnessPriorKey
  {
    std::string protein;
    int start = -1;
    int end = -1;
    bool operator<(const WitnessPriorKey& o) const { return std::tie(protein, start, end) < std::tie(o.protein, o.start, o.end); }
  };

  /// The engine's running counts behind n_ms2_prior and n_ms3_prior (FLASHIda, under analysis_mutex_): MS2
  /// identifications applied per key, and MS3 received per key they were placed under at planning.
  struct OPENMS_DLLAPI WitnessPriorCounts
  {
    std::map<WitnessPriorKey, int> ms2;
    std::map<WitnessPriorKey, int> ms3;
    std::map<int, WitnessPriorKey> planned_ms3;   ///< scan id -> key, for a planned MS3 placed under a key; until received
  };

  /// The last MS1 whose scan_results row was written (any kind: a survey, a Monitor or a Window MS1, any CV), as the
  /// row prints it: the time, and each mass with its per-charge intensities summed.
  struct OPENMS_DLLAPI WitnessLastMs1
  {
    uint64_t resolve_ts = 0;                       ///< steady-clock ms; 0 = no MS1 resolved yet
    std::string tracking_id;
    std::vector<std::pair<double, double>> masses; ///< (mass, summed intensity), each as the row prints it
  };

  /// x_ms1_plan: the largest summed intensity of a mass of @p ms1 within 10 ppm of @p mono_mass at an offset of 0 or
  /// +-1.00335 Da; 0 when none matches.
  OPENMS_DLLAPI double witnessMs1Plan(const WitnessLastMs1& ms1, double mono_mass);

  /// The parent-level inputs of one plan.
  struct OPENMS_DLLAPI WitnessPlanContext
  {
    WitnessRaw parent;            ///< z_ms2, mass_ms2, x_prec_*, id_score, x_ms1_plan, x_ms2_*, x_plan_age, x_rt_frac
    bool has_key = false;
    WitnessPriorKey key;
    std::string parent_letters;   ///< the parent's identified stretch, letters only (placement)
    int prior_ms2 = 0;            ///< MS2 identifications applied under key before this plan; the parent is added
    int prior_ms3 = 0;            ///< MS3 received under key before this plan
    std::string parent_id;        ///< the parent MS2's tracking id
    const WitnessPriorCounts* counts = nullptr;   ///< the engine's running counts, read when the key is known
  };

  /// One scored candidate charge, for the [MS3-W] line and, when admitted, the command it becomes.
  struct OPENMS_DLLAPI WitnessRecord
  {
    std::string ion_type;
    int ion_index = 0;
    int charge = 0;
    WitnessRaw raw;
    double w = -1.0;              ///< -1: untested
    int uncovered = -1;           ///< uncovered bonds inside the fragment at planning; -1 unknown
    std::string verdict;          ///< admitted | refused | untested | capped | withdrawn | kept
    std::string reason;           ///< below_threshold | missing_input=<name> | cleavage_cap | budget
    bool placed = false;          ///< placed on the parent's stretch (its MS3 is counted under the key when received)
    WitnessPriorKey key;
  };

  struct OPENMS_DLLAPI WitnessGate
  {
    const WitnessModel* model = nullptr;
    double threshold = 0.0;
    bool rank_by_w = true;        ///< rank_by yield: the walk by W
    std::string gate_name = "characterization";
    WitnessPlanContext ctx;
    double cid3 = std::numeric_limits<double>::quiet_NaN();   ///< the MS3 stage's collision energy (first variant)
    int ms3_per_target = 1;       ///< commands per admitted charge (CE variants under MS3 exploration)
    /// ce1 of a raw-mass target (exhaustive, discriminating): the winner MS2's collision energy, the stage 0 it replays.
    double default_ce1 = std::numeric_limits<double>::quiet_NaN();
    double queue_ahead = 0.0;     ///< q0: MS3 queued at planning, the proteoform's own unsent ones excluded
    int admitted = 0;             ///< MS3 commands admitted so far in this plan (the walk's k)
    /// A target fragment's ProForma string, without the fixed-modification prefix; "" when it cannot be rendered.
    std::function<std::string(const std::string& ion_type, int ion_index)> fragment_proforma;
    /// Uncovered bonds inside a target fragment (logged only); -1 unknown. Set by the planner that has the model.
    std::function<int(const std::string& ion_type, int ion_index)> uncovered_bonds;
    std::vector<WitnessRecord> records;

    /// One candidate charge scored at rank @p rank: its raw values, W (or untested with the first missing input), its
    /// verdict against the threshold. Records nothing.
    WitnessRecord score(const std::string& ion_type, int ion_index, int charge, double ce1,
                        const FragmentAnalysis::FragmentScores& s1, int rank) const;
  };

  /// The [MS3-W] line for @p records: @p stage plan | rescore, the parent-level inputs of @p g once, each record's
  /// candidate group. @p extra (" ms1=<id>") is appended after the stage. Empty when @p records is empty.
  OPENMS_DLLAPI std::string witnessLine(const WitnessGate& g, const char* stage, const std::string& extra,
                                        const std::vector<WitnessRecord>& records);

  /// What a queued W-scored MS3 is re-scored with at an MS1 (ScanCommandQueue::registerWitnessed): its record (the
  /// raw values it was admitted at), the gate's threshold and name, its plan's parent and precursor, the model.
  struct OPENMS_DLLAPI WitnessRegistration
  {
    std::shared_ptr<const WitnessRecord> record;
    std::shared_ptr<const WitnessGate> gate;   ///< the plan's gate, for the threshold and the [MS3-W] parent fields
    int group = 0;                             ///< the commands withdrawn together (a target's CE variants); 0 = alone
    int precursor_id = 0;
  };
} // namespace OpenMS
