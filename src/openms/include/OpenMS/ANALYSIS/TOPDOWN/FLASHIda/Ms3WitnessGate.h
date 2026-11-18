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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms1Buffer.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms2Params.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3WitnessModel.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include <array>
#include <chrono>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
  /**
    @brief The W gate of one MS3 plan (ADR-0079, package 17 of the ms3_esm W gate).

    W is computed once, when the plan is made, for every eligible candidate of the planner's pool; the candidates are
    ranked by W and the plan takes, from the top, those with an uncovered bond and W >= threshold until max_targets is
    filled (the cleavage cap as before); they are sent in that order and the queue is first in, first out. Nothing is
    re-scored afterwards.

    Every input is the value the engine's logs show for it, read back at the precision they print it at (LogFormat.h):
    the candidate's own cells from its command composed before the decision (ScanCommandQueue::composeMS3), the parent
    MS2's values as its scan_results.tsv and identification.tsv rows print them, the last MS1 row's masses; and the eight
    peak features of the precursor from the MS1 buffer (Ms1Buffer.h), at full precision.
  */

  /// @p v as a default-formatted ostream prints it (6 significant digits), read back.
  OPENMS_DLLAPI double asLoggedSig6(double v);
  /// @p v as std::fixed / setprecision(4) prints it, read back.
  OPENMS_DLLAPI double asLoggedFixed4(double v);

  /// The residue letters of a ProForma string: the leading <...> global modification removed, every [...] removed, and
  /// every other non-letter (parentheses of a ranged shift) dropped.
  OPENMS_DLLAPI std::string proformaLetters(const std::string& proforma);
  /// L of package 17: the capital letters of a ProForma string after removing its leading <...> and every [...].
  OPENMS_DLLAPI int proformaResidueCount(const std::string& proforma);
  /// 1 if a fragment's ProForma string, without its leading <...>, writes a mass shift ('['), else 0.
  OPENMS_DLLAPI bool proformaHasShift(const std::string& proforma);

  /// Stage @p stage (0 = the first ';'-group, -1 = the last) of a logged command cell, its first ','-value (the anchor of
  /// co-isolation notches), read back; NaN when empty.
  OPENMS_DLLAPI double loggedStage(const std::string& cell, int stage);

  /// The last MS1 whose scan_results row was written (any kind: a survey, a Monitor or a Window MS1, any CV), as the
  /// row prints it: the time, and each mass with its per-charge intensities summed.
  struct OPENMS_DLLAPI WitnessLastMs1
  {
    uint64_t resolve_ts = 0;                       ///< engine-clock ms; 0 = no MS1 resolved yet
    std::string tracking_id;
    std::vector<std::pair<double, double>> masses; ///< (mass, summed intensity), each as the row prints it
  };

  /// x_ms1_plan: the largest summed intensity of a mass of @p ms1 within 10 ppm of @p mono_mass at an offset of 0 or
  /// +-1.00335 Da; 0 when none matches.
  OPENMS_DLLAPI double witnessMs1Plan(const WitnessLastMs1& ms1, double mono_mass);

  /// The parent MS2's values its log rows print (scan_results.tsv, the first identification.tsv MS2 row): the parent
  /// side of the 29 raw inputs.
  struct OPENMS_DLLAPI WitnessParentRow
  {
    double id_score = std::numeric_limits<double>::quiet_NaN();          ///< identification_score, fixed 4
    double ms2_ions = std::numeric_limits<double>::quiet_NaN();          ///< entries of ms2_fragments (an empty list: 1)
    double mass_count = std::numeric_limits<double>::quiet_NaN();
    double tic_coverage = std::numeric_limits<double>::quiet_NaN();
    double precursor_share = std::numeric_limits<double>::quiet_NaN();
    double fragment_count = std::numeric_limits<double>::quiet_NaN();
    double tag_count = std::numeric_limits<double>::quiet_NaN();
  };

  /// The plan-level part of the gate: what every candidate of one plan shares.
  struct OPENMS_DLLAPI WitnessPlanContext
  {
    std::string parent_id;                 ///< the parent MS2's tracking id (ms3_gate.tsv plan_id)
    std::string survey_id;                 ///< the tracking id the parent MS2 was picked from (its parent_tracking_id)
    std::string parent_path = "regular";   ///< regular | sweep_winner
    WitnessParentRow parent;
    double x_ms1_plan = std::numeric_limits<double>::quiet_NaN();   ///< NaN: no MS1 resolved at or before the plan
    std::string last_ms1_id;
    uint64_t last_ms1_ts = 0;
    uint64_t plan_ts = 0;                  ///< t: the engine clock when the plan was made (ms)
    std::chrono::steady_clock::time_point eval_start;
    double mass_ms2 = std::numeric_limits<double>::quiet_NaN();     ///< M: stage 0 of mono_mass as logged
    PeakFeatures pk3{};
    // The buffer selection of the plan, as logged.
    std::string buffer_kind_from = "none";
    double scan_range_lo = 500.0, scan_range_hi = 1800.0;
    uint64_t buffer_surveys = 0, buffer_centroids = 0;
    uint64_t buffer_first_resolve_ts = 0, buffer_last_resolve_ts = 0;
    std::string buffer_last_survey_id;
  };

  /// The uncovered and covered bonds of a target fragment (bond j between residues j and j+1 of the fragment, 1..L-1).
  struct OPENMS_DLLAPI WitnessBonds
  {
    int n_uncovered = -1;                  ///< -1: unknown
    std::vector<int> uncovered, covered;
  };

  /// One row of ms3_gate.tsv: a candidate target of the plan and one of its MS3 commands (a CE sweep's variant).
  struct OPENMS_DLLAPI WitnessRecord
  {
    int candidate = 0;                     ///< 1-based, in the planner's order
    int variant = 0;                       ///< 0-based CE variant
    int planner_rank = 0;
    int w_rank = 0;                        ///< 1-based rank by W among the plan's scored targets; 0 = not ranked
    std::string ion_type;
    int ion_index = 0;
    ScanCommand cmd{};                     ///< the composed command (scan_id 0 until sent)
    std::string ms3_proteoform;
    std::array<double, kW17RawCount> raw{};
    std::array<float, kW17InputCount> x{};
    double w = std::numeric_limits<double>::quiet_NaN();
    double target_w = std::numeric_limits<double>::quiet_NaN();   ///< the largest W over the target's variants
    std::string failed_input;              ///< the first input the gate could not fill; "" = none
    WitnessBonds bonds;
    std::string decision;                  ///< send | discard
    std::string reason;                    ///< sent | below_threshold | no_uncovered_bond | failed | not_chosen | capped | not_eligible:<gate>
    int in_sent = 0;                       ///< the cleavage was sent for this proteoform before (governed plans)
    int room_before = -1;                  ///< the cleavage cap's room when the target was decided; -1 = not consulted
    std::string tracking_id;               ///< set when this command is sent
  };

  /// What a candidate's commands are composed from: the planner's target descriptors.
  struct OPENMS_DLLAPI WitnessCandidate
  {
    std::string ion_type;
    int ion_index = 0;
    int charge = 0;
    double frag_mz = 0.0;
    double iso_width = 0.0;
    FragmentAnalysis::FragmentScores stage1_scores;
    Ms2Params stage0_params;
  };

  struct OPENMS_DLLAPI WitnessGate
  {
    const WitnessModel* model = nullptr;
    double threshold = 0.0;
    bool apply = true;                     ///< false: shadow mode -- decided and logged, the plan sent as without the gate
    bool rank_by_w = true;                 ///< rank_by yield: the plan is ranked by W; intensity: the planner's order, W tests
    int max_targets = 0;
    WitnessPlanContext ctx;
    /// The commands a candidate would be (one per CE variant), composed without a tracking id, with their ProForma.
    std::function<std::vector<std::pair<ScanCommand, std::string>>(const WitnessCandidate&)> compose;
    /// The target's bonds from the planner's coverage model.
    std::function<WitnessBonds(const std::string& ion_type, int ion_index)> bonds;
    std::vector<WitnessRecord> records;
    int64_t eval_us = 0;                   ///< the plan's evaluation time, set when it is decided
    int candidates = 0;                    ///< targets scored

    /// Score one candidate target: compose its commands, fill one record per command (its raw values from the command's
    /// logged cells, the plan context and the peak features; W; or the first input it could not fill). Returns the
    /// index of its first record.
    size_t scoreTarget(const WitnessCandidate& c, int planner_rank);
  };
} // namespace OpenMS
