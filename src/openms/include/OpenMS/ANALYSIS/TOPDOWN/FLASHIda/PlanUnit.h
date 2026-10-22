// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include <map>
#include <vector>

namespace OpenMS
{
  /// @brief The plain data of survey re-planning: what a survey's selection wrote, and the plan it hands the queue.
  ///
  /// A survey's MS2 decisions stay provisional while they wait in FLASHIda's queue; a newer survey at the same CV
  /// re-plans them, and only the drain's dequeue makes them permanent. That needs acquisition memory to be written
  /// as a JOURNAL rather than in place, so the writes of a species can be applied when it commits and discarded
  /// when it is dropped. These types are that journal and the unit the queue keeps it in.
  ///
  /// Header-only and free of `Config`, in the shape of `CandidateAdmission.h`: every value arrives already decided.

  /// One write PrecursorSelection used to make in place on admission (ADR-0028 / dynamic exclusion), now journaled.
  struct MemoryStamp
  {
    enum class Kind
    {
      AllMassRt,
      MassQscore,
      TqMassRt,
      TqMzRt,
      AuthoredAcquired,
      LegacyTargetChargeErase
    };
    Kind   kind = Kind::AllMassRt;
    int    key = 0;            ///< nominal mass, or integer m/z for TqMzRt, or legacy target mass bin
    int    charge = 0;         ///< AuthoredAcquired / LegacyTargetChargeErase only; 0 = species-level
    double value = 0.0;        ///< rt (seconds) or score
    double legacy_mass = 0.0;  ///< LegacyTargetChargeErase: the exact target_mass_charge_map_ key
    int    seq = 0;            ///< write order within the survey; replay order (mass_qscore_map_ is an overwrite)
  };

  /// Everything one survey's selection wrote for one species (keyed by nominal mass; ADR-0036 siblings fold in).
  struct SpeciesStamps
  {
    int nominal_mass = 0;
    std::vector<MemoryStamp> species_level;             ///< applied on the species' FIRST commit
    std::map<int, std::vector<MemoryStamp>> per_charge; ///< charge -> stamps applied when THAT charge's command commits
  };

  /// What selection returns besides selectedPeakGroups(): per-species stamps + survey-level ones owned by no unit.
  struct SelectionJournal
  {
    std::map<int, SpeciesStamps> species; ///< nominal -> stamps, admitted species only
    std::vector<MemoryStamp> survey_level; ///< stamps of candidates that passed exclusion but were refused
  };

  /// Plain = R roster; Quant = has a Q group; Sweep = exploration groups.
  enum class UnitRole
  {
    Plain,
    Quant,
    Sweep
  };

  /// One species of one survey's plan, at one CV.
  struct PlanUnit
  {
    int unit_id = 0;
    double faims_cv = 0.0;
    int nominal_mass = 0;
    int decision_survey_id = 0;       ///< tracking id of the MS1 that (last) chose it
    UnitRole role = UnitRole::Plain;
    bool separate = false;            ///< precursor_charges == separate (remainder may be re-planned)
    std::vector<ScanCommand> commands; ///< the unit's MS2s / sweep variants, in push order
    std::vector<double> purity;        ///< index-parallel to commands: windowPurity on the decision survey
    std::vector<int> quant_group_ids;  ///< Quantification groups to close if dropped
    std::vector<int> sweep_group_ids;  ///< Exploration groups to drop if dropped
    SpeciesStamps stamps;
    // state, owned by the queue:
    bool committed = false;
    std::vector<int> sent_charges;
    int commands_left = 0;            ///< set on insert; a unit with none left is erased
  };

  /// Admission key the reconcile step tests a new unit against (mirrors CandidateAdmission's exclusion + qscore gates).
  struct UnitAdmissionKey
  {
    int nominal_mass = 0;
    int integer_mz = 0;     ///< anchor's rounded centre m/z
    int anchor_charge = 0;
    double score = 0.0;
    bool has_authored = false;
  };

  /// Would the committed stamps @p committed have refused @p k at admission?
  ///
  /// MUST mirror CandidateAdmission.h admitCandidate's exclusion gate and its qscore bar (the `Barred` and
  /// `ScoreNotBetter` refusals). Pinned against admitCandidate by a bidirectional drift-guard table in
  /// PlanUnit_test: a change to either gate that the other does not follow turns a row red.
  inline bool wouldBar(const std::vector<MemoryStamp>& committed, const UnitAdmissionKey& k)
  {
    using K = MemoryStamp::Kind;
    bool mass_barred = false, mz_barred = false, anchor_spent = false;
    const double* qbar = nullptr;
    for (const auto& s : committed)
    {
      if (s.kind == K::TqMassRt && s.key == k.nominal_mass) mass_barred = true;
      if (s.kind == K::TqMzRt && s.key == k.integer_mz) mz_barred = true;
      if (s.kind == K::AuthoredAcquired && s.key == k.nominal_mass && s.charge == k.anchor_charge) anchor_spent = true;
      if (s.kind == K::MassQscore && s.key == k.nominal_mass) qbar = &s.value; // last write wins, as the map's overwrite
    }
    if (k.has_authored) return anchor_spent; // an authored species reads only its per-charge bar
    if (mass_barred || mz_barred) return true;
    return qbar != nullptr && k.score < *qbar; // ScoreNotBetter
  }
} // namespace OpenMS
