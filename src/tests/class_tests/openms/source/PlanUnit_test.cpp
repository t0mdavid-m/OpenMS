// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/CandidateAdmission.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PlanUnit.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>

#include "FLASHIda_TestAccess.h"
#include "FLASHIda_TestHelpers.h"  // loadTsvScans

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;

// wouldBar is the reconcile gate of survey re-planning: it asks whether stamps a unit committed in the
// meantime would have refused a candidate at admission. It is a SECOND statement of admitCandidate's
// exclusion gate and qscore bar, so it is pinned against the real gate rather than against a table of
// expected answers alone -- the drift guard below fails if either side changes without the other.

namespace
{
  const double PROTON = 1.00727646677;
  const double ISO_GAP = 1.00235;

  // SurveyMemory is private to PrecursorSelection; its type is reached through the test seam's
  // accessor rather than by adding a second one, and its members are public.
  using Memory = decltype(FLASHIdaTestAccess::surveyMemory(std::declval<const FLASHIda&>()));

  /// A synthetic MS1 PeakGroup resolving exactly @p charges, three peaks each (one peak would make
  /// getMzRange report the charge as unresolved; see FLASHIda_CandidateAdmission_test).
  PeakGroup makeEnvelope(double mono_mass, const std::vector<int>& charges, double qscore, float snr)
  {
    PeakGroup pg(charges.front(), charges.back(), true);
    pg.setMonoisotopicMass(mono_mass);
    for (int z : charges)
    {
      for (int iso = 0; iso < 3; ++iso)
      {
        FLASHHelperClasses::LogMzPeak lp;
        lp.mz = (mono_mass + z * PROTON + iso * ISO_GAP) / z;
        lp.abs_charge = z;
        pg.push_back(lp);
      }
      pg.setChargeSNR(z, snr);
      pg.setChargeIsotopeCosine(z, 0.99f);
    }
    pg.setQscore(qscore);
    pg.setRepAbsCharge(charges.front());
    return pg;
  }

  MemoryStamp stamp(MemoryStamp::Kind kind, int key, double value, int charge = 0)
  {
    MemoryStamp s;
    s.kind = kind;
    s.key = key;
    s.charge = charge;
    s.value = value;
    return s;
  }

  /// Replays @p stamps into acquisition memory in order, the way PrecursorSelection wrote them in place.
  Memory memoryFrom(const std::vector<MemoryStamp>& stamps)
  {
    using K = MemoryStamp::Kind;
    Memory m;
    for (const auto& s : stamps)
    {
      switch (s.kind)
      {
        case K::AllMassRt: m.all_mass_rt_map_[s.key] = s.value; break;
        case K::MassQscore: m.mass_qscore_map_[s.key] = s.value; break;
        case K::TqMassRt: m.tqscore_exceeding_mass_rt_map_[s.key] = s.value; break;
        case K::TqMzRt: m.tqscore_exceeding_mz_rt_map_[s.key] = s.value; break;
        case K::AuthoredAcquired: m.authored_acquired_rt_map_[{s.key, s.charge}] = s.value; break;
        case K::LegacyTargetChargeErase: break;   // not acquisition memory; no admission gate reads it
      }
    }
    return m;
  }

  /// The real gate, with its context looked up from @p m exactly as PrecursorSelection's selection loop
  /// does it in an exclusion-applying phase. Every threshold is zero, so only exclusion and the
  /// qscore bar can refuse.
  bool admitCandidateRefuses(const Memory& m, const UnitAdmissionKey& k)
  {
    PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, k.score, 10.0f);

    double bar_value = 0.0;
    bool   has_bar = false;
    {
      auto bar_it = m.mass_qscore_map_.find(k.nominal_mass);
      if (bar_it != m.mass_qscore_map_.end()) { bar_value = bar_it->second; has_bar = true; }
    }

    AdmissionContext ctx;
    ctx.has_authored = k.has_authored;
    if (k.has_authored) ctx.authored = {12, 13, 14};
    ctx.anchor_spent = m.authored_acquired_rt_map_.count({k.nominal_mass, k.anchor_charge}) > 0;
    ctx.mass_barred = m.tqscore_exceeding_mass_rt_map_.find(k.nominal_mass) != m.tqscore_exceeding_mass_rt_map_.end();
    ctx.mz_barred = m.tqscore_exceeding_mz_rt_map_.find(k.integer_mz) != m.tqscore_exceeding_mz_rt_map_.end();
    ctx.qscore_bar = has_bar ? &bar_value : nullptr;
    ctx.where = "test";

    AdmissionVerdict v = admitCandidate(pg, k.anchor_charge, k.score, ctx);
    return v.reason == AdmissionReason::Barred || v.reason == AdmissionReason::ScoreNotBetter;
  }

  UnitAdmissionKey key(double score, bool authored = false)
  {
    UnitAdmissionKey k;
    k.nominal_mass = 12351;
    k.integer_mz = 951;
    k.anchor_charge = 13;
    k.score = score;
    k.has_authored = authored;
    return k;
  }

  /// A DDA survey config for PrecursorSelection driven directly: max_precursors 3, and the arming threshold
  /// as the one varying value, so one run arms the mass bars and the other leaves only the qscore bar.
  std::string realSurveyConfig(double tqscore_threshold)
  {
    return std::string(R"JSON({
  "deconvolution": {
    "score_threshold": 0.0,
    "tqscore_threshold": )JSON") + std::to_string(tqscore_threshold) + R"JSON(,
    "min_charge": 4,
    "max_charge": 50,
    "min_mass": 500,
    "max_mass": 50000,
    "tol": [10, 10, 10]
  },
  "flashtnt": { "min_length": 3, "max_length": 8 },
  "faims": { "cv_values": [], "max_cv_skip": 0 },
  "scheduling": {
    "cycle_time": { "enabled": false, "value_ms": 60000 },
    "scan_timeout": { "enabled": true, "value_ms": 30000 },
    "agc_interval_seconds": 30
  },
  "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "consider_all_charges": false,
    "strict_inclusion": false,
    "tie_threshold": 0.1,
    "rank_by": "qscore",
    "max_precursors": 3,
    "precursor_charges": "single"
  },
  "characterization": { "mode": "off", "max_targets": 10 },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000,
             "agc_target": 800000, "max_it": 246 },
    "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }
  },
  "tagging": {},
  "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 }
}
)JSON";
  }
}

START_TEST(PlanUnit, "$Id$")

/////////////////////////////////////////////////////////////

// The bidirectional drift guard. Each row is checked twice: against its expected answer, so the table
// cannot agree vacuously, and against admitCandidate over the same stamps replayed into memory, so a
// change to either gate that the other does not follow turns a row red.
START_SECTION(wouldBar_matches_admitCandidate_on_the_gate_table)
{
  using K = MemoryStamp::Kind;
  struct Row
  {
    std::string name;
    std::vector<MemoryStamp> stamps;
    UnitAdmissionKey k;
    bool expect_bar;
  };

  const std::vector<Row> rows = {
    {"none", {}, key(0.90), false},
    {"AllMassRt carries no bar", {stamp(K::AllMassRt, 12351, 60.0)}, key(0.90), false},
    {"TqMass only", {stamp(K::TqMassRt, 12351, 60.0)}, key(0.90), true},
    {"TqMass on another mass", {stamp(K::TqMassRt, 12000, 60.0)}, key(0.90), false},
    {"TqMz only", {stamp(K::TqMzRt, 951, 60.0)}, key(0.90), true},
    {"TqMass and TqMz", {stamp(K::TqMassRt, 12351, 60.0), stamp(K::TqMzRt, 951, 60.0)}, key(0.90), true},
    {"authored, anchor charge spent", {stamp(K::AuthoredAcquired, 12351, 60.0, 13)}, key(0.90, true), true},
    {"authored, another charge spent", {stamp(K::AuthoredAcquired, 12351, 60.0, 12)}, key(0.90, true), false},
    {"authored ignores both mass bars",
     {stamp(K::TqMassRt, 12351, 60.0), stamp(K::TqMzRt, 951, 60.0)}, key(0.90, true), false},
    {"authored skips the qscore bar", {stamp(K::MassQscore, 12351, 0.95)}, key(0.90, true), false},
    {"qbar above score", {stamp(K::MassQscore, 12351, 0.95)}, key(0.90), true},
    {"qbar below score", {stamp(K::MassQscore, 12351, 0.85)}, key(0.90), false},
    {"qbar equal to score", {stamp(K::MassQscore, 12351, 0.90)}, key(0.90), false},
    {"two qbars, the lower last", {stamp(K::MassQscore, 12351, 0.95), stamp(K::MassQscore, 12351, 0.85)}, key(0.90), false},
    {"two qbars, the higher last", {stamp(K::MassQscore, 12351, 0.85), stamp(K::MassQscore, 12351, 0.95)}, key(0.90), true},
  };

  for (const auto& r : rows)
  {
    const bool bars = wouldBar(r.stamps, r.k);
    const bool real = admitCandidateRefuses(memoryFrom(r.stamps), r.k);
    TEST_EQUAL(r.name + (bars ? ": bars" : ": admits"), r.name + (r.expect_bar ? ": bars" : ": admits"))
    TEST_EQUAL(r.name + (bars ? ": bars" : ": admits"), r.name + (real ? ": bars" : ": admits"))
  }
}
END_SECTION

START_SECTION(wouldBar_empty_never_bars)
{
  TEST_EQUAL(wouldBar({}, key(0.0)), false)
  TEST_EQUAL(wouldBar({}, key(0.90)), false)
  TEST_EQUAL(wouldBar({}, key(0.0, true)), false)
  TEST_EQUAL(wouldBar({}, UnitAdmissionKey{}), false)
}
END_SECTION

// The same guard on a REAL survey, against the whole selection loop rather than admitCandidate alone.
// Survey N admits its species and journals their stamps; those stamps are applied the way a commit will
// apply them (roll the survey back, replay); survey N+1 sees the SAME spectrum at the same retention
// time. A species admitted at N must be refused at N+1 exactly when wouldBar says its stamps bar it.
//
// Run twice so both answers occur: with the mass bars armed every species is barred, and with the
// arming threshold above every qscore only the qscore bar remains, which an equal score does not fail.
START_SECTION(wouldBar_agrees_with_selection_on_a_real_survey)
{
  const auto scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_ecoli_rich.txt");
  ABORT_IF(scans.empty())
  const auto& s = scans.front();

  int barred = 0, readmitted = 0;
  for (double tq : {0.3, 1.0})
  {
    const std::string cfg_json = realSurveyConfig(tq);
    Config cfg{cfg_json};
    Deconvolution deconv(cfg, cfg.toleranceList());
    PrecursorSelection ps(cfg, deconv);

    const int n = ps.filterAndRank(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 1, 0.0);
    TEST_EQUAL(n > 0, true)

    // The admission key of every species admitted at N. The score is the one selection recorded (the
    // MassQscore stamp), not triggerScores(), which is a float and would round across the bar.
    std::vector<UnitAdmissionKey> keys;
    for (int i = 0; i < n; ++i)
    {
      const PeakGroup& pg = ps.selectedPeakGroups()[i];
      UnitAdmissionKey k;
      k.nominal_mass = SpectralDeconvolution::getNominalMass(pg.getMonoMass());
      k.anchor_charge = ps.triggerCharges()[i];
      const auto [lo, hi] = pg.getMzRange(k.anchor_charge);
      k.integer_mz = (int)round((lo + hi) / 2.0);
      const auto sp = ps.journal().species.find(k.nominal_mass);
      ABORT_IF(sp == ps.journal().species.end())
      bool has_score = false;
      for (const auto& st : sp->second.species_level)
      {
        if (st.kind == MemoryStamp::Kind::MassQscore)
        {
          k.score = st.value;
          has_score = true;
        }
      }
      TEST_EQUAL(has_score, true)
      keys.push_back(k);
    }

    // Every stamp of survey N, in write order, applied onto its entry state.
    std::vector<MemoryStamp> stamps = ps.journal().survey_level;
    for (const auto& kv : ps.journal().species)
    {
      const auto f = PrecursorSelection::flatten(kv.second);
      stamps.insert(stamps.end(), f.begin(), f.end());
    }
    std::sort(stamps.begin(), stamps.end(), [](const MemoryStamp& a, const MemoryStamp& b) { return a.seq < b.seq; });
    ps.rollbackAcquisitionMemory();
    ps.applyStamps(stamps);

    const int n2 = ps.filterAndRank(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 1, 0.0);
    std::set<int> admitted_again;
    for (int i = 0; i < n2; ++i)
      admitted_again.insert(SpectralDeconvolution::getNominalMass(ps.selectedPeakGroups()[i].getMonoMass()));

    for (const auto& k : keys)
    {
      const bool bars = wouldBar(stamps, k);
      const bool refused = admitted_again.count(k.nominal_mass) == 0;
      (bars ? barred : readmitted) += 1;
      const std::string name = "tq=" + std::to_string(tq) + " nominal=" + std::to_string(k.nominal_mass);
      TEST_EQUAL(name + (bars ? ": bars" : ": admits"), name + (refused ? ": bars" : ": admits"))
    }
  }
  // Both directions were exercised, so agreement is not vacuous on either side.
  TEST_EQUAL(barred > 0, true)
  TEST_EQUAL(readmitted > 0, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
