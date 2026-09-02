// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/NotchSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <algorithm>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{
  /// @brief The per-candidate acquisition decision, lifted out of PrecursorSelection's selection loop.
  ///
  /// Answers, for one MS1 PeakGroup: may it be acquired, at which anchor charge, at which set of
  /// charges, and does it cost a `max_precursors` slot.
  ///
  /// It lives in a header of free functions over plain data — the shape `NotchSelection.h` already
  /// uses, and for the same two reasons. First, the selection loop it came from is ~450 lines whose
  /// `break`/`continue` targets are load-bearing (see the comment at PrecursorSelection.cpp's
  /// `charges_to_process` loop), so a decision that can be read and tested apart from that control
  /// flow is worth more than one that is merely commented. Second, `PrecursorSelection` and
  /// `ScanCommandQueue::buildMS2` must agree on the acquisition charge set or the wire disagrees with
  /// every record of itself — sharing the extraction rather than the policy is what makes them agree
  /// by construction.
  ///
  /// Deliberately free of `Config`, exactly as `NotchSelection.h` is: every knob arrives as a plain
  /// value the caller has already read, and every acquisition-memory lookup arrives as an answer the
  /// caller has already performed. That is what makes these functions drivable from a test with a
  /// synthetic PeakGroup and no engine, no spectrum and no config document.
  ///
  /// @note Charge sets here are ACQUISITION charges — what a scan isolates. The distinction from the
  /// row's AUTHORED charge set (a restriction, ADR-0028) and from the INTENDED charge set (what
  /// "acquired in full" is measured against, ADR-0036) is the one in `CONTEXT.md`; the names are used
  /// in that sense throughout.

  /// What one survey has already done to one species, keyed by nominal mass by the caller.
  ///
  /// A **split envelope** (ADR-0036) is one Precursor arriving as several PeakGroups, each carrying
  /// part of its charge range — which happens only for inclusion targets, because the overlap
  /// collapse that merges near-identical features is deliberately skipped for targeted ones. This
  /// record is what lets the second PeakGroup of a species be told apart from the first, and it is
  /// call-local: a survey starts with an empty ledger.
  struct SpeciesSurveyRecord
  {
    /// Charges that already-SELECTED PeakGroups of this species RESOLVED — not merely acquired.
    /// The distinction matters: a sibling earns a scan only by supplying something no earlier
    /// PeakGroup could, and "could" is about resolution, not about which charge won the anchor.
    std::vector<int> resolved;

    /// Charges this species has actually ISOLATED so far this survey. Subtracted from a sibling's
    /// acquisition set so the scans PARTITION the envelope rather than overlapping it.
    std::vector<int> acquired;

    /// The row's AUTHORED charge set, empty for an unrestricted row. Carried so the caller can say,
    /// after the candidate loop has ended, which named charges the survey still owed — see the
    /// budget-exhausted diagnostic in `PrecursorSelection`. Not read by any decision here.
    std::vector<int> authored;

    /// The monoisotopic mass of the PeakGroup that opened this record, for that same diagnostic.
    /// The ledger keys on the NOMINAL mass, which is an integer bin and not a mass to print.
    double mono_mass = 0.0;
  };

  /// nominal mass -> what this survey has done to it. Call-local; empty at the start of each survey.
  using SurveyLedger = std::map<int, SpeciesSurveyRecord>;

  /// Why a candidate was admitted or refused. Doubles as the `[CHARGE-SET]` reason vocabulary, so a
  /// refusal can say which gate it failed rather than being reported as a gate it never reached.
  enum class AdmissionReason
  {
    Admitted = 0,        ///< acquired; the anchor and the acquisition charge set are in the verdict
    NoAnchorResolved,    ///< no candidate charge survived: every one absent, spent, or below the floor
    IntendedSetComplete, ///< a SIBLING that adds nothing: every charge it offers was already resolved
    BelowMinCharge,      ///< the anchor is under the level's charge floor
    ScoreBelowThreshold, ///< score < qscore_threshold
    SnrBelowThreshold,   ///< the anchor's own envelope does not rise above noise
    Barred,              ///< dynamic exclusion: a tqscore_exceeding_* bar, or a spent authored charge
    ScoreNotBetter,      ///< the species has been acquired at a score this survey does not beat
    BelowMinChargeStates ///< the acquisition charge set is smaller than min_charge_states (ADR-0040)
  };

  /// The anchor charge of one candidate, with the score that will be logged and excluded on.
  ///
  /// Split from the admission decision below because the real dependency order forces it: the
  /// geometry, the inclusion-target match and the within-survey m/z dedup all sit between the two and
  /// all need the anchor. Collapsing them into one call would mean either recomputing the anchor or
  /// moving those three, and both are changes of behaviour dressed as refactors.
  struct AnchorChoice
  {
    int    charge = -1;                                  ///< < 0 means "no anchor" — see @p reason
    double score = 0.0;                                  ///< the anchor's OWN score, never a sibling's
    AdmissionReason reason = AdmissionReason::NoAnchorResolved;

    bool resolved() const { return charge >= 0; }
  };

  /// Everything the anchor pick needs, as plain values the caller has already looked up.
  struct AnchorContext
  {
    /// The row's AUTHORED CHARGE SET (ADR-0028). Empty means the row named no charge, which is
    /// "no opinion" rather than "every charge": the anchor then comes from the envelope.
    std::vector<int> authored;

    /// `authored_acquired_rt_map_` for this nominal mass — the named charges already spent inside the
    /// retention-time window. Never consulted when @c authored is empty, since exclusion is
    /// mass-keyed for every species without an authored set.
    std::vector<int> spent_charges;

    /// What this survey has already done to this species, or nullptr if this is its first PeakGroup.
    ///
    /// Non-null makes this candidate a SIBLING, and a sibling is judged on what it can ADD rather
    /// than on its score (ADR-0036). Before that ADR this was a bare "already acquired this survey"
    /// flag that refused every sibling outright, which is what lost the charges only a sibling could
    /// resolve — under `single` for three surveys, and permanently if the peak stopped eluting
    /// first.
    const SpeciesSurveyRecord* seen = nullptr;

    /// Neither `fan_out` nor an SNR gate here, deliberately. The acquisition mode decides how many
    /// charges one SCAN isolates — that is `AdmissionContext::fan_out`'s job further down — and it
    /// has no say in what a species is OWED, because only an authored row owes charges at all
    /// (ADR-0036, amended 2026-08-28). SNR likewise only ever gated the notches, never the named
    /// set: under an authored row SNR ORDERS the anchor candidates and does not admit them.
    int    min_charge = 0;           ///< the level's charge floor; 0 = no floor
    bool consider_all_charges = false;  ///< rank across the envelope rather than the representative charge
  };

  /// Everything the admission decision needs, past the anchor.
  struct AdmissionContext
  {
    bool has_authored = false;       ///< the row named charges; mirrors `!AnchorContext::authored.empty()`
    std::vector<int> authored;

    /// The named charges already spent inside the retention-time window, for the NOTCH filter.
    /// Read unconditionally, matching the anchor pick.
    std::vector<int> spent_charges;

    /// Is the ANCHOR's own charge already spent? Separate from @c spent_charges above, and not
    /// derivable from it, because the two are gated differently in the loop this came from: the
    /// anchor pick and the notch filter consult the per-charge map unconditionally, while the
    /// exclusion bar consults it only in the selection phases that apply exclusion at all. Deriving
    /// one from the other would quietly apply exclusion in a phase that does not.
    bool anchor_spent = false;

    /// Dynamic exclusion, already looked up: the nominal-mass bar and the integer-m/z bar. They are
    /// ORed, and both are armed together, so exempting one alone changes nothing whenever the anchor
    /// charge repeats.
    bool mass_barred = false;
    bool mz_barred   = false;

    /// `mass_qscore_map_` for this nominal mass, or nullptr if the species has not been acquired
    /// inside the retention-time window. The best score it has been acquired at.
    const double* qscore_bar = nullptr;

    /// As AnchorContext::seen — non-null makes this candidate a sibling. It changes two things here:
    /// the sibling's acquisition set is reduced by what the species already isolated (so the scans
    /// partition rather than overlap), and the sibling is exempt from the qscore bar entirely,
    /// because its admission was earned by completing the intended set rather than by scoring well.
    const SpeciesSurveyRecord* seen = nullptr;

    bool   fan_out = false;          ///< precursor_charges != Single: the scan takes the whole envelope
    int    min_charge = 0;
    double snr_threshold = 0.0;      ///< the NOTCH gate; the anchor's own gate is @c anchor_snr_threshold
    double anchor_snr_threshold = 0.0;   ///< zeroed by the caller for a matched target (the SNR waiver)
    double qscore_threshold = 0.0;       ///< likewise zeroed for a matched target
    double tqscore_threshold = 0.0;
    int    max_notches = MAX_NOTCHES_PER_STAGE;
    int    min_charge_states = 1;    ///< floor on the ACQUIRED set size; 1 = no floor. Siblings exempt.
    std::string where;               ///< short label for selectNotches' clamp message
  };

  /// What the caller should do with this candidate.
  struct AdmissionVerdict
  {
    bool admit = false;

    /// Costs no `max_precursors` slot. Always false today; ADR-0036 turns it on for a PeakGroup that
    /// is completing a species another PeakGroup of the same species already started.
    bool is_sibling = false;

    /// Should this acquisition arm the two `tqscore_exceeding_*` bars? Separated from the admit
    /// decision because the arming side and the reading side are independent choices (ADR-0037).
    bool arms_bars = false;

    /// Did this candidate get past dynamic exclusion, whatever the qscore ledger then decided?
    ///
    /// Exposed because the loop records the species' acquisition-memory timestamp at that point,
    /// BEFORE the ledger can refuse — so a candidate the ledger turns away still refreshes the
    /// record that keeps its score alive through expiry. Collapsing this into @c admit would move
    /// that write, which is a change to when a mass's stored score expires rather than a
    /// simplification.
    bool passed_exclusion = false;

    /// The score to record against the species. Never lowers a bar the species already holds.
    double record_score = 0.0;

    /// Every charge this acquisition isolates: the anchor first, then its notches in the order
    /// `selectNotches` returned them. Size one under `single`.
    std::vector<int> acquisition_charges;

    AdmissionReason reason = AdmissionReason::NoAnchorResolved;
  };

  namespace CandidateAdmissionDetail
  {
    inline bool contains(const std::vector<int>& v, int x)
    {
      return std::find(v.begin(), v.end(), x) != v.end();
    }
  }

  /// The INTENDED CHARGE SET this PeakGroup can contribute (ADR-0036): what "acquired in full" is
  /// measured against, restricted to charges this PeakGroup actually resolved.
  ///
  /// **Only a row that NAMES charges has one.** The named charges it resolved, whatever the
  /// acquisition mode — a named charge is the method asking for that charge specifically, so
  /// `single` still owes them all and just pays one per survey.
  ///
  /// An UNRESTRICTED row has no intended set, and that is a scope decision rather than an oversight
  /// (ADR-0036, amended 2026-08-28). A `-1` row expresses no opinion about charge, so nothing about
  /// it is ever incomplete: the best-scoring mass is what the method asked for and acquiring it once
  /// satisfies the request. Returning `{}` here is what makes `completesSomething` false for every
  /// sibling of such a row, and therefore what keeps an unrestricted row byte-identical in `single`,
  /// `separate` and `multiplexed` alike. The acquisition mode is deliberately not consulted: it
  /// decides how many charges one SCAN isolates, not which charges the Precursor is owed.
  ///
  /// Charges below the level's floor are excluded: an intended set that named one would let a
  /// sibling be admitted to acquire something the floor forbids.
  inline std::vector<int> intendedChargeSet(const PeakGroup& pg, const std::vector<int>& authored,
                                            int min_charge)
  {
    std::vector<int> out;
    if (authored.empty()) return out;

    for (int c : authored)
    {
      if (min_charge > 0 && c < min_charge) continue;
      auto [lo, hi] = pg.getMzRange(c);
      if (hi > lo) out.push_back(c);
    }
    return out;
  }

  namespace CandidateAdmissionDetail
  {
    /// Is this candidate worth a scan given what the species already has (ADR-0036)?
    ///
    /// The FIRST PeakGroup of a species always is — there is nothing to compare against. A SIBLING
    /// earns one only by offering an intended-set member that no already-selected PeakGroup of the
    /// species RESOLVED. That predicate does three things at once, and all three are wanted:
    ///
    ///   - it completes a split envelope, which is the whole point;
    ///   - it refuses a sibling that merely repeats what the first PeakGroup could already see, so
    ///     the duplicate near-identical peak groups a targeted survey produces cost nothing;
    ///   - it preserves ADR-0028's deferral, because "resolved" is not "acquired": under `single`
    ///     the first PeakGroup resolves all of its named charges while acquiring one, so a sibling
    ///     offering the same charges adds nothing and the rest still arrive one survey at a time.
    ///
    /// An UNRESTRICTED row has an empty intended set, so the loop below never runs and every sibling
    /// of one is refused — which is exactly the pre-ADR-0036 behaviour, reached without a mode test.
    inline bool completesSomething(const PeakGroup& pg, const std::vector<int>& authored,
                                   int min_charge, const SpeciesSurveyRecord* seen)
    {
      if (seen == nullptr) { return true; }
      for (int c : intendedChargeSet(pg, authored, min_charge))
        if (!contains(seen->resolved, c)) { return true; }
      return false;
    }
  }

  /// Pick the anchor charge of one candidate.
  ///
  /// Under an AUTHORED charge set the anchor is drawn from the SET rather than from the envelope, so
  /// it has to be resolved before anything else looks at the PeakGroup. SNR ORDERS the named charges
  /// here; it does not gate them — a matched target's anchor is exempt from the SNR threshold, which
  /// is most of why one authors a target, and the notch side keeps the gate.
  ///
  /// With no authored set the anchor is the envelope's own: the best-qscore charge under
  /// `consider_all_charges`, otherwise the representative charge.
  inline AnchorChoice pickAnchor(const PeakGroup& pg, const AnchorContext& ctx)
  {
    AnchorChoice out;

    if (!ctx.authored.empty())
    {
      int   charge = -1;
      float best_snr = -1.0f;
      for (int c : ctx.authored)
      {
        // The set can only SUBTRACT, so it must not smuggle a charge past the level floor.
        if (ctx.min_charge > 0 && c < ctx.min_charge) continue;
        auto [lo, hi] = pg.getMzRange(c);
        if (hi <= lo) continue;  // never resolved -- there is no MEASURED window to isolate
        if (CandidateAdmissionDetail::contains(ctx.spent_charges, c)) continue;  // already spent
        const float snr = pg.getChargeSNR(c);
        if (snr > best_snr) { best_snr = snr; charge = c; }
      }
      // Every named charge spent, absent, or below the floor: the species is finished. There is
      // deliberately no fallback to an unnamed charge -- that would be the set ADDING.
      if (charge < 0) { out.reason = AdmissionReason::NoAnchorResolved; return out; }
      if (!CandidateAdmissionDetail::completesSomething(pg, ctx.authored, ctx.min_charge, ctx.seen))
      {
        out.reason = AdmissionReason::IntendedSetComplete;
        return out;
      }

      // The anchor's OWN qscore, not a sibling's. The pair used to be bound together and then the
      // charge alone reassigned by inclusion matching, so a scan could be logged and excluded on the
      // score of a charge it never isolated. getAllQscores() is empty for a PeakGroup that never went
      // through updateQscore, and the fallback below is the path a synthetic one takes.
      const std::unordered_map<int, float> per_charge_qscores = pg.getAllQscores();
      auto qs_it = per_charge_qscores.find(charge);
      out.charge = charge;
      out.score  = qs_it != per_charge_qscores.end() ? (double)qs_it->second : pg.getQscore();
      out.reason = AdmissionReason::Admitted;
      return out;
    }

    // No authored set, so no intended charge set and nothing a sibling could complete: this refuses
    // every sibling of an unrestricted row, in every acquisition mode. It is the pre-ADR-0036 guard
    // ("this species was already acquired this survey") reached through the same predicate the
    // authored path uses, rather than through a second rule that could drift from it.
    if (!CandidateAdmissionDetail::completesSomething(pg, ctx.authored, ctx.min_charge, ctx.seen))
    {
      out.reason = AdmissionReason::IntendedSetComplete;
      return out;
    }

    if (ctx.consider_all_charges)
    {
      out.charge = pg.getBestQScoreCharge();
      out.score  = pg.getBestQScore();
    }
    else
    {
      out.charge = pg.getRepAbsCharge();
      out.score  = pg.getQscore();
    }
    out.reason = AdmissionReason::Admitted;
    return out;
  }

  /// The ACQUISITION CHARGE SET of one candidate: every charge this scan will isolate.
  ///
  /// Both non-single modes take the SAME set and differ only in scan count (ADR-0016) — multiplexed
  /// co-isolates it as notches in one scan, separate emits one scan per member — so deriving both
  /// from this one call is what keeps that true rather than aspirational. It is the same
  /// `peakGroupNotchCandidates` + `selectNotches` pair `buildMS2` uses, which is what makes the set
  /// recorded as acquired identical by construction to the set actually isolated.
  ///
  /// An authored charge set FILTERS the candidates and never extends them; `selectNotches` still
  /// applies the SNR gate afterwards, so a named-but-weak charge is refused as a notch even though
  /// the anchor is exempt.
  inline std::vector<int> acquisitionChargeSet(const PeakGroup& pg, int anchor_charge,
                                               const AdmissionContext& ctx)
  {
    std::vector<int> out {anchor_charge};
    if (!ctx.fan_out) return out;

    std::vector<NotchCandidate> notch_candidates = peakGroupNotchCandidates(pg, optimal_window_margin_);
    if (ctx.has_authored)
    {
      notch_candidates.erase(
          std::remove_if(notch_candidates.begin(), notch_candidates.end(),
                         [&ctx](const NotchCandidate& n) {
                           return !CandidateAdmissionDetail::contains(ctx.authored, n.charge)
                                  || CandidateAdmissionDetail::contains(ctx.spent_charges, n.charge)
                                  || (ctx.min_charge > 0 && n.charge < ctx.min_charge);
                         }),
          notch_candidates.end());
    }
    for (const NotchCandidate& n :
         selectNotches(notch_candidates, anchor_charge, ctx.snr_threshold, ctx.max_notches, ctx.where))
      out.push_back(n.charge);

    // A SIBLING isolates only what the species has not already isolated this survey (ADR-0036), so
    // the scans PARTITION the envelope instead of overlapping it. Re-isolating a charge a previous
    // scan of the same survey already took splits its ions across two fills and makes both weaker
    // than the single isolation would have been -- the same reasoning ADR-0016's SNR gate rests on.
    // The anchor is never dropped: it was chosen from charges that are still owed.
    if (ctx.seen != nullptr)
    {
      out.erase(std::remove_if(out.begin() + 1, out.end(),
                               [&ctx](int c) {
                                 return CandidateAdmissionDetail::contains(ctx.seen->acquired, c);
                               }),
                out.end());
    }
    return out;
  }

  /// The admission decision proper, past the anchor and past the caller's target matching.
  ///
  /// Gate order is the loop's own and is not incidental: the score gate precedes the SNR gate, both
  /// precede dynamic exclusion, and the qscore ledger is consulted last — after the acquisition
  /// charge set has been computed, because an authored species records per charge and skips the
  /// mass-keyed ledger entirely.
  inline AdmissionVerdict admitCandidate(const PeakGroup& pg, int anchor_charge, double score,
                                         const AdmissionContext& ctx)
  {
    AdmissionVerdict v;

    if (ctx.min_charge > 0 && anchor_charge < ctx.min_charge)
    {
      v.reason = AdmissionReason::BelowMinCharge;
      return v;
    }
    if (score < ctx.qscore_threshold)
    {
      v.reason = AdmissionReason::ScoreBelowThreshold;
      return v;
    }
    if (pg.getChargeSNR(anchor_charge) < ctx.anchor_snr_threshold)
    {
      v.reason = AdmissionReason::SnrBelowThreshold;
      return v;
    }

    v.is_sibling = ctx.seen != nullptr;

    // Dynamic exclusion. An authored species is keyed per (mass, charge) so `single` can come back
    // for the charges the row named; everything else is keyed on the mass and on the integer m/z of
    // the anchor's isolation centre, ORed.
    //
    // A species that matched an inclusion row does not READ either mass-keyed bar (ADR-0037). It
    // still arms them below, once, so every other species sees exactly what it sees today; what
    // changes is that a target is no longer barred from its own re-acquisition for the whole
    // retention-time window regardless of how much better a later survey resolves it. With
    // tqscore_threshold at 0.1 against observed qscores of 0.48-0.98, the first acquisition armed
    // those bars every time and the qscore ledger below was unreachable.
    if (ctx.has_authored)
    {
      if (ctx.anchor_spent)
      {
        v.reason = AdmissionReason::Barred;
        return v;
      }
    }
    else if (ctx.mass_barred || ctx.mz_barred)
    {
      v.reason = AdmissionReason::Barred;
      return v;
    }

    // Past exclusion. The caller stamps acquisition memory here, before the ledger below can refuse.
    v.passed_exclusion = true;

    v.acquisition_charges = acquisitionChargeSet(pg, anchor_charge, ctx);

    // ADR-0040. A species whose set is too small to vote on is not worth the scans, and refusing
    // HERE -- rather than at the consumer -- is what frees its max_precursors slot for one that is.
    //
    // The sibling exemption is load-bearing, not defensive: a sibling's set is REDUCED by what the
    // species already isolated (ADR-0036), so the PeakGroup supplying the last charge of a split
    // envelope presents a set of size 1. Without this clause the floor refuses precisely the charge
    // ADR-0036 exists to recover, and split-envelope completion regresses silently.
    //
    // Placed after passed_exclusion deliberately: all_mass_rt_map_ has no readers, it only carries
    // mass_qscore_map_ entries through expiry, and a species never acquired has no qscore entry --
    // so operator[] default-inserts a 0.0 bar, which bars nothing. Moving the gate earlier would
    // reorder exclusion stamping, which is a behaviour change dressed as a refactor.
    if (ctx.min_charge_states > 1 && !v.is_sibling
        && static_cast<int>(v.acquisition_charges.size()) < ctx.min_charge_states)
    {
      v.reason = AdmissionReason::BelowMinChargeStates;
      return v;   // v.admit stays false -> no slot spent
    }

    // The qscore ledger is mass-keyed and therefore skipped entirely for an authored species: it
    // would retire the MASS after its first charge and leave the other named charges unreachable,
    // and its "previously acquired at a higher score" guard would refuse the second named charge on
    // the next survey, since charges are ranked descending and every later one scores lower.
    if (!ctx.has_authored)
    {
      // "Previously acquired at a higher qscore" -- the ratchet, and it is older than any of this
      // (`4457abeb4c change algorithm to max qscore`). Whether it governs a species at all is a
      // CONFIGURATION question, not a rule to be special-cased here: the arming test three lines
      // below decides it, keyed on `tqscore_threshold`. Above a species' qscore the bars never arm
      // and this comparison runs every survey; below it, the species is barred outright further up
      // and never reaches here. ADR-0037 proposed to override that for matched inclusion targets
      // and was WITHDRAWN -- it removed the choice rather than adding a capability, and bought
      // re-acquisitions on qscore improvements of +0.0002.
      //
      // No sibling exemption: siblings exist only for authored rows now (ADR-0036, amended
      // 2026-08-28), and an authored row never reaches this branch at all.
      if (ctx.qscore_bar != nullptr && score < *ctx.qscore_bar)
      {
        v.reason = AdmissionReason::ScoreNotBetter;
        return v;
      }

      v.record_score = score;
      v.arms_bars    = v.record_score > ctx.tqscore_threshold;
    }

    v.admit  = true;
    v.reason = AdmissionReason::Admitted;
    return v;
  }
}  // namespace OpenMS
