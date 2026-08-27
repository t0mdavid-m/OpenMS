// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/CandidateAdmission.h>

#include <string>
#include <vector>

using namespace OpenMS;

// The per-candidate acquisition decision, driven directly with synthetic PeakGroups.
//
// This suite exists because the behaviour it guards is UNREACHABLE from a spectrum fixture. The
// defect ADR-0036 addresses needs one species' charge envelope to arrive as several PeakGroups, and
// no committed fixture does that: ms1_cytc.txt resolves cytC at 13 charges in ONE PeakGroup (see
// FLASHIda_ChargeModes_test's CM-05), and the inclusion goldens' duplicate peak groups at 12351.30x
// are weak leftovers rather than envelope halves. An end-to-end golden therefore cannot see either
// the defect or its fix, which is exactly why the decision was lifted out of the selection loop.
//
// THIS FILE PINS THE EXTRACTION, NOT THE NEW BEHAVIOUR. Every expectation below is the selection
// loop's behaviour as it stands, so a green run here plus byte-identical log goldens is the proof
// that lifting the decision out changed nothing. The cases that pin ADR-0036/0037 -- sibling
// admission, the partitioned acquisition set, the budget exemption, the bar exemption for matched
// targets, and the selector/buildMS2 drift guard -- arrive with the behaviour change itself, where
// their expectations are written once rather than flipped.

namespace
{
  const double PROTON = 1.00727646677;
  const double ISO_GAP = 1.00235;

  /// A synthetic MS1 PeakGroup resolving exactly @p charges.
  ///
  /// Two details are load-bearing and both are silent when wrong:
  ///
  /// - THREE peaks per charge, not one. getMzRange returns {mz, mz} for a single peak, and every
  ///   consumer reads hi <= lo as "never resolved" -- so a one-peak-per-charge helper would make
  ///   every assertion below pass vacuously.
  /// - The constructor's charge range bounds the per-charge SNR vector AND getMzRange, which returns
  ///   {-1, -10} outside it. That is what lets one PeakGroup genuinely not resolve another's
  ///   charges, which is the whole scenario under test.
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

  AnchorContext anchorCtx(std::vector<int> authored = {}, std::vector<int> spent = {})
  {
    AnchorContext c;
    c.authored = std::move(authored);
    c.spent_charges = std::move(spent);
    return c;
  }

  /// Defaults that admit: no bars, no floor, thresholds at zero.
  AdmissionContext admitCtx(std::vector<int> authored = {}, std::vector<int> spent = {})
  {
    AdmissionContext c;
    c.has_authored = !authored.empty();
    c.authored = std::move(authored);
    c.spent_charges = std::move(spent);
    c.where = "test";
    return c;
  }
}

START_TEST(FLASHIda_CandidateAdmission, "$Id$")

/////////////////////////////////////////////////////////////

// CA-B1: with no authored charge set the anchor is the envelope's own representative charge, and the
// score is the PeakGroup's scalar qscore. This is the path every non-targeted species takes.
START_SECTION(anchor_without_an_authored_set_is_the_representative_charge)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.94, 10.0f);
  pg.setRepAbsCharge(13);

  AnchorChoice a = pickAnchor(pg, anchorCtx());

  TEST_EQUAL(a.resolved(), true)
  TEST_EQUAL(a.charge, 13)
  TEST_REAL_SIMILAR(a.score, 0.94)
}
END_SECTION

// CA-B2: consider_all_charges reads the PER-CHARGE qscore map, which only updateQscore populates.
// A PeakGroup that never went through it has no per-charge scores at all, so this ranks over an
// empty map and yields no anchor. Pinned rather than worked around: it is the honest boundary of
// what a synthetic PeakGroup can express, and a future helper that starts populating qscores_ should
// have to notice this expectation change rather than silently strengthen four other tests.
START_SECTION(consider_all_charges_ranks_over_the_per_charge_qscore_map)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.94, 10.0f);
  AnchorContext ctx = anchorCtx();
  ctx.consider_all_charges = true;

  AnchorChoice a = pickAnchor(pg, ctx);

  TEST_EQUAL(pg.getAllQscores().empty(), true)   // the premise, stated rather than assumed
  TEST_EQUAL(a.resolved(), false)
}
END_SECTION

// CA-B3: under an AUTHORED charge set the anchor is drawn from the SET, by SNR, and overrides the
// representative charge. SNR orders here and does not gate -- a matched target's anchor is exempt
// from the SNR threshold, which is most of why one authors a target.
START_SECTION(authored_set_moves_the_anchor_to_its_highest_snr_member)
{
  PeakGroup pg = makeEnvelope(12351.39, {10, 11, 12, 13, 14}, 0.94, 4.0f);
  pg.setRepAbsCharge(12);        // the envelope's own pick, which the authored set must override
  pg.setChargeSNR(13, 40.0f);    // the strongest NAMED charge
  pg.setChargeSNR(11, 90.0f);    // stronger still, but not named -- must not win

  AnchorChoice a = pickAnchor(pg, anchorCtx({10, 13}));

  TEST_EQUAL(a.charge, 13)
  TEST_REAL_SIMILAR(a.score, 0.94)   // getAllQscores() is empty, so the anchor falls back to getQscore()
}
END_SECTION

// CA-B4: a named charge the survey never resolved acquires nothing. An isolation window has to be
// MEASURED, so there is nothing to isolate and the set cannot conjure it.
START_SECTION(a_named_charge_that_was_never_resolved_is_skipped)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.94, 10.0f);

  AnchorChoice named_only_outside = pickAnchor(pg, anchorCtx({16, 17}));
  TEST_EQUAL(named_only_outside.resolved(), false)
  TEST_EQUAL(named_only_outside.reason == AdmissionReason::NoAnchorResolved, true)

  AnchorChoice one_inside = pickAnchor(pg, anchorCtx({13, 16, 17}));
  TEST_EQUAL(one_inside.charge, 13)
}
END_SECTION

// CA-B5: when every named charge has been spent inside the retention-time window the species is
// finished. There is deliberately no fallback to an unnamed charge -- that would be the set ADDING.
START_SECTION(a_fully_spent_authored_set_yields_no_anchor)
{
  PeakGroup pg = makeEnvelope(12351.39, {10, 13, 16}, 0.94, 10.0f);

  AnchorChoice a = pickAnchor(pg, anchorCtx({10, 13, 16}, {10, 13, 16}));

  TEST_EQUAL(a.resolved(), false)
  TEST_EQUAL(a.reason == AdmissionReason::NoAnchorResolved, true)

  AnchorChoice partial = pickAnchor(pg, anchorCtx({10, 13, 16}, {10, 16}));
  TEST_EQUAL(partial.charge, 13)   // the one still owed
}
END_SECTION

// CA-B6: the authored set restricts and must not smuggle a charge past the level's floor.
START_SECTION(an_authored_charge_below_the_floor_is_not_smuggled_in)
{
  PeakGroup pg = makeEnvelope(12351.39, {8, 9, 10, 13}, 0.94, 10.0f);
  pg.setChargeSNR(8, 99.0f);   // the strongest named charge, but below the floor

  AnchorContext ctx = anchorCtx({8, 13});
  ctx.min_charge = 10;

  TEST_EQUAL(pickAnchor(pg, ctx).charge, 13)
}
END_SECTION

// CA-B7: a SIBLING is judged on what it can ADD, not refused for existing (ADR-0036).
//
// The reported defect, reduced to its two peak groups: a row naming 12..18, one PeakGroup resolving
// z12-14 and another z15-18 within 0.1 Da. Before this rule the second was refused outright and
// z15-18 were deferred until the first PeakGroup's charges were spent -- three surveys under
// `single` -- and lost entirely if the peak stopped eluting first.
START_SECTION(a_sibling_that_completes_the_intended_set_is_admitted)
{
  const std::vector<int> named {12, 13, 14, 15, 16, 17, 18};
  PeakGroup pg_b = makeEnvelope(12351.33, {15, 16, 17, 18}, 0.72, 10.0f);

  SpeciesSurveyRecord seen;             // what PG_A, resolving z12-14, left behind
  seen.resolved = {12, 13, 14};
  seen.acquired = {12, 13, 14};
  seen.best_score = 0.94;

  AnchorContext ctx = anchorCtx(named);
  ctx.seen = &seen;
  ctx.fan_out = true;

  AnchorChoice a = pickAnchor(pg_b, ctx);
  TEST_EQUAL(a.resolved(), true)
  TEST_EQUAL(a.reason == AdmissionReason::Admitted, true)
  TEST_EQUAL(a.charge >= 15 && a.charge <= 18, true)   // drawn from what only PG_B could resolve
}
END_SECTION

// CA-B7b: the other half of the same rule, and the one that keeps it from becoming "admit
// everything". A sibling resolving only charges an earlier PeakGroup ALSO resolved adds nothing and
// is refused -- which is what makes the duplicate near-identical peak groups a targeted survey
// produces cost nothing, and what preserves ADR-0028's one-named-charge-per-survey deferral.
START_SECTION(a_sibling_that_adds_nothing_is_refused)
{
  const std::vector<int> named {12, 13, 14, 15, 16, 17, 18};
  PeakGroup dup = makeEnvelope(12351.40, {12, 13, 14}, 0.71, 10.0f);

  SpeciesSurveyRecord seen;
  seen.resolved = {12, 13, 14};
  seen.acquired = {13};               // `single`: only ONE was taken, but all three were RESOLVED
  seen.best_score = 0.94;

  AnchorContext ctx = anchorCtx(named);
  ctx.seen = &seen;

  AnchorChoice a = pickAnchor(dup, ctx);
  TEST_EQUAL(a.resolved(), false)
  TEST_EQUAL(a.reason == AdmissionReason::IntendedSetComplete, true)
}
END_SECTION

// CA-B7c: an UNRESTRICTED row splits the same way, and the rule has to reach it -- the defect is
// identical there, only the guard that used to kill the sibling differs. Under `single` the intended
// set is the anchor alone, so no sibling is ever admitted and "the best-qscore mass is enough"
// stands; under multiplexed/separate the intended set is the signal-bearing envelope and a sibling
// completing it is admitted.
START_SECTION(an_unrestricted_row_admits_siblings_only_when_it_asked_for_several_charges)
{
  PeakGroup pg_b = makeEnvelope(12351.33, {15, 16, 17, 18}, 0.72, 10.0f);

  SpeciesSurveyRecord seen;
  seen.resolved = {12, 13, 14};
  seen.acquired = {12, 13, 14};

  AnchorContext fan = anchorCtx();      // authored = {} -- an unrestricted row
  fan.seen = &seen;
  fan.fan_out = true;
  TEST_EQUAL(pickAnchor(pg_b, fan).reason == AdmissionReason::Admitted, true)

  AnchorContext single = anchorCtx();
  single.seen = &seen;
  single.fan_out = false;
  TEST_EQUAL(pickAnchor(pg_b, single).reason == AdmissionReason::IntendedSetComplete, true)
}
END_SECTION

// CA-B8: the two per-candidate gates, and the order they fire in. Score precedes SNR; both precede
// dynamic exclusion. The order matters for the reason a refusal reports, which is the only thing
// telling an operator which gate to loosen.
START_SECTION(score_and_snr_gates_fire_in_order)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.50, 2.0f);

  AdmissionContext low_score = admitCtx();
  low_score.qscore_threshold = 0.80;
  low_score.anchor_snr_threshold = 100.0;   // would also fail, but must not be the reported reason
  TEST_EQUAL(admitCandidate(pg, 13, 0.50, low_score).reason == AdmissionReason::ScoreBelowThreshold, true)

  AdmissionContext low_snr = admitCtx();
  low_snr.anchor_snr_threshold = 100.0;
  AdmissionVerdict v = admitCandidate(pg, 13, 0.50, low_snr);
  TEST_EQUAL(v.admit, false)
  TEST_EQUAL(v.reason == AdmissionReason::SnrBelowThreshold, true)

  TEST_EQUAL(admitCandidate(pg, 13, 0.50, admitCtx()).admit, true)   // neither gate armed
}
END_SECTION

// CA-B9: dynamic exclusion is ORed over the nominal-mass bar and the integer-m/z bar. Either alone
// refuses -- which is why exempting only one of them changes nothing whenever the anchor charge
// repeats, and the reason ADR-0037 had to name both.
START_SECTION(either_exclusion_bar_alone_refuses_an_unauthored_species)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.94, 10.0f);

  AdmissionContext mass_only = admitCtx();
  mass_only.mass_barred = true;
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, mass_only).reason == AdmissionReason::Barred, true)

  AdmissionContext mz_only = admitCtx();
  mz_only.mz_barred = true;
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, mz_only).reason == AdmissionReason::Barred, true)
}
END_SECTION

// CA-B10: an authored species is excluded per (mass, charge) and ignores the mass-keyed bars
// entirely -- otherwise `single` could never come back for the second and third named charge.
START_SECTION(an_authored_species_is_excluded_per_charge_not_by_the_mass_bars)
{
  PeakGroup pg = makeEnvelope(12351.39, {10, 13, 16}, 0.94, 10.0f);

  AdmissionContext barred = admitCtx({10, 13, 16});
  barred.mass_barred = true;
  barred.mz_barred = true;
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, barred).admit, true)   // the mass bars do not reach it

  AdmissionContext spent = admitCtx({10, 13, 16});
  spent.anchor_spent = true;
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, spent).reason == AdmissionReason::Barred, true)

  // anchor_spent and spent_charges are NOT interchangeable, and this pins why. The loop consults the
  // per-charge map unconditionally for the anchor pick and the notch filter, but only in the
  // exclusion-applying selection phases for the bar. So a populated spent_charges must NOT bar on
  // its own -- deriving the bar from it would apply exclusion in a phase that does not.
  AdmissionContext notch_filter_only = admitCtx({10, 13, 16}, {13});
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, notch_filter_only).admit, true)
}
END_SECTION

// CA-B11: the qscore ledger. A species already acquired at a higher score is refused; one this
// survey resolves at least as well is re-acquired and raises the bar.
START_SECTION(the_qscore_bar_admits_only_an_equal_or_better_survey)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.80, 10.0f);

  double bar = 0.90;
  AdmissionContext worse = admitCtx();
  worse.qscore_bar = &bar;
  TEST_EQUAL(admitCandidate(pg, 13, 0.80, worse).reason == AdmissionReason::ScoreNotBetter, true)

  AdmissionContext better = admitCtx();
  better.qscore_bar = &bar;
  AdmissionVerdict v = admitCandidate(pg, 13, 0.95, better);
  TEST_EQUAL(v.admit, true)
  TEST_REAL_SIMILAR(v.record_score, 0.95)

  AdmissionContext equal = admitCtx();
  equal.qscore_bar = &bar;
  TEST_EQUAL(admitCandidate(pg, 13, 0.90, equal).admit, true)   // `<` is strict: equal re-acquires

  // A ledger refusal still PASSED exclusion, and the caller has to be able to tell the two apart:
  // the species' acquisition-memory timestamp is stamped once exclusion is cleared, before the
  // ledger runs, so a candidate the ledger turns away still refreshes the record that keeps its
  // stored score alive through expiry. A refusal at the bars stamps nothing.
  AdmissionVerdict ledger_refusal = admitCandidate(pg, 13, 0.80, worse);
  TEST_EQUAL(ledger_refusal.admit, false)
  TEST_EQUAL(ledger_refusal.passed_exclusion, true)

  AdmissionContext bar_refusal = admitCtx();
  bar_refusal.mass_barred = true;
  TEST_EQUAL(admitCandidate(pg, 13, 0.80, bar_refusal).passed_exclusion, false)
}
END_SECTION

// CA-B12: arming the exclusion bars is a separate decision from admitting, and it keys on the score
// the species is recorded at rather than on the raw score. An authored species arms nothing.
START_SECTION(the_exclusion_bars_are_armed_by_the_recorded_score)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.94, 10.0f);

  AdmissionContext above = admitCtx();
  above.tqscore_threshold = 0.10;
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, above).arms_bars, true)

  AdmissionContext below = admitCtx();
  below.tqscore_threshold = 0.99;
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, below).arms_bars, false)

  AdmissionContext authored = admitCtx({12, 13, 14});
  authored.tqscore_threshold = 0.10;
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, authored).arms_bars, false)
}
END_SECTION

// CA-B13: the acquisition charge set. Size one under `single`; the whole SNR-positive envelope when
// the mode asks for several charges, anchor first.
START_SECTION(the_acquisition_charge_set_is_the_anchor_or_the_envelope)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.94, 10.0f);

  AdmissionContext single = admitCtx();
  std::vector<int> one = acquisitionChargeSet(pg, 13, single);
  TEST_EQUAL(one.size(), 1)
  TEST_EQUAL(one[0], 13)

  AdmissionContext fan = admitCtx();
  fan.fan_out = true;
  std::vector<int> many = acquisitionChargeSet(pg, 13, fan);
  TEST_EQUAL(many.size(), 3)
  TEST_EQUAL(many[0], 13)   // the anchor leads; selectNotches drops it from its own output
}
END_SECTION

// CA-B14: an authored set FILTERS the notch candidates and never extends them, and a spent charge is
// filtered out too -- so a co-isolating scan cannot re-isolate a charge an earlier survey took.
START_SECTION(an_authored_set_filters_the_notches_and_never_extends_them)
{
  PeakGroup pg = makeEnvelope(12351.39, {10, 11, 12, 13, 14, 15, 16}, 0.94, 10.0f);

  AdmissionContext fan = admitCtx({10, 13, 16});
  fan.fan_out = true;
  std::vector<int> got = acquisitionChargeSet(pg, 13, fan);
  TEST_EQUAL(got.size(), 3)
  for (int z : got) TEST_EQUAL(z == 10 || z == 13 || z == 16, true)

  AdmissionContext with_spent = admitCtx({10, 13, 16}, {10});
  with_spent.fan_out = true;
  std::vector<int> fewer = acquisitionChargeSet(pg, 13, with_spent);
  TEST_EQUAL(fewer.size(), 2)
  for (int z : fewer) TEST_EQUAL(z == 13 || z == 16, true)
}
END_SECTION

// CA-B15: the two scans PARTITION the envelope (ADR-0036). A sibling drops charges the species has
// already isolated this survey, so its ions are not split across two fills -- the same reasoning
// ADR-0016's SNR gate rests on. The anchor is never dropped: it was picked from charges still owed.
START_SECTION(a_siblings_scan_is_reduced_by_what_the_species_already_isolated)
{
  PeakGroup pg_b = makeEnvelope(12351.33, {13, 15, 16, 17, 18}, 0.72, 10.0f);

  SpeciesSurveyRecord seen;
  seen.resolved = {12, 13, 14};
  seen.acquired = {12, 13, 14};        // z13 already isolated by the earlier PeakGroup

  AdmissionContext ctx = admitCtx();
  ctx.fan_out = true;
  ctx.seen = &seen;

  std::vector<int> got = acquisitionChargeSet(pg_b, 16, ctx);
  TEST_EQUAL(got[0], 16)                                                   // the anchor leads
  for (int z : got) TEST_EQUAL(z == 13, false)                             // z13 not re-isolated
  TEST_EQUAL(got.size(), 4)                                                // 16 + 15,17,18
}
END_SECTION

// CA-B16: an admitted sibling costs no max_precursors slot (ADR-0036). The predicate that admitted
// it just declared it the same species as the PeakGroup that already paid, so charging again would
// contradict the admission in the same breath -- and would let one split target eat a whole survey.
START_SECTION(an_admitted_sibling_is_not_charged_a_budget_slot)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.94, 10.0f);

  TEST_EQUAL(admitCandidate(pg, 13, 0.94, admitCtx()).is_sibling, false)   // first of its species

  SpeciesSurveyRecord seen;
  seen.resolved = {15, 16};
  AdmissionContext sib = admitCtx();
  sib.seen = &seen;
  AdmissionVerdict v = admitCandidate(pg, 13, 0.94, sib);
  TEST_EQUAL(v.admit, true)
  TEST_EQUAL(v.is_sibling, true)
}
END_SECTION

// CA-B17: a matched inclusion target is not barred by dynamic exclusion (ADR-0037). Without this the
// qscore ledger below is unreachable: production sets tqscore_threshold to 0.1 against observed
// qscores of 0.48-0.98, so the first acquisition arms both bars every time.
START_SECTION(a_matched_target_does_not_read_the_exclusion_bars)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.94, 10.0f);

  AdmissionContext target = admitCtx();
  target.target_matched = true;
  target.mass_barred = true;
  target.mz_barred = true;
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, target).admit, true)

  // ...and the exemption is scoped to targets. An unmatched species stays barred exactly as today,
  // which is what keeps every non-inclusion run byte-identical.
  AdmissionContext other = admitCtx();
  other.mass_barred = true;
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, other).reason == AdmissionReason::Barred, true)
}
END_SECTION

// CA-B18: completion beats the score bar outright (ADR-0037, M2-wide). This is not decoration:
// siblings rank below the first PeakGroup by construction, so a score test would turn every one of
// them away and a split envelope could never be completed however permissive the other guards were.
START_SECTION(an_admitted_sibling_is_exempt_from_the_score_bar)
{
  PeakGroup pg = makeEnvelope(12351.33, {15, 16, 17, 18}, 0.55, 10.0f);
  double bar = 0.60;                    // an EARLIER survey's bar, higher than this sibling scores

  SpeciesSurveyRecord seen;
  seen.resolved = {12, 13, 14};
  seen.best_score = 0.72;

  AdmissionContext sib = admitCtx();
  sib.qscore_bar = &bar;
  sib.seen = &seen;
  TEST_EQUAL(admitCandidate(pg, 16, 0.55, sib).admit, true)

  // A NON-sibling scoring the same is still refused -- the exemption is about completion, not about
  // being lenient.
  AdmissionContext first = admitCtx();
  first.qscore_bar = &bar;
  TEST_EQUAL(admitCandidate(pg, 16, 0.55, first).reason == AdmissionReason::ScoreNotBetter, true)
}
END_SECTION

// CA-B19: the bar only ever RISES (W-max). A sibling admitted on charge grounds may score well below
// what the species was already acquired at; letting it write its own score would lower the bar the
// next survey faces and turn the ratchet into a random walk.
START_SECTION(a_low_scoring_sibling_cannot_lower_the_bar)
{
  PeakGroup pg = makeEnvelope(12351.33, {15, 16, 17, 18}, 0.55, 10.0f);
  double bar = 0.60;

  SpeciesSurveyRecord seen;
  seen.resolved = {12, 13, 14};
  seen.best_score = 0.72;

  AdmissionContext sib = admitCtx();
  sib.qscore_bar = &bar;
  sib.seen = &seen;

  AdmissionVerdict v = admitCandidate(pg, 16, 0.55, sib);
  TEST_EQUAL(v.admit, true)
  TEST_REAL_SIMILAR(v.record_score, 0.72)   // the survey's best, not this sibling's 0.55
}
END_SECTION

// CA-B20: a matched target arms the bars ONCE and never refreshes them (ADR-0037) -- refreshing
// would suppress its bin and window neighbours for its whole elution, longer than today.
//
// The second half is the one that matters for "non-target behaviour is unchanged": a species that is
// NOT a target still arms on a later, better acquisition. Its first acquisition scoring below
// tqscore_threshold never armed the bars, so nothing ever blocked it from coming back -- and gating
// arming on "acquired before" would leave such a mass permanently unbarred.
START_SECTION(only_a_matched_target_stops_re_arming_the_bars)
{
  PeakGroup pg = makeEnvelope(12351.39, {12, 13, 14}, 0.94, 10.0f);
  double bar = 0.50;

  AdmissionContext first_target = admitCtx();
  first_target.target_matched = true;
  first_target.tqscore_threshold = 0.10;
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, first_target).arms_bars, true)

  AdmissionContext again_target = admitCtx();
  again_target.target_matched = true;
  again_target.tqscore_threshold = 0.10;
  again_target.qscore_bar = &bar;                       // acquired earlier in this rt window
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, again_target).arms_bars, false)

  AdmissionContext again_other = admitCtx();
  again_other.tqscore_threshold = 0.10;
  again_other.qscore_bar = &bar;                        // same shape, but NOT a target
  TEST_EQUAL(admitCandidate(pg, 13, 0.94, again_other).arms_bars, true)
}
END_SECTION

/////////////////////////////////////////////////////////////

END_TEST
