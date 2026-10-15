// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Quantification.h>

#include <string>
#include <vector>

using namespace OpenMS;

///////////////////////////////////////////////////////////////////////////////
//
// Quantification::measure -- the VERDICT logic (ADR-0038).
//
// This file exists because the pre-ADR-0038 suite would have passed with the whole differential
// test replaced by `return true;`. Both C++ fixtures set fold_change_threshold to 0.01, which makes
// `(fc > t) || (1/fc > t)` unconditionally true, and there was no negative case anywhere. Driving
// the full engine cannot fix that cheaply either: the interesting cases need spectra with
// SPECIFIC channels empty, which is awkward to arrange through an acquisition and trivial here.
//
// measure() is a pure function of (spectrum, config), so every case below is a synthesized TMT6
// reporter region plus a config, and no acquisition is involved.
//
///////////////////////////////////////////////////////////////////////////////

namespace
{
  // TMT6 reporter centres, from TMTSixPlexQuantitationMethod. Channel order == the ordinals
  // conditions resolve to, which is the order `intensities` below is given in.
  const double kTmt6[6] = {126.127725, 127.124760, 128.134433, 129.131468, 130.141141, 131.138176};

  /// A spectrum carrying exactly the six reporter peaks at the given intensities, plus one
  /// non-reporter peak so the spectrum is not degenerate. An intensity of 0 omits the peak
  /// entirely -- which is what an absent channel looks like to the extractor.
  struct Spec
  {
    std::vector<double> mzs, ints;
    explicit Spec(const std::vector<double>& intensities)
    {
      for (size_t i = 0; i < 6; ++i)
        if (intensities[i] > 0) { mzs.push_back(kTmt6[i]); ints.push_back(intensities[i]); }
      mzs.push_back(700.5); ints.push_back(1234.0);   // an ordinary fragment ion
    }
  };

  /// A loadable config with quantification on. `conditions` is spliced verbatim so each section can
  /// choose its own grouping; `extra` appends further quantification keys.
  std::string quantJson(const std::string& conditions, const std::string& extra = "")
  {
    return std::string(R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": "qscore", "max_precursors": 5 },
      "characterization": { "mode": "off" },
      "quantification": { "enabled": true, "labelling": "tmt6plex", "reporter_mz_tol": 0.002,
                          "conditions": )") + conditions + extra + R"( },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
        "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 },
        "ms2_quant": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 30 }
      }
    })";
  }

  // The default 3-vs-3 grouping, i.e. what the retired positional split used to hardcode.
  const char* kThreeVsThree =
      R"([ { "name": "a", "channels": ["126","127","128"] },
            { "name": "b", "channels": ["129","130","131"] } ])";

  // An all-zero matrix is the IDENTITY -- IsobaricQuantitationMethod skips "0.0"/"-1"/"NA" and
  // leaves self_contribution at 100 -- so it is the documented way to turn correction off.
  //
  // Every section that asserts an exact mean or ratio uses it, and that is not a convenience: with
  // the stock matrix the reported channels are CORRECTED, so they are deliberately not the raw
  // intensities fed in, and asserting raw numbers against them would be asserting that correction
  // does nothing. The one section that wants correction on is the one that tests correction.
  const char* kNoCorrection =
      R"(, "correction_matrix": ["0.0/0.0/0.0/0.0","0.0/0.0/0.0/0.0","0.0/0.0/0.0/0.0",
                                 "0.0/0.0/0.0/0.0","0.0/0.0/0.0/0.0","0.0/0.0/0.0/0.0"])";

  /// threshold 1.4, correction off -- the arithmetic default for this file.
  std::string thresh14() { return std::string(R"(, "fold_change_threshold": 1.4)") + kNoCorrection; }
  /// threshold 1.4 with correction left at the method's STOCK matrix.
  std::string thresh14Corrected() { return R"(, "fold_change_threshold": 1.4)"; }

  Quantification::Result run(const Config& cfg, const std::vector<double>& intensities)
  {
    Quantification q(cfg);
    Spec s(intensities);
    return q.measure(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), 1.0, 2, "test");
  }
}

START_TEST(Quantification, "$Id$")

// ---------------------------------------------------------------------------------------------
START_SECTION(differential_and_not_differential_are_both_reachable)
{
  // THE test the old suite lacked. A threshold of 1.4 with a 4:1 ratio is differential; the same
  // threshold with a 1:1 ratio is not. Replace measure()'s body with `return true;` and the second
  // half fails -- which is exactly what nothing asserted before ADR-0038.
  Config cfg(quantJson(kThreeVsThree, thresh14()));

  auto up = run(cfg, {4000, 4000, 4000, 1000, 1000, 1000});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(up.verdict)), "differential")
  TEST_REAL_SIMILAR(up.fold_change, 4.0)
  TEST_REAL_SIMILAR(up.condition_means[0], 4000.0)
  TEST_REAL_SIMILAR(up.condition_means[1], 1000.0)
  // ADR-0056: the p-value is computed under every test. Exact repeats have no variance, so it is 0
  // for different means; the fold change of equal channels is the same geometric or arithmetic.
  TEST_REAL_SIMILAR(up.p_value, 0.0)

  // Symmetric: the test is (fc > t || 1/fc > t), so the direction must not matter.
  auto down = run(cfg, {1000, 1000, 1000, 4000, 4000, 4000});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(down.verdict)), "differential")
  TEST_REAL_SIMILAR(down.fold_change, 0.25)

  auto flat = run(cfg, {1000, 1000, 1000, 1000, 1000, 1000});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(flat.verdict)), "not_differential")
  TEST_REAL_SIMILAR(flat.fold_change, 1.0)
  TEST_REAL_SIMILAR(flat.p_value, 1.0)

  // Just inside the threshold, so the comparison is strict rather than approximate.
  auto just_under = run(cfg, {1300, 1300, 1300, 1000, 1000, 1000});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(just_under.verdict)), "not_differential")
}
END_SECTION

START_SECTION(a_wholly_absent_condition_is_differential_with_no_finite_ratio)
{
  // A species present in one condition and absent in the other is the strongest result the
  // experiment can produce. The pre-ADR-0038 gate REJECTED it -- any empty channel killed the
  // spectrum -- which is the hole the (unreachable) only_one_condition flag existed to patch.
  Config cfg(quantJson(kThreeVsThree, thresh14()));

  auto b_absent = run(cfg, {5000, 5000, 5000, 0, 0, 0});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(b_absent.verdict)), "differential")
  TEST_REAL_SIMILAR(b_absent.fold_change, -1.0)   // no finite ratio
  TEST_REAL_SIMILAR(b_absent.condition_means[0], 5000.0)
  TEST_REAL_SIMILAR(b_absent.condition_means[1], 0.0)

  auto a_absent = run(cfg, {0, 0, 0, 5000, 5000, 5000});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(a_absent.verdict)), "differential")
  TEST_REAL_SIMILAR(a_absent.fold_change, -1.0)

  // Both absent is no signal, not a maximal one.
  auto none = run(cfg, {0, 0, 0, 0, 0, 0});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(none.verdict)), "incomplete_channels")

  // -1 therefore means two things, and condition_means is what separates them: a NOT-MEASURED row
  // carries no means at all, while this one carries two. That pairing is the logged contract.
  TEST_EQUAL(b_absent.condition_means.size(), 2)
}
END_SECTION

START_SECTION(a_partially_empty_condition_is_incomplete_not_differential)
{
  // Distinct from the wholly-absent case on purpose. One zero folded into a mean biases the ratio
  // downward, so there is no honest number to report -- reporting it as differential would invent a
  // result out of a measurement failure.
  Config cfg(quantJson(kThreeVsThree, thresh14()));
  auto r = run(cfg, {5000, 5000, 0, 1000, 1000, 1000});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(r.verdict)), "incomplete_channels")
  TEST_REAL_SIMILAR(r.fold_change, -1.0)
}
END_SECTION

START_SECTION(a_wholly_absent_condition_against_a_partial_one_is_incomplete)
{
  // ADR-0038, amended 2026-10-15: absence is read only against a FULLY measured condition. The two
  // rules above were stated without a precedence, and checking the wholly-absent one first let a
  // single reporter in one of three replicates, against an empty condition, read as the strongest
  // result the experiment can produce -- and buy an identification scan (seen on hardware).
  Config cfg(quantJson(kThreeVsThree, thresh14()));

  auto b_absent_a_partial = run(cfg, {0, 1489, 0, 0, 0, 0});  // the hardware row
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(b_absent_a_partial.verdict)), "incomplete_channels")
  TEST_REAL_SIMILAR(b_absent_a_partial.fold_change, -1.0)

  // The mirror: wholly_absent[0] is a separate arm of the same test.
  auto a_absent_b_partial = run(cfg, {0, 0, 0, 0, 5000, 0});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(a_absent_b_partial.verdict)), "incomplete_channels")

  // Two of three is still partial.
  auto two_of_three = run(cfg, {5000, 5000, 0, 0, 0, 0});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(two_of_three.verdict)), "incomplete_channels")
}
END_SECTION

START_SECTION(a_channel_in_no_condition_is_ignored_by_the_gate_but_still_logged)
{
  // Four samples in six-plex chemistry, or an 11-plex bridge channel. The old gate read EVERY
  // channel in the scheme, so a kit run below capacity rejected every spectrum for its whole run.
  // Only channels named in a condition are gated now; the rest are reported and otherwise ignored.
  Config cfg(quantJson(R"([ { "name": "a", "channels": ["126","127"] },
                            { "name": "b", "channels": ["129","130"] } ])",
                       thresh14()));

  // 128 and 131 are assigned to nothing and are empty. This must still measure cleanly.
  auto r = run(cfg, {4000, 4000, 0, 1000, 1000, 0});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(r.verdict)), "differential")
  TEST_REAL_SIMILAR(r.fold_change, 4.0)

  // ...and every channel is reported regardless, including the unassigned empty ones: gating and
  // logging are deliberately different questions.
  TEST_EQUAL(r.channels.size(), 6)
  TEST_REAL_SIMILAR(r.channels[2], 0.0)
  TEST_REAL_SIMILAR(r.channels[5], 0.0)
}
END_SECTION

START_SECTION(conditions_need_not_be_contiguous_or_equal_sized)
{
  // The retired positional split could only ever express first-three / last-three. Naming channels
  // makes an interleaved or asymmetric design expressible, which is the point of the reshape.
  Config cfg(quantJson(R"([ { "name": "a", "channels": ["126","128","130"] },
                            { "name": "b", "channels": ["127","129","131"] } ])",
                       thresh14()));
  // Interleaved: a = {4000,4000,4000}, b = {1000,1000,1000} despite alternating positions.
  auto r = run(cfg, {4000, 1000, 4000, 1000, 4000, 1000});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(r.verdict)), "differential")
  TEST_REAL_SIMILAR(r.fold_change, 4.0)

  // Asymmetric sizes: means, not sums, so 2-vs-4 is still a fair ratio.
  Config asym(quantJson(R"([ { "name": "a", "channels": ["126","127"] },
                             { "name": "b", "channels": ["128","129","130","131"] } ])",
                        thresh14()));
  auto r2 = run(asym, {4000, 4000, 1000, 1000, 1000, 1000});
  TEST_REAL_SIMILAR(r2.condition_means[0], 4000.0)
  TEST_REAL_SIMILAR(r2.condition_means[1], 1000.0)
  TEST_REAL_SIMILAR(r2.fold_change, 4.0)
}
END_SECTION

START_SECTION(the_ratio_direction_follows_the_conditions_ARRAY_order)
{
  // The array order IS the numerator. Authored as an object, nlohmann's std::map would sort the two
  // alphabetically and the numerator would depend on what the author named their groups -- the trap
  // Config.cpp already documents for additional_ms2 dispatch order.
  const char* ab = R"([ { "name": "zzz", "channels": ["126","127","128"] },
                        { "name": "aaa", "channels": ["129","130","131"] } ])";
  const char* ba = R"([ { "name": "aaa", "channels": ["129","130","131"] },
                        { "name": "zzz", "channels": ["126","127","128"] } ])";
  const std::vector<double> ints {4000, 4000, 4000, 1000, 1000, 1000};

  // Names are deliberately anti-alphabetical to the array order: if anything sorted, these two
  // would agree instead of being reciprocals.
  auto first  = run(Config(quantJson(ab, thresh14())), ints);
  auto second = run(Config(quantJson(ba, thresh14())), ints);
  TEST_REAL_SIMILAR(first.fold_change, 4.0)
  TEST_REAL_SIMILAR(second.fold_change, 0.25)
}
END_SECTION

START_SECTION(isotope_correction_is_applied_and_an_all_zero_matrix_disables_it)
{
  // Correction has never been applied on this path -- IsobaricQuantifier was never constructed --
  // and it matters for the N/C schemes ADR-0038 enables, where neighbouring channels sit 6.32 mDa
  // apart and cross-talk mixes the two conditions together.
  //
  // An all-zero matrix is the identity (IsobaricQuantitationMethod.cpp skips "0.0"/"-1"/"NA" and
  // leaves self_contribution at 100), so it is the documented way to turn correction off. Comparing
  // stock against identity is what proves correction is running at all.
  // Intensities chosen to be all-different, so a correction that redistributes between neighbours
  // cannot coincidentally land back on the raw values.
  const std::vector<double> ints {4000, 3000, 2000, 1500, 1000, 800};

  // The ONLY place in this file that leaves the stock matrix in place -- everywhere else asserts
  // exact arithmetic and therefore turns correction off (see kNoCorrection).
  auto stock = run(Config(quantJson(kThreeVsThree, thresh14Corrected())), ints);
  auto off   = run(Config(quantJson(kThreeVsThree, thresh14())), ints);

  TEST_EQUAL(stock.channels.size(), 6)
  TEST_EQUAL(off.channels.size(), 6)
  // Uncorrected channels are the raw intensities.
  TEST_REAL_SIMILAR(off.channels[0], 4000.0)
  TEST_REAL_SIMILAR(off.channels[5], 800.0)
  // Corrected ones are not -- if this ever becomes equal, correction silently stopped running.
  TEST_EQUAL(stock.channels[0] != off.channels[0] || stock.channels[1] != off.channels[1], true)
}
END_SECTION

START_SECTION(an_unusable_reporter_region_is_reported_not_guessed)
{
  // A spectrum with no reporter ions at all. This is what an activation that cannot release the
  // TMT reporter produces -- the ETD misconfiguration ADR-0038 was written about. It must be
  // REPORTED, because with the activation and tolerance guards deliberately declined, quant_verdict
  // is the only route by which such a run is ever diagnosed.
  Config cfg(quantJson(kThreeVsThree, thresh14()));
  auto r = run(cfg, {0, 0, 0, 0, 0, 0});
  TEST_EQUAL(r.verdict != Quantification::Result::Verdict::Differential, true)
  TEST_EQUAL(r.verdict != Quantification::Result::Verdict::NotDifferential, true)
  TEST_REAL_SIMILAR(r.fold_change, -1.0)
}
END_SECTION

START_SECTION(labelling_vocabulary_matches_the_shipped_schemes)
{
  // Borrowed from TopDownIsobaricQuantification's own setValidStrings MINUS "none". Drift here
  // would let a config name a scheme the factory cannot build, which measure() can only report as
  // extraction_failed -- a silent quant-dead run.
  const auto& names = Quantification::labellingNames();
  TEST_EQUAL(names.size(), 7)
  for (const auto& n : names)
    TEST_EQUAL(Quantification::channelNames(n).empty(), false)

  TEST_EQUAL(Quantification::channelNames("tmt6plex").size(), 6)
  TEST_EQUAL(Quantification::channelNames("tmt10plex").size(), 10)
  TEST_EQUAL(Quantification::channelNames("tmt18plex").size(), 18)
  TEST_EQUAL(Quantification::channelNames("itraq4plex").size(), 4)
  // "none" is NOT a labelling here: quantification.enabled is the switch.
  TEST_EQUAL(Quantification::channelNames("none").empty(), true)
}
END_SECTION

///////////////////////////////////////////////////////////////////////////////
//
// ADR-0056 -- the authored test, the geometric ratio, the empty condition, and the group aggregate.
//
// Every expected number below was computed with scipy.stats.ttest_ind(equal_var=False) on log2
// channels -- the offline pipeline's call -- not with this engine, so the literals are independent
// of the implementation they check. The hardware vectors are LOGGED channels of the 2026-10-13
// `separate` R1 run (already isotope-corrected there), so they are fed with correction OFF.
//
///////////////////////////////////////////////////////////////////////////////

// Real Q scans, 2026-10-13 separate R1 (glucose 126-128 / acetate 129-131).
const std::vector<double> kScan13301 {322572, 328683, 255722, 144385, 147727, 150777};       // 7838 Da z9
const std::vector<double> kScan9107  {243326, 303897, 297053, 332567, 343227, 385884};       // 10393 Da z13
const std::vector<double> kScan9108  {259410, 306478, 296230, 349575, 355020, 374754};       // 10393 Da z12
const std::vector<double> kScan13373 {302887, 291906, 233010, 126518, 137301, 151060};       // 18295 Da z23
const std::vector<double> kScan13375 {204293, 204299, 165680, 85729.7, 83777, 95913.6};      // 18295 Da z22
// The synthetic fixture's RAW reporter peaks (ms2_quant_tmt.txt).
const std::vector<double> kSynthetic {3866.58, 9540.50, 1857.11, 13153.43, 9235.59, 7387.50};
// kScan13373 with the acetate condition wholly absent, and with one acetate channel missing.
const std::vector<double> kAcetateAbsent  {302887, 291906, 233010, 0, 0, 0};
const std::vector<double> kAcetatePartial {302887, 291906, 233010, 126518, 0, 151060};
const std::vector<double> kNoReporters    {0, 0, 0, 0, 0, 0};

auto withTest = [](const std::string& test, const std::string& more = std::string()) {
  return std::string(R"(, "test": ")") + test + R"(")" + more + kNoCorrection;
};

// A separate-mode config, so a min_usable_scans above one loads (ADR-0056 decision 7).
auto groupJson = [](const std::string& extra) {
  return std::string(R"({
    "deconvolution": { "tol": [10, 10, 10] },
    "precursor_selection": { "rank_by": "qscore", "max_precursors": 5, "precursor_charges": "separate",
                             "min_charge_states": 2, "max_charge_states": 3 },
    "characterization": { "mode": "off" },
    "quantification": { "enabled": true, "labelling": "tmt6plex", "reporter_mz_tol": 0.002,
                        "conditions": )") + kThreeVsThree + extra + kNoCorrection + R"( },
    "ms_settings": {
      "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
      "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 },
      "ms2_quant": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 30 }
    }
  })";
};

// One group member's command: tracking id, charge, window_snr, intensity.
auto memberCmd = [](int tracking_id, int charge, double window_snr, double intensity) {
  ScanCommand c {};
  c.scan_id = tracking_id;
  c.num_stages = 1;
  c.stages[0].charge_state = charge;
  c.window_snr = window_snr;
  c.precursor_intensity = intensity;
  return c;
};

using V = Quantification::Result::Verdict;

START_SECTION(welch_p_is_the_pipelines_number)
{
  Config cfg(quantJson(kThreeVsThree, withTest("welch")));
  auto r = run(cfg, kScan13301);
  TEST_EQUAL(r.complete, true)
  TEST_REAL_SIMILAR(r.fold_change, 2.03525259)
  TEST_REAL_SIMILAR(r.p_value, 0.01120125581)
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(r.verdict)), "differential")
}
END_SECTION

START_SECTION(the_three_tests_decide_differently_on_the_same_channels)
{
  // The synthetic fixture: a large ratio (1/FC 2.36) on scattered channels (p 0.204).
  TEST_EQUAL(run(Config(quantJson(kThreeVsThree, withTest("fold_change"))), kSynthetic).verdict == V::Differential, true)
  TEST_EQUAL(run(Config(quantJson(kThreeVsThree, withTest("welch"))), kSynthetic).verdict == V::NotDifferential, true)
  TEST_EQUAL(run(Config(quantJson(kThreeVsThree, withTest("welch_and_fold_change"))), kSynthetic).verdict
                 == V::NotDifferential, true)
  // A real scan the other way round: a small ratio (1/FC 1.25) on tight channels (p 0.032).
  TEST_EQUAL(run(Config(quantJson(kThreeVsThree, withTest("fold_change"))), kScan9108).verdict == V::NotDifferential, true)
  TEST_EQUAL(run(Config(quantJson(kThreeVsThree, withTest("welch"))), kScan9108).verdict == V::Differential, true)
  TEST_EQUAL(run(Config(quantJson(kThreeVsThree, withTest("welch_and_fold_change"))), kScan9108).verdict
                 == V::NotDifferential, true)
  // Both numbers are computed under every test; only the verdict differs.
  auto fc = run(Config(quantJson(kThreeVsThree, withTest("fold_change"))), kSynthetic);
  TEST_REAL_SIMILAR(fc.fold_change, 0.4242073633)
  TEST_REAL_SIMILAR(fc.p_value, 0.2043606735)
}
END_SECTION

START_SECTION(condition_means_are_geometric_and_the_ratio_is_their_quotient)
{
  // A skewed triple separates the means: arithmetic 34000 vs geometric 4641.6.
  auto r = run(Config(quantJson(kThreeVsThree, withTest("fold_change"))), {1000, 1000, 100000, 10000, 10000, 10000});
  TEST_REAL_SIMILAR(r.condition_means[0], 4641.588834)
  TEST_REAL_SIMILAR(r.condition_means[1], 10000.0)
  TEST_REAL_SIMILAR(r.fold_change, r.condition_means[0] / r.condition_means[1])
  TEST_REAL_SIMILAR(r.fold_change, 0.4641588834)
  // Direction is the ratio's: an arithmetic mean would have called condition a enriched here.
  TEST_EQUAL(r.condition_means[0] < r.condition_means[1], true)
  TEST_EQUAL(r.verdict == V::Differential, true)
}
END_SECTION

START_SECTION(an_empty_condition_is_authored)
{
  for (const std::string test : {"fold_change", "welch", "welch_and_fold_change"})
  {
    auto diff = run(Config(quantJson(kThreeVsThree, withTest(test, R"(, "empty_condition": "differential")"))),
                    kAcetateAbsent);
    TEST_EQUAL(diff.verdict == V::Differential, true)
    TEST_EQUAL(diff.complete, false)
    TEST_REAL_SIMILAR(diff.fold_change, -1.0)
    TEST_REAL_SIMILAR(diff.p_value, -1.0)
    TEST_EQUAL(diff.condition_means[0] > 0.0, true)
    TEST_REAL_SIMILAR(diff.condition_means[1], 0.0)

    auto inc = run(Config(quantJson(kThreeVsThree, withTest(test, R"(, "empty_condition": "incomplete")"))),
                   kAcetateAbsent);
    TEST_EQUAL(inc.verdict == V::IncompleteChannels, true)

    // A PARTIAL condition is never an empty condition, under either setting.
    for (const std::string ec : {"differential", "incomplete"})
    {
      auto part = run(Config(quantJson(kThreeVsThree, withTest(test, R"(, "empty_condition": ")" + ec + R"(")"))),
                      kAcetatePartial);
      TEST_EQUAL(part.verdict == V::IncompleteChannels, true)
    }
  }
}
END_SECTION

START_SECTION(zero_variance_takes_p_0_or_1)
{
  // Exact repeats have no variance: t is undefined, and scipy's answer is p 0 for different means
  // and p 1 for equal ones. The value used is the value logged.
  auto apart = run(Config(quantJson(kThreeVsThree, withTest("welch"))), {4000, 4000, 4000, 1000, 1000, 1000});
  TEST_REAL_SIMILAR(apart.p_value, 0.0)
  TEST_EQUAL(apart.verdict == V::Differential, true)
  auto same = run(Config(quantJson(kThreeVsThree, withTest("welch"))), {1000, 1000, 1000, 1000, 1000, 1000});
  TEST_REAL_SIMILAR(same.p_value, 1.0)
  TEST_EQUAL(same.verdict == V::NotDifferential, true)
}
END_SECTION

START_SECTION(p_is_logged_under_fold_change_too_and_minus_one_without_two_channels)
{
  auto three = run(Config(quantJson(kThreeVsThree, withTest("fold_change"))), kScan13301);
  TEST_REAL_SIMILAR(three.p_value, 0.01120125581)
  // One channel in a condition: legal under fold_change, no variance, no p.
  const char* one_v_three = R"([ { "name": "a", "channels": ["126"] }, { "name": "b", "channels": ["129","130","131"] } ])";
  auto one = run(Config(quantJson(one_v_three, withTest("fold_change"))), kScan13301);
  TEST_REAL_SIMILAR(one.p_value, -1.0)
  TEST_EQUAL(one.complete, true)
  TEST_EQUAL(one.verdict == V::Differential, true)
}
END_SECTION

START_SECTION(a_group_is_decided_on_the_median_of_shares_not_by_vote)
{
  // Two charge states of one real species. Alone, z13 is NOT a Welch hit (p 0.061) and z12 IS
  // (p 0.032): a vote would have tied 1-1. The aggregate is a hit (p 0.043).
  Config welch(groupJson(R"(, "test": "welch")"));
  Quantification q(welch);
  q.addMember(7, memberCmd(101, 13, 2.0, 1e6), 501);
  q.addMember(7, memberCmd(102, 12, 1.0, 1e6), 502);
  TEST_EQUAL(q.deposit(101, run(welch, kScan9107)), 7)
  TEST_EQUAL(q.isComplete(7), false)
  TEST_EQUAL(q.deposit(102, run(welch, kScan9108)), 7)
  TEST_EQUAL(q.isComplete(7), true)
  auto g = q.decide(7);
  TEST_EQUAL(g.result.verdict == V::Differential, true)
  TEST_REAL_SIMILAR(g.result.fold_change, 0.7950028757)
  TEST_REAL_SIMILAR(g.result.p_value, 0.04286976721)
  TEST_EQUAL(g.usable, 2)
  TEST_EQUAL(g.members, 2)
  TEST_EQUAL(g.has_winner, true)
  TEST_EQUAL(g.winner_charge, 13)        // the higher window_snr
  TEST_EQUAL(g.winner_precursor_id, 501)
  TEST_EQUAL(g.result.channels.size(), 6)
  TEST_REAL_SIMILAR(g.result.channels[0], 0.1306408534)
  TEST_REAL_SIMILAR(g.result.channels[5], 0.1977442951)

  // Same group under welch_and_fold_change: a hit by p, not by ratio (1/FC 1.26 < 1.4).
  Config both(groupJson(R"(, "test": "welch_and_fold_change")"));
  Quantification q2(both);
  q2.addMember(8, memberCmd(201, 13, 2.0, 1e6), 601);
  q2.addMember(8, memberCmd(202, 12, 1.0, 1e6), 602);
  q2.deposit(201, run(both, kScan9107));
  q2.deposit(202, run(both, kScan9108));
  auto g2 = q2.decide(8);
  TEST_EQUAL(g2.result.verdict == V::NotDifferential, true)
  TEST_EQUAL(g2.has_winner, true)
}
END_SECTION

START_SECTION(aggregate_sum_mean_median_and_scan_normalize)
{
  auto decideTwo = [&](const std::string& extra) {
    Config cfg(groupJson(extra));
    Quantification q(cfg);
    q.addMember(1, memberCmd(11, 23, 1.8, 6.2e6), 1);
    q.addMember(1, memberCmd(12, 22, 1.8, 6.0e6), 2);
    q.deposit(11, run(cfg, kScan13373));
    q.deposit(12, run(cfg, kScan13375));
    return q.decide(1);
  };
  // Even count: the median is the mean of the two middle values (= the mean of two shares).
  auto med = decideTwo(R"(, "aggregate": "median")");
  TEST_REAL_SIMILAR(med.result.fold_change, 2.069697842)
  TEST_REAL_SIMILAR(med.result.p_value, 0.002788670769)
  TEST_REAL_SIMILAR(med.result.channels[0], 0.2435158024)
  // window_snr ties: intensity breaks it, as before.
  TEST_EQUAL(med.winner_charge, 23)

  // Without normalisation, sum and mean give the SAME ratio and p -- they differ by a constant per
  // channel -- so only the logged channels can tell a `sum` that runs `mean` from a real one.
  auto sum = decideTwo(R"(, "scan_normalize": false, "aggregate": "sum")");
  auto mean = decideTwo(R"(, "scan_normalize": false, "aggregate": "mean")");
  TEST_REAL_SIMILAR(sum.result.fold_change, 2.053398224)
  TEST_REAL_SIMILAR(sum.result.p_value, 0.002989884053)
  TEST_REAL_SIMILAR(mean.result.fold_change, 2.053398224)
  TEST_REAL_SIMILAR(sum.result.channels[0], 507180.0)
  TEST_REAL_SIMILAR(mean.result.channels[0], 253590.0)

  // Three members: the median picks the MIDDLE share per channel, never a mean.
  Config cfg(groupJson(R"(, "aggregate": "median")"));
  Quantification q(cfg);
  q.addMember(2, memberCmd(21, 9, 1.0, 1e6), 1);
  q.addMember(2, memberCmd(22, 23, 1.0, 1e6), 1);
  q.addMember(2, memberCmd(23, 22, 1.0, 1e6), 1);
  q.deposit(21, run(cfg, kScan13301));
  q.deposit(22, run(cfg, kScan13373));
  q.deposit(23, run(cfg, kScan13375));
  auto three = q.decide(2);
  TEST_EQUAL(three.usable, 3)
  TEST_REAL_SIMILAR(three.result.channels[0], 0.2432950737)
  TEST_REAL_SIMILAR(three.result.channels[2], 0.1894425076)
  TEST_REAL_SIMILAR(three.result.fold_change, 2.063504791)
  TEST_REAL_SIMILAR(three.result.p_value, 0.006533441894)
}
END_SECTION

START_SECTION(complete_members_decide_and_absence_only_when_none_is_complete)
{
  // A complete member and an absent one: the complete one IS the group (pipeline's COMPLETE_ONLY).
  Config cfg(groupJson(""));
  Quantification q(cfg);
  q.addMember(1, memberCmd(11, 23, 1.0, 1e6), 1);
  q.addMember(1, memberCmd(12, 22, 5.0, 1e6), 2);   // the absent one has the better window
  q.deposit(11, run(cfg, kScan13373));
  q.deposit(12, run(cfg, kAcetateAbsent));
  auto g = q.decide(1);
  TEST_EQUAL(g.usable, 1)
  TEST_EQUAL(g.members, 2)
  TEST_REAL_SIMILAR(g.result.fold_change, 1.98750199)
  TEST_EQUAL(g.winner_charge, 23)                     // usable members only, whatever its window

  // No complete member, two absent against full: the empty condition decides.
  auto twoAbsent = [&](const std::string& extra) {
    Config c(groupJson(extra));
    Quantification qq(c);
    qq.addMember(3, memberCmd(31, 23, 1.0, 1e6), 1);
    qq.addMember(3, memberCmd(32, 22, 4.0, 1e6), 2);
    qq.deposit(31, run(c, kAcetateAbsent));
    qq.deposit(32, run(c, kAcetateAbsent));
    return qq.decide(3);
  };
  auto diff = twoAbsent(R"(, "empty_condition": "differential")");
  TEST_EQUAL(diff.result.verdict == V::Differential, true)
  TEST_REAL_SIMILAR(diff.result.fold_change, -1.0)
  TEST_EQUAL(diff.has_winner, true)
  TEST_EQUAL(diff.winner_charge, 22)
  TEST_EQUAL(diff.usable, 0)
  auto inc = twoAbsent(R"(, "empty_condition": "incomplete")");
  TEST_EQUAL(inc.result.verdict == V::IncompleteChannels, true)
  TEST_EQUAL(inc.has_winner, true)                    // undecided, still a window: identify all buys it
}
END_SECTION

START_SECTION(the_floor_applies_to_both_paths)
{
  auto decide = [&](const std::string& extra, const std::vector<double>& a, const std::vector<double>& b) {
    Config c(groupJson(extra));
    Quantification q(c);
    q.addMember(1, memberCmd(11, 15, 2.0, 1e6), 1);
    q.addMember(1, memberCmd(12, 16, 1.0, 1e6), 2);
    q.deposit(11, run(c, a));
    q.deposit(12, run(c, b));
    return q.decide(1);
  };
  const std::string floor2 = R"(, "min_usable_scans": 2)";
  // A: one complete + one absent, floor 2 -> not enough usable members.
  auto A = decide(floor2, kScan13373, kAcetateAbsent);
  TEST_EQUAL(A.result.verdict == V::IncompleteChannels, true)
  TEST_EQUAL(A.usable, 1)
  TEST_EQUAL(A.has_winner, true)
  // B: two absent against full, floor 2 -> the absence path is met.
  auto B = decide(floor2, kAcetateAbsent, kAcetateAbsent);
  TEST_EQUAL(B.result.verdict == V::Differential, true)
  // C: one absent + one with no reporters, floor 2 -> the absence path is NOT met.
  auto C = decide(floor2, kAcetateAbsent, kNoReporters);
  TEST_EQUAL(C.result.verdict == V::IncompleteChannels, true)
  // Floor 1: A is decided by its complete member, C by its absent one.
  auto A1 = decide("", kScan13373, kAcetateAbsent);
  TEST_EQUAL(A1.result.verdict == V::Differential, true)
  TEST_REAL_SIMILAR(A1.result.fold_change, 1.98750199)
  auto C1 = decide("", kAcetateAbsent, kNoReporters);
  TEST_EQUAL(C1.result.verdict == V::Differential, true)
  TEST_EQUAL(C1.winner_charge, 15)
}
END_SECTION

START_SECTION(an_undecided_group_still_names_its_identification_window)
{
  // Neither member is usable nor absent-against-full: incomplete, and the window is still named,
  // so `identify: all` can buy it (ADR-0056 decision 9).
  Config cfg(groupJson(""));
  Quantification q(cfg);
  q.addMember(1, memberCmd(11, 15, 1.0, 1e6), 1);
  q.addMember(1, memberCmd(12, 16, 3.0, 1e6), 2);
  q.deposit(11, run(cfg, kAcetatePartial));
  q.deposit(12, run(cfg, kNoReporters));
  auto g = q.decide(1);
  TEST_EQUAL(g.result.verdict == V::IncompleteChannels, true)
  TEST_EQUAL(g.has_winner, true)
  TEST_EQUAL(g.winner_charge, 16)
  TEST_EQUAL(g.usable, 0)
}
END_SECTION

START_SECTION(an_unusable_member_is_absent_from_the_aggregate_not_opposed_to_it)
{
  Config cfg(groupJson(""));
  Quantification q(cfg);
  q.addMember(1, memberCmd(11, 15, 1.0, 1e6), 1);
  q.addMember(1, memberCmd(12, 16, 9.0, 1e6), 2);
  q.deposit(11, run(cfg, kScan13301));
  q.deposit(12, run(cfg, kNoReporters));
  auto g = q.decide(1);
  TEST_EQUAL(g.usable, 1)
  TEST_EQUAL(g.members, 2)
  TEST_REAL_SIMILAR(g.result.fold_change, 2.03525259)
  TEST_REAL_SIMILAR(g.result.p_value, 0.01120125581)
  TEST_EQUAL(g.winner_charge, 15)
  TEST_EQUAL(g.usable_charges.size(), 1)
  TEST_EQUAL(g.usable_charges[0], 15)
}
END_SECTION

START_SECTION(a_sub_count_channel_is_present_and_dominates_a_geometric_mean)
{
  // The presence floor stays 1e-3 (owner, 2026-10-15): 0.01 counts is PRESENT, and under log2 it
  // pulls condition a's mean from 6667 (arithmetic) to 100. Pinned so it is a known quantity.
  auto r = run(Config(quantJson(kThreeVsThree, withTest("fold_change"))), {10000, 10000, 0.01, 10000, 10000, 10000});
  TEST_EQUAL(r.complete, true)
  TEST_REAL_SIMILAR(r.condition_means[0], 100.0)
  TEST_REAL_SIMILAR(r.fold_change, 0.01)
  TEST_EQUAL(r.verdict == V::Differential, true)
}
END_SECTION

START_SECTION(a_partial_condition_logs_zero_means)
{
  // A geometric mean of a condition with a 0 channel is undefined; it is logged as 0, never NaN.
  auto r = run(Config(quantJson(kThreeVsThree, withTest("fold_change"))), {5000, 5000, 0, 1000, 1000, 1000});
  TEST_EQUAL(r.verdict == V::IncompleteChannels, true)
  TEST_REAL_SIMILAR(r.condition_means[0], 0.0)
  TEST_REAL_SIMILAR(r.condition_means[1], 1000.0)
  TEST_REAL_SIMILAR(r.fold_change, -1.0)
  TEST_REAL_SIMILAR(r.p_value, -1.0)
}
END_SECTION

END_TEST
