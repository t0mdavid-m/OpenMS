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

  // Symmetric: the test is (fc > t || 1/fc > t), so the direction must not matter.
  auto down = run(cfg, {1000, 1000, 1000, 4000, 4000, 4000});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(down.verdict)), "differential")
  TEST_REAL_SIMILAR(down.fold_change, 0.25)

  auto flat = run(cfg, {1000, 1000, 1000, 1000, 1000, 1000});
  TEST_STRING_EQUAL(std::string(Quantification::verdictName(flat.verdict)), "not_differential")
  TEST_REAL_SIMILAR(flat.fold_change, 1.0)

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

END_TEST
