// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// ADR-0069: the predicted bond yield, E = p x (L - 1), the MS3 yield gate decides by (Ms3YieldModel.h). The formula
// against the specification's own test vectors, the clipping of every input to its training range, the gate's >=, the
// untested target with no ion index, and the drift guard that makes the scan_commands row's predicted_bond_yield the
// decision value: the inputs read back off a built MS3 equal the planner's inputs, field by field, for a single MS3,
// an MS3 sweep's production MS3 (built from a variant, Exploration.cpp) and a quantification 'Q' MS3.

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3YieldModel.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include <cmath>
#include <cstring>
#include <limits>
#include <string>

using namespace OpenMS;

namespace
{
  // The specification's test vectors (MS3_gate_spec.md section 4, core-5).
  YieldInputs vec(double i_frag, double i_prec, double z, double m_frag, double m_prec, double iso_cos, int length)
  {
    YieldInputs in;
    in.i_frag = i_frag;
    in.i_prec = i_prec;
    in.z_prec = z;
    in.m_frag = m_frag;
    in.m_prec = m_prec;
    in.iso_cos = iso_cos;
    in.length = length;
    return in;
  }

  double etaOf(const YieldPrediction& pred) { return std::log(pred.p / (1.0 - pred.p)); }

  // A mid-range input every clipping section varies one field of.
  YieldInputs mid() { return vec(1.0e4, 5.0e5, 10, 4000.0, 10000.0, 0.95, 30); }

  // No runtime block / no log files.
  const char* kConfig = R"({
  "deconvolution": { "min_charge": 1, "max_charge": 50, "min_mass": 100, "max_mass": 50000, "tol": [10, 10, 10] },
  "precursor_selection": { "rank_by": "qscore", "max_precursors": 3 },
  "identification": { "reference": { "sequence": "PEPTCDEK" } },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000 },
    "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 30 }
  }
})";

  // The dispatching MS2: every stage-0 input set to a distinctive value.
  ScanCommand ms2Ctx()
  {
    ScanCommand ctx{};
    ctx.scan_id = 1;
    ctx.msn_level = 2;
    ctx.num_stages = 1;
    ctx.mono_mass = 15321.4567;
    ctx.precursor_intensity = 818339.0;
    ctx.iso_cos = 0.91;
    ctx.stages[0].precursor_mz = 1179.58;
    ctx.stages[0].isolation_width = 2.0;
    ctx.stages[0].charge_state = 13;
    ctx.stages[0].collision_energy = 30;
    std::strncpy(ctx.stages[0].activation_type, "HCD", sizeof(ctx.stages[0].activation_type) - 1);
    return ctx;
  }

  FragmentAnalysis::FragmentScores fragScores()
  {
    FragmentAnalysis::FragmentScores s;
    s.mono_mass = 13194.8123;
    s.iso_cos = 0.9886;
    s.precursor_intensity = 13194.8;
    s.qscore = 0.7;
    return s;
  }

  bool sameInputs(const YieldInputs& a, const YieldInputs& b)
  {
    return a.i_frag == b.i_frag && a.i_prec == b.i_prec && a.z_prec == b.z_prec && a.m_frag == b.m_frag
           && a.m_prec == b.m_prec && a.iso_cos == b.iso_cos && a.length == b.length;
  }
} // namespace

START_TEST(FLASHIda_Ms3YieldModel, "$Id$")

// T1: the specification's three rows, to the digits it prints.
START_SECTION(spec_test_vectors)
{
  const YieldCoefficients c;   // core-5
  TOLERANCE_ABSOLUTE(6e-5)   // eta and E are printed to 4 decimals; E = p x (L - 1) then pins p to ~2e-6
  {
    const YieldPrediction r = predictBondYield(c, vec(13194.8, 818339, 4, 3444.8950, 3721.0067, 0.9886, 30));
    TEST_EQUAL(r.tested, true)
    TEST_REAL_SIMILAR(etaOf(r), -4.0082)
    TEST_REAL_SIMILAR(r.p, 0.017842)
    TEST_REAL_SIMILAR(r.e, 0.5174)
  }
  {
    const YieldPrediction r = predictBondYield(c, vec(1306.96, 66399.5, 8, 1989.2173, 6153.3083, 0.967746, 18));
    TEST_REAL_SIMILAR(etaOf(r), -5.2709)
    TEST_REAL_SIMILAR(r.p, 0.005113)
    TEST_REAL_SIMILAR(r.e, 0.0869)
  }
  {
    const YieldPrediction r = predictBondYield(c, vec(200000, 5000000, 6, 2000, 8000, 0.99, 18));
    TEST_REAL_SIMILAR(etaOf(r), -0.1408)
    TEST_REAL_SIMILAR(r.p, 0.464852)
    TEST_REAL_SIMILAR(r.e, 7.9025)
  }
}
END_SECTION

// T2: every input beyond its training range predicts exactly what the range's bound predicts.
START_SECTION(every_input_clips_to_its_training_range)
{
  using R = YieldClipRanges;
  const YieldCoefficients c;
  auto pAt = [&c](YieldInputs in) { return predictBondYield(c, in).p; };
  auto with = [](double YieldInputs::*field, double v) {
    YieldInputs in = mid();
    in.*field = v;
    return in;
  };
  // i_frag, i_prec, z_prec, iso_cos: 10x outside each bound equals the bound.
  TEST_EQUAL(pAt(with(&YieldInputs::i_frag, R::kFragIntensityMin / 10)), pAt(with(&YieldInputs::i_frag, R::kFragIntensityMin)))
  TEST_EQUAL(pAt(with(&YieldInputs::i_frag, R::kFragIntensityMax * 10)), pAt(with(&YieldInputs::i_frag, R::kFragIntensityMax)))
  TEST_EQUAL(pAt(with(&YieldInputs::i_prec, R::kPrecIntensityMin / 10)), pAt(with(&YieldInputs::i_prec, R::kPrecIntensityMin)))
  TEST_EQUAL(pAt(with(&YieldInputs::i_prec, R::kPrecIntensityMax * 10)), pAt(with(&YieldInputs::i_prec, R::kPrecIntensityMax)))
  TEST_EQUAL(pAt(with(&YieldInputs::z_prec, R::kPrecChargeMin / 10)), pAt(with(&YieldInputs::z_prec, R::kPrecChargeMin)))
  TEST_EQUAL(pAt(with(&YieldInputs::z_prec, R::kPrecChargeMax * 10)), pAt(with(&YieldInputs::z_prec, R::kPrecChargeMax)))
  TEST_EQUAL(pAt(with(&YieldInputs::iso_cos, R::kIsoCosMin / 10)), pAt(with(&YieldInputs::iso_cos, R::kIsoCosMin)))
  TEST_EQUAL(pAt(with(&YieldInputs::iso_cos, R::kIsoCosMax * 10)), pAt(with(&YieldInputs::iso_cos, R::kIsoCosMax)))
  // The mass ratio, through m_frag over a fixed m_prec of 10000.
  TEST_EQUAL(pAt(with(&YieldInputs::m_frag, 10000.0 * R::kMassRatioMin / 10)), pAt(with(&YieldInputs::m_frag, 10000.0 * R::kMassRatioMin)))
  TEST_EQUAL(pAt(with(&YieldInputs::m_frag, 10000.0 * R::kMassRatioMax * 10)), pAt(with(&YieldInputs::m_frag, 10000.0 * R::kMassRatioMax)))
  // Inside the range the input moves p: the clip is a clip, not a constant.
  TEST_NOT_EQUAL(pAt(with(&YieldInputs::iso_cos, 0.90)), pAt(with(&YieldInputs::iso_cos, 0.99)))
}
END_SECTION

// T3: an unset score is 0 in the engine; it clips to its range's lower bound, so p stays finite.
START_SECTION(an_unset_zero_clips_to_the_lower_bound)
{
  using R = YieldClipRanges;
  const YieldCoefficients c;
  const YieldPrediction zero = predictBondYield(c, vec(0, 0, 0, 4000, 10000, 0, 30));
  const YieldPrediction low = predictBondYield(c, vec(R::kFragIntensityMin, R::kPrecIntensityMin, R::kPrecChargeMin, 4000, 10000,
                                                     R::kIsoCosMin, 30));
  TEST_EQUAL(std::isfinite(zero.p), true)
  TEST_EQUAL(zero.p, low.p)
  TEST_EQUAL(zero.e, low.e)
}
END_SECTION

// T4: a precursor mass of 0 would divide by zero; the ratio takes its upper bound instead.
START_SECTION(a_zero_precursor_mass_clips_the_ratio_to_one)
{
  const YieldCoefficients c;
  YieldInputs z = mid();
  z.m_prec = 0.0;
  YieldInputs one = mid();
  one.m_frag = one.m_prec;   // ratio exactly 1
  TEST_EQUAL(std::isfinite(predictBondYield(c, z).p), true)
  TEST_EQUAL(predictBondYield(c, z).p, predictBondYield(c, one).p)
  YieldInputs nan = mid();
  nan.m_frag = std::numeric_limits<double>::quiet_NaN();
  TEST_EQUAL(predictBondYield(c, nan).p, predictBondYield(c, one).p)
}
END_SECTION

// T5: the gate admits at E >= threshold, exactly.
START_SECTION(threshold_compares_ge)
{
  const YieldCoefficients c;
  const YieldPrediction pred = predictBondYield(c, mid());
  YieldGateConfig at;
  at.threshold = pred.e;
  TEST_EQUAL(yieldAdmits(at, pred), true)
  YieldGateConfig above;
  above.threshold = std::nextafter(pred.e, std::numeric_limits<double>::infinity());
  TEST_EQUAL(yieldAdmits(above, pred), false)
}
END_SECTION

// T6: a target with no ion index has no L, so the gate does not test it and admits it.
START_SECTION(no_ion_index_is_untested)
{
  const YieldCoefficients c;
  YieldInputs in = mid();
  in.length = 0;
  const YieldPrediction pred = predictBondYield(c, in);
  TEST_EQUAL(pred.tested, false)
  TEST_EQUAL(pred.e, -1.0)
  YieldGateConfig strict;
  strict.threshold = 1.0e9;
  TEST_EQUAL(yieldAdmits(strict, pred), true)
}
END_SECTION

// T7: the drift guard. What the row reads back off a built MS3 is what the planner decided with, field by field.
START_SECTION(row_inputs_equal_plan_inputs)
{
  Config cfg{std::string(kConfig)};
  ScanCommandQueue queue(cfg);
  const ScanCommand ctx = ms2Ctx();
  const FragmentAnalysis::FragmentScores s1 = fragScores();
  ScanConfig ms3;
  ms3.activation = "CID";
  ms3.collision_energy = 35;
  const YieldInputs plan = yieldInputsOf(ctx, s1, 139);   // a three-digit ion at a two-digit charge: the description budget

  // A single MS3, as dispatchFromModel builds it.
  const ScanCommand single = queue.buildMS3(ctx, ms3, 1100.25, 12, 2.0, ctx.scan_id, 'y', 139, 1, s1);
  TEST_EQUAL(sameInputs(yieldInputsOf(single), plan), true)
  TEST_EQUAL(ms3IonOf(single).first, std::string("y"))
  TEST_EQUAL(ms3IonOf(single).second, 139)

  // An MS3 sweep's production MS3: built FROM the winning variant, with the variant's *_s1 fields as its scores
  // (Exploration.cpp, the follow-up after a level-3 sweep).
  FragmentAnalysis::FragmentScores wfs;
  wfs.mono_mass = single.mono_mass_s1;
  wfs.iso_cos = single.iso_cos_s1;
  wfs.precursor_intensity = single.precursor_intensity_s1;
  const ScanCommand production = queue.buildMS3(single, ms3, 1100.25, 12, 2.0, single.scan_id, 'y', 139, 1, wfs);
  TEST_EQUAL(sameInputs(yieldInputsOf(production), plan), true)

  // A quantification MS3, as dispatchQuant builds it.
  const ScanCommand q = queue.buildMS3(ctx, ms3, 1100.25, 12, 2.0, ctx.scan_id, 'y', 139, 1, s1, nullptr, {}, nullptr, 'Q');
  TEST_EQUAL(sameInputs(yieldInputsOf(q), plan), true)

  // The same E both ways.
  const YieldCoefficients c;
  TEST_EQUAL(predictBondYield(c, yieldInputsOf(single)).e, predictBondYield(c, plan).e)

  // A raw mass with no ion: no ion index on the command, so the row's E is untested, as the planner's was.
  const ScanCommand raw = queue.buildMS3(ctx, ms3, 1100.25, 12, 2.0, ctx.scan_id, '\0', 0, 1, s1);
  TEST_EQUAL(yieldInputsOf(raw).length, 0)
  TEST_EQUAL(predictBondYield(c, yieldInputsOf(raw)).tested, false)
}
END_SECTION

END_TEST
