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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/YieldGateSettings.h>

#include <string>
#include <utility>

namespace OpenMS
{
  /**
    @brief The predicted bond yield of an MS3 (ADR-0069): pure functions, no state.

    E = p * (L - 1), the backbone bonds inside the MS3's target fragment its own fragments are expected to witness; p is
    a logistic model of five values known when the MS3 is planned, each clipped to its training range first, and L is
    the target's ion index. The planner decides with yieldInputsOf(ms_ctx, scores, ion_index); the scan_commands row
    logs yieldInputsOf(ms3), the same values read back off the built command, through the same predictBondYield --
    so the logged E is the decision E (FLASHIda_Ms3YieldModel_test row_inputs_equal_plan_inputs).
  */
  /// The five model inputs and the fragment length, as read: raw, unclipped (the events line logs them so).
  struct OPENMS_DLLAPI YieldInputs
  {
    double i_frag = 0;   ///< the fragment's intensity at its charge in the source MS2 (stage-1 precursor_intensity)
    double i_prec = 0;   ///< the source MS2 precursor's intensity (stage-0 precursor_intensity)
    double z_prec = 0;   ///< the source MS2 precursor's charge (stage-0 charge)
    double m_frag = 0;   ///< the fragment's monoisotopic mass (stage-1 mono_mass)
    double m_prec = 0;   ///< the source MS2 precursor's monoisotopic mass (stage-0 mono_mass)
    double iso_cos = 0;  ///< the fragment's isotope cosine (stage-1 iso_cos)
    int length = 0;      ///< L, the target's ion index; <= 0 = a mass that matches no ion, which the gate does not test
  };

  struct OPENMS_DLLAPI YieldPrediction
  {
    bool tested = false;  ///< false iff the target has no ion index (then p and e are -1)
    double p = -1;        ///< the predicted chance that one bond inside the fragment is witnessed
    double e = -1;        ///< the predicted bond yield, p * (L - 1)
  };

  /// The training range each input is clipped to before the model reads it (MS3_gate_spec.md section 2).
  struct OPENMS_DLLAPI YieldClipRanges
  {
    static constexpr double kFragIntensityMin = 372.0, kFragIntensityMax = 5.0e6;
    static constexpr double kPrecIntensityMin = 6.8e3, kPrecIntensityMax = 8.4e7;
    static constexpr double kPrecChargeMin = 4.0, kPrecChargeMax = 49.0;
    static constexpr double kMassRatioMin = 0.011, kMassRatioMax = 1.0;
    static constexpr double kIsoCosMin = 0.85, kIsoCosMax = 1.0005;
  };

  /// The prediction for @p in. A non-finite input (or a mass ratio with a non-positive precursor mass) clips to its
  /// range's upper bound for the ratio and its lower bound otherwise, so the result is always finite.
  OPENMS_DLLAPI YieldPrediction predictBondYield(const YieldCoefficients& c, const YieldInputs& in);

  /// Whether @p gate admits an MS3 predicted @p pred: an untested one always, a tested one when e >= threshold.
  OPENMS_DLLAPI bool yieldAdmits(const YieldGateConfig& gate, const YieldPrediction& pred);

  /// The planner's inputs: the dispatching MS2 command @p ms_ctx (stage 0 of the MS3 buildMS3 will make from it), the
  /// fragment charge's stage-1 scores @p s1 (buildMS3's frag_scores) and the target's @p ion_index.
  OPENMS_DLLAPI YieldInputs yieldInputsOf(const ScanCommand& ms_ctx, const FragmentAnalysis::FragmentScores& s1, int ion_index);

  /// The row's inputs: the same values read back off a built MS3 command.
  OPENMS_DLLAPI YieldInputs yieldInputsOf(const ScanCommand& ms3);

  /// The fragment ion an MS3's scan_description names after its last '@' (("y", 39) for "...@12y39"); ("", 0) for
  /// any other command or an MS3 whose target is a raw mass. The one parser IdaLogger's ion columns and
  /// yieldInputsOf(ms3) share.
  OPENMS_DLLAPI std::pair<std::string, int> ms3IonOf(const ScanCommand& cmd);

  /// The " e=... p=... threshold=... i_frag=... ... length=..." tail every [MS3-YIELD] line carries.
  OPENMS_DLLAPI std::string yieldEventFields(const YieldGateConfig& gate, const YieldInputs& in, const YieldPrediction& pred);

  /// Write one [MS3-YIELD] line to the event log: @p gate_name ("characterization" | "quantification"), the precursor,
  /// the ion and charge tested, @p verdict ("refused" | "untested"; an untested line also carries reason=no_ion_index)
  /// and yieldEventFields. Admitted MS3 write no line: their E is the scan_commands row's predicted_bond_yield.
  OPENMS_DLLAPI void logYieldVerdict(const char* gate_name, int precursor_id, const std::string& ion_type, int ion_index,
                                     int charge, const char* verdict, const YieldGateConfig& gate, const YieldInputs& in,
                                     const YieldPrediction& pred);
} // namespace OpenMS
