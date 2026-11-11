// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

/**
  @brief The settings of one MS3 yield gate (ADR-0069), split out of Ms3YieldModel.h so Config.h can hold them by value.

  A header of its own for the reason RefineSettings.h is one: the model's input builders need ScanCommand.h and
  FragmentAnalysis.h, and FragmentAnalysis.h includes Config.h. This file includes nothing from FLASHIda.
*/
namespace OpenMS
{
  /// The six coefficients of the predicted-bond-yield model, defaulted to the core-5 fit (MS3_gate_spec.md, 2026-11-09:
  /// two Caco-2 histone live-ID runs, CID 35, 240k, 4 microscans, AGC 5e5, max_it 502 ms). The transforms are fixed in
  /// the model: eta = intercept + fragment_intensity * log10(I_frag) + precursor_intensity * log10(I_prec)
  /// + precursor_charge * z_prec + mass_ratio * (M_frag / M_prec) + fragment_iso_cos * iso_cos_frag.
  struct OPENMS_DLLAPI YieldCoefficients
  {
    double intercept = -31.632174;
    double fragment_intensity = 1.529105;
    double precursor_intensity = 0.474888;
    double precursor_charge = -0.105683;
    double mass_ratio = -2.766660;
    double fragment_iso_cos = 21.747494;
  };

  /// identification.characterization.yield_gate / quantification.ms3_yield_gate: an MS3 is admitted when its predicted
  /// bond yield E = p * (L - 1) is at least `threshold` expected bonds.
  /// yield_gate.rank_by: the order a gated plan walks its candidates in. Yield: by predicted bond yield, highest
  /// first, a target with no ion index after every tested one. Intensity: strongest first with the yield as a
  /// pass/fail test only -- the walk before the key, kept for runs that collect training data.
  enum class YieldRank
  {
    Yield = 0,
    Intensity
  };

  struct OPENMS_DLLAPI YieldGateConfig
  {
    double threshold = 1.0;
    YieldCoefficients coefficients;
    YieldRank rank_by = YieldRank::Yield;
  };
} // namespace OpenMS
