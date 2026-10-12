// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <cstdint>

/**
  @brief The settings of one Refine run (ADR-0049), split out of Refine.h so Config.h can hold them by value.

  A header of its own because of an include cycle, not for tidiness: Refine.h needs FragmentAnalysis.h
  (complete types in its declarations), and FragmentAnalysis.h includes Config.h. Config.h including
  Refine.h would therefore close a loop that #pragma once breaks in whichever order a TU happens to
  include them -- one of the two headers would see the other's types undeclared. This file includes
  nothing from FLASHIda, so both Config.h and Refine.h can include it.
*/
namespace OpenMS
{
  /// The settings of one Refine run. All but tolerance_ppm are characterization.refine.*; tolerance_ppm is level(2)'s.
  struct OPENMS_DLLAPI RefineSettings
  {
    double fdr = 0.2;            ///< characterization.refine.fdr
    int shuffles = 300;          ///< characterization.refine.shuffles
    double grid_da = 0.01;       ///< characterization.refine.grid_da
    double max_shift_da = 700.0; ///< characterization.refine.max_shift_da  (grid runs -max..+max)
    int min_ions = 3;            ///< characterization.refine.min_ions      (smallest cut tried)
    double tolerance_ppm = 10.0; ///< level(2).tolerance_ppm
    int max_levels = 20;         ///< characterization.refine.max_levels    (the cap on the kept levels)
    int jump_penalty = 2;        ///< characterization.refine.jump_penalty  (ions a jump costs: S = ions - penalty x jumps)
    static constexpr uint64_t kSeed = 12345ULL; ///< the shuffle stream's seed (the prototype's)
    static constexpr double kSame = 0.004;      ///< two shifts within this are one (prototype SAME)
  };
} // namespace OpenMS
