// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <cstdint>
#include <string>

namespace OpenMS
{
  /**
   * @brief Carries acquisition and exploration metadata for one deconvolved spectrum.
   *
   * Populated by the exploration engine (Phase 7). Zero overhead when not populated
   * because DeconvolvedSpectrum stores this in std::optional.
   *
   * Serialized to mzML via MSSpectrum::setMetaValue() in DeconvolvedSpectrum::toSpectrum().
   */
  struct OptimizationMetadata
  {
    int group_id = 0;
    int variant_index = -1;
    int total_variants = 0;
    bool is_best_variant = false;
    int rank = 0;
    int msn_level_optimized = 0;
    int exploration_metric = 0;  ///< ExplorationMetric enum cast to int
    int parent_tracking_id = 0;
    double collision_energy = 0;
    double isolation_width = 0;
    std::string activation_type;
    double precursor_mass = 0;
    int precursor_charge = 0;
    double fragmentation_quality_score = -1;
    float tic_coverage = 0;
    int fragment_count = 0;
    uint64_t start_ms = 0;
    uint64_t complete_ms = 0;
    int exploration_scans = 0;
  };
} // namespace OpenMS
