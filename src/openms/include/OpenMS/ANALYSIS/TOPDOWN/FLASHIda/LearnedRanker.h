// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Müller $
// $Authors: Tom Müller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <onnxruntime_cxx_api.h>

#include <memory>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
   * @brief Per-scan global context feeding the LearnedRanker alongside PeakGroups.
   *
   * Mirrors the 8-float "global" slice of simulation_project/src/sim/feature_schema.json.
   * Fields must match that schema in count and order.
   */
  struct OPENMS_DLLAPI LearnedRankerGlobalContext
  {
    float rt                      = 0.0f;
    float faims_cv                = 0.0f;
    float queue_depth             = 0.0f;
    float elapsed_frac            = 0.0f;
    float recent_id_density       = 0.0f;
    float time_since_last_id      = 0.0f;
    float mass_exclusion_map_size = 0.0f;
    float targets_remaining       = 0.0f;
  };

  /**
   * @brief ONNX-backed learned ranker for FLASHIda precursor selection.
   *
   * Loads a trained policy from disk and scores candidate PeakGroups using the
   * same 16-float per-candidate + 8-float global feature layout the Python
   * training pipeline uses.
   *
   * Fail-loud policy: construction throws if the ONNX model cannot be loaded;
   * scoring throws if onnxruntime raises an exception. No fallback to other
   * SelectionMetrics.
   */
  class OPENMS_DLLAPI LearnedRanker
  {
  public:
    /// Load an ONNX model from @p model_path. Throws std::runtime_error on failure.
    explicit LearnedRanker(const std::string& model_path);

    // Non-copyable, movable (Ort::Session is not copyable)
    LearnedRanker(const LearnedRanker&) = delete;
    LearnedRanker& operator=(const LearnedRanker&) = delete;
    LearnedRanker(LearnedRanker&&) = default;
    LearnedRanker& operator=(LearnedRanker&&) = default;

    ~LearnedRanker() = default;

    /**
     * @brief Score each candidate PeakGroup. Throws on inference failure.
     *
     * @param candidates Deconvolved peak groups to score.
     * @param ctx Global context features for this scan.
     * @return Per-candidate scores (higher = more promising), same size as @p candidates.
     */
    std::vector<float> score(const std::vector<PeakGroup>& candidates,
                             const LearnedRankerGlobalContext& ctx) const;

  private:
    /// Build the (N, 16) candidate feature tensor from the peak-group vector.
    void buildCandidateTensor_(const std::vector<PeakGroup>& candidates,
                               std::vector<float>& out) const;

    std::unique_ptr<Ort::Env> env_;
    mutable std::unique_ptr<Ort::Session> session_;
    std::string model_path_;

    static constexpr int CANDIDATE_FEATURE_DIM = 16;
    static constexpr int GLOBAL_FEATURE_DIM    = 8;
  };
} // namespace OpenMS
