// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Müller $
// $Authors: Tom Müller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/LearnedRanker.h>

#include <array>
#include <cmath>
#include <stdexcept>

namespace OpenMS
{
  LearnedRanker::LearnedRanker(const std::string& model_path)
    : model_path_(model_path)
  {
    try
    {
      env_ = std::make_unique<Ort::Env>(ORT_LOGGING_LEVEL_WARNING, "LearnedRanker");
      Ort::SessionOptions opts;
      opts.SetIntraOpNumThreads(1);
      opts.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
      session_ = std::make_unique<Ort::Session>(*env_, model_path.c_str(), opts);
    }
    catch (const Ort::Exception& e)
    {
      throw std::runtime_error(
        std::string("LearnedRanker: failed to load ONNX model '") + model_path + "': " + e.what());
    }
  }

  std::vector<float> LearnedRanker::score(const std::vector<PeakGroup>& candidates,
                                          const LearnedRankerGlobalContext& ctx) const
  {
    const int N = static_cast<int>(candidates.size());
    if (N == 0) return {};

    std::vector<float> cand_data;
    cand_data.reserve(static_cast<size_t>(N) * CANDIDATE_FEATURE_DIM);
    buildCandidateTensor_(candidates, cand_data);

    std::array<float, GLOBAL_FEATURE_DIM> global_data{
      ctx.rt, ctx.faims_cv, ctx.queue_depth, ctx.elapsed_frac,
      ctx.recent_id_density, ctx.time_since_last_id,
      ctx.mass_exclusion_map_size, ctx.targets_remaining};

    std::vector<uint8_t> mask_data(static_cast<size_t>(N), 1);

    std::array<int64_t, 3> cand_shape{1, N, CANDIDATE_FEATURE_DIM};
    std::array<int64_t, 2> global_shape{1, GLOBAL_FEATURE_DIM};
    std::array<int64_t, 2> mask_shape{1, N};

    auto mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);

    Ort::Value cand_tensor = Ort::Value::CreateTensor<float>(
      mem, cand_data.data(), cand_data.size(), cand_shape.data(), cand_shape.size());
    Ort::Value global_tensor = Ort::Value::CreateTensor<float>(
      mem, global_data.data(), global_data.size(), global_shape.data(), global_shape.size());
    Ort::Value mask_tensor = Ort::Value::CreateTensor<bool>(
      mem, reinterpret_cast<bool*>(mask_data.data()), mask_data.size(),
      mask_shape.data(), mask_shape.size());

    const char* input_names[]  = {"candidates", "global_ctx", "mask"};
    const char* output_names[] = {"scores", "values"};
    std::array<Ort::Value, 3> inputs{std::move(cand_tensor), std::move(global_tensor), std::move(mask_tensor)};

    std::vector<Ort::Value> outputs;
    try
    {
      outputs = session_->Run(Ort::RunOptions{nullptr}, input_names, inputs.data(), 3,
                              output_names, 2);
    }
    catch (const Ort::Exception& e)
    {
      throw std::runtime_error(std::string("LearnedRanker: inference failed: ") + e.what());
    }

    float* scores_ptr = outputs[0].GetTensorMutableData<float>();
    std::vector<float> result(scores_ptr, scores_ptr + N);
    return result;
  }

  void LearnedRanker::buildCandidateTensor_(const std::vector<PeakGroup>& candidates,
                                            std::vector<float>& out) const
  {
    // Feature order MUST match simulation_project/src/sim/feature_schema.json v1.0.0.
    // Per-candidate layout (16 floats):
    //   0: mono_mass
    //   1: log_intensity
    //   2: isotope_cosine         — scan-level isotope pattern cosine
    //   3: charge_cosine          — charge-specific isotope cosine (getChargeIsotopeCosine)
    //   4: snr
    //   5: charge_score
    //   6: qscore
    //   7: idscore_best           — not yet exposed on PeakGroup; Plan 5B placeholder (0.0)
    //   8: ppm_error
    //   9: abs_charge
    //  10: charge_range_width
    //  11: target_decoy
    //  12: is_targeted
    //  13: mono_mass_nearest_excluded_ppm       — Plan 5B placeholder (0.0)
    //  14: rt_since_last_match_for_nearest_mass — Plan 5B placeholder (0.0)
    //  15: priority_hint                        — Plan 5B placeholder (0.0)
    out.clear();
    for (const auto& pg : candidates)
    {
      const int rep_charge = pg.getRepAbsCharge();

      const float mono_mass    = static_cast<float>(pg.getMonoMass());
      const float intensity    = static_cast<float>(pg.getIntensity());
      const float log_intensity = (intensity > 0.0f) ? std::log(intensity) : 0.0f;
      const float iso_cos      = pg.getIsotopeCosine();
      const float ch_iso_cos   = pg.getChargeIsotopeCosine(rep_charge);
      const float snr          = pg.getSNR();
      const float ch_score     = pg.getChargeScore();
      const float qscore       = static_cast<float>(pg.getQscore());
      const float idscore_best = 0.0f;  // not yet exposed on PeakGroup; Plan 5B
      const float ppm_err      = pg.getAvgPPMError();
      const float abs_charge   = static_cast<float>(rep_charge);

      // getAbsChargeRange() returns std::tuple<int, int>
      const auto charge_range       = pg.getAbsChargeRange();
      const int  cmin               = std::get<0>(charge_range);
      const int  cmax               = std::get<1>(charge_range);
      const float charge_range_width = static_cast<float>(cmax - cmin);

      // getTargetDecoyType() returns PeakGroup::TargetDecoyType (unscoped enum)
      const float target_decoy = static_cast<float>(static_cast<int>(pg.getTargetDecoyType()));
      const float is_targeted  = pg.isTargeted() ? 1.0f : 0.0f;

      // Exclusion-map state — Plan 5B will wire these through:
      const float mono_mass_nearest_excluded_ppm       = 0.0f;
      const float rt_since_last_match_for_nearest_mass = 0.0f;
      const float priority_hint                        = 0.0f;

      out.push_back(mono_mass);
      out.push_back(log_intensity);
      out.push_back(iso_cos);
      out.push_back(ch_iso_cos);
      out.push_back(snr);
      out.push_back(ch_score);
      out.push_back(qscore);
      out.push_back(idscore_best);
      out.push_back(ppm_err);
      out.push_back(abs_charge);
      out.push_back(charge_range_width);
      out.push_back(target_decoy);
      out.push_back(is_targeted);
      out.push_back(mono_mass_nearest_excluded_ppm);
      out.push_back(rt_since_last_match_for_nearest_mass);
      out.push_back(priority_hint);
    }
  }
} // namespace OpenMS
