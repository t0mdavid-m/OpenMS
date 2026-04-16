// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Exploration.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>

namespace OpenMS
{

  Exploration::Exploration(const Config& config, FragmentAnalysis& fragments)
    : config_(config), fragments_(fragments)
  {
    DoubleList expl_tol;
    for (const auto& [lvl, cfg] : config.levels())
      expl_tol.push_back(cfg.exploration_tolerance_ppm);
    exploration_deconv_ = std::make_unique<Deconvolution>(config, expl_tol);
  }

  std::vector<Exploration::VariantParams> Exploration::buildVariants_(
      const MSLevelConfig& cfg, const ScanConfig& base_config) const
  {
    std::vector<VariantParams> variants;

    // Determine which activations to sweep
    std::vector<std::string> acts = cfg.activations;
    if (acts.empty())
      acts.push_back(base_config.activation);

    // Build CE and RT value arrays
    std::vector<double> ces;
    for (double ce = cfg.ce_min; ce <= cfg.ce_max + 1e-9; ce += cfg.ce_step)
      ces.push_back(ce);

    std::vector<double> rts;
    if (cfg.rt_max > cfg.rt_min)
    {
      for (double rt = cfg.rt_min; rt <= cfg.rt_max + 1e-9; rt += cfg.rt_step)
        rts.push_back(rt);
    }

    for (const auto& act : acts)
    {
      bool sweep_ce = (act == "HCD" || act == "CID" || act == "EThcD");
      bool sweep_rt = (act == "ETD" || act == "EThcD");

      if (sweep_ce && sweep_rt)
      {
        // Cross-product (EThcD)
        for (double ce : ces)
          for (double rt : rts)
            variants.push_back({act, ce, rt});
      }
      else if (sweep_ce)
      {
        // CE-only (HCD/CID)
        for (double ce : ces)
          variants.push_back({act, ce, base_config.reaction_time});
      }
      else if (sweep_rt)
      {
        // RT-only (ETD)
        for (double rt : rts)
          variants.push_back({act, static_cast<double>(base_config.collision_energy), rt});
      }
      else
      {
        // Unknown activation — single variant with base config values
        variants.push_back({act, static_cast<double>(base_config.collision_energy), base_config.reaction_time});
      }
    }

    return variants;
  }

  std::vector<ScanCommand> Exploration::initiate(int msn_level, const PeakGroup& pg, int charge,
      double faims_cv, ScanCommandQueue& queue, const ScanCommand* ms_ctx,
      char ion_type, int frag_index,
      const MS3FragmentMatcher::ProteoformContext& proto_ctx)
  {
    std::vector<ScanCommand> commands;

    const auto& cfg = config_.level(msn_level);
    if (!config_.hasExploration(msn_level)) return commands;

    // Build base ScanConfig from the level's primary scan config, then apply overrides
    ScanConfig base_config = cfg.scans[0];
    base_config.applyOverrides(cfg.overrides);

    auto variant_params = buildVariants_(cfg, base_config);
    if (variant_params.empty()) return commands;

    // Prepend baseline variant for RemainingPrecursor metric
    bool needs_baseline = (cfg.exploration == ExplorationMetric::RemainingPrecursor);
    if (needs_baseline)
      variant_params.insert(variant_params.begin(),
                            {base_config.activation, 0.0, base_config.reaction_time});

    // Compute precursor_mz and isolation_width from PeakGroup
    auto [mz1, mz2] = pg.getMzRange(charge);
    double precursor_mz = (mz1 + mz2) / 2.0;
    double isolation_width = std::max(mz2 - mz1, 2.0);
    double precursor_mass = pg.getMonoMass();

    ExplorationGroup group;
    group.group_id = next_group_id_++;
    group.msn_level = msn_level;
    group.exploration_metric = cfg.exploration;
    group.precursor_mz = precursor_mz;
    group.precursor_mass = precursor_mass;
    group.isolation_width = isolation_width;
    group.precursor_charge = charge;
    group.precursor_pg = pg;
    group.faims_cv = faims_cv;
    group.start_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    // Capture originating MS2 context for buildMS3 (stage 0)
    if (ms_ctx != nullptr) group.originating_cmd = *ms_ctx;
    group.fragment_ion_type = ion_type;
    group.fragment_ion_index = frag_index;
    group.proteoform_ctx = proto_ctx;

    for (int i = 0; i < static_cast<int>(variant_params.size()); ++i)
    {
      const auto& vp = variant_params[i];
      ExplorationVariant v;
      v.variant_index = (needs_baseline && i == 0) ? -1 : (needs_baseline ? i - 1 : i);
      v.collision_energy = vp.collision_energy;
      v.reaction_time = vp.reaction_time;
      v.is_baseline = (needs_baseline && i == 0);
      v.activation_type = vp.activation;

      ScanConfig variant_config = base_config;
      variant_config.collision_energy = static_cast<int>(vp.collision_energy);
      variant_config.reaction_time = vp.reaction_time;
      variant_config.activation = vp.activation;

      int expl_priority = (msn_level >= 3) ? 1 : 2;
      ScanCommand cmd;
      if (msn_level >= 3 && ms_ctx != nullptr)
      {
        cmd = queue.buildMS3(*ms_ctx, variant_config,
                             precursor_mz, charge, isolation_width,
                             ion_type, frag_index, expl_priority);
      }
      else
      {
        cmd = queue.buildMS2(pg, charge, variant_config, expl_priority);
      }
      cmd.faims_cv = faims_cv;

      int id_int = cmd.scan_id;
      std::string id_str = ScanCommandQueue::encode(id_int);
      v.tracking_id = id_str;
      v.cmd = cmd;
      if (msn_level >= 3 && ion_type != '\0' && frag_index > 0)
      {
        double frag_mass_kda = precursor_mz * charge / 1000.0;
        std::snprintf(cmd.scan_description, 16, "%sE%.1f@%d%c%d",
                     id_str.c_str(), frag_mass_kda, charge, ion_type, frag_index);
      }
      else
      {
        std::snprintf(cmd.scan_description, 16, "%sE%.1f@%d",
                     id_str.c_str(), precursor_mass / 1000.0, charge);
      }

      group.variants.push_back(v);
      variant_tracking_map_[id_int] = {group.group_id, i};
      commands.push_back(cmd);

      std::cout << "[TRACK-CREATE] id=" << id_str
                << " ms_level=" << msn_level << " type=exploration"
                << " activation=" << vp.activation
                << " CE=" << vp.collision_energy
                << " RT=" << vp.reaction_time << std::endl;
    }

    active_groups_[group.group_id] = std::move(group);
    return commands;
  }

  bool Exploration::isExplorationVariant(int tracking_id) const
  {
    return variant_tracking_map_.find(tracking_id) != variant_tracking_map_.end();
  }

  Exploration::FeedResultInfo Exploration::feedResult(int tracking_id,
      const double* mzs, const double* ints, int length,
      double rt, ScanCommandQueue& queue)
  {
    // Look up the group to get correct precursor context for deconvolution
    auto vit = variant_tracking_map_.find(tracking_id);
    if (vit == variant_tracking_map_.end()) return {};

    int group_id = vit->second.group_id;
    auto git = active_groups_.find(group_id);
    if (git == active_groups_.end()) return {};
    const ExplorationGroup& group = git->second;

    // Deconvolve with correct precursor context from the group
    DeconvolvedSpectrum ms2_deconv(tracking_id);
    if (mzs != nullptr && ints != nullptr && length > 0)
    {
      exploration_deconv_->deconvolveMSn(mzs, ints, length, rt,
                                         group.precursor_mass, group.precursor_charge);
      ms2_deconv = exploration_deconv_->storedMS2();
    }

    return feedResultImpl_(tracking_id, ms2_deconv, mzs, ints, length, rt, queue);
  }

  Exploration::FeedResultInfo Exploration::feedResultForTest(int tracking_id,
      const DeconvolvedSpectrum& ms2_deconv, double rt, ScanCommandQueue& queue)
  {
    return feedResultImpl_(tracking_id, ms2_deconv, nullptr, nullptr, 0, rt, queue);
  }

  Exploration::FeedResultInfo Exploration::feedResultImpl_(int tracking_id,
      const DeconvolvedSpectrum& ms2_deconv,
      const double* mzs, const double* ints, int length,
      double rt, ScanCommandQueue& queue)
  {
    (void)rt;
    FeedResultInfo info;

    // Look up the variant reference
    auto vit = variant_tracking_map_.find(tracking_id);
    if (vit == variant_tracking_map_.end()) return info;

    int group_id = vit->second.group_id;
    int variant_index = vit->second.variant_index;
    variant_tracking_map_.erase(vit);

    // Find the group
    auto git = active_groups_.find(group_id);
    if (git == active_groups_.end()) return info;
    ExplorationGroup& group = git->second;

    if (variant_index < 0 || variant_index >= static_cast<int>(group.variants.size())) return info;
    ExplorationVariant& v = group.variants[variant_index];
    if (v.received) return info;

    v.result = ms2_deconv;
    double remaining_ratio = -1.0;
    FragmentMatchResult frag{};
    v.score = computeExplorationScore_(group.exploration_metric, ms2_deconv, group, mzs, ints, length, &remaining_ratio, &frag, v.activation_type);
    v.tic_coverage = computeTICCoverage_(ms2_deconv);
    v.fragment_count = static_cast<int>(frag.count);
    v.received = true;

    // Handle baseline for RemainingPrecursor
    if (v.is_baseline)
    {
      double iso_half = group.isolation_width / 2.0;
      double mz_low = group.precursor_mz - iso_half;
      double mz_high = group.precursor_mz + iso_half;
      double baseline_sum = 0.0;
      if (mzs != nullptr && ints != nullptr)
      {
        for (int bi = 0; bi < length; ++bi)
        {
          if (mzs[bi] >= mz_low && mzs[bi] <= mz_high)
            baseline_sum += ints[bi];
        }
      }
      group.baseline_intensity = baseline_sum;
      group.has_baseline = true;
      v.score = 0.0;  // baseline score not meaningful
    }

    // Count real (non-baseline) variants for metadata
    int real_variant_count = 0;
    for (const auto& vr : group.variants)
      if (!vr.is_baseline) ++real_variant_count;

    // Populate per-variant metadata in the return info
    info.group_id = group.group_id;
    info.variant_index = v.variant_index;
    info.total_variants = real_variant_count;
    info.collision_energy = v.collision_energy;
    info.activation_type = v.activation_type;
    info.reaction_time = v.reaction_time;
    info.score = v.score;
    info.tic_coverage = v.tic_coverage;
    info.fragment_count = v.fragment_count;
    info.exploration_metric = static_cast<int>(group.exploration_metric);
    info.matched_protein = frag.matched_protein;
    info.proteoform_sequence = frag.proteoform_sequence;
    info.remaining_ratio = remaining_ratio;
    std::string parent_enc = ScanCommandQueue::encode(group.originating_cmd.scan_id);
    std::strncpy(info.parent_scan_id, parent_enc.c_str(), 3);
    info.parent_scan_id[3] = '\0';

    // Populate identification context from the group
    info.ms2_context.proteoform_sequence = (group.proteoform_ctx.region_start >= 0)
        ? config_.targeting().protein_sequence.substr(group.proteoform_ctx.region_start,
            group.proteoform_ctx.region_end - group.proteoform_ctx.region_start)
        : "";
    info.ms2_context.start_pos = group.proteoform_ctx.region_start;
    info.ms2_context.end_pos = group.proteoform_ctx.region_end;
    info.ms2_context.ptm_sites = group.proteoform_ctx.ptm_sites;
    info.ms2_context.fragment_ion_type = group.fragment_ion_type;
    info.ms2_context.fragment_ion_index = group.fragment_ion_index;
    if (group.originating_cmd.num_stages > 0)
    {
      info.ms2_context.ms1_precursor_mass = group.originating_cmd.mono_mass;
      info.ms2_context.ms1_precursor_mz = group.originating_cmd.stages[0].precursor_mz;
      info.ms2_context.ms1_precursor_charge = group.originating_cmd.stages[0].charge_state;
    }
    const auto& vcmd = group.variants[variant_index].cmd;
    if (vcmd.num_stages >= 2)
    {
      info.ms2_context.fragment_mz = vcmd.stages[1].precursor_mz;
      info.ms2_context.fragment_charge = vcmd.stages[1].charge_state;
    }
    info.ms2_context.fragment_mass = group.precursor_mass;
    info.identification_result = group.variants[variant_index].identification_result;

    auto& meta = v.result.getOrCreateOptimizationMetadata();
    meta.group_id = group.group_id;
    meta.variant_index = v.variant_index;
    meta.total_variants = real_variant_count;
    meta.is_best_variant = false;
    meta.msn_level_optimized = group.msn_level;
    meta.exploration_metric = static_cast<int>(group.exploration_metric);
    meta.collision_energy = v.collision_energy;
    meta.activation_type = v.activation_type;
    meta.precursor_mass = group.precursor_mass;
    meta.precursor_charge = group.precursor_charge;
    meta.fragmentation_quality_score = v.score;
    meta.tic_coverage = v.tic_coverage;
    meta.fragment_count = v.fragment_count;
    meta.start_ms = group.start_ms;
    meta.complete_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());
    meta.exploration_scans = real_variant_count;

    bool all_received = std::all_of(group.variants.begin(), group.variants.end(),
                                    [](const ExplorationVariant& x){ return x.received; });
    if (!all_received) return info;

    // MS3 FragmentCount: batch re-score with calibrated subsequence matching
    if (group.exploration_metric == ExplorationMetric::FragmentCount && group.msn_level >= 3)
    {
      std::vector<const DeconvolvedSpectrum*> variant_spectra;
      for (auto& var : group.variants)
        variant_spectra.push_back(var.received ? &var.result : nullptr);

      std::vector<MS3FragmentMatcher::MatchResult> detailed_results;
      auto calibrated_scores = MS3FragmentMatcher::calibrateAndScore(
        variant_spectra,
        config_.targeting().protein_sequence,
        group.proteoform_ctx,
        group.fragment_ion_type,
        group.fragment_ion_index,
        MS3FragmentMatcher::LOOSE_TOLERANCE_PPM,
        config_.level(group.msn_level).tolerance_ppm,
        &detailed_results);

      for (size_t vi = 0; vi < calibrated_scores.size(); ++vi)
      {
        group.variants[vi].score = calibrated_scores[vi];
        group.variants[vi].fragment_count = static_cast<int>(calibrated_scores[vi]);
        if (vi < detailed_results.size())
          group.variants[vi].identification_result = detailed_results[vi];
      }

      // Re-populate identification_result after batch re-scoring updated it
      info.identification_result = group.variants[variant_index].identification_result;
    }

    int best_idx = -1;
    double best_score = -1.0;
    for (int i = 0; i < static_cast<int>(group.variants.size()); ++i)
    {
      if (group.variants[i].is_baseline) continue;  // skip baseline
      if (group.variants[i].score > best_score)
      {
        best_score = group.variants[i].score;
        best_idx = i;
      }
    }
    if (best_idx < 0) return info;  // shouldn't happen
    group.winner_index = best_idx;
    group.complete = true;
    group.variants[best_idx].result.getOrCreateOptimizationMetadata().is_best_variant = true;

    std::cout << "[EXPL-WINNER] group=" << group.group_id
              << " winner_idx=" << best_idx
              << " activation=" << group.variants[best_idx].activation_type
              << " CE=" << group.variants[best_idx].collision_energy
              << " RT=" << group.variants[best_idx].reaction_time
              << " score=" << best_score << std::endl;

    const auto& level_config = config_.level(group.msn_level);
    if (!level_config.overrides.empty())
    {
      ScanConfig prod_config = level_config.scans[0];
      prod_config.collision_energy = static_cast<int>(group.variants[best_idx].collision_energy);
      prod_config.reaction_time = group.variants[best_idx].reaction_time;
      prod_config.activation = group.variants[best_idx].activation_type;

      ScanCommand prod_cmd;
      if (group.msn_level >= 3)
      {
        prod_cmd = queue.buildMS3(group.variants[best_idx].cmd, prod_config,
                                   group.precursor_mz, group.precursor_charge,
                                   group.isolation_width,
                                   group.fragment_ion_type, group.fragment_ion_index, 1);
      }
      else
      {
        prod_cmd = queue.buildMS2(group.precursor_pg, group.precursor_charge, prod_config, 2);
      }
      prod_cmd.faims_cv = group.faims_cv;

      std::string prod_id = ScanCommandQueue::encode(prod_cmd.scan_id);
      std::cout << "[TRACK-CREATE] id=" << prod_id
                << " ms_level=" << group.msn_level << " type=production"
                << std::endl;

      info.commands.push_back(prod_cmd);
    }
    else if (group.msn_level < 3)
    {
      std::cout << "exploration call site" << std::endl;
      auto next_nlr = initiateNextLevel(group.msn_level,
          group.variants[best_idx].result, group.faims_cv, queue,
          &group.variants[best_idx].cmd);
      info.commands.insert(info.commands.end(), next_nlr.commands.begin(), next_nlr.commands.end());
    }

    for (const auto& vr : group.variants)
      variant_tracking_map_.erase(queue.decode(vr.tracking_id));

    active_groups_.erase(git);
    return info;
  }

  Exploration::NextLevelResult Exploration::initiateNextLevel(int msn_level,
      const DeconvolvedSpectrum& result, double faims_cv, ScanCommandQueue& queue,
      const ScanCommand* ms_ctx)
  {
    NextLevelResult nlr;

    int next_level = msn_level + 1;
    const auto& this_cfg = config_.level(msn_level);
    const auto& next_cfg = config_.level(next_level);
    if (this_cfg.selection == SelectionMetric::None) return nlr;

    const auto& seq = config_.targeting().protein_sequence;
    int num_targets = this_cfg.max_targets;

    // Get fragment ion targets via FragmentAnalysis
    DeconvolvedSpectrum result_copy = result;
    std::vector<double> masses(num_targets), qscores(num_targets);
    std::vector<double> wstarts(num_targets), wends(num_targets);
    std::vector<int> charges(num_targets);
    std::vector<char> ion_types(num_targets, '\0');
    std::vector<int> frag_indices(num_targets, 0);

    int found = 0;

    switch (next_cfg.selection)
    {
      case SelectionMetric::Intensity:
      case SelectionMetric::QScore:
        found = fragments_.getTopFragmentMatches(seq, num_targets,
            masses.data(), qscores.data(), charges.data(),
            wstarts.data(), wends.data(),
            ion_types.data(), frag_indices.data(), result_copy);
        break;
      case SelectionMetric::TerminalFragments:
        found = fragments_.getTerminalFragmentIons(seq, num_targets,
            masses.data(), qscores.data(), charges.data(),
            wstarts.data(), wends.data(),
            ion_types.data(), frag_indices.data(), result_copy);
        break;
      case SelectionMetric::AmbiguityResolution:
        found = fragments_.getAmbiguityEnclosingIons(seq, num_targets,
            masses.data(), qscores.data(), charges.data(),
            wstarts.data(), wends.data(),
            ion_types.data(), frag_indices.data(), result_copy);
        break;
      default:
        break;
    }

    // Cache proteoform context for MS3 subsequence scoring
    MS3FragmentMatcher::ProteoformContext proto_ctx;
    if (next_level >= 3)
    {
      const auto& pinfo = fragments_.getLastProteoformInfo();
      proto_ctx.region_start = pinfo.region_start;
      proto_ctx.region_end = pinfo.region_end;
      proto_ctx.ptm_sites = pinfo.ptm_sites;
      // If no truncation detected, use full protein sequence bounds
      if (proto_ctx.region_start < 0)
        proto_ctx.region_start = 0;
      if (proto_ctx.region_end < 0)
        proto_ctx.region_end = static_cast<int>(seq.size());
    }

    // Populate fragment matching metadata
    nlr.fragment_count = found;
    if (!seq.empty() && found > 0)
    {
      nlr.matched_protein = config_.targeting().fasta_file;
      nlr.proteoform_sequence = seq;
      // TIC coverage: sum of matched qscores / total spectrum intensity
      double total_tic = 0.0;
      for (Size i = 0; i < result.size(); ++i)
        total_tic += result[i].getIntensity();
      double matched_intensity = 0.0;
      for (int i = 0; i < found; ++i)
        matched_intensity += qscores[i];
      if (total_tic > 0.0)
        nlr.tic_coverage = static_cast<float>(matched_intensity / total_tic);
    }

    num_targets = std::min(num_targets, found);

    // Build commands for each selected fragment target
    ScanConfig next_scan_config = next_cfg.scans[0];

    int charge_floor = this_cfg.min_charge;

    if (config_.hasExploration(next_level))
    {
      // Recursive exploration at next level
      for (int ti = 0; ti < num_targets; ++ti)
      {
        int abs_charge = std::abs(charges[ti]);
        if (charge_floor > 0 && abs_charge < charge_floor)
          continue;
        PeakGroup frag_pg(abs_charge, abs_charge, true);
        frag_pg.setMonoisotopicMass(masses[ti]);
        FLASHHelperClasses::LogMzPeak lp_lo;
        lp_lo.mz = wstarts[ti];
        lp_lo.abs_charge = abs_charge;
        frag_pg.push_back(lp_lo);
        FLASHHelperClasses::LogMzPeak lp_hi;
        lp_hi.mz = wends[ti];
        lp_hi.abs_charge = abs_charge;
        frag_pg.push_back(lp_hi);

        auto sub_cmds = initiate(next_level, frag_pg, std::abs(charges[ti]), faims_cv, queue, ms_ctx,
                                 ion_types[ti], frag_indices[ti], proto_ctx);
        nlr.commands.insert(nlr.commands.end(), sub_cmds.begin(), sub_cmds.end());

        for (size_t sci = 0; sci < sub_cmds.size(); ++sci)
        {
          MS2Context mc;
          mc.proteoform_sequence = nlr.proteoform_sequence;
          mc.start_pos = proto_ctx.region_start;
          mc.end_pos = proto_ctx.region_end;
          mc.ptm_sites = proto_ctx.ptm_sites;
          if (ms_ctx != nullptr)
          {
            mc.ms1_precursor_mass = ms_ctx->mono_mass;
            mc.ms1_precursor_mz = ms_ctx->stages[0].precursor_mz;
            mc.ms1_precursor_charge = ms_ctx->stages[0].charge_state;
          }
          mc.fragment_ion_type = ion_types[ti];
          mc.fragment_ion_index = frag_indices[ti];
          mc.fragment_mass = masses[ti];
          mc.fragment_mz = (wstarts[ti] + wends[ti]) / 2.0;
          mc.fragment_charge = std::abs(charges[ti]);
          nlr.ms3_contexts.push_back(mc);
        }
      }
    }
    else
    {
      // Direct command building for each fragment target
      for (int ti = 0; ti < num_targets; ++ti)
      {
        double frag_mz = (wstarts[ti] + wends[ti]) / 2.0;
        int frag_charge = std::abs(charges[ti]);
        if (charge_floor > 0 && frag_charge < charge_floor)
          continue;
        double iso_width = wends[ti] - wstarts[ti];

        ScanCommand cmd;
        if (next_level >= 3 && ms_ctx != nullptr)
        {
          // MS3: proper two-stage command via buildMS3 with config CE/activation
          cmd = queue.buildMS3(*ms_ctx, next_scan_config, frag_mz, frag_charge, iso_width,
                               ion_types[ti], frag_indices[ti], 1);
        }
        else
        {
          // MS2: single-stage command via buildMS2
          PeakGroup frag_pg(frag_charge, frag_charge, true);
          frag_pg.setMonoisotopicMass(masses[ti]);
          FLASHHelperClasses::LogMzPeak lp;
          lp.mz = frag_mz;
          lp.abs_charge = frag_charge;
          frag_pg.push_back(lp);

          cmd = queue.buildMS2(frag_pg, frag_charge, next_scan_config, 1);
          cmd.msn_level = next_level;
        }
        cmd.faims_cv = faims_cv;

        std::string id_str = ScanCommandQueue::encode(cmd.scan_id);
        std::cout << "[TRACK-CREATE] id=" << id_str
                  << " ms_level=" << next_level << " type=next_level"
                  << std::endl;

        nlr.commands.push_back(cmd);

        if (next_level >= 3 && ms_ctx != nullptr)
        {
          MS2Context mc;
          mc.proteoform_sequence = nlr.proteoform_sequence;
          mc.start_pos = proto_ctx.region_start;
          mc.end_pos = proto_ctx.region_end;
          mc.ptm_sites = proto_ctx.ptm_sites;
          mc.ms1_precursor_mass = ms_ctx->mono_mass;
          mc.ms1_precursor_mz = ms_ctx->stages[0].precursor_mz;
          mc.ms1_precursor_charge = ms_ctx->stages[0].charge_state;
          mc.fragment_ion_type = ion_types[ti];
          mc.fragment_ion_index = frag_indices[ti];
          mc.fragment_mass = masses[ti];
          mc.fragment_mz = frag_mz;
          mc.fragment_charge = frag_charge;
          nlr.ms3_contexts.push_back(mc);
        }
      }
    }

    return nlr;
  }

  int Exploration::activeGroupCount() const
  {
    return static_cast<int>(active_groups_.size());
  }

  Exploration::ExplorationGroup Exploration::getGroup(int group_id) const
  {
    return active_groups_.at(group_id);
  }

  double Exploration::computeExplorationScore_(ExplorationMetric metric,
      const DeconvolvedSpectrum& spec,
      const ExplorationGroup& group,
      const double* mzs, const double* ints, int length,
      double* out_remaining_ratio, FragmentMatchResult* out_frag,
      const std::string& activation_type) const
  {
    FragmentMatchResult fmr;
    switch (metric)
    {
      case ExplorationMetric::MassCount:
        fmr = computeFragmentMatch_(spec, group.msn_level, activation_type);
        if (out_frag) *out_frag = fmr;
        return computeMassCount_(spec);
      case ExplorationMetric::RemainingPrecursor:
        fmr = computeFragmentMatch_(spec, group.msn_level, activation_type);
        if (out_frag) *out_frag = fmr;
        return computeRemainingPrecursorScore_(group, mzs, ints, length, out_remaining_ratio);
      case ExplorationMetric::FragmentCount:
      {
        fmr = computeFragmentMatch_(spec, group.msn_level, activation_type);
        if (out_frag) *out_frag = fmr;
        return fmr.count;
      }
      default:
        fmr = computeFragmentMatch_(spec, group.msn_level, activation_type);
        if (out_frag) *out_frag = fmr;
        return computeMassCount_(spec);
    }
  }

  double Exploration::computeMassCount_(const DeconvolvedSpectrum& spec) const
  {
    return static_cast<double>(spec.size());
  }

  double Exploration::computeRemainingPrecursorScore_(const ExplorationGroup& group,
      const double* mzs, const double* ints, int length, double* out_ratio) const
  {
    if (out_ratio) *out_ratio = -1.0;

    if (length <= 0 || mzs == nullptr || ints == nullptr)
      return 0.0;

    double iso_half = group.isolation_width / 2.0;
    double mz_low = group.precursor_mz - iso_half;
    double mz_high = group.precursor_mz + iso_half;

    double remaining_intensity = 0.0;
    for (int i = 0; i < length; ++i)
    {
      if (mzs[i] >= mz_low && mzs[i] <= mz_high)
        remaining_intensity += ints[i];
    }

    // Reference: baseline isolation-window intensity (CE=0 scan)
    double reference;
    if (group.has_baseline)
    {
      reference = group.baseline_intensity;
      if (reference <= 0.0)
        return -1.0;  // Baseline failed: ratio = 1/1, score = 0
    }
    else
    {
      return -2.0;  // Baseline not yet received: ratio = 1/1, score = 0
    }

    double ratio = remaining_intensity / reference;
    if (out_ratio) *out_ratio = ratio;
    double target = config_.level(group.msn_level).remaining_precursor_target;
    double deviation = std::abs(ratio - target);
    double score = 1.0 - deviation;
    if (score < 0.0) score = -3.0;
    return score;
  }

  Exploration::FragmentMatchResult Exploration::computeFragmentMatch_(const DeconvolvedSpectrum& spec, int msn_level,
                                                                      const std::string& activation_type) const
  {
    FragmentMatchResult result;
    const auto& seq = config_.targeting().protein_sequence;
    if (seq.empty() || spec.empty())
      return result;

    DeconvolvedSpectrum spec_copy = spec;

    const int max_matches = 100;
    std::vector<double> masses(max_matches), qscores(max_matches);
    std::vector<double> wstarts(max_matches), wends(max_matches);
    std::vector<int> charges(max_matches);
    std::vector<char> ion_types(max_matches, '\0');
    std::vector<int> frag_indices(max_matches, 0);

    int count = fragments_.getTopFragmentMatches(
        seq, max_matches,
        masses.data(), qscores.data(), charges.data(),
        wstarts.data(), wends.data(),
        ion_types.data(), frag_indices.data(),
        spec_copy, activation_type,
        config_.level(msn_level).exploration_tolerance_ppm);

    result.count = static_cast<double>(count);
    if (count > 0)
    {
      result.matched_protein = config_.targeting().fasta_file;
      result.proteoform_sequence = seq;
    }
    return result;
  }

  float Exploration::computeTICCoverage_(const DeconvolvedSpectrum& spec) const
  {
    if (spec.empty()) return 0.0f;
    float total = 0.0f;
    for (Size i = 0; i < spec.size(); ++i)
      total += static_cast<float>(spec[i].getIntensity());
    return total > 0.0f ? 1.0f : 0.0f;
  }

} // namespace OpenMS
