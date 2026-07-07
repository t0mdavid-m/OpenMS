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
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>

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
    exploration_deconv_ = std::make_unique<Deconvolution>(config, config.explorationToleranceList());
  }

  std::vector<Exploration::VariantParams> Exploration::buildVariants_(
      const MSLevelConfig& cfg, const ScanConfig& base_config) const
  {
    std::vector<VariantParams> variants;

    // Determine which activations to sweep
    std::vector<std::string> acts = cfg.activations;
    if (acts.empty())
      acts.push_back(base_config.activation);

    // Build collision energy and reaction time value arrays
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
      ScanCommandQueue& queue, const MSSpectrum* source_spectrum, const ScanCommand* ms_ctx,
      char ion_type, int frag_index,
      const MS3FragmentMatcher::ProteoformContext& proto_ctx,
      const FragmentAnalysis::FragmentScores& frag_scores,
      const Ms2Params* stage0_params)
  {
    std::vector<ScanCommand> commands;

    const auto& cfg = config_.level(msn_level);
    if (!config_.hasExploration(msn_level)) return commands;

    // Build base ScanConfig from the level's primary scan config, then apply overrides
    ScanConfig base_config = cfg.scans[0];
    base_config.applyOverrides(cfg.overrides);

    auto variant_params = buildVariants_(cfg, base_config);
    if (variant_params.empty()) return commands;

    // Prepend a CE-0 baseline (variant_index -1) to EVERY exploration metric. It is skipped in winner
    // selection (is_baseline) and excluded from the ProteoformTracker feed, so it never wins or pools;
    // it fires one extra pre-scan per group and measures the un-fragmented reference. needs_baseline is
    // kept (always true now) because the per-variant index/is_baseline logic below reads it.
    const bool needs_baseline = true;
    variant_params.insert(variant_params.begin(), {base_config.activation, 0.0, 0.0});

    // Compute precursor_mz and isolation_width from PeakGroup. The 2.0 Th minimum applies to MS3+
    // fragment isolation only (buildMS3 stage[1] also re-floors); an MS2 group scores over its
    // natural precursor window, matching the actual (unfloored) commanded MS2 isolation.
    auto [mz1, mz2] = pg.getMzRange(charge);
    double precursor_mz = (mz1 + mz2) / 2.0;
    double isolation_width = (msn_level >= 3) ? std::max(mz2 - mz1, 2.0) : (mz2 - mz1);
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
    group.faims_cv = (ms_ctx != nullptr) ? ms_ctx->faims_cv : 0.0;  // Item 1: CV travels via the context
    group.start_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

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
                             precursor_mz, charge, isolation_width, ms_ctx->scan_id,
                             ion_type, frag_index, expl_priority, frag_scores,
                             stage0_params, proto_ctx);
      }
      else
      {
        cmd = queue.buildMS2(pg, charge, variant_config, expl_priority, ms_ctx ? ms_ctx->scan_id : 0);
      }
      cmd.faims_cv = group.faims_cv;  // CV inferred from the MS1 context (group.faims_cv <- ms_ctx->faims_cv)

      // Stamp window-SNR over the ACTUAL fragmentation window: the last stage (MS2 -> stage[0],
      // MS3 -> stage[1]); the signal is that stage's selected-precursor intensity (stage[1] uses _s1).
      if (source_spectrum != nullptr && cmd.num_stages >= 1)
      {
        const int si = cmd.num_stages - 1;
        const double signal = (cmd.num_stages >= 2) ? cmd.precursor_intensity_s1 : cmd.precursor_intensity;
        const double half = cmd.stages[si].isolation_width / 2.0;
        cmd.window_snr = FragmentAnalysis::windowSnr(*source_spectrum,
            cmd.stages[si].precursor_mz - half, cmd.stages[si].precursor_mz + half, signal);
      }

      int id_int = cmd.scan_id;
      std::string id_str = ScanCommandQueue::encode(id_int);
      v.tracking_id = id_str;
      v.cmd = cmd;
      if (msn_level >= 3 && ion_type != '\0' && frag_index > 0)
      {
        double frag_mass_kda = precursor_mz * charge / 1000.0;
        std::string mass_tok = ScanCommandQueue::formatMassToken(frag_mass_kda, charge, ion_type, frag_index);
        std::snprintf(cmd.scan_description, 16, "%sE%s@%d%c%d",
                     id_str.c_str(), mass_tok.c_str(), charge, ion_type, frag_index);
      }
      else
      {
        std::string mass_tok = ScanCommandQueue::formatMassToken(precursor_mass / 1000.0, charge, '\0', 0);
        std::snprintf(cmd.scan_description, 16, "%sE%s@%d",
                     id_str.c_str(), mass_tok.c_str(), charge);
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
      double rt, ScanCommandQueue& queue, ProteoformTracker* tracker, int precursor_id)
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

    return feedResultImpl_(tracking_id, ms2_deconv, mzs, ints, length, queue, tracker, precursor_id);
  }

  Exploration::FeedResultInfo Exploration::feedResultImpl_(int tracking_id,
      const DeconvolvedSpectrum& msn_deconv,
      const double* mzs, const double* ints, int length,
      ScanCommandQueue& queue, ProteoformTracker* tracker, int precursor_id)
  {
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

    v.result = msn_deconv;
    double remaining_ratio = -1.0;
    FragmentAnalysis::ProteoformMatch frag{};
    v.score = computeExplorationScore_(group.exploration_metric, msn_deconv, group, mzs, ints, length, &remaining_ratio, &frag, v.activation_type);
    v.tic_coverage = computeTICCoverage_(msn_deconv);
    v.fragment_count = frag.total_match_count;
    v.received = true;

    // Handle baseline (present for EVERY exploration group as of #18 baseline-on-all)
    if (v.is_baseline)
    {
      group.baseline_intensity = precursorWindowIntensity_(group, mzs, ints, length);
      group.has_baseline = true;
      v.score = 0.0;

      // Empty baseline window => no CE variant can be scored. Abort the group:
      // cancel its still-queued / in-flight child scans (no follow-up production scan),
      // and account for the cancelled children so the all-received path below cleans up.
      if (group.exploration_metric == ExplorationMetric::RemainingPrecursor &&
          group.baseline_intensity <= 0.0)
      {
        group.baseline_failed = true;

        std::vector<int> child_ids;
        for (const auto& cv : group.variants)
          if (!cv.is_baseline) child_ids.push_back(queue.decode(cv.tracking_id));

        std::vector<int> removed = queue.cancelByScanIds(child_ids);

        // Cancelled children will never return — mark them received (score 0) so the group
        // can complete, and erase their routing so any already-dispatched (in-flight) result
        // that returns later is a harmless no-op (it won't re-create the erased group).
        for (int cid : removed)
        {
          auto cit = variant_tracking_map_.find(cid);
          if (cit == variant_tracking_map_.end()) continue;
          int gidx = cit->second.variant_index; 
          if (gidx >= 0 && gidx < static_cast<int>(group.variants.size()))
          {
            group.variants[gidx].received = true;
            group.variants[gidx].score = 0.0;
          }
          variant_tracking_map_.erase(cit);
        }
      }
    }

    // remaining_ratio is metric-independent: every exploration group now carries a CE-0 baseline (#18), so
    // ratio the variant's own isolation-window intensity against the baseline (same window-sum semantics as
    // computeRemainingPrecursorScore_). The un-fragmented baseline variant ratios against itself -> 1.0; a
    // fragment scan's surviving precursor is depleted -> < 1. Empty/late baseline -> stays -1.0 (N/A). The
    // baseline dequeues FIFO-first, so has_baseline is already set when its CE siblings are scored.
    if (group.has_baseline && group.baseline_intensity > 0.0)
    {
      remaining_ratio = v.is_baseline
                          ? 1.0
                          : precursorWindowIntensity_(group, mzs, ints, length) / group.baseline_intensity;
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
    // @Claude Please implement this more generically. At some point we might want to add MS4 or whatever
    if (group.msn_level >= 3 && v.cmd.num_stages >= 2)
    {
      info.stage0_collision_energy = v.cmd.stages[0].collision_energy;
      info.stage0_activation_type  = std::string(v.cmd.stages[0].activation_type);
      info.stage0_reaction_time    = v.cmd.stages[0].reaction_time;
    }
    info.score = v.score;
    info.tic_coverage = v.tic_coverage;
    info.fragment_count = v.fragment_count;
    info.exploration_metric = static_cast<int>(group.exploration_metric);
    info.matched_protein = frag.matched_protein;
    info.proteoform_sequence = frag.proteoform_sequence;
    info.identification_result = frag;

    // Stage this MS2 variant into the ProteoformTracker for pooling (MS3 evidence folds via foldMs3,
    // not here). The CE-0 baseline is excluded: it carries no fragmentation evidence and must never
    // pollute the pooled model.
    if (tracker != nullptr && group.msn_level == 2 && !v.is_baseline)
    {
      tracker->feedScan(precursor_id,
                        2,
                        Ms2Params{v.collision_energy, v.activation_type, v.reaction_time},
                        tracking_id,
                        v.result,
                        frag,
                        frag.score,
                        v.cmd);
    }
    
    info.remaining_ratio = remaining_ratio;
    std::string parent_enc = ScanCommandQueue::encode(group.originating_cmd.scan_id);
    std::strncpy(info.parent_scan_id, parent_enc.c_str(), 3);
    info.parent_scan_id[3] = '\0';
    
    // @Claude solve inline and make generic for ms3
    auto buildMS2ContextForVariant = [&](int group_variant_index)
    {
      MS2Context ctx;
      ctx.proteoform_sequence = (group.proteoform_ctx.region_start >= 0)
          ? config_.characterization().protein_sequence.substr(group.proteoform_ctx.region_start,
              group.proteoform_ctx.region_end - group.proteoform_ctx.region_start)
          : "";
      ctx.start_pos = group.proteoform_ctx.region_start;
      ctx.end_pos = group.proteoform_ctx.region_end;
      ctx.ptm_sites = group.proteoform_ctx.ptm_sites;
      ctx.fragment_ion_type = group.fragment_ion_type;
      ctx.fragment_ion_index = group.fragment_ion_index;
      if (group.originating_cmd.num_stages > 0)
      {
        ctx.ms1_precursor_mass = group.originating_cmd.mono_mass;
        ctx.ms1_precursor_mz = group.originating_cmd.stages[0].precursor_mz;
        ctx.ms1_precursor_charge = group.originating_cmd.stages[0].charge_state;
      }

      const auto& variant_cmd = group.variants[group_variant_index].cmd;
      if (variant_cmd.num_stages >= 2)
      {
        ctx.fragment_mz = variant_cmd.stages[1].precursor_mz;
        ctx.fragment_charge = variant_cmd.stages[1].charge_state;
        ctx.ms2_isolation_width = group.originating_cmd.stages[0].isolation_width;
        ctx.ms2_charge_intensity = group.originating_cmd.precursor_intensity;
        ctx.ms2_window_snr = group.originating_cmd.window_snr;
        ctx.ms3_isolation_width = variant_cmd.stages[1].isolation_width;
        ctx.ms3_charge_intensity = variant_cmd.precursor_intensity_s1;
        ctx.ms3_window_snr = variant_cmd.window_snr;
      }
      else if (variant_cmd.num_stages >= 1)
      {
        ctx.ms2_isolation_width = variant_cmd.stages[0].isolation_width;
        ctx.ms2_charge_intensity = variant_cmd.precursor_intensity;
        ctx.ms2_window_snr = variant_cmd.window_snr;
        ctx.ms1_precursor_mass = group.precursor_mass;
        ctx.ms1_precursor_mz = group.precursor_mz;
        ctx.ms1_precursor_charge = group.precursor_charge;
      }
      ctx.fragment_mass = group.precursor_mass;
      return ctx;
    };

    info.ms2_context = buildMS2ContextForVariant(variant_index);

    if (group.msn_level >= 3 && !info.ms2_context.proteoform_sequence.empty())
    {
      info.proteoform_sequence = info.ms2_context.proteoform_sequence;
      info.matched_protein = config_.targeting().fasta_file.empty()
          ? config_.characterization().protein_sequence : config_.targeting().fasta_file;
    }
    
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

      std::vector<FragmentAnalysis::ProteoformMatch> detailed_results;
      auto calibrated_scores = MS3FragmentMatcher::calibrateAndScore(
        variant_spectra,
        config_.characterization().protein_sequence,
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
      info.score = group.variants[variant_index].score;

      // Return calibrated identification rows for other variants in the completed group.
      for (size_t vi = 0; vi < group.variants.size(); ++vi)
      {
        if (static_cast<int>(vi) == variant_index) continue;
        const auto& gv = group.variants[vi];
        if (gv.identification_result.fragments.empty()) continue;

        FeedResultInfo::IdentificationRowInfo row;
        row.tracking_id = gv.tracking_id;
        row.identification_result = gv.identification_result;
        row.ms2_context = buildMS2ContextForVariant(static_cast<int>(vi));
        row.tic_coverage = gv.tic_coverage;   // this variant's OWN tic (per-scan)
        info.additional_identification_rows.push_back(row);
      }
    }

    // Select the winner — unless the group was aborted (empty baseline) or no variant scored.
    int best_idx = -1;
    double best_score = -1.0;
    if (!group.baseline_failed)
    {
      for (int i = 0; i < static_cast<int>(group.variants.size()); ++i)
      {
        if (group.variants[i].is_baseline) continue;  // skip baseline
        if (group.variants[i].score > best_score)
        {
          best_score = group.variants[i].score;
          best_idx = i;
        }
      }
    }

    // Finalize the proteoform model before the next-level dispatch below
    if (tracker != nullptr && group.msn_level == 2)
    {
      tracker->finalizeMS2(precursor_id);
    }

    if (group.baseline_failed || best_idx < 0)
    {
      // Empty-baseline abort (or, defensively, no scorable variant): no winner, no follow-up scan.
      // Children were already cancelled in the baseline hook. Falls through to the cleanup below
      // (this replaces the old leaking `if (best_idx < 0) return info;`).
      group.complete = true;
      std::cout << "[EXPL-ABORT] group=" << group.group_id
                << " reason=" << (group.baseline_failed ? "empty-baseline" : "no-winner")
                << " winner=none" << std::endl;
    }
    else
    {
      group.winner_index = best_idx;
      group.complete = true;
      group.variants[best_idx].result.getOrCreateOptimizationMetadata().is_best_variant = true;
      info.winner_tracking_id = group.variants[best_idx].tracking_id;

      std::cout << "[EXPL-WINNER] group=" << group.group_id
                << " winner_idx=" << best_idx
                << " activation=" << group.variants[best_idx].activation_type
                << " CE=" << group.variants[best_idx].collision_energy
                << " RT=" << group.variants[best_idx].reaction_time
                << " score=" << best_score << std::endl;

      const auto& level_config = config_.level(group.msn_level);

      // Send production scan after exploration
      if (!level_config.overrides.empty())
      {
        ScanConfig prod_config = level_config.scans[0];
        prod_config.collision_energy = static_cast<int>(group.variants[best_idx].collision_energy);
        prod_config.reaction_time = group.variants[best_idx].reaction_time;
        prod_config.activation = group.variants[best_idx].activation_type;

        ScanCommand prod_cmd;
        if (group.msn_level >= 3)
        {
          // Reconstruct the winning variant's stage-1 fragment scores for the production MS3 (buildMS3
          // keeps its explicit frag_scores param; the other callers pass genuine scores, so do not fold
          // this into buildMS3 — an ms2-command would carry zero-initialised *_s1 fields).
          const ScanCommand& wcmd = group.variants[best_idx].cmd;
          FragmentAnalysis::FragmentScores wfs;
          wfs.mono_mass = wcmd.mono_mass_s1;                 
          wfs.qscore = wcmd.qscore_s1;
          wfs.charge_cos = wcmd.charge_cos_s1;               
          wfs.charge_snr = wcmd.charge_snr_s1;
          wfs.iso_cos = wcmd.iso_cos_s1;                     
          wfs.snr = wcmd.snr_s1;
          wfs.charge_score = wcmd.charge_score_s1;           
          wfs.ppm_error = wcmd.ppm_error_s1;
          wfs.precursor_intensity = wcmd.precursor_intensity_s1; 
          wfs.peakgroup_intensity = wcmd.peakgroup_intensity_s1;
          prod_cmd = queue.buildMS3(group.variants[best_idx].cmd, prod_config,
                                     group.precursor_mz, group.precursor_charge,
                                     group.isolation_width, group.variants[best_idx].cmd.scan_id,
                                     group.fragment_ion_type, group.fragment_ion_index, 1, wfs,
                                     nullptr, group.proteoform_ctx);
        }
        else
        {
          prod_cmd = queue.buildMS2(group.precursor_pg, group.precursor_charge, prod_config, 2, 0);
        }

        prod_cmd.faims_cv = group.faims_cv;
        if (best_idx >= 0 && best_idx < static_cast<int>(group.variants.size()))
        {
          const double inherited = group.variants[best_idx].cmd.window_snr;
          if (inherited >= 0.0) prod_cmd.window_snr = inherited;
        }

        std::string prod_id = ScanCommandQueue::encode(prod_cmd.scan_id);
        std::cout << "[TRACK-CREATE] id=" << prod_id
                  << " ms_level=" << group.msn_level << " type=production"
                  << std::endl;

        info.commands.push_back(prod_cmd);
        if (group.msn_level >= 3)
          info.ms3_context_cache.emplace_back(prod_cmd.scan_id, buildMS2ContextForVariant(best_idx));
      }
      // inititate next level - @Claude make interface more generic. Should be for last level with selection criterion
      else if (group.msn_level < 3)
      {
        auto next_nlr = initiateNextLevel(group.msn_level,
            group.variants[best_idx].result, group.faims_cv, queue,
            &group.variants[best_idx].cmd, tracker, precursor_id);  // P5: same precursor_id model key
        info.commands.insert(info.commands.end(), next_nlr.commands.begin(), next_nlr.commands.end());
      }
      else  // overrides empty && msn_level >= 3: MS3 exploration with no production scan -> fold the winning variant now @Claude should be more generic and happen after all winning variants or if a production scan is found.
      {
        if (tracker != nullptr && group.fragment_ion_type != '\0')
        {
          const ExplorationVariant& w = group.variants[best_idx];
          if (!w.identification_result.fragments.empty())
          {
            Ms2Params stage0{w.cmd.stages[0].collision_energy,
                             std::string(w.cmd.stages[0].activation_type),
                             w.cmd.stages[0].reaction_time};
            tracker->feedScan(precursor_id, /*ms_level=*/3, stage0,
                              /*scan_id=*/queue.decode(w.tracking_id), w.result,
                              w.identification_result, w.identification_result.score, w.cmd);
            const std::string trig = std::string(1, group.fragment_ion_type)
                                   + std::to_string(group.fragment_ion_index);
            tracker->foldMs3(precursor_id, trig, w.tracking_id);
          }
        }
      }
    }

    // F1 fallback: stamp any command whose builder left the parent empty with the group's originating id
    // (only fill empties — never overwrite a builder-set parent, e.g. buildMS3's immediate MS2 parent).
    // Relocated here from the former processScan feedResult push-loops so callers receive fully-parented
    // commands and no longer post-stamp.
    for (auto& c : info.commands)
    {
      if (c.parent_scan_id[0] == '\0')
        std::strncpy(c.parent_scan_id, info.parent_scan_id, 4);
    }

    // Pre-encode the pushed children's tracking ids here (relocated from the former processScan
    // expl_children loops). encode() is the same pure base-94 mapping, applied in info.commands order,
    // so scan_results.tsv's child_ids column is byte-identical.
    info.child_ids.clear();
    for (const auto& c : info.commands)
      info.child_ids.push_back(ScanCommandQueue::encode(c.scan_id));

    for (const auto& vr : group.variants)
      variant_tracking_map_.erase(queue.decode(vr.tracking_id));

    // (Model finalize was hoisted above, before the next-level MS3 dispatch — 9b.)

    active_groups_.erase(git);
    return info;
  }

  Exploration::NextLevelResult Exploration::initiateNextLevel(int msn_level,
      const DeconvolvedSpectrum& result, double faims_cv, ScanCommandQueue& queue,
      const ScanCommand* ms_ctx, ProteoformTracker* tracker, int precursor_id)
  {
    NextLevelResult nlr;

    int next_level = msn_level + 1;
    const auto& this_cfg = config_.level(msn_level);
    const auto& next_cfg = config_.level(next_level);
    if (this_cfg.selection == SelectionMetric::None) return nlr;
    if (this_cfg.scans.empty() || next_cfg.scans.empty()
        || next_cfg.selection == SelectionMetric::None) return nlr;

    const auto& seq = config_.characterization().protein_sequence;
    int num_targets = this_cfg.max_targets;

    // Extract activation type from the scan command that produced this result
    std::string scan_activation = (ms_ctx != nullptr && ms_ctx->num_stages > 0)
        ? std::string(ms_ctx->stages[0].activation_type)
        : config_.level(msn_level).scans[0].activation;

    // Get fragment ion targets via FragmentAnalysis
    // @Claude should completely move to proteoform tracker
    DeconvolvedSpectrum result_copy = result;
    std::vector<double> masses(num_targets), qscores(num_targets);
    std::vector<double> wstarts(num_targets), wends(num_targets);
    std::vector<int> charges(num_targets);
    std::vector<char> ion_types(num_targets, '\0');
    std::vector<int> frag_indices(num_targets, 0);
    std::vector<FragmentAnalysis::FragmentScores> frag_scores(num_targets);  // stage-1 (fragment) scoring for MS3 cmd rows

    int found = 0;
    FragmentAnalysis::ProteoformMatch frag_result;

    // @Claude this logic should move into ProteoformTracker
    switch (next_cfg.selection)
    {
      case SelectionMetric::Intensity:
      case SelectionMetric::QScore:
        found = fragments_.getTopFragmentMatches(seq, num_targets,
            masses.data(), qscores.data(), charges.data(),
            wstarts.data(), wends.data(),
            ion_types.data(), frag_indices.data(), result_copy, frag_result, scan_activation, 0.0, frag_scores.data());
        break;
      case SelectionMetric::TerminalFragments:
        found = fragments_.getTerminalFragmentIons(seq, num_targets,
            masses.data(), qscores.data(), charges.data(),
            wstarts.data(), wends.data(),
            ion_types.data(), frag_indices.data(), result_copy, frag_result, scan_activation, 0.0, frag_scores.data());
        break;
      case SelectionMetric::AmbiguityResolution:
        found = fragments_.getAmbiguityEnclosingIons(seq, num_targets,
            masses.data(), qscores.data(), charges.data(),
            wstarts.data(), wends.data(),
            ion_types.data(), frag_indices.data(), result_copy, frag_result, scan_activation, 0.0, frag_scores.data());
        break;
      default:
        break;
    }

    // Cache proteoform context for MS3 subsequence scoring
    // @Claude this logic should move into ProteoformTracker
    MS3FragmentMatcher::ProteoformContext proto_ctx;
    if (next_level >= 3)
    {
      proto_ctx.region_start = frag_result.region_start;
      proto_ctx.region_end = frag_result.region_end;
      proto_ctx.ptm_sites = frag_result.ptm_sites;
      // If no truncation detected, use full protein sequence bounds
      if (proto_ctx.region_start < 0)
        proto_ctx.region_start = 0;
      if (proto_ctx.region_end < 0)
        proto_ctx.region_end = static_cast<int>(seq.size());
    }

    // Populate fragment matching metadata
    // @Claude this logic should move into ProteoformTracker
    nlr.fragment_count = frag_result.total_match_count;
    nlr.proteoform_match = frag_result;
    if (!seq.empty() && found > 0)
    {
      nlr.matched_protein = config_.targeting().fasta_file.empty()
          ? config_.characterization().protein_sequence : config_.targeting().fasta_file;
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

    if (tracker != nullptr && next_level >= 3 && ms_ctx != nullptr && ms_ctx->num_stages > 0 && found > 0)
    {
      // @Claude bail early if the ProteoformModel does not exist.
      const ProteoformModel* existing = tracker->getModel(precursor_id);
      if (existing == nullptr || existing->proteoform_sequence.empty())
      {
        Ms2Params p;
        p.collision_energy = ms_ctx->stages[0].collision_energy;
        p.activation_type  = std::string(ms_ctx->stages[0].activation_type);
        p.reaction_time    = ms_ctx->stages[0].reaction_time;
        tracker->feedScan(precursor_id, /*ms_level=*/2, p, ms_ctx->scan_id, result, frag_result,
                          frag_result.score, *ms_ctx);
        tracker->finalizeMS2(precursor_id);
      }
    }

    const ProteoformModel* model = (tracker != nullptr && next_level >= 3 && ms_ctx != nullptr)
        ? tracker->getModel(precursor_id) : nullptr;
    if (model != nullptr && !model->proteoform_sequence.empty())
    {
      auto targets = tracker->planNextScans(precursor_id);

      // @Claude this should be bound by num_tartgets
      for (const Ms3Target& target : targets)
      {
        const int frag_charge = std::abs(target.frag_charge);
        if (charge_floor > 0 && frag_charge < charge_floor) continue;
        const char ion_type = target.ion_type.empty() ? '\0' : target.ion_type[0];

        if (config_.hasExploration(next_level))
        {
          const double half = target.iso_width / 2.0;
          const double wstart = target.frag_mz - half;
          const double wend = target.frag_mz + half;
          PeakGroup frag_pg(frag_charge, frag_charge, true);
          frag_pg.setMonoisotopicMass(target.frag_mass);
          FLASHHelperClasses::LogMzPeak lp_lo;
          lp_lo.mz = wstart;
          lp_lo.abs_charge = frag_charge;
          frag_pg.push_back(lp_lo);
          FLASHHelperClasses::LogMzPeak lp_hi;
          lp_hi.mz = wend;
          lp_hi.abs_charge = frag_charge;
          frag_pg.push_back(lp_hi);

          auto sub_cmds = initiate(next_level, frag_pg, frag_charge, queue, &result.getOriginalSpectrum(), ms_ctx,
                                   ion_type, target.ion_index, proto_ctx, target.stage1_scores, &target.stage0_params);
          nlr.commands.insert(nlr.commands.end(), sub_cmds.begin(), sub_cmds.end());

          for (size_t sci = 0; sci < sub_cmds.size(); ++sci)
          {
            MS2Context mc;
            mc.proteoform_sequence = nlr.proteoform_sequence;
            mc.start_pos = proto_ctx.region_start;
            mc.end_pos = proto_ctx.region_end;
            mc.ptm_sites = proto_ctx.ptm_sites;
            mc.ms1_precursor_mass = ms_ctx->mono_mass;
            mc.ms1_precursor_mz = ms_ctx->stages[0].precursor_mz;
            mc.ms1_precursor_charge = ms_ctx->stages[0].charge_state;
            mc.fragment_ion_type = ion_type;
            mc.fragment_ion_index = target.ion_index;
            mc.fragment_mass = target.frag_mass;
            mc.fragment_mz = target.frag_mz;
            mc.fragment_charge = frag_charge;
            mc.ms2_isolation_width = ms_ctx->stages[0].isolation_width;
            mc.ms2_charge_intensity = ms_ctx->precursor_intensity;
            mc.ms2_window_snr = ms_ctx->window_snr;
            if (sub_cmds[sci].num_stages >= 2)
            {
              const ScanCommand& sub = sub_cmds[sci];
              mc.ms3_isolation_width = sub.stages[1].isolation_width;
              mc.ms3_charge_intensity = sub.precursor_intensity_s1;
              mc.ms3_window_snr = sub.window_snr;
            }
            nlr.ms3_contexts.push_back(mc);
          }
        }
        else
        {
          // Single MS3: build directly from the target descriptors, applying the model's stage[0].
          ScanCommand cmd = queue.buildMS3(*ms_ctx, next_scan_config, target.frag_mz, frag_charge,
                                           target.iso_width, ms_ctx->scan_id,
                                           ion_type, target.ion_index, 1, target.stage1_scores, &target.stage0_params,
                                           proto_ctx);
          cmd.faims_cv = faims_cv;

          // Item 2: stamp the MS3 fragment-window SNR onto the command before push (same inputs/formula
          // as the legacy direct path — computed over the MS2 source spectrum).
          if (cmd.num_stages >= 2)
          {
            const double half3 = cmd.stages[1].isolation_width / 2.0;
            cmd.window_snr = FragmentAnalysis::windowSnr(result.getOriginalSpectrum(),
                cmd.stages[1].precursor_mz - half3, cmd.stages[1].precursor_mz + half3,
                cmd.precursor_intensity_s1);
          }

          std::string id_str = ScanCommandQueue::encode(cmd.scan_id);
          std::cout << "[TRACK-CREATE] id=" << id_str
                    << " ms_level=" << next_level << " type=next_level"
                    << std::endl;

          nlr.commands.push_back(cmd);

          MS2Context mc;
          mc.proteoform_sequence = nlr.proteoform_sequence;
          mc.start_pos = proto_ctx.region_start;
          mc.end_pos = proto_ctx.region_end;
          mc.ptm_sites = proto_ctx.ptm_sites;
          mc.ms1_precursor_mass = ms_ctx->mono_mass;
          mc.ms1_precursor_mz = ms_ctx->stages[0].precursor_mz;
          mc.ms1_precursor_charge = ms_ctx->stages[0].charge_state;
          mc.fragment_ion_type = ion_type;
          mc.fragment_ion_index = target.ion_index;
          mc.fragment_mass = target.frag_mass;
          mc.fragment_mz = target.frag_mz;
          mc.fragment_charge = frag_charge;
          mc.ms2_isolation_width = ms_ctx->stages[0].isolation_width;
          mc.ms2_charge_intensity = ms_ctx->precursor_intensity;
          mc.ms2_window_snr = ms_ctx->window_snr;
          if (cmd.num_stages >= 2)
          {
            mc.ms3_isolation_width = cmd.stages[1].isolation_width;
            mc.ms3_charge_intensity = cmd.precursor_intensity_s1;
            mc.ms3_window_snr = cmd.window_snr;
          }
          nlr.ms3_contexts.push_back(mc);
        }
      }

      return nlr;
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
      double* out_remaining_ratio, FragmentAnalysis::ProteoformMatch* out_frag,
      const std::string& activation_type) const
  {
    FragmentAnalysis::ProteoformMatch fmr;
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
        return static_cast<double>(fmr.total_match_count);
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

  double Exploration::precursorWindowIntensity_(const ExplorationGroup& group,
      const double* mzs, const double* ints, int length) const
  {
    if (length <= 0 || mzs == nullptr || ints == nullptr)
      return 0.0;

    double iso_half = group.isolation_width / 2.0;
    double mz_low = group.precursor_mz - iso_half;
    double mz_high = group.precursor_mz + iso_half;

    double sum = 0.0;
    for (int i = 0; i < length; ++i)
    {
      if (mzs[i] >= mz_low && mzs[i] <= mz_high)
        sum += ints[i];
    }
    return sum;
  }

  double Exploration::computeRemainingPrecursorScore_(const ExplorationGroup& group,
      const double* mzs, const double* ints, int length, double* out_ratio) const
  {
    if (out_ratio) *out_ratio = -1.0;

    if (length <= 0 || mzs == nullptr || ints == nullptr)
      return 0.0;

    double remaining_intensity = precursorWindowIntensity_(group, mzs, ints, length);

    // Reference: baseline isolation-window intensity (CE=0 scan)
    double reference;
    if (group.has_baseline)
    {
      reference = group.baseline_intensity;
      if (reference <= 0.0)
        return 0.0;  // Baseline failed (no in-window signal): score 0; out_ratio stays -1.0 (N/A)
    }
    else
    {
      return 0.0;  // Baseline not yet received: score 0; out_ratio stays -1.0 (N/A)
    }

    double ratio = remaining_intensity / reference;
    if (out_ratio) *out_ratio = ratio;
    double target = config_.level(group.msn_level).remaining_precursor_target;
    double deviation = std::abs(ratio - target);
    double score = 1.0 - deviation;
    if (score < 0.0) score = 0.0;  // floor: score is always in [0, 1]
    return score;
  }

  FragmentAnalysis::ProteoformMatch Exploration::computeFragmentMatch_(const DeconvolvedSpectrum& spec, int msn_level,
                                                                      const std::string& activation_type) const
  {
    FragmentAnalysis::ProteoformMatch result;
    const auto& seq = config_.characterization().protein_sequence;
    if (seq.empty() || spec.empty())
      return result;

    DeconvolvedSpectrum spec_copy = spec;

    const int max_matches = 100;
    std::vector<double> masses(max_matches), qscores(max_matches);
    std::vector<double> wstarts(max_matches), wends(max_matches);
    std::vector<int> charges(max_matches);
    std::vector<char> ion_types(max_matches, '\0');
    std::vector<int> frag_indices(max_matches, 0);

    fragments_.getTopFragmentMatches(
        seq, max_matches,
        masses.data(), qscores.data(), charges.data(),
        wstarts.data(), wends.data(),
        ion_types.data(), frag_indices.data(),
        spec_copy, result, activation_type,
        config_.level(msn_level).exploration_tolerance_ppm);

    if (result.total_match_count > 0)
    {
      result.matched_protein = config_.targeting().fasta_file.empty()
          ? config_.characterization().protein_sequence : config_.targeting().fasta_file;
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
