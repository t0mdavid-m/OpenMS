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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <iostream>

namespace OpenMS
{

  Exploration::Exploration(const Config& config, Deconvolution& deconv, FragmentAnalysis& fragments)
    : config_(config), deconv_(deconv), fragments_(fragments)
  {
  }

  std::vector<double> Exploration::buildCEVariants_(double ce_min, double ce_max, double ce_step) const
  {
    std::vector<double> ces;
    for (double ce = ce_min; ce <= ce_max + 1e-9; ce += ce_step)
      ces.push_back(ce);
    return ces;
  }

  std::vector<ScanCommand> Exploration::initiate(int msn_level, const PeakGroup& pg, int charge,
      double faims_cv, ScanCommandQueue& queue, const ScanCommand* ms_ctx)
  {
    std::vector<ScanCommand> commands;

    const auto& cfg = config_.level(msn_level);
    if (!config_.hasExploration(msn_level)) return commands;

    std::vector<double> ces = buildCEVariants_(cfg.ce_min, cfg.ce_max, cfg.ce_step);
    if (ces.empty()) return commands;

    // Prepend CE=0 baseline for RemainingPrecursor metric
    bool needs_baseline = (cfg.exploration == ExplorationMetric::RemainingPrecursor);
    if (needs_baseline)
      ces.insert(ces.begin(), 0.0);

    // Compute precursor_mz and isolation_width from PeakGroup
    auto [mz1, mz2] = pg.getMzRange(charge);
    double precursor_mz = (mz1 + mz2) / 2.0;
    double isolation_width = mz2 - mz1;
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

    // Build base ScanConfig from the level's primary scan config, then apply overrides
    ScanConfig base_config = cfg.scans[0];
    base_config.applyOverrides(cfg.overrides);

    for (int i = 0; i < static_cast<int>(ces.size()); ++i)
    {
      ExplorationVariant v;
      v.variant_index = (needs_baseline && i == 0) ? -1 : (needs_baseline ? i - 1 : i);
      v.collision_energy = ces[i];
      v.is_baseline = (needs_baseline && i == 0);
      v.activation_type = base_config.activation;

      ScanConfig variant_config = base_config;
      variant_config.collision_energy = static_cast<int>(ces[i]);

      int expl_priority = (msn_level >= 3) ? 1 : 2;  // MS3 variants = p1, MS2 variants = p2
      ScanCommand cmd = queue.buildMS2(pg, charge, variant_config, expl_priority);
      cmd.faims_cv = faims_cv;

      int id_int = cmd.scan_id;
      std::string id_str = ScanCommandQueue::encode(id_int);
      v.tracking_id = id_str;
      std::snprintf(cmd.scan_description, 16, "%sE%.1f@%d",
                   id_str.c_str(), precursor_mass / 1000.0, charge);

      group.variants.push_back(v);
      variant_tracking_map_[id_int] = {group.group_id, i};
      commands.push_back(cmd);

      std::cout << "[TRACK-CREATE] id=" << id_str
                << " ms_level=" << msn_level << " type=exploration"
                << " CE=" << ces[i] << std::endl;
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
      deconv_.deconvolveMSn(mzs, ints, length, rt,
                            group.precursor_mass, group.precursor_charge);
      ms2_deconv = deconv_.storedMS2();
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
    v.score = computeExplorationScore_(group.exploration_metric, ms2_deconv, group, mzs, ints, length);
    v.tic_coverage = computeTICCoverage_(ms2_deconv);
    v.fragment_count = static_cast<int>(ms2_deconv.size());
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
    info.score = v.score;
    info.tic_coverage = v.tic_coverage;
    info.fragment_count = v.fragment_count;
    info.exploration_metric = static_cast<int>(group.exploration_metric);

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
              << " CE=" << group.variants[best_idx].collision_energy
              << " score=" << best_score << std::endl;

    const auto& level_config = config_.level(group.msn_level);
    if (!level_config.overrides.empty())
    {
      ScanConfig prod_config = level_config.scans[0];
      prod_config.collision_energy = static_cast<int>(group.variants[best_idx].collision_energy);
      prod_config.activation = group.variants[best_idx].activation_type;

      ScanCommand prod_cmd = queue.buildMS2(group.precursor_pg, group.precursor_charge, prod_config);
      prod_cmd.faims_cv = group.faims_cv;
      prod_cmd.priority = 2;

      std::string prod_id = ScanCommandQueue::encode(prod_cmd.scan_id);
      std::cout << "[TRACK-CREATE] id=" << prod_id
                << " ms_level=" << group.msn_level << " type=production"
                << std::endl;

      info.commands.push_back(prod_cmd);
    }
    else
    {
      auto next_nlr = initiateNextLevel(group.msn_level,
          group.variants[best_idx].result, group.faims_cv, queue,
          &group.originating_cmd);
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

    if (config_.hasExploration(next_level))
    {
      // Recursive exploration at next level
      for (int ti = 0; ti < num_targets; ++ti)
      {
        PeakGroup frag_pg(std::abs(charges[ti]), std::abs(charges[ti]), true);
        frag_pg.setMonoisotopicMass(masses[ti]);
        FLASHHelperClasses::LogMzPeak lp;
        lp.mz = (wstarts[ti] + wends[ti]) / 2.0;
        lp.abs_charge = std::abs(charges[ti]);
        frag_pg.push_back(lp);

        auto sub_cmds = initiate(next_level, frag_pg, std::abs(charges[ti]), faims_cv, queue, ms_ctx);
        nlr.commands.insert(nlr.commands.end(), sub_cmds.begin(), sub_cmds.end());
      }
    }
    else
    {
      // Direct command building for each fragment target
      for (int ti = 0; ti < num_targets; ++ti)
      {
        double frag_mz = (wstarts[ti] + wends[ti]) / 2.0;
        int frag_charge = std::abs(charges[ti]);
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
      const double* mzs, const double* ints, int length) const
  {
    switch (metric)
    {
      case ExplorationMetric::MassCount:
        return computeMassCount_(spec);
      case ExplorationMetric::RemainingPrecursor:
        return computeRemainingPrecursorScore_(group, mzs, ints, length);
      case ExplorationMetric::FragmentCount:
        return computeFragmentCount_(spec);
      default:
        return computeMassCount_(spec);
    }
  }

  double Exploration::computeMassCount_(const DeconvolvedSpectrum& spec) const
  {
    return static_cast<double>(spec.size());
  }

  double Exploration::computeRemainingPrecursorScore_(const ExplorationGroup& group,
      const double* mzs, const double* ints, int length) const
  {
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
        return 0.0;  // Baseline failed: ratio = 1/1, score = 0
    }
    else
    {
      return 0.0;  // Baseline not yet received: ratio = 1/1, score = 0
    }

    double ratio = remaining_intensity / reference;
    double score = 1.0 - ratio;
    if (score < 0.0) score = 0.0;
    if (score > 1.0) score = 1.0;
    return score;
  }

  double Exploration::computeFragmentCount_(const DeconvolvedSpectrum& spec) const
  {
    const auto& seq = config_.targeting().protein_sequence;
    if (seq.empty() || spec.empty())
      return 0.0;

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
        spec_copy);

    return static_cast<double>(count);
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
