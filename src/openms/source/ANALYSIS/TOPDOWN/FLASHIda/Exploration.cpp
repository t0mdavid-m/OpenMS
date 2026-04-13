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
      double faims_cv, ScanCommandQueue& queue)
  {
    std::vector<ScanCommand> commands;

    const auto& cfg = config_.level(msn_level);
    if (!config_.hasExploration(msn_level)) return commands;

    std::vector<double> ces = buildCEVariants_(cfg.ce_min, cfg.ce_max, cfg.ce_step);
    if (ces.empty()) return commands;

    // Compute precursor_mz and isolation_width from PeakGroup
    auto [mz1, mz2] = pg.getMzRange(charge);
    double precursor_mz = (mz1 + mz2) / 2.0;
    double isolation_width = mz2 - mz1;
    double precursor_mass = pg.getMonoMass();
    if (precursor_mass <= 0) precursor_mass = precursor_mz;  // defensive fallback

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

    // Build base ScanConfig from the level's primary scan config, then apply overrides
    ScanConfig base_config = cfg.scans[0];
    base_config.applyOverrides(cfg.overrides);

    for (int i = 0; i < static_cast<int>(ces.size()); ++i)
    {
      ExplorationVariant v;
      v.variant_index = i;
      v.collision_energy = ces[i];
      v.activation_type = cfg.exploration_activation;

      ScanConfig variant_config = base_config;
      variant_config.collision_energy = static_cast<int>(ces[i]);
      variant_config.activation = cfg.exploration_activation;

      ScanCommand cmd = queue.buildMS2(pg, charge, variant_config, 0);
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

  std::vector<ScanCommand> Exploration::feedResult(int tracking_id,
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

  std::vector<ScanCommand> Exploration::feedResultForTest(int tracking_id,
      const DeconvolvedSpectrum& ms2_deconv, double rt, ScanCommandQueue& queue)
  {
    return feedResultImpl_(tracking_id, ms2_deconv, nullptr, nullptr, 0, rt, queue);
  }

  std::vector<ScanCommand> Exploration::feedResultImpl_(int tracking_id,
      const DeconvolvedSpectrum& ms2_deconv,
      const double* mzs, const double* ints, int length,
      double rt, ScanCommandQueue& queue)
  {
    (void)rt;
    std::vector<ScanCommand> commands;

    // Look up the variant reference
    auto vit = variant_tracking_map_.find(tracking_id);
    if (vit == variant_tracking_map_.end()) return commands;

    int group_id = vit->second.group_id;
    int variant_index = vit->second.variant_index;
    variant_tracking_map_.erase(vit);

    // Find the group
    auto git = active_groups_.find(group_id);
    if (git == active_groups_.end()) return commands;
    ExplorationGroup& group = git->second;

    if (variant_index < 0 || variant_index >= static_cast<int>(group.variants.size())) return commands;
    ExplorationVariant& v = group.variants[variant_index];
    if (v.received) return commands;

    v.result = ms2_deconv;
    v.score = computeExplorationScore_(group.exploration_metric, ms2_deconv, group, mzs, ints, length);
    v.tic_coverage = computeTICCoverage_(ms2_deconv);
    v.fragment_count = static_cast<int>(ms2_deconv.size());
    v.received = true;

    auto& meta = v.result.getOrCreateOptimizationMetadata();
    meta.group_id = group.group_id;
    meta.variant_index = variant_index;
    meta.total_variants = static_cast<int>(group.variants.size());
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
    meta.exploration_scans = static_cast<int>(group.variants.size());

    bool all_received = std::all_of(group.variants.begin(), group.variants.end(),
                                    [](const ExplorationVariant& x){ return x.received; });
    if (!all_received) return commands;

    int best_idx = 0;
    double best_score = group.variants[0].score;
    for (int i = 1; i < static_cast<int>(group.variants.size()); ++i)
    {
      if (group.variants[i].score > best_score)
      {
        best_score = group.variants[i].score;
        best_idx = i;
      }
    }
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
      prod_cmd.priority = 1;

      std::string prod_id = ScanCommandQueue::encode(prod_cmd.scan_id);
      std::cout << "[TRACK-CREATE] id=" << prod_id
                << " ms_level=" << group.msn_level << " type=production"
                << std::endl;

      commands.push_back(prod_cmd);
    }
    else
    {
      auto next_cmds = initiateNextLevel(group.msn_level,
          group.variants[best_idx].result, group.faims_cv, queue);
      commands.insert(commands.end(), next_cmds.begin(), next_cmds.end());
    }

    for (const auto& vr : group.variants)
      variant_tracking_map_.erase(queue.decode(vr.tracking_id));

    active_groups_.erase(git);
    return commands;
  }

  std::vector<ScanCommand> Exploration::initiateNextLevel(int msn_level,
      const DeconvolvedSpectrum& result, double faims_cv, ScanCommandQueue& queue)
  {
    std::vector<ScanCommand> commands;

    int next_level = msn_level + 1;
    const auto& this_cfg = config_.level(msn_level);
    const auto& next_cfg = config_.level(next_level);
    if (this_cfg.selection == SelectionMetric::None) return commands;

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

        auto sub_cmds = initiate(next_level, frag_pg, std::abs(charges[ti]), faims_cv, queue);
        commands.insert(commands.end(), sub_cmds.begin(), sub_cmds.end());
      }
    }
    else
    {
      // Direct command building for each fragment target
      for (int ti = 0; ti < num_targets; ++ti)
      {
        // Build as PeakGroup for the buildMS2 factory
        PeakGroup frag_pg(std::abs(charges[ti]), std::abs(charges[ti]), true);
        frag_pg.setMonoisotopicMass(masses[ti]);
        FLASHHelperClasses::LogMzPeak lp;
        lp.mz = (wstarts[ti] + wends[ti]) / 2.0;
        lp.abs_charge = std::abs(charges[ti]);
        frag_pg.push_back(lp);

        ScanCommand cmd = queue.buildMS2(frag_pg, std::abs(charges[ti]), next_scan_config, 1);
        cmd.msn_level = next_level;
        cmd.faims_cv = faims_cv;

        std::string id_str = ScanCommandQueue::encode(cmd.scan_id);
        std::cout << "[TRACK-CREATE] id=" << id_str
                  << " ms_level=" << next_level << " type=next_level"
                  << std::endl;

        commands.push_back(cmd);
      }
    }

    return commands;
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

    // Sum intensity within the precursor isolation window
    double iso_half = group.isolation_width / 2.0;
    double mz_low = group.precursor_mz - iso_half;
    double mz_high = group.precursor_mz + iso_half;

    double remaining_intensity = 0.0;
    for (int i = 0; i < length; ++i)
    {
      if (mzs[i] >= mz_low && mzs[i] <= mz_high)
        remaining_intensity += ints[i];
    }

    // Reference: charge-specific intensity from the precursor PeakGroup
    double reference = group.precursor_pg.getChargeIntensity(
        std::abs(group.precursor_charge));
    if (reference <= 0.0)
      reference = group.precursor_pg.getIntensity();
    if (reference <= 0.0)
      return 0.0;

    double ratio = remaining_intensity / reference;
    // Score = 1 - ratio, clamped to [0, 1]. Higher = less remaining = better fragmentation.
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
