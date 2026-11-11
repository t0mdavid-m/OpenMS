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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/AnnotatorMetrics.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/EventLog.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3YieldModel.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/NotchSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <algorithm>
#include <utility>
#include <sstream>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <optional>

namespace OpenMS
{
  namespace
  {
    /// The chosen protein of a matcher that has not identified anything yet (ADR-0061 push 4b): none, so
    /// IdentificationConfig::proteinOr / accessionOr fall back to the inline reference. A named empty string
    /// rather than a temporary, so the reference proteinOr returns can never bind to a dead argument.
    const std::string kNoChosenProtein;
  }

  Exploration::Exploration(const Config& config, FragmentAnalysis& fragments)
    : config_(config), fragments_(fragments)
  {
    exploration_deconv_ = std::make_unique<Deconvolution>(config, config.explorationToleranceList());
  }

  std::vector<int> Exploration::sweepPermutation(int n, uint64_t seed)
  {
    std::vector<int> perm;
    if (n <= 0) return perm;
    perm.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) perm[static_cast<size_t>(i)] = i;
    uint64_t state = seed;
    auto next = [&state]() {   // splitmix64
      uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
      z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
      z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
      return z ^ (z >> 31);
    };
    for (int i = n - 1; i > 0; --i)
    {
      const int j = static_cast<int>(next() % static_cast<uint64_t>(i + 1));   // modulo bias is irrelevant at sweep sizes
      std::swap(perm[static_cast<size_t>(i)], perm[static_cast<size_t>(j)]);
    }
    return perm;
  }

  std::vector<Exploration::VariantParams> Exploration::buildVariants_(
      const MSLevelConfig& cfg, const ScanConfig& base_config) const
  {
    std::vector<VariantParams> variants;

    // Determine which activations to sweep
    std::vector<std::string> acts = cfg.activations;
    if (acts.empty())
      acts.push_back(base_config.activation);

    // Grid points from min by step, then the top edge if the range does not divide (ADR-0047 decision 8).
    // Index arithmetic, not accumulation: `ce += step` drifts, and the old `ce <= max + 1e-9` loop dropped
    // the top edge of every range the step does not divide (15..40/2 stopped at 39) -- the decimation
    // study whose result this engine implements rests on both edges being acquired.
    auto grid = [](double lo, double hi, double step) {
      std::vector<double> g;
      const int n = static_cast<int>(std::floor((hi - lo) / step + 1e-9)) + 1;
      for (int i = 0; i < n; ++i) g.push_back(lo + i * step);
      if (g.empty() || g.back() < hi - 1e-9) g.push_back(hi);
      return g;
    };
    std::vector<double> ces = grid(cfg.ce_min, cfg.ce_max, cfg.ce_step);
    std::vector<double> rts;
    if (cfg.rt_max > cfg.rt_min) rts = grid(cfg.rt_min, cfg.rt_max, cfg.rt_step);

    for (const auto& act : acts)
    {
      // ONE definition of which axis an activation varies, shared with Config's validator
      // (Config.h). These used to be re-inlined here as string literals, so the engine held two
      // copies of the same rule with nothing pinning them together.
      const bool sweep_ce = needsCollisionEnergy(act);   // HCD | CID | EThcD
      const bool sweep_rt = needsReactionTime(act);      // ETD | EThcD

      // No reference scan heads the block (ADR-0048 superseded ADR-0029): the precursor share is
      // referenced against the scan's own spectrum, so a sweep is exactly its grid points. A grid
      // authored to start at CE 0 acquires that point as an ordinary variant.
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
      const Ms2Params* stage0_params,
      const std::vector<NotchCandidate>* stage1_notches,
      const std::vector<int>* allowed_charges)
  {
    std::vector<ScanCommand> commands;

    const auto& cfg = config_.level(msn_level);
    if (!config_.hasExploration(msn_level)) return commands;

    // Build base ScanConfig from the level's primary scan config, then apply overrides
    ScanConfig base_config = cfg.scans[0];
    base_config.applyOverrides(cfg.overrides);

    auto variant_params = buildVariants_(cfg, base_config);
    if (variant_params.empty()) return commands;

    // NOTHING IS PREPENDED HERE, and nothing is synthesized in buildVariants_ either: a group is exactly
    // its grid points (ADR-0048). The CE-0 reference scan this used to insert -- first here, then per
    // swept activation in buildVariants_ (ADR-0029) -- is gone; the precursor share is referenced
    // against the scan's own total, so no variant needs a denominator scan.

    // Compute precursor_mz and isolation_width from PeakGroup. The interval a metric SUMS must be the
    // interval the instrument was told to ISOLATE (ADR-0026 decision 2). buildMS2 widens the measured
    // range by optimal_window_margin_ per side before commanding it (ScanCommandQueue.cpp:291-295), so
    // an MS2 group's window is (mz2 - mz1) + 2 * margin -- 0.8 Th WIDER than the bare PeakGroup span
    // this used to take. The 2.0 Th minimum applies to MS3+ fragment isolation only (buildMS3 stage[1]
    // re-floors identically at ScanCommandQueue.cpp:451).
    //
    // MS2-ONLY, AND THE MARGIN GOES BEFORE THE FLOOR:
    //   - buildMS3 applies NO margin, and Exploration passes its own isolation_width straight into
    //     buildMS3. Widening at MS3 would therefore RE-CREATE the same disagreement in the other
    //     direction *and* change the commanded MS3 width, moving three exploration golden modes.
    //   - FLASHIda_exploration_test.cpp:2926 pins TEST_REAL_SIMILAR(group.isolation_width, 2.0) on the
    //     MS3 floor. Applying the margin AFTER the floor gives 2.8 and breaks it; applying it BEFORE
    //     gives max(0.8, 2.0) == 2.0 -- harmless but pointless. MS2-only is the correct scope.
    //   - The CENTRE needs no correction: buildMS2 computes (mz1 + mz2) / 2 BEFORE applying the margin
    //     (ScanCommandQueue.cpp:292), so the two sides already agree on precursor_mz.
    //
    // The margin is also what keeps a charge resolved at a SINGLE isotope (mz2 == mz1) from handing the
    // metric a zero-width window to sum over. Nothing downstream turns this width into a scan range any
    // more -- a pre-scan reads out its level's configured range (ADR-0044).
    auto [mz1, mz2] = pg.getMzRange(charge);
    double precursor_mz = (mz1 + mz2) / 2.0;
    double isolation_width = (msn_level >= 3)
                               ? std::max(mz2 - mz1, 2.0)
                               : (mz2 - mz1) + 2.0 * optimal_window_margin_;
    double precursor_mass = pg.getMonoMass();

    ExplorationGroup group;
    group.group_id = next_group_id_++;
    group.msn_level = msn_level;

    // exploration.order: under "random" this sweep acquires its whole grid -- every activation and both axes -- in its own
    // order, drawn from the group id so it is reproducible; grid_positions keeps each variant's place in the grid for the
    // tie-break. Under "ascending" the acquisition order is the grid order, as before.
    std::vector<int> grid_positions(variant_params.size());
    for (size_t gi = 0; gi < grid_positions.size(); ++gi) grid_positions[gi] = static_cast<int>(gi);
    if (cfg.sweep_order == SweepOrder::Random && variant_params.size() > 1)
    {
      const std::vector<int> perm = sweepPermutation(static_cast<int>(variant_params.size()),
                                                     kSweepSeed ^ static_cast<uint64_t>(group.group_id));
      std::vector<VariantParams> shuffled;
      shuffled.reserve(variant_params.size());
      for (const int p : perm) shuffled.push_back(variant_params[static_cast<size_t>(p)]);
      variant_params = std::move(shuffled);
      grid_positions = perm;
    }
    // An unassigned exhaustive target (unknown ion class) cannot be scored by a READING metric --
    // the matcher refuses to project without an ion frame, so every variant scores 0 and the winner
    // is a coin flip. RemainingPrecursor scores from isolation-window intensity alone, needs no
    // proteoform, and being a MEASURING metric it also earns the ADR-0020 close-out production scan
    // with no further edit. ADR-0023 decision 11.
    //
    // msn_level >= 3 is LOAD-BEARING: initiate()'s ion_type parameter defaults to '\0' and the MS2
    // call site (FLASHIda.cpp) passes nothing, so an ungated test would force this metric onto EVERY
    // MS2 exploration group and move four committed golden modes.
    //
    // Keyed on the ion CLASS, never on the 'u' literal: every projection site refuses on the class
    // (MS3FragmentMatcher::extractSubsequence, calibrateAndScore), the sentinel's value is TU-local
    // to ProteoformTracker.cpp's anonymous namespace and invisible from here, and any FUTURE
    // unprojectable class must take this branch too.
    //
    // '\0' IS EXEMPT, AND THAT EXEMPTION IS THE POINT. It means "no ion identity was RECORDED", not
    // "an ion identity the matcher refuses" -- the same asymmetry planExhaustive_ spells out for an
    // empty activation_type ("NOT EMPTY *AND* NOT CAPABLE => refuse"; failing closed on the empty
    // value is the quieter bug). Production never reaches MS3 with it: Ms3Target::ion_type is
    // non-empty at both assignment sites and the conversion in initiateNextLevel is the only
    // string->char site. Only unit tests calling initiate(3, ...) without an ion pass it, and
    // forcing the metric there turned two MS3 sections into empty-baseline aborts -- one reported as
    // a score of 0, the other as a SegFault that hid 44 of the file's 48 sections.
    const bool has_ion_identity  = (ion_type != '\0');
    const bool unprojectable_ion = has_ion_identity && !MS3FragmentMatcher::isKnownIonClass(ion_type);
    // Every metric but kinetic_model: an unassigned mass cannot be READ (no ion identity to match), so it
    // is measured by remaining_precursor instead (ADR-0023 d.11). kinetic_model already measures -- it
    // fits the precursor share, which a 'u' target carries like any other -- and forcing it onto
    // remaining_precursor would silently drop the authored k2 for every unassigned mass of an exhaustive
    // run, the stomped-key shape ADR-0047 d.6 forbids. So it keeps its metric.
    group.exploration_metric = (msn_level >= 3 && unprojectable_ion && cfg.exploration != ExplorationMetric::KineticModel)
                                 ? ExplorationMetric::RemainingPrecursor
                                 : cfg.exploration;
    group.precursor_mz = precursor_mz;
    group.precursor_mass = precursor_mass;
    group.isolation_width = isolation_width;
    group.precursor_charge = charge;
    group.precursor_pg = pg;
    if (allowed_charges != nullptr) group.allowed_charges = *allowed_charges;
    group.faims_cv = (ms_ctx != nullptr) ? ms_ctx->faims_cv : 0.0;  // Item 1: CV travels via the context
    group.start_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    if (ms_ctx != nullptr) group.originating_cmd = *ms_ctx;
    group.fragment_ion_type = ion_type;
    group.fragment_ion_index = frag_index;
    group.proteoform_ctx = proto_ctx;
    // Captured for the post-sweep production scan, which rebuilds from the winning VARIANT rather
    // than from the Ms3Target and so has no other route to the notch set (ADR-0016).
    if (stage1_notches != nullptr) group.stage1_notches = *stage1_notches;

    // -1 is the baseline marker in the logged variant_index column, and there may now be one per
    // swept activation rather than exactly one per group. Real variants keep a single running
    // 0-based counter across the whole group, so a single-activation sweep -- every committed
    // config -- still numbers 0,1,2,... exactly as the old `i - 1` did.
    for (int i = 0; i < static_cast<int>(variant_params.size()); ++i)
    {
      const auto& vp = variant_params[i];
      ExplorationVariant v;
      v.variant_index = i;   // one running counter; no reference scan to skip (ADR-0048); the acquisition position
      v.grid_position = grid_positions[static_cast<size_t>(i)];
      v.collision_energy = vp.collision_energy;
      v.reaction_time = vp.reaction_time;
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
                             stage0_params, proto_ctx, stage1_notches);
      }
      else
      {
        cmd = queue.buildMS2(pg, charge, variant_config, expl_priority, ms_ctx ? ms_ctx->scan_id : 0,
                             &group.allowed_charges);
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

      EventLine() << "[TRACK-CREATE] id=" << id_str
                << " ms_level=" << msn_level << " type=exploration"
                << " activation=" << vp.activation
                << " CE=" << vp.collision_energy
                << " RT=" << vp.reaction_time;
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
    const VariantInput in = prepareVariant(tracking_id, mzs, ints, length, rt);
    if (!in.known) return {};
    return applyVariant(in, nullptr, false, queue, tracker, precursor_id);
  }

  Exploration::VariantInput Exploration::prepareVariant(int tracking_id, const double* mzs, const double* ints,
                                                        int length, double rt)
  {
    VariantInput in;
    in.tracking_id = tracking_id;

    // Look up the group to get correct precursor context for deconvolution
    auto vit = variant_tracking_map_.find(tracking_id);
    if (vit == variant_tracking_map_.end()) return in;

    int group_id = vit->second.group_id;
    auto git = active_groups_.find(group_id);
    if (git == active_groups_.end()) return in;
    const ExplorationGroup& group = git->second;

    // Charge ceiling = the highest charge THIS VARIANT ACTUALLY ISOLATED, not the group's anchor
    // (ADR-0016). With a co-isolated charge set every member is genuinely present in the isolation, so
    // a fragment of the highest member may itself carry that charge; capping at the anchor silently
    // discards those fragments before matching -- and since MS3 targets are chosen from the matched
    // fragment list, that thins the target set with nothing anywhere reading wrong.
    //
    // This mirrors the two regular return paths, which already call maxIsolatedCharge
    // (FLASHIda.cpp:330 for MS2, :526 for MS3); the exploration path was the one that did not, and it
    // is the only deconvolveMSn in this file -- it serves MS2 *and* MS3 variants. The LAST cascade
    // stage is the one that fragmented, so stage 0 for an MS2 variant and stage 1 for an MS3 variant.
    //
    // The mass is deliberately NOT touched: every co-isolated member is the same neutral mass, which
    // is exactly why the spectrum is not chimeric. Only the charge ceiling was ever charge-specific.
    // With no notches this resolves to abs(anchor charge), i.e. the previous value, so every existing
    // golden stays byte-identical.
    int precursor_charge = group.precursor_charge;
    // ADR-0045: a trap pre-scan is measured, never identified. Decided on the variant's OWN command --
    // the engine's record of what it asked for -- never on anything the returning scan carries
    // (ADR-0042). An out-of-range slot leaves it false, the same fallthrough precursor_charge takes.
    bool measured_only = false;
    const int variant_slot = vit->second.variant_index;   // array index into group.variants,
                                                          // NOT ExplorationVariant::variant_index
                                                          // (which is -1 for the CE-0 baseline)
    if (variant_slot >= 0 && variant_slot < static_cast<int>(group.variants.size()))
    {
      const ScanCommand& vcmd = group.variants[variant_slot].cmd;
      if (vcmd.num_stages > 0) precursor_charge = maxIsolatedCharge(vcmd, vcmd.num_stages - 1);
      measured_only = isTrapAnalyzer(vcmd.analyzer);
      in.activation_type = group.variants[variant_slot].activation_type;
    }

    // Deconvolve with correct precursor context from the group -- unless the variant is a trap
    // pre-scan, which is never deconvolved: unit-resolution peaks have no business in a ppm-tolerance
    // identification, and the only thing read from such a scan is the raw window sum.
    in.deconv = DeconvolvedSpectrum(tracking_id);
    if (!measured_only && mzs != nullptr && ints != nullptr && length > 0)
    {
      exploration_deconv_->deconvolveMSn(mzs, ints, length, rt,
                                         group.precursor_mass, precursor_charge);
      in.deconv = exploration_deconv_->storedMS2();
      in.deconvolved = true;
    }
    last_result_deconvolved_ = !measured_only;

    in.known = true;
    in.msn_level = group.msn_level;
    in.measured_only = measured_only;
    // The precursor SHARE is read for every variant under every metric (ADR-0048 d.1): it is what the
    // share-based metrics score, what the Kinetic Model fits (ADR-0047), and what scan_results logs.
    // No reference scan exists; the share is the window's fraction of the scan's own signal. Read here,
    // where the raw arrays still exist; applyVariant stores it on the variant.
    in.precursor_share = precursorShare_(group, mzs, ints, length);
    return in;
  }

  bool Exploration::variantNeedsSearch(const VariantInput& in) const
  {
    // An empty spectrum has nothing to search; a trap variant is never matched (ADR-0045); an MS3 variant is not
    // an MS2 search (ADR-0061 d.2) and stays on the analysis thread. Any reference is searched, a FASTA as well as
    // an inline sequence (ADR-0061 push 4c).
    return in.known && in.msn_level == 2 && !in.measured_only && in.deconvolved && !in.deconv.empty()
           && config_.identification().enabled();
  }

  FragmentAnalysis::ProteoformMatch Exploration::searchVariant(const VariantInput& in) const
  {
    // The worker set in.prsm (ADR-0061 push 4c); without it, the inline reference, as before the search.
    return computeFragmentMatch_(in.deconv, in.msn_level, in.activation_type,
                                 config_.identification().proteinOr(kNoChosenProtein),
                                 config_.identification().accessionOr(kNoChosenProtein), in.prsm ? &*in.prsm : nullptr);
  }

  Exploration::FeedResultInfo Exploration::feedResultImpl_(int tracking_id,
      const DeconvolvedSpectrum& msn_deconv,
      const double* mzs, const double* ints, int length,
      ScanCommandQueue& queue, ProteoformTracker* tracker, int precursor_id, bool measured_only)
  {
    // The deconvolution-bypass entry (ExplorationTestAccess): the caller supplies the spectrum. The share
    // still needs the group's window, so it is read here exactly as prepareVariant reads it.
    VariantInput in;
    in.tracking_id = tracking_id;
    in.measured_only = measured_only;
    in.deconv = msn_deconv;
    auto rit = variant_tracking_map_.find(tracking_id);
    if (rit != variant_tracking_map_.end())
    {
      auto git = active_groups_.find(rit->second.group_id);
      if (git != active_groups_.end()) in.precursor_share = precursorShare_(git->second, mzs, ints, length);
    }
    return applyVariant(in, nullptr, false, queue, tracker, precursor_id);
  }

  Exploration::FeedResultInfo Exploration::applyVariant(const VariantInput& in,
      const FragmentAnalysis::ProteoformMatch* searched, bool dropped,
      ScanCommandQueue& queue, ProteoformTracker* tracker, int precursor_id)
  {
    const int tracking_id = in.tracking_id;
    const DeconvolvedSpectrum& msn_deconv = in.deconv;
    const bool measured_only = in.measured_only;
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
    v.measured_only = measured_only;
    v.precursor_share = in.precursor_share;   // read by prepareVariant (ADR-0048 d.1)
    // The variant keeps its search for the cascade (ADR-0061 push 4c). A dropped search keeps an empty one: it
    // identifies nothing and gives characterization nothing (ADR-0062 d.5), at the next level too.
    if (dropped)
      v.prsm = FragmentAnalysis::SearchedPrSM{};
    else
      v.prsm = in.prsm;
    FragmentAnalysis::ProteoformMatch frag{};
    const FragmentAnalysis::ProteoformMatch no_match{};
    if (measured_only)
    {
      // Window and total only (ADR-0045). The metric is share-based by config -- Config::validate
      // refuses a trap sweep under a counting metric -- so the -1.0 arm is unreachable from config and
      // fails closed: a variant with nothing to count must not win at 0 (winner selection seeds -1.0,
      // strictly greater). Under the Kinetic Model no variant is scored at all: the group is decided
      // by the fit at completion.
      v.score = (group.exploration_metric == ExplorationMetric::RemainingPrecursor)
                  ? computeRemainingPrecursorScore_(group, v.precursor_share)
                  : -1.0;
    }
    else if (dropped)
    {
      // ADR-0062 d.5: the backlog dropped this variant's search, so it has no match. A count it never took
      // must not compete, so fragment_count leaves it unscored (-1.0, which never wins); every other metric
      // scores from what was read under the lock, as for any variant. It is not pooled (below).
      v.score = (group.exploration_metric == ExplorationMetric::FragmentCount || isAnnotatorMetric(group.exploration_metric))
                  ? -1.0   // ADR-0074: an annotator metric reads the PrSM the dropped search never made
                  : computeExplorationScore_(group.exploration_metric, msn_deconv, group, v.precursor_share,
                                             &frag, v.activation_type, &no_match);
    }
    else
    {
      v.score = computeExplorationScore_(group.exploration_metric, msn_deconv, group, v.precursor_share,
                                         &frag, v.activation_type, searched);
    }
    v.tic_coverage = measured_only ? 0.0f : computeTICCoverage_(msn_deconv);
    v.fragment_count = frag.total_match_count;
    v.received = true;

    const int real_variant_count = static_cast<int>(group.variants.size());

    // Populate per-variant metadata in the return info
    info.group.group_id = group.group_id;
    info.group.variant_index = v.variant_index;
    info.group.total_variants = real_variant_count;
    info.fragmentation.collision_energy = v.collision_energy;
    info.fragmentation.activation_type = v.activation_type;
    info.fragmentation.reaction_time = v.reaction_time;
    // Two-level only by construction: FragmentationInfo carries exactly one stage0_* triplet, so this
    // records "the isolation stage before the fragmentation stage" and nothing deeper. Generalizing to
    // MS4+ means a per-stage vector here AND in the scan_results CE/activation/reaction columns, which
    // are emitted as a fixed "stage0;stage1" pair (IdaLogger) and golden-locked in that shape. Deliberately
    // not built for: there is no MS4 path to validate it against. The other wall is the `msn_level < 3`
    // cascade test in feedResultImpl_.
    if (group.msn_level >= 3 && v.cmd.num_stages >= 2)
    {
      info.fragmentation.stage0_collision_energy = v.cmd.stages[0].collision_energy;
      info.fragmentation.stage0_activation_type  = std::string(v.cmd.stages[0].activation_type);
      info.fragmentation.stage0_reaction_time    = v.cmd.stages[0].reaction_time;
    }
    info.metric.score = v.score;
    info.metric.tic_coverage = v.tic_coverage;
    info.metric.fragment_count = v.fragment_count;
    info.metric.exploration_metric = static_cast<int>(group.exploration_metric);
    info.identification.matched_protein = frag.matched_protein;
    info.identification.proteoform_sequence = frag.proteoform_sequence;
    // MS2 only. At MS3 `frag` is the whole-protein MS2 matcher's match of a sub-fragment spectrum: it has no
    // equivalent ion and no adjusted mass, so as an identification row it read "0;0;0" / "0.0000". An MS3
    // variant's identification is the calibrated re-score's alone (below, fragment_count); a measuring metric
    // identifies no pre-scan (ADR-0020), and a variant of a dropped sweep is never scored.
    if (group.msn_level == 2) info.identification.result = frag;

    // Stage this MS2 variant into the ProteoformTracker for pooling (MS3 evidence folds via foldMs3,
    // not here). A measured-only (trap) variant is never pooled (ADR-0045), nor one whose search the
    // backlog dropped (ADR-0062 d.5).
    if (tracker != nullptr && group.msn_level == 2 && !v.measured_only && !dropped)
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
    
    info.metric.precursor_share = v.precursor_share;
    std::string parent_enc = ScanCommandQueue::encode(group.originating_cmd.scan_id);
    std::strncpy(info.parent_scan_id, parent_enc.c_str(), 3);
    info.parent_scan_id[3] = '\0';
    
    // Renders the group's MS2Context for one variant. Reads group.proteoform_ctx (the TRIGGERING scan's
    // render context, not the winner) and the variant's own command, so it is correct at any point in
    // feedResult; it is a lambda rather than a member because it closes over `group`. Note the name is
    // literal: the value describes the MS2 that produced the group, even when the group is at MS3.
    auto buildMS2ContextForVariant = [&](int group_variant_index)
    {
      MS2Context ctx;
      // With no protein chosen and none inline (a FASTA reference, an unidentified Precursor) there is nothing to cut.
      const std::string& render_protein = config_.identification().proteinOr(group.proteoform_ctx.protein_sequence);
      ctx.proteoform_sequence = (group.proteoform_ctx.region_start >= 0 && !render_protein.empty())
          ? render_protein.substr(group.proteoform_ctx.region_start,
              group.proteoform_ctx.region_end - group.proteoform_ctx.region_start)
          : "";
      ctx.protein_sequence = group.proteoform_ctx.protein_sequence;
      ctx.protein_accession = group.proteoform_ctx.protein_accession;
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

    info.identification.ms2_context = buildMS2ContextForVariant(variant_index);

    if (group.msn_level >= 3 && !info.identification.ms2_context.proteoform_sequence.empty())
    {
      info.identification.proteoform_sequence = info.identification.ms2_context.proteoform_sequence;
      info.identification.matched_protein = config_.identification().accessionOr(info.identification.ms2_context.protein_accession);
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

    // MS3 FragmentCount: batch re-score with calibrated subsequence matching. An annotator metric (ADR-0074 d.3, d.6)
    // runs the same batch for the variants' identification -- it is a reading metric, the winner is its own
    // evidence -- and is then scored by the annotator below, not by the match count.
    const bool annotator_metric = isAnnotatorMetric(group.exploration_metric);
    if ((group.exploration_metric == ExplorationMetric::FragmentCount || annotator_metric) && group.msn_level >= 3)
    {
      std::vector<const DeconvolvedSpectrum*> variant_spectra;
      for (auto& var : group.variants)
        // A measured-only variant has no spectrum to score; nullptr, so it is excluded rather than scored 0
        // (unreachable from config -- a trap sweep cannot be FragmentCount -- but not from a direct initiate).
        variant_spectra.push_back(var.received && !var.measured_only ? &var.result : nullptr);

      // Score the completed exploration group against the LIVE WINNER ("Full"), not the render context.
      // group.proteoform_ctx stays the triggering-scan RENDER context (drives buildMS3 + the MS2Context);
      // tracker is null in unit tests -> fall back to the render context (the only one available then).
      const MS3FragmentMatcher::ProteoformContext score_ctx = (tracker != nullptr)
          ? tracker->buildWinnerProteoformContext(precursor_id)
          : group.proteoform_ctx;
      std::vector<FragmentAnalysis::ProteoformMatch> detailed_results;
      auto calibrated_scores = ProteoformTracker::scoreCalibratedVariants(
        variant_spectra,
        config_.identification().proteinOr(score_ctx.protein_sequence),
        score_ctx,
        group.fragment_ion_type,
        group.fragment_ion_index,
        MS3FragmentMatcher::LOOSE_TOLERANCE_PPM,
        config_.level(group.msn_level).tolerance_ppm,
        &detailed_results,
        config_.identification().fixed_mods);

      for (size_t vi = 0; vi < calibrated_scores.size(); ++vi)
      {
        if (!annotator_metric) group.variants[vi].score = calibrated_scores[vi];
        group.variants[vi].fragment_count = static_cast<int>(calibrated_scores[vi]);
        if (vi < detailed_results.size())
          group.variants[vi].identification_result = detailed_results[vi];
      }

      if (annotator_metric)
      {
        // ADR-0074 d.3-d.5, d.7: every variant's raw peaks annotated against the target fragment as a chain of its
        // own (a b target carries its water loss), with the ion types of the variant's own MS3 activation. No PrSM
        // (a c/z/a/x target, or no identified protein): every variant keeps -1.0 and the group ends no-winner.
        FLASHAnnotatorAlgorithm::PrSM target;
        if (AnnotatorMetrics::ms3TargetPrSM(config_.identification().proteinOr(score_ctx.protein_sequence), score_ctx,
                                            group.fragment_ion_type, group.fragment_ion_index, group.precursor_charge, target))
        {
          FLASHAnnotatorAlgorithm annotator;
          std::string annotator_activation;
          bool annotator_set = false;
          for (auto& var : group.variants)
          {
            if (!var.received || var.measured_only) continue;
            if (!annotator_set || var.activation_type != annotator_activation)
            {
              annotator.setParameters(AnnotatorMetrics::parameters(FragmentAnalysis::getIonTypesForFragmentationMethod(var.activation_type),
                                                                   config_.identification().flashtnt.fixed_mod,
                                                                   config_.level(group.msn_level).tolerance_ppm));
              annotator_activation = var.activation_type;
              annotator_set = true;
            }
            const AnnotatorMetrics::Values values = AnnotatorMetrics::annotate(annotator, var.result.getOriginalSpectrum(), target);
            var.identification_result.raw_score = values.raw_score;
            var.identification_result.bond_coverage = values.bond_coverage;
            var.identification_result.explained_signal = values.explained_signal;
            var.score = group.exploration_metric == ExplorationMetric::RawScore ? values.raw_score
                        : group.exploration_metric == ExplorationMetric::BondCoverage ? values.bond_coverage
                        : values.explained_signal;
          }
        }
      }

      // Re-populate identification_result after batch re-scoring updated it
      info.identification.result = group.variants[variant_index].identification_result;
      info.metric.score = group.variants[variant_index].score;

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
        row.flash_extender_score = gv.identification_result.score;   // C: self-ID row carries its own score
        info.identification.additional_rows.push_back(row);
      }
    }

    // The group's decision. Two shapes (ADR-0047):
    //   * every metric but kinetic_model picks a WINNER -- strict `> best_score` from -1.0, so a variant
    //     that was not scored (-1.0: no share could be read, or a trap pre-scan under a counting metric)
    //     can never win; if nothing was scored, best_idx stays -1 and the no-winner arm handles it;
    //   * kinetic_model fits the Kinetic Model to the shares and takes the MODELED MAXIMUM -- an energy no
    //     pre-scan need have run at, so no variant wins and a follow-up is always acquired.
    const bool kinetic = (group.exploration_metric == ExplorationMetric::KineticModel);
    int best_idx = -1;
    double best_score = -1.0;
    if (kinetic)
    {
      std::vector<KineticModel::SharePoint> pts;
      for (const auto& var : group.variants) pts.push_back({var.collision_energy, var.precursor_share});
      group.kinetic_fit = KineticModel::fitGompertz(pts);
      const auto& lvl = config_.level(group.msn_level);
      // One activation per kinetic sweep (Config::validate), so the first variant's activation is the sweep's.
      const auto k2_it = lvl.kinetic_k2.find(group.variants.front().activation_type);
      group.kinetic_k2 = (k2_it != lvl.kinetic_k2.end()) ? k2_it->second : -1.0;
      group.kinetic_max = KineticModel::modeledMaximum(group.kinetic_fit, group.kinetic_k2,
                                                       static_cast<int>(std::lround(lvl.ce_max)));
      info.kinetic.fitted = group.kinetic_fit.ok;
      info.kinetic.fit_points = group.kinetic_fit.n_points;
      info.kinetic.p0 = group.kinetic_fit.p0;
      info.kinetic.ce_max_s = group.kinetic_fit.ce_max_s;
      info.kinetic.width = group.kinetic_fit.width;
      info.kinetic.rss = group.kinetic_fit.rss;
      info.kinetic.k2 = group.kinetic_k2;
      info.kinetic.modeled_ce = group.kinetic_max.ce;
      info.kinetic.modeled_b = group.kinetic_max.b;
      info.kinetic.at_lower_edge = group.kinetic_max.at_lower_edge;
      info.kinetic.at_upper_edge = group.kinetic_max.at_upper_edge;
      info.kinetic.energy_source = group.kinetic_max.ce <= 0 ? "none"
                                 : group.kinetic_max.at_lower_edge ? "modeled_edge_lower"
                                 : group.kinetic_max.at_upper_edge ? "modeled_edge_upper" : "modeled";
    }
    else
    {
      // A tie goes to the earlier GRID position (exploration.order): under "ascending" that is the earlier-acquired
      // variant, as before; under "random" the same scores pick the same winner whatever order they were acquired in.
      // best_idx >= 0 keeps an unscored variant (-1) from ever winning, as the strict > from -1 always has.
      for (int i = 0; i < static_cast<int>(group.variants.size()); ++i)
      {
        const auto& var = group.variants[i];
        if (var.score > best_score
            || (best_idx >= 0 && var.score == best_score && var.grid_position < group.variants[best_idx].grid_position))
        {
          best_score = var.score;
          best_idx = i;
        }
      }
      info.kinetic.energy_source = best_idx >= 0 ? "measured" : "none";
    }

    // Finalize the proteoform model before the next-level dispatch below
    if (tracker != nullptr && group.msn_level == 2)
    {
      tracker->finalizeMS2(precursor_id);
    }

    const bool decided = kinetic ? (group.kinetic_max.ce > 0) : (best_idx >= 0);
    if (!decided)
    {
      // No decision: no follow-up scan, and at MS3 the target is not acquired at this level (ADR-0047 d.7:
      // no fallback to the authored energy -- an energy the model could not justify is not acquired in its
      // name). Falls through to the cleanup below.
      group.complete = true;
      EventLine() << "[EXPL-ABORT] group=" << group.group_id
                << " reason=" << (kinetic ? "no-fit" : "no-winner") << " winner=none";
    }
    else
    {
      // ref_idx: the variant whose command stands in for "the sweep's command" when the follow-up is
      // rebuilt (stage-0 context, fragment scores, window_snr). Under kinetic_model no variant won, so the
      // first one serves -- every variant of a group shares that context.
      const int ref_idx = kinetic ? 0 : best_idx;
      const int chosen_ce = kinetic ? group.kinetic_max.ce
                                    : static_cast<int>(group.variants[best_idx].collision_energy);
      group.complete = true;
      if (!kinetic)
      {
        group.winner_index = best_idx;
        group.variants[best_idx].result.getOrCreateOptimizationMetadata().is_best_variant = true;
        info.group.winner_tracking_id = group.variants[best_idx].tracking_id;
      }
      else
      {
        EventLine() << "[EXPL-MODELED] group=" << group.group_id
                  << " activation=" << group.variants[0].activation_type
                  << " p0=" << group.kinetic_fit.p0 << " ce_max_s=" << group.kinetic_fit.ce_max_s
                  << " w=" << group.kinetic_fit.width << " k2=" << group.kinetic_k2
                  << " points=" << group.kinetic_fit.n_points << " rss=" << group.kinetic_fit.rss
                  << " CE=" << chosen_ce << " B=" << group.kinetic_max.b
                  << (group.kinetic_max.at_lower_edge ? " edge=lower" : group.kinetic_max.at_upper_edge ? " edge=upper" : "");
      }

      // C: emit identification.tsv rows for NON-WINNER variants that did NOT self-identify (empty own FLASHTnT
      // match) but whose masses re-matched the WINNER ladder in finalizeMS2 (mapNonWinnerMs2_) — they fold into
      // the pool (contributing_scan_ids) yet would otherwise have no id row. Render the winner proteoform + the
      // scan's re-matched fragments, flash_extender_score = -1 (no own ID = the distinguisher). Not under the
      // Kinetic Model: no exploration winner exists there; the follow-up is identified on the ordinary path.
      if (!kinetic && tracker != nullptr && group.msn_level == 2)
      {
        for (size_t vi = 0; vi < group.variants.size(); ++vi)
        {
          if (static_cast<int>(vi) == best_idx) continue;               // winner already has its own self-ID row
          const auto& gv = group.variants[vi];
          if (gv.measured_only) continue;                               // a trap pre-scan is never pooled (ADR-0045)
          if (!gv.identification_result.fragments.empty()) continue;    // self-identified -> already logged
          const FragmentAnalysis::ProteoformMatch* rem =
            tracker->getRematchedNonWinnerMatch(precursor_id, gv.cmd.scan_id);
          if (rem == nullptr) continue;                                 // no winner-ladder re-match -> no row
          FeedResultInfo::IdentificationRowInfo row;
          row.tracking_id = gv.tracking_id;
          row.identification_result = *rem;                             // winner proteoform + re-matched fragments
          row.ms2_context = buildMS2ContextForVariant(static_cast<int>(vi));
          row.tic_coverage = gv.tic_coverage;
          row.flash_extender_score = -1.0;                             // no own ID
          info.identification.additional_rows.push_back(row);
        }
      }

      if (!kinetic)
        EventLine() << "[EXPL-WINNER] group=" << group.group_id
                  << " winner_idx=" << best_idx
                  << " activation=" << group.variants[best_idx].activation_type
                  << " CE=" << group.variants[best_idx].collision_energy
                  << " RT=" << group.variants[best_idx].reaction_time
                  << " score=" << best_score;

      const auto& level_config = config_.level(group.msn_level);

      // Send production scan after exploration. Two independent reasons to re-acquire (ADR-0020):
      //
      //  1. overrides non-empty -- the pre-scans ran at DEGRADED settings (that is what overrides are
      //     for), so no variant was acquired at production fidelity. Applies at every level.
      //  2. a MEASURING metric at MS3 -- MassCount/RemainingPrecursor score from bulk signal and never
      //     match fragments, so the whole group leaves NO identification even though every variant was
      //     acquired at production settings. FragmentCount is exempt: it matches each variant to score
      //     it, so the batch re-score below already populated identification_result and the inline fold
      //     at the tail of this function carries the evidence.
      //
      // Deliberately MS3-only: an MS2 variant is matched under every metric (computeFragmentMatch_ runs
      // the whole-protein matcher, which is the right one at MS2), so an MS2 sweep already has its
      // evidence and cascades to MS3 instead of re-acquiring itself.
      //
      //  3. the Kinetic Model, at EVERY level (ADR-0047 d.5): its chosen energy is computed, and in general
      //     no pre-scan ran at it, so there is nothing to cascade from. One rule, one code path -- a
      //     pre-scan that happens to sit at the modeled energy is not reused.
      const bool measuring_ms3_sweep = group.msn_level >= 3 && isMeasuringMetric(group.exploration_metric);
      if (kinetic || !level_config.overrides.empty() || measuring_ms3_sweep)
      {
        ScanConfig prod_config = level_config.scans[0];
        prod_config.collision_energy = chosen_ce;
        prod_config.reaction_time = group.variants[ref_idx].reaction_time;
        prod_config.activation = group.variants[ref_idx].activation_type;

        ScanCommand prod_cmd;
        if (group.msn_level >= 3)
        {
          // Annotation seed = the finalized model (matches pooled_identification), not the MS2 group's
          // triggering render context. Self-guards on a finalized winner (empty ctx -> region_start < 0).
          if (tracker != nullptr)
          {
            MS3FragmentMatcher::ProteoformContext win_ctx = tracker->buildWinnerProteoformContext(precursor_id);
            if (win_ctx.region_start >= 0) group.proteoform_ctx = win_ctx;
          }
          // Reconstruct the winning variant's stage-1 fragment scores for the production MS3 (buildMS3
          // keeps its explicit frag_scores param; the other callers pass genuine scores, so do not fold
          // this into buildMS3 — an ms2-command would carry zero-initialised *_s1 fields).
          const ScanCommand& wcmd = group.variants[ref_idx].cmd;
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
          // group.stage1_notches, not the winning variant's own notches: this rebuilds the fragment
          // stage from the group's descriptors, and the variant command's notch slots would be read
          // as stage-0's inherited set by buildMS3's own inheritance step.
          prod_cmd = queue.buildMS3(group.variants[ref_idx].cmd, prod_config,
                                     group.precursor_mz, group.precursor_charge,
                                     group.isolation_width, group.variants[ref_idx].cmd.scan_id,
                                     group.fragment_ion_type, group.fragment_ion_index, 1, wfs,
                                     nullptr, group.proteoform_ctx, &group.stage1_notches);
        }
        else
        {
          prod_cmd = queue.buildMS2(group.precursor_pg, group.precursor_charge, prod_config, 2, 0,
                                    &group.allowed_charges);
        }

        prod_cmd.faims_cv = group.faims_cv;
        {
          const double inherited = group.variants[ref_idx].cmd.window_snr;
          if (inherited >= 0.0) prod_cmd.window_snr = inherited;
        }

        std::string prod_id = ScanCommandQueue::encode(prod_cmd.scan_id);
        EventLine() << "[TRACK-CREATE] id=" << prod_id
                  << " ms_level=" << group.msn_level << " type=production";

        info.commands.push_back(prod_cmd);
        if (group.msn_level >= 3)
          info.ms2_context_cache.emplace_back(prod_cmd.scan_id, buildMS2ContextForVariant(ref_idx));
      }
      // Cascade to the next level. The `< 3` is the MS4 wall: it says "MS2 groups cascade, MS3 groups are
      // terminal" where the general rule would be "cascade while the next level has a selection criterion"
      // -- which initiateNextLevel already checks for itself (next_cfg.selection == None returns empty).
      //
      // That self-check is no longer a CHEAP early return. Since identification was decoupled from the MS3
      // switch, initiateNextLevel matches fragments before it decides it has no commands to build, so with
      // characterization.mode off this call runs a tagger + extender pass whose result is discarded here
      // (only .commands and .ms3_contexts are read below). It is redundant rather than wrong -- every
      // variant of this group was already identified per-variant by computeFragmentMatch_. The condition is
      // deliberately NOT duplicated at this call site: one owner of the dispatch decision beats two that can
      // drift apart, which is how acquisition geometry once ended up sourced from an exclusion flag.
      // Dropping the literal would make this level-agnostic here, but the levels beneath it are not:
      // config_ only ever materializes {1,2,3}, buildMS3 hardcodes a 2-stage command, and the logging
      // emits a fixed stage0;stage1 pair. Not generalized -- no MS4 path exists to validate it.
      else if (group.msn_level < 3)
      {
        // The winner's own search (ADR-0061 push 4c): the next level is planned from the PrSM its worker found, so
        // no proteome search runs here, under the analysis lock. No search ran (the inline path): the inline matcher.
        const ExplorationVariant& winner = group.variants[best_idx];
        auto next_nlr = initiateNextLevel(group.msn_level,
            winner.result, group.faims_cv, queue,
            &winner.cmd, tracker, precursor_id,  // P5: same precursor_id model key
            winner.prsm ? &*winner.prsm : nullptr);
        info.commands.insert(info.commands.end(), next_nlr.commands.begin(), next_nlr.commands.end());
        // Carry the per-command MS2 contexts out ALONGSIDE the commands. Without this the MS3s this
        // branch dispatches are acquired but never identified: a fixed-CE MS3 returns on the REGULAR
        // MS3 path, which resolves its parent MS2 context out of FLASHIda::ms2_context_cache_ and,
        // on a miss, silently skips the whole identification block -- no identification.tsv row, no
        // tracker feedScan/foldMs3, so the MS3 never reaches pooled_identification either. The
        // regular MS2->MS3 path already seeds that cache from NextLevelResult::ms3_contexts
        // (FLASHIda.cpp); mirror it here so an MS2-exploration winner behaves identically.
        // ms3_contexts is index-parallel to next_nlr.commands, but info.commands may already hold
        // other entries, so key by scan_id (what the cache is keyed on) rather than by position.
        for (size_t i = 0; i < next_nlr.commands.size() && i < next_nlr.ms3_contexts.size(); ++i)
        {
          if (next_nlr.commands[i].msn_level >= 3)
            info.ms2_context_cache.emplace_back(next_nlr.commands[i].scan_id, next_nlr.ms3_contexts[i]);
        }
        // ADR-0066: a governed dispatch returns its plan instead of commands; the caller swaps it in.
        info.governed = std::move(next_nlr.governed);
      }
      // overrides empty && msn_level >= 3 && a READING metric (FragmentCount): an MS3 exploration group
      // with no production scan. Since ADR-0020 this is the ONLY way to reach here at MS3 -- a measuring
      // metric now takes the production branch above -- and it is exactly the case where the group already
      // holds its own evidence, because the FragmentCount batch re-score populated identification_result.
      // Nothing else will fold that into the trajectory (the production-MS3 branch folds via the returning
      // scan, and the cascade branch has no MS3 to fold), so the winning variant is folded here, inline, at
      // group completion. Known asymmetry, not generalized: folding is per-group, so a precursor whose
      // fragments are explored in several groups folds each winner as its group closes rather than once at
      // the end. That only matters if a later group should be able to revise an earlier fold.
      else
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
      const ScanCommand* ms_ctx, ProteoformTracker* tracker, int precursor_id,
      const FragmentAnalysis::SearchedPrSM* searched)
  {
    return applyNextLevel(searchNextLevel(msn_level, result, ms_ctx, searched), msn_level, result, faims_cv, queue, ms_ctx,
                          tracker, precursor_id);
  }

  // The measuring half (ADR-0062): reads only its arguments and the immutable Config, so a search worker
  // runs it with no lock held. The FragmentAnalysis it drives holds nothing but the Config.
  Exploration::NextLevelSearch Exploration::searchNextLevel(int msn_level, const DeconvolvedSpectrum& result,
                                                            const ScanCommand* ms_ctx,
                                                            const FragmentAnalysis::SearchedPrSM* searched) const
  {
    NextLevelSearch nlr;

    int next_level = msn_level + 1;
    const auto& this_cfg = config_.level(msn_level);
    const auto& next_cfg = config_.level(next_level);
    // The protein matched against. A searched MS2 (ADR-0061 push 4c): the one the proteome search chose, empty when it
    // chose none. Otherwise the single-protein matcher's fallback, the inline reference (ADR-0061 push 4b).
    const std::string& seq = (searched != nullptr) ? searched->protein_sequence
                                                   : config_.identification().proteinOr(kNoChosenProtein);

    // Identification is gated by the REFERENCE, not by characterization.mode. The four conditions that
    // used to bail here gate MS3 DISPATCH, and they are re-asserted below, immediately above the
    // next_cfg.scans[0] read they protect. Bailing on them HERE cost every "mode": "off" run its whole
    // identification payload -- no sequence tags, no fragment matches, no identification.tsv row at all --
    // because this function owns identification as well as MS3 target selection, and the former sits above
    // the latter. With no reference there is nothing to match against, so that is the honest gate. Not the
    // configured sequence (ADR-0061 push 4c): a FASTA reference has none, and an MS2 with no PrSM (or matched
    // against no protein) is an unidentified MS2 whose model may still dispatch, not a skipped one.
    if (!config_.identification().enabled()) return nlr;
    nlr.ran = true;

    int num_targets = this_cfg.max_targets;

    // Extract activation type from the scan command that produced this result.
    // The fallback indexes this_cfg.scans[0], which the entry gate above used to guard. On both reachable
    // paths it is still non-empty -- validate() requires ms_settings.ms2 whenever level 1 selects, and an
    // MS3 group can only exist once dispatch has run, which requires ms_settings.ms3 -- but that is a
    // three-step reachability argument standing where an explicit guard used to be, so keep it explicit.
    // Same value on every path that reaches it; this only refuses to read past the end of an empty list.
    std::string scan_activation = (ms_ctx != nullptr && ms_ctx->num_stages > 0)
        ? std::string(ms_ctx->stages[0].activation_type)
        : (this_cfg.scans.empty() ? std::string() : this_cfg.scans[0].activation);

    // Get fragment ion targets. The SELECTION already lives on ProteoformTracker
    // (selectNextLevelTargets, below); what remains here is the out-parameter marshalling it inherited
    // from the C-style signature. Seven parallel arrays go in and only two come back out meaningfully --
    // frag_result and qscores (summed into tic_coverage). masses/charges/wstarts/wends/ion_types/
    // frag_indices/frag_scores are filled and never read again, because the real target descriptors
    // arrive later via tracker->planNextScans. Collapsing them is safe but touches the shared signature.
    DeconvolvedSpectrum result_copy = result;
    std::vector<double> masses(num_targets), qscores(num_targets);
    std::vector<double> wstarts(num_targets), wends(num_targets);
    std::vector<int> charges(num_targets);
    std::vector<char> ion_types(num_targets, '\0');
    std::vector<int> frag_indices(num_targets, 0);
    std::vector<FragmentAnalysis::FragmentScores> frag_scores(num_targets);  // stage-1 (fragment) scoring for MS3 cmd rows

    int found = 0;
    FragmentAnalysis::ProteoformMatch& frag_result = nlr.frag_result;

    // #46: the metric->matcher selection switch now lives in ProteoformTracker::selectNextLevelTargets
    // (byte-identical: same matcher, same args). Behavior fixes (bail-early / num_targets bound) are deferred.
    // selectNextLevelTargets chooses its matcher by switching on the NEXT level's metric, and its default
    // arm returns 0 without calling one. Config projects characterization.mode onto that metric, so under
    // "mode": "off" it is None -- the switch falls to default and matches nothing. Relaxing the entry gate
    // above therefore buys nothing on its own; the matcher has to be named explicitly. getTopFragmentMatches
    // is the only matcher reachable from method.json since ADR-0014 deleted the chooser, and it is what
    // every MS3-enabled config already resolves to, so this is a no-op wherever the metric is set.
    const SelectionMetric matcher =
        (next_cfg.selection == SelectionMetric::None) ? SelectionMetric::QScore : next_cfg.selection;

    found = ProteoformTracker::selectNextLevelTargets(fragments_, matcher, seq, num_targets,
        masses.data(), qscores.data(), charges.data(),
        wstarts.data(), wends.data(),
        ion_types.data(), frag_indices.data(), result_copy, frag_result, scan_activation, 0.0, frag_scores.data(), searched);
    // The chosen protein travels with the match (ADR-0061 push 4b): a match that identified something carries
    // the accession of the protein it was matched against. Overwrites Refine's "input_protein" placeholder. A
    // searched match already carries its PrSM's accession.
    if (searched == nullptr && !frag_result.proteoform_sequence.empty())
      frag_result.matched_protein = config_.identification().accessionOr(kNoChosenProtein);

    // RENDER context (build time): the TRIGGERING scan's frag_result, NOT the (unfinalized) winner. At
    // MS2/next-level build time the winner is not yet finalized, so buildWinnerProteoformContext would
    // return an empty context -> fragmentProForma renders "" (empty/collapsed ms3_proteoform). This one
    // source flows to group.proteoform_ctx (:215), the three buildMS3 ms3_proteoform renders and the
    // cached MS2Context. Winner SCORING is decoupled: the returning MS3 (FLASHIda.cpp) and the completed
    // exploration group (below) score against a fresh buildWinnerProteoformContext instead.
    MS3FragmentMatcher::ProteoformContext& proto_ctx = nlr.proto_ctx;
    if (next_level >= 3)
    {
      proto_ctx.region_start = frag_result.region_start;
      proto_ctx.region_end = frag_result.region_end;
      proto_ctx.ptm_sites = frag_result.ptm_sites;
      if (proto_ctx.region_start < 0) proto_ctx.region_start = 0;
      if (proto_ctx.region_end < 0) proto_ctx.region_end = static_cast<int>(config_.identification().proteinOr(seq).size());
      // The protein the match chose (empty when it identified nothing: readers then take the fallback).
      proto_ctx.protein_sequence = frag_result.proteoform_sequence;
      proto_ctx.protein_accession = frag_result.matched_protein;
    }

    // Populate fragment matching metadata. This is report assembly, not identification -- the match was
    // produced by ProteoformTracker::selectNextLevelTargets above; these lines only shape it for the
    // caller's log rows. Moving it onto the tracker would give it a caller-facing return type it does not
    // otherwise need, so it stays here by choice, not by oversight.
    nlr.fragment_count = frag_result.total_match_count;
    if (!seq.empty() && found > 0)
    {
      nlr.matched_protein = config_.identification().accessionOr(frag_result.matched_protein);
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

    return nlr;
  }

  // The acquiring half (ADR-0062): everything below the seam, on the analysis thread under its lock.
  Exploration::NextLevelResult Exploration::applyNextLevel(const NextLevelSearch& search, int msn_level,
      const DeconvolvedSpectrum& result, double faims_cv, ScanCommandQueue& queue,
      const ScanCommand* ms_ctx, ProteoformTracker* tracker, int precursor_id)
  {
    NextLevelResult nlr;
    if (!search.ran) return nlr;

    const int next_level = msn_level + 1;
    const auto& this_cfg = config_.level(msn_level);
    const auto& next_cfg = config_.level(next_level);
    const FragmentAnalysis::ProteoformMatch& frag_result = search.frag_result;
    const MS3FragmentMatcher::ProteoformContext& proto_ctx = search.proto_ctx;
    nlr.fragment_count = search.fragment_count;
    nlr.proteoform_match = frag_result;
    nlr.matched_protein = search.matched_protein;
    nlr.proteoform_sequence = search.proteoform_sequence;
    nlr.tic_coverage = search.tic_coverage;

    // MS3 DISPATCH GATE -- the four conditions lifted off the entry gate at the top of this function.
    //
    // It sits HERE, and not one line lower, because the very next statement indexes next_cfg.scans[0]
    // UNGUARDED. Config::validate requires an ms_settings.ms3 block only when characterization.mode != off,
    // and its own comment calls the reverse case "the direction that segfaults". That validation covered
    // this read only for as long as mode == off could never reach this function -- which is exactly what
    // the relaxed entry gate changes. A "mode": "off" config with no ms_settings.ms3 block now arrives
    // here, so the guard has to travel with the read rather than stay at the door.
    //
    // Everything ABOVE this line is measurement: fragment matching plus the report assembly that shapes it
    // for the caller's log rows. Everything BELOW acquires. That is the whole seam -- identification is
    // gated by the sequence, dispatch by the mode.
    if (this_cfg.selection == SelectionMetric::None || next_cfg.selection == SelectionMetric::None
        || this_cfg.scans.empty() || next_cfg.scans.empty())
      return nlr;

    // Build commands for each selected fragment target
    ScanConfig next_scan_config = next_cfg.scans[0];
    int charge_floor = this_cfg.min_charge;

    // Gate the pooled baseline (feedScan + finalizeMS2) on the actual match, NOT on found>0. found is the
    // MS3-target count; under AmbiguityResolution an unambiguous ID selects zero targets (found==0) yet still
    // carries a real proteoform match — keying on the match keeps that ID in pooled_identification.tsv. For
    // intensity/qscore metrics found>0 iff fragments non-empty, so this is equivalent for shipped configs.
    // Exhaustive builds a model WITHOUT a match too (ADR-0055): its pool is the winner scan's masses, and for an
    // unidentified Precursor the winner scan is the one this dispatch comes from -- the regular MS2, the MS2
    // exploration metric's winner (the cascade passes variants[best_idx]), or the production follow-up MS2.
    const bool identified = !frag_result.proteoform_sequence.empty() && !frag_result.fragments.empty();
    const bool exhaustive = config_.characterization().objective == CharacterizationObjective::Exhaustive;
    if (tracker != nullptr && next_level >= 3 && ms_ctx != nullptr && ms_ctx->num_stages > 0
        && (identified || exhaustive))
    {
      // NOT a bail: this block RUNS when no model exists. It is the one-shot baseline for the
      // non-exploration MS2->MS3 path, where nothing upstream has fed the tracker -- the MS2 exploration
      // path finalizes its own model in feedResult, so `existing` is already populated there and this is
      // skipped. Turning it into an early-out would leave the regular path with no model, and
      // planNextScans below returns zero targets without one: no MS3 at all, silently.
      // An unidentified model is re-pointed at every new dispatch scan; dispatched_nominal_masses survives
      // finalize, so a species already fragmented for this Precursor is never re-targeted.
      ensureModel_(precursor_id, result, frag_result, identified, *ms_ctx, *tracker);
    }

    // The acquiring tail -- planNextScans over the finalized model, the per-target MS3 builds and the
    // [MS3-DISPATCH] marker -- lives in dispatchFromModel so the refine MS3 return path can re-dispatch
    // through the same code after a decision. Only commands and their contexts come back from it.
    NextLevelResult dispatched = dispatchFromModel(precursor_id, result, faims_cv, queue, ms_ctx, tracker, proto_ctx,
                                                   next_scan_config, charge_floor, next_level, nlr.proteoform_sequence);
    nlr.commands = std::move(dispatched.commands);
    nlr.ms3_contexts = std::move(dispatched.ms3_contexts);
    nlr.governed = std::move(dispatched.governed);   // ADR-0066: engaged only when per_proteoform governs the dispatch
    return nlr;
  }

  void Exploration::ensureModel_(int precursor_id, const DeconvolvedSpectrum& result,
                                 const FragmentAnalysis::ProteoformMatch& frag_result, bool identified,
                                 const ScanCommand& ms_ctx, ProteoformTracker& tracker)
  {
    const ProteoformModel* existing = tracker.getModel(precursor_id);
    if (existing == nullptr || existing->proteoform_sequence.empty())
    {
      Ms2Params p;
      p.collision_energy = ms_ctx.stages[0].collision_energy;
      p.activation_type  = std::string(ms_ctx.stages[0].activation_type);
      p.reaction_time    = ms_ctx.stages[0].reaction_time;
      tracker.feedScan(precursor_id, /*ms_level=*/2, p, ms_ctx.scan_id, result, frag_result,
                       frag_result.score, ms_ctx);
      tracker.finalizeMS2(precursor_id, identified ? std::nullopt : std::optional<int>(ms_ctx.scan_id));
    }
  }

  Exploration::QuantDispatch Exploration::dispatchQuant(int precursor_id, const FragmentAnalysis::ProteoformMatch& frag_result,
                                                        const DeconvolvedSpectrum& result, double faims_cv, ScanCommandQueue& queue,
                                                        const ScanCommand& ms_ctx, ProteoformTracker& tracker,
                                                        const ScanConfig& quant_scan, int n)
  {
    QuantDispatch out;
    const bool identified = !frag_result.proteoform_sequence.empty() && !frag_result.fragments.empty();
    if (!identified || ms_ctx.num_stages <= 0 || n <= 0) return out;
    ensureModel_(precursor_id, result, frag_result, identified, ms_ctx, tracker);

    bool has_cysteine = false;
    const std::vector<Ms3Target> candidates = tracker.labelledFragments(precursor_id, 'C', &has_cysteine);
    if (!has_cysteine)
    {
      out.outcome = "no_cysteine";
      EventLine() << "[QUANT-MS3] precursor_id=" << precursor_id << " targets=0 reason=no_cysteine";
      return out;
    }
    const MS3FragmentMatcher::ProteoformContext proto_ctx = tracker.buildWinnerProteoformContext(precursor_id);
    const MSSpectrum& source = result.getOriginalSpectrum();
    int eligible = 0;
    // Window SNR > 1 (ADR-0064 d.4, fixed): the fragment's own signal outweighs the rest of the window it would be
    // isolated in, read on this MS2 with buildMS3's floored width -- the formula the characterization MS3 logs.
    auto windowSnrOf = [&source](const Ms3Target& t) {
      const double width = std::max(t.iso_width, 2.0);
      return FragmentAnalysis::windowSnr(source, t.frag_mz - width / 2.0, t.frag_mz + width / 2.0,
                                         t.stage1_scores.precursor_intensity);
    };
    const std::optional<YieldGateConfig>& yield_gate = config_.quantification().ms3_yield_gate;
    bool yield_refused = false;   // ADR-0069: a fragment had an SNR-positive charge and the yield gate refused every one
    // ADR-0069: the fragment's charges strongest first, each through the SNR test and then the yield test; the first to
    // pass both buys the MS3. None: the fragment is refused. @p any_snr: it had an SNR-positive charge.
    auto walkCharges = [&](const Ms3Target& candidate, Ms3Target& t, double& snr, double& e, bool& any_snr) {
      for (const Ms3Target& c : tracker.labelledFragmentCharges(precursor_id, candidate.ion_type, candidate.ion_index))
      {
        const double s = windowSnrOf(c);
        if (!(s > 1.0)) continue;
        any_snr = true;
        const YieldInputs in = yieldInputsOf(ms_ctx, c.stage1_scores, c.ion_index);
        const YieldPrediction pred = predictBondYield(yield_gate->coefficients, in);
        if (yieldAdmits(*yield_gate, pred))
        {
          t = c;
          snr = s;
          e = pred.e;
          return true;
        }
        logYieldVerdict("quantification", precursor_id, c.ion_type, c.ion_index, std::abs(c.frag_charge), "refused",
                        *yield_gate, in, pred);
      }
      return false;
    };
    auto buy = [&](const Ms3Target& t, double snr) {
      const char ion_type = t.ion_type.empty() ? '\0' : t.ion_type[0];
      // Stage 0 is the identification scan's own isolation and activation (ADR-0003): no per-ion override.
      ScanCommand cmd = queue.buildMS3(ms_ctx, quant_scan, t.frag_mz, std::abs(t.frag_charge), t.iso_width, ms_ctx.scan_id,
                                       ion_type, t.ion_index, 1, t.stage1_scores, nullptr, proto_ctx, nullptr, 'Q');
      cmd.faims_cv = faims_cv;
      cmd.window_snr = snr;
      out.commands.push_back(cmd);
    };
    if (yield_gate && yield_gate->rank_by == YieldRank::Yield)
    {
      // ADR-0069 (amended): every labelled fragment is walked first, then the buyable ones are taken by predicted
      // bond yield, highest first (stable: equal yields keep the strongest-first order), up to n.
      struct Buyable
      {
        Ms3Target target;
        double snr = 0.0;
        double e = 0.0;
      };
      std::vector<Buyable> buyable;
      for (const Ms3Target& candidate : candidates)
      {
        Buyable b;
        b.target = candidate;
        bool any_snr = false;
        const bool bought = walkCharges(candidate, b.target, b.snr, b.e, any_snr);
        if (any_snr) ++eligible;   // eligible counts fragments with an SNR-positive charge, as without the gate
        if (bought) buyable.push_back(b);
        else yield_refused = yield_refused || any_snr;
      }
      std::stable_sort(buyable.begin(), buyable.end(), [](const Buyable& a, const Buyable& b) { return a.e > b.e; });
      for (const Buyable& b : buyable)
      {
        if (static_cast<int>(out.commands.size()) >= n) break;
        buy(b.target, b.snr);
      }
    }
    else
    {
      for (const Ms3Target& candidate : candidates)
      {
        if (static_cast<int>(out.commands.size()) >= n) break;
        Ms3Target t = candidate;
        double snr = 0.0;
        if (!yield_gate)
        {
          snr = windowSnrOf(t);
          if (!(snr > 1.0)) continue;
          ++eligible;
        }
        else
        {
          // rank_by intensity: strongest fragment first; a refused fragment hands its place to the next.
          bool any_snr = false;
          double e = 0.0;
          const bool bought = walkCharges(candidate, t, snr, e, any_snr);
          if (any_snr) ++eligible;   // eligible counts fragments with an SNR-positive charge, as without the gate
          if (!bought)
          {
            yield_refused = yield_refused || any_snr;
            continue;
          }
        }
        buy(t, snr);
      }
    }
    // Nothing bought: yield_refused when the quantification yield gate refused an eligible fragment (ADR-0069), else no
    // labelled fragment was eligible.
    out.outcome = !out.commands.empty() ? "bought" : (yield_refused ? "yield_refused" : "no_labelled_fragment");
    {
      EventLine ev;
      ev << "[QUANT-MS3] precursor_id=" << precursor_id << " targets=" << out.commands.size()
         << " labelled=" << candidates.size() << " eligible=" << eligible << " budget=" << n;
      if (out.commands.empty()) ev << " reason=" << out.outcome;
    }
    return out;
  }

  Exploration::NextLevelResult Exploration::dispatchFromModel(int precursor_id, const DeconvolvedSpectrum& result,
      double faims_cv, ScanCommandQueue& queue, const ScanCommand* ms_ctx, ProteoformTracker* tracker,
      MS3FragmentMatcher::ProteoformContext proto_ctx, const ScanConfig& next_scan_config, int charge_floor,
      int next_level, const std::string& proteoform_sequence)
  {
    NextLevelResult nlr;
    nlr.proteoform_sequence = proteoform_sequence;

    const ProteoformModel* model = (tracker != nullptr && next_level >= 3 && ms_ctx != nullptr)
        ? tracker->getModel(precursor_id) : nullptr;
    // Exhaustive dispatches for an unidentified model too (ADR-0055); planExhaustive_ refuses one with no pool.
    if (model != nullptr && (!model->proteoform_sequence.empty()
                             || config_.characterization().objective == CharacterizationObjective::Exhaustive))
    {
      // Render/annotation seed = the finalized ProteoformTracker model (the SAME object
      // pooled_identification emits), NOT the exploration-metric winner's frag_result. Under a CE sweep
      // the metric winner (e.g. mass_count) can carry a fused blind-mod decomposition while the id-best
      // winner is split; seeding the render from the model keeps scan_commands ms3_proteoform and the
      // identification.tsv MS3 leaf == pooled. buildWinnerProteoformContext returns a region_start < 0
      // context when there is no finalized winner, so this is a no-op there (and where they coincide).
      MS3FragmentMatcher::ProteoformContext win_ctx = tracker->buildWinnerProteoformContext(precursor_id);
      if (win_ctx.region_start >= 0) proto_ctx = win_ctx;

      // ADR-0066: characterization.per_proteoform governs the dispatch of an IDENTIFIED model (an unidentified exhaustive
      // model has no proteoform key and stays ungoverned). The plan is made against the proteoform's registry snapshot
      // through the tracker's gate, and comes back as items for the caller's replaceMs3Plan instead of commands.
      // Ungoverned, everything below runs exactly as before.
      std::optional<GovernedMs3> gov;
      Ms3RegistrySnapshot snap;
      Ms3PlanGate gate;
      if (config_.characterization().per_proteoform && !model->proteoform_sequence.empty())
      {
        gov.emplace();
        const int P = static_cast<int>(model->proteoform_sequence.size());
        gov->key.accession = model->matched_protein;
        gov->key.region_start = (model->region_start < 0) ? 0 : model->region_start;   // normalised as planNextScans does
        gov->key.region_end = (model->region_end < 0) ? P : model->region_end;
        gov->key.commanded_mass = ms_ctx->mono_mass;                                    // the commanded precursor mass
        gov->ms2_scan_id = ms_ctx->scan_id;
        snap = queue.ms3Snapshot(gov->key);
        gov->snapshot_generation = snap.generation;
        gate.snapshot = &snap;
        gate.spectrum = &result.getOriginalSpectrum();   // the window SNR below is computed over the same spectrum
        gate.charge_floor = charge_floor;
        gate.max_targets = config_.characterization().max_targets;
      }
      // ADR-0069: the yield gate, governed or not. Ungoverned, the gate has no snapshot and no spectrum, so the plan gains
      // only the isolatable / charge-floor tests in the charge walk (the skips below then never fire) and the yield test.
      const bool yield_gated = config_.characterization().yield_gate.has_value();
      if (yield_gated)
      {
        gate.yield = &*config_.characterization().yield_gate;
        gate.ms_ctx = ms_ctx;   // stage 0 of every MS3 buildMS3 makes below
        if (!gov)
        {
          gate.charge_floor = charge_floor;
          gate.max_targets = config_.characterization().max_targets;
        }
      }

      auto targets = tracker->planNextScans(precursor_id, (gov || yield_gated) ? &gate : nullptr);
      // characterization.order: under "random" the plan's MS3 are sent in one reproducible permutation, drawn from the
      // dispatching scan's id. Only the send order changes: the plan above (its targets, slots and stamps) is made.
      if (config_.characterization().send_order == Ms3SendOrder::Random && targets.size() > 1)
      {
        const uint64_t seed = kSweepSeed ^ static_cast<uint64_t>(ms_ctx != nullptr ? ms_ctx->scan_id : precursor_id);
        const std::vector<int> perm = sweepPermutation(static_cast<int>(targets.size()), seed);
        std::vector<Ms3Target> shuffled;
        shuffled.reserve(targets.size());
        for (const int p : perm) shuffled.push_back(targets[static_cast<size_t>(p)]);
        targets = std::move(shuffled);
      }

      // A raw-mass plan (exhaustive, or a discriminating round still undecided after planning) stamped each target's
      // nominal mass; its items carry it so a withdrawal can give it back. A mapped plan stamps nothing.
      const CharacterizationObjective objective = config_.characterization().objective;
      const bool raw_mass_plan = objective == CharacterizationObjective::Exhaustive
                                 || (objective == CharacterizationObjective::Discriminating && model->refine.has_value()
                                     && !model->refine->decided);
      // One governed item: what the target acquires (cleavage, the anchor and stage-1 notch charges at their MS2
      // intensities), who planned it, and the stamp it set. @p cmds: one command, or one sweep group's variants.
      auto makeItem = [&](const Ms3Target& target, int frag_charge, const std::vector<ScanCommand>& cmds) {
        Ms3PlanItem item;
        item.cleavage = ProteoformTracker::ms3CleavageOf(target);
        item.charges.emplace_back(frag_charge, cmds.front().precursor_intensity_s1);
        for (const NotchCandidate& n : target.notches) item.charges.emplace_back(std::abs(n.charge), n.intensity);
        item.commands = cmds;
        item.owner_precursor_id = precursor_id;
        item.fragment_mass = target.frag_mass;
        item.nominal = raw_mass_plan ? SpectralDeconvolution::getNominalMass(target.frag_mass) : 0;
        const auto v = variant_tracking_map_.find(cmds.front().scan_id);
        if (v != variant_tracking_map_.end()) item.sweep_group_id = v->second.group_id;
        return item;
      };

      // Already bounded: planNextScans stops adding once targets.size() == config_.level(2).max_targets,
      // which is the same value num_targets was read from above (this_cfg == level 2). The two skips below
      // (charge floor, unisolatable frag_mz) can only yield fewer targets, never more. NB with MS3
      // exploration on, each target still fans out into one command PER CE VARIANT -- the budget bounds
      // targets, not commands.
      //
      // Counters for the [MS3-DISPATCH] marker below (ADR-0023 D-e). Counted HERE rather than on the
      // tracker's [MS3-PLAN] line because buildVariants_ is private to this class and the tracker
      // holds no reference to it -- and because the two skips below mean the tracker's target count
      // is an upper bound, not what was dispatched.
      int dispatched_targets = 0;
      int variants_per_target = 0;
      for (const Ms3Target& target : targets)
      {
        const int frag_charge = std::abs(target.frag_charge);
        if (charge_floor > 0 && frag_charge < charge_floor) continue;
        // A tracker match with no backing peak carries frag_mz/frag_charge 0 -- ProteoformTracker
        // initialises matched_mz/matched_charge to 0 and deliberately keeps the observation when no
        // peak is found ("never drop a matched fragment for lack of a peak"). Such a fragment is a
        // real match but cannot be isolated, so it is not a dispatchable MS3 target. Skipping here
        // costs no tracking id and writes no log row; the charge_floor guard above does not cover it
        // because characterization.min_fragment_charge defaults to 0, leaving that test inert.
        if (target.frag_mz <= 0.0 || frag_charge <= 0) continue;
        const char ion_type = target.ion_type.empty() ? '\0' : target.ion_type[0];
        ++dispatched_targets;

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
                                   ion_type, target.ion_index, proto_ctx, target.stage1_scores, &target.stage0_params,
                                   &target.notches);
          if (gov)
          {
            if (!sub_cmds.empty()) gov->items.push_back(makeItem(target, frag_charge, sub_cmds));   // one group, one item
          }
          else
            nlr.commands.insert(nlr.commands.end(), sub_cmds.begin(), sub_cmds.end());
          variants_per_target = static_cast<int>(sub_cmds.size());

          for (size_t sci = 0; sci < sub_cmds.size(); ++sci)
          {
            MS2Context mc;
            mc.proteoform_sequence = nlr.proteoform_sequence;
            mc.start_pos = proto_ctx.region_start;
            mc.end_pos = proto_ctx.region_end;
            mc.ptm_sites = proto_ctx.ptm_sites;
            mc.protein_sequence = proto_ctx.protein_sequence;
            mc.protein_accession = proto_ctx.protein_accession;
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
            if (gov)
              gov->contexts.emplace_back(sub_cmds[sci].scan_id, mc);
            else
              nlr.ms3_contexts.push_back(mc);
          }
        }
        else
        {
          // Single MS3: build directly from the target descriptors, applying the model's stage[0].
          ScanCommand cmd = queue.buildMS3(*ms_ctx, next_scan_config, target.frag_mz, frag_charge,
                                           target.iso_width, ms_ctx->scan_id,
                                           ion_type, target.ion_index, 1, target.stage1_scores, &target.stage0_params,
                                           proto_ctx, &target.notches);
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
          EventLine() << "[TRACK-CREATE] id=" << id_str
                    << " ms_level=" << next_level << " type=next_level";

          if (gov)
            gov->items.push_back(makeItem(target, frag_charge, {cmd}));   // one command, one item (separate: per charge)
          else
            nlr.commands.push_back(cmd);
          variants_per_target = 1;  // no CE sweep at this level: one command per target

          MS2Context mc;
          mc.proteoform_sequence = nlr.proteoform_sequence;
          mc.start_pos = proto_ctx.region_start;
          mc.end_pos = proto_ctx.region_end;
          mc.ptm_sites = proto_ctx.ptm_sites;
          mc.protein_sequence = proto_ctx.protein_sequence;
          mc.protein_accession = proto_ctx.protein_accession;
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
          if (gov)
            gov->contexts.emplace_back(cmd.scan_id, mc);
          else
            nlr.ms3_contexts.push_back(mc);
        }
      }

      // The dispatch half of the plan marker (ADR-0023 D-e). commands is the figure that actually
      // predicts instrument cost: max_targets bounds TARGETS, and with an MS3 CE sweep on each
      // target fans out into one command per variant -- so a budget of 20 can be ~100 scans.
      // Emitted even at zero targets, so a grep for [MS3-DISPATCH] pairs with every [MS3-PLAN]
      // line. stdout only, like every engine marker.
      size_t dispatched_commands = nlr.commands.size();
      if (gov)
        for (const Ms3PlanItem& item : gov->items) dispatched_commands += item.commands.size();
      EventLine() << "[MS3-DISPATCH] precursor_id=" << precursor_id
                << " targets=" << dispatched_targets
                << " variants=" << variants_per_target
                << " commands=" << dispatched_commands;
      if (dispatched_targets > 0 && !nlr.commands.empty())
      {
        // ADR-0065: the operator's line for this MS3 batch, from the commands [MS3-DISPATCH] counts.
        std::ostringstream op;
        op << '#' << ScanCommandQueue::encode(ms_ctx->scan_id) << "  " << dispatched_targets << " targets ";
        std::vector<std::string> ions;
        for (const MS2Context& mc : nlr.ms3_contexts)
        {
          const std::string ion = std::string(1, mc.fragment_ion_type) + std::to_string(mc.fragment_ion_index);
          if (std::find(ions.begin(), ions.end(), ion) == ions.end()) ions.push_back(ion);
        }
        for (const std::string& ion : ions) op << ' ' << ion;
        const ScanCommand& first = nlr.commands.front();
        if (first.num_stages >= 2) op << "  " << first.stages[1].activation_type << ' ' << first.stages[1].collision_energy;
        op << "  (" << nlr.commands.size() << " scans)";
        EventLog::say(OperatorKind::Ms3, op.str());
      }

      nlr.governed = std::move(gov);
      return nlr;
    }

    return nlr;
  }

  int Exploration::activeGroupCount() const
  {
    return static_cast<int>(active_groups_.size());
  }

  int Exploration::activeGroupCount(int msn_level) const
  {
    int n = 0;
    for (const auto& kv : active_groups_)
      if (kv.second.msn_level == msn_level) ++n;
    return n;
  }

  bool Exploration::dropGroup(int group_id)
  {
    auto git = active_groups_.find(group_id);
    if (git == active_groups_.end()) return false;
    for (const auto& v : git->second.variants)
      if (v.received) return false;
    for (const auto& v : git->second.variants)
      variant_tracking_map_.erase(v.cmd.scan_id);
    active_groups_.erase(git);
    EventLine() << "[EXPL-DROP] group=" << group_id;
    return true;
  }

  int Exploration::lastGroupId() const
  {
    return next_group_id_ - 1;
  }

  Exploration::ExplorationGroup Exploration::getGroup(int group_id) const
  {
    return active_groups_.at(group_id);
  }

  double Exploration::computeExplorationScore_(ExplorationMetric metric,
      const DeconvolvedSpectrum& spec,
      const ExplorationGroup& group, double share,
      FragmentAnalysis::ProteoformMatch* out_frag,
      const std::string& activation_type,
      const FragmentAnalysis::ProteoformMatch* searched) const
  {
    // Every Orbitrap variant is matched under every metric at MS2 (the whole-protein matcher is the
    // right one there); at MS3 the FragmentCount batch re-score is the matcher that counts. A search
    // worker may already have run the same matcher on the same spectrum (ADR-0062).
    // Inline, an MS3 variant is matched against its group's CHOSEN protein (ADR-0061 push 4c; none chosen: nothing),
    // an MS2 variant no search ran on against the inline reference.
    const bool ms3 = group.msn_level >= 3;
    FragmentAnalysis::ProteoformMatch fmr = (searched != nullptr) ? *searched
        : computeFragmentMatch_(spec, group.msn_level, activation_type,
                                ms3 ? group.proteoform_ctx.protein_sequence : config_.identification().proteinOr(kNoChosenProtein),
                                config_.identification().accessionOr(ms3 ? group.proteoform_ctx.protein_accession : kNoChosenProtein));
    if (out_frag) *out_frag = fmr;
    switch (metric)
    {
      case ExplorationMetric::MassCount:
        return computeMassCount_(spec);
      case ExplorationMetric::RemainingPrecursor:
        return computeRemainingPrecursorScore_(group, share);
      case ExplorationMetric::FragmentCount:
        return static_cast<double>(fmr.total_match_count);
      case ExplorationMetric::KineticModel:
        // ADR-0047: no pre-scan is scored; the group is decided by the fit at completion. -1.0 is "not
        // scored", so the winner loop never crowns one.
        return -1.0;
      // ADR-0074: at MS2 the values the search's annotator put on the variant's own PrSM (-1.0 without one, so a
      // variant that identified nothing cannot win); at MS3 -1.0 here, scored by the annotator when the group completes.
      case ExplorationMetric::RawScore:
        return ms3 ? -1.0 : fmr.raw_score;
      case ExplorationMetric::BondCoverage:
        return ms3 ? -1.0 : fmr.bond_coverage;
      case ExplorationMetric::ExplainedSignal:
        return ms3 ? -1.0 : fmr.explained_signal;
      default:
        return computeMassCount_(spec);
    }
  }

  double Exploration::computeMassCount_(const DeconvolvedSpectrum& spec) const
  {
    return static_cast<double>(spec.size());
  }

  double Exploration::precursorShare_(const ExplorationGroup& group,
      const double* mzs, const double* ints, int length) const
  {
    const double iso_half = group.isolation_width / 2.0;
    return KineticModel::precursorShare(mzs, ints, length, group.precursor_mz - iso_half, group.precursor_mz + iso_half);
  }

  double Exploration::computeRemainingPrecursorScore_(const ExplorationGroup& group, double share) const
  {
    // No measurement (no arrays, or a scan with no signal at all): -1.0, NOT 0.0, and the distinction is
    // the whole bar -- winner selection seeds best_score at -1.0 and takes `score > best_score`, so a 0.0
    // here would WIN an unmeasured variant at zero, and at MS3 a measuring metric would then fire a
    // production scan at that coin-flip CE. A share of exactly 0 (signal elsewhere, none in the window)
    // IS a measurement and scores like any other.
    if (share < 0.0) return -1.0;
    const double target = config_.level(group.msn_level).remaining_precursor_target;
    double score = 1.0 - std::abs(share - target);
    if (score < 0.0) score = 0.0;  // floor: score is always in [0, 1]
    return score;
  }

  FragmentAnalysis::ProteoformMatch Exploration::computeFragmentMatch_(const DeconvolvedSpectrum& spec, int msn_level,
                                                                      const std::string& activation_type,
                                                                      const std::string& protein_sequence,
                                                                      const std::string& protein_accession,
                                                                      const FragmentAnalysis::SearchedPrSM* searched) const
  {
    FragmentAnalysis::ProteoformMatch result;
    // A searched MS2 is matched on the protein its PrSM lies in (ADR-0061 push 4c); the single-protein matcher on the
    // protein its caller names, and on nothing when that is empty.
    const std::string& seq = (searched != nullptr) ? searched->protein_sequence : protein_sequence;
    if (spec.empty() || (searched == nullptr && seq.empty()))
      return result;

    DeconvolvedSpectrum spec_copy = spec;

    const int max_matches = 100;
    std::vector<double> masses(max_matches), qscores(max_matches);
    std::vector<double> wstarts(max_matches), wends(max_matches);
    std::vector<int> charges(max_matches);
    std::vector<char> ion_types(max_matches, '\0');
    std::vector<int> frag_indices(max_matches, 0);

    ProteoformTracker::identifyExplorationFragments(
        fragments_, seq, max_matches,
        masses.data(), qscores.data(), charges.data(),
        wstarts.data(), wends.data(),
        ion_types.data(), frag_indices.data(),
        spec_copy, result, activation_type,
        config_.level(msn_level).exploration_tolerance_ppm, nullptr, searched);

    // The accession of the protein matched against (ADR-0061 push 4b), on every match that identified something.
    // Overwrites Refine's "input_protein" placeholder. A searched match already carries its PrSM's accession.
    if (searched == nullptr && (result.total_match_count > 0 || !result.proteoform_sequence.empty()))
    {
      result.matched_protein = protein_accession;
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
