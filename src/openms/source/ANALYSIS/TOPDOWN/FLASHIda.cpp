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
// $Maintainer: Tom David Müller$
// $Authors: Kyowon Jeong, Tom David Müller $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <functional>
#include <map>
#include <sstream>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{

/// constructor
FLASHIda::FLASHIda(char* arg) :
    config_(std::string(arg)),
    logger_(config_),
    queue_(config_),
    deconv_(config_, config_.toleranceList()),
    fragments_(config_),
    selection_(config_, deconv_),
    quant_(config_),
    faims_(config_),
    exploration_(config_, fragments_),
    tracker_(config_, logger_)
{
  #ifdef _OPENMP
    omp_set_num_threads(4);
  #endif

    engine_start_time_ = std::chrono::steady_clock::now();

    // Initialize FAIMS CV atomic for getNextScanCommand reads.
    // relaxed is safe: no other thread can observe this object before construction completes.
    current_faims_cv_.store(faims_.isEnabled() ? faims_.currentCV() : 0.0, std::memory_order_relaxed);
  }

  FLASHIda::~FLASHIda() = default;

  int FLASHIda::processScan(const double* mzs, const double* ints, int length,
                             double rt_min, int ms_level, const char* scan_description,
                             double faims_cv)
  {
    std::lock_guard<std::mutex> lock(analysis_mutex_);

    // Tracking ID extraction
    std::string desc_str = scan_description ? std::string(scan_description) : "";
    if (desc_str.size() < 3) return 0;
    std::string id_str = desc_str.substr(0, 3);
    int tracking_id = queue_.decode(id_str);

    // Skip AGC scans — calibration-only, no data to process
    if (desc_str.size() >= 4 && desc_str[3] == 'A')
    {
      queue_.resolvePending(tracking_id);
      return 0;
    }

    // Retrieve timestamps from pending map
    uint64_t enqueue_ts = 0;
    uint64_t dequeue_ts = 0;
    {
      auto peeked = queue_.peekPending(tracking_id);
      if (peeked.has_value())
      {
        enqueue_ts = peeked->enqueue_timestamp_ms;
        dequeue_ts = peeked->dequeue_timestamp_ms;
      }
      else
      {
        std::cout << "[TRACK-RESOLVE] id=" << id_str << " status=not_found" << std::endl;
        return 0;
      }
    }

    // Stamp received timestamp (instrument → C++ handoff)
    uint64_t received_ts = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()
      ).count()
    );

    auto resolved = queue_.resolvePending(tracking_id);
    if (!resolved.has_value())
    {
      // peekPending above normally guarantees this; reaching here means an upstream race.
      std::cout << "[TRACK-RESOLVE] id=" << id_str << " status=context_lost_race" << std::endl;
      return 0;
    }
    ScanCommand ctx = resolved.value();

    // Context-support gate, validated once before branching. required_stages = the stage count this
    // level needs in its queued command; ms_level>3 is unsupported (returns 0).
    const int required_stages = (ms_level == 1) ? 0 : (ms_level == 2) ? 1 : (ms_level == 3) ? 2 : -1;
    if (required_stages < 0 || ctx.msn_level != ms_level || ctx.num_stages < required_stages)
    {
      std::cout << "[TRACK-RESOLVE] id=" << id_str << " status=context_unsupported"
                << " ms_level=" << ms_level << " ctx_msn=" << ctx.msn_level
                << " num_stages=" << ctx.num_stages << std::endl;
      return 0;
    }

    // Locals shared by all levels: the follow-up-scan flag and the child-id accumulator.
    const bool is_follow_up_scan =
        std::strlen(ctx.scan_description) >= 4 &&
        (ctx.scan_description[3] == 'F' || ctx.scan_description[3] == 'C');
    std::vector<std::string> child_ids;

    // Carriers for the file-backed TSV rows (scan_results + identification): each branch fills these,
    // and they are written once after the branch dispatch below. The IDA log, the [TRACK-*] stdout
    // lines, and all side effects (queue pushes, atomics, cache erase) stay inline in each branch.
    int return_code = 0;
    bool has_scan_row = false;
    IdaLogger::ScanRowDescriptor scan_row{};
    std::vector<IdaLogger::IdRowDescriptor> id_rows;

    // scan_results fields that are identical in every branch — set once here. Each branch fills its
    // own branch-specific fields; the rest keep their scan_row{} sentinel defaults.
    scan_row.tracking_id = id_str;
    scan_row.ms_level    = ms_level;   // gate guarantees ms_level in {1,2,3}
    scan_row.rt          = rt_min;
    scan_row.enqueue_ts  = enqueue_ts;
    scan_row.dequeue_ts  = dequeue_ts;
    scan_row.received_ts = received_ts;

    if (ms_level == 1)
    {
      // ===== MS1: deconvolve -> select top-N -> push MS2 (or CE sweep) -> FAIMS cycle =====
      // Selection=none: skip MS1 precursor selection entirely
      if (config_.level(1).selection == SelectionMetric::None)
        return 0;

      int n = selection_.filterAndRank(mzs, ints, length, rt_min, 1, faims_cv);
      const auto& selected = selection_.selectedPeakGroups();
      const auto& sel_charges = selection_.triggerCharges();
      int commands_pushed = 0;
      std::vector<ScanCommand> ms2_commands;
      // Item 2: the survey MS1 spectrum every MS2 isolation-window SNR is measured against. Hoisted so both
      // the exploration path (passed into initiate) and the normal path (inline) can stamp window_snr.
      const MSSpectrum& raw_ms1 = selection_.deconvolvedMS1().getOriginalSpectrum();

      if (config_.hasExploration(2))
      {
        // Exploration path: initiate CE sweep variants INSTEAD of regular MS2
        for (int i = 0; i < n; i++)
        {
          ScanCommand ms1_ctx{};
          ms1_ctx.scan_id  = tracking_id;
          ms1_ctx.faims_cv = ctx.faims_cv;  // Item 1: CV travels via the context (the resolved MS1 carries it)
          // initiate stamps each command's parent_scan_id from ms1_ctx.scan_id (= tracking_id), reads the CV
          // off the context, and measures each variant's window_snr over raw_ms1.
          auto cmds = exploration_.initiate(2, selected[i], sel_charges[i], queue_, &raw_ms1, &ms1_ctx);
          for (auto& c : cmds)
          {
            queue_.push(c);
            ms2_commands.push_back(c);
            commands_pushed++;
          }
        }
      }
      else
      {
        // Normal path: push MS2 for each precursor, for each scan config
        for (int i = 0; i < n; i++)
        {
          for (const auto& sc : config_.level(2).scans)
          {
            ScanConfig ms2_config = sc;
            ScanCommand cmd = queue_.buildMS2(selected[i], sel_charges[i], ms2_config, 2, tracking_id);
            cmd.faims_cv = faims_cv;  // MS2 carries parent MS1's CV
            // Item 2: stamp the isolation-window SNR before pushing (was the standalone post-push loop +
            // queue window-SNR map). Same formula/inputs as before.
            if (cmd.num_stages > 0)
            {
              const double half = cmd.stages[0].isolation_width / 2.0;
              cmd.window_snr = FragmentAnalysis::windowSnr(raw_ms1, cmd.stages[0].precursor_mz - half,
                                                           cmd.stages[0].precursor_mz + half, cmd.precursor_intensity);
            }
            queue_.push(cmd);
            ms2_commands.push_back(cmd);
            commands_pushed++;
          }
        }
      }

      // IDA log entry (MS1 only).
      logger_.writeIDALogEntry(rt_min, tracking_id, id_str, ms2_commands, selection_.deconvolvedMS1());

      for (const auto& c : ms2_commands)
        child_ids.push_back(ScanCommandQueue::encode(c.scan_id));
      int ms1_mass_count = static_cast<int>(selection_.deconvolvedMS1().size());
      scan_row.mass_count      = ms1_mass_count;
      scan_row.commands_pushed = commands_pushed;
      scan_row.child_ids       = child_ids;
      scan_row.tag_count       = 0;   // MS1 logs 0 (overrides the -1 default)
      scan_row.deconv_spectrum = &selection_.deconvolvedMS1();
      scan_row.fragment_count  = 0;   // MS1 logs 0 (overrides the -1 default)
      has_scan_row             = true;

      // FAIMS CV cycling: update skip policy, advance to next CV, push MS1
      if (faims_.isEnabled())
      {
        double current_cv = faims_.currentCV();
        faims_.updateSkip(current_cv, commands_pushed);

        double next_cv = faims_.advanceToNextCV();
        ScanCommand ms1 = queue_.makeMS1();
        ms1.faims_cv = next_cv;
        ms1.scan_id = queue_.nextTrackingId();
        ms1.priority = 0;  // priority 0 to send before pending MS2s

        std::string next_ms1_id = ScanCommandQueue::encode(ms1.scan_id);
        std::snprintf(ms1.scan_description, 16, "%sS", next_ms1_id.c_str());

        queue_.push(ms1);
        std::cout << "[TRACK-CREATE] id=" << next_ms1_id << " ms_level=1 type=cv_transition cv=" << next_cv << std::endl;
      }

      // Update atomics for lock-free reads by getNextScanCommand
      exploration_active_.store(exploration_.activeGroupCount() > 0, std::memory_order_release);
      current_faims_cv_.store(faims_.isEnabled() ? faims_.currentCV() : 0.0, std::memory_order_release);

      return_code = commands_pushed;
    }
    else if (ms_level == 2)
    {
      int commands_pushed = 0;

      // ===== MS2 (exploration variant): score variant -> push children -> fill row =====
      if (exploration_.isExplorationVariant(tracking_id))
      {
        auto expl_result = exploration_.feedResult(tracking_id, mzs, ints, length, rt_min, queue_, tracker_);
        for (auto& c : expl_result.commands)
          queue_.push(c);  // parent_scan_id already stamped by feedResult; children pre-encoded in expl_result.child_ids

        int expl_mass_count = exploration_.explorationDeconvMassCount();
        const DeconvolvedSpectrum* expl_spec = exploration_.explorationDeconvSpectrum();
        std::string parent_id(expl_result.parent_scan_id);

        // Fill scan_results row in canonical order (copy expl_result.* by value — it dies at block end).
        // DIVERGENCE vs MS3-expl: MS2 sets tag_count + fragment_count; MS3-expl leaves both at -1.
        scan_row.mass_count          = expl_mass_count;
        scan_row.commands_pushed     = static_cast<int>(expl_result.commands.size());
        scan_row.child_ids           = expl_result.child_ids;
        scan_row.tag_count           = expl_result.identification_result.tag_count;
        scan_row.matched_protein     = expl_result.matched_protein;
        scan_row.proteoform_sequence = FragmentAnalysis::toProForma(expl_result.proteoform_sequence, expl_result.identification_result.ptm_sites);
        scan_row.deconv_spectrum     = expl_spec;
        scan_row.parent_tracking_id  = parent_id;
        scan_row.tic_coverage        = expl_result.tic_coverage;
        scan_row.fragment_count      = expl_result.fragment_count;
        scan_row.exploration_group_id = expl_result.group_id;
        scan_row.exploration_metric  = expl_result.exploration_metric;
        scan_row.variant_index       = expl_result.variant_index;
        scan_row.total_variants      = expl_result.total_variants;
        scan_row.collision_energy    = std::to_string(expl_result.collision_energy);
        scan_row.exploration_score   = expl_result.score;
        scan_row.remaining_ratio     = expl_result.remaining_ratio;
        scan_row.activation_type     = expl_result.activation_type;
        scan_row.reaction_time       = std::to_string(expl_result.reaction_time);
        scan_row.winner_tracking_id  = expl_result.winner_tracking_id;  // "" except the group-completing row
        has_scan_row = true;

        // MS2 identification row (copy by value before expl_result goes out of scope).
        // DIVERGENCE vs MS3-expl: MS2 gate requires a non-empty proteoform_sequence AND fragments.
        if (!expl_result.proteoform_sequence.empty() && !expl_result.identification_result.fragments.empty())
          id_rows.push_back({id_str, 2, 'E', expl_result.ms2_context, expl_result.identification_result});

        exploration_active_.store(exploration_.activeGroupCount() > 0, std::memory_order_release);
        return_code = static_cast<int>(expl_result.commands.size());
      }
      else
      {
        // ===== MS2 (regular): deconvolve -> targeting + follow-ups -> MS3 targeting -> row =====
        double precursor_mass = ctx.mono_mass;
        int precursor_charge = ctx.stages[0].charge_state;
        deconv_.deconvolveMSn(mzs, ints, length, rt_min, precursor_mass, precursor_charge);

        // Tag-based targeting
        bool tags_found = false;
        int tags_count = 0;   // real sequence-tag count (logged as tag_count, and cached for MS3 rows)
        if (selection_.hasTargetProteinDatabase() && precursor_mass > 0)
        {
          std::string ms2_activation = std::string(ctx.stages[0].activation_type);
          tags_count = selection_.processMS2ForTagBasedTargeting(precursor_mass, ms2_activation);
          tags_found = tags_count > 0;
        }

        // Quantification follow-up (independent of tags)
        if (config_.quantification().enabled && !config_.quantification().follow_up_scan.activation.empty()
            && !is_follow_up_scan)
        {
          if (quant_.isDifferentiallyAbundant(mzs, ints, length, rt_min, 2, "ms2_quant",
                                              config_.quantification().reporter_mz_tol, config_.quantification().fold_change_threshold, false))
          {
            auto followup = queue_.buildFollowUp(ctx, config_.quantification().follow_up_scan, 'F');
            queue_.push(followup);
            child_ids.push_back(ScanCommandQueue::encode(followup.scan_id));
            commands_pushed++;
          }
        }

        // Conditional MS2 follow-up -- only when tags detected AND explicitly enabled
        if (config_.targeting().conditional_ms2_enabled && tags_found && !is_follow_up_scan)
        {
          auto cond = queue_.buildFollowUp(ctx, config_.targeting().tagging_follow_up_scan, 'C');
          queue_.push(cond);
          child_ids.push_back(ScanCommandQueue::encode(cond.scan_id));
          commands_pushed++;
        }

        // MS3 targeting via selection_strategy
        Exploration::NextLevelResult ms3_targeting;
        if (config_.level(2).selection != SelectionMetric::None)
        {
          // initiateNextLevel stamps each command's parent with the MS2 id (ctx.scan_id) at creation.
          ms3_targeting = exploration_.initiateNextLevel(2, deconv_.storedMS2(), ctx.faims_cv, queue_, &ctx);
          for (auto& c : ms3_targeting.commands)
          {
            queue_.push(c);
            child_ids.push_back(ScanCommandQueue::encode(c.scan_id));
            commands_pushed++;
          }

          // Cache MS2 context for each MS3 command (for non-exploration identification.tsv)
          for (size_t ci = 0; ci < ms3_targeting.commands.size(); ++ci)
          {
            if (ms3_targeting.commands[ci].msn_level >= 3 && ci < ms3_targeting.ms3_contexts.size())
            {
              // Carry the parent MS2's identification tag count (real when a proteoform matched) to its
              // MS3 rows — not the FASTA-DB-gated tags_count, which is 0 unless a tag-targeting DB is loaded.
              ms3_targeting.ms3_contexts[ci].tag_count = (!ms3_targeting.proteoform_match.fragments.empty())
                  ? ms3_targeting.proteoform_match.tag_count : tags_count;
              ms2_context_cache_[ms3_targeting.commands[ci].scan_id] = ms3_targeting.ms3_contexts[ci];
            }
          }
        }

        // Write MS2 identification row if proteoform was matched
        if (!ms3_targeting.proteoform_sequence.empty() && !ms3_targeting.proteoform_match.fragments.empty())
        {
          Exploration::MS2Context ms2_ctx;
          ms2_ctx.proteoform_sequence = ms3_targeting.proteoform_match.proteoform_sequence;
          ms2_ctx.start_pos = ms3_targeting.proteoform_match.region_start;
          ms2_ctx.end_pos = ms3_targeting.proteoform_match.region_end;
          ms2_ctx.ptm_sites = ms3_targeting.proteoform_match.ptm_sites;
          ms2_ctx.ms1_precursor_mass = ctx.mono_mass;
          ms2_ctx.ms1_precursor_mz = ctx.stages[0].precursor_mz;
          ms2_ctx.ms1_precursor_charge = ctx.stages[0].charge_state;
          // MS2 isolation-window reporting from the resolved command (window-SNR was recorded at MS1
          // time, keyed by this scan_id). The MS3 triplet stays 0 on MS2 rows.
          if (ctx.num_stages > 0)
          {
            ms2_ctx.ms2_isolation_width = ctx.stages[0].isolation_width;
            ms2_ctx.ms2_charge_intensity = ctx.precursor_intensity;
            ms2_ctx.ms2_window_snr = ctx.window_snr;  // Item 2: SNR travels on the command (was the queue map)
          }
          id_rows.push_back({id_str, 2, 'R', ms2_ctx, ms3_targeting.proteoform_match});  // copy by value (ms2_ctx is local)
        }

        int ms2_mass_count = deconv_.hasStoredMS2() ? static_cast<int>(deconv_.storedMS2().size()) : 0;
        // Log the identification tagger's real tag count when a proteoform matched; otherwise fall back
        // to tags_count (the FASTA tag-targeting case, and plain-DDA selection==None which is legitimately 0).
        int tag_count = (!ms3_targeting.proteoform_match.fragments.empty())
            ? ms3_targeting.proteoform_match.tag_count : tags_count;
        const DeconvolvedSpectrum* ms2_spec = deconv_.hasStoredMS2() ? &deconv_.storedMS2() : nullptr;
        std::string parent_id(ctx.parent_scan_id);

        // Log the actual MS2 scan-config CE/activation/reaction (single stage) the engine used.
        std::string ms2_collision_energy = ctx.num_stages > 0 ? std::to_string(ctx.stages[0].collision_energy) : "0";
        std::string ms2_activation_type = ctx.num_stages > 0
            ? std::string(ctx.stages[0].activation_type)
            : config_.level(2).scans[0].activation;
        std::string ms2_reaction_time = ctx.num_stages > 0 ? std::to_string(ctx.stages[0].reaction_time) : "0";
        // Fill scan_results row; exploration fields stay at their sentinel defaults.
        scan_row.mass_count         = ms2_mass_count;
        scan_row.commands_pushed    = commands_pushed;
        scan_row.child_ids          = child_ids;
        scan_row.tag_count          = tag_count;
        scan_row.matched_protein    = ms3_targeting.matched_protein;
        scan_row.proteoform_sequence = FragmentAnalysis::toProForma(ms3_targeting.proteoform_sequence, ms3_targeting.proteoform_match.ptm_sites);
        scan_row.deconv_spectrum    = ms2_spec;
        scan_row.parent_tracking_id = parent_id;
        scan_row.tic_coverage       = ms3_targeting.tic_coverage;
        scan_row.fragment_count     = ms3_targeting.fragment_count;
        scan_row.collision_energy   = ms2_collision_energy;
        scan_row.activation_type    = ms2_activation_type;
        scan_row.reaction_time      = ms2_reaction_time;
        has_scan_row = true;

        std::cout << "[TRACK-RESOLVE] id=" << id_str
                  << " rt=" << rt_min
                  << " commands=" << commands_pushed
                  << std::endl;

        // Update atomic for lock-free reads by getNextScanCommand
        exploration_active_.store(exploration_.activeGroupCount() > 0, std::memory_order_release);

        return_code = commands_pushed;
      }
    }
    else if (ms_level == 3)
    {
      // ===== MS3 (exploration variant): score variant -> push children -> fill row =====
      if (exploration_.isExplorationVariant(tracking_id))
      {
        auto expl_result = exploration_.feedResult(tracking_id, mzs, ints, length, rt_min, queue_, tracker_);
        for (auto& c : expl_result.commands)
          queue_.push(c);  // parent_scan_id already stamped by feedResult; children pre-encoded in expl_result.child_ids

        int expl_mass_count = exploration_.explorationDeconvMassCount();
        const DeconvolvedSpectrum* expl_spec = exploration_.explorationDeconvSpectrum();
        std::string parent_id(expl_result.parent_scan_id);

        // Fill scan_results row in canonical order (copy expl_result.* by value — it dies at block end).
        // DIVERGENCE vs MS2-expl: MS3-expl leaves tag_count + fragment_count at the -1 sentinel.
        scan_row.mass_count          = expl_mass_count;
        scan_row.commands_pushed     = static_cast<int>(expl_result.commands.size());
        scan_row.child_ids           = expl_result.child_ids;
        scan_row.matched_protein     = expl_result.matched_protein;
        scan_row.proteoform_sequence = FragmentAnalysis::toProForma(expl_result.proteoform_sequence, expl_result.identification_result.ptm_sites);
        scan_row.deconv_spectrum     = expl_spec;
        scan_row.parent_tracking_id  = parent_id;
        scan_row.tic_coverage        = expl_result.tic_coverage;
        scan_row.exploration_group_id = expl_result.group_id;
        scan_row.exploration_metric  = expl_result.exploration_metric;
        scan_row.variant_index       = expl_result.variant_index;
        scan_row.total_variants      = expl_result.total_variants;
        scan_row.collision_energy    = std::to_string(expl_result.collision_energy);
        scan_row.exploration_score   = expl_result.score;
        scan_row.remaining_ratio     = expl_result.remaining_ratio;
        scan_row.activation_type     = expl_result.activation_type;
        scan_row.reaction_time       = std::to_string(expl_result.reaction_time);
        scan_row.winner_tracking_id  = expl_result.winner_tracking_id;  // "" except the group-completing row
        has_scan_row = true;

        // Identification rows (copy by value before expl_result goes out of scope): primary + winners.
        // DIVERGENCE vs MS2-expl: MS3 gate is fragments-only, and MS3 also emits the winner-batch rows.
        if (!expl_result.identification_result.fragments.empty())
          id_rows.push_back({id_str, 3, 'E', expl_result.ms2_context, expl_result.identification_result});
        for (const auto& row : expl_result.additional_identification_rows)
        {
          if (!row.identification_result.fragments.empty())
            id_rows.push_back({row.tracking_id, 3, 'E', row.ms2_context, row.identification_result});
        }

        exploration_active_.store(exploration_.activeGroupCount() > 0, std::memory_order_release);
        return_code = static_cast<int>(expl_result.commands.size());
      }
      else
      {
        // ===== MS3 (regular): reuse ctx -> deconvolve -> fragment-match -> row =====
        // Reuse the resolved ctx from the top (gate guarantees num_stages >= 2).
        // Do NOT re-resolve here — the pending entry is already gone.
        double precursor_mass = 0.0;
        int precursor_charge = 0;
        if (resolved.has_value() && resolved->num_stages >= 2)
        {
          precursor_charge = resolved->stages[1].charge_state;
          // Pair the fragment charge with the fragment mass (mono_mass_s1), not the MS2-precursor mass.
          // A consistent (mass,charge) precursor caps MS3 sub-fragment charges to the parent (fragZ <= parentZ).
          precursor_mass = resolved->mono_mass_s1;
        }

        int ms3_mass_count = 0;
        if (mzs != nullptr && ints != nullptr && length > 0)
        {
          deconv_.deconvolveMSn(mzs, ints, length, rt_min, precursor_mass, precursor_charge);
          ms3_mass_count = deconv_.hasStoredMS2()
              ? static_cast<int>(deconv_.storedMS2().size()) : 0;
        }

        const DeconvolvedSpectrum* ms3_spec = deconv_.hasStoredMS2() ? &deconv_.storedMS2() : nullptr;

        // Results-row values that must outlive the inner identification block (cached_ms2_ctx/ms3_matches
        // are scoped there, and cache_it is erased below). Defaults below are overwritten on a match.
        std::string ms3_proteoform, ms3_matched_protein;
        int ms3_frag_count = 0, ms3_tag_count = 0;
        float ms3_tic = 0.0f;

        // Identification: look up MS2 context from cache and run fragment matching
        {
          auto cache_it = ms2_context_cache_.find(tracking_id);
          if (cache_it != ms2_context_cache_.end() && ms3_spec != nullptr && !ms3_spec->empty())
          {
            const auto& cached_ms2_ctx = cache_it->second;
            MS3FragmentMatcher::ProteoformContext proto_ctx;
            proto_ctx.region_start = cached_ms2_ctx.start_pos;
            proto_ctx.region_end = cached_ms2_ctx.end_pos;
            proto_ctx.ptm_sites = cached_ms2_ctx.ptm_sites;

            std::vector<const DeconvolvedSpectrum*> spectra = {ms3_spec};
            std::vector<FragmentAnalysis::ProteoformMatch> ms3_matches;
            MS3FragmentMatcher::calibrateAndScore(
              spectra,
              config_.targeting().protein_sequence,
              proto_ctx,
              cached_ms2_ctx.fragment_ion_type,
              cached_ms2_ctx.fragment_ion_index,
              MS3FragmentMatcher::LOOSE_TOLERANCE_PPM,
              config_.level(3).tolerance_ppm,
              &ms3_matches);

            if (!ms3_matches.empty() && !ms3_matches[0].fragments.empty())
            {
              // Copy cached_ms2_ctx (a reference into ms2_context_cache_, erased below) and ms3_matches[0]
              // (inner-scope) BY VALUE into the id-row now — the bottom write must not read them after the erase.
              id_rows.push_back({id_str, 3, 'R', cached_ms2_ctx, ms3_matches[0]});
              // Capture the decision values the engine used so the results row agrees with this match.
              // Render the proteoform with its discovered PTMs via the same renderer the id row uses;
              // cached_ms2_ctx.ptm_sites is the cached parent-MS2 PTM set (toProForma returns the bare sequence when empty).
              ms3_proteoform = FragmentAnalysis::toProForma(cached_ms2_ctx.proteoform_sequence, cached_ms2_ctx.ptm_sites);
              ms3_frag_count = ms3_matches[0].total_match_count;
              ms3_tag_count = cached_ms2_ctx.tag_count;
              ms3_matched_protein = config_.targeting().fasta_file.empty()
                  ? config_.targeting().protein_sequence : config_.targeting().fasta_file;
              if (ms3_mass_count > 0)
              {
                float r = static_cast<float>(ms3_matches[0].total_match_count) / static_cast<float>(ms3_mass_count);
                ms3_tic = r > 1.0f ? 1.0f : r;  // matched-fragment coverage of the MS3 deconvolution
              }
            }

            ms2_context_cache_.erase(cache_it);
          }
        }

        std::string parent_id;
        if (resolved.has_value())
          parent_id = std::string(resolved->parent_scan_id);

        // 2-stage CE/activation/reaction = MS2 isolation stage ; MS3 fragmentation stage.
        std::string ms3_collision_energy = "0", ms3_activation_type = "", ms3_reaction_time = "0";
        if (resolved.has_value() && resolved->num_stages >= 2)
        {
          ms3_collision_energy = std::to_string(resolved->stages[0].collision_energy) + ";" + std::to_string(resolved->stages[1].collision_energy);
          ms3_activation_type  = std::string(resolved->stages[0].activation_type) + ";" + std::string(resolved->stages[1].activation_type);
          ms3_reaction_time    = std::to_string(resolved->stages[0].reaction_time) + ";" + std::to_string(resolved->stages[1].reaction_time);
        }
        // MS3 rows log fragment_count & tag_count as the -1 sentinel:
        //   - fragment_count: MS3 matching is finalized only in the calibrated round; the matched count
        //     lives in identification.tsv (ms3_fragments), not here.
        //   - tag_count: tagging is an MS2 feature, not used for fragment-based MS3 id.
        // Fill scan_results row; commands_pushed, tag_count, fragment_count keep their sentinel defaults
        // (0/-1/-1); child_ids stays empty; exploration fields stay at defaults; tic_coverage stays real.
        scan_row.mass_count         = ms3_mass_count;
        scan_row.matched_protein    = ms3_matched_protein;
        scan_row.proteoform_sequence = ms3_proteoform;
        scan_row.deconv_spectrum    = ms3_spec;
        scan_row.parent_tracking_id = parent_id;
        scan_row.tic_coverage       = ms3_tic;
        scan_row.collision_energy   = ms3_collision_energy;
        scan_row.activation_type    = ms3_activation_type;
        scan_row.reaction_time      = ms3_reaction_time;
        has_scan_row = true;

        return_code = 0;
      }
    }

    // ===== Write the file-backed TSV rows once: scan_results row, then identification rows. =====
    // (IDA log + [TRACK-*] stdout stay inline in their branches.) Per-file row order is observable;
    // cross-file scan<->id order is NOT (separate TSVs, per-file golden compare in LogGoldenComparer).
    // Do not add a cross-file interleaving assertion without revisiting this.
    if (has_scan_row)
      logger_.writeScanResultRow(scan_row);
    for (const auto& r : id_rows)
      logger_.writeIdentificationRow(r);
    return return_code;
  }

  int FLASHIda::getNextScanCommand(ScanCommand& out)
  {
    // No analysis_mutex_ acquired — queue methods lock internally, exploration/FAIMS via atomics
    double faims_cv = faims_.isEnabled() ? current_faims_cv_.load(std::memory_order_acquire) : 0.0;

    // Step 1: AGC scan if needed
    if (queue_.needsAGC())
    {
      out = queue_.makeAGC();
      out.faims_cv = faims_cv;
      out.scan_id = queue_.nextTrackingId();
      queue_.recordAGCTime();

      uint64_t now_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::steady_clock::now().time_since_epoch()).count());
      out.enqueue_timestamp_ms = now_ms;
      out.dequeue_timestamp_ms = now_ms;

      // Scan description: {3-char ID}A
      std::string id_str = ScanCommandQueue::encode(out.scan_id);
      std::snprintf(out.scan_description, 16, "%sA", id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=agc" << std::endl;
      queue_.registerPending(out);
      {
        std::lock_guard<std::mutex> lk(analysis_mutex_);
        logger_.writeScanCommandRow(out);
      }
      return 1;
    }

    // Step 2: Cycle time -- force MS1 if too long since last survey scan
    // Suppressed while any exploration group is active.
    // Queued at priority 0 (not returned immediately) so it goes through normal dequeue.
    if (config_.scheduling().cycle_time_enabled
        && !exploration_active_.load(std::memory_order_acquire)
        && queue_.msSinceLastMS1() > static_cast<uint64_t>(config_.scheduling().cycle_time_ms)
    )
    {
      ScanCommand ms1_cmd = queue_.makeMS1();
      ms1_cmd.faims_cv = faims_cv;
      ms1_cmd.scan_id = queue_.nextTrackingId();
      ms1_cmd.priority = 0;
      ms1_cmd.scan_description[3] = 'C';
      queue_.recordMS1Time();

      std::string id_str = ScanCommandQueue::encode(ms1_cmd.scan_id);
      std::snprintf(ms1_cmd.scan_description, 16, "%sS", id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=cycle_time" << std::endl;
      queue_.push(ms1_cmd);
      // Fall through to Step 3 (cleanup) and Step 4 (dequeue)
    }

    // Step 3: Cleanup expired commands
    queue_.cleanupExpired();

    // Step 4: Dequeue by priority (0 = highest -> 3 = lowest)
    auto dequeued = queue_.dequeue();
    if (dequeued.has_value())
    {
      out = dequeued.value();
      if (out.msn_level == 1 && out.is_agc == 0)
        queue_.recordMS1Time();
      // faims_cv already set at creation time (MS2 -> parent CV, CV-transition MS1 -> next CV)
      {
        std::lock_guard<std::mutex> lk(analysis_mutex_);
        logger_.writeScanCommandRow(out);
      }
      return 1;
    }

    // Step 5: Idle cycle -- queue empty, keep the instrument busy with AGC + MS1
    // Create an AGC command (returned immediately) and push an MS1 at priority 3
    // into the queue as a fallback scan (lowest priority, behind follow-ups/MS3/MS2).
    {
      // 5a: AGC
      ScanCommand agc_cmd = queue_.makeAGC();
      agc_cmd.faims_cv = faims_cv;
      agc_cmd.scan_id = queue_.nextTrackingId();
      queue_.recordAGCTime();

      uint64_t now_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::steady_clock::now().time_since_epoch()).count());
      agc_cmd.enqueue_timestamp_ms = now_ms;
      agc_cmd.dequeue_timestamp_ms = now_ms;

      std::string agc_id_str = ScanCommandQueue::encode(agc_cmd.scan_id);
      std::snprintf(agc_cmd.scan_description, 16, "%sA", agc_id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << agc_id_str << " ms_level=1 type=idle_agc" << std::endl;

      // 5b: MS1 -- use default priority 3 (lowest, behind follow-ups/MS3/MS2)
      ScanCommand ms1_cmd = queue_.makeMS1();
      ms1_cmd.faims_cv = faims_cv;
      ms1_cmd.scan_id = queue_.nextTrackingId();
      ms1_cmd.priority = 3;

      std::string ms1_id_str = ScanCommandQueue::encode(ms1_cmd.scan_id);
      std::snprintf(ms1_cmd.scan_description, 16, "%sS", ms1_id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << ms1_id_str << " ms_level=1 type=idle_ms1" << std::endl;

      // Push MS1 into priority-3 queue for next dequeue call
      queue_.push(ms1_cmd);

      out = agc_cmd;
      queue_.registerPending(out);
      {
        std::lock_guard<std::mutex> lk(analysis_mutex_);
        logger_.writeScanCommandRow(out);
      }
      return 1;
    }
  }

  int FLASHIda::getNextTrackingId()
  {
    return queue_.nextTrackingId();
  }

} // namespace OpenMS
