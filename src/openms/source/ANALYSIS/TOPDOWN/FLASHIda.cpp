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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/NotchSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanRole.h>
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
  namespace
  {
    /// Does this quantification verdict buy the identification scan, under the authored objective?
    /// (ADR-0039, amending ADR-0038's hardcoded `verdict == Differential`.)
    ///
    /// This is the WHOLE of the acquisition decision quantification makes. It deliberately does not
    /// touch the verdict itself: quant_verdict reports the MEASUREMENT, this decides the PURCHASE.
    /// Keeping those apart is what lets the objective change without moving a single log column.
    bool quantBuysIdentification(const Quantification::Result& q, const QuantConfig& cfg)
    {
      using V = Quantification::Result::Verdict;
      switch (cfg.identify)
      {
        case QuantIdentify::None:
          return false;
        case QuantIdentify::All:
          return true;
        case QuantIdentify::Quantified:
          return q.verdict == V::Differential || q.verdict == V::NotDifferential;
        case QuantIdentify::Differential:
          if (q.verdict != V::Differential) { return false; }
          if (cfg.enriched_in < 0) { return true; }  // "either" -- ADR-0038's symmetric behaviour
          // condition_means, NEVER fold_change. A wholly-absent condition is Differential with
          // fold_change == -1.0 (Quantification.cpp) -- a SENTINEL, not a ratio -- so a
          // `fold_change > threshold` test here would silently drop every species present in one
          // condition and absent in the other, which is the strongest result the experiment can
          // produce. Pinned by processScan_enriched_in_survives_a_wholly_absent_condition.
          //
          // A plain `>` is sufficient and correct: this line is reached only when the verdict is
          // already Differential, so the threshold is satisfied in one direction and the larger
          // mean is the enriched condition.
          return q.condition_means[cfg.enriched_in] > q.condition_means[1 - cfg.enriched_in];
      }
      return false;  // unreachable for a valid enumerator
    }

    /// consensus_charges: every balloting charge, with the agreeing ones starred (ADR-0040).
    ///
    /// ';'-joined, matching every other list column. The star is what makes a dissent legible at a
    /// glance -- `11*;14*;17` says which charge was discarded and therefore which population the
    /// identification scan was NOT taken from. Abstentions do not appear: they cast no ballot, so
    /// listing them would read as dissent.
    std::string formatConsensusCharges_(const Quantification::Consensus& c)
    {
      std::string out;
      for (size_t i = 0; i < c.ballot_charges.size(); ++i)
      {
        if (!out.empty()) out += ";";
        out += std::to_string(c.ballot_charges[i]);
        if (i < c.ballot_agreed.size() && c.ballot_agreed[i]) out += "*";
      }
      return out;
    }
  } // namespace

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

    // The monitor scan's own deconvolution engine, built only if a level asks for one (ADR-0042).
    // In the ctor BODY rather than the init list deliberately: it is conditional, and a null pointer
    // is what makes the feature free when off. Same toleranceList() as deconv_, so a monitor reading
    // is numerically comparable with a survey's.
    if (config_.level(2).monitor_ms1_enabled || config_.level(3).monitor_ms1_enabled)
      monitor_deconv_ = std::make_unique<Deconvolution>(config_, config_.toleranceList());

    // Initialize FAIMS CV atomic for getNextScanCommand reads.
    // relaxed is safe: no other thread can observe this object before construction completes.
    current_faims_cv_.store(faims_.isEnabled() ? faims_.currentCV() : 0.0, std::memory_order_relaxed);
  }

  FLASHIda::~FLASHIda() = default;

  /**
   * @brief The MS1 SURVEY arm: deconvolve -> select top-N -> push MS2 (or CE sweep) -> FAIMS cycle.
   *
   * Extracted VERBATIM from processScan's `ms_level == 1` branch (ADR-0042). Nothing here changed;
   * the extraction exists so that the branch can fork into a deciding arm and an observing one with
   * NOTHING AFTER THE FORK -- which is what makes a side effect added to this function unreachable
   * from a monitor scan by construction rather than by a guard someone has to remember.
   *
   * The `selection == None` early return deliberately stays in processScan rather than moving here:
   * it returns without writing a scan_results row, and that asymmetry is documented at the bottom
   * of processScan.
   *
   * @return the number of commands pushed, which becomes processScan's return_code.
   */
  void FLASHIda::observeMonitoringMS1_(Deconvolution& d, const Config& cfg,
                                       const double* mzs, const double* ints, int length,
                                       double rt, double faims_cv,
                                       IdaLogger::ScanRowDescriptor& results_row)
  {
    // The WHOLE observing arm. Read the `static` note on the declaration before adding a line here:
    // nothing in this function can reach an acquisition decision, and that is enforced by the
    // absence of `this` rather than by anyone remembering.
    d.deconvolveMS1(mzs, ints, length, rt, faims_cv);

    // Order the masses the way a survey's row orders them, so a monitor row and a survey row can be
    // read down the same column. This is the RANKING sort only -- the inclusion-priority pass that
    // follows it in PrecursorSelection is a selection decision and is deliberately not called.
    sortByLevelMetric(d.deconvolvedMS1(), cfg, 1);

    results_row.mass_count      = static_cast<int>(d.deconvolvedMS1().size());
    // Engine-owned and outlives processScan (monitor_deconv_ is a FLASHIda member), which is the
    // invariant ScanRowDescriptor::deconv_spectrum documents.
    results_row.deconv_spectrum = &d.deconvolvedMS1();

    // commands_pushed stays 0, child_ids stays empty, and every other field keeps its sentinel --
    // all of which is the truth about a scan that decided nothing.
    //
    // NO ida.log entry, deliberately. That file is the record of acquisition DECISIONS and the
    // FLASHDeconv coupling file; a monitor scan makes none and must not couple. Its entry would read
    // "- 0 targets", which BOTH readers (IdaLogger::parseFLASHIdaLog and PrecursorSelection's
    // target_log_files loader) skip by design -- so it would carry nothing to any consumer while
    // adding noise to the one stream with an outside reader.
  }

  int FLASHIda::runSurveyMS1_(const double* mzs, const double* ints, int length, double rt_min,
                              double faims_cv, int parent_tracking_id, const ScanCommand& parent_ctx,
                              int instrument_scan_number,
                              IdaLogger::ScanRowDescriptor& results_row,
                              std::vector<std::string>& child_ids)
  {
    int n = selection_.filterAndRank(mzs, ints, length, rt_min, 1, faims_cv);
    const auto& selected = selection_.selectedPeakGroups();
    const auto& sel_charges = selection_.triggerCharges();
    int commands_pushed = 0;
    std::vector<ScanCommand> ms2_commands;
    // Index-parallel to ms2_commands; points into selection_.selectedPeakGroups(), a member vector
    // that outlives the writeIDALogEntry call below.
    std::vector<const PeakGroup*> ms2_sources;
    const MSSpectrum& raw_ms1 = selection_.deconvolvedMS1().getOriginalSpectrum();

    // Shared MS2 push: stamp the precursor_id, enqueue, and remember the command for the IDA log.
    // Equal-priority commands drain FIFO, so this preserves scan_commands row order and child_ids.
    // @p src is the PeakGroup this command was built from. Recorded index-parallel to
    // ms2_commands so ida.log's ChargeRange can report the species' real charge envelope without
    // looking the mass back up in the deconvolved spectrum -- several PeakGroups routinely share
    // one mass in a survey (ADR-0036) and each carries a different charge subset, so a lookup
    // would pick one of them arbitrarily. See ADR-0035 decision 6.
    auto push_ms2_command = [&](ScanCommand& c, int precursor_id, const PeakGroup& src)
    {
      stampAndPush_(c, precursor_id);
      ms2_commands.push_back(c);
      ms2_sources.push_back(&src);
      commands_pushed++;
    };

    if (config_.hasExploration(2))
    {
      // Exploration path: initiate CE sweep variants INSTEAD of regular MS2
      for (int i = 0; i < n; i++)
      {
        // One precursor_id per MS1-selected precursor; every CE-sweep variant of this precursor
        // shares it (they are the same MS1 selection -> same charge -> same model).
        const int precursor_id = allocPrecursorId_();
        ScanCommand ms1_ctx{};
        ms1_ctx.scan_id  = parent_tracking_id;
        ms1_ctx.faims_cv = parent_ctx.faims_cv;  // CV travels via the resolved MS1 context
        // Same index i as the production dispatch below, so a species carrying an authored
        // charge set has its CE-sweep variants isolate the charges its production scan would.
        const std::vector<int>* allowed_i =
            i < (int)selection_.triggerAuthoredCharges().size() ? &selection_.triggerAuthoredCharges()[i] : nullptr;
        auto cmds = exploration_.initiate(2, selected[i], sel_charges[i], queue_, &raw_ms1, &ms1_ctx,
                                          '\0', 0, {}, {}, nullptr, nullptr, allowed_i);
        for (auto& c : cmds)
          push_ms2_command(c, precursor_id, selected[i]);
      }
    }
    else
    {
      // Normal path: push MS2 for each precursor, for each scan config
      //
      // ADR-0040: the QUANTIFICATION GROUP is minted here, one per SPECIES per survey, and keyed
      // on NOMINAL MASS rather than on loop adjacency. Adjacency happens to hold today -- the emit
      // loop pushes a species' charges consecutively -- but ADR-0036 siblings put two PeakGroups
      // of one species in this list at slightly different masses (12351.39 / 12351.33 in the
      // committed fixtures), and those are ONE reporter population. Nominal mass is the ~1 Da bin
      // every other acquisition map already keys on, so it groups them and adjacency would not.
      std::map<int, int> quant_group_by_nominal;   // survey-local: nominal mass -> group_id
      for (int i = 0; i < n; i++)
      {
        // One precursor_id per MS1-selected precursor; all of its MS2 scan-config commands share it.
        const int precursor_id = allocPrecursorId_();
        // With quantification on, Config puts the quantification scan in the roster's PRIMARY
        // slot and holds ms_settings.ms2 back as the scan a differential verdict buys, so
        // scans[0] is the 'Q' (ADR-0038). Everything else on the roster is an ordinary 'R'.
        const bool quant_owns_primary =
            config_.quantification().enabled && config_.quantification().has_quant_scan;
        const auto& roster = config_.level(2).scans;
        for (size_t si = 0; si < roster.size(); ++si)
        {
          ScanConfig ms2_config = roster[si];
          const char marker = (quant_owns_primary && si == 0) ? 'Q' : 'R';
          const std::vector<int>* allowed_i =
              i < (int)selection_.triggerAuthoredCharges().size() ? &selection_.triggerAuthoredCharges()[i] : nullptr;
          ScanCommand cmd = queue_.buildMS2(selected[i], sel_charges[i], ms2_config, 2, parent_tracking_id, allowed_i, marker);
          // MS2 inherits the parent MS1's *processing* CV (the faims_cv arg), NOT parent_ctx.faims_cv:
          // in a FAIMS-skip run the MS1 command's stored CV differs from the CV it was processed at.
          cmd.faims_cv = faims_cv;
          if (cmd.num_stages > 0)
          {
            const double half = cmd.stages[0].isolation_width / 2.0;
            cmd.window_snr = FragmentAnalysis::windowSnr(raw_ms1, cmd.stages[0].precursor_mz - half,
                                                         cmd.stages[0].precursor_mz + half, cmd.precursor_intensity);
          }
          push_ms2_command(cmd, precursor_id, selected[i]);

          // Register the 'Q' with its species' group. AFTER the push, because addMember indexes on
          // cmd.scan_id and reads the window_snr/intensity set above -- both are final only here.
          if (marker == 'Q')
          {
            const int nominal = SpectralDeconvolution::getNominalMass(selected[i].getMonoMass());
            auto git = quant_group_by_nominal.find(nominal);
            if (git == quant_group_by_nominal.end())
              git = quant_group_by_nominal.emplace(nominal, allocQuantGroupId_()).first;
            quant_.addMember(git->second, cmd, precursor_id);
          }
        }
      }
    }

    // IDA log entry (MS1 only).
    logger_.writeIDALogEntry(rt_min, parent_tracking_id, instrument_scan_number,
                             ms2_commands, ms2_sources, selection_.deconvolvedMS1());

    for (const auto& c : ms2_commands)
      child_ids.push_back(ScanCommandQueue::encode(c.scan_id));
    int ms1_mass_count = static_cast<int>(selection_.deconvolvedMS1().size());
    results_row.mass_count      = ms1_mass_count;
    results_row.commands_pushed = commands_pushed;
    results_row.child_ids       = child_ids;
    results_row.deconv_spectrum = &selection_.deconvolvedMS1();

    // FAIMS CV cycling: update skip policy, advance to next CV, push MS1.
    // isCycling(), not isEnabled(): with a single configured CV the "next" CV is always the
    // current one, so this would push a redundant priority-0 MS1 after every MS1 and double the
    // survey rate for no gain. A fixed-CV run still gets its CV -- it travels on every command
    // via current_faims_cv_, it just never transitions (ADR-0012).
    if (faims_.isCycling())
    {
      double current_cv = faims_.currentCV();
      faims_.updateSkip(current_cv, commands_pushed);

      double next_cv = faims_.advanceToNextCV();
      ScanCommand ms1 = queue_.makeMS1();
      ms1.faims_cv = next_cv;
      ms1.scan_id = queue_.nextTrackingId();
      ms1.priority = 0;  // priority 0 to send before pending MS2s

      std::string next_ms1_id = ScanCommandQueue::encode(ms1.scan_id);
      ScanCommandQueue::writeMs1Description(ms1, ScanRole::Survey);

      queue_.push(ms1);
      std::cout << "[TRACK-CREATE] id=" << next_ms1_id << " ms_level=1 type=cv_transition cv=" << next_cv << std::endl;
    }

    // Update atomics for lock-free reads by getNextScanCommand
    publishExplorationState_();
    current_faims_cv_.store(faims_.isEnabled() ? faims_.currentCV() : 0.0, std::memory_order_release);

    return commands_pushed;
  }

  int FLASHIda::processScan(const double* mzs, const double* ints, int length,
                             double rt_min, int ms_level, const char* scan_description,
                             double faims_cv, int instrument_scan_number)
  {
    std::lock_guard<std::mutex> lock(analysis_mutex_);

    // Tracking ID extraction (perspective: scan processed is parent of potential children)
    std::string parent_desc_str = scan_description ? std::string(scan_description) : "";
    if (parent_desc_str.size() < 3) return 0;
    std::string parent_id_str = parent_desc_str.substr(0, 3);
    int parent_tracking_id = queue_.decode(parent_id_str);

    // THE ECHO's role, and the ONLY gate entitled to read it (ADR-0042).
    //
    // Two copies of the role exist and they are not interchangeable. This one is decoded from the
    // description the INSTRUMENT handed back, because this gate has to run before resolvePending()
    // below and there is no engine record to consult yet. Every DECISION gate reads `cmd_role`,
    // decoded from the engine's own queued command, because that copy cannot be truncated, spoofed
    // or lost by the trailer round trip. ADR-0008/0035 separate the identity channels but never said
    // which copy of the ROLE wins; it is the engine's, everywhere except right here.
    const std::optional<ScanRole> echo_role =
        roleFromDescription(parent_desc_str.c_str(), parent_desc_str.size());

    // Skip AGC scans — calibration-only, no data to process.
    // `Agc` is the only row in kScanRoleTraits with observes == false, and an undecodable
    // description resolves as observing, so this is byte-identical to the `desc[3] == 'A'` test it
    // replaces — while also being correct for any future peakless role.
    if (!roleObserves(echo_role))
    {
      queue_.resolvePending(parent_tracking_id);
      return 0;
    }

    // Retrieve timestamps from pending map
    uint64_t enqueue_ts = 0;
    uint64_t dequeue_ts = 0;
    {
      auto peeked = queue_.peekPending(parent_tracking_id);
      if (peeked.has_value())
      {
        enqueue_ts = peeked->enqueue_timestamp_ms;
        dequeue_ts = peeked->dequeue_timestamp_ms;
      }
      else
      {
        std::cout << "[TRACK-RESOLVE] id=" << parent_id_str << " status=not_found" << std::endl;
        return 0;
      }
    }

    // Stamp received timestamp (instrument → C++ handoff)
    uint64_t received_ts = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()
      ).count()
    );

    auto resolved = queue_.resolvePending(parent_tracking_id);
    if (!resolved.has_value())
    {
      // peekPending above normally guarantees this; reaching here means an upstream race.
      std::cout << "[TRACK-RESOLVE] id=" << parent_id_str << " status=context_lost_race" << std::endl;
      return 0;
    }
    ScanCommand parent_ctx = resolved.value();

    // Context-support gate, validated once before branching. required_stages = the stage count this
    // level needs in its queued command; ms_level>3 is unsupported (returns 0).
    const int required_stages = (ms_level == 1) ? 0 : (ms_level == 2) ? 1 : (ms_level == 3) ? 2 : -1;
    if (required_stages < 0 || parent_ctx.msn_level != ms_level || parent_ctx.num_stages < required_stages)
    {
      std::cout << "[TRACK-RESOLVE] id=" << parent_id_str << " status=context_unsupported"
                << " ms_level=" << ms_level << " ctx_msn=" << parent_ctx.msn_level
                << " num_stages=" << parent_ctx.num_stages << std::endl;
      return 0;
    }

    // Locals shared by all levels: the two scan-role flags and the child-id accumulator.
    //
    // 'F' is gone with ADR-0038. The two roles it used to conflate are now separate marks:
    //   'C'  a tagging follow-up  -- bought by tag detection, must not buy another (depth one)
    //   'Q'  the QUANTIFICATION scan -- rostered, primary, and the only scan ever measured
    // The identification scan a differential verdict buys is an ordinary 'R', which is what keeps
    // 'R' meaning "identification MS2" in every mode. It is safe for it to be indistinguishable
    // from a rostered MS2 precisely because the screen gate below keys on 'Q' rather than on
    // "is not a follow-up": an 'R' is never re-screened, so it can never buy a second scan.
    // THE ENGINE RECORD's role -- the copy every DECISION gate reads, and the one the MS1 fork below
    // dispatches on. Not `echo_role`: see the note at its declaration. Comparing an engaged optional
    // to a ScanRole is exactly the old `len >= 4 && desc[3] == 'C'` test, since an undecodable
    // description yields nullopt and nullopt equals no role.
    const std::optional<ScanRole> cmd_role = roleOf(parent_ctx);
    const bool is_follow_up_scan = (cmd_role == ScanRole::FollowUp);
    const bool is_quant_scan     = (cmd_role == ScanRole::Quantification);
    std::vector<std::string> child_ids;

    // Carriers for the file-backed TSV rows (scan_results + identification): each branch fills these,
    // and they are written once after the branch dispatch below. The IDA log, the [TRACK-*] stdout
    // lines, and all side effects (queue pushes, atomics, cache erase) stay inline in each branch.
    int return_code = 0;
    IdaLogger::ScanRowDescriptor results_row{};
    std::vector<IdaLogger::IdRowDescriptor> id_rows;

    // scan_results fields that are identical in every branch — set once here. Each branch fills its
    // own branch-specific fields; the rest keep their results_row{} sentinel defaults.
    results_row.tracking_id = parent_id_str;
    results_row.ms_level    = ms_level;   // gate guarantees ms_level in {1,2,3}
    results_row.rt          = rt_min;
    results_row.enqueue_ts  = enqueue_ts;
    results_row.dequeue_ts  = dequeue_ts;
    results_row.received_ts = received_ts;

    if (ms_level == 1)
    {
      // ===== THE FORK (ADR-0042) =====
      // An MS1 either DECIDES or it OBSERVES, and the two arms share nothing below this point.
      //
      // The observing arm is tested FIRST and deliberately sits above the selection == None gate:
      // observation is not gated on a selection metric, and that gate returns before any row is
      // written. It reads cmd_role -- the ENGINE's record of what it asked for -- not the atomic, so
      // a monitor scan that returns after its sweep has already ended is still recognised as one.
      // Re-reading exploration_active_ here would make behaviour depend on when the instrument got
      // round to answering.
      if (!roleDecides(cmd_role))
      {
        if (monitor_deconv_)
          observeMonitoringMS1_(*monitor_deconv_, config_, mzs, ints, length, rt_min, faims_cv,
                                results_row);
        // return_code stays 0, which the host already treats as an ordinary gate rejection. No
        // bridge or C# contract moves.
      }
      // ===== MS1 SURVEY: deconvolve -> select top-N -> push MS2 (or CE sweep) -> FAIMS cycle =====
      // Selection=none: skip MS1 precursor selection entirely
      else if (config_.level(1).selection == SelectionMetric::None)
      {
        return 0;
      }
      else
      {
        return_code = runSurveyMS1_(mzs, ints, length, rt_min, faims_cv,
                                    parent_tracking_id, parent_ctx, instrument_scan_number,
                                    results_row, child_ids);
      }
      // NOTHING BELOW THIS LINE, and that is the seam. A monitor scan cannot acquire a decision by
      // someone adding one here, because there is no here -- any new MS1 side effect has to go
      // inside one arm or the other, which forces the author to say which. The shared tail below
      // (writeScanResultRow) runs for both arms and is the only thing they have in common.
    }
    else if (ms_level == 2)
    {
      int commands_pushed = 0;

      // Precursor ID for grounding protoeoform model
      const int precursor_id = precursorIdForTracking_(parent_tracking_id);

      // ===== MS2 (exploration variant): score variant -> push children -> fill row =====
      if (exploration_.isExplorationVariant(parent_tracking_id))
      {
        auto expl_result = exploration_.feedResult(parent_tracking_id, mzs, ints, length, rt_min, queue_, &tracker_, precursor_id);
        for (auto& c : expl_result.commands)
        {
          stampAndPush_(c, precursor_id);
        }
        // Seed the parent-MS2 context of every MS3 this group dispatched, so each one identifies when it
        // returns on the regular MS3 path. Empty unless the group cascaded to MS3.
        for (const auto& [ms3_scan_id, parent_ms2_context] : expl_result.ms2_context_cache)
          ms2_context_cache_[ms3_scan_id] = parent_ms2_context;

        int expl_mass_count = exploration_.explorationDeconvMassCount();
        const DeconvolvedSpectrum* expl_spec = exploration_.explorationDeconvSpectrum();
        std::string parent_id(expl_result.parent_scan_id);

        // Identification yield for an MS2 exploration VARIANT. Every MS2 variant is identified under
        // every metric (computeFragmentMatch_'s whole-protein matcher is correct at MS2), so these are
        // real readings, gated only on a sequence being configured.
        if (!config_.characterization().protein_sequence.empty())
        {
          results_row.tag_count      = expl_result.identification.result.tag_count;
          results_row.fragment_count = expl_result.identification.result.total_match_count;
          results_row.tic_coverage   = expl_result.metric.tic_coverage;
        }
        // Fill scan_results row in canonical order (copy expl_result.* by value — it dies at block end).
        results_row.mass_count          = expl_mass_count;
        results_row.commands_pushed     = static_cast<int>(expl_result.commands.size());
        results_row.child_ids           = expl_result.child_ids;
        results_row.deconv_spectrum     = expl_spec;
        results_row.parent_tracking_id  = parent_id;
        results_row.exploration_group_id = expl_result.group.group_id;
        results_row.exploration_metric  = expl_result.metric.exploration_metric;
        results_row.variant_index       = expl_result.group.variant_index;
        results_row.total_variants      = expl_result.group.total_variants;
        results_row.collision_energy    = std::to_string(expl_result.fragmentation.collision_energy);
        results_row.exploration_score   = expl_result.metric.score;
        results_row.remaining_ratio     = expl_result.metric.remaining_ratio;
        results_row.activation_type     = expl_result.fragmentation.activation_type;
        results_row.reaction_time       = std::to_string(expl_result.fragmentation.reaction_time);
        results_row.winner_tracking_id  = expl_result.group.winner_tracking_id;  // "" except the group-completing row

        // MS2 identification row (copy by value before expl_result goes out of scope).
        if (!expl_result.identification.proteoform_sequence.empty() && !expl_result.identification.result.fragments.empty())
          id_rows.push_back({parent_id_str, 2, 'E', expl_result.identification.ms2_context, expl_result.identification.result, precursor_id, expl_result.metric.tic_coverage, expl_result.identification.result.score});
        // C: winner-re-matched non-winner identification rows (MS2 'E'), each rendering the winner proteoform;
        // flash_extender_score = -1 marks "no own ID". Populated by Exploration after MS2 finalize.
        for (const auto& row : expl_result.identification.additional_rows)
          if (!row.identification_result.fragments.empty())
            id_rows.push_back({row.tracking_id, 2, 'E', row.ms2_context, row.identification_result, precursor_id, row.tic_coverage, row.flash_extender_score});

        publishExplorationState_();
        return_code = static_cast<int>(expl_result.commands.size());
      }
      else
      {
        // ===== MS2 (regular): deconvolve -> targeting + follow-ups -> MS3 targeting -> row =====
        double precursor_mass = parent_ctx.mono_mass;
        // The HIGHEST charge this scan actually isolated, not the anchor's (ADR-0016). With a
        // co-isolated charge set every member is genuinely present in the isolation, so a fragment of
        // the highest member may itself carry that charge; capping at the anchor would discard real
        // fragments. Identical to stages[0].charge_state whenever there are no notches.
        int precursor_charge = maxIsolatedCharge(parent_ctx, 0);
        deconv_.deconvolveMSn(mzs, ints, length, rt_min, precursor_mass, precursor_charge);

        // Tag-based targeting. processMS2ForTagBasedTargeting returns how many tags matched the FASTA
        // target database, but nothing consumes that as a COUNT -- it is a gate, and the only question
        // asked of it is whether a target protein was detected. It is also 0 whenever tags existed but
        // matched no DB protein, which is why it is not the number any log reports; the identification
        // tag count on ProteoformMatch is, and that one is taken before any protein is consulted.
        // Not on a quantification scan (ADR-0038). Its activation and energy were chosen to release
        // a reporter ion, not to fragment informatively, so tags read off it would be read from the
        // wrong spectrum -- and a 'C' follow-up bought here would be a third scan on the precursor
        // before the identification scan has even run. Skipping the detection also skips the work.
        bool tags_found = false;
        if (selection_.hasTargetProteinDatabase() && precursor_mass > 0 && !is_quant_scan)
        {
          std::string ms2_activation = std::string(parent_ctx.stages[0].activation_type);
          tags_found = selection_.processMS2ForTagBasedTargeting(precursor_mass, ms2_activation) > 0;
        }

        // --- Quantification: measure THIS scan, and buy the identification scan (ADR-0038) -------
        // Gated on the scan's own mark, not on "is not a follow-up". Only a 'Q' is measured, so the
        // identification scan it buys -- an ordinary 'R' -- is never re-screened and depth stays one
        // without the marker having to distinguish a follow-up from a primary.
        if (is_quant_scan && config_.quantification().enabled)
        {
          const auto q = quant_.measure(mzs, ints, length, rt_min, 2, "ms2_quant");

          // Every field is logged, so what scan_results reports IS what the engine decided on.
          results_row.quant_verdict = Quantification::verdictName(q.verdict);
          results_row.quant_fold_change = q.fold_change;
          results_row.quant_channels = q.channels;
          results_row.quant_condition_means = {q.condition_means[0], q.condition_means[1]};

          // ADR-0040. The buy moved from PER-SCAN to PER-GROUP. Under precursor_charges "single"
          // the group has one member, completes on this first return, and the consensus of one
          // measurement IS that measurement -- so this path is byte-identical to ADR-0039's.
          const int gid = quant_.deposit(parent_tracking_id, q);
          if (gid > 0 && quant_.isComplete(gid))
          {
            const auto cons = quant_.decide(gid);

            // The consensus belongs to the GROUP, not to any one scan, so it is written on
            // whichever 'Q' completed it and left empty everywhere else. consensus_id_charge is
            // carried explicitly so a reader never has to infer that this row is not the row the
            // decision was about.
            results_row.consensus_verdict     = Quantification::verdictName(cons.result.verdict);
            results_row.consensus_fold_change = cons.result.fold_change;
            results_row.consensus_agreement   = std::to_string(cons.agreeing) + "/"
                                                + std::to_string(cons.total_ballots);
            results_row.consensus_charges     = formatConsensusCharges_(cons);
            results_row.consensus_id_charge   = cons.winner_charge;

            // ADR-0039's predicate, UNCHANGED: it receives the synthesized consensus Result, so
            // `identify` and `enriched_in` keep working -- including the condition_means-not-
            // fold_change trap, because the aggregate carries real means.
            if (cons.has_winner && quantBuysIdentification(cons.result, config_.quantification()))
            {
              // BOTH the geometry and the precursor_id come from the WINNER, never from the scan
              // currently returning. The group completes on its LAST member, which may well be a
              // dissenter: building from the winner's command while stamping the dissenter's
              // precursor_id would fragment one charge and file the identification -- and every MS3
              // target downstream of it -- under another charge's proteoform model.
              auto ident = queue_.buildFollowUp(cons.winner_cmd,
                                                config_.quantification().identification_scan, 'R');
              stampAndPush_(ident, cons.winner_precursor_id);
              child_ids.push_back(ScanCommandQueue::encode(ident.scan_id));
              commands_pushed++;
            }
            quant_.closeGroup(gid);
          }
        }

        // Conditional MS2 follow-up -- only when tags detected AND explicitly enabled
        if (config_.targeting().conditional_ms2_enabled && tags_found && !is_follow_up_scan)
        {
          auto cond = queue_.buildFollowUp(parent_ctx, config_.targeting().tagging_follow_up_scan, 'C');
          stampAndPush_(cond, precursor_id);  // P5: follow-up inherits the MS2's precursor_id
          child_ids.push_back(ScanCommandQueue::encode(cond.scan_id));
          commands_pushed++;
        }

        // Identification + MS3 targeting, gated by the SEQUENCE rather than by characterization.mode.
        // initiateNextLevel owns both, and it decides internally which of the two it performs: with mode
        // off it matches fragments and returns no commands, with mode on it also dispatches MS3. Gating
        // the CALL on level(2).selection (i.e. on mode) meant an MS2 in a mode: off run was never
        // identified at all -- no sequence tags, no fragments, no identification.tsv row.
        // The identification-row write below sits OUTSIDE this block deliberately and reads
        // ms3_targeting, which stays default-constructed when the call is skipped.
        // Also skipped on a quantification scan (ADR-0038): MS3 targets picked off a reporter-release
        // spectrum would be selected from a spectrum optimised for something else. Identification and
        // MS3 both belong to the 'R' identification scan, which is the one chosen to fragment
        // informatively -- and which this Q scan may be about to buy.
        Exploration::NextLevelResult ms3_targeting;
        if (!config_.characterization().protein_sequence.empty() && !is_quant_scan)
        {
          // initiateNextLevel stamps each command's parent with the MS2 id (parent_ctx.scan_id) at creation.
          // It serves BOTH MS3 shapes, chosen by config_.hasExploration(3): with MS3 exploration off it
          // builds one fixed-CE MS3 per target straight from ms_settings.ms3; with it on it builds a CE
          // sweep. The `exploration_` receiver names the owning object, not a precondition on the scan.
          ms3_targeting = exploration_.initiateNextLevel(2, deconv_.storedMS2(), parent_ctx.faims_cv, queue_, &parent_ctx, &tracker_, precursor_id);
          for (auto& c : ms3_targeting.commands)
          {
            stampAndPush_(c, precursor_id);  // P5: MS3 children inherit the parent MS2's precursor_id
            child_ids.push_back(ScanCommandQueue::encode(c.scan_id));
            commands_pushed++;
          }

          // Cache MS2 context for each MS3 command (for non-exploration identification.tsv)
          for (size_t ci = 0; ci < ms3_targeting.commands.size(); ++ci)
          {
            if (ms3_targeting.commands[ci].msn_level >= 3 && ci < ms3_targeting.ms3_contexts.size())
            {
              ms2_context_cache_[ms3_targeting.commands[ci].scan_id] = ms3_targeting.ms3_contexts[ci];
            }
          }
        }

        // Write MS2 identification row if proteoform was matched.
        // Gate on the actual match (proteoform_match), NOT on ms3_targeting.proteoform_sequence — the
        // latter is only set when found>0 (MS3 targets were selected), so an unambiguous ID that yields
        // zero MS3 targets would otherwise be dropped from identification.tsv.
        if (!ms3_targeting.proteoform_match.proteoform_sequence.empty() && !ms3_targeting.proteoform_match.fragments.empty())
        {
          Exploration::MS2Context ms2_ctx;
          ms2_ctx.proteoform_sequence = ms3_targeting.proteoform_match.proteoform_sequence;
          ms2_ctx.start_pos = ms3_targeting.proteoform_match.region_start;
          ms2_ctx.end_pos = ms3_targeting.proteoform_match.region_end;
          ms2_ctx.ptm_sites = ms3_targeting.proteoform_match.ptm_sites;
          ms2_ctx.ms1_precursor_mass = parent_ctx.mono_mass;
          ms2_ctx.ms1_precursor_mz = parent_ctx.stages[0].precursor_mz;
          ms2_ctx.ms1_precursor_charge = parent_ctx.stages[0].charge_state;
          // MS2 isolation-window reporting from the resolved command (window-SNR was recorded at MS1
          // time, keyed by this scan_id). The MS3 triplet stays 0 on MS2 rows.
          if (parent_ctx.num_stages > 0)
          {
            ms2_ctx.ms2_isolation_width = parent_ctx.stages[0].isolation_width;
            ms2_ctx.ms2_charge_intensity = parent_ctx.precursor_intensity;
            ms2_ctx.ms2_window_snr = parent_ctx.window_snr;  // Item 2: SNR travels on the command (was the queue map)
          }
          id_rows.push_back({parent_id_str, 2, 'R', ms2_ctx, ms3_targeting.proteoform_match, precursor_id, ms3_targeting.tic_coverage, ms3_targeting.proteoform_match.score});  // copy by value (ms2_ctx is local)
        }

        int ms2_mass_count = deconv_.hasStoredMS2() ? static_cast<int>(deconv_.storedMS2().size()) : 0;
        const DeconvolvedSpectrum* ms2_spec = deconv_.hasStoredMS2() ? &deconv_.storedMS2() : nullptr;
        std::string parent_id(parent_ctx.parent_scan_id);

        // Log the actual MS2 scan-config CE/activation/reaction (single stage) the engine used.
        std::string ms2_collision_energy = parent_ctx.num_stages > 0 ? std::to_string(parent_ctx.stages[0].collision_energy) : "0";
        std::string ms2_activation_type = parent_ctx.num_stages > 0
            ? std::string(parent_ctx.stages[0].activation_type)
            : config_.level(2).scans[0].activation;
        std::string ms2_reaction_time = parent_ctx.num_stages > 0 ? std::to_string(parent_ctx.stages[0].reaction_time) : "0";
        // Identification yield. The DISTINCTION the sentinels exist for: -1 means the tagger never ran
        // on this spectrum, 0 means it ran and read nothing. Identification runs exactly when a protein
        // sequence is configured, so that is the condition -- NOT whether it happened to match anything,
        // which is the reading 0 carries. ms3_targeting.tic_coverage is written only inside its own
        // found>0 guard and otherwise stays 0.0f, so it needs the same treatment or a spectrum that
        // matched nothing would be indistinguishable from one never examined.
        if (!config_.characterization().protein_sequence.empty())
        {
          results_row.tag_count      = ms3_targeting.proteoform_match.tag_count;
          results_row.fragment_count = ms3_targeting.proteoform_match.total_match_count;
          results_row.tic_coverage   = ms3_targeting.tic_coverage;
        }
        // Fill scan_results row; exploration fields stay at their sentinel defaults.
        results_row.mass_count         = ms2_mass_count;
        results_row.commands_pushed    = commands_pushed;
        results_row.child_ids          = child_ids;
        results_row.deconv_spectrum    = ms2_spec;
        results_row.parent_tracking_id = parent_id;
        results_row.collision_energy   = ms2_collision_energy;
        results_row.activation_type    = ms2_activation_type;
        results_row.reaction_time      = ms2_reaction_time;

        std::cout << "[TRACK-RESOLVE] id=" << parent_id_str
                  << " rt=" << rt_min
                  << " commands=" << commands_pushed
                  << std::endl;

        // Update atomic for lock-free reads by getNextScanCommand
        publishExplorationState_();

        return_code = commands_pushed;
      }
    }
    else if (ms_level == 3)
    {
      // Precursor ID for grounding protoeoform model
      const int precursor_id = precursorIdForTracking_(parent_tracking_id);

      // ===== MS3 (exploration variant): score variant -> push children -> fill row =====
      if (exploration_.isExplorationVariant(parent_tracking_id))
      {
        auto expl_result = exploration_.feedResult(parent_tracking_id, mzs, ints, length, rt_min, queue_, &tracker_, precursor_id);
        for (auto& c : expl_result.commands)
        {
          stampAndPush_(c, precursor_id);
        }
        // Seed the parent-MS2 context of the production MS3 a completed MS3 group re-acquires (overrides
        // set); that scan is not a variant, so it returns on the regular MS3 path and needs the lookup.
        for (const auto& [ms3_scan_id, parent_ms2_context] : expl_result.ms2_context_cache)
          ms2_context_cache_[ms3_scan_id] = parent_ms2_context;

        // Retire this variant's own cache entry (a different key from the ones just seeded above, which
        // belong to newly built commands). An MS3 exploration variant is dispatched through the same
        // initiateNextLevel path as a regular MS3, so it was given a parent-MS2 context — but it returns
        // HERE, and this branch never reads the cache. Without the erase the entry lives for the rest of
        // the run: nothing else is keyed on it and only the regular MS3 branch erases. Matches the
        // established one-shot policy (a context is consumed by the scan it was minted for, so a
        // duplicate or late MS3 yields no row).
        ms2_context_cache_.erase(parent_tracking_id);

        int expl_mass_count = exploration_.explorationDeconvMassCount();
        const DeconvolvedSpectrum* expl_spec = exploration_.explorationDeconvSpectrum();
        std::string parent_id(expl_result.parent_scan_id);

        // MS3 exploration variant. tag_count/fragment_count stay at -1: an MS3 spectrum is scored against
        // a known ladder rather than re-tagged for identity, so a count here would describe the
        // sub-fragment spectrum, not the precursor. tic_coverage IS real on MS3 and is reported.
        results_row.tic_coverage        = expl_result.metric.tic_coverage;
        // Logging
        results_row.mass_count          = expl_mass_count;
        results_row.commands_pushed     = static_cast<int>(expl_result.commands.size());
        results_row.child_ids           = expl_result.child_ids;
        results_row.deconv_spectrum     = expl_spec;
        results_row.parent_tracking_id  = parent_id;
        results_row.exploration_group_id = expl_result.group.group_id;
        results_row.exploration_metric  = expl_result.metric.exploration_metric;
        results_row.variant_index       = expl_result.group.variant_index;
        results_row.total_variants      = expl_result.group.total_variants;
        // MS3-exploration rows always log the 2-stage "ms2;ms3" CE/activation/reaction (stage[0] MS2
        // isolation ; stage[1] MS3 fragmentation). stage0_* is always populated on the MS3 branch
        // (buildMS3 emits 2 stages with the real MS2 params), so the single-stage form is unreachable.
        results_row.collision_energy = std::to_string(expl_result.fragmentation.stage0_collision_energy) + ";" + std::to_string(expl_result.fragmentation.collision_energy);
        results_row.activation_type  = expl_result.fragmentation.stage0_activation_type + ";" + expl_result.fragmentation.activation_type;
        results_row.reaction_time    = std::to_string(expl_result.fragmentation.stage0_reaction_time) + ";" + std::to_string(expl_result.fragmentation.reaction_time);
        results_row.exploration_score   = expl_result.metric.score;
        results_row.remaining_ratio     = expl_result.metric.remaining_ratio;
        results_row.winner_tracking_id  = expl_result.group.winner_tracking_id;

        // Identification rows
        if (!expl_result.identification.result.fragments.empty())
          id_rows.push_back({parent_id_str, 3, 'E', expl_result.identification.ms2_context, expl_result.identification.result, precursor_id, expl_result.metric.tic_coverage, expl_result.identification.result.score});
        for (const auto& row : expl_result.identification.additional_rows)
        {
          // Winner-batch rows carry each variant's OWN tic (row.tic_coverage), not the group-completing
          // scan's, so every MS3-'E' identification row reports the tic of the scan it describes.
          if (!row.identification_result.fragments.empty())
            id_rows.push_back({row.tracking_id, 3, 'E', row.ms2_context, row.identification_result, precursor_id, row.tic_coverage, row.flash_extender_score});
        }

        publishExplorationState_();
        return_code = static_cast<int>(expl_result.commands.size());
      }
      else
      {
        // ===== MS3 (regular): reuse parent_ctx -> deconvolve -> fragment-match -> row =====
        // Reuse the resolved parent_ctx from the top (gate guarantees num_stages >= 2).
        // Do NOT re-resolve here — the pending entry is already gone.
        // The context-support gate (top of processScan) guarantees parent_ctx.num_stages >= 2 here.
        // Highest charge isolated at the MS3 stage, so a co-isolated fragment charge set keeps its
        // sub-fragments rather than having them capped at the anchor (ADR-0016). Equals
        // stages[1].charge_state when there are no notches.
        const int precursor_charge = maxIsolatedCharge(parent_ctx, 1);
        // Pair the fragment charge with the fragment mass (mono_mass_s1), not the MS2-precursor mass.
        // A consistent (mass,charge) precursor caps MS3 sub-fragment charges to the parent (fragZ <= parentZ).
        const double precursor_mass = parent_ctx.mono_mass_s1;

        int ms3_mass_count = 0;
        if (mzs != nullptr && ints != nullptr && length > 0)
        {
          deconv_.deconvolveMSn(mzs, ints, length, rt_min, precursor_mass, precursor_charge);
          ms3_mass_count = deconv_.hasStoredMS2()
              ? static_cast<int>(deconv_.storedMS2().size()) : 0;
        }

        const DeconvolvedSpectrum* ms3_spec = deconv_.hasStoredMS2() ? &deconv_.storedMS2() : nullptr;

        // ms3_tic must outlive the inner identification block (ms3_matches is scoped there, cache_it is erased
        // below); it feeds the identification row (moved from scan_results).
        float ms3_tic = 0.0f;

        // Identification: look up MS2 context from cache and run fragment matching
        {
          auto cache_it = ms2_context_cache_.find(parent_tracking_id);
          if (cache_it != ms2_context_cache_.end() && ms3_spec != nullptr && !ms3_spec->empty())
          {
            const auto& cached_ms2_ctx = cache_it->second;
            // Score the returning MS3 against the LIVE WINNER proteoform held by the tracker
            // (ADR-0002), not this triggering MS2 scan's cached context. fragment_ion_type/index still
            // come from cached_ms2_ctx below. Empty context (no finalized winner) => MS3 matches nothing.
            MS3FragmentMatcher::ProteoformContext proto_ctx = tracker_.buildWinnerProteoformContext(precursor_id);

            std::vector<const DeconvolvedSpectrum*> spectra = {ms3_spec};
            std::vector<FragmentAnalysis::ProteoformMatch> ms3_matches;
            ProteoformTracker::scoreCalibratedVariants(
              spectra,
              config_.characterization().protein_sequence,
              proto_ctx,
              cached_ms2_ctx.fragment_ion_type,
              cached_ms2_ctx.fragment_ion_index,
              MS3FragmentMatcher::LOOSE_TOLERANCE_PPM,
              config_.level(3).tolerance_ppm,
              &ms3_matches);
            
            // logging if fragment ions are found.
            if (!ms3_matches.empty() && !ms3_matches[0].fragments.empty())
            {
              if (ms3_mass_count > 0)
              {
                float r = static_cast<float>(ms3_matches[0].total_match_count) / static_cast<float>(ms3_mass_count);
                ms3_tic = r > 1.0f ? 1.0f : r;
              }
              id_rows.push_back({parent_id_str, 3, 'R', cached_ms2_ctx, ms3_matches[0], precursor_id, ms3_tic, ms3_matches[0].score});
              // parent_ctx, not resolved-> : `resolved` is unconditionally non-empty here (the top of
              // processScan returns on an empty resolve) and parent_ctx is its value. num_stages >= 2 is
              // likewise guaranteed by the context-support gate for ms_level 3, and ms3_spec != nullptr by
              // the enclosing if. Those three guards were always true, and reading them the other way --
              // as if a stage-less MS3 could reach this point -- is what they cost.
              {
                Ms2Params parent_ms2_params{parent_ctx.stages[0].collision_energy,
                                            std::string(parent_ctx.stages[0].activation_type),
                                            parent_ctx.stages[0].reaction_time};
                tracker_.feedScan(precursor_id, 3, parent_ms2_params, parent_tracking_id, *ms3_spec,
                                  ms3_matches[0], ms3_matches[0].score, parent_ctx);
                if (cached_ms2_ctx.fragment_ion_type != '\0')
                {
                  const std::string trig = std::string(1, cached_ms2_ctx.fragment_ion_type)
                                         + std::to_string(cached_ms2_ctx.fragment_ion_index);
                  tracker_.foldMs3(precursor_id, trig, parent_id_str);
                }
              }
            }

            ms2_context_cache_.erase(cache_it);
          }
        }

        std::string parent_id(parent_ctx.parent_scan_id);

        // 2-stage CE/activation/reaction = MS2 isolation stage ; MS3 fragmentation stage. Always the
        // two-stage form: the context-support gate rejects an ms_level-3 scan whose command carries
        // fewer than 2 stages, so the single-stage "0"/""/"0" fallback this used to keep was dead.
        std::string ms3_collision_energy = std::to_string(parent_ctx.stages[0].collision_energy) + ";" + std::to_string(parent_ctx.stages[1].collision_energy);
        std::string ms3_activation_type  = std::string(parent_ctx.stages[0].activation_type) + ";" + std::string(parent_ctx.stages[1].activation_type);
        std::string ms3_reaction_time    = std::to_string(parent_ctx.stages[0].reaction_time) + ";" + std::to_string(parent_ctx.stages[1].reaction_time);
        // Fill scan_results row (pure acquisition event): identification payload moved off scan_results —
        // the MS3 fragment proteoform is on scan_commands.ms3_proteoform; matched_protein/tic on identification.
        // Regular MS3. tag_count/fragment_count stay at -1 for the same reason as the exploration MS3
        // row above: no tagger runs on an MS3 spectrum. ms3_tic is this scan's real matched-fragment
        // ratio and is declared outside the identification block precisely so it can be read here.
        results_row.tic_coverage       = ms3_tic;
        // commands_pushed keeps its 0 default; child_ids empty; exploration fields at defaults.
        results_row.mass_count         = ms3_mass_count;
        results_row.deconv_spectrum    = ms3_spec;
        results_row.parent_tracking_id = parent_id;
        results_row.collision_energy   = ms3_collision_energy;
        results_row.activation_type    = ms3_activation_type;
        results_row.reaction_time      = ms3_reaction_time;

        return_code = 0;
      }
    }

    // ===== Write the file-backed TSV rows once: scan_results row, then identification rows. =====
    // (IDA log + [TRACK-*] stdout stay inline in their branches.) Per-file row order is observable;
    // cross-file scan<->id order is NOT (separate TSVs, per-file golden compare in LogGoldenComparer).
    // Do not add a cross-file interleaving assertion without revisiting this.
    // Every scan that reaches here was dispatched into a branch that filled results_row; the only
    // level that leaves without a row (MS1 selection==None) returned above, so write unconditionally.
    logger_.writeScanResultRow(results_row);
    for (const auto& r : id_rows)
      logger_.writeIdentificationRow(r);
    return return_code;
  }

  int FLASHIda::getNextScanCommand(ScanCommand& out)
  {
    // This runs on the INSTRUMENT EVENT THREAD and acquires analysis_mutex_ nowhere. That is the
    // contract, not an accident: processScan holds that lock across a whole MS1 deconvolution, so any
    // acquisition on this path parks the instrument for the duration. Pinned by
    // FLASHIda_ProcessScan_test::drain_completes_while_analysis_mutex_held, which covers all three
    // emitting paths separately because each used to take the lock at its own site.
    //
    // How each piece of cross-thread state is reached instead:
    //   queue_                    -- ScanCommandQueue locks queue_mutex_ internally
    //   precursor_id_by_tracking_ -- precursorIdForTracking_ locks precursor_map_mutex_ internally
    //   logger_                   -- IdaLogger serialises each stream itself
    //   exploration_active_, current_faims_cv_ -- atomics, release/acquire
    //   faims_.isEnabled()        -- NOT synchronised. Safe only because FAIMS::enabled_ is assigned
    //                                once in the FAIMS constructor and never again. That immutability
    //                                is the drain's last unsynchronised cross-thread read and is
    //                                pinned nowhere; a "disable FAIMS mid-run" feature would make it
    //                                a real race.
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
      ScanCommandQueue::writeMs1Description(out, ScanRole::Agc);

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=agc" << std::endl;
      queue_.registerPending(out);
      // No analysis_mutex_. Taking it here would park the instrument event thread behind a whole
      // deconvolution to write one row -- see the note on analysis_mutex_ in FLASHIda.h. Serialising
      // a log stream is IdaLogger's job, not the analysis lock's.
      logger_.writeScanCommandRow(out);
      return 1;
    }

    const bool sweeping = exploration_active_.load(std::memory_order_acquire);

    // Step 2: Cycle time -- force MS1 if too long since last survey scan
    // Suppressed while any exploration group is active.
    // Queued at priority 0 (not returned immediately) so it goes through normal dequeue.
    if (config_.scheduling().cycle_time_enabled
        && !sweeping
        && queue_.msSinceLastMS1() > static_cast<uint64_t>(config_.scheduling().cycle_time_ms)
    )
    {
      ScanCommand ms1_cmd = queue_.makeMS1();
      ms1_cmd.faims_cv = faims_cv;
      ms1_cmd.scan_id = queue_.nextTrackingId();
      ms1_cmd.priority = 0;
      // (A `scan_description[3] = 'C'` used to sit here and was DEAD: the snprintf four lines below
      // overwrites the whole buffer, so a cycle-time MS1 has always gone out as 'S' -- correctly, it
      // is a survey. Removed with ADR-0038 rather than left to look like a third marker site.)
      queue_.recordMS1Time();

      std::string id_str = ScanCommandQueue::encode(ms1_cmd.scan_id);
      ScanCommandQueue::writeMs1Description(ms1_cmd, ScanRole::Survey);

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=cycle_time" << std::endl;
      queue_.push(ms1_cmd);
      // Fall through to Step 3 (cleanup) and Step 4 (dequeue)
    }

    // Step 2b: MONITOR scan -- a periodic MS1 that lets the operator watch the source while a sweep
    // is running, and that drives NO acquisition decision (ADR-0042).
    //
    // NOT an `else` of Step 2, though the two are behaviourally complementary (Step 2 fires only
    // when nothing is sweeping, Step 2b only when something is). takeMonitorDue must be reached on
    // EVERY drain, because its `sweeping == false` path is what re-arms a level for its next
    // episode. Put this in an else and the arm survives the end of a sweep, so the next sweep
    // silently loses its anchor scan -- visible only as a missing first data point.
    //
    // Both levels are always asked, even once one has answered yes, so each consumes its own clock:
    // it is ONE MS1 and it serves both levels' purposes, so neither should stay due afterwards.
    {
      int due_level = 0;
      for (int i = 0; i < 2; ++i)
      {
        const auto& lvl_cfg = config_.level(i + 2);
        const bool lvl_sweeping = lvl_cfg.monitor_ms1_enabled
                                  && exploration_active_lvl_[i].load(std::memory_order_acquire) > 0;
        if (queue_.takeMonitorDue(i, lvl_sweeping, lvl_cfg.monitor_ms1_interval_ms) && due_level == 0)
          due_level = i + 2;
      }

      if (due_level > 0)
      {
        // ms_settings.ms1 VERBATIM, so a monitor reading is directly comparable with the run's real
        // surveys -- that comparability is the point, and it is why there is no cheaper dedicated
        // scan config. At production settings (240k, 4 microscans) this costs ~2 s.
        ScanCommand mon_cmd = queue_.makeMS1();
        mon_cmd.faims_cv = faims_cv;   // the CV the sweep froze it at; a monitor scan never advances it
        mon_cmd.scan_id  = queue_.nextTrackingId();
        // Priority 0 is MANDATORY, not a preference. Priority 3 is ADR-0031's idle-survey sentinel
        // that five role-blind drain loops break on; 1 and 2 are the sweep's own lanes, so a monitor
        // scan would queue behind the very sweep it is meant to interrupt and the authored interval
        // would never govern. 0 is the lane the cycle-time and CV-transition MS1s already use.
        mon_cmd.priority = 0;
        // NO recordMS1Time() here -- see the drain tail. A monitor scan is not a survey.
        ScanCommandQueue::writeMs1Description(mon_cmd, ScanRole::Monitor);

        std::cout << "[TRACK-CREATE] id=" << ScanCommandQueue::encode(mon_cmd.scan_id)
                  << " ms_level=1 type=monitor lvl=" << due_level << std::endl;
        queue_.push(mon_cmd);   // pushed, not returned: inherits Step 4's registerPending, the
                                // timestamps and the scan_commands.tsv row
      }
    }

    // Step 3: Cleanup expired commands
    queue_.cleanupExpired();

    // Step 4: Dequeue by priority (0 = highest -> 3 = lowest)
    auto dequeued = queue_.dequeue();

    if (!dequeued.has_value())
    {
      // Step 5: Idle survey -- the queue is drained, so emit a fallback survey MS1 at the lowest
      // priority and hand it straight back, so the instrument is never left without a command.
      //
      // There is deliberately NO AGC prescan on this path (ADR-0031). A prescan is a scheduled
      // flux measurement, emitted by Step 1 when scheduling.agc_interval_seconds has elapsed and
      // at no other time; it is not filler for an idle instrument. Fabricating one here also
      // called recordAGCTime(), which reset that very timer -- so Step 1 could only ever fire in a
      // run whose queue stayed non-empty for a whole interval, and the authored interval never
      // governed the real cadence.
      ScanCommand ms1_cmd = queue_.makeMS1();
      ms1_cmd.faims_cv = faims_cv;
      ms1_cmd.scan_id = queue_.nextTrackingId();
      ms1_cmd.priority = 3;

      std::string ms1_id_str = ScanCommandQueue::encode(ms1_cmd.scan_id);
      ScanCommandQueue::writeMs1Description(ms1_cmd, ScanRole::Survey);

      std::cout << "[TRACK-CREATE] id=" << ms1_id_str << " ms_level=1 type=idle_ms1" << std::endl;

      dequeued = queue_.pushAndDequeue(ms1_cmd);

      // Re-entering dequeue() rather than returning ms1_cmd directly is what lets the idle survey
      // inherit the whole Step 4 tail below -- recordMS1Time(), the enqueue/dequeue timestamps,
      // pending_scan_map_ registration, and the log row carrying precursor_id. A bypass would get
      // none of those and would silently stop resetting the cycle-time clock.
      //
      // It is also the correct read under a concurrent processScan: if that thread pushed
      // higher-priority work BEFORE our atomic push-and-take, that command wins and ours waits at
      // priority 3. Taking under the SAME lock as the push is what makes this safe -- see
      // ScanCommandQueue::pushAndDequeue().
      //
      // The fallback below is now UNREACHABLE BY CONSTRUCTION: we pushed at a clamped priority
      // under the lock we then dequeue under, so dequeue_ always finds something, and
      // cancelByScanIds can no longer interleave either. It is retained only as a belt-and-braces
      // backstop for the hard invariant that getNextScanCommand must never return 0 for an empty
      // queue. Do NOT "restore" the old push()-then-dequeue() split to make it reachable again:
      // that gap let a concurrent drainer steal this survey, after which this thread logged its
      // local copy of a command the thief had ALREADY logged -- a duplicate tracking id, and one
      // fewer id allocated, because the thief never built a command of its own.
      if (!dequeued.has_value()) dequeued = ms1_cmd;
    }

    out = dequeued.value();
    // The added roleDecides conjunct keeps a MONITOR scan from resetting the SURVEY clock (ADR-0042).
    // It satisfies `msn_level == 1 && is_agc == 0` exactly as a survey does, so without this it would
    // stand in for one -- and that is character-for-character the ADR-0031 bug, where the idle path
    // called recordAGCTime() and the authored AGC interval silently stopped governing the cadence.
    // Here the victim would be scheduling.cycle_time: monitoring a long sweep would starve the very
    // survey the cycle time exists to force. Inert for every command that exists today, since Monitor
    // is the only non-AGC role with decides == false.
    if (out.msn_level == 1 && out.is_agc == 0 && roleDecides(roleOf(out)))
      queue_.recordMS1Time();
    // faims_cv already set at creation time (MS2 -> parent CV, CV-transition MS1 -> next CV)
    //
    // No analysis_mutex_ -- this is the one drain site that reads shared state, and it now takes
    // only the leaf locks it actually needs. precursorIdForTracking_ locks precursor_map_mutex_
    // itself; takeMS3Proteoform locks queue_mutex_ itself. Neither can be held for longer than a
    // hash lookup, whereas analysis_mutex_ is held across a whole deconvolution.
    //
    // The VALUE read here is already final and would be even without the map lock: the write
    // happens before its queue_.push(), and the dequeue() above acquired the same queue_mutex_
    // that push released. The lock exists for the container -- a concurrent insert can rehash the
    // table while find() is walking it. See FLASHIda.h.
    //
    // Take the wide MS3 fragment ProForma stashed at buildMS3 time (empty for non-MS3 commands) so the
    // scan_commands.tsv row carries the fired MS3 target's fragment identity.
    std::string ms3pf = (out.msn_level == 3) ? queue_.takeMS3Proteoform(out.scan_id) : std::string();
    logger_.writeScanCommandRow(out, precursorIdForTracking_(out.scan_id), ms3pf);
    return 1;
  }

  int FLASHIda::getNextTrackingId()
  {
    return queue_.nextTrackingId();
  }

} // namespace OpenMS
