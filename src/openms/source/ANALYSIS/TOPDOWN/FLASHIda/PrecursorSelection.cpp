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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/CandidateAdmission.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/NotchSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>

namespace OpenMS
{

  PrecursorSelection::PrecursorSelection(const Config& config, Deconvolution& deconv) :
      config_(config),
      deconv_(deconv)
  {
    // --- Load target log files ---
    const auto& targeting = config_.targeting();
    std::stringstream log_files_msg{};
    for (const auto& log_file : targeting.target_log_files)
    {
      log_files_msg << log_file << " ";
      std::ifstream instream(log_file);
      if (instream.good())
      {
        String line;
        double rt = .0, mass, qscore;
        int charge;
        while (std::getline(instream, line))
        {
          // " 0 targets" WITH THE LEADING SPACE — the same anchor as IdaLogger::parseFLASHIdaLog,
          // and for the same reason: the bare substring also matches "10 targets" / "20 targets".
          // Here the consequence is different and worse than a mis-keyed map. A dropped header leaves
          // `rt` at the PREVIOUS header's value (or 0.0 at the head of the file), so those masses are
          // filed in target_mass_rt_map_ against the wrong retention time — and mode 3 discards
          // rt == 0 entries outright.
          if (line.find(" 0 targets") != line.npos)
            continue;
          if (line.hasPrefix("MS1"))
          {
            Size st = line.find("RT ") + 3;
            Size ed = line.find('(') - 2;
            String n = line.substr(st, ed - st + 1);
            rt = atof(n.c_str());
          }
          if (line.hasPrefix("Mass"))
          {
            Size st = 5;
            Size ed = line.find('\t');
            String n = line.substr(st, ed - st + 1);
            mass = atof(n.c_str());

            st = line.find("Score=") + 6;
            ed = line.find('\t', st);
            n = line.substr(st, ed - st + 1);
            qscore = atof(n.c_str());

            st = line.find("Z=") + 2;
            ed = line.find('\t', st);
            n = line.substr(st, ed - st + 1);
            charge = n.toInt();

            if (targeting.mode == 1 || targeting.mode == 2)
            {
              if (target_mass_rt_map_.find(mass) == target_mass_rt_map_.end())
                target_mass_rt_map_[mass] = std::vector<double>();
              target_mass_rt_map_[mass].push_back(rt * 60.0);
              if (target_mass_qscore_map_.find(mass) == target_mass_qscore_map_.end())
                target_mass_qscore_map_[mass] = std::vector<double>();
              target_mass_qscore_map_[mass].push_back(qscore);
              if (target_mass_charge_map_.find(mass) == target_mass_charge_map_.end())
                target_mass_charge_map_[mass] = std::vector<int>();
              target_mass_charge_map_[mass].push_back(charge);
            }
          }
          else if (line.hasPrefix("AllMass"))
          {
            if (targeting.mode == 3)
            {
              Size st = 8;
              Size ed = line.size();
              String n = line.substr(st, ed - st + 1);
              std::stringstream tmp_stream(n);
              String str;
              std::vector<double> results;
              while (getline(tmp_stream, str, ' '))
                results.push_back(atof(str.c_str()));
              if (exclusion_rt_masses_map_.find(rt * 60.0) == exclusion_rt_masses_map_.end())
                exclusion_rt_masses_map_[rt * 60.0] = std::vector<double>();
              for (double m : results)
                exclusion_rt_masses_map_[rt * 60.0].push_back(m);
            }
          }
        }
        instream.close();
      }
    }

    if (targeting.mode == 1)
      std::cout << log_files_msg.str() << "file(s) is(are) used for inclusion mode\n";
    else if (targeting.mode == 2)
      std::cout << log_files_msg.str() << "file(s) is(are) used for in-depth mode\n";
    else if (targeting.mode == 3)
      std::cout << log_files_msg.str() << "file(s) is(are) used for exclusion mode\n";

    // Parse TSV inclusion list files
    if (targeting.mode == 1 && !targeting.inclusion_list_file.empty())
    {
      parseInclusionListTSV_(targeting.inclusion_list_file);
      std::cout << inclusion_targets_.size() << " targets loaded from TSV inclusion list\n";
    }

    // Load FASTA files for tag-based targeting
    if (!targeting.fasta_file.empty())
    {
      std::vector<FASTAFile::FASTAEntry> entries;
      FASTAFile().load(targeting.fasta_file, entries);
      target_protein_database_.insert(target_protein_database_.end(), entries.begin(), entries.end());
      std::cout << target_protein_database_.size() << " protein entries loaded for tag-based targeting\n";
    }

    if (!target_protein_database_.empty())
    {
      std::cout << "Tag-based targeting: min_tag_length=" << targeting.min_tag_length
                << ", max_tag_length=" << targeting.max_tag_length
                << ", tolerance=" << targeting.tag_matching_tolerance_ppm << " ppm"
                << ", max_flanking_mass_diff=" << targeting.max_flanking_mass_diff << " Da\n";
    }

    // Load PTM TSV files
    if (!targeting.ptm_list_file.empty())
      parseTargetPTMsTSV_(targeting.ptm_list_file);
    if (!target_ptms_.empty())
    {
      std::cout << target_ptms_.size() << " PTM modifications loaded for target expansion (max "
                << targeting.max_total_ptm_count << " total per proteoform)\n";
    }
  }

  int PrecursorSelection::filterAndRank(const double* mzs, const double* ints, int length,
                                        double rt, int ms_level, double faims_cv)
  {
    if (ms_level != 1)
    {
      return 0;
    }

    target_masses_.clear();
    excluded_masses_.clear();
    if (config_.targeting().mode == 1)  // Unified inclusion mode (merged mode 4 functionality)
    {
      active_targets_.clear();
      target_priority_map_.clear();

      if (!inclusion_targets_.empty())  // TSV-based targets
      {
        for (const auto& target : inclusion_targets_)
        {
          if (target.isActiveAt(rt))
          {
            active_targets_.push_back(&target);
            target_masses_.push_back(target.mass);

            // Build priority map for tie-breaking (use highest priority for each nominal mass)
            int nominal = SpectralDeconvolution::getNominalMass(target.mass);
            if (target_priority_map_.find(nominal) == target_priority_map_.end()
                || target.priority > target_priority_map_[nominal])
            {
              target_priority_map_[nominal] = target.priority;
            }
          }
        }
      }
      else  // Legacy .log/.out targets (existing behavior)
      {
        for (const auto& [mass, rts] : target_mass_rt_map_)
        {
          target_masses_.push_back(mass);
        }
      }

      std::sort(target_masses_.begin(), target_masses_.end());
      deconv_.engine().setTargetMasses(target_masses_, false);
    }
    else if (config_.targeting().mode == 3)
    {
      for (const auto& [prt, masses] : exclusion_rt_masses_map_)
      {
        if (std::abs(rt - prt) >= config_.targeting().rt_window && prt != 0) continue;
        for (double mass : masses)
        {
          excluded_masses_.push_back(mass);
        }
      }
      std::sort(excluded_masses_.begin(), excluded_masses_.end());
    }

    selected_peak_groups_.clear();

    // Deconvolve MS1 spectrum (result stored in deconv_.deconvolvedMS1())
    deconv_.deconvolveMS1(mzs, ints, length, rt, faims_cv);
    // per spec deconvolution
    filterPeakGroupsUsingMassExclusion_(ms_level, rt);
    return (int)selected_peak_groups_.size();
  }

  void PrecursorSelection::filterPeakGroupsUsingMassExclusion_(const int ms_level, const double rt)
  {
    if (config_.level(ms_level).selection == SelectionMetric::Intensity)
    {
      deconv_.deconvolvedMS1().sortByIntensity();
    }
    else
    {
      if (config_.targeting().consider_all_charges) {
        deconv_.deconvolvedMS1().sortByQScoreAllCharges();
      }
      else {
        deconv_.deconvolvedMS1().sortByQscore();
      }
    }

    // Apply priority tie-breaking when TSV targets are loaded
    if (config_.targeting().mode == 1 && !inclusion_targets_.empty())
    {
      std::stable_sort(deconv_.deconvolvedMS1().begin(), deconv_.deconvolvedMS1().end(),
        [this](const PeakGroup& a, const PeakGroup& b) {
          if (std::abs(a.getQscore() - b.getQscore()) < config_.targeting().tie_threshold)
          {
            int nom_a = SpectralDeconvolution::getNominalMass(a.getMonoMass());
            int nom_b = SpectralDeconvolution::getNominalMass(b.getMonoMass());
            int pri_a = target_priority_map_.count(nom_a) ? target_priority_map_.at(nom_a) : 0;
            int pri_b = target_priority_map_.count(nom_b) ? target_priority_map_.at(nom_b) : 0;
            return pri_a > pri_b;  // Higher priority first
          }
          return a.getQscore() > b.getQscore();
        });
    }

    Size mass_count = (Size)config_.level(ms_level).max_targets;
    trigger_charges_.clear();
    trigger_scores_.clear();
    trigger_charges_.reserve(mass_count);
    trigger_scores_.reserve(mass_count);
    trigger_left_isolation_mzs_.clear();
    trigger_left_isolation_mzs_.reserve(mass_count);
    trigger_right_isolation_mzs_.clear();
    trigger_right_isolation_mzs_.reserve(mass_count);
    trigger_ids_.clear();
    trigger_ids_.reserve(mass_count);
    trigger_authored_charges_.clear();
    trigger_authored_charges_.reserve(mass_count);
    std::vector<int>* charges = nullptr;

    selected_peak_groups_.reserve(mass_count);
    // Nominal masses that already spent one of their AUTHORED charges in THIS survey (ADR-0028).
    // Without it, the phase-1 pass of non-strict inclusion would come back for a second named
    // charge in the same survey and break `single`'s "one anchor per species per survey"
    // invariant (see the break at the bottom of the emit loop). Only bites under `single`:
    // separate/multiplexed spend the whole set in one pass anyway.
    std::set<int> authored_acquired_this_survey;

    std::set<double> current_selected_mzs;    // current selected mzs
    std::set<double> current_selected_masses; // current selected mzs

    std::unordered_map<int, double> new_mz_rt_map_;
    std::unordered_map<int, double> new_mass_rt_map_;
    std::unordered_map<int, double> new_all_mass_rt_map_;
    std::unordered_map<int, double> new_mass_score_map_;
    std::unordered_map<int, double> t_mass_score_map_;

    // IN-DEPTH mode (targeting == "in_depth"). NOT exclusion -- this comment said "exclusion mode"
    // for a long time while guarding mode 2, which is the opposite branch.
    //
    // Builds a de-prioritization product: t_mass_score_map_[nominal] *= (1 - qscore) over every
    // target-log observation in the RT window. A mass observed often, or observed well, drives the
    // product toward 0, and iteration 0 below then skips it while 1 - product exceeds
    // tqscore_threshold. It is SOFT: iteration 1 has no such guard and back-fills anything skipped,
    // so the effect is invisible unless the slot budget is contended.
    // TODO: Update IDScore bla, currently only qscore
    if (config_.targeting().mode == 2)
    {
      for (const auto& [mass, rts] : target_mass_rt_map_)
      {
        int nominal_mass = SpectralDeconvolution::getNominalMass(mass);
        auto qscores = target_mass_qscore_map_[mass];
        for (uint i = 0; i < rts.size(); i++)
        {
          double prt = rts[i];
          double qscore = qscores[i];
          if (std::abs(rt - prt) < config_.targeting().rt_window)
          {
            auto inter = t_mass_score_map_.find(nominal_mass);
            if (inter == t_mass_score_map_.end()) { t_mass_score_map_[nominal_mass] = 1 - qscore; }
            else { t_mass_score_map_[nominal_mass] *= 1 - qscore; }
          }
        }
      }
    }

    // remove expired entries for tqscore_exceeding_mz_rt_map_
    for (const auto& [m, r] : tqscore_exceeding_mz_rt_map_)
    {
      if (rt - r > config_.targeting().rt_window) { continue; }
      new_mz_rt_map_[m] = r;
    }
    new_mz_rt_map_.swap(tqscore_exceeding_mz_rt_map_);
    std::unordered_map<int, double>().swap(new_mz_rt_map_);

    // remove expired entries for tqscore_exceeding_mass_rt_map_
    for (const auto& [m, r] : tqscore_exceeding_mass_rt_map_)
    {
      if (rt - r > config_.targeting().rt_window) { continue; }
      new_mass_rt_map_[m] = r;
    }
    new_mass_rt_map_.swap(tqscore_exceeding_mass_rt_map_);
    std::unordered_map<int, double>().swap(new_mass_rt_map_);

    // remove expired entries for authored_acquired_rt_map_ (ADR-0028). Same keep-if-fresh shape:
    // a named charge becomes re-acquirable once its acquisition falls out of the RT window, which
    // is what the mass-keyed maps above have always done for masses.
    {
      std::map<std::pair<int, int>, double> fresh;
      for (const auto& [key, r] : authored_acquired_rt_map_)
      {
        if (rt - r > config_.targeting().rt_window) { continue; }
        fresh[key] = r;
      }
      authored_acquired_rt_map_.swap(fresh);
    }

    // remove expired entries for all_mass_rt_map_, mass_qscore_map_
    for (const auto& item : all_mass_rt_map_)
    {
      if (rt - item.second > config_.targeting().rt_window) { continue; }
      new_all_mass_rt_map_[item.first] = item.second;
      new_mass_score_map_[item.first] = mass_qscore_map_[item.first];
    }
    new_all_mass_rt_map_.swap(all_mass_rt_map_);
    std::unordered_map<int, double>().swap(new_all_mass_rt_map_);
    new_mass_score_map_.swap(mass_qscore_map_);
    std::unordered_map<int, double>().swap(new_mass_score_map_);

    const int selection_phase_start = 0;
    const int selection_phase_end = 2; // inclusive
    // When selection_phase == 0, consider only the masses whose tqscore did not exceed total qscore threshold.
    // when selection_phase == 1, consider all other masses for selection but the same m/z is avoided
    // when selection_phase == 2, consider all.
    // for target inclusive masses, qscore precursor snr threshold is not applied.
    // In all phase, for target exclusive mode, all the exclusive masses are excluded. For target inclusive mode, only the target masses are considered.

    for (int iteration = config_.targeting().mode == 2 ? 0 : 1; iteration < 2; iteration++)
    // for mass exclusion, first collect masses with exclusion list. Then collect without exclusion. This works the best
    {
      for (int selection_phase = selection_phase_start; selection_phase <= selection_phase_end; selection_phase++)
      {
        // Phase 0: targets (for inclusion mode) or tqscore-filtered masses
        // Phase 1: non-targets (only if non-strict inclusion mode and targets exist)
        if (selection_phase > 0)
        {
          // Allow phase 1 only for non-strict inclusion mode with active targets
          if (!(config_.targeting().mode == 1 && !config_.targeting().strict_inclusion && target_masses_.size() > 0 && selection_phase == 1))
          {
            break;
          }
        }

        // The budget counts SPECIES, not acquisitions. It matters only for
        // precursor_charges: "separate", which emits one scan per charge state of ONE species and must
        // cost one slot for the whole envelope -- exactly what "multiplexed" costs for the same envelope
        // in a single scan. Counting acquisitions instead would make max_precursors: 3 spend everything
        // on the first species the moment it had three charges. In every other mode a PeakGroup pushes at
        // most once, so this equals selected_peak_groups_.size() and nothing changes.
        Size species_selected = 0;

        // Iterate over candidates (sorted by qscore)
        for (const auto& pg : deconv_.deconvolvedMS1())
        {
          // dont acquire the same mass multiple times
          if (species_selected >= mass_count) { break; }
          const size_t pushed_before_pg = selected_peak_groups_.size();

          // The ANCHOR charge -- one per PeakGroup. This was a multi-element list only under
          // charge_based_exclusion, which is gone (ADR-0021); how many charges the scan then
          // ISOLATES is decided by precursor_charges further down, from the PeakGroup itself.
          //
          // Deliberately still a one-element loop rather than straight-line code: the body below
          // uses `break` and `continue` throughout, and both would silently retarget the enclosing
          // per-PeakGroup loop if this one were removed -- e.g. the qscore-threshold `break`, which
          // means "stop considering this species" and would become "stop selecting entirely".
          // The AUTHORED CHARGE SET in force for this species (ADR-0028). Resolved from the mass
          // alone, so it is available BEFORE the anchor is picked -- which it has to be, because
          // under an authored set the anchor is drawn from the set rather than from the envelope.
          const double mass_pg = pg.getMonoMass();
          const double delta_pg = 2 * config_.level(1).tolerance_ppm * mass_pg * 1e-6;
          const int nominal_pg = SpectralDeconvolution::getNominalMass(mass_pg);
          const std::vector<int> authored =
              config_.targeting().mode == 1 ? authoredChargesFor_(mass_pg, delta_pg) : std::vector<int>();
          const bool has_authored = !authored.empty();

          // The ANCHOR pick now lives in CandidateAdmission.h, unchanged in what it decides: drawn
          // from the authored SET when the row names charges (SNR orders the named charges, it does
          // not gate them), otherwise the envelope's own representative or best-qscore charge.
          //
          // The guard is `reason`, NOT `resolved()`, and the two differ on purpose. Only the authored
          // path ever refused a candidate for having no anchor; the unauthored path has always taken
          // the envelope's charge as-is, including a representative charge the deconvolution never
          // assigned. Guarding on `resolved()` would tighten that, which is a behaviour change and
          // does not belong in a lift-and-shift.
          const std::vector<int> spent_charges = spentAuthoredCharges_(nominal_pg);

          AnchorContext anchor_ctx;
          anchor_ctx.authored = authored;
          anchor_ctx.spent_charges = spent_charges;
          anchor_ctx.already_this_survey = authored_acquired_this_survey.count(nominal_pg) > 0;
          anchor_ctx.min_charge = config_.level(ms_level).min_charge;
          anchor_ctx.consider_all_charges = config_.targeting().consider_all_charges;

          const AnchorChoice anchor = pickAnchor(pg, anchor_ctx);
          if (anchor.reason != AdmissionReason::Admitted) { continue; }

          struct ChargeCandidate { int charge; double score; };
          std::vector<ChargeCandidate> charges_to_process {{anchor.charge, anchor.score}};

          for (const auto& cc : charges_to_process)
          {
            // No budget check for the additional charges of THIS species under "separate": the slot was
            // already paid for by the per-PeakGroup guard above, and the mode's contract is the whole
            // SNR-positive envelope. The set is bounded by the observed charge range and the SNR gate.
            if (config_.targeting().precursor_charges != ChargeAcquisitionMode::Separate
                && selected_peak_groups_.size() >= mass_count) { break; }
            int charge = cc.charge;
            double score = cc.score;

            // Per-level charge filter: ms1.min_charge controls what MS1 picks.
            //
            // Kept HERE, ahead of the geometry, rather than relying on admitCandidate's own copy of
            // the same rule. Moving it down would let a below-floor candidate reach the legacy
            // charge-list matching below, which assigns the `charges` pointer -- and that pointer
            // outlives this iteration, so a later candidate would erase from a list an earlier,
            // refused one had selected. admitCandidate re-checks the floor because it must be
            // decidable on its own; that copy is unreachable from here and is meant to be.
            if (config_.level(ms_level).min_charge > 0 && charge < config_.level(ms_level).min_charge)
              continue;

            double mass = pg.getMonoMass();

            auto [mz1, mz2] = pg.getMzRange(charge);

            double center_mz = (mz1 + mz2) / 2.0;

            mz1 -= optimal_window_margin_;
            mz2 += optimal_window_margin_;
            int integer_mz = (int)round(center_mz);

            int nominal_mass = SpectralDeconvolution::getNominalMass(mass);
            bool target_matched = false;
            double snr_threshold = config_.targeting().snr_threshold;
            double qscore_threshold = config_.targeting().qscore_threshold;
            double tqscore_factor_for_exclusion = 1.0;



            // Only triggered in exclusion mode
            if (iteration == 0)
            {
              auto inter = t_mass_score_map_.find(nominal_mass);
              if (inter != t_mass_score_map_.end()) { tqscore_factor_for_exclusion = t_mass_score_map_[nominal_mass]; }
              if (1 - tqscore_factor_for_exclusion > config_.targeting().tqscore_threshold) { continue; }
            }

            // Inclusion mode (unified: supports TSV-based targets and legacy .log/.out targets)
            if (config_.targeting().mode == 1 && target_masses_.size() > 0)
            {
              double delta = 2 * config_.level(1).tolerance_ppm * mass * 1e-6;
              auto ub = std::upper_bound(target_masses_.begin(), target_masses_.end(), mass + delta);

              while (!target_matched)
              {
                if (ub != target_masses_.end())
                {
                  if (std::abs(*ub - mass) < delta) // target is detected.
                  {
                    // Check charge matching for both TSV and legacy modes
                    if (!inclusion_targets_.empty())  // TSV mode
                    {
                      // A pure MASS question now. The mass matched at the binary search above, and
                      // any charge restriction the matching rows carry was already applied when the
                      // anchor was picked (authoredChargesFor_ + the authored anchor loop), so
                      // there is nothing left to reassign here.
                      //
                      // This replaces a loop over active_targets_ that took the first row whose
                      // charge fell inside the envelope WITHOUT checking that row's mass -- so with
                      // several rows live at one retention time the isolated charge could be
                      // supplied by an unrelated target, and which row won depended on std::sort's
                      // unspecified order among equal masses. Invisible in CI only because every
                      // committed list writes -1, which hit the `charge < 0` arm on the first row.
                      target_matched = true;
                    }
                    else if (!target_mass_charge_map_.empty())  // Legacy .log mode with charge
                    {
                      auto it = target_mass_charge_map_.find(*ub);
                      if (it != target_mass_charge_map_.end())
                      {
                        charges = &it->second;
                        if (std::find(charges->begin(), charges->end(), charge) != charges->end())
                        {
                          target_matched = true;
                        }
                      }
                    }
                    else  // Legacy mode without charge (mass-only matching)
                    {
                      target_matched = true;
                    }

                    if (!target_matched)
                    {
                      if (ub == target_masses_.begin()) { break; }
                      ub--;
                      continue;
                    }
                  }
                  if (mass - *ub > delta) { break; }
                }
                if (ub == target_masses_.begin()) { break; }
                ub--;
              }

              if (target_matched)
              {
                snr_threshold = 0.0;
                qscore_threshold = 0.0; // stop exclusion for targets. todo tqscore lowest first? charge change.
              }
              else if (selection_phase == 0)
              {
                // Phase 0: only targets (strict behavior)
                continue;
              }
              // Phase 1: non-targets proceed with default thresholds
            }
            else if (config_.targeting().mode == 1 && config_.targeting().strict_inclusion)
            {
              // Strict inclusion mode with no active targets - skip all candidates
              continue;
            }
            // EXCLUSION mode (targeting == "exclusion_masses"). NOT deep -- this comment said
            // "deep mode" while guarding mode 3. Unlike the soft de-prioritization at mode 2, this
            // is a HARD skip: a matched mass is never selected, whatever the slot budget.
            else if (config_.targeting().mode == 3 && excluded_masses_.size() > 0)
            {
              bool to_exclude = false;
              double delta = 2 * config_.level(1).tolerance_ppm * mass * 1e-6;
              auto ub = std::upper_bound(excluded_masses_.begin(), excluded_masses_.end(), mass + delta);

              while (! to_exclude)
              {
                if (ub != excluded_masses_.end())
                {
                  if (std::abs(*ub - mass) < delta) // target is detected.
                  {
                    to_exclude = true;
                  }
                  if (mass - *ub > delta) { break; }
                }
                if (ub == excluded_masses_.begin()) { break; }
                ub--;
              }

              if (to_exclude) { continue; }
            }

            if (current_selected_mzs.find(center_mz) != current_selected_mzs.end()) // mz has been triggered
            {
              if (selection_phase < selection_phase_end) { continue; }
              if (! target_matched && current_selected_masses.find(pg.getMonoMass()) == current_selected_masses.end()) // but mass is different
              {
                continue;
              }
            }

            // The per-candidate gates and dynamic exclusion now live in CandidateAdmission.h. Same
            // decisions and the same order among themselves: score, then the anchor's OWN SNR, then
            // exclusion -- per-charge for an authored species, ORed over the nominal-mass and
            // integer-m/z bars for everyone else. Exclusion applies only in the selection phases
            // that apply it at all, so the three lookups are gated here rather than inside.
            //
            // Two deliberate non-equivalences, both invisible in behaviour:
            //
            // - These sit AFTER the within-survey m/z dedup above, where the score and SNR gates used
            //   to sit before it. Both are side-effect-free skips of the same candidate, so which one
            //   reports first is unobservable; leaving the dedup textually in place keeps the diff
            //   legible.
            // - The score gate was a `break` out of charges_to_process and is now a refusal.
            //   That loop has exactly one element -- one anchor per PeakGroup (ADR-0021) -- so
            //   `break` and `continue` reach the same next statement, the species count below, which
            //   does not fire because nothing was pushed.
            //
            // spent_charges and anchor_spent answer DIFFERENT questions and neither derives from the
            // other: the notch filter reads the per-charge record unconditionally, the exclusion bar
            // only in the exclusion-applying phases.
            double bar_value = 0.0;
            bool   has_bar = false;
            {
              auto bar_it = mass_qscore_map_.find(nominal_mass);
              if (bar_it != mass_qscore_map_.end()) { bar_value = bar_it->second; has_bar = true; }
            }

            AdmissionContext adm_ctx;
            adm_ctx.target_matched = target_matched;
            adm_ctx.has_authored = has_authored;
            adm_ctx.authored = authored;
            adm_ctx.spent_charges = spent_charges;
            adm_ctx.anchor_spent = selection_phase < selection_phase_end
                                   && authored_acquired_rt_map_.count({nominal_mass, charge}) > 0;
            adm_ctx.mass_barred = selection_phase < selection_phase_end
                                  && tqscore_exceeding_mass_rt_map_.find(nominal_mass) != tqscore_exceeding_mass_rt_map_.end();
            adm_ctx.mz_barred = selection_phase < selection_phase_end
                                && tqscore_exceeding_mz_rt_map_.find(integer_mz) != tqscore_exceeding_mz_rt_map_.end();
            adm_ctx.qscore_bar = has_bar ? &bar_value : nullptr;
            adm_ctx.fan_out = config_.targeting().precursor_charges != ChargeAcquisitionMode::Single;
            adm_ctx.min_charge = config_.level(ms_level).min_charge;
            adm_ctx.snr_threshold = config_.targeting().snr_threshold;
            adm_ctx.anchor_snr_threshold = snr_threshold;
            adm_ctx.qscore_threshold = qscore_threshold;
            adm_ctx.tqscore_threshold = config_.targeting().tqscore_threshold;
            adm_ctx.max_notches = MAX_NOTCHES_PER_STAGE;
            adm_ctx.where = "MS2 z=" + std::to_string(charge);

            const AdmissionVerdict verdict = admitCandidate(pg, charge, score, adm_ctx);

            // Acquisition memory is stamped once exclusion is cleared and BEFORE the qscore ledger
            // can refuse, so a candidate the ledger turns away still refreshes the record that keeps
            // its stored score alive through expiry. That ordering is the whole reason the verdict
            // reports passing exclusion separately from being admitted.
            if (verdict.passed_exclusion) { all_mass_rt_map_[nominal_mass] = rt; }
            if (! verdict.admit) { continue; }

            // The species' ACQUISITION CHARGE SET -- every charge this survey will fragment, not just
            // the anchor. Both non-single modes acquire the SAME set and differ only in scan count
            // (ADR-0016): multiplexed co-isolates it as notches in one scan, separate emits one scan
            // per member. Deriving both from this one call is what keeps that true; it is also how the
            // MS3 side does it (ProteoformTracker::planNextScans).
            //
            // Recording only the anchor would leave its siblings eligible and the next survey's
            // fallback would land on a charge already fragmented (ADR-0018). The set comes from the
            // same peakGroupNotchCandidates + selectNotches pair buildMS2 uses, so what is recorded as
            // acquired is by construction what gets isolated. Computed once: the qscore accumulation,
            // the RT map and the emit loop below must all record the same set.
            //
            // The anchor is `charge` as pickAnchor resolved it -- under an authored charge set that
            // is the highest-SNR NAMED charge, not the envelope's representative one. selectNotches
            // drops the anchor from its output, so a stale anchor would emit it twice.
            //
            // Computed inside admitCandidate, from the same peakGroupNotchCandidates + selectNotches
            // pair buildMS2 uses. Binding a reference here rather than copying keeps the single
            // source visible: the qscore accumulation, the per-charge record and the emit loop below
            // must all record the same set, and now they cannot each derive their own.
            const std::vector<int>& acquired_charges = verdict.acquisition_charges;
            std::vector<int> authored_dropped;  // named, observed, but refused -- reported below

            if (has_authored)
            {
              // Per-charge exclusion (ADR-0028): record what was actually isolated, charge by
              // charge, and write NO mass-keyed entry -- see the mass_qscore_map_ block below.
              for (int c : acquired_charges) { authored_acquired_rt_map_[{nominal_mass, c}] = rt; }
              authored_acquired_this_survey.insert(nominal_mass);

              // Say what was named and refused. A silent drop reads as "we isolated everything you
              // asked for" when we did not -- the same failure shape [NOTCH-CLAMP] exists for.
              std::stringstream charge_set_msg;
              // rt identifies the survey. Without it, lines from different surveys -- and, in a
              // ctest log, from different engine runs -- are indistinguishable, which is most of
              // what made the CM-04 scoping failure slow to attribute.
              charge_set_msg << "[CHARGE-SET] rt=" << rt << " m=" << mass << " authored=";
              for (size_t ai = 0; ai < authored.size(); ++ai)
                charge_set_msg << (ai ? ";" : "") << authored[ai];
              charge_set_msg << " anchor=" << charge << " acquired=";
              for (size_t ai = 0; ai < acquired_charges.size(); ++ai)
                charge_set_msg << (ai ? "," : "") << acquired_charges[ai];
              // "unacquired", not "dropped": under `single` the other named charges are DEFERRED to
              // later surveys, not refused, and per-charge exclusion is what will deliver them.
              // Calling that a drop would read as a failure every single survey.
              for (int c : authored)
              {
                if (std::find(acquired_charges.begin(), acquired_charges.end(), c) != acquired_charges.end())
                  continue;
                auto [lo, hi] = pg.getMzRange(c);
                authored_dropped.push_back(c);
                charge_set_msg << (authored_dropped.size() == 1 ? " unacquired=" : ",") << c << "(";
                if (hi <= lo) charge_set_msg << "not resolved";
                else if (config_.level(ms_level).min_charge > 0 && c < config_.level(ms_level).min_charge)
                  charge_set_msg << "below min_charge " << config_.level(ms_level).min_charge;
                // Ahead of the SNR arm deliberately: a spent charge has a perfectly good SNR and
                // would otherwise be reported as having failed a gate it never reached.
                else if (authored_acquired_rt_map_.count({nominal_mass, c}) > 0)
                  charge_set_msg << "already acquired this rt window";
                else if (config_.targeting().precursor_charges == ChargeAcquisitionMode::Single)
                  charge_set_msg << "deferred, precursor_charges single";
                else charge_set_msg << "snr " << pg.getChargeSNR(c) << " < " << config_.targeting().snr_threshold;
                charge_set_msg << ")";
              }
              std::cout << charge_set_msg.str() << std::endl;
            }

            // Compute total qscore. Mass-keyed. Runs ONCE per species even when the scan isolates
            // several of its charges -- see the emit loop below for why that matters.
            //
            // Skipped entirely for an authored charge set (ADR-0028), and both halves have to go.
            // The tqscore_exceeding_* writes would retire the MASS after its first charge, leaving
            // the other named charges unreachable; and the "previously acquired with higher
            // qscore" guard is mass-keyed too, so on the next survey the second charge -- ranked
            // lower by construction -- would `continue` against the first charge's score. Per-charge
            // state was already recorded above.
            // The ledger DECISION -- whether this survey beats the score the species was acquired at
            // -- was made inside admitCandidate and has already turned this candidate away if it did
            // not. What is left here is the WRITE, which belongs to the class that owns the maps.
            if (! has_authored)
            {
              mass_qscore_map_[nominal_mass] = verdict.record_score;
              if (verdict.arms_bars)
              {
                tqscore_exceeding_mass_rt_map_[nominal_mass] = rt;
                tqscore_exceeding_mz_rt_map_[integer_mz] = rt;
              }
            }

            // For legacy .log mode with charge targeting, remove this charge from list
            if (config_.targeting().mode == 1 && charges != nullptr) {
              auto it = std::find(charges->begin(), charges->end(), charge);
              if (it != charges->end()) {
                charges->erase(it);
              }
            }

            // Store acquisition. "separate" emits the acquisition charge set as one record PER charge
            // -- each its own Precursor with its own precursor_id (Config.h:82) -- while single and
            // multiplexed emit exactly one, byte-identically to before.
            //
            // The mass-level bookkeeping above stays deliberately OUTSIDE this loop. It ran once for
            // the species and must not run again per charge: the tqscore_exceeding_mass_rt_map_ guard
            // and the "previously acquired with higher qscore" guard in the mass_qscore_map_ block are
            // both keyed on nominal_mass, so a second pass would `continue` on every sibling -- the
            // first writes the mass on the anchor, and charges are ranked descending so every later
            // one scores lower. That is exactly why the fan-out could not simply be a walk of
            // charges_to_process, and why "separate" was inert whenever charge_based_exclusion --
            // the only thing that made that list multi-charge -- was off, i.e. in every shipped config.
            const bool emit_per_charge =
                config_.targeting().precursor_charges == ChargeAcquisitionMode::Separate;

            for (size_t ei = 0, emit_count = emit_per_charge ? acquired_charges.size() : 1;
                 ei < emit_count; ++ei)
            {
              const int emit_charge = emit_per_charge ? acquired_charges[ei] : charge;

              // The anchor keeps the geometry computed above: inclusion-mode target matching may have
              // recomputed mz1/mz2/center_mz for a matched charge (:486-489), and recomputing here
              // would silently discard that. Siblings measure their own window the same way, margin
              // included, so every emitted record carries MEASURED geometry.
              double e_mz1 = mz1, e_mz2 = mz2, e_center_mz = center_mz;
              if (emit_charge != charge)
              {
                std::tie(e_mz1, e_mz2) = pg.getMzRange(emit_charge);
                e_center_mz = (e_mz1 + e_mz2) / 2.0;
                e_mz1 -= optimal_window_margin_;
                e_mz2 += optimal_window_margin_;
              }

              id_mass_map_[window_id_] = nominal_mass;
              id_mz_map_[window_id_] = (int)round(e_center_mz);
              id_qscore_map_[window_id_] = score;
              id_charge_map_[window_id_] = emit_charge;
              // The set the geometry writer must not exceed. Under an authored charge set this is
              // exactly what was recorded above; otherwise empty, meaning unrestricted, which
              // reproduces the envelope-wide behaviour every non-targeted species has always had.
              trigger_authored_charges_.push_back(has_authored ? acquired_charges : std::vector<int>());

              trigger_ids_.push_back(window_id_);
              window_id_++;

              selected_peak_groups_.push_back(pg);
              trigger_charges_.push_back(emit_charge);
              trigger_scores_.push_back(score);

              trigger_left_isolation_mzs_.push_back(e_mz1);
              trigger_right_isolation_mzs_.push_back(e_mz2);
              current_selected_masses.insert(pg.getMonoMass());
              current_selected_mzs.insert(e_center_mz);
            }

            // ONE anchor per PeakGroup per survey. How many charges the survey actually FRAGMENTS is
            // decided by precursor_charges at the emit loop above, from the PeakGroup's own envelope
            // -- never here.
            //
            // This distinction is the whole history of this loop. It once walked a multi-charge list
            // built by the charge_based_exclusion developer flag, which meant one proteoform could
            // consume as many max_precursors slots as it had charges (with max_precursors 3 and three
            // charges, P2 and P3 were never fragmented). Adding this break fixed that, but the
            // "separate" mode was then implemented as a suppressed break over the same list -- so
            // acquisition GEOMETRY silently depended on an exclusion-KEYING flag, and "separate"
            // equalled "single" wherever that flag sat at its default, i.e. everywhere.
            //
            // The flag is now gone (ADR-0021) and the break is unconditional. Asking for several
            // charge states is precursor_charges: "separate" or "multiplexed", and nothing else.
            break;
          }  // end for charges_to_process

          // One species consumed, however many acquisitions it produced.
          if (selected_peak_groups_.size() > pushed_before_pg) ++species_selected;
        }
      }
    }
  }

  void PrecursorSelection::removeFromExclusionList(int id)
  {
    // Check if id is valid
    if (id >= window_id_) { return; }

    // Obtain information needed for removal
    int nominal_mass = id_mass_map_[id];
    int integer_mz = id_mz_map_[id];
    double qscore = id_qscore_map_[id];

    // Remove from mass exclusion
    if (tqscore_exceeding_mass_rt_map_.find(nominal_mass) != tqscore_exceeding_mass_rt_map_.end())
    {
      tqscore_exceeding_mass_rt_map_.erase(nominal_mass);
    }

    // Remove from mz exclusion
    if (tqscore_exceeding_mz_rt_map_.find(integer_mz) != tqscore_exceeding_mz_rt_map_.end()) { tqscore_exceeding_mz_rt_map_.erase(integer_mz); }

    // Remove the per-charge entry an authored charge set would have written (ADR-0028). Keyed on
    // the charge this very window isolated, so the named charge becomes selectable again rather
    // than the whole set.
    authored_acquired_rt_map_.erase({nominal_mass, id_charge_map_[id]});

    // Remove qscore from further calculations
    if (mass_qscore_map_.find(nominal_mass) != mass_qscore_map_.end()) { mass_qscore_map_[nominal_mass] /= 1 - qscore; }
  }

  std::vector<int> PrecursorSelection::authoredChargesFor_(double mass, double delta) const
  {
    std::vector<int> out;
    // active_targets_ is already RT-filtered (filterAndRank, :190-207), so "active" needs no
    // second RT test here -- which is exactly what keeps two rows for one mass in DIFFERENT RT
    // windows independent while unioning two rows in the SAME window.
    for (const InclusionTarget* t : active_targets_)
    {
      if (std::abs(t->mass - mass) >= delta) continue;
      // An unrestricted row wins the union outright. Intersecting instead would let a bare `-1`
      // row silently narrow a neighbouring row's set, and "no charge named" has to mean "no
      // opinion", not "every charge".
      if (t->charges.empty()) return {};
      out.insert(out.end(), t->charges.begin(), t->charges.end());
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
  }

  std::vector<int> PrecursorSelection::spentAuthoredCharges_(int nominal_mass) const
  {
    // authored_acquired_rt_map_ is keyed on the pair, and std::map orders pairs lexicographically,
    // so one nominal mass's charges are a contiguous range. Walking it beats querying per charge:
    // the anchor pick, the notch filter and the exclusion bar all want the same set, and gathering
    // it once is what keeps them consulting identical state.
    std::vector<int> out;
    for (auto it = authored_acquired_rt_map_.lower_bound({nominal_mass, std::numeric_limits<int>::min()});
         it != authored_acquired_rt_map_.end() && it->first.first == nominal_mass; ++it)
    {
      out.push_back(it->first.second);
    }
    return out;
  }

  void PrecursorSelection::parseInclusionListTSV_(const String& filename)
  {
    std::ifstream instream(filename);
    if (!instream.good())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    String line;
    bool header_skipped = false;

    while (std::getline(instream, line))
    {
      line = line.trim();
      if (line.empty() || line.hasPrefix("#")) continue;

      StringList cells;
      line.split('\t', cells);

      // Skip header row (check if first cell looks like "mass" header)
      if (!header_skipped && cells.size() > 0 && cells[0].toLower() == "mass")
      {
        header_skipped = true;
        continue;
      }

      // Expect 5 columns: mass, charge, rt_start, rt_end, priority
      if (cells.size() < 5)
      {
        OPENMS_LOG_WARN << "Skipping malformed line (expected 5 columns): " << line << std::endl;
        continue;
      }

      InclusionTarget target;
      target.mass = cells[0].toDouble();

      // The AUTHORED CHARGE SET (ADR-0028). "-1" or empty leaves the row unrestricted and yields
      // an EMPTY set; a single charge is just a set of one; "10;13;16" names three.
      //
      // ';' rather than ',' deliberately: this file is tab-separated so either would parse, but a
      // spreadsheet in a comma-decimal locale writes 12351,3 and a comma list would then be
      // indistinguishable from a mangled number. The wire grammar's ',' axis is a different thing
      // in a different file and does not reach here.
      String charge_str = cells[1].trim();
      if (!charge_str.empty() && charge_str != "-1")
      {
        // String::split pushes the whole string and returns false when the delimiter is absent, so
        // a lone "12" arrives here as a one-element list and needs no special case. It does NOT
        // trim on that path, hence the trim below.
        StringList charge_parts;
        charge_str.split(';', charge_parts);
        for (const String& part : charge_parts)
        {
          const String tok = String(part).trim();
          if (tok.empty()) continue;
          const int z = tok.toInt();
          // Loud rather than silent, per the strict-schema stance (ADR-0007): a charge that cannot
          // be isolated is a typo, and a row silently reduced to "unrestricted" would acquire the
          // whole envelope while reading as though it named three charges.
          if (z <= 0)
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "inclusion list: charge must be a positive integer, a ';'-separated list of them, or -1",
              tok);
          }
          target.charges.push_back(z);
        }
        std::sort(target.charges.begin(), target.charges.end());
        target.charges.erase(std::unique(target.charges.begin(), target.charges.end()),
                             target.charges.end());
      }
      // RT values in minutes, convert to seconds
      target.rt_start = cells[2].toDouble() * 60.0;
      target.rt_end = cells[3].toDouble() * 60.0;
      target.priority = cells[4].toInt();

      inclusion_targets_.push_back(target);
    }

    // Sort by mass for efficient binary search later
    std::sort(inclusion_targets_.begin(), inclusion_targets_.end(),
      [](const InclusionTarget& a, const InclusionTarget& b) { return a.mass < b.mass; });
  }

  void PrecursorSelection::parseTargetPTMsTSV_(const String& filename)
  {
    std::ifstream instream(filename);
    if (!instream.good())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    String line;
    bool header_skipped = false;

    while (std::getline(instream, line))
    {
      line = line.trim();
      if (line.empty() || line.hasPrefix("#")) continue;

      StringList cells;
      line.split('\t', cells);

      // Skip header row (check if first cell looks like "name" header)
      if (!header_skipped && cells.size() > 0 && cells[0].toLower() == "name")
      {
        header_skipped = true;
        continue;
      }

      // Expect 3 columns: name, mass, max_count
      if (cells.size() < 3)
      {
        OPENMS_LOG_WARN << "Skipping malformed PTM line (expected 3 columns): " << line << std::endl;
        continue;
      }

      TargetPTM ptm;
      ptm.name = cells[0].trim();
      ptm.mass = cells[1].toDouble();
      ptm.max_count = cells[2].toInt();

      target_ptms_.push_back(ptm);
    }
  }

  std::vector<double> PrecursorSelection::generatePTMCombinations_(double base_mass,
                                                                   const std::vector<TargetPTM>& ptms) const
  {
    std::vector<double> result;
    // Unmodified base mass intentionally excluded - only PTM-modified forms

    if (ptms.empty())
    {
      return result;
    }

    // Use iterative approach to generate all combinations
    // For each PTM, generate 0 to max_count occurrences, respecting global max
    std::function<void(Size, double, int)> generate;
    generate = [&](Size ptm_idx, double current_mass, int total_count) {
      if (ptm_idx >= ptms.size())
      {
        if (total_count > 0)  // Only include modified forms
        {
          result.push_back(current_mass);
        }
        return;
      }

      const TargetPTM& ptm = ptms[ptm_idx];
      // Generate both negative and positive mass shifts: -max_count to +max_count
      for (int count = -ptm.max_count; count <= ptm.max_count; ++count)
      {
        int new_total = total_count + std::abs(count);  // Track absolute count for pruning
        if (new_total > config_.targeting().max_total_ptm_count) continue;  // Skip if exceeds global max
        double new_mass = current_mass + count * ptm.mass;
        generate(ptm_idx + 1, new_mass, new_total);
      }
    };

    // Generate all combinations starting from base mass
    generate(0, base_mass, 0);

    // Remove duplicates (masses within 0.01 Da are considered equal)
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end(),
      [](double a, double b) { return std::abs(a - b) < 0.01; }), result.end());

    return result;
  }

  void PrecursorSelection::addDynamicTargets_(const std::vector<double>& masses,
                                              double rt,
                                              int priority)
  {
    // Add new targets to inclusion list
    // Note: target_masses_ and deconv_.engine().setTargetMasses() are managed by filterAndRank()
    // which clears and rebuilds them from inclusion_targets_ on each MS1 scan.
    // We only need to add to inclusion_targets_ here.
    for (double mass : masses)
    {
      InclusionTarget target;
      target.mass = mass;
      // No authored charge set: a tag-expanded target is a MASS hypothesis, and nothing in the tag
      // match says anything about which charge states of it will be worth isolating. Leaving
      // `charges` empty keeps acquisition unrestricted, exactly as the old `-1` did.
      target.rt_start = rt;
      target.rt_end = rt + config_.targeting().rt_window;
      target.priority = priority;

      inclusion_targets_.push_back(target);
    }

    // Re-sort inclusion_targets_ by mass for efficient lookup
    std::sort(inclusion_targets_.begin(), inclusion_targets_.end(),
      [](const InclusionTarget& a, const InclusionTarget& b) { return a.mass < b.mass; });

    std::cout << "Added " << masses.size() << " dynamic target masses (RT window: "
              << rt << "-" << (rt + config_.targeting().rt_window) << "s)\n";
    std::cout << "  Masses: ";
    for (Size i = 0; i < masses.size(); ++i)
    {
      if (i > 0) std::cout << ", ";
      std::cout << std::fixed << std::setprecision(4) << masses[i];
    }
    std::cout << std::endl;
  }

  int PrecursorSelection::processMS2ForTagBasedTargeting(double precursor_mass, const std::string& activation_type)
  {
    // Early exit if tag-based targeting not enabled
    if (target_protein_database_.empty())
    {
      return false;
    }

    // Require deconvolveMSn() to be called first
    if (!deconv_.hasStoredMS2())
    {
      return false;
    }

    // Use stored MS2 deconvolution
    DeconvolvedSpectrum dspec = deconv_.storedMS2();
    if (dspec.empty())
    {
      return false;
    }

    // Sort deconvolved spectrum by mass
    dspec.sort();

    // Create and configure tagger with tag length parameters
    FLASHTaggerAlgorithm tagger;
    tagger.setParameters(
      FragmentAnalysis::buildTaggerParam(config_, FragmentAnalysis::getIonTypesForFragmentationMethod(activation_type)));

    // Run tag generation
    tagger.run(dspec, config_.targeting().tag_matching_tolerance_ppm);

    // Get the generated tags
    std::vector<FLASHHelperClasses::Tag> tags;
    tagger.fillTags(tags);

    if (tags.empty())
    {
      return false;
    }

    // Prepare protein sequences for matching (replace I with L for matching)
    std::vector<String> protein_seqs;
    protein_seqs.reserve(target_protein_database_.size());
    for (const auto& fe : target_protein_database_)
    {
      String seq = fe.sequence;
      std::replace(seq.begin(), seq.end(), 'I', 'L');
      protein_seqs.push_back(seq);
    }

    // Match tags against target protein database
    std::set<Size> matched_protein_indices;
    for (const auto& tag : tags)
    {
      std::cout << tag.getSequence() << std::endl;
      for (Size protein_idx = 0; protein_idx < protein_seqs.size(); ++protein_idx)
      {
        const String& pseq = protein_seqs[protein_idx];

        std::vector<int> positions;
        std::vector<double> flanking_mass_diffs;

        FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(
          positions,
          flanking_mass_diffs,
          config_.targeting().max_flanking_mass_diff,
          pseq,
          tag);

        if (!positions.empty())
        {
          matched_protein_indices.insert(protein_idx);
        }
      }
    }

    if (matched_protein_indices.empty())
    {
      return false;
    }

    // Target protein detected! Expand target masses with PTM combinations
    // Use precursor mass from iAPI as base (not theoretical protein mass)
    std::vector<double> new_targets;
    if (precursor_mass > 0)
    {
      std::vector<double> ptm_masses = generatePTMCombinations_(precursor_mass, target_ptms_);
      new_targets.insert(new_targets.end(), ptm_masses.begin(), ptm_masses.end());
    }

    // Remove duplicates and already-tracked masses
    std::sort(new_targets.begin(), new_targets.end());
    new_targets.erase(std::unique(new_targets.begin(), new_targets.end(),
      [](double a, double b) { return std::abs(a - b) < 0.01; }), new_targets.end());

    // Filter out masses already in expanded set
    std::vector<double> truly_new_targets;
    for (double mass : new_targets)
    {
      int nominal = SpectralDeconvolution::getNominalMass(mass);
      bool already_tracked = false;
      for (double tracked : expanded_target_masses_)
      {
        if (std::abs(SpectralDeconvolution::getNominalMass(tracked) - nominal) <= 1)
        {
          already_tracked = true;
          break;
        }
      }
      if (!already_tracked)
      {
        truly_new_targets.push_back(mass);
        expanded_target_masses_.insert(mass);
      }
    }

    if (truly_new_targets.empty())
    {
      return static_cast<int>(tags.size()); // target detected (all masses already tracked); tag_count = #tags
    }

    // Add to dynamic inclusion list
    addDynamicTargets_(truly_new_targets, deconv_.storedMS2RT(), 100); // High priority

    return static_cast<int>(tags.size());
  }

} // namespace OpenMS
