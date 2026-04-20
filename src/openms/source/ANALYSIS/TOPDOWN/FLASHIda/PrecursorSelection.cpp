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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>

namespace OpenMS
{

  /// optimal window margin
  inline const double optimal_window_margin_ = .4;

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
          if (line.find("0 targets") != line.npos)
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
    // IDScore replaces QScore but not intensity
    if (config_.targeting().use_idscore)
    {
      if (config_.targeting().consider_all_charges && config_.targeting().hcd_energy < 0) {
        deconv_.deconvolvedMS1().sortByIDScoreAllCharges();
      }
      else if (config_.targeting().consider_all_charges) {
        deconv_.deconvolvedMS1().sortByIDScoreAllCharges(config_.targeting().hcd_energy);
      }
      else if (config_.targeting().hcd_energy < 0) {
        deconv_.deconvolvedMS1().sortByIDScoreRepresentative();
      }
      else {
        deconv_.deconvolvedMS1().sortByIDScoreRepresentative(config_.targeting().hcd_energy);
      }
    }
    else if (config_.level(ms_level).selection == SelectionMetric::Intensity)
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
    trigger_hcds_.clear();
    trigger_scores_.clear();
    trigger_charges_.reserve(mass_count);
    trigger_hcds_.reserve(mass_count);
    trigger_scores_.reserve(mass_count);
    trigger_left_isolation_mzs_.clear();
    trigger_left_isolation_mzs_.reserve(mass_count);
    trigger_right_isolation_mzs_.clear();
    trigger_right_isolation_mzs_.reserve(mass_count);
    trigger_ids_.clear();
    trigger_ids_.reserve(mass_count);
    std::vector<int>* charges = nullptr;

    selected_peak_groups_.reserve(mass_count);
    std::set<double> current_selected_mzs;    // current selected mzs
    std::set<double> current_selected_masses; // current selected mzs

    std::unordered_map<int, double> new_mz_rt_map_;
    std::unordered_map<int, double> new_mass_rt_map_;
    std::unordered_map<int, double> new_all_mass_rt_map_;
    std::unordered_map<int, double> new_mass_score_map_;
    std::unordered_map<int, double> t_mass_score_map_;

    // exclusion mode
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

        // Iterate over candidates (sorted by qscore)
        for (const auto& pg : deconv_.deconvolvedMS1())
        {
          // dont acquire the same mass multiple times
          if (selected_peak_groups_.size() >= mass_count) { break; }

          struct ChargeCandidate { int charge; double score; int hcd; };
          std::vector<ChargeCandidate> charges_to_process;

          if (config_.targeting().charge_based_exclusion)
          {
            auto [min_c, max_c] = pg.getAbsChargeRange();
            const auto& all_qs = pg.getAllQscores();
            for (int c = min_c; c <= max_c; ++c)
            {
              if (all_qs.count(c) == 0) { continue; }
              double s;
              int h = config_.targeting().hcd_energy;
              if (config_.targeting().use_idscore && config_.targeting().hcd_energy < 0)
              {
                s = pg.getBestIDScoreForCharge(c);
                h = pg.getBestHCDForCharge(c);
              }
              else if (config_.targeting().use_idscore)
              {
                s = pg.getIDScoreForChargeAndHCD(c, config_.targeting().hcd_energy);
              }
              else
              {
                s = all_qs.at(c);
              }
              charges_to_process.push_back({c, s, h});
            }
            std::sort(charges_to_process.begin(), charges_to_process.end(),
                      [](const ChargeCandidate& a, const ChargeCandidate& b) { return a.score > b.score; });
          }
          else
          {
            int charge;
            double score;
            int hcd = config_.targeting().hcd_energy;

            if (config_.targeting().use_idscore && config_.targeting().consider_all_charges && config_.targeting().hcd_energy < 0) {
              charge = pg.getBestIDScoreCharge();
              score = pg.getBestIDScore();
              hcd = pg.getBestIDScoreHCD();
            }
            else if (config_.targeting().use_idscore && config_.targeting().consider_all_charges) {
              charge = pg.getBestIDScoreChargeForHCD(config_.targeting().hcd_energy);
              score = pg.getBestIDScoreForHCD(config_.targeting().hcd_energy);
            }
            else if (config_.targeting().use_idscore && !config_.targeting().consider_all_charges && config_.targeting().hcd_energy < 0) {
              charge = pg.getRepAbsCharge();
              score = pg.getBestIDScoreForCharge(charge);
              hcd = pg.getBestHCDForCharge(charge);
            }
            else if (config_.targeting().use_idscore && !config_.targeting().consider_all_charges) {
              charge = pg.getRepAbsCharge();
              score = pg.getIDScoreForChargeAndHCD(charge, config_.targeting().hcd_energy);
            }
            else if (!config_.targeting().use_idscore && config_.targeting().consider_all_charges) {
              charge = pg.getBestQScoreCharge();
              score = pg.getBestQScore();
            }
            else {
              charge = pg.getRepAbsCharge();
              score = pg.getQscore();
            }
            charges_to_process.push_back({charge, score, hcd});
          }

          for (const auto& cc : charges_to_process)
          {
            if (selected_peak_groups_.size() >= mass_count) { break; }
            int charge = cc.charge;
            double score = cc.score;
            int hcd = cc.hcd;

            // Per-level charge filter: ms1.min_charge controls what MS1 picks
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
                    if (!inclusion_targets_.empty())  // TSV mode: check active_targets_ with charge
                    {
                      for (const auto* t : active_targets_)
                      {
                        if (t->charge < 0) {
                          target_matched = true;
                          break;
                        }
                        auto [min_charge, max_charge] = pg.getAbsChargeRange();
                        if ((t->charge >= min_charge) && (t->charge <= max_charge))
                        {
                          // Update with matched charge
                          charge = t->charge;
                          std::tie(mz1, mz2) = pg.getMzRange(charge);
                          center_mz = (mz1 + mz2) / 2.0;
                          mz1 -= optimal_window_margin_;
                          mz2 += optimal_window_margin_;
                          integer_mz = (int)round(center_mz);

                          target_matched = true;
                          break;
                        }
                      }
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
            // deep mode
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

            if (score < qscore_threshold) { break; }

            // TODO: Check
            if (pg.getChargeSNR(charge) < snr_threshold) { continue; }

            if (current_selected_mzs.find(center_mz) != current_selected_mzs.end()) // mz has been triggered
            {
              if (selection_phase < selection_phase_end) { continue; }
              if (! target_matched && current_selected_masses.find(pg.getMonoMass()) == current_selected_masses.end()) // but mass is different
              {
                continue;
              }
            }

            if (config_.targeting().charge_based_exclusion
                && tqscore_exceeding_mass_charge_set_.count({nominal_mass, charge}) > 0)
            {
              continue;
            }

            // selection phase 0, skip masses over tqscore threshold
            if (selection_phase < selection_phase_end - 1)
            {
              if (tqscore_exceeding_mass_rt_map_.find(nominal_mass) != tqscore_exceeding_mass_rt_map_.end()
                  || tqscore_exceeding_mz_rt_map_.find(integer_mz) != tqscore_exceeding_mz_rt_map_.end())
              {
                continue;
              }
            }

            // save mass acquisition
            all_mass_rt_map_[nominal_mass] = rt;

            if (config_.targeting().charge_based_exclusion)
            {
              // Per-(mass, charge) accumulation. No mass-level writes — the mass is never globally excluded.
              const auto key = std::make_pair(nominal_mass, charge);
              if (!config_.targeting().use_idscore) {
                auto inter = mass_charge_qscore_map_.find(key);
                if (inter == mass_charge_qscore_map_.end())
                {
                  mass_charge_qscore_map_[key] = score;
                }
                else {
                  mass_charge_qscore_map_[key] = std::max(inter->second, score);
                }
                if (mass_charge_qscore_map_[key] > config_.targeting().tqscore_threshold)
                {
                  tqscore_exceeding_mass_charge_set_.insert(key);
                }
              }
              else {
                auto inter = mass_charge_qscore_map_.find(key);
                if (inter == mass_charge_qscore_map_.end()) { mass_charge_qscore_map_[key] = 1 - score; }
                else { mass_charge_qscore_map_[key] *= 1 - score; }
                if (1 - mass_charge_qscore_map_[key] * tqscore_factor_for_exclusion > config_.targeting().tqscore_threshold)
                {
                  tqscore_exceeding_mass_charge_set_.insert(key);
                }
              }
            }
            else if (!config_.targeting().use_idscore) {
              // Compute total qscore
              auto inter = mass_qscore_map_.find(nominal_mass);
              if (inter == mass_qscore_map_.end())
              {
                mass_qscore_map_[nominal_mass] = score;
              }
              else {
                // If mass has previously been acquired with higher qscore, skip
                if (score < mass_qscore_map_[nominal_mass]) {
                  continue;
                }
                mass_qscore_map_[nominal_mass] = score;
              }

              // Add to exclusion list if neccessary
              if (mass_qscore_map_[nominal_mass] > config_.targeting().tqscore_threshold)
              {
                tqscore_exceeding_mass_rt_map_[nominal_mass] = rt;
                tqscore_exceeding_mz_rt_map_[integer_mz] = rt;
              }
            }
            else {
              // Compute total qscore
              auto inter = mass_qscore_map_.find(nominal_mass);
              if (inter == mass_qscore_map_.end()) { mass_qscore_map_[nominal_mass] = 1 - score; }
              else { mass_qscore_map_[nominal_mass] *= 1 - score; }

              // Add to exclusion list if neccessary
              if (1 - mass_qscore_map_[nominal_mass] * tqscore_factor_for_exclusion > config_.targeting().tqscore_threshold)
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

            // Store acquisition
            id_mass_map_[window_id_] = nominal_mass;
            id_mz_map_[window_id_] = integer_mz;
            id_qscore_map_[window_id_] = score;
            id_charge_map_[window_id_] = charge;
            trigger_ids_.push_back(window_id_);
            window_id_++;

            selected_peak_groups_.push_back(pg);
            trigger_charges_.push_back(charge);
            trigger_hcds_.push_back(hcd);
            trigger_scores_.push_back(score);

            trigger_left_isolation_mzs_.push_back(mz1);
            trigger_right_isolation_mzs_.push_back(mz2);
            current_selected_masses.insert(pg.getMonoMass());
            current_selected_mzs.insert(center_mz);
          }  // end for charges_to_process
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

    // Remove qscore from further calculations
    if (mass_qscore_map_.find(nominal_mass) != mass_qscore_map_.end()) { mass_qscore_map_[nominal_mass] /= 1 - qscore; }
  }

  void PrecursorSelection::getIsolationWindows(double* wstart,
                                               double* wend,
                                               double* qscores,
                                               int* charges,
                                               int* min_charges,
                                               int* max_charges,
                                               double* mono_masses,
                                               double* charge_cos,
                                               double* charge_snrs,
                                               double* iso_cos,
                                               double* snrs,
                                               double* charge_scores,
                                               double* ppm_errors,
                                               double* precursor_intensities,
                                               double* peakgroup_intensities,
                                               int* hcds,
                                               int* ids)
  {
    for (Size i = 0; i < selected_peak_groups_.size(); i++)
    {
      if (trigger_charges_[i] == 0) { continue; }
      auto peakgroup = selected_peak_groups_[i];
      charges[i] = trigger_charges_[i];
      auto cr = peakgroup.getAbsChargeRange();
      min_charges[i] = std::get<0>(cr);
      max_charges[i] = std::get<1>(cr);

      wstart[i] = trigger_left_isolation_mzs_[i];
      wend[i] = trigger_right_isolation_mzs_[i];

      qscores[i] = trigger_scores_[i];
      mono_masses[i] = peakgroup.getMonoMass();
      charge_cos[i] = peakgroup.getChargeIsotopeCosine(charges[i]);
      charge_snrs[i] = peakgroup.getChargeSNR(charges[i]);
      iso_cos[i] = peakgroup.getIsotopeCosine();
      snrs[i] = peakgroup.getSNR();
      charge_scores[i] = peakgroup.getChargeScore();
      ppm_errors[i] = peakgroup.getAvgPPMError();
      peakgroup_intensities[i] = peakgroup.getIntensity();
      precursor_intensities[i] = peakgroup.getChargeIntensity(charges[i]);
      hcds[i] = trigger_hcds_[i];
      ids[i] = trigger_ids_[i];
    }
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
      // Charge: -1 or empty means "any charge"
      String charge_str = cells[1].trim();
      target.charge = (charge_str.empty() || charge_str == "-1") ? -1 : cells[1].toInt();
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
      target.charge = -1;  // Any charge
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

  bool PrecursorSelection::processMS2ForTagBasedTargeting(double precursor_mass, const std::string& activation_type)
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
    Param tagger_param = tagger.getDefaults();
    tagger_param.setValue("min_length", config_.targeting().min_tag_length);
    tagger_param.setValue("max_length", config_.targeting().max_tag_length);
    tagger_param.setValue("ion_type", FragmentAnalysis::getIonTypesForFragmentationMethod(activation_type));
    tagger.setParameters(tagger_param);

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
      return true; // Target protein was detected, but all masses already tracked
    }

    // Add to dynamic inclusion list
    addDynamicTargets_(truly_new_targets, deconv_.storedMS2RT(), 100); // High priority

    return true;
  }

} // namespace OpenMS
