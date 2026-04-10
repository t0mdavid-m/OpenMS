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
    queue_(config_),
    deconv_(config_),
    fragments_(config_),
    selection_(config_, deconv_),
    quant_(config_),
    faims_(config_),
    exploration_(config_)
{
  #ifdef _OPENMP
    omp_set_num_threads(4);
  #endif
  }

  // Fragment analysis: C-pointer overloads delegate to fragments_ + deconv_.storedMS2()

  int FLASHIda::getBestMS2Masses(int n, double* masses, double* qscores, int* charges,
                                 double* window_starts, double* window_ends)
  {
    if (!deconv_.hasStoredMS2()) return 0;
    return fragments_.getBestMS2Masses(n, masses, qscores, charges, window_starts, window_ends,
                                       deconv_.storedMS2());
  }

  int FLASHIda::getTopFragmentMatches(const String& protein_sequence, int n,
                                      double* masses, double* qscores, int* charges,
                                      double* window_starts, double* window_ends,
                                      char* ion_types, int* fragment_indices,
                                      const String& fragmentation_method)
  {
    if (!deconv_.hasStoredMS2()) return 0;
    return fragments_.getTopFragmentMatches(protein_sequence, n, masses, qscores, charges,
                                            window_starts, window_ends, ion_types, fragment_indices,
                                            deconv_.storedMS2(), fragmentation_method);
  }

  int FLASHIda::getTerminalFragmentIons(const String& protein_sequence, int n,
                                        double* masses, double* qscores, int* charges,
                                        double* window_starts, double* window_ends,
                                        char* ion_types, int* fragment_indices,
                                        const String& fragmentation_method)
  {
    if (!deconv_.hasStoredMS2()) return 0;
    return fragments_.getTerminalFragmentIons(protein_sequence, n, masses, qscores, charges,
                                              window_starts, window_ends, ion_types, fragment_indices,
                                              deconv_.storedMS2(), fragmentation_method);
  }

  int FLASHIda::getAmbiguityEnclosingIons(const String& protein_sequence, int n,
                                          double* masses, double* qscores, int* charges,
                                          double* window_starts, double* window_ends,
                                          char* ion_types, int* fragment_indices,
                                          const String& fragmentation_method)
  {
    if (!deconv_.hasStoredMS2()) return 0;
    return fragments_.getAmbiguityEnclosingIons(protein_sequence, n, masses, qscores, charges,
                                                window_starts, window_ends, ion_types, fragment_indices,
                                                deconv_.storedMS2(), fragmentation_method);
  }

  int FLASHIda::getConfigInt(const std::string& key) const
  {
    return config_.getInt(key);
  }

  double FLASHIda::getConfigDouble(const std::string& key) const
  {
    return config_.getDouble(key);
  }

  std::map<int, std::vector<std::vector<float>>> FLASHIda::parseFLASHIdaLog(const String& in_log_file)
  {
    std::map<int, std::vector<std::vector<float>>>
      precursor_map_for_real_time_acquisition; // ms1 scan -> mass, charge ,score, mz range, precursor int, mass int, color

    if (in_log_file.empty()) { return precursor_map_for_real_time_acquisition; }

    std::ifstream f(in_log_file.c_str());
    if (! f.good())
    {
      std::cout << "FLASHIda log file " << in_log_file << " is NOT found. FLASHIda support is not active." << std::endl;
      return precursor_map_for_real_time_acquisition;
    }


    std::cout << "FLASHIda log file used: " << in_log_file << std::endl;
    std::ifstream instream(in_log_file);
    if (instream.good())
    {
      String line;
      int scan;
      float mass, charge, w1, w2, qscore, pint, mint, z1, z2;
      float features[6];
      while (std::getline(instream, line))
      {
        if (line.find("0 targets") != line.npos) { continue; }
        if (line.hasPrefix("MS1"))
        {
          Size st = line.find("MS1 Scan# ") + 10;
          Size ed = line.find(' ', st);
          String n = line.substr(st, ed);
          scan = atoi(n.c_str());
          precursor_map_for_real_time_acquisition[scan]
            = std::vector<std::vector<float>>(); //// ms1 scan -> mass, charge ,score, mz range, precursor int, mass int, color
        }
        if (line.hasPrefix("Mass"))
        {
          Size st = 5;
          Size ed = line.find('\t');
          String n = line.substr(st, ed);
          mass = (float)atof(n.c_str());

          st = line.find("Z=") + 2;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          charge = (float)atof(n.c_str());

          st = line.find("Score=") + 6;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          qscore = (float)atof(n.c_str());

          st = line.find("[") + 1;
          ed = line.find('-', st);
          n = line.substr(st, ed);
          w1 = (float)atof(n.c_str());

          st = line.find('-', ed) + 1;
          ed = line.find(']', st);
          n = line.substr(st, ed);
          w2 = (float)atof(n.c_str());

          st = line.find("PrecursorIntensity=", ed) + 19;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          pint = (float)atof(n.c_str());

          st = line.find("PrecursorMassIntensity=", ed) + 23;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          mint = (float)atof(n.c_str());

          st = line.find("Features=", ed) + 9;
          // ed = line.find(' ', st);

          st = line.find('[', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[0] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[1] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[2] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[3] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[4] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(']', st);
          n = line.substr(st, ed);
          features[5] = (float)atof(n.c_str());

          st = line.find("ChargeRange=[", ed) + 13;
          ed = line.find('-', st);
          n = line.substr(st, ed);
          z1 = (float)atof(n.c_str());

          st = line.find("-", ed) + 1;
          ed = line.find(']', st);
          n = line.substr(st, ed);
          z2 = (float)atof(n.c_str());
          std::vector<float> e(15);
          e[0] = mass;
          e[1] = charge;
          e[2] = qscore;
          e[3] = w1;
          e[4] = w2;
          e[5] = pint;
          e[6] = mint;
          e[7] = z1;
          e[8] = z2;
          for (int i = 9; i < 15; i++)
          {
            e[i] = features[i - 9];
          }
          precursor_map_for_real_time_acquisition[scan].push_back(e);
        }
      }
      instream.close();
    }
    else { std::cout << in_log_file << " not found\n"; }
    int mass_cntr = 0;
    for (auto& v : precursor_map_for_real_time_acquisition)
    {
      std::sort(v.second.begin(), v.second.end(), [](const std::vector<float>& left, const std::vector<float>& right) { return left[0] < right[0]; });
      mass_cntr += v.second.size();
    }

    std::cout << "Used precursor size : " << precursor_map_for_real_time_acquisition.size() << " precursor masses : " << mass_cntr << std::endl;

    return precursor_map_for_real_time_acquisition;
  }

  // ---- FAIMS CV adaptive skip state machine ----
  //
  // Behavioral audit of C# ScanScheduler (reference implementation):
  //
  // Q1. CV advance trigger: updateCV() called after every MS1 deconvolution with
  //     the CV of that scan and precursor count. getFAIMSMS1Scan() does actual cycling.
  // Q2. Threshold comparison: precursors < MassThreshold (strictly less than). Default=15.
  // Q3. Skip counter: CVSkipAmount[pos] doubles (min=1, cap=MaxCVSkip).
  //     CVSkipCount[pos] resets to 0 in BOTH branches of updateCV().
  //     In getFAIMSMS1Scan(), CVSkipCount[i]++ when skipping.
  // Q4. Forced cycle: when CVSkipAmount >= MaxCVSkip, skip amount stops growing.
  //     CV is visited when CVSkipCount >= CVSkipAmount (counter exhausted).
  // Q5. CV transition: getFAIMSMS1Scan() enqueues AGC+MS1 for next CV.
  // Q6. Call order: deconvolve → create MS2s → updateCV(cv, precursors) → getFAIMSMS1Scan(true).
  // Q7. Cycling: forward only (currentCV++), wraps at end. Increment-first (index 0 skipped on first call).
  // Q8. Single CV: Flash.cs uses UnifiedScanProcessor, ScanScheduler with UseFAIMS=false.
  // Q9. Non-FAIMS: faims_cv = 0.0 on all commands.

  int FLASHIda::processScan(const double* mzs, const double* ints, int length,
                             double rt_min, int ms_level, const char* scan_description,
                             double faims_cv)
  {
    std::lock_guard<std::mutex> lock(mutex_);

    if (ms_level == 1)
    {
      queue_.recordMS1Time();

      // Selection=none: skip MS1 precursor selection entirely
      if (config_.level(1).selection == SelectionMetric::None)
        return 0;

      // MS1 path: deconvolve, score, filter, select top-N, push MS2 commands
      double parent_cv = config_.faims().enabled ? faims_cv : 0.0;

      int n = selection_.filterAndRank(mzs, ints, length, rt_min, 1, nullptr);
      const auto& selected = selection_.selectedPeakGroups();
      const auto& sel_charges = selection_.triggerCharges();
      const auto& sel_hcds = selection_.triggerHcds();
      int commands_pushed = 0;
      for (int i = 0; i < n; i++)
      {
        ScanCommand cmd = queue_.buildMS2(selected[i], sel_charges[i], sel_hcds[i]);
        cmd.faims_cv = parent_cv;  // MS2 carries parent MS1's CV
        queue_.push(cmd);
        commands_pushed++;
      }

      // Initiate exploration for selected precursors if MS2 exploration is enabled
      if (config_.hasExploration(2))
      {
        for (int i = 0; i < n; i++)
        {
          auto [mz1, mz2] = selected[i].getMzRange(sel_charges[i]);
          double center_mz = (mz1 + mz2) / 2.0;
          auto cmds = exploration_.initiate(2, center_mz,
              selected[i].getMonoMass(),
              sel_charges[i], parent_cv, queue_);
          for (auto& c : cmds) queue_.push(c);
        }
      }

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

        std::string id_str = ScanCommandQueue::encode(ms1.scan_id);
        std::snprintf(ms1.scan_description, 16, "%sS", id_str.c_str());

        auto now = std::chrono::steady_clock::now();
        ms1.enqueue_timestamp_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()).count();

        queue_.push(ms1);
        std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=cv_transition cv=" << next_cv << std::endl;
      }

      return commands_pushed;
    }
    else if (ms_level == 2)
    {
      return processMS2Path_(mzs, ints, length, rt_min, scan_description);
    }
    return 0;
  }

  std::vector<FLASHIda::MS3Target> FLASHIda::selectMS3Targets_()
  {
    std::vector<MS3Target> targets;
    if (config_.targeting().ms3_mode == 0 || !deconv_.hasStoredMS2())
      return targets;

    const int n = config_.level(3).max_targets;
    std::vector<double> masses(n), qscores(n), wstarts(n), wends(n);
    std::vector<int> charges(n);
    std::vector<char> ion_types(n, '\0');
    std::vector<int> frag_indices(n, 0);

    int count = 0;

    if (config_.targeting().ms3_mode == 1 || config_.targeting().ms3_mode == 2)
    {
      // Mode 1 (SourceCID) and Mode 2 (SPS): Use getBestMS2Masses
      count = fragments_.getBestMS2Masses(n, masses.data(), qscores.data(), charges.data(),
                                          wstarts.data(), wends.data(), deconv_.storedMS2());
    }
    else if (config_.targeting().ms3_mode == 3 && !config_.targeting().protein_sequence.empty())
    {
      // Mode 3 (HCD-triggered): Use getTopFragmentMatches
      count = fragments_.getTopFragmentMatches(config_.targeting().protein_sequence, n, masses.data(), qscores.data(),
                                               charges.data(), wstarts.data(), wends.data(),
                                               ion_types.data(), frag_indices.data(), deconv_.storedMS2(), "HCD");
    }
    else if (config_.targeting().ms3_mode == 4 && !config_.targeting().protein_sequence.empty())
    {
      // Mode 4 (EThcD-triggered): Use getTerminalFragmentIons
      count = fragments_.getTerminalFragmentIons(config_.targeting().protein_sequence, n, masses.data(), qscores.data(),
                                                  charges.data(), wstarts.data(), wends.data(),
                                                  ion_types.data(), frag_indices.data(), deconv_.storedMS2(), "EThcD");
    }

    for (int i = 0; i < count; i++)
    {
      double center_mz = (wstarts[i] + wends[i]) / 2.0;
      double iso_width = wends[i] - wstarts[i];
      targets.push_back({center_mz, charges[i], iso_width, ion_types[i], frag_indices[i]});
    }
    return targets;
  }

  int FLASHIda::processMS2Path_(const double* mzs, const double* ints, int length,
                                 double rt_min, const char* scan_desc)
  {
    int commands_pushed = 0;

    // Step 1: Decode tracking ID from scan_description -- fixed position chars 0-2
    std::string desc_str = scan_desc ? scan_desc : "";
    if (desc_str.size() < 3)
      return 0;

    std::string id_str = desc_str.substr(0, 3);
    int tracking_id = queue_.decode(id_str);

    // Check if this is an exploration variant (before pending scan lookup)
    if (exploration_.isExplorationVariant(tracking_id))
    {
      // Deconvolve the MS2 result for scoring
      DeconvolvedSpectrum ms2_deconv(tracking_id);
      if (mzs != nullptr && ints != nullptr && length > 0)
      {
        deconv_.deconvolveMS2(mzs, ints, length, rt_min, 0.0, 0);
        ms2_deconv = deconv_.storedMS2();
      }

      auto cmds = exploration_.feedResult(tracking_id, ms2_deconv, rt_min, queue_);
      for (auto& c : cmds) queue_.push(c);
      return commands_pushed;
    }

    // Step 2: Look up pending scan context via queue
    auto resolved = queue_.resolvePending(tracking_id);
    if (!resolved.has_value())
    {
      std::cout << "[TRACK-RESOLVE] id=" << id_str << " status=not_found" << std::endl;
      return 0;
    }
    ScanCommand ctx = resolved.value();

    // Step 3: Deconvolve MS2 with precursor context
    double precursor_mass = 0;
    int precursor_charge = 0;
    if (ctx.num_stages > 0)
    {
      precursor_charge = ctx.stages[0].charge_state;
      precursor_mass = ctx.stages[0].precursor_mz * precursor_charge
                       - precursor_charge * FLASHHelperClasses::getChargeMass(true);
    }

    deconv_.deconvolveMS2(mzs, ints, length, rt_min, precursor_mass, precursor_charge);

    // Step 4: Route by mode
    // Tag-based targeting
    bool tags_found = false;
    if (selection_.hasTargetProteinDatabase() && precursor_mass > 0)
    {
      tags_found = selection_.processMS2ForTagBasedTargeting(precursor_mass);
    }

    // Quantification follow-up (independent of tags)
    if (config_.quantification().enabled && config_.level(2).scans.size() >= 2)
    {
      if (quant_.isDifferentiallyAbundant(mzs, ints, length, rt_min, 2, "ms2_quant",
                                          config_.quantification().reporter_mz_tol, config_.quantification().fold_change_threshold, false))
      {
        queue_.push(queue_.buildFollowUpMS2(ctx));
        commands_pushed++;
      }
    }

    // Conditional MS2 follow-up -- only when tags detected AND explicitly enabled
    if (config_.targeting().conditional_ms2_enabled && config_.level(2).scans.size() >= 2 && tags_found)
    {
      queue_.push(queue_.buildConditionalFollowUp(ctx));
      commands_pushed++;
    }

    // Step 5: MS3 targeting -- uses config levels when exploration is configured,
    // falls back to legacy MS3 targeting path for standard MS3 targeting
    if (config_.hasExploration(3))
    {
      // MS3 exploration: create exploration groups for top fragments
      auto cmds = exploration_.initiateNextLevel(2, deconv_.storedMS2(), ctx.faims_cv, queue_);
      for (auto& c : cmds) queue_.push(c);
    }
    else if (config_.targeting().ms3_enabled && config_.targeting().ms3_mode > 0)
    {
      // Legacy MS3 targeting (non-exploration)
      auto ms3_targets = selectMS3Targets_();
      for (const auto& t : ms3_targets)
      {
        ScanCommand ms3_cmd = queue_.buildMS3(ctx, t.center_mz, t.charge, t.iso_width,
                                               t.ion_type, t.frag_index);
        queue_.push(ms3_cmd);
        commands_pushed++;
      }
    }
    else if (config_.level(3).selection != SelectionMetric::None && !(config_.targeting().ms3_mode > 0))
    {
      // New selection_strategy MS3 targeting (no exploration, not legacy)
      auto cmds = exploration_.initiateNextLevel(2, deconv_.storedMS2(), ctx.faims_cv, queue_);
      for (auto& c : cmds) queue_.push(c);
    }

    std::cout << "[TRACK-RESOLVE] id=" << id_str
              << " rt=" << rt_min
              << " commands=" << commands_pushed
              << std::endl;

    return commands_pushed;
  }

  int FLASHIda::getNextScanCommand(ScanCommand& out)
  {
    std::lock_guard<std::mutex> lock(mutex_);

    // Step 1: AGC scan if needed
    if (queue_.needsAGC())
    {
      out = queue_.makeAGC();
      out.faims_cv = faims_.isEnabled() ? faims_.currentCV() : 0.0;
      out.scan_id = queue_.nextTrackingId();
      queue_.recordAGCTime();

      // Scan description: {3-char ID}A
      std::string id_str = ScanCommandQueue::encode(out.scan_id);
      std::snprintf(out.scan_description, 16, "%sA", id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=agc" << std::endl;
      return 1;
    }

    // Step 2: Cycle time -- force MS1 if too long since last survey scan
    // Suppressed while any exploration group is active
    bool exploration_active = exploration_.activeGroupCount() > 0;
    if (config_.scheduling().cycle_time_enabled && !exploration_active
        && queue_.msSinceLastMS1() > static_cast<uint64_t>(config_.scheduling().cycle_time_ms))
    {
      out = queue_.makeMS1();
      out.faims_cv = faims_.isEnabled() ? faims_.currentCV() : 0.0;
      out.scan_id = queue_.nextTrackingId();
      queue_.recordMS1Time();

      std::string id_str = ScanCommandQueue::encode(out.scan_id);
      std::snprintf(out.scan_description, 16, "%sS", id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=cycle_time" << std::endl;
      return 1;
    }

    // Step 3: Cleanup expired commands
    queue_.cleanupExpired();

    // Step 4: Dequeue by priority (0 = highest -> 3 = lowest)
    auto dequeued = queue_.dequeue();
    if (dequeued.has_value())
    {
      out = dequeued.value();
      // faims_cv already set at creation time (MS2 -> parent CV, CV-transition MS1 -> next CV)
      return 1;
    }

    // Step 5: Idle cycle -- queue empty, keep the instrument busy with AGC + MS1
    // Create an AGC command (returned immediately) and push an MS1 at priority 0
    // into the queue so the next dequeue returns it before any MS2s (priority 1+).
    {
      // 5a: AGC
      ScanCommand agc_cmd = queue_.makeAGC();
      agc_cmd.faims_cv = faims_.isEnabled() ? faims_.currentCV() : 0.0;
      agc_cmd.scan_id = queue_.nextTrackingId();
      queue_.recordAGCTime();

      std::string agc_id_str = ScanCommandQueue::encode(agc_cmd.scan_id);
      std::snprintf(agc_cmd.scan_description, 16, "%sA", agc_id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << agc_id_str << " ms_level=1 type=idle_agc" << std::endl;

      // 5b: MS1 -- override priority to 0 (makeMS1 defaults to 3)
      ScanCommand ms1_cmd = queue_.makeMS1();
      ms1_cmd.faims_cv = faims_.isEnabled() ? faims_.currentCV() : 0.0;
      ms1_cmd.scan_id = queue_.nextTrackingId();
      ms1_cmd.priority = 0;
      queue_.recordMS1Time();

      std::string ms1_id_str = ScanCommandQueue::encode(ms1_cmd.scan_id);
      std::snprintf(ms1_cmd.scan_description, 16, "%sS", ms1_id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << ms1_id_str << " ms_level=1 type=idle_ms1" << std::endl;

      // Push MS1 into priority-0 queue for next dequeue call
      queue_.push(ms1_cmd);

      out = agc_cmd;
      return 1;
    }
  }

  int FLASHIda::getNextTrackingId()
  {
    return queue_.nextTrackingId();
  }

} // namespace OpenMS
