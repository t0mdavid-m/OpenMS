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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string>

namespace OpenMS
{

  /// optimal window margin (same constant used in FLASHIda.cpp)
  static const double optimal_window_margin_ = .4;

  // Static member definition: all 94 printable ASCII characters (0x21-0x7E)
  const std::string ScanCommandQueue::tracking_alphabet_ = "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

  ScanCommandQueue::ScanCommandQueue(const Config& config) :
      config_(config)
  {
    last_agc_time_ = std::chrono::steady_clock::now();
  }

  // --- Tracking ID encoding/decoding ---

  std::string ScanCommandQueue::encode(int value)
  {
    const int base = static_cast<int>(tracking_alphabet_.size());
    char buf[4] = {tracking_alphabet_[0], tracking_alphabet_[0], tracking_alphabet_[0], '\0'};
    for (int i = 2; i >= 0; --i)
    {
      buf[i] = tracking_alphabet_[value % base];
      value /= base;
    }
    return std::string(buf);
  }

  int ScanCommandQueue::decode(const std::string& s) const
  {
    const int base = static_cast<int>(tracking_alphabet_.size());
    int value = 0;
    for (char c : s)
    {
      auto pos = tracking_alphabet_.find(c);
      if (pos == std::string::npos) return -1;
      value = value * base + static_cast<int>(pos);
    }
    return value;
  }

  std::string ScanCommandQueue::formatMassToken(double mass_kda, int charge, char ion_type, int ion_index)
  {
    // 15-char scan_description budget: id(3) + marker(1) + massToken(this) + '@'(1) + charge + ion.
    // Pack the most decimals that still leave the trailing ion intact (so the cap never truncates the ion).
    const int charge_digits = static_cast<int>(std::to_string(charge < 0 ? -charge : charge).size());
    const int ion_part = (ion_type != '\0' && ion_index > 0)
                           ? 1 + static_cast<int>(std::to_string(ion_index).size()) : 0;
    int token_budget = 15 - 5 - charge_digits - ion_part;  // 5 = id(3) + marker(1) + '@'(1)
    if (token_budget < 2) token_budget = 2;                // always allow at least "Nk"
    char buf[24];
    for (int d = 6; d >= 0; --d)
    {
      std::snprintf(buf, sizeof(buf), "%.*fk", d, mass_kda);
      if (static_cast<int>(std::strlen(buf)) <= token_budget) return std::string(buf);
    }
    std::snprintf(buf, sizeof(buf), "%.0fk", mass_kda);
    return std::string(buf);
  }

  int ScanCommandQueue::nextTrackingIdInt_()
  {
    const int base = static_cast<int>(tracking_alphabet_.size());
    const int max_id = base * base * base - 1;
    int id = tracking_id_counter_++;
    if (tracking_id_counter_ > max_id)
      tracking_id_counter_ = 0;
    return id;
  }

  int ScanCommandQueue::nextTrackingId()
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    return nextTrackingIdInt_();
  }

  // --- AGC and timing ---

  bool ScanCommandQueue::needsAGC() const
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_agc_time_).count();
    return static_cast<uint64_t>(elapsed) > config_.scheduling().agc_interval_ms;
  }

  uint64_t ScanCommandQueue::msSinceLastMS1() const
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    auto now = std::chrono::steady_clock::now();
    return static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(now - last_ms1_time_).count());
  }

  void ScanCommandQueue::recordMS1Time()
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    last_ms1_time_ = std::chrono::steady_clock::now();
  }

  void ScanCommandQueue::recordAGCTime()
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    last_agc_time_ = std::chrono::steady_clock::now();
  }

  // --- Command builders ---

  ScanCommand ScanCommandQueue::makeMS1() const
  {
    ScanCommand cmd{};
    cmd.msn_level = 1;
    cmd.priority = 3; // lowest priority -- MS1 is the fallback
    cmd.is_agc = 0;
    cmd.num_stages = 0;
    cmd.orbitrap_resolution = config_.level(1).scans[0].resolution;
    cmd.agc_target = config_.level(1).scans[0].agc_target;
    cmd.first_mass = config_.level(1).scans[0].first_mass;
    cmd.last_mass = config_.level(1).scans[0].last_mass;
    cmd.max_it = config_.level(1).scans[0].max_it;
    cmd.microscans = config_.level(1).scans[0].microscans;
    cmd.rf_lens = config_.level(1).scans[0].rf_lens;
    cmd.source_cid = config_.level(1).scans[0].source_cid;
    cmd.source_cid_scaling = config_.level(1).scans[0].source_cid_scaling;
    std::strncpy(cmd.data_type, config_.level(1).scans[0].data_type.c_str(), sizeof(cmd.data_type) - 1);
    cmd.data_type[sizeof(cmd.data_type) - 1] = '\0';
    std::strncpy(cmd.scan_rate, config_.level(1).scans[0].scan_rate.c_str(), sizeof(cmd.scan_rate) - 1);
    cmd.scan_rate[sizeof(cmd.scan_rate) - 1] = '\0';

    // Copy analyzer string safely
    std::strncpy(cmd.analyzer, config_.level(1).scans[0].analyzer.c_str(), sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';

    std::strncpy(cmd.scan_description, "S", sizeof(cmd.scan_description) - 1);
    cmd.scan_description[sizeof(cmd.scan_description) - 1] = '\0';

    return cmd;
  }

  ScanCommand ScanCommandQueue::makeAGC() const
  {
    ScanCommand cmd{};
    cmd.msn_level = 1;
    cmd.priority = 0; // highest priority -- AGC is time-critical
    cmd.is_agc = 1;
    cmd.num_stages = 0;
    cmd.orbitrap_resolution = 0;
    cmd.agc_target = 30000;
    cmd.first_mass = config_.level(1).scans[0].first_mass;
    cmd.last_mass = config_.level(1).scans[0].last_mass;
    cmd.max_it = 1;
    std::strncpy(cmd.scan_rate, "Turbo", sizeof(cmd.scan_rate) - 1);
    std::strncpy(cmd.analyzer, "IonTrap", sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    std::strncpy(cmd.scan_description, "A", sizeof(cmd.scan_description) - 1);
    cmd.scan_description[sizeof(cmd.scan_description) - 1] = '\0';

    return cmd;
  }

  ScanCommand ScanCommandQueue::buildMS2(const PeakGroup& pg, int charge, const ScanConfig& scan_config, int priority, int parent_scan_id)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    ScanCommand cmd{};
    int id = nextTrackingIdInt_();
    cmd.scan_id = id;
    cmd.msn_level = 2;
    cmd.priority = priority;
    cmd.is_agc = 0;
    cmd.num_stages = 1;

    // Resolution and analyzer from the provided ScanConfig
    cmd.orbitrap_resolution = scan_config.resolution;
    std::strncpy(cmd.analyzer, scan_config.analyzer.c_str(), sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';

    // Mass range from MS2 config (0 = unset, Thermo API inherits from method default)
    cmd.first_mass = scan_config.first_mass;
    cmd.last_mass = scan_config.last_mass;

    // Scan properties from MS2 config
    cmd.agc_target = config_.level(2).scans[0].agc_target;
    cmd.max_it = config_.level(2).scans[0].max_it;

    // New scan parameters from MS2 config
    cmd.microscans = scan_config.microscans;
    cmd.rf_lens = scan_config.rf_lens;
    cmd.source_cid = scan_config.source_cid;
    cmd.source_cid_scaling = scan_config.source_cid_scaling;
    std::strncpy(cmd.data_type, scan_config.data_type.c_str(), sizeof(cmd.data_type) - 1);
    cmd.data_type[sizeof(cmd.data_type) - 1] = '\0';
    std::strncpy(cmd.scan_rate, scan_config.scan_rate.c_str(), sizeof(cmd.scan_rate) - 1);
    cmd.scan_rate[sizeof(cmd.scan_rate) - 1] = '\0';

    // Isolation window from peak group m/z range
    auto [mz1, mz2] = pg.getMzRange(charge);
    double center_mz = (mz1 + mz2) / 2.0;
    mz1 -= optimal_window_margin_;
    mz2 += optimal_window_margin_;
    double iso_width = mz2 - mz1;

    cmd.stages[0].precursor_mz = center_mz;
    cmd.stages[0].isolation_width = iso_width;
    cmd.stages[0].charge_state = charge;

    // CE and activation from scan_config
    cmd.stages[0].collision_energy = static_cast<double>(scan_config.collision_energy);
    std::strncpy(cmd.stages[0].activation_type, scan_config.activation.c_str(),
                 sizeof(cmd.stages[0].activation_type) - 1);
    cmd.stages[0].activation_type[sizeof(cmd.stages[0].activation_type) - 1] = '\0';
    cmd.stages[0].reaction_time = scan_config.reaction_time;
    cmd.stages[0].reagent_max_it = scan_config.reagent_max_it;
    cmd.stages[0].reagent_agc_target = scan_config.reagent_agc_target;

    // Scan description: {3-char ID}R{mass_kDa}k@{charge}  (adaptive-precision mass token)
    std::string id_str = encode(id);
    std::string mass_tok = formatMassToken(pg.getMonoMass() / 1000.0, charge, '\0', 0);
    std::snprintf(cmd.scan_description, 16, "%sR%s@%d", id_str.c_str(), mass_tok.c_str(), charge);

    // Precursor scoring data for diagnostic TSV output
    cmd.qscore = pg.getQscore();
    cmd.mono_mass = pg.getMonoMass();
    cmd.charge_cos = pg.getChargeIsotopeCosine(std::abs(charge));
    cmd.charge_snr = pg.getChargeSNR(std::abs(charge));
    cmd.iso_cos = pg.getIsotopeCosine();
    cmd.snr = pg.getSNR();
    cmd.charge_score = pg.getChargeScore();
    cmd.ppm_error = pg.getAvgPPMError();
    cmd.precursor_intensity = pg.getChargeIntensity(std::abs(charge));
    cmd.peakgroup_intensity = pg.getIntensity();
    cmd.hcd_energy = scan_config.collision_energy;
    cmd.pad2 = 0;

    // Parent tracking ID
    if (parent_scan_id > 0)
    {
      std::string parent_enc = encode(parent_scan_id);
      std::strncpy(cmd.parent_scan_id, parent_enc.c_str(), 3);
      cmd.parent_scan_id[3] = '\0';
    }

    std::cout << "[TRACK-CREATE] id=" << id_str
              << " ms_level=2"
              << " mz=" << center_mz
              << " z=" << charge
              << " mass=" << pg.getMonoMass()
              << std::endl;

    return cmd;
  }

  ScanCommand ScanCommandQueue::buildMS3(const ScanCommand& ms2_ctx, const ScanConfig& ms3_config,
                                          double frag_mz, int frag_charge, double iso_width,
                                          char ion_type, int frag_index, int priority,
                                          const FragmentAnalysis::FragmentScores& frag_scores)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    ScanCommand cmd{};
    int id = nextTrackingIdInt_();
    cmd.scan_id = id;
    cmd.msn_level = 3;
    cmd.priority = priority;
    cmd.is_agc = 0;
    cmd.num_stages = 2;

    // Mass range from MS3 config (0 = unset, Thermo API inherits from method default)
    cmd.first_mass = ms3_config.first_mass;
    cmd.last_mass = ms3_config.last_mass;

    cmd.max_it = ms3_config.max_it;
    cmd.agc_target = ms3_config.agc_target;
    cmd.orbitrap_resolution = ms3_config.resolution;
    std::strncpy(cmd.analyzer, ms3_config.analyzer.c_str(), sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';

    // New scan parameters from MS3 config
    cmd.microscans = ms3_config.microscans;
    cmd.rf_lens = ms3_config.rf_lens;
    cmd.source_cid = ms3_config.source_cid;
    cmd.source_cid_scaling = ms3_config.source_cid_scaling;
    std::strncpy(cmd.data_type, ms3_config.data_type.c_str(), sizeof(cmd.data_type) - 1);
    cmd.data_type[sizeof(cmd.data_type) - 1] = '\0';
    std::strncpy(cmd.scan_rate, ms3_config.scan_rate.c_str(), sizeof(cmd.scan_rate) - 1);
    cmd.scan_rate[sizeof(cmd.scan_rate) - 1] = '\0';

    // Parent tracking ID from MS2 context
    std::string parent_enc = encode(ms2_ctx.scan_id);
    std::strncpy(cmd.parent_scan_id, parent_enc.c_str(), 3);
    cmd.parent_scan_id[3] = '\0';

    // Stage 0: MS2 precursor (from MS2 context)
    cmd.stages[0] = ms2_ctx.stages[0];

    // Stage 1: Fragment target
    cmd.stages[1].precursor_mz = frag_mz;
    cmd.stages[1].isolation_width = std::max(iso_width, 2.0);
    cmd.stages[1].charge_state = frag_charge;
    cmd.stages[1].collision_energy = static_cast<double>(ms3_config.collision_energy);
    std::strncpy(cmd.stages[1].activation_type, ms3_config.activation.c_str(), sizeof(cmd.stages[1].activation_type) - 1);
    cmd.stages[1].activation_type[sizeof(cmd.stages[1].activation_type) - 1] = '\0';
    cmd.stages[1].reaction_time = ms3_config.reaction_time;
    cmd.stages[1].reagent_max_it = ms3_config.reagent_max_it;
    cmd.stages[1].reagent_agc_target = ms3_config.reagent_agc_target;

    // Description: {3-char ID}R{frag_mass_kDa}k@{frag_charge}[{ion_type}{frag_index}]  (adaptive-precision)
    std::string id_str = encode(id);
    double frag_mass_kda = frag_mz * frag_charge / 1000.0;
    if (ion_type != '\0' && frag_index > 0)
    {
      std::string mass_tok = formatMassToken(frag_mass_kda, frag_charge, ion_type, frag_index);
      std::snprintf(cmd.scan_description, 16, "%sR%s@%d%c%d",
                    id_str.c_str(), mass_tok.c_str(), frag_charge, ion_type, frag_index);
    }
    else
    {
      std::string mass_tok = formatMassToken(frag_mass_kda, frag_charge, '\0', 0);
      std::snprintf(cmd.scan_description, 16, "%sR%s@%d",
                    id_str.c_str(), mass_tok.c_str(), frag_charge);
    }

    std::cout << "[TRACK-CREATE] id=" << id_str
              << " ms_level=3"
              << " frag_mz=" << frag_mz
              << std::endl;

    // 2-stage scoring for the MS3 commands row: stage-0 = the MS2 precursor (carried from ms2_ctx),
    // stage-1 = the selected fragment precursor (frag_scores, re-computed from its PeakGroup).
    cmd.qscore               = ms2_ctx.qscore;               cmd.qscore_s1               = frag_scores.qscore;
    cmd.mono_mass            = ms2_ctx.mono_mass;            cmd.mono_mass_s1            = frag_scores.mono_mass;
    cmd.charge_cos           = ms2_ctx.charge_cos;           cmd.charge_cos_s1           = frag_scores.charge_cos;
    cmd.charge_snr           = ms2_ctx.charge_snr;           cmd.charge_snr_s1           = frag_scores.charge_snr;
    cmd.iso_cos              = ms2_ctx.iso_cos;              cmd.iso_cos_s1              = frag_scores.iso_cos;
    cmd.snr                  = ms2_ctx.snr;                  cmd.snr_s1                  = frag_scores.snr;
    cmd.charge_score         = ms2_ctx.charge_score;         cmd.charge_score_s1         = frag_scores.charge_score;
    cmd.ppm_error            = ms2_ctx.ppm_error;            cmd.ppm_error_s1            = frag_scores.ppm_error;
    cmd.precursor_intensity  = ms2_ctx.precursor_intensity;  cmd.precursor_intensity_s1  = frag_scores.precursor_intensity;
    cmd.peakgroup_intensity  = ms2_ctx.peakgroup_intensity;  cmd.peakgroup_intensity_s1  = frag_scores.peakgroup_intensity;
    cmd.hcd_energy           = ms2_ctx.hcd_energy;           cmd.hcd_energy_s1           = static_cast<int32_t>(ms3_config.collision_energy);

    return cmd;
  }

  ScanCommand ScanCommandQueue::buildFollowUp(const ScanCommand& ctx,
      const ScanConfig& follow_up_config, char suffix, int priority)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    ScanCommand cmd = ctx;
    cmd.scan_id = nextTrackingIdInt_();
    cmd.priority = priority;

    // Follow-up's logical parent is the triggering MS2 (ctx) itself, not the grandparent MS1 that
    // ctx.parent_scan_id points to (#8). Mirror buildMS3 and encode ctx's own scan_id as the parent.
    std::string parent_enc = encode(ctx.scan_id);
    std::strncpy(cmd.parent_scan_id, parent_enc.c_str(), 3);
    cmd.parent_scan_id[3] = '\0';

    std::strncpy(cmd.analyzer, follow_up_config.analyzer.c_str(), sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    cmd.orbitrap_resolution = follow_up_config.resolution;
    cmd.stages[0].collision_energy = static_cast<double>(follow_up_config.collision_energy);
    std::strncpy(cmd.stages[0].activation_type, follow_up_config.activation.c_str(),
                 sizeof(cmd.stages[0].activation_type) - 1);
    cmd.stages[0].activation_type[sizeof(cmd.stages[0].activation_type) - 1] = '\0';

    std::string id_str = encode(cmd.scan_id);
    std::string mass_tok = formatMassToken(cmd.mono_mass / 1000.0, cmd.stages[0].charge_state, '\0', 0);
    std::snprintf(cmd.scan_description, 16, "%s%c%s@%d",
                  id_str.c_str(), suffix, mass_tok.c_str(), cmd.stages[0].charge_state);

    std::cout << "[TRACK-CREATE] id=" << id_str
              << " ms_level=2 type=followup_" << suffix
              << std::endl;

    return cmd;
  }

  // --- Queue operations ---

  void ScanCommandQueue::push(ScanCommand cmd)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    cmd.enqueue_timestamp_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());
    int p = std::clamp(cmd.priority, 0, 3);
    queues_[p].push_back(cmd);
  }

  std::optional<ScanCommand> ScanCommandQueue::dequeue()
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    for (int p = 0; p < 4; ++p)
    {
      if (!queues_[p].empty())
      {
        ScanCommand cmd = queues_[p].front();
        queues_[p].pop_front();
        cmd.dequeue_timestamp_ms = static_cast<uint64_t>(
          std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now().time_since_epoch()).count());
        pending_scan_map_[cmd.scan_id] = cmd;
        return cmd;
      }
    }
    return std::nullopt;
  }

  void ScanCommandQueue::registerPending(const ScanCommand& cmd)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    pending_scan_map_[cmd.scan_id] = cmd;
  }

  std::optional<ScanCommand> ScanCommandQueue::resolvePending(int id)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    auto it = pending_scan_map_.find(id);
    if (it == pending_scan_map_.end())
      return std::nullopt;
    ScanCommand cmd = it->second;
    pending_scan_map_.erase(it);
    return cmd;
  }

  void ScanCommandQueue::setWindowSnr(int scan_id, double value)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    window_snr_map_[scan_id] = value;
  }

  double ScanCommandQueue::windowSnr(int scan_id) const
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    auto it = window_snr_map_.find(scan_id);
    return it == window_snr_map_.end() ? -1.0 : it->second;
  }

  void ScanCommandQueue::cleanupExpired()
  {
    if (!config_.scheduling().timeout_enabled)
      return;

    std::lock_guard<std::mutex> lock(queue_mutex_);
    auto now_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    for (int p = 0; p < 4; ++p)
    {
      auto it = queues_[p].begin();
      while (it != queues_[p].end())
      {
        if (it->enqueue_timestamp_ms > 0 &&
            (now_ms - it->enqueue_timestamp_ms) > static_cast<uint64_t>(config_.scheduling().timeout_ms))
        {
          std::string id_str = encode(it->scan_id);
          std::cout << "[TRACK-EXPIRE] id=" << id_str
                    << " age_ms=" << (now_ms - it->enqueue_timestamp_ms)
                    << " ms_level=" << it->msn_level
                    << std::endl;
          it = queues_[p].erase(it);
        }
        else
        {
          ++it;
        }
      }
    }
  }

  std::vector<int> ScanCommandQueue::cancelByScanIds(const std::vector<int>& scan_ids)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    std::vector<int> removed;

    auto is_target = [&scan_ids](int id) {
      return std::find(scan_ids.begin(), scan_ids.end(), id) != scan_ids.end();
    };

    // Remove matching commands still waiting in the priority queues.
    for (int p = 0; p < 4; ++p)
    {
      auto it = queues_[p].begin();
      while (it != queues_[p].end())
      {
        if (is_target(it->scan_id))
        {
          removed.push_back(it->scan_id);
          std::cout << "[TRACK-CANCEL] id=" << encode(it->scan_id)
                    << " ms_level=" << it->msn_level << " state=queued" << std::endl;
          it = queues_[p].erase(it);
        }
        else
        {
          ++it;
        }
      }
    }

    // Remove matching in-flight commands from the pending map; their late
    // results are dropped once the caller erases the matching tracking entry.
    auto pit = pending_scan_map_.begin();
    while (pit != pending_scan_map_.end())
    {
      if (is_target(pit->first))
      {
        removed.push_back(pit->first);
        std::cout << "[TRACK-CANCEL] id=" << encode(pit->first)
                  << " ms_level=" << pit->second.msn_level << " state=in-flight" << std::endl;
        pit = pending_scan_map_.erase(pit);
      }
      else
      {
        ++pit;
      }
    }

    return removed;
  }

  // --- Introspection ---

  size_t ScanCommandQueue::pendingScanMapSize() const
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    return pending_scan_map_.size();
  }

  size_t ScanCommandQueue::queueSize(int priority) const
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    if (priority < 0 || priority > 3)
      return 0;
    return queues_[priority].size();
  }

} // namespace OpenMS
