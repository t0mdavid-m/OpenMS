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
#include <cstring>
#include <iostream>

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
    cmd.agc_target = config_.level(1).scans[0].agc_target;
    cmd.first_mass = config_.level(1).scans[0].first_mass;
    cmd.last_mass = config_.level(1).scans[0].last_mass;
    cmd.max_it = config_.level(1).scans[0].max_it;
    std::strncpy(cmd.analyzer, "IonTrap", sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    std::strncpy(cmd.scan_description, "A", sizeof(cmd.scan_description) - 1);
    cmd.scan_description[sizeof(cmd.scan_description) - 1] = '\0';

    return cmd;
  }

  ScanCommand ScanCommandQueue::buildMS2(const PeakGroup& pg, int charge, const ScanConfig& scan_config, int priority)
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

    // Instrument defaults from MS1 config
    cmd.agc_target = config_.level(1).scans[0].agc_target;
    cmd.first_mass = config_.level(1).scans[0].first_mass;
    cmd.last_mass = config_.level(1).scans[0].last_mass;
    cmd.max_it = config_.level(1).scans[0].max_it;

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

    // Scan description: {3-char ID}R{mass_kDa:.1f}@{charge}
    std::string id_str = encode(id);
    std::snprintf(cmd.scan_description, 16, "%sR%.1f@%d", id_str.c_str(), pg.getMonoMass() / 1000.0, charge);

    // Timestamp
    cmd.enqueue_timestamp_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

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

    // Store in pending map for MS2 tracking resolution
    pending_scan_map_[id] = cmd;

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
                                          char ion_type, int frag_index, int priority)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    ScanCommand cmd{};
    int id = nextTrackingIdInt_();
    cmd.scan_id = id;
    cmd.msn_level = 3;
    cmd.priority = priority;
    cmd.is_agc = 0;
    cmd.num_stages = 2;

    // Copy analyzer/resolution from MS2 context
    cmd.orbitrap_resolution = ms2_ctx.orbitrap_resolution;
    std::strncpy(cmd.analyzer, ms2_ctx.analyzer, sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    cmd.agc_target = ms2_ctx.agc_target;
    cmd.first_mass = ms2_ctx.first_mass;
    cmd.last_mass = ms2_ctx.last_mass;
    cmd.max_it = ms2_ctx.max_it;

    // Stage 0: MS2 precursor (from MS2 context)
    cmd.stages[0] = ms2_ctx.stages[0];

    // Stage 1: Fragment target
    cmd.stages[1].precursor_mz = frag_mz;
    cmd.stages[1].isolation_width = iso_width;
    cmd.stages[1].charge_state = frag_charge;
    cmd.stages[1].collision_energy = static_cast<double>(ms3_config.collision_energy);
    std::strncpy(cmd.stages[1].activation_type, ms3_config.activation.c_str(), sizeof(cmd.stages[1].activation_type) - 1);
    cmd.stages[1].activation_type[sizeof(cmd.stages[1].activation_type) - 1] = '\0';

    // Description: {3-char ID}R{frag_mass_kDa:.1f}@{frag_charge}[{ion_type}{frag_index}]
    std::string id_str = encode(id);
    double frag_mass_kda = frag_mz * frag_charge / 1000.0;
    if (ion_type != '\0' && frag_index > 0)
      std::snprintf(cmd.scan_description, 16, "%sR%.1f@%d%c%d",
                    id_str.c_str(), frag_mass_kda, frag_charge, ion_type, frag_index);
    else
      std::snprintf(cmd.scan_description, 16, "%sR%.1f@%d",
                    id_str.c_str(), frag_mass_kda, frag_charge);

    cmd.enqueue_timestamp_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    // D6: Store in pending map for future MS3 resolution
    pending_scan_map_[id] = cmd;

    std::cout << "[TRACK-CREATE] id=" << id_str
              << " ms_level=3"
              << " frag_mz=" << frag_mz
              << std::endl;

    return cmd;
  }

  ScanCommand ScanCommandQueue::buildFollowUp(const ScanCommand& ctx,
      const ScanConfig& follow_up_config, char suffix, int priority)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    ScanCommand cmd = ctx;
    cmd.scan_id = nextTrackingIdInt_();
    cmd.priority = priority;

    std::strncpy(cmd.analyzer, follow_up_config.analyzer.c_str(), sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    cmd.orbitrap_resolution = follow_up_config.resolution;
    cmd.stages[0].collision_energy = static_cast<double>(follow_up_config.collision_energy);
    std::strncpy(cmd.stages[0].activation_type, follow_up_config.activation.c_str(),
                 sizeof(cmd.stages[0].activation_type) - 1);
    cmd.stages[0].activation_type[sizeof(cmd.stages[0].activation_type) - 1] = '\0';

    std::string id_str = encode(cmd.scan_id);
    std::snprintf(cmd.scan_description, 16, "%s%c%.1f@%d",
                  id_str.c_str(), suffix, cmd.mono_mass / 1000.0, cmd.stages[0].charge_state);

    cmd.enqueue_timestamp_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    pending_scan_map_[cmd.scan_id] = cmd;

    std::cout << "[TRACK-CREATE] id=" << id_str
              << " ms_level=2 type=followup_" << suffix
              << std::endl;

    return cmd;
  }

  // --- Queue operations ---

  void ScanCommandQueue::push(ScanCommand cmd)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
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
        return cmd;
      }
    }
    return std::nullopt;
  }

  void ScanCommandQueue::registerPending(int id, ScanCommand cmd)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    pending_scan_map_[id] = cmd;
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

  void ScanCommandQueue::cleanupExpired()
  {
    // config_ is const ref, immutable after construction — no lock needed for this check
    if (!config_.scheduling().timeout_enabled)
      return;

    std::lock_guard<std::mutex> lock(queue_mutex_);
    auto now_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    auto it = pending_scan_map_.begin();
    while (it != pending_scan_map_.end())
    {
      if (it->second.enqueue_timestamp_ms > 0 &&
          (now_ms - it->second.enqueue_timestamp_ms) > static_cast<uint64_t>(config_.scheduling().timeout_ms))
      {
        std::string id_str = encode(it->first);
        std::cout << "[TRACK-EXPIRE] id=" << id_str
                  << " age_ms=" << (now_ms - it->second.enqueue_timestamp_ms)
                  << std::endl;
        it = pending_scan_map_.erase(it);
      }
      else
      {
        ++it;
      }
    }
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

  // --- Overrides ---

  void ScanCommandQueue::applyOverrides(ScanCommand& cmd,
      const std::unordered_map<std::string, std::string>& overrides) const
  {
    for (auto ov_it = overrides.begin(); ov_it != overrides.end(); ++ov_it)
    {
      const auto& okey = ov_it->first;
      const auto& oval = ov_it->second;
      if (okey == "agc_target") cmd.agc_target = static_cast<int32_t>(std::stod(oval));
      else if (okey == "max_injection_time_ms") cmd.max_it = std::stod(oval);
      else if (okey == "isolation_width") cmd.stages[0].isolation_width = std::stod(oval);
    }
  }

} // namespace OpenMS
