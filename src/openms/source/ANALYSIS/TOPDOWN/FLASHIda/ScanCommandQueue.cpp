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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/NotchSelection.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <set>
#include <string>

namespace OpenMS
{

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

  void ScanCommandQueue::writeMs1Description(ScanCommand& cmd, ScanRole role)
  {
    // The 16 is the pre-existing scan_description budget the four migrated call sites all used, NOT
    // sizeof(scan_description) (256). Keeping it is what makes this byte-identical to them.
    const std::string id_str = encode(cmd.scan_id);
    std::snprintf(cmd.scan_description, 16, "%s%c", id_str.c_str(), roleMarker(role));
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
      tracking_id_counter_ = 1;  // wrap to 1, NOT 0 -- 0 is the reserved "no parent" sentinel.
                                 // See the counter's declaration in ScanCommandQueue.h.
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

  bool ScanCommandQueue::takeMonitorDue(int idx, bool sweeping, double interval_ms)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    if (idx < 0 || idx > 1) return false;

    // Not sweeping (or not configured): RE-ARM. This is the whole reason the caller must reach here
    // on every drain rather than only while a sweep is running -- the reset is what makes the next
    // episode start with an anchor scan instead of inheriting the last one's elapsed clock.
    if (!sweeping)
    {
      monitor_armed_[idx] = false;
      return false;
    }

    const auto now = std::chrono::steady_clock::now();

    // EPISODE ANCHOR: this level has just started sweeping, so emit immediately and start the clock.
    // Deliberately independent of interval_ms -- a sweep shorter than one interval still gets the
    // one reading that says what the source looked like while it ran.
    if (!monitor_armed_[idx])
    {
      monitor_armed_[idx] = true;
      last_monitor_time_[idx] = now;
      return true;
    }

    const auto elapsed =
      std::chrono::duration_cast<std::chrono::milliseconds>(now - last_monitor_time_[idx]).count();
    if (static_cast<double>(elapsed) > interval_ms)
    {
      last_monitor_time_[idx] = now;
      return true;
    }
    return false;
  }

  void ScanCommandQueue::backdateMonitorClock(int idx, double ms)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    if (idx < 0 || idx > 1) return;
    last_monitor_time_[idx] -= std::chrono::milliseconds(static_cast<long long>(ms));
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

  /// Recompute the log-only energy mirrors from the stages that own the truth.
  ///
  /// hcd_energy / hcd_energy_s1 duplicate stages[0]/stages[1].collision_energy purely so the log
  /// writers can emit them; nothing else in the engine reads them, and ScanFactory builds the
  /// instrument request from stages[] alone. Because they are a duplicate, they can go stale --
  /// buildMS3 refreshed stages[0].collision_energy from the tracker's per-ion stage0_params but
  /// left the mirror pointing at the MS2 context's energy.
  ///
  /// This runs in the BUILDERS, not only in push(), and the distinction is load-bearing: two
  /// independent readers hold copies of a command. push() takes its argument by value, so
  /// normalising there reaches the queued copy (which scan_commands.tsv reads at dequeue) but not
  /// the caller's copy -- and FLASHIda.cpp keeps exactly such a copy of every MS2 command to feed
  /// writeIDALogEntry. Setting the mirror at construction precedes both readers, so no command
  /// ever exists with a stale value.
  ///
  /// lround rather than a truncating cast mirrors ScanFactory's (int)Math.Round, so the logged
  /// value equals the value actually sent to the instrument.
  static void syncEnergyMirrors_(ScanCommand& cmd)
  {
    cmd.hcd_energy    = cmd.num_stages > 0 ? static_cast<int32_t>(std::lround(cmd.stages[0].collision_energy)) : 0;
    cmd.hcd_energy_s1 = cmd.num_stages > 1 ? static_cast<int32_t>(std::lround(cmd.stages[1].collision_energy)) : 0;
  }

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
    cmd.faims_enabled = config_.faims().enabled ? 1 : 0;
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

    // The AGC prescan is split by physics, not by "whatever makeMS1 copies" (ADR-0011).
    //
    // SOURCE-REGION comes from the survey's config: these are upstream of the analyzer and decide
    // WHICH ions arrive, so a flux estimate taken through a different source configuration does not
    // describe the scans it is used to gain-correct. ScanFactory now emits this group
    // unconditionally, so leaving them at 0 here would actively command RF lens 0 rather than
    // omitting the key -- this is load-bearing, not cosmetic. Pre-port set exactly these three on
    // the AGC scan from MS1 (Flash.cs@cd0d086:282-284).
    cmd.rf_lens = config_.level(1).scans[0].rf_lens;
    cmd.source_cid = config_.level(1).scans[0].source_cid;
    cmd.source_cid_scaling = config_.level(1).scans[0].source_cid_scaling;
    cmd.faims_enabled = config_.faims().enabled ? 1 : 0;

    // ANALYZER-SIDE stays hardcoded, joining agc_target/max_it/analyzer/scan_rate above: this is a
    // fast ion-trap prescan, not a small copy of the Orbitrap survey. microscans in particular must
    // NOT come from config -- the shipped method.json sets ms1.microscans = 4, which would quadruple
    // a priority-0 scan with a 1 ms max IT. Both values are pre-port (Flash.cs@cd0d086:279,285).
    cmd.microscans = 1;
    std::strncpy(cmd.data_type, "Profile", sizeof(cmd.data_type) - 1);
    cmd.data_type[sizeof(cmd.data_type) - 1] = '\0';

    std::strncpy(cmd.scan_rate, "Turbo", sizeof(cmd.scan_rate) - 1);
    std::strncpy(cmd.analyzer, "IonTrap", sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    std::strncpy(cmd.scan_description, "A", sizeof(cmd.scan_description) - 1);
    cmd.scan_description[sizeof(cmd.scan_description) - 1] = '\0';

    return cmd;
  }

  ScanCommand ScanCommandQueue::buildMS2(const PeakGroup& pg, int charge, const ScanConfig& scan_config, int priority, int parent_scan_id,
                                         const std::vector<int>* acquisition_set, char marker)
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

    // Scan properties from the PASSED scan_config, like every other field in this function.
    // These two used to read config_.level(2).scans[0] directly, which silently overrode the
    // argument: FLASHIda.cpp loops over every level-2 scan config, so ms_settings.ms2[1..N]
    // acquired at ms2[0]'s AGC target and max IT -- and an exploration override of either key
    // was inert at level 2. buildMS3 and buildFollowUp never had this; buildMS2 was the outlier.
    // A ScanConfig fully determines its scan's instrument parameters (ADR-0009).
    cmd.agc_target = scan_config.agc_target;
    cmd.max_it = scan_config.max_it;

    // New scan parameters from MS2 config
    cmd.microscans = scan_config.microscans;
    cmd.rf_lens = scan_config.rf_lens;
    cmd.source_cid = scan_config.source_cid;
    cmd.source_cid_scaling = scan_config.source_cid_scaling;
    cmd.faims_enabled = config_.faims().enabled ? 1 : 0;
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

    // Co-isolation notches: the other charge states of THIS SAME PeakGroup, when the acquisition mode
    // asks for them (ADR-0016). Geometry comes from getMzRange(z) per charge -- the same measured
    // window the anchor uses, never derived from theory.
    //
    // WHICH charges is the caller's decision, not this function's: @p acquisition_set is the set
    // admission recorded (anchor first, then the rest in charge_selection order, already gated,
    // authored-filtered, sibling-reduced and capped at max_charge_states). This walk is the one that
    // reaches the instrument, so it writes exactly that set. It used to re-run selectNotches here with
    // MAX_NOTCHES_PER_STAGE, so max_charge_states never reached the wire for an unrestricted row
    // (ADR-0059: the quant_multiplexed golden, cap 3, isolated ten charges) -- the wire disagreeing
    // with every record of itself, as an unfiltered authored set once did (acquired=16,13,10 logged
    // while the command carried all ten windows).
    //
    // With no set supplied (a direct caller with no admission behind it) the same policy is applied
    // here from config, so there is still exactly one ranking and one cap.
    if (config_.targeting().precursor_charges == ChargeAcquisitionMode::Multiplexed)
    {
      std::vector<NotchCandidate> notch_candidates = peakGroupNotchCandidates(pg, optimal_window_margin_);
      std::vector<NotchCandidate> notches;
      if (acquisition_set != nullptr && !acquisition_set->empty())
      {
        for (int c : *acquisition_set)
        {
          if (c == charge) continue;   // the anchor is the stage itself
          for (const NotchCandidate& n : notch_candidates)
            if (n.charge == c) { notches.push_back(n); break; }
        }
      }
      else
      {
        const int max_charge_states = config_.targeting().max_charge_states;
        const int max_notches = max_charge_states > 0 ? std::min(max_charge_states - 1, MAX_NOTCHES_PER_STAGE)
                                                      : MAX_NOTCHES_PER_STAGE;
        notches = selectNotches(notch_candidates, charge, config_.targeting().snr_threshold, max_notches,
                                config_.targeting().charge_selection, "MS2 z=" + std::to_string(charge));
      }
      writeNotchesForStage(cmd, 0, notches);
    }

    // Scan description: {3-char ID}R{mass_kDa}k@{charge}  (adaptive-precision mass token)
    std::string id_str = encode(id);
    std::string mass_tok = formatMassToken(pg.getMonoMass() / 1000.0, charge, '\0', 0);
    // `marker` is 'R' for an identification MS2 and 'Q' for the quantification scan (ADR-0038).
    // Both are built by this one function because they are the same KIND of command -- a rostered
    // MS2 of a selected precursor -- differing only in role, which is exactly what the marker says.
    std::snprintf(cmd.scan_description, 16, "%s%c%s@%d", id_str.c_str(), marker, mass_tok.c_str(), charge);

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
    cmd.pad2 = 0;   // hcd_energy is derived from stages[0] by syncEnergyMirrors_ below

    // Parent tracking ID. `> 0` is exactly right, but only because nextTrackingIdInt_ never issues
    // 0 -- it is the reserved "no parent / root" sentinel FLASHIda.cpp passes literally for a root
    // MS2. That reservation is the whole contract: while the counter DID issue 0, a survey holding
    // it produced MS2s with an empty parent_tracking_id, indistinguishable from having no parent.
    // The AGC prescan hid this by absorbing 0 on every fresh engine until ADR-0031 deleted it.
    //
    // buildMS3 and buildFollowUp encode their parent unconditionally, which reads as inconsistent
    // and is not: with 0 unreachable, every parent they are handed is a real id, so a guard there
    // could never fire. Do not "harmonise" by removing this one -- it is what keeps a root MS2's
    // parent empty rather than "!!!".
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

    syncEnergyMirrors_(cmd);
    return cmd;
  }

  ScanCommand ScanCommandQueue::buildMS3(const ScanCommand& ms2_ctx, const ScanConfig& ms3_config,
                                          double frag_mz, int frag_charge, double iso_width, int parent_scan_id,
                                          char ion_type, int frag_index, int priority,
                                          const FragmentAnalysis::FragmentScores& frag_scores,
                                          const Ms2Params* stage0_params,
                                          const MS3FragmentMatcher::ProteoformContext& proto_ctx,
                                          const std::vector<NotchCandidate>* stage1_notches)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    ScanCommand cmd{};
    int id = nextTrackingIdInt_();
    cmd.scan_id = id;
    // Render the wide MS3 fragment ProForma here (moved off the callers) and stash it for the
    // scan_commands.tsv ms3_proteoform column (drained by takeMS3Proteoform). Inline write —
    // queue_mutex_ is already held and is non-recursive.
    const std::string ms3_proteoform = MS3FragmentMatcher::fragmentProForma(
        config_.characterization().protein_sequence, proto_ctx, ion_type, frag_index);
    if (!ms3_proteoform.empty()) ms3_cmd_proteoform_[id] = ms3_proteoform;
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
    cmd.faims_enabled = config_.faims().enabled ? 1 : 0;
    std::strncpy(cmd.data_type, ms3_config.data_type.c_str(), sizeof(cmd.data_type) - 1);
    cmd.data_type[sizeof(cmd.data_type) - 1] = '\0';
    std::strncpy(cmd.scan_rate, ms3_config.scan_rate.c_str(), sizeof(cmd.scan_rate) - 1);
    cmd.scan_rate[sizeof(cmd.scan_rate) - 1] = '\0';

    // Parent tracking ID (mandatory param; caller passes the immediate parent, normally ms2_ctx.scan_id)
    std::string parent_enc = encode(parent_scan_id);
    std::strncpy(cmd.parent_scan_id, parent_enc.c_str(), 3);
    cmd.parent_scan_id[3] = '\0';

    // Stage 0: MS2 precursor (from MS2 context)
    cmd.stages[0] = ms2_ctx.stages[0];

    // Stage 0 override: if caller supplies per-ion optimised MS2 params, apply them now.
    // Precursor mz/isolation_width/charge_state are NOT overridden — those come from ms2_ctx
    // and describe the physical isolation window, which is fixed by the original MS2.
    if (stage0_params != nullptr)
    {
      cmd.stages[0].collision_energy = stage0_params->collision_energy;
      std::strncpy(cmd.stages[0].activation_type, stage0_params->activation_type.c_str(),
                   sizeof(cmd.stages[0].activation_type) - 1);
      cmd.stages[0].activation_type[sizeof(cmd.stages[0].activation_type) - 1] = '\0';
      cmd.stages[0].reaction_time = stage0_params->reaction_time;
    }

    // Stage 0's notches are INHERITED from the MS2 context, not re-selected (ADR-0016 decision 6):
    // if the MS2 co-isolated a charge set, the replay that regenerates the fragment must isolate the
    // same set, or it produces a fraction of the fragment signal the MS3 then tries to work with.
    // No SNR gate and no re-ranking here -- these already passed both when the MS2 was built, and
    // re-doing either would need the MS1 PeakGroup, which this builder does not have.
    //
    // Stage 0 and stage 1 own DISJOINT notch blocks, so this can never starve the fragment stage and
    // write order does not matter. It used to: both stages drew from one pool of
    // MAX_ISOLATION_STAGES - num_stages == 8 slots, stage 0 was written first and took all of them
    // whenever the parent MS2 was a full multiplex, and stage 1 -- the whole point of an MS3 -- was
    // left with zero (observed as notches/stage=[8,0] in the multiplexed_ms2 golden mode).
    {
      auto inherited_range = notchesForStage(ms2_ctx, 0);
      if (inherited_range.second > 0)
      {
        std::vector<NotchCandidate> inherited;
        inherited.reserve(inherited_range.second);
        for (int i = 0; i < inherited_range.second; ++i)
          inherited.push_back({inherited_range.first[i].charge_state,
                               inherited_range.first[i].precursor_mz,
                               inherited_range.first[i].isolation_width, 0.0, 0.0});
        writeNotchesForStage(cmd, 0, inherited);
      }
    }

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

    // Stage 1's notches: the other charge states of this fragment, chosen by the tracker, which has
    // already gated and ranked them. This stage's own MAX_NOTCHES_PER_STAGE block, independent of
    // stage 0's, so truncation here means the tracker handed over more than a 10-plex -- which its
    // own selectNotches clamp should have prevented. Report it rather than trimming in silence.
    if (stage1_notches != nullptr && !stage1_notches->empty())
    {
      const int wrote = writeNotchesForStage(cmd, 1, *stage1_notches);
      if (wrote < static_cast<int>(stage1_notches->size()))
        std::cout << "[NOTCH-CLAMP] MS3 stage1 kept=" << wrote
                  << " dropped=" << (static_cast<int>(stage1_notches->size()) - wrote)
                  << " (over this stage's own cap of " << MAX_NOTCHES_PER_STAGE << ")" << std::endl;
    }

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
    // hcd_energy / hcd_energy_s1 are derived from stages[0]/stages[1] below. Copying
    // ms2_ctx.hcd_energy here was the stale-mirror bug: the stage0_params override above
    // refreshes stages[0].collision_energy but left this pointing at the MS2's energy.

    syncEnergyMirrors_(cmd);
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

    // A ScanConfig fully determines the scan's instrument parameters (ADR-0009), so EVERY one of
    // them comes from follow_up_config below. `cmd = ctx` above is deliberate and must stay: it
    // carries what makes this a follow-up *of that precursor* -- the targeting (mono_mass,
    // precursor_mz, isolation_width, charge_state), the precursor scoring fields, and faims_cv.
    // Rebuilding from ScanCommand{} instead would zero faims_cv, so the follow-up would be acquired
    // at a different compensation voltage and sample a different ion population entirely.
    std::strncpy(cmd.analyzer, follow_up_config.analyzer.c_str(), sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    cmd.orbitrap_resolution = follow_up_config.resolution;
    cmd.agc_target = follow_up_config.agc_target;
    cmd.max_it = follow_up_config.max_it;
    cmd.first_mass = follow_up_config.first_mass;
    cmd.last_mass = follow_up_config.last_mass;
    cmd.microscans = follow_up_config.microscans;
    cmd.rf_lens = follow_up_config.rf_lens;
    cmd.source_cid = follow_up_config.source_cid;
    cmd.source_cid_scaling = follow_up_config.source_cid_scaling;
    cmd.faims_enabled = config_.faims().enabled ? 1 : 0;
    std::strncpy(cmd.data_type, follow_up_config.data_type.c_str(), sizeof(cmd.data_type) - 1);
    cmd.data_type[sizeof(cmd.data_type) - 1] = '\0';
    std::strncpy(cmd.scan_rate, follow_up_config.scan_rate.c_str(), sizeof(cmd.scan_rate) - 1);
    cmd.scan_rate[sizeof(cmd.scan_rate) - 1] = '\0';
    // hcd_energy is derived from stages[0] by syncEnergyMirrors_ below; the ctx copy this started from no longer
    // leaks through, since stages[0].collision_energy below is the follow-up's own energy.

    cmd.stages[0].collision_energy = static_cast<double>(follow_up_config.collision_energy);
    std::strncpy(cmd.stages[0].activation_type, follow_up_config.activation.c_str(),
                 sizeof(cmd.stages[0].activation_type) - 1);
    cmd.stages[0].activation_type[sizeof(cmd.stages[0].activation_type) - 1] = '\0';
    // Activation-coupled parameters MUST travel with the activation being overridden above.
    // Leaving them inherited is how an HCD follow-up ended up carrying the parent ETD's ion-ion
    // reaction settings (reaction_time/reagent_*), which reached the instrument via ScanFactory.
    cmd.stages[0].reaction_time = follow_up_config.reaction_time;
    cmd.stages[0].reagent_max_it = follow_up_config.reagent_max_it;
    cmd.stages[0].reagent_agc_target = follow_up_config.reagent_agc_target;

    std::string id_str = encode(cmd.scan_id);
    std::string mass_tok = formatMassToken(cmd.mono_mass / 1000.0, cmd.stages[0].charge_state, '\0', 0);
    std::snprintf(cmd.scan_description, 16, "%s%c%s@%d",
                  id_str.c_str(), suffix, mass_tok.c_str(), cmd.stages[0].charge_state);

    std::cout << "[TRACK-CREATE] id=" << id_str
              << " ms_level=2 type=followup_" << suffix
              << std::endl;

    syncEnergyMirrors_(cmd);
    return cmd;
  }

  // --- Queue operations ---

  void ScanCommandQueue::push(ScanCommand cmd)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    push_(cmd);
  }

  void ScanCommandQueue::push_(ScanCommand cmd)
  {
    // Backstop only -- the builders already did this at construction, which is what makes the
    // value correct for every reader rather than just for the queue (see syncEnergyMirrors_).
    // Same helper on the same input, so it cannot disagree with them; it exists so a future
    // command that reaches the queue without going through a builder still holds the invariant.
    syncEnergyMirrors_(cmd);

    // Causal invariant of survey re-planning: a command is only ever derived from one that went out. A parent
    // still registered in the plan registry was never dequeued, so its child could survive the parent's drop.
    const char* const parent_begin = cmd.parent_scan_id;
    const std::string parent(parent_begin, std::find(parent_begin, parent_begin + sizeof(cmd.parent_scan_id), '\0'));
    if (!parent.empty() && unit_of_.count(decode(parent)) > 0)
    {
      std::cout << "[REPLAN-INVARIANT] child=" << encode(cmd.scan_id) << " parent=" << parent << std::endl;
    }

    cmd.enqueue_timestamp_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());
    int p = std::clamp(cmd.priority, 0, 3);
    queues_[p].push_back(cmd);
  }

  std::optional<ScanCommand> ScanCommandQueue::dequeue()
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    return dequeue_();
  }

  std::optional<ScanCommand> ScanCommandQueue::dequeue_()
  {
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

        // Commit at dequeue: the first command of a unit to go out commits its species, and each command
        // commits its own charge's stamps. Both move to the inbox the analysis thread folds, in commit order.
        auto u = unit_of_.find(cmd.scan_id);
        if (u != unit_of_.end())
        {
          const int unit_id = u->second;
          PlanUnit& unit = units_.at(unit_id);
          const int z = cmd.num_stages > 0 ? std::abs(cmd.stages[0].charge_state) : 0;
          if (!unit.committed)
          {
            unit.committed = true;
            for (const auto& s : unit.stamps.species_level)
            {
              inbox_.push_back(s);
              committed_since_fold_.push_back(s);
            }
          }
          auto pc = unit.stamps.per_charge.find(z);
          if (pc != unit.stamps.per_charge.end())
          {
            for (const auto& s : pc->second)
            {
              inbox_.push_back(s);
              committed_since_fold_.push_back(s);
            }
            unit.stamps.per_charge.erase(pc);
          }
          unit.sent_charges.push_back(z);
          unit_of_.erase(u);
          if (--unit.commands_left == 0) units_.erase(unit_id);
        }
        adjustable_.erase(cmd.scan_id);
        return cmd;
      }
    }
    return std::nullopt;
  }

  void ScanCommandQueue::forgetQueued_(int scan_id)
  {
    adjustable_.erase(scan_id);
    decision_info_.erase(scan_id);
    auto u = unit_of_.find(scan_id);
    if (u == unit_of_.end()) return;
    auto it = units_.find(u->second);
    unit_of_.erase(u);
    if (it != units_.end() && --it->second.commands_left == 0) units_.erase(it); // an uncommitted unit's stamps go with it
  }

  std::optional<ScanCommand> ScanCommandQueue::pushAndDequeue(ScanCommand cmd)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    push_(cmd);
    // Cannot be nullopt: push_ clamps into [0,3] and we hold the lock across both, so the queue we
    // just wrote to is still non-empty here. That is the whole point -- see the header.
    return dequeue_();
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
          ms3_cmd_proteoform_.erase(it->scan_id);  // drop any stashed MS3 proteoform (no-op for non-MS3)
          forgetQueued_(it->scan_id);              // never acquired: leaves the plan registry (no-op if unregistered)
          it = queues_[p].erase(it);
        }
        else
        {
          ++it;
        }
      }
    }
  }

  std::string ScanCommandQueue::takeMS3Proteoform(int scan_id)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    auto it = ms3_cmd_proteoform_.find(scan_id);
    if (it == ms3_cmd_proteoform_.end()) return std::string();
    std::string pf = std::move(it->second);
    ms3_cmd_proteoform_.erase(it);
    return pf;
  }

  // --- Survey re-planning: the plan registry ---

  namespace
  {
    /// FAIMS CVs are authored values copied, never computed, so equality is exact up to representation.
    bool sameCv(double a, double b)
    {
      return std::abs(a - b) < 1e-6;
    }

    /// The admission key of a unit's anchor (its first command), as selection keyed the species: the nominal mass,
    /// the integer centre m/z buildMS2 isolated, the anchor charge, the score selection recorded (the unit's own
    /// MassQscore stamp, else the command's qscore) and whether an authored charge set governs it.
    UnitAdmissionKey anchorKey(const PlanUnit& u)
    {
      const ScanCommand& a = u.commands.front();
      UnitAdmissionKey k;
      k.nominal_mass = u.nominal_mass;
      k.anchor_charge = a.num_stages > 0 ? std::abs(a.stages[0].charge_state) : 0;
      k.integer_mz = static_cast<int>(std::round(a.stages[0].precursor_mz));
      k.score = a.qscore;
      for (const auto& s : u.stamps.species_level)
      {
        if (s.kind == MemoryStamp::Kind::MassQscore) k.score = s.value;
      }
      for (const auto& [z, stamps] : u.stamps.per_charge)
      {
        for (const auto& s : stamps)
        {
          if (s.kind == MemoryStamp::Kind::AuthoredAcquired) k.has_authored = true;
        }
      }
      return k;
    }
  } // namespace

  ReplanResult ScanCommandQueue::replacePlan(double cv, std::vector<PlanUnit> units, const std::vector<AdjustRequest>& adjust,
                                             const std::vector<int>& replace_remainder_of)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    ReplanResult r;

    // 1. Which units at this CV lose their queued commands: every uncommitted one, and the committed ones whose
    //    remainder the new plan re-derives. Every other committed unit is kept as it is.
    std::set<int> leaving;
    for (const auto& [id, u] : units_)
    {
      if (!sameCv(u.faims_cv, cv)) continue;
      const bool remainder = std::find(replace_remainder_of.begin(), replace_remainder_of.end(), id) != replace_remainder_of.end();
      if (!u.committed || remainder)
      {
        leaving.insert(id);
        if (!u.committed)
        {
          r.dropped_quant_groups.insert(r.dropped_quant_groups.end(), u.quant_group_ids.begin(), u.quant_group_ids.end());
          r.dropped_sweep_groups.insert(r.dropped_sweep_groups.end(), u.sweep_group_ids.begin(), u.sweep_group_ids.end());
        }
      }
      else
      {
        ++r.kept;
      }
    }
    if (!leaving.empty())
    {
      for (int p = 0; p < 4; ++p)
      {
        auto it = queues_[p].begin();
        while (it != queues_[p].end())
        {
          auto u = unit_of_.find(it->scan_id);
          if (u != unit_of_.end() && leaving.count(u->second) > 0)
          {
            r.dropped_scan_ids.push_back(it->scan_id);
            decision_info_.erase(it->scan_id);
            ms3_cmd_proteoform_.erase(it->scan_id);
            unit_of_.erase(u);
            it = queues_[p].erase(it);
          }
          else
          {
            ++it;
          }
        }
      }
      for (int id : leaving) units_.erase(id); // an uncommitted unit's stamps go with it: never acquired, never excluded
    }

    // 2. Insert the new plan. A unit whose species committed since the last fold would have been refused had the
    //    survey seen that commit; it is reconciled away instead. A unit arriving committed is a remainder.
    for (auto& unit : units)
    {
      if (unit.unit_id <= 0 || units_.count(unit.unit_id) > 0)
      {
        do
        {
          ++unit_id_counter_;
        } while (units_.count(unit_id_counter_) > 0);
        unit.unit_id = unit_id_counter_;
      }
      if (unit.commands.empty()) continue;
      if (!unit.committed && wouldBar(committed_since_fold_, anchorKey(unit)))
      {
        r.reconciled_unit_ids.push_back(unit.unit_id);
        continue;
      }
      unit.faims_cv = cv;
      unit.commands_left = static_cast<int>(unit.commands.size());
      for (size_t k = 0; k < unit.commands.size(); ++k)
      {
        const ScanCommand& cmd = unit.commands[k];
        decision_info_[cmd.scan_id] = DecisionInfo{unit.decision_survey_id, k < unit.purity.size() ? unit.purity[k] : -1.0};
        push_(cmd);
        unit_of_[cmd.scan_id] = unit.unit_id;
      }
      ++r.added;
      const int id = unit.unit_id;
      units_[id] = std::move(unit);
    }

    // 3. Rewrite the adjustable commands in place: same scan id, same queue position, a fresh enqueue time.
    const auto now_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now().time_since_epoch()).count());
    for (const auto& a : adjust)
    {
      if (adjustable_.count(a.scan_id) == 0) continue;
      for (int p = 0; p < 4; ++p)
      {
        auto it = std::find_if(queues_[p].begin(), queues_[p].end(), [&a](const ScanCommand& c) { return c.scan_id == a.scan_id; });
        if (it == queues_[p].end()) continue;
        ScanCommand cmd = a.rewritten;
        cmd.scan_id = a.scan_id;
        cmd.priority = it->priority;
        syncEnergyMirrors_(cmd);
        cmd.enqueue_timestamp_ms = now_ms;
        *it = cmd;
        decision_info_[a.scan_id] = DecisionInfo{a.decision_survey_id, a.isolation_purity};
        ++r.adjusted;
        break;
      }
    }
    return r;
  }

  void ScanCommandQueue::registerAdjustable(int scan_id, std::vector<int> resolved_charges, int q_survey_id)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    adjustable_[scan_id] = AdjustableR{std::move(resolved_charges), q_survey_id};
    decision_info_[scan_id] = DecisionInfo{q_survey_id, -1.0};
  }

  ScanCommandQueue::Waiting ScanCommandQueue::waitingAt(double cv) const
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    Waiting w;
    for (const auto& [id, u] : units_)
    {
      if (!u.committed || !sameCv(u.faims_cv, cv)) continue;
      PlanUnit copy = u;
      copy.commands.clear();
      copy.purity.clear();
      for (size_t k = 0; k < u.commands.size(); ++k)
      {
        if (unit_of_.count(u.commands[k].scan_id) == 0) continue; // already sent
        copy.commands.push_back(u.commands[k]);
        copy.purity.push_back(k < u.purity.size() ? u.purity[k] : -1.0);
      }
      if (!copy.commands.empty()) w.committed_units.push_back(std::move(copy));
    }
    for (int p = 0; p < 4; ++p)
    {
      for (const auto& c : queues_[p])
      {
        auto a = adjustable_.find(c.scan_id);
        if (a != adjustable_.end()) w.adjustable.emplace_back(c, a->second.resolved_charges);
      }
    }
    return w;
  }

  std::vector<MemoryStamp> ScanCommandQueue::takeCommitted()
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    std::vector<MemoryStamp> out = std::move(inbox_);
    inbox_.clear();
    committed_since_fold_.clear();
    return out;
  }

  std::optional<DecisionInfo> ScanCommandQueue::takeDecisionInfo(int scan_id)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    auto it = decision_info_.find(scan_id);
    if (it == decision_info_.end()) return std::nullopt;
    DecisionInfo info = it->second;
    decision_info_.erase(it);
    return info;
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
          ms3_cmd_proteoform_.erase(it->scan_id);  // drop any stashed MS3 proteoform (no-op for non-MS3)
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
        ms3_cmd_proteoform_.erase(pit->first);  // drop any stashed MS3 proteoform (no-op for non-MS3)
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
