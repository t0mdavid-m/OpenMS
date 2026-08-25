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
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <cstdint>
#include <cstring>
#include <utility>

namespace OpenMS
{
  /// Maximum number of cascade (MSⁿ) isolation stages per scan command.
  ///
  /// This is the instrument's documented ';'-axis limit: PrecursorMass, IsolationWidth, ChargeStates,
  /// ActivationType and CollisionEnergy each accept "a maximum of 10 values" as ';'-separated groups.
  /// The engine only ever builds 0, 1 or 2 stages.
  static constexpr int MAX_ISOLATION_STAGES = 10;

  /// Maximum co-isolation windows the instrument accepts per fragmentation stage, and hence the
  /// maximum notches per stage (a 10-plex is the anchor plus 9 notches).
  ///
  /// 10 is MSXTargets' own limit — "AGC target values for MSX windows, in m/z order … a maximum of
  /// 10 values" — which is the per-scan MSX window count, matching the Q Exactive method editor's
  /// "MSX count … 1 (no spectral multiplexing) to 10 fillings". It is a DIFFERENT ten from
  /// MAX_ISOLATION_STAGES: that one counts cascade stages on the ';' axis, this one counts parallel
  /// windows on the ',' axis. Conflating them capped an MS3's fragment stage at 10 - num_stages and,
  /// worse, made stage 0 and stage 1 compete for one shared pool.
  static constexpr int MAX_NOTCHES_PER_STAGE = 9;

  /// Notch capacity of a scan command: MAX_NOTCHES_PER_STAGE for each of the two stages that can
  /// carry notches. Stage k owns the FIXED block [k * MAX_NOTCHES_PER_STAGE, +MAX_NOTCHES_PER_STAGE),
  /// so the two stages never contend and write order is irrelevant.
  static constexpr int MAX_NOTCHES = 2 * MAX_NOTCHES_PER_STAGE;

  /// One co-isolation notch: a parallel isolation window inside a single fragmentation stage.
  ///
  /// Geometry only, and that is structural rather than conventional: every notch of a stage fires
  /// into the SAME fragmentation event, so there is no per-notch collision energy, activation type or
  /// reaction time to carry — the wire has exactly one of each per ';'-group. Scores are likewise
  /// absent; the *_s1 scoring block is two-stage and reports the anchor's (ADR-0016).
  ///
  /// Layout: 2 doubles (16) + 2 int32 (8) = 24 bytes.
  struct OPENMS_DLLAPI Notch
  {
    double precursor_mz;    ///< Centre m/z of this window
    double isolation_width; ///< Full width of this window (Da), measured, never derived from theory
    int32_t charge_state;   ///< The charge state this window isolates
    int32_t pad_;           ///< Explicit padding; Notch is 8-aligned because of the doubles
  };
  static_assert(sizeof(Notch) == 24, "Notch must be 24 bytes for P/Invoke");

  /// Blittable struct representing one isolation stage in a multi-stage MS2/MS3 scan.
  /// Layout: 5 doubles (40) + 2 int32 (8) + char[32] (32) = 80 bytes.
  struct OPENMS_DLLAPI IsolationStage
  {
    double precursor_mz;         ///< Center m/z for isolation window
    double isolation_width;      ///< Full width of isolation window (Da)
    double collision_energy;     ///< Normalized collision energy (%)
    /// Ion/ion reaction time (ms). 0 = not used, and that is an ABI sentinel ONLY -- it means "this
    /// stage has no ion-ion reaction", never "react for zero milliseconds". The instrument REJECTS a
    /// reaction time of 0 on an ETD-family stage, so an un-fragmented ETD reference is commanded at
    /// Config.h's MIN_REACTION_TIME_MS (0.03), not at 0. Whether the key is emitted at all is decided
    /// by the stage's ACTIVATION rather than by this value (ADR-0030).
    double reaction_time;
    double reagent_max_it;       ///< Reagent max injection time (ms), 0 = not used
    int32_t reagent_agc_target;  ///< Reagent AGC target, 0 = not used
    int32_t charge_state;        ///< Precursor charge state (0 = unknown)
    char activation_type[32];    ///< Activation method (e.g., "HCD", "ETD", "EThcD")
  };
  static_assert(sizeof(IsolationStage) == 80, "IsolationStage must be 80 bytes for P/Invoke");

  /// Blittable struct representing a complete scan command for the instrument.
  /// Layout: 1248 (existing) + 8 (dequeue_timestamp_ms) + 8 (microscans+pad3) + 24 (rf_lens+source_cid+source_cid_scaling)
  ///       + 64 (data_type+scan_rate) + 4 (parent_scan_id) + 84 (stage-1 scoring) + 608 (reserved) = 2048.
  /// The trailing 608 covers window_snr (8 @1440) + faims_enabled (4 @1448) + stage0_notch_count
  /// (4 @1452) + stage1_notch_count (4 @1456) + pad4 (4 @1460) + notches (432 @1464)
  /// + reserved_ (152 @1896):
  /// new fields are carved OUT of reserved_ so the 2048 total and every existing offset stay fixed.
  struct OPENMS_DLLAPI ScanCommand
  {
    int32_t scan_id;             ///< Unique tracking ID (encoded as 3-char string in scan description)
    int32_t msn_level;           ///< MS level: 1 = MS1, 2 = MS2, 3 = MS3
    int32_t priority;            ///< Queue priority: 0 = highest, 3 = lowest
    int32_t is_agc;              ///< 1 if this is an AGC calibration scan, 0 otherwise
    int32_t num_stages;          ///< Number of valid isolation stages (0 for MS1)
    int32_t orbitrap_resolution; ///< Orbitrap resolution (e.g., 120000)
    int32_t agc_target;          ///< AGC target value
    int32_t pad1;                ///< Padding for 8-byte alignment
    double first_mass;           ///< Scan range start (m/z)
    double last_mass;            ///< Scan range end (m/z)
    double max_it;               ///< Maximum injection time (ms)
    char analyzer[32];           ///< Mass analyzer name (e.g., "Orbitrap", "IonTrap")
    char scan_description[256];  ///< Human-readable description (includes tracking ID)
    IsolationStage stages[MAX_ISOLATION_STAGES]; ///< Isolation stages array
    uint64_t enqueue_timestamp_ms;  ///< Timestamp when command was enqueued (steady_clock ms)
    uint64_t dequeue_timestamp_ms;  ///< Timestamp when command was dequeued/sent to instrument (steady_clock ms)

    // Precursor scoring data (populated by buildMS2Command_ for diagnostic output)
    double qscore;                  ///< Quality score from PeakGroup::getQscore()
    double mono_mass;               ///< Monoisotopic mass from PeakGroup::getMonoMass()
    double charge_cos;              ///< Charge isotope cosine from PeakGroup::getChargeIsotopeCosine(charge)
    double charge_snr;              ///< Charge SNR from PeakGroup::getChargeSNR(charge)
    double iso_cos;                 ///< Isotope cosine from PeakGroup::getIsotopeCosine()
    double snr;                     ///< Total SNR from PeakGroup::getSNR()
    double charge_score;            ///< Charge score from PeakGroup::getChargeScore()
    double ppm_error;               ///< Avg PPM error from PeakGroup::getAvgPPMError()
    double precursor_intensity;     ///< Per-charge intensity from PeakGroup::getChargeIntensity(charge)
    double peakgroup_intensity;     ///< Total peak group intensity from PeakGroup::getIntensity()
    int32_t hcd_energy;             ///< HCD collision energy
    int32_t pad2;                   ///< Alignment padding
    double faims_cv;                ///< FAIMS compensation voltage (0.0 if non-FAIMS)
    int32_t microscans;            ///< Number of microscans (0 = use method default)
    int32_t pad3;                  ///< Alignment padding
    double rf_lens;                ///< RF lens voltage (0 = use method default)
    double source_cid;             ///< Source CID energy (0 = use method default)
    double source_cid_scaling;     ///< Source CID scaling factor (0 = use method default)
    char data_type[32];            ///< Data type (e.g., "Centroid", "Profile"; empty = method default)
    char scan_rate[32];            ///< Scan rate (e.g., "Normal", "Turbo"; empty = method default)
    char parent_scan_id[4];        ///< Parent scan's encoded tracking ID (3 chars + null; empty for MS1)

    // --- Stage-1 (MS3 fragment) precursor scoring. Stage-0 lives in the top-level scalars above
    //     (qscore..hcd_energy); for MS3 commands those carry the MS2 precursor, these carry the fragment.
    //     int32 is FIRST so it absorbs the 4-aligned slot at offset 1356 and the doubles land 8-aligned
    //     at 1360 with ZERO implicit padding (reordering shifts every offset — see ScanCommandLayout tests).
    int32_t hcd_energy_s1;          ///< @1356 MS3 fragment-stage collision energy
    double  mono_mass_s1;           ///< @1360 fragment PeakGroup::getMonoMass()
    double  qscore_s1;              ///< @1368 fragment PeakGroup::getQscore()
    double  charge_cos_s1;          ///< @1376 fragment getChargeIsotopeCosine(frag_charge)
    double  charge_snr_s1;          ///< @1384 fragment getChargeSNR(frag_charge)
    double  iso_cos_s1;             ///< @1392 fragment getIsotopeCosine()
    double  snr_s1;                 ///< @1400 fragment getSNR()
    double  charge_score_s1;        ///< @1408 fragment getChargeScore()
    double  ppm_error_s1;           ///< @1416 fragment getAvgPPMError()
    double  precursor_intensity_s1; ///< @1424 fragment getChargeIntensity(frag_charge)
    double  peakgroup_intensity_s1; ///< @1432 fragment getIntensity()
    double  window_snr = -1.0;     ///< @1440 isolation-window SNR (FragmentAnalysis::windowSnr) of THIS command's
                                   ///<       commanded window; -1.0 = not computed. Set at build time and travels
                                   ///<       with the command (replaces the former ScanCommandQueue window-SNR map).

    /// @1448 1 if FAIMS is in use for this run, 0 otherwise. Carved from reserved_ (ADR-0012).
    ///
    /// faims_cv alone cannot express this: 0.0 is both "no FAIMS" and a legitimate compensation
    /// voltage, so ScanFactory had to infer intent from |cv| > 0.001 and could therefore only ever
    /// say "FAIMS on" -- never "FAIMS off", which is an ACTIVE instruction to the instrument, not
    /// the absence of one. int32_t rather than bool to match is_agc and avoid P/Invoke marshalling
    /// ambiguity.
    int32_t faims_enabled;

    /// @1452 / @1456 Number of EXTRA co-isolation notches at cascade stage 0 / stage 1 (ADR-0017,
    /// amended by ADR-0019: notches have their own array and their own per-stage cap).
    ///
    /// A notch is a parallel isolation window within ONE fragmentation stage, as against a stage,
    /// which is a further fragmentation performed in sequence. All notches of a stage fire into the
    /// same fragmentation event and are read out as one spectrum. Every notch here holds a different
    /// charge state of the SAME species, so a multi-notch spectrum is not chimeric.
    ///
    /// Each count is at most MAX_NOTCHES_PER_STAGE, so either stage can be a full 10-plex
    /// independently of the other. The descriptors live in notches[] below, in fixed per-stage
    /// blocks; use notchesForStage() rather than open-coding the offset.
    ///
    /// num_stages does NOT count notches, and that is load-bearing: processScan's context gate,
    /// syncEnergyMirrors_, Exploration's `si = num_stages - 1` and ScanFactory's clamp all key on
    /// it, and stay correct only while it means cascade depth.
    int32_t stage0_notch_count;
    int32_t stage1_notch_count;

    /// @1460 Explicit padding so notches[] lands 8-aligned at 1464. Without it the compiler inserts
    /// the same 4 bytes implicitly and the C# mirror has nothing to line up against.
    int32_t pad4;

    /// @1464 Co-isolation notch descriptors, 18 * 24 = 432 bytes, carved from reserved_ (588 -> 152).
    ///
    /// Stage k owns the fixed block [k * MAX_NOTCHES_PER_STAGE, + MAX_NOTCHES_PER_STAGE). Fixed
    /// blocks rather than a packed shared pool, deliberately: they let each stage reach a full
    /// 10-plex, they make write order irrelevant, and they remove the failure mode where an MS3
    /// inheriting its parent's precursor charge set left its own fragment stage with zero slots.
    ///
    /// They no longer live in the unused tail of stages[]: that arrangement (ADR-0017) capped the
    /// total at MAX_ISOLATION_STAGES - num_stages, i.e. 8 shared for an MS3, on the mistaken reading
    /// that the instrument's 10-value PrecursorMass limit was a joint budget for stages and notches.
    /// It is not — that 10 counts ';' groups, while the ',' windows have their own limit of 10 per
    /// group (MSXTargets). A notch also needs only 3 of IsolationStage's 8 fields, so a dedicated
    /// 24-byte record buys 18 slots for less than the 800 bytes stages[20] would have cost.
    Notch notches[MAX_NOTCHES];

    char reserved_[152];           ///< Reserved for future fields (consume from here, never change total size)
  };
  static_assert(sizeof(ScanCommand) == 2048, "ScanCommand must be 2048 bytes for P/Invoke");

  /// The co-isolation notches belonging to cascade stage @p k (0 or 1), as {first, count}.
  ///
  /// Stage k's notches occupy the FIXED block starting at k * MAX_NOTCHES_PER_STAGE, so the offset
  /// depends on nothing but k — not on num_stages, not on the other stage's count. Still prefer this
  /// accessor to open-coding `k * MAX_NOTCHES_PER_STAGE`: it is the one place the block layout is
  /// stated, and it applies the bounds guard.
  ///
  /// An empty range is returned as {nullptr, 0}, so a caller can loop without a count check.
  inline std::pair<const Notch*, int> notchesForStage(const ScanCommand& cmd, int k)
  {
    if (k < 0 || k > 1) return {nullptr, 0};
    const int count = (k == 0) ? cmd.stage0_notch_count : cmd.stage1_notch_count;
    if (count <= 0) return {nullptr, 0};
    const int begin = k * MAX_NOTCHES_PER_STAGE;
    if (count > MAX_NOTCHES_PER_STAGE || begin + count > MAX_NOTCHES) return {nullptr, 0};
    return {&cmd.notches[begin], count};
  }

  /// Total isolation windows the instrument is asked for at cascade stage @p k: the stage itself
  /// plus its notches. This is the per-stage element count on the wire, where the stage's values are
  /// ','-joined within the stage's ';'-separated group.
  inline int windowsAtStage(const ScanCommand& cmd, int k)
  {
    if (k < 0 || k >= cmd.num_stages) return 0;
    return 1 + notchesForStage(cmd, k).second;
  }

} // namespace OpenMS
