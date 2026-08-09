// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

// Verifies the ScanCommand / IsolationStage struct layout — the load-bearing 2048-byte
// C++/C# ABI contract. Expected offsets mirror the C# side in
// FlashIDA/src/Flash.Tests/ScanCommandLayoutTests.cs (P3_U03); keep the two in lockstep.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/NotchSelection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include <cstddef>
#include <cstring>
#include <string>
#include <vector>

using namespace OpenMS;

START_TEST(ScanCommandLayout, "$Id$")

START_SECTION([EXTRA] ScanCommand and IsolationStage struct layout)
{
  // Struct sizes (also guarded at compile time by static_assert in ScanCommand.h)
  TEST_EQUAL(sizeof(ScanCommand), 2048)
  TEST_EQUAL(sizeof(IsolationStage), 80)

  // ScanCommand field offsets
  TEST_EQUAL(offsetof(ScanCommand, scan_id), 0)
  TEST_EQUAL(offsetof(ScanCommand, msn_level), 4)
  TEST_EQUAL(offsetof(ScanCommand, priority), 8)
  TEST_EQUAL(offsetof(ScanCommand, is_agc), 12)
  TEST_EQUAL(offsetof(ScanCommand, num_stages), 16)
  TEST_EQUAL(offsetof(ScanCommand, orbitrap_resolution), 20)
  TEST_EQUAL(offsetof(ScanCommand, agc_target), 24)
  TEST_EQUAL(offsetof(ScanCommand, pad1), 28)
  TEST_EQUAL(offsetof(ScanCommand, first_mass), 32)
  TEST_EQUAL(offsetof(ScanCommand, last_mass), 40)
  TEST_EQUAL(offsetof(ScanCommand, max_it), 48)
  TEST_EQUAL(offsetof(ScanCommand, analyzer), 56)
  TEST_EQUAL(offsetof(ScanCommand, scan_description), 88)
  TEST_EQUAL(offsetof(ScanCommand, stages), 344)
  TEST_EQUAL(offsetof(ScanCommand, enqueue_timestamp_ms), 1144)
  TEST_EQUAL(offsetof(ScanCommand, dequeue_timestamp_ms), 1152)
  TEST_EQUAL(offsetof(ScanCommand, qscore), 1160)
  TEST_EQUAL(offsetof(ScanCommand, mono_mass), 1168)
  TEST_EQUAL(offsetof(ScanCommand, charge_cos), 1176)
  TEST_EQUAL(offsetof(ScanCommand, charge_snr), 1184)
  TEST_EQUAL(offsetof(ScanCommand, iso_cos), 1192)
  TEST_EQUAL(offsetof(ScanCommand, snr), 1200)
  TEST_EQUAL(offsetof(ScanCommand, charge_score), 1208)
  TEST_EQUAL(offsetof(ScanCommand, ppm_error), 1216)
  TEST_EQUAL(offsetof(ScanCommand, precursor_intensity), 1224)
  TEST_EQUAL(offsetof(ScanCommand, peakgroup_intensity), 1232)
  TEST_EQUAL(offsetof(ScanCommand, hcd_energy), 1240)
  TEST_EQUAL(offsetof(ScanCommand, pad2), 1244)
  TEST_EQUAL(offsetof(ScanCommand, faims_cv), 1248)
  TEST_EQUAL(offsetof(ScanCommand, microscans), 1256)
  TEST_EQUAL(offsetof(ScanCommand, pad3), 1260)
  TEST_EQUAL(offsetof(ScanCommand, rf_lens), 1264)
  TEST_EQUAL(offsetof(ScanCommand, source_cid), 1272)
  TEST_EQUAL(offsetof(ScanCommand, source_cid_scaling), 1280)
  TEST_EQUAL(offsetof(ScanCommand, data_type), 1288)
  TEST_EQUAL(offsetof(ScanCommand, scan_rate), 1320)
  TEST_EQUAL(offsetof(ScanCommand, parent_scan_id), 1352)
  TEST_EQUAL(offsetof(ScanCommand, hcd_energy_s1), 1356)
  TEST_EQUAL(offsetof(ScanCommand, mono_mass_s1), 1360)
  TEST_EQUAL(offsetof(ScanCommand, qscore_s1), 1368)
  TEST_EQUAL(offsetof(ScanCommand, charge_cos_s1), 1376)
  TEST_EQUAL(offsetof(ScanCommand, charge_snr_s1), 1384)
  TEST_EQUAL(offsetof(ScanCommand, iso_cos_s1), 1392)
  TEST_EQUAL(offsetof(ScanCommand, snr_s1), 1400)
  TEST_EQUAL(offsetof(ScanCommand, charge_score_s1), 1408)
  TEST_EQUAL(offsetof(ScanCommand, ppm_error_s1), 1416)
  TEST_EQUAL(offsetof(ScanCommand, precursor_intensity_s1), 1424)
  TEST_EQUAL(offsetof(ScanCommand, peakgroup_intensity_s1), 1432)
  TEST_EQUAL(offsetof(ScanCommand, window_snr), 1440)
  // Carved out of reserved_, which has now moved 1448 -> 1460 and shrunk 600 -> 588 across two
  // changes (faims_enabled, ADR-0012; then the two notch counts, ADR-0017). Every offset above is
  // unchanged and sizeof stays 2048; that is the whole point of consuming from the tail.
  TEST_EQUAL(offsetof(ScanCommand, faims_enabled), 1448)
  TEST_EQUAL(offsetof(ScanCommand, stage0_notch_count), 1452)
  TEST_EQUAL(offsetof(ScanCommand, stage1_notch_count), 1456)
  TEST_EQUAL(offsetof(ScanCommand, reserved_), 1460)
  // Both notch counts are int32 and reserved_ began 4-mod-8, so the pair lands with no padding.
  // A double carved here would have needed an explicit pad first, like pad1/pad2/pad3 above.
  TEST_EQUAL(sizeof(ScanCommand) - offsetof(ScanCommand, reserved_), 588)

  // IsolationStage field offsets
  TEST_EQUAL(offsetof(IsolationStage, precursor_mz), 0)
  TEST_EQUAL(offsetof(IsolationStage, isolation_width), 8)
  TEST_EQUAL(offsetof(IsolationStage, collision_energy), 16)
  TEST_EQUAL(offsetof(IsolationStage, reaction_time), 24)
  TEST_EQUAL(offsetof(IsolationStage, reagent_max_it), 32)
  TEST_EQUAL(offsetof(IsolationStage, reagent_agc_target), 40)
  TEST_EQUAL(offsetof(IsolationStage, charge_state), 44)
  TEST_EQUAL(offsetof(IsolationStage, activation_type), 48)
}
END_SECTION

// The notch packing (ADR-0017) is the one implicit part of the layout: stage-1's notches begin at
// num_stages + stage0_notch_count, NOT at a fixed slot. Reading them from a fixed slot is the
// mistake notchesForStage() exists to prevent, so it is pinned here rather than left to call sites.
START_SECTION(notchesForStage_packing)
{
  // MS2, two notches at stage 0. stages[0] is the cascade stage; notches follow at [1] and [2].
  ScanCommand ms2 {};
  ms2.num_stages = 1;
  ms2.stage0_notch_count = 2;
  ms2.stages[0].charge_state = 17;   // anchor
  ms2.stages[1].charge_state = 16;   // notch 1
  ms2.stages[2].charge_state = 15;   // notch 2

  auto n0 = notchesForStage(ms2, 0);
  TEST_EQUAL(n0.second, 2)
  TEST_EQUAL(n0.first[0].charge_state, 16)
  TEST_EQUAL(n0.first[1].charge_state, 15)
  TEST_EQUAL(windowsAtStage(ms2, 0), 3)
  // An MS2 has no cascade stage 1, so neither notches nor windows exist there.
  TEST_EQUAL(notchesForStage(ms2, 1).second, 0)
  TEST_EQUAL(windowsAtStage(ms2, 1), 0)

  // MS3 with notches at BOTH stages -- the case a fixed-slot offset gets wrong.
  ScanCommand ms3 {};
  ms3.num_stages = 2;
  ms3.stage0_notch_count = 2;
  ms3.stage1_notch_count = 3;
  ms3.stages[0].charge_state = 17;   // cascade stage 0 anchor
  ms3.stages[1].charge_state = 4;    // cascade stage 1 anchor
  ms3.stages[2].charge_state = 16;   // stage-0 notch 1
  ms3.stages[3].charge_state = 15;   // stage-0 notch 2
  ms3.stages[4].charge_state = 5;    // stage-1 notch 1
  ms3.stages[5].charge_state = 3;    // stage-1 notch 2
  ms3.stages[6].charge_state = 2;    // stage-1 notch 3

  auto s0 = notchesForStage(ms3, 0);
  TEST_EQUAL(s0.second, 2)
  TEST_EQUAL(s0.first[0].charge_state, 16)
  TEST_EQUAL(s0.first[1].charge_state, 15)
  auto s1 = notchesForStage(ms3, 1);
  TEST_EQUAL(s1.second, 3)
  TEST_EQUAL(s1.first[0].charge_state, 5)
  TEST_EQUAL(s1.first[1].charge_state, 3)
  TEST_EQUAL(s1.first[2].charge_state, 2)
  TEST_EQUAL(windowsAtStage(ms3, 0), 3)
  TEST_EQUAL(windowsAtStage(ms3, 1), 4)

  // No notches: an empty range, so a caller may loop without a count check.
  ScanCommand plain {};
  plain.num_stages = 1;
  TEST_EQUAL(notchesForStage(plain, 0).second, 0)
  TEST_EQUAL(notchesForStage(plain, 0).first == nullptr, true)
  TEST_EQUAL(windowsAtStage(plain, 0), 1)

  // Overflow is refused rather than read past stages[]: 1 cascade + 10 notches exceeds 10 slots.
  ScanCommand over {};
  over.num_stages = 1;
  over.stage0_notch_count = 10;
  TEST_EQUAL(notchesForStage(over, 0).second, 0)

  // A stage index outside {0,1} has no notch counter and yields nothing.
  TEST_EQUAL(notchesForStage(ms3, 2).second, 0)
}
END_SECTION

// selectNotches is the one place the co-isolation policy lives, called from three sites, so its
// rules are pinned here rather than inferred from any one caller (ADR-0016).
START_SECTION(selectNotches_policy)
{
  // charge, mz, width, snr. z=17 is the anchor; z=9 is below noise.
  const std::vector<NotchCandidate> cands = {
    {17, 1000.5, 3.2, 8.0},
    {16,  938.2, 3.0, 5.0},
    {15,  883.9, 2.9, 9.0},
    { 9,  512.1, 2.0, 0.4},
  };

  // The anchor is the cascade stage itself, never a notch; the sub-threshold charge is dropped;
  // what survives is ordered by DESCENDING SNR so a clamp keeps the strongest.
  {
    auto out = selectNotches(cands, /*anchor=*/17, /*snr_threshold=*/1.0, /*max_notches=*/9);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].charge, 15)   // snr 9.0
    TEST_EQUAL(out[1].charge, 16)   // snr 5.0
  }

  // The clamp keeps the top N by SNR. (It also reports the drop on stdout; that a value was dropped
  // is never silent, which is the point.)
  {
    auto out = selectNotches(cands, 17, 1.0, 1);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].charge, 15)
  }

  // max_notches 0 -- an MS3 whose stage-0 notches already consumed every slot -- yields nothing
  // rather than overflowing.
  TEST_EQUAL(selectNotches(cands, 17, 1.0, 0).size(), 0)

  // A higher threshold gates more out; a threshold above every candidate leaves nothing.
  TEST_EQUAL(selectNotches(cands, 17, 6.0, 9).size(), 1)    // only z=15 (snr 9.0)
  TEST_EQUAL(selectNotches(cands, 17, 100.0, 9).size(), 0)

  // Anchoring elsewhere moves which charge is excluded, and the anchor's own SNR is irrelevant.
  {
    auto out = selectNotches(cands, /*anchor=*/15, 1.0, 9);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].charge, 17)   // snr 8.0
    TEST_EQUAL(out[1].charge, 16)   // snr 5.0
  }

  // Not isolatable is not the same as below noise: a candidate with no geometry is dropped whatever
  // its SNR, because a window needs a centre and a width.
  {
    const std::vector<NotchCandidate> broken = {
      {16, 0.0,   3.0, 99.0},   // no m/z
      {15, 883.9, 0.0, 99.0},   // no width
      { 0, 512.1, 2.0, 99.0},   // no charge
      {14, 820.0, 2.5,  4.0},   // fine
    };
    auto out = selectNotches(broken, 17, 1.0, 9);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].charge, 14)
  }

  // Empty in, empty out.
  TEST_EQUAL(selectNotches({}, 17, 1.0, 9).size(), 0)
}
END_SECTION

// writeNotchesForStage is the only writer of the packed slots, and its ordering requirement (stage 0
// before stage 1) is what keeps notchesForStage able to find them.
START_SECTION(writeNotchesForStage_roundTrip)
{
  ScanCommand cmd {};
  cmd.num_stages = 2;
  cmd.stages[0].charge_state = 17;
  cmd.stages[0].collision_energy = 30;
  std::strncpy(cmd.stages[0].activation_type, "HCD", sizeof(cmd.stages[0].activation_type) - 1);
  cmd.stages[1].charge_state = 4;
  cmd.stages[1].collision_energy = 25;
  std::strncpy(cmd.stages[1].activation_type, "CID", sizeof(cmd.stages[1].activation_type) - 1);

  TEST_EQUAL(writeNotchesForStage(cmd, 0, {{16, 938.2, 3.0, 5.0}, {15, 883.9, 2.9, 9.0}}), 2)
  TEST_EQUAL(writeNotchesForStage(cmd, 1, {{5, 1001.2, 2.0, 4.0}}), 1)

  auto s0 = notchesForStage(cmd, 0);
  auto s1 = notchesForStage(cmd, 1);
  TEST_EQUAL(s0.second, 2)
  TEST_EQUAL(s1.second, 1)
  TEST_EQUAL(s0.first[0].charge_state, 16)
  TEST_EQUAL(s1.first[0].charge_state, 5)

  // A notch inherits its stage's fragmentation settings: all notches of a stage fire into the SAME
  // event, and the wire carries collision energy and activation per stage, not per notch.
  TEST_REAL_SIMILAR(s0.first[0].collision_energy, 30.0)
  TEST_EQUAL(std::string(s0.first[0].activation_type), "HCD")
  TEST_REAL_SIMILAR(s1.first[0].collision_energy, 25.0)
  TEST_EQUAL(std::string(s1.first[0].activation_type), "CID")

  // Wire arity: one ';'-group per cascade stage, each 1 + its notch count.
  TEST_EQUAL(windowsAtStage(cmd, 0), 3)
  TEST_EQUAL(windowsAtStage(cmd, 1), 2)

  // The ceiling a returning spectrum is deconvolved against is the highest charge ISOLATED, not the
  // anchor's -- capping at the anchor would discard fragments of the higher members.
  TEST_EQUAL(maxIsolatedCharge(cmd, 0), 17)   // anchor 17 still the highest here
  TEST_EQUAL(maxIsolatedCharge(cmd, 1), 5)    // notch z=5 exceeds the anchor's 4
}
END_SECTION

END_TEST
