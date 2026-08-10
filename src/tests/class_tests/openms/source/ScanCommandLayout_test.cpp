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
  // Carved out of reserved_, which has now moved 1448 -> 1896 and shrunk 600 -> 152 across three
  // changes (faims_enabled, ADR-0012; the two notch counts, ADR-0017; the notch array, ADR-0019).
  // Every offset above is unchanged and sizeof stays 2048; that is the whole point of consuming from
  // the tail.
  TEST_EQUAL(offsetof(ScanCommand, faims_enabled), 1448)
  TEST_EQUAL(offsetof(ScanCommand, stage0_notch_count), 1452)
  TEST_EQUAL(offsetof(ScanCommand, stage1_notch_count), 1456)
  // pad4 is EXPLICIT: notches is 8-aligned (two doubles) and 1460 is 4-mod-8, so without it the
  // compiler inserts the same four bytes silently and the C# mirror has nothing to line up against.
  TEST_EQUAL(offsetof(ScanCommand, pad4), 1460)
  TEST_EQUAL(offsetof(ScanCommand, notches), 1464)
  TEST_EQUAL(sizeof(ScanCommand::notches), 432)   // 18 * 24
  TEST_EQUAL(offsetof(ScanCommand, reserved_), 1896)
  TEST_EQUAL(sizeof(ScanCommand) - offsetof(ScanCommand, reserved_), 152)

  // The two tens are DIFFERENT limits and must not be collapsed into one: MAX_ISOLATION_STAGES is
  // the instrument's ';'-axis cap (cascade depth), MAX_NOTCHES_PER_STAGE + 1 is its ','-axis cap
  // (MSX windows per fragmentation stage). Reading the first as a joint budget is what capped an
  // MS3's two stages at 8 windows between them.
  TEST_EQUAL(MAX_ISOLATION_STAGES, 10)
  TEST_EQUAL(MAX_NOTCHES_PER_STAGE, 9)
  TEST_EQUAL(MAX_NOTCHES, 18)

  // IsolationStage field offsets
  TEST_EQUAL(offsetof(IsolationStage, precursor_mz), 0)
  TEST_EQUAL(offsetof(IsolationStage, isolation_width), 8)
  TEST_EQUAL(offsetof(IsolationStage, collision_energy), 16)
  TEST_EQUAL(offsetof(IsolationStage, reaction_time), 24)
  TEST_EQUAL(offsetof(IsolationStage, reagent_max_it), 32)
  TEST_EQUAL(offsetof(IsolationStage, reagent_agc_target), 40)
  TEST_EQUAL(offsetof(IsolationStage, charge_state), 44)
  TEST_EQUAL(offsetof(IsolationStage, activation_type), 48)

  // Notch field offsets -- geometry only, 24 bytes, no per-notch CE or activation by construction.
  TEST_EQUAL(sizeof(Notch), 24)
  TEST_EQUAL(offsetof(Notch, precursor_mz), 0)
  TEST_EQUAL(offsetof(Notch, isolation_width), 8)
  TEST_EQUAL(offsetof(Notch, charge_state), 16)
  TEST_EQUAL(offsetof(Notch, pad_), 20)
}
END_SECTION

// Stage k's notches live in the FIXED block [k * MAX_NOTCHES_PER_STAGE, +MAX_NOTCHES_PER_STAGE) of
// notches[] (ADR-0019). Fixed blocks are what let both stages reach a full 10-plex and what make
// write order irrelevant; the previous packed-shared arrangement had stage-1's offset depend on
// stage-0's count, so writing out of order aliased them.
START_SECTION(notchesForStage_packing)
{
  // MS2, two notches at stage 0.
  ScanCommand ms2 {};
  ms2.num_stages = 1;
  ms2.stage0_notch_count = 2;
  ms2.stages[0].charge_state = 17;   // anchor lives in stages[], not notches[]
  ms2.notches[0].charge_state = 16;
  ms2.notches[1].charge_state = 15;

  auto n0 = notchesForStage(ms2, 0);
  TEST_EQUAL(n0.second, 2)
  TEST_EQUAL(n0.first[0].charge_state, 16)
  TEST_EQUAL(n0.first[1].charge_state, 15)
  TEST_EQUAL(windowsAtStage(ms2, 0), 3)
  // An MS2 has no cascade stage 1, so no windows there. (Its notch BLOCK still exists and is simply
  // left at count 0 -- the block is addressed by k, not by num_stages.)
  TEST_EQUAL(notchesForStage(ms2, 1).second, 0)
  TEST_EQUAL(windowsAtStage(ms2, 1), 0)

  // MS3 with notches at BOTH stages. Stage 1's block starts at slot 9 regardless of stage 0's count.
  ScanCommand ms3 {};
  ms3.num_stages = 2;
  ms3.stage0_notch_count = 2;
  ms3.stage1_notch_count = 3;
  ms3.stages[0].charge_state = 17;   // cascade stage 0 anchor
  ms3.stages[1].charge_state = 4;    // cascade stage 1 anchor
  ms3.notches[0].charge_state = 16;  // stage-0 notch 1
  ms3.notches[1].charge_state = 15;  // stage-0 notch 2
  ms3.notches[MAX_NOTCHES_PER_STAGE + 0].charge_state = 5;   // stage-1 notch 1
  ms3.notches[MAX_NOTCHES_PER_STAGE + 1].charge_state = 3;   // stage-1 notch 2
  ms3.notches[MAX_NOTCHES_PER_STAGE + 2].charge_state = 2;   // stage-1 notch 3

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

  // Stage 0's block cannot reach into stage 1's: a stage-0 count of 9 stops one slot short of
  // notches[9], which is stage 1's first. This is the guarantee the shared pool did not give.
  ScanCommand full {};
  full.num_stages = 2;
  full.stage0_notch_count = MAX_NOTCHES_PER_STAGE;
  full.stage1_notch_count = MAX_NOTCHES_PER_STAGE;
  TEST_EQUAL(notchesForStage(full, 0).second, 9)
  TEST_EQUAL(notchesForStage(full, 1).second, 9)
  TEST_EQUAL(notchesForStage(full, 0).first + MAX_NOTCHES_PER_STAGE
             == notchesForStage(full, 1).first, true)
  // Both stages at a full 10-plex simultaneously -- 20 isolation windows in one MS3 command, which
  // the shared 8-slot pool could not express at all.
  TEST_EQUAL(windowsAtStage(full, 0), 10)
  TEST_EQUAL(windowsAtStage(full, 1), 10)

  // No notches: an empty range, so a caller may loop without a count check.
  ScanCommand plain {};
  plain.num_stages = 1;
  TEST_EQUAL(notchesForStage(plain, 0).second, 0)
  TEST_EQUAL(notchesForStage(plain, 0).first == nullptr, true)
  TEST_EQUAL(windowsAtStage(plain, 0), 1)

  // A count past this stage's own cap is refused rather than read into the next stage's block.
  ScanCommand over {};
  over.num_stages = 1;
  over.stage0_notch_count = MAX_NOTCHES_PER_STAGE + 1;
  TEST_EQUAL(notchesForStage(over, 0).second, 0)

  // A stage index outside {0,1} has no notch counter and yields nothing.
  TEST_EQUAL(notchesForStage(ms3, 2).second, 0)
  TEST_EQUAL(notchesForStage(ms3, -1).second, 0)
}
END_SECTION

// selectNotches is the one place the co-isolation policy lives, called from three sites, so its
// rules are pinned here rather than inferred from any one caller (ADR-0016).
START_SECTION(selectNotches_policy)
{
  // charge, mz, width, snr, intensity. z=17 is the anchor; z=9 is below noise. SNR and intensity are
  // deliberately ANTI-CORRELATED here (z=15 is the cleanest, z=16 the most abundant) so the two
  // possible ranking rules give different answers and the test can tell them apart -- with a
  // correlated fixture an SNR-ordered implementation would pass unnoticed.
  const std::vector<NotchCandidate> cands = {
    {17, 1000.5, 3.2, 8.0, 5.0e6},
    {16,  938.2, 3.0, 5.0, 9.0e6},
    {15,  883.9, 2.9, 9.0, 2.0e6},
    { 9,  512.1, 2.0, 0.4, 8.0e6},   // abundant but below noise: intensity does not buy admission
  };

  // The anchor is the cascade stage itself, never a notch; the sub-threshold charge is dropped even
  // though it is the second most abundant; what survives is ordered by DESCENDING INTENSITY, because
  // a co-isolated fill exists to harvest ion current and SNR only decides admission.
  {
    auto out = selectNotches(cands, /*anchor=*/17, /*snr_threshold=*/1.0, /*max_notches=*/9);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].charge, 16)   // intensity 9.0e6, despite the lower SNR
    TEST_EQUAL(out[1].charge, 15)   // intensity 2.0e6, despite the highest SNR
  }

  // The clamp keeps the top N by intensity. (It also reports the drop on stdout; that a value was
  // dropped is never silent, which is the point.)
  {
    auto out = selectNotches(cands, 17, 1.0, 1);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].charge, 16)
  }

  // Equal intensities fall back to ascending charge, so the output is a total order and reproducible
  // rather than dependent on the input container's iteration order (ms2_by_charge is unordered).
  {
    const std::vector<NotchCandidate> tied = {
      {12, 700.0, 2.0, 4.0, 1.0e6},
      {14, 800.0, 2.0, 4.0, 1.0e6},
      {13, 750.0, 2.0, 4.0, 1.0e6},
    };
    auto out = selectNotches(tied, 17, 1.0, 9);
    TEST_EQUAL(out.size(), 3)
    TEST_EQUAL(out[0].charge, 12)
    TEST_EQUAL(out[1].charge, 13)
    TEST_EQUAL(out[2].charge, 14)
  }

  // A full 10-plex: nine notches survive alongside the anchor, which is the per-stage ceiling.
  {
    std::vector<NotchCandidate> wide;
    for (int z = 2; z <= 20; ++z) wide.push_back({z, 500.0 + z, 2.0, 4.0, 1.0e6 * z});
    auto out = selectNotches(wide, /*anchor=*/20, 1.0, MAX_NOTCHES_PER_STAGE);
    TEST_EQUAL(out.size(), 9)
    TEST_EQUAL(out[0].charge, 19)   // most intense of the non-anchors
    TEST_EQUAL(out[8].charge, 11)   // ninth; z<=10 dropped as least intense
  }

  // max_notches 0 -- an MS3 whose stage-0 notches already consumed every slot -- yields nothing
  // rather than overflowing.
  TEST_EQUAL(selectNotches(cands, 17, 1.0, 0).size(), 0)

  // A higher threshold gates more out; a threshold above every candidate leaves nothing.
  TEST_EQUAL(selectNotches(cands, 17, 6.0, 9).size(), 1)    // only z=15 (snr 9.0)
  TEST_EQUAL(selectNotches(cands, 17, 100.0, 9).size(), 0)

  // Anchoring elsewhere moves which charge is excluded, and the anchor's own scores are irrelevant.
  {
    auto out = selectNotches(cands, /*anchor=*/15, 1.0, 9);
    TEST_EQUAL(out.size(), 2)
    TEST_EQUAL(out[0].charge, 16)   // intensity 9.0e6
    TEST_EQUAL(out[1].charge, 17)   // intensity 5.0e6
  }

  // Not isolatable is not the same as below noise: a candidate with no geometry is dropped whatever
  // its SNR or intensity, because a window needs a centre and a width.
  {
    const std::vector<NotchCandidate> broken = {
      {16, 0.0,   3.0, 99.0, 9.0e6},   // no m/z
      {15, 883.9, 0.0, 99.0, 9.0e6},   // no width
      { 0, 512.1, 2.0, 99.0, 9.0e6},   // no charge
      {14, 820.0, 2.5,  4.0, 1.0e6},   // fine
    };
    auto out = selectNotches(broken, 17, 1.0, 9);
    TEST_EQUAL(out.size(), 1)
    TEST_EQUAL(out[0].charge, 14)
  }

  // Empty in, empty out.
  TEST_EQUAL(selectNotches({}, 17, 1.0, 9).size(), 0)
}
END_SECTION

// writeNotchesForStage is the only writer of the notch blocks. With fixed per-stage blocks it has no
// ordering requirement left, which this pins by writing stage 1 FIRST -- under the old packed layout
// that put stage-1's notches where stage-0's were expected.
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

  // Stage 1 written BEFORE stage 0, deliberately.
  TEST_EQUAL(writeNotchesForStage(cmd, 1, {{5, 1001.2, 2.0, 4.0, 3.0e6}}), 1)
  TEST_EQUAL(writeNotchesForStage(cmd, 0, {{16, 938.2, 3.0, 5.0, 9.0e6},
                                           {15, 883.9, 2.9, 9.0, 2.0e6}}), 2)

  auto s0 = notchesForStage(cmd, 0);
  auto s1 = notchesForStage(cmd, 1);
  TEST_EQUAL(s0.second, 2)
  TEST_EQUAL(s1.second, 1)
  TEST_EQUAL(s0.first[0].charge_state, 16)
  TEST_EQUAL(s0.first[1].charge_state, 15)
  TEST_EQUAL(s1.first[0].charge_state, 5)
  TEST_REAL_SIMILAR(s0.first[0].precursor_mz, 938.2)
  TEST_REAL_SIMILAR(s0.first[0].isolation_width, 3.0)
  TEST_REAL_SIMILAR(s1.first[0].precursor_mz, 1001.2)

  // Stage 0's write did not disturb stage 1's, and vice versa: disjoint blocks.
  TEST_EQUAL(cmd.stage0_notch_count, 2)
  TEST_EQUAL(cmd.stage1_notch_count, 1)

  // Wire arity: one ';'-group per cascade stage, each 1 + its notch count.
  TEST_EQUAL(windowsAtStage(cmd, 0), 3)
  TEST_EQUAL(windowsAtStage(cmd, 1), 2)

  // The ceiling a returning spectrum is deconvolved against is the highest charge ISOLATED, not the
  // anchor's -- capping at the anchor would discard fragments of the higher members.
  TEST_EQUAL(maxIsolatedCharge(cmd, 0), 17)   // anchor 17 still the highest here
  TEST_EQUAL(maxIsolatedCharge(cmd, 1), 5)    // notch z=5 exceeds the anchor's 4

  // Over-long input is trimmed to this stage's own cap and reported, never spilled into the next
  // stage's block. The caller's selectNotches clamp normally prevents reaching this.
  {
    ScanCommand big {};
    big.num_stages = 2;
    std::vector<NotchCandidate> many;
    for (int z = 2; z <= 20; ++z) many.push_back({z, 500.0 + z, 2.0, 4.0, 1.0e6});
    TEST_EQUAL(writeNotchesForStage(big, 0, many), MAX_NOTCHES_PER_STAGE)
    TEST_EQUAL(big.stage1_notch_count, 0)                  // untouched by stage 0's overflow
    TEST_EQUAL(notchesForStage(big, 1).second, 0)
  }

  // Both stages at a full 10-plex: 20 windows in one command, the ceiling the user asked for.
  {
    ScanCommand full {};
    full.num_stages = 2;
    std::vector<NotchCandidate> nine;
    for (int z = 2; z <= 10; ++z) nine.push_back({z, 500.0 + z, 2.0, 4.0, 1.0e6});
    TEST_EQUAL(writeNotchesForStage(full, 0, nine), 9)
    TEST_EQUAL(writeNotchesForStage(full, 1, nine), 9)
    TEST_EQUAL(windowsAtStage(full, 0), 10)
    TEST_EQUAL(windowsAtStage(full, 1), 10)
  }
}
END_SECTION

END_TEST
