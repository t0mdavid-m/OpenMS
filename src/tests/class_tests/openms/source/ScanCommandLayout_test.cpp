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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include <cstddef>

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
  // Carved out of reserved_, which moved 1448 -> 1452 and shrank 600 -> 596. Every offset above is
  // unchanged and sizeof stays 2048; that is the whole point of consuming from the tail.
  TEST_EQUAL(offsetof(ScanCommand, faims_enabled), 1448)
  TEST_EQUAL(offsetof(ScanCommand, reserved_), 1452)

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

END_TEST
