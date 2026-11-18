// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <cstdint>

namespace OpenMS
{
  /**
    @brief The clock of the W gate's time base (ADR-0079; package 17, SPEC 2): a survey's resolve_ts and an MS3 plan's
    time t, in ms of the steady clock.

    One function for both, so the MS1 buffer's horizon and the peak features compare times of one clock. A test advances
    it (FLASHIdaTestAccess::advanceGateClock) to place surveys and plans minutes apart; production never does, and then
    nowMs() is the steady clock exactly as every other engine timestamp. Process-wide.
  */
  struct OPENMS_DLLAPI EngineClock
  {
    /// steady_clock ms since its epoch, plus the test offset (0 in production).
    static uint64_t nowMs();
    /// Test only: move nowMs() forward by @p ms.
    static void advanceForTest(int64_t ms);
    /// Test only: back to the steady clock.
    static void resetForTest();
  };
} // namespace OpenMS
