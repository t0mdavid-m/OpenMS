// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#pragma once

#include <chrono>
#include <cstdio>
#include <cstdlib>

namespace OpenMS
{
  /**
    @brief Opt-in wall-clock phase marks of FLASHDeconv and FLASHTnT, for runtime analysis.

    With the environment variable FLASH_TRACE_PHASES set, each call writes "[phase] <seconds> <name>" to stderr, the seconds of
    std::steady_clock (CLOCK_MONOTONIC on Linux), so that the marks line up with a sampler of the same clock. Without it, a call reads
    one cached flag and returns. Call it from the main thread only (the marks of parallel workers would interleave).
  */
  inline void flashPhase(const char* name)
  {
    static const bool on = std::getenv("FLASH_TRACE_PHASES") != nullptr;
    if (! on) return;
    const double t = std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
    std::fprintf(stderr, "[phase] %.3f %s\n", t, name);
  }
} // namespace OpenMS
