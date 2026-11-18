// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/EngineClock.h>

#include <atomic>
#include <chrono>

namespace OpenMS
{
  namespace
  {
    std::atomic<int64_t> g_test_offset_ms{0};
  }

  uint64_t EngineClock::nowMs()
  {
    const int64_t steady = std::chrono::duration_cast<std::chrono::milliseconds>(
                               std::chrono::steady_clock::now().time_since_epoch()).count();
    return static_cast<uint64_t>(steady + g_test_offset_ms.load(std::memory_order_relaxed));
  }

  void EngineClock::advanceForTest(int64_t ms) { g_test_offset_ms.fetch_add(ms, std::memory_order_relaxed); }

  void EngineClock::resetForTest() { g_test_offset_ms.store(0, std::memory_order_relaxed); }
} // namespace OpenMS
