// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// Test-only accessors for FLASHIda / Exploration private state.
//
// These structs are the SINGLE home for test scaffolding that needs private access. The production
// headers (FLASHIda.h / Exploration.h) expose only `friend struct FLASHIdaTestAccess;` /
// `friend struct ExplorationTestAccess;` — no *ForTest methods leak into the shipped API.
//
// They are NAMED structs at namespace scope (NOT an anonymous namespace): a `friend` declaration names
// one concrete type, so the befriended type must be identical across every test translation unit. The
// methods are `static inline`, so including this header in multiple test TUs does not violate the ODR.
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Exploration.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include <atomic>
#include <cstddef>

namespace OpenMS
{
  /// Reaches into FLASHIda's private queue_ / exploration_active_ for the test suite.
  /// Replaces the former FLASHIda::pushCommandForTest / getQueueForTest / getQueueSizeForTest /
  /// explorationActive methods.
  struct FLASHIdaTestAccess
  {
    static ScanCommandQueue& queue(FLASHIda& f) { return f.queue_; }
    static void push(FLASHIda& f, ScanCommand cmd) { f.queue_.push(cmd); }
    static size_t queueSize(const FLASHIda& f, int priority) { return f.queue_.queueSize(priority); }
    static bool explorationActive(const FLASHIda& f)
    {
      // The exact atomic the cycle-time injection gate consults (set on group formation, cleared on
      // completion), so a test can assert suppression against the engine's real decision variable.
      return f.exploration_active_.load(std::memory_order_acquire);
    }
  };

  /// Reaches into Exploration's private feedResultImpl_ / getGroup for the test suite.
  /// Replaces the former Exploration::feedResultForTest and the public Exploration::getGroup.
  struct ExplorationTestAccess
  {
    /// Feed a pre-deconvolved result, bypassing deconvolution (was Exploration::feedResultForTest).
    static Exploration::FeedResultInfo feedResult(Exploration& e, int tracking_id,
                                                  const DeconvolvedSpectrum& ms2_deconv,
                                                  double rt, ScanCommandQueue& queue)
    {
      return e.feedResultImpl_(tracking_id, ms2_deconv, nullptr, nullptr, 0, rt, queue);
    }

    /// Get exploration group by ID (caller must ensure group exists).
    static Exploration::ExplorationGroup group(const Exploration& e, int group_id)
    {
      return e.getGroup(group_id);
    }
  };
} // namespace OpenMS
