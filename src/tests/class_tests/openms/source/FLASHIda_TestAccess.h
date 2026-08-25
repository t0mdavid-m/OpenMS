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
#include <mutex>

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
    /// Number of pending (MS3 scan_id -> parent MS2 context) entries. An MS3 that returns on the
    /// REGULAR path identifies only if its id is in here, so this is the observable for "the engine
    /// dispatched an MS3 it will be able to identify". Entries are erased on that return, so read it
    /// while the MS3 commands are still outstanding (drive with an EMPTY ms3_ion_map).
    static size_t ms2ContextCacheSize(const FLASHIda& f) { return f.ms2_context_cache_.size(); }
    static bool explorationActive(const FLASHIda& f)
    {
      // The exact atomic the cycle-time injection gate consults (set on group formation, cleared on
      // completion), so a test can assert suppression against the engine's real decision variable.
      return f.exploration_active_.load(std::memory_order_acquire);
    }

    /// The analysis mutex itself, so a test can HOLD it and observe who blocks.
    ///
    /// processScan takes this function-scoped across its whole body (FLASHIda.cpp:84), including the
    /// MS1 deconvolution. Holding it here therefore simulates "a deconvolution is in flight" exactly,
    /// with no timing race: whatever the drain does, it does it against a genuinely held lock.
    ///
    /// This is the observable for the ONE property that cannot be checked any other way -- that
    /// getNextScanCommand acquires nothing processScan holds. It is scoped by mutex IDENTITY, which is
    /// both its strength (no threshold, no flake) and its limit: it says nothing about the drain
    /// blocking on some OTHER lock, which is why the logger locks are per-stream rather than shared.
    static std::mutex& analysisMutex(FLASHIda& f) { return f.analysis_mutex_; }

    /// tracking_id -> precursor_id, or 0 if untracked. Forwards to the private accessor rather than
    /// reading the map, so the test exercises the same path the drain does.
    static int precursorIdFor(const FLASHIda& f, int tracking_id) { return f.precursorIdForTracking_(tracking_id); }
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
      (void)rt;  // feedResultImpl_ no longer takes rt (this bypass path never used it)
      return e.feedResultImpl_(tracking_id, ms2_deconv, nullptr, nullptr, 0, queue);
    }

    /// Get exploration group by ID (caller must ensure group exists).
    static Exploration::ExplorationGroup group(const Exploration& e, int group_id)
    {
      return e.getGroup(group_id);
    }

    /// Has ANY un-fragmented reference returned for this group yet?
    /// Replaces the former scalar `ExplorationGroup::has_baseline`, which could not survive the
    /// move to one reference per swept activation (ADR-0029).
    static bool hasAnyBaseline(const Exploration::ExplorationGroup& g)
    {
      return !g.baseline_intensity.empty();
    }

    /// The reference intensity of a single-activation group, or -1.0 if none has returned.
    ///
    /// For fixtures that sweep exactly ONE activation, which is every fixture that predates
    /// ADR-0029 -- it asserts the same quantity the old scalar `baseline_intensity` held, without
    /// each call site having to name that fixture's activation string. Tests that specifically
    /// exercise the per-activation keying assert on `baseline_intensity.at("HCD")` directly, so
    /// nothing here hides the key.
    static double soleBaseline(const Exploration::ExplorationGroup& g)
    {
      if (g.baseline_intensity.size() != 1) return -1.0;
      return g.baseline_intensity.begin()->second;
    }
  };
} // namespace OpenMS
