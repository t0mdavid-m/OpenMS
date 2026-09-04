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

    // ---- ADR-0042: the monitor scan's two observables --------------------------------------------

    /// The engine's ACQUISITION MEMORY -- everything PrecursorSelection remembers across surveys.
    ///
    /// This is the read side of the neutrality drift-guard. It returns the whole struct BY VALUE and
    /// compares with the compiler-generated operator==, so a container added to SurveyMemory later is
    /// covered the moment it is declared, with no edit here and none in the test. Do not replace this
    /// with a hand-picked tuple of members -- that is precisely the checklist step the struct exists
    /// to delete.
    static PrecursorSelection::SurveyMemory surveyMemory(const FLASHIda& f) { return f.selection_.memory_; }

    /// Next precursor id the engine would mint. Part of the neutrality fingerprint: a monitor scan
    /// must not consume one, and an id counter that moved is the cheapest evidence that it selected.
    static int nextPrecursorId(const FLASHIda& f) { return f.next_precursor_id_; }

    /// The FAIMS wheel's current CV. A monitor scan must not advance it (advanceToNextCV lives on the
    /// survey arm), so this is the third leg of the fingerprint.
    ///
    /// The isEnabled() guard is NOT optional and is copied from the engine's own two call sites
    /// (FLASHIda.cpp's ctor and publishExplorationState_, both `faims_.isEnabled() ? currentCV() : 0.0`).
    /// FAIMS::currentCV() is a bare `cv_values_[current_cv_index_]` with no empty check, so on a
    /// config with `"cv_values": []` -- which is every non-FAIMS config, i.e. most of them -- calling
    /// it unguarded reads off the end of an empty vector and SEGFAULTs. Do not "simplify" this.
    static double faimsCurrentCV(const FLASHIda& f)
    {
      return f.faims_.isEnabled() ? f.faims_.currentCV() : 0.0;
    }

    /// Back-date a level's monitor clock so a WALL-CLOCK cadence becomes deterministic in a test.
    /// @param idx 0 = level 2, 1 = level 3.
    ///
    /// This exists because the authored cadence is in milliseconds of real time, which no test can
    /// wait for. Without it the interval branch of takeMonitorDue would ship unexercised -- exactly
    /// the hole ADR-0031 documents for agc_interval_seconds, which all 43 committed configs pin at
    /// 9999999 so that its interval logic has never once been executed by a test.
    static void backdateMonitorClock(FLASHIda& f, int idx, double ms) { f.queue_.backdateMonitorClock(idx, ms); }

    /// Is the monitor scan's private Deconvolution engine present? Null unless a level enables
    /// monitor_ms1, which is what makes the feature free when off.
    static bool hasMonitorDeconv(const FLASHIda& f) { return f.monitor_deconv_ != nullptr; }

    /// Retention time stamped on each of the two MS1 deconvolution results. This is the ONLY
    /// observable that separates the two MS1 arms unconditionally.
    ///
    /// Everything else is contingent on the SPECTRUM. A survey writes acquisition memory and mints
    /// ids only if it actually selects; and peak-group COUNT is no better, because a spectrum can
    /// deconvolve to nothing while the survey arm still ran in full -- which is exactly what
    /// ms1_standard's first scan does. But both arms call deconvolveMS1 unconditionally, and that
    /// stamps the fed rt onto the stored spectrum whatever it contains. The arms deconvolve into
    /// DIFFERENT engines by construction (survey -> deconv_, monitor -> its own monitor_deconv_),
    /// so feeding a distinctive rt and asking which engine now carries it says which arm ran.
    static double surveyDeconvRT(const FLASHIda& f)
    {
      return f.deconv_.deconvolvedMS1().getOriginalSpectrum().getRT();
    }
    static double monitorDeconvRT(const FLASHIda& f)
    {
      return f.monitor_deconv_ ? f.monitor_deconv_->deconvolvedMS1().getOriginalSpectrum().getRT() : -1.0;
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
