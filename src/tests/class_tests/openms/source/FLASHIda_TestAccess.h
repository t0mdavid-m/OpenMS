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

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <mutex>
#include <utility>
#include <vector>

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

    // ---- The target-log loader (PrecursorSelection's constructor, targeting.target_log_files) --------

    /// What the loader filed for exclusion_masses (mode 3): survey rt in SECONDS -> that survey's
    /// AllMass= masses. Taken from a PrecursorSelection, not a FLASHIda, so a test can build one from a
    /// Config alone -- the loader needs no spectrum.
    static std::map<double, std::vector<double>> exclusionRtMasses(const PrecursorSelection& s)
    {
      return s.exclusion_rt_masses_map_;
    }

    /// What the loader filed for inclusion / in_depth (modes 1, 2): reference mass -> the rts (SECONDS)
    /// of its Mass= lines.
    static std::map<double, std::vector<double>> targetMassRts(const PrecursorSelection& s)
    {
      return s.target_mass_rt_map_;
    }

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

    /// A PrecursorSelection's acquisition memory, WRITABLE. For a member that no engine path reaches --
    /// removeFromExclusionList has had no caller since its C export was removed -- the only way to
    /// exercise it is to lay down the record a survey would have written and call it directly.
    static PrecursorSelection::SurveyMemory& selectionMemory(PrecursorSelection& s) { return s.memory_; }

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

    // ---- One survey's selection on chosen PeakGroups (PrecursorSelection) -------------------------

    /// Run one MS1 survey's selection at @p rt (SECONDS) on @p peak_groups in place of a deconvolved
    /// spectrum, and return how many acquisitions it selected. filterAndRank runs first on an EMPTY
    /// spectrum: that rebuilds the survey's target list exactly as a real survey does and deconvolves
    /// nothing. The PeakGroups then go where the deconvolution leaves its result, and the selection
    /// proper runs on them -- so order, masses and charges are the test's, not a fixture's.
    static int selectFrom(PrecursorSelection& s, const std::vector<PeakGroup>& peak_groups, double rt)
    {
      s.filterAndRank(nullptr, nullptr, 0, rt, 1, 0.0);
      DeconvolvedSpectrum& ms1 = s.deconv_.deconvolvedMS1();
      ms1.clear();
      for (const PeakGroup& pg : peak_groups) ms1.push_back(pg);
      s.selected_peak_groups_.clear();
      s.filterPeakGroupsUsingMassExclusion_(1, rt);
      return (int)s.selected_peak_groups_.size();
    }

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

    // ---- Precursor selection on a caller-built survey --------------------------------------------

    /// Run ONE survey's precursor selection on a deconvolved MS1 the caller built, skipping the
    /// deconvolution. Returns one (monoisotopic mass, isolated charge) pair per acquisition the survey
    /// emitted, in emit order -- under "separate" a species contributes one pair per charge.
    ///
    /// It is filterAndRank minus deconvolveMS1: the selection is cleared, the engine's deconvolved
    /// MS1 is replaced, and the same filterPeakGroupsUsingMassExclusion_ runs with the engine's own
    /// config and acquisition memory. It exists because a synthetic PeakGroup can put more
    /// resolvable species into one survey than any committed fixture does (ms1_standard yields <= 1
    /// selectable mass per scan, ms1_cytc one species), and the slot budget only matters once there
    /// are more species than slots.
    ///
    /// @p target_masses stands in for the list filterAndRank builds for inclusion mode from target
    /// logs (mass-only matching); leave it empty for every other mode. With it non-empty and
    /// strict_inclusion false, the survey runs both the target phase and the non-target phase.
    static std::vector<std::pair<double, int>> selectFromDeconvolved(FLASHIda& f, const DeconvolvedSpectrum& ms1,
                                                                     double rt,
                                                                     std::vector<double> target_masses = {})
    {
      PrecursorSelection& sel = f.selection_;
      std::sort(target_masses.begin(), target_masses.end());
      sel.target_masses_ = std::move(target_masses);
      sel.selected_peak_groups_.clear();
      sel.deconv_.deconvolvedMS1() = ms1;
      sel.filterPeakGroupsUsingMassExclusion_(1, rt);

      std::vector<std::pair<double, int>> out;
      const std::vector<int>& charges = sel.triggerCharges();
      for (size_t i = 0; i < sel.selected_peak_groups_.size() && i < charges.size(); ++i)
        out.emplace_back(sel.selected_peak_groups_[i].getMonoMass(), charges[i]);
      return out;
    }

    // ---- FASTA tag-based targeting: what an inclusion-mode survey matches against -----------------

    /// Add tag-expanded targets exactly as a tag detection does (PrecursorSelection::addDynamicTargets_).
    static void addDynamicTargets(PrecursorSelection& s, const std::vector<double>& masses, double rt, int priority)
    {
      s.addDynamicTargets_(masses, rt, priority);
    }

    /// The sorted target masses of an inclusion-mode survey at @p rt (seconds) -- the list filterAndRank
    /// hands the deconvolution engine and matches every candidate against.
    static std::vector<double> inclusionTargetMasses(PrecursorSelection& s, double rt)
    {
      s.collectInclusionTargets_(rt);
      return s.target_masses_;
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
    // hasAnyBaseline / soleBaseline are gone with the baseline itself (ADR-0048): the precursor share
    // is referenced against the scan's own spectrum, so no group holds a reference intensity.
  };
} // namespace OpenMS
