// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace OpenMS
{
  /// One candidate co-isolation window: a charge state of the species being acquired, with the
  /// measured geometry needed to isolate it and the SNR that decides whether it is worth isolating.
  ///
  /// mz and width are MEASURED (PeakGroup::getMzRange for a precursor, the stored observation for a
  /// fragment), never derived from theory — an isolation window must match what was actually seen.
  struct OPENMS_DLLAPI NotchCandidate
  {
    int charge = 0;
    double mz = 0;
    double width = 0;
    double snr = 0;
  };

  /// Pick the co-isolation notches for one acquisition charge set (ADR-0016).
  ///
  /// Policy, in one place because three call sites need the same answer and a second copy is how the
  /// exclusion side and the geometry side drift apart:
  ///   1. drop the anchor (it is the cascade stage itself, not a notch),
  ///   2. drop any charge whose own envelope does not rise above noise (snr < @p snr_threshold) —
  ///      a charge contributing no signal still consumes part of the scan's ion budget, and under
  ///      MSX the budget is shared with equal injection time per notch, so a junk notch actively
  ///      costs the good ones,
  ///   3. order by descending SNR, so a clamp keeps the strongest,
  ///   4. clamp to @p max_notches and SAY what was dropped — a silent truncation reads as "we
  ///      isolated the whole envelope" when we did not.
  ///
  /// @param candidates every observed charge of the species, anchor included
  /// @param anchor_charge the charge the cascade stage itself isolates
  /// @param snr_threshold the same targeting().snr_threshold the selector already applies
  /// @param max_notches MAX_ISOLATION_STAGES - num_stages, which is also the instrument's own
  ///        10-value PrecursorMass limit
  /// @param where short label for the drop message (e.g. "MS2 z=17 m=12358.3")
  inline std::vector<NotchCandidate> selectNotches(std::vector<NotchCandidate> candidates,
                                                   int anchor_charge, double snr_threshold,
                                                   int max_notches, const std::string& where = "")
  {
    std::vector<NotchCandidate> out;
    if (max_notches <= 0) return out;

    for (const NotchCandidate& c : candidates)
    {
      if (c.charge == anchor_charge) continue;          // the stage itself, not a notch
      if (c.charge <= 0 || c.mz <= 0 || c.width <= 0) continue;  // not isolatable
      if (c.snr < snr_threshold) continue;              // below noise
      out.push_back(c);
    }

    std::sort(out.begin(), out.end(),
              [](const NotchCandidate& a, const NotchCandidate& b) { return a.snr > b.snr; });

    if (static_cast<int>(out.size()) > max_notches)
    {
      const int dropped = static_cast<int>(out.size()) - max_notches;
      out.resize(max_notches);
      std::cout << "[NOTCH-CLAMP] " << where << " kept=" << max_notches << " dropped=" << dropped
                << " (lowest-SNR charge states; cap is MAX_ISOLATION_STAGES - num_stages, which is"
                   " also the instrument's 10-value PrecursorMass limit)" << std::endl;
    }
    return out;
  }

  /// Every observed charge of one MS1 PeakGroup, as notch candidates with MEASURED geometry.
  ///
  /// Shared by PrecursorSelection (which records the acquired set for charge-keyed exclusion) and
  /// ScanCommandQueue::buildMS2 (which writes the isolation geometry). Sharing the extraction, not
  /// just the policy, is what makes the set recorded as acquired identical BY CONSTRUCTION to the set
  /// actually isolated — two independent walks of the same PeakGroup is exactly how the exclusion
  /// side and the geometry side would drift apart.
  ///
  /// @p margin is ScanCommandQueue's optimal_window_margin_, applied symmetrically as it is for the
  /// anchor stage, so a notch's width is computed the same way the cascade stage's is.
  inline std::vector<NotchCandidate> peakGroupNotchCandidates(const PeakGroup& pg, double margin)
  {
    std::vector<NotchCandidate> cands;
    auto [min_c, max_c] = pg.getAbsChargeRange();
    for (int c = min_c; c <= max_c; ++c)
    {
      auto [lo, hi] = pg.getMzRange(c);
      if (hi <= lo) continue;   // charge not present in this envelope
      cands.push_back({c, (lo + hi) / 2.0, (hi + margin) - (lo - margin), pg.getChargeSNR(c)});
    }
    return cands;
  }

  /// Write @p notches into @p cmd's stage slots for cascade stage @p k and set that stage's count.
  ///
  /// Notches are packed at stages[num_stages ...], stage-0's first then stage-1's, so stage 1 must
  /// be written AFTER stage 0 — writing them out of order would place stage-1's notches where
  /// notchesForStage() looks for stage-0's. Returns how many were written (may be fewer than
  /// requested if the slots run out, which selectNotches' clamp normally prevents).
  ///
  /// Collision energy and activation are copied from the stage the notches belong to: all notches of
  /// a stage fire into the SAME fragmentation event, and the wire has no per-notch slot for them.
  inline int writeNotchesForStage(ScanCommand& cmd, int k, const std::vector<NotchCandidate>& notches)
  {
    if (k < 0 || k >= cmd.num_stages) return 0;
    const int begin = cmd.num_stages + ((k == 1) ? cmd.stage0_notch_count : 0);
    int written = 0;
    for (const NotchCandidate& n : notches)
    {
      const int slot = begin + written;
      if (slot >= MAX_ISOLATION_STAGES) break;
      cmd.stages[slot] = cmd.stages[k];        // inherit CE / activation / reagent settings
      cmd.stages[slot].precursor_mz = n.mz;
      cmd.stages[slot].isolation_width = n.width;
      cmd.stages[slot].charge_state = n.charge;
      ++written;
    }
    if (k == 0) cmd.stage0_notch_count = written;
    else        cmd.stage1_notch_count = written;
    return written;
  }

  /// The highest charge a scan actually isolates at cascade stage @p k: the anchor or any notch.
  ///
  /// This is the deconvolution charge ceiling for a returning spectrum. A fragment of the z=17
  /// member of a co-isolated set may itself carry z=17, so bounding by the anchor alone would
  /// discard real fragments; bounding by the maximum keeps the "fragment <= precursor" guarantee
  /// valid over the set rather than over one charge (ADR-0016).
  inline int maxIsolatedCharge(const ScanCommand& cmd, int k)
  {
    if (k < 0 || k >= cmd.num_stages) return 0;
    int hi = std::abs(cmd.stages[k].charge_state);
    auto notches = notchesForStage(cmd, k);
    for (int i = 0; i < notches.second; ++i) hi = std::max(hi, std::abs(notches.first[i].charge_state));
    return hi;
  }

} // namespace OpenMS
