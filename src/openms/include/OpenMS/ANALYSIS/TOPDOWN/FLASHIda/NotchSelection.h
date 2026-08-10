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
  /// measured geometry needed to isolate it, the SNR that decides whether it is worth isolating at
  /// all, and the intensity that decides which survive a clamp.
  ///
  /// mz and width are MEASURED (PeakGroup::getMzRange for a precursor, the stored observation for a
  /// fragment), never derived from theory — an isolation window must match what was actually seen.
  ///
  /// snr and intensity have distinct jobs and are not interchangeable: SNR is the ADMISSION test
  /// ("is there signal here at all"), intensity is the RANKING ("which charges carry the ion current
  /// worth spending a fill on"). A weak-but-clean charge passes the gate and still loses the slot.
  struct OPENMS_DLLAPI NotchCandidate
  {
    int charge = 0;
    double mz = 0;
    double width = 0;
    double snr = 0;
    double intensity = 0;
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
  ///   3. order by descending INTENSITY, so a clamp keeps the charges carrying the most ion current,
  ///   4. clamp to @p max_notches and SAY what was dropped — a silent truncation reads as "we
  ///      isolated the whole envelope" when we did not.
  ///
  /// Ordering by intensity rather than SNR is deliberate. SNR is a purity measure: a charge sitting
  /// in a clean part of the spectrum can out-score a far more abundant one in a crowded region, and
  /// under a clamp that trades away most of the envelope's ion current for tidiness. What a
  /// co-isolated fill is *for* is harvesting that current, so intensity ranks and SNR only admits.
  ///
  /// The anchor is always kept — it is the stage itself and the scan's identity channel (ADR-0008) —
  /// so the isolated set is the anchor plus the @p max_notches most intense of the rest. With an
  /// envelope of 10 or fewer charges that is every SNR-positive charge, and the order is moot.
  ///
  /// @param candidates every observed charge of the species, anchor included
  /// @param anchor_charge the charge the cascade stage itself isolates
  /// @param snr_threshold the same targeting().snr_threshold the selector already applies
  /// @param max_notches MAX_NOTCHES_PER_STAGE — this stage's own budget, not shared with the other
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

    // Descending intensity, charge as the tiebreak so the order is total and the output reproducible
    // (two charges of one envelope can carry byte-equal intensities).
    std::sort(out.begin(), out.end(),
              [](const NotchCandidate& a, const NotchCandidate& b) {
                if (a.intensity != b.intensity) return a.intensity > b.intensity;
                return a.charge < b.charge;
              });

    if (static_cast<int>(out.size()) > max_notches)
    {
      const int dropped = static_cast<int>(out.size()) - max_notches;
      out.resize(max_notches);
      std::cout << "[NOTCH-CLAMP] " << where << " kept=" << max_notches << " dropped=" << dropped
                << " (least intense charge states; cap is MAX_NOTCHES_PER_STAGE, the instrument's"
                   " 10 MSX windows per fragmentation stage minus the anchor)" << std::endl;
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
      cands.push_back({c, (lo + hi) / 2.0, (hi + margin) - (lo - margin), pg.getChargeSNR(c),
                       static_cast<double>(pg.getChargeIntensity(c))});
    }
    return cands;
  }

  /// Write @p notches into cascade stage @p k's notch block and set that stage's count.
  ///
  /// Stage k's block is fixed at [k * MAX_NOTCHES_PER_STAGE, + MAX_NOTCHES_PER_STAGE), so the two
  /// stages can be written in either order and neither can consume the other's slots. Returns how
  /// many were written — fewer than requested only if the caller passed more than a stage's cap,
  /// which selectNotches' clamp normally prevents.
  ///
  /// No collision energy or activation is written: a notch is geometry only, because every notch of a
  /// stage fires into the SAME fragmentation event and the wire carries one CE and one activation per
  /// ';' group. That used to be a copy from stages[k] into a spare stage slot; the Notch record makes
  /// it structural instead, so there is no per-notch CE to drift from the stage's.
  inline int writeNotchesForStage(ScanCommand& cmd, int k, const std::vector<NotchCandidate>& notches)
  {
    if (k < 0 || k > 1) return 0;
    const int begin = k * MAX_NOTCHES_PER_STAGE;
    int written = 0;
    for (const NotchCandidate& n : notches)
    {
      if (written >= MAX_NOTCHES_PER_STAGE) break;
      Notch& slot = cmd.notches[begin + written];
      slot.precursor_mz = n.mz;
      slot.isolation_width = n.width;
      slot.charge_state = n.charge;
      slot.pad_ = 0;
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
