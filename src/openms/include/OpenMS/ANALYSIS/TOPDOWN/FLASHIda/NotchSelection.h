// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>
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
  /// SNR is the ADMISSION test ("is there signal here at all") and always stays one. Which of the
  /// admitted charges rank first is authored (ChargeSelection, ADR-0059): snr, qscore, intensity or
  /// purity. qscore is the charge's own per-charge qscore; purity is its isolation purity on the
  /// survey the candidates were built from, -1 when no spectrum was supplied.
  struct OPENMS_DLLAPI NotchCandidate
  {
    int charge = 0;
    double mz = 0;
    double width = 0;
    double snr = 0;
    double intensity = 0;
    double qscore = 0;
    double purity = -1.0;
  };

  /// The ranking key of @p c under @p rule. Higher ranks first.
  inline double chargeRankKey(const NotchCandidate& c, ChargeSelection rule)
  {
    switch (rule)
    {
      case ChargeSelection::Snr:       return c.snr;
      case ChargeSelection::Qscore:    return c.qscore;
      case ChargeSelection::Intensity: return c.intensity;
      case ChargeSelection::Purity:    return c.purity;
    }
    return c.intensity;
  }

  /// Order @p candidates by @p rule: the rule's key descending, then intensity descending, then the
  /// lower charge -- total, so the output is reproducible when two charges tie on the key.
  inline void rankCharges(std::vector<NotchCandidate>& candidates, ChargeSelection rule)
  {
    std::sort(candidates.begin(), candidates.end(),
              [rule](const NotchCandidate& a, const NotchCandidate& b) {
                const double ka = chargeRankKey(a, rule);
                const double kb = chargeRankKey(b, rule);
                if (ka != kb) return ka > kb;
                if (a.intensity != b.intensity) return a.intensity > b.intensity;
                return a.charge < b.charge;
              });
  }

  /// Pick the co-isolation notches for one acquisition charge set (ADR-0016).
  ///
  /// Policy, in one place because three call sites need the same answer and a second copy is how the
  /// exclusion side and the geometry side drift apart:
  ///   1. drop the anchor (it is the cascade stage itself, not a notch),
  ///   2. drop any charge whose own envelope does not rise above noise (snr < @p snr_threshold) —
  ///      a charge contributing no signal still consumes part of the scan's ion budget, and under
  ///      MSX the budget is shared with equal injection time per notch, so a junk notch actively
  ///      costs the good ones,
  ///   3. order by @p rule (rankCharges), so a clamp keeps the charges the rule ranks first,
  ///   4. clamp to @p max_notches and SAY what was dropped — a silent truncation reads as "we
  ///      isolated the whole envelope" when we did not.
  ///
  /// The rule is authored for a precursor (precursor_selection.charge_selection, ADR-0059) and is
  /// the same one that chose the anchor, so the anchor is the set's first member by construction.
  /// MS3 passes Intensity: a co-isolated fragment fill is for harvesting ion current, and the
  /// authored key is a precursor_selection key.
  ///
  /// The anchor is always kept — it is the stage itself and the scan's identity channel (ADR-0008) —
  /// so the isolated set is the anchor plus the @p max_notches first-ranked of the rest. With an
  /// envelope of 10 or fewer charges and no authored cap that is every SNR-positive charge.
  ///
  /// @param candidates every observed charge of the species, anchor included
  /// @param anchor_charge the charge the cascade stage itself isolates
  /// @param snr_threshold the same targeting().snr_threshold the selector already applies
  /// @param max_notches this stage's own budget (at most MAX_NOTCHES_PER_STAGE)
  /// @param rule the ranking that decides order and which charges survive the clamp
  /// @param where short label for the drop message (e.g. "MS2 z=17 m=12358.3")
  inline std::vector<NotchCandidate> selectNotches(std::vector<NotchCandidate> candidates,
                                                   int anchor_charge, double snr_threshold,
                                                   int max_notches, ChargeSelection rule,
                                                   const std::string& where = "")
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

    rankCharges(out, rule);

    if (static_cast<int>(out.size()) > max_notches)
    {
      const int dropped = static_cast<int>(out.size()) - max_notches;
      out.resize(max_notches);
      std::cout << "[NOTCH-CLAMP] " << where << " kept=" << max_notches << " dropped=" << dropped
                << " (last-ranked charge states; cap is max_charge_states - 1, at most MAX_NOTCHES_PER_STAGE,"
                   " the instrument's 10 MSX windows per fragmentation stage minus the anchor)" << std::endl;
    }
    return out;
  }

  /// The per-side widening applied to a MEASURED m/z range before it is commanded as an isolation
  /// window, so the envelope is transmitted whole rather than clipped at its edges: the species does
  /// not stop at the outermost isotope the deconvolution happened to resolve, and an unpadded window
  /// shaves the wings off the fill.
  ///
  /// Applied SYMMETRICALLY at every site — a window built from a measured [mz1, mz2] is
  /// (mz2 - mz1) + 2 * this value, i.e. 0.8 Th wider than the span, never 0.4.
  ///
  /// It lives here because peakGroupNotchCandidates below already takes it as @p margin and every
  /// other consumer must agree with it. It previously existed as THREE independent TU-local copies —
  /// ScanCommandQueue.cpp, PrecursorSelection.cpp, and FragmentAnalysis.cpp's anonymous namespace —
  /// all .4, with nothing keeping them in step. Exploration used none of the three: it took the bare
  /// PeakGroup span (Exploration.cpp:164-168), so the window an exploration metric SUMMED came out
  /// 0.8 Th narrower than the one buildMS2 COMMANDED for the very same precursor, and every score
  /// was computed over an interval the instrument was never asked to isolate. ADR-0026 decision 2
  /// closed that gap: the copies were collapsed onto this symbol and Exploration was given it too.
  ///
  /// The invariant that holds now, at MS2: the interval a metric sums and the interval the
  /// instrument was told to isolate are the SAME interval. A fourth copy, or a consumer reaching for
  /// the bare span again, reopens the gap silently — such scores stay perfectly plausible, they just
  /// describe a window that was never acquired.
  ///
  /// The margin is also the window FLOOR. A charge resolved at a single isotope has mz2 == mz1, so
  /// without it the commanded isolation_width is 0 and first_mass == last_mass == precursor_mz —
  /// both positive, which FlashIDA/src/Flash/ScanFactory.cs:245-247 sends to the instrument as a
  /// zero-width "DefineMZRange".
  inline constexpr double optimal_window_margin_ = 0.4;

  /// Every observed charge of one MS1 PeakGroup, as notch candidates with MEASURED geometry.
  ///
  /// Shared by PrecursorSelection (which records the acquired set for charge-keyed exclusion) and
  /// ScanCommandQueue::buildMS2 (which writes the isolation geometry). Sharing the extraction, not
  /// just the policy, is what makes the set recorded as acquired identical BY CONSTRUCTION to the set
  /// actually isolated — two independent walks of the same PeakGroup is exactly how the exclusion
  /// side and the geometry side would drift apart.
  ///
  /// @p margin is optimal_window_margin_ above, applied symmetrically as it is for the anchor stage,
  /// so a notch's width is computed the same way the cascade stage's is.
  ///
  /// @p raw is the survey's raw spectrum. When given, each candidate's isolation purity is measured
  /// over exactly the window it would be commanded with (FragmentAnalysis::windowPurity); without it
  /// purity stays -1, which only the Purity rule reads. qscore is the charge's own per-charge qscore,
  /// falling back to the PeakGroup's qscore for one that never went through updateQscore.
  inline std::vector<NotchCandidate> peakGroupNotchCandidates(const PeakGroup& pg, double margin,
                                                              const MSSpectrum* raw = nullptr)
  {
    std::vector<NotchCandidate> cands;
    const std::unordered_map<int, float> qscores = pg.getAllQscores();
    auto [min_c, max_c] = pg.getAbsChargeRange();
    for (int c = min_c; c <= max_c; ++c)
    {
      auto [lo, hi] = pg.getMzRange(c);
      if (hi <= lo) continue;   // charge not present in this envelope
      NotchCandidate n;
      n.charge = c;
      n.mz = (lo + hi) / 2.0;
      n.width = (hi + margin) - (lo - margin);
      n.snr = pg.getChargeSNR(c);
      n.intensity = static_cast<double>(pg.getChargeIntensity(c));
      const auto qs = qscores.find(c);
      n.qscore = qs != qscores.end() ? static_cast<double>(qs->second) : pg.getQscore();
      if (raw != nullptr)
        n.purity = FragmentAnalysis::windowPurity(*raw, n.mz - n.width / 2.0, n.mz + n.width / 2.0, n.intensity);
      cands.push_back(n);
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
