// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <deque>
#include <string>
#include <vector>

/// One mass offset per run and MS level (ADR-0081). Plain state and pure rules, free of Config, so the engine and
/// FLASHIda_MassOffset_test drive exactly the same code. The rules are the analysis' references:
///   - the estimate: tools/v24/off24_ref.py, reading A (prior N(start, 0.6^2 ppm^2), every raw error of a confident hit
///     one observation N(offset, sigma^2), nothing forgotten);
///   - the matching inputs and the fallback vote: tools/v23/off_ref_impl.py (wide_errors, vote), with the v24 restart.
namespace OpenMS
{
  /// The spreads of one raw error, per level, measured on the 12 E. coli runs (results/v24/off24/choices.json); Tom's
  /// 0.58 and 1.63 ppm at the precision the analysis' test vectors were produced with.
  constexpr double kMassOffsetSigmaMs3Ppm = 0.5787868575652204;
  constexpr double kMassOffsetSigmaMs2Ppm = 1.6305928194273822;
  /// Reading A: the prior's width in ppm.
  constexpr double kMassOffsetPriorWidthPpm = 0.6;
  /// A confident MS3 matches at least this many distinct bonds; a confident MS2 at least this many fragments.
  constexpr int kMassOffsetConfidentMs3Bonds = 3;
  constexpr int kMassOffsetConfidentMs2Fragments = 5;

  /// The run mass offset of one MS level (reading A). The state is the start, a count and a sum; the estimate is the
  /// posterior mean (start / w0 + sum / sigma^2) / (1 / w0 + n / sigma^2), w0 = 0.6^2.
  class OPENMS_DLLAPI MassOffsetState
  {
  public:
    MassOffsetState() = default;
    MassOffsetState(double start_ppm, double sigma_ppm);

    /// The offset a scan resolved now is matched at.
    double estimate() const;
    /// Every raw error (ppm) of one confident hit, each one observation.
    void feed(const std::vector<double>& raw_errors);
    /// After a jump of the fallback vote: the prior is centred on the vote's centre and every observation is dropped
    /// (they were taken at a wrong offset).
    void restart(double centre_ppm);

    double start() const { return start_; }
    long count() const { return n_; }

  private:
    double start_ = 0.0;
    double sigma_ = kMassOffsetSigmaMs3Ppm;
    long n_ = 0;
    double sum_ = 0.0;
  };

  /// The outcome of one vote. `result` is one of "too_few_errors", "too_few_scans", "runner_up", "confirms", "jump".
  struct MassOffsetVoteResult
  {
    std::string result;
    double centre = 0.0;
    int errors = 0;      ///< errors within 3 ppm of the densest window's centre
    int scans = 0;       ///< scans with an error within 3 ppm of the centre
    int runner_up = 0;   ///< the best window 8 ppm or more away
    bool has_centre = false;
  };

  /// The fallback vote of the MS3 level (v23 recipe; MS3 only, armed by an authored MS3 start value). It counts the MS3
  /// with at least 5 masses since the last confident one and keeps the wide errors of the last 20 such MS3; when the
  /// counter reaches 20 it votes over them and, on a clear winner more than 3 ppm from the estimate, jumps.
  class OPENMS_DLLAPI MassOffsetVote
  {
  public:
    static constexpr int kMinMasses = 5;
    static constexpr int kCounterLimit = 20;
    static constexpr int kRingScans = 20;
    static constexpr double kWindowPpm = 300.0;
    static constexpr int kMinErrors = 16;
    static constexpr int kMinScans = 5;
    static constexpr double kNearPpm = 3.0;
    static constexpr double kRunnerUpPpm = 8.0;
    static constexpr double kJumpPpm = 3.0;

    /// The raw error (ppm) of every mass against the closest ladder mass, kept if within +-300 ppm
    /// (off_ref_impl.wide_errors). @p ladder_sorted ascending.
    static std::vector<double> wideErrors(const std::vector<double>& masses, const std::vector<double>& ladder_sorted);

    /// The vote over error lists (off_ref_impl.vote), judged against @p estimate_ppm. Never fails.
    static MassOffsetVoteResult vote(const std::deque<std::vector<double>>& ring, double estimate_ppm);

    /// One resolved MS3: @p n_masses its deconvolved masses, @p confident whether it fed the state, @p wide its
    /// wideErrors. Returns true when a vote was taken (its outcome in @p out).
    bool step(int n_masses, bool confident, const std::vector<double>& wide, double estimate_ppm,
              MassOffsetVoteResult& out);

    int counter() const { return miss_; }

  private:
    int miss_ = 0;
    std::deque<std::vector<double>> ring_;
  };

  /// The run's two offsets (ADR-0081): one state per MS level, never shared, and the MS3 vote. Built from the method's
  /// deconvolution.mass_offset_start_ppm; the vote is armed only when the method authored an MS3 start value.
  struct OPENMS_DLLAPI MassOffsets
  {
    MassOffsetState ms2;
    MassOffsetState ms3;
    MassOffsetVote vote;
    bool vote_armed = false;

    MassOffsets() = default;
    MassOffsets(double ms2_start_ppm, double ms3_start_ppm, bool ms3_vote_armed);

    /// One scored MS3 (scored at @p used_ppm, the MS3 estimate when it was matched): a confident one (>= 3 distinct bonds)
    /// feeds its raw errors; the vote, when armed, is stepped and may jump. Writes [MASS-OFFSET] (and [MASS-OFFSET-VOTE]
    /// when a vote is taken) for an MS3 with masses. Returns whether it was confident.
    bool onMs3(const std::string& tracking_id, double used_ppm, int n_masses, int matched, int bonds,
               const std::vector<double>& raw_errors, const std::vector<double>& wide_errors);
    /// One MS2 identified at the identification threshold (matched at @p used_ppm): with >= 5 matched fragments it feeds
    /// their raw errors and writes [MASS-OFFSET]. Returns whether it was confident.
    bool onMs2(const std::string& tracking_id, double used_ppm, int matched, const std::vector<double>& raw_errors);
    /// [MASS-OFFSET-START], one line per level.
    void logStart() const;
  };
} // namespace OpenMS
