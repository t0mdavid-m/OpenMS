// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MassOffset.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/EventLog.h>

#include <algorithm>
#include <cmath>
#include <iomanip>

namespace OpenMS
{
  namespace
  {
    double medianOf(std::vector<double> v)
    {
      std::sort(v.begin(), v.end());
      const size_t n = v.size();
      const size_t m = n / 2;
      return (n % 2) ? v[m] : (v[m - 1] + v[m]) / 2.0;
    }
  } // namespace

  // --- MassOffsetState (reading A) -----------------------------------------------------------------

  MassOffsetState::MassOffsetState(double start_ppm, double sigma_ppm) : start_(start_ppm), sigma_(sigma_ppm) {}

  double MassOffsetState::estimate() const
  {
    if (n_ == 0) return start_;
    // The reference's arithmetic, in its order (off24_ref.OffsetState._mean): the prior precision is 1 / 0.6^2.
    const double prior_precision = 1.0 / (kMassOffsetPriorWidthPpm * kMassOffsetPriorWidthPpm);
    const double w = 1.0 / (sigma_ * sigma_);
    return (start_ * prior_precision + sum_ * w) / (prior_precision + static_cast<double>(n_) * w);
  }

  void MassOffsetState::feed(const std::vector<double>& raw_errors)
  {
    for (const double e : raw_errors)
    {
      sum_ += e;   // a running sum in the order fed, as the reference keeps it
      ++n_;
    }
  }

  void MassOffsetState::restart(double centre_ppm)
  {
    start_ = centre_ppm;
    n_ = 0;
    sum_ = 0.0;
  }

  // --- MassOffsetVote (v23 recipe) -----------------------------------------------------------------

  std::vector<double> MassOffsetVote::wideErrors(const std::vector<double>& masses, const std::vector<double>& ladder)
  {
    std::vector<double> out;
    if (masses.empty() || ladder.empty()) return out;
    const long n = static_cast<long>(ladder.size());
    for (const double x : masses)
    {
      // j = min(max(bisect_left(tm, x), 1), n - 1); the neighbours j - 1 and j (a one-entry ladder: both are entry 0,
      // as Python's index -1 is).
      long j = static_cast<long>(std::lower_bound(ladder.begin(), ladder.end(), x) - ladder.begin());
      j = std::min(std::max(j, 1L), n - 1);
      const double a = ladder[static_cast<size_t>(j >= 1 ? j - 1 : n - 1)];
      const double b = ladder[static_cast<size_t>(j)];
      const double lo = (x - a) / a * 1e6;
      const double hi = (x - b) / b * 1e6;
      const double e = std::abs(lo) <= std::abs(hi) ? lo : hi;
      if (std::abs(e) <= kWindowPpm) out.push_back(e);
    }
    return out;
  }

  MassOffsetVoteResult MassOffsetVote::vote(const std::deque<std::vector<double>>& ring, double estimate_ppm)
  {
    MassOffsetVoteResult r;
    r.result = "too_few_errors";
    std::vector<double> all;
    for (const auto& lst : ring) all.insert(all.end(), lst.begin(), lst.end());
    if (static_cast<int>(all.size()) < kMinErrors) return r;

    // 1 ppm bins on [-300, 300] (600 bins; 300 falls in the last), then the 597 windows of four bins.
    // The bin is found by comparing with the integer edges (bisect_right - 1), never by flooring x + 300, which rounds
    // a tiny negative error into the bin above.
    constexpr int kBins = 600;
    std::vector<double> edges(kBins + 1);
    for (int k = 0; k <= kBins; ++k) edges[static_cast<size_t>(k)] = -300.0 + static_cast<double>(k);
    std::vector<int> counts(kBins, 0);
    for (const double x : all)
    {
      if (x < edges.front() || x > edges.back()) continue;
      const int i = (x == edges.back()) ? kBins - 1
                                        : static_cast<int>(std::upper_bound(edges.begin(), edges.end(), x) - edges.begin()) - 1;
      counts[static_cast<size_t>(i)] += 1;
    }
    std::vector<int> h(kBins - 3, 0);
    for (int j = 0; j < kBins - 3; ++j) h[static_cast<size_t>(j)] = counts[j] + counts[j + 1] + counts[j + 2] + counts[j + 3];
    int jbest = 0;
    for (int j = 1; j < kBins - 3; ++j)
      if (h[static_cast<size_t>(j)] > h[static_cast<size_t>(jbest)]) jbest = j;   // the first maximum
    const double c = -300.0 + static_cast<double>(jbest) + 2.0;

    std::vector<double> near;
    for (const double x : all)
      if (std::abs(x - c) <= kNearPpm) near.push_back(x);
    if (static_cast<int>(near.size()) < kMinErrors) return r;

    r.has_centre = true;
    r.centre = medianOf(near);
    r.errors = static_cast<int>(near.size());
    for (const auto& lst : ring)
      for (const double e : lst)
        if (std::abs(e - r.centre) <= kNearPpm)
        {
          ++r.scans;
          break;
        }
    for (int k = 0; k < kBins - 3; ++k)
      if (std::abs(-300.0 + static_cast<double>(k) + 2.0 - c) >= kRunnerUpPpm) r.runner_up = std::max(r.runner_up, h[static_cast<size_t>(k)]);

    if (r.scans < kMinScans) r.result = "too_few_scans";
    else if (r.errors < 2 * r.runner_up) r.result = "runner_up";
    else if (std::abs(r.centre - estimate_ppm) <= kJumpPpm) r.result = "confirms";
    else r.result = "jump";
    return r;
  }

  bool MassOffsetVote::step(int n_masses, bool confident, const std::vector<double>& wide, double estimate_ppm,
                            MassOffsetVoteResult& out)
  {
    if (confident) miss_ = 0;
    else if (n_masses >= kMinMasses) ++miss_;
    if (n_masses >= kMinMasses)
    {
      ring_.push_back(wide);
      if (static_cast<int>(ring_.size()) > kRingScans) ring_.pop_front();
    }
    if (miss_ >= kCounterLimit && static_cast<int>(ring_.size()) >= kRingScans)
    {
      miss_ = 0;
      out = vote(ring_, estimate_ppm);
      return true;
    }
    return false;
  }

  // --- MassOffsets ---------------------------------------------------------------------------------

  MassOffsets::MassOffsets(double ms2_start_ppm, double ms3_start_ppm, bool ms3_vote_armed)
    : ms2(ms2_start_ppm, kMassOffsetSigmaMs2Ppm), ms3(ms3_start_ppm, kMassOffsetSigmaMs3Ppm), vote_armed(ms3_vote_armed)
  {
  }

  namespace
  {
    void writeErrors(EventLine& ev, const std::vector<double>& errors)
    {
      ev << " raw_errors=";
      if (errors.empty()) ev << "none";
      for (size_t i = 0; i < errors.size(); ++i) ev << (i > 0 ? ";" : "") << std::fixed << std::setprecision(4) << errors[i];
    }
  } // namespace

  bool MassOffsets::onMs3(const std::string& tracking_id, double used_ppm, int n_masses, int matched, int bonds,
                          const std::vector<double>& raw_errors, const std::vector<double>& wide_errors)
  {
    const bool confident = bonds >= kMassOffsetConfidentMs3Bonds;
    if (confident) ms3.feed(raw_errors);
    MassOffsetVoteResult vr;
    const bool voted = vote_armed && vote.step(n_masses, confident, wide_errors, used_ppm, vr);
    if (voted && vr.result == "jump") ms3.restart(vr.centre);
    if (n_masses > 0)
    {
      EventLine ev;
      ev << "[MASS-OFFSET] level=3 tid=" << tracking_id << " used_ppm=" << std::fixed << std::setprecision(4) << used_ppm
         << " masses=" << n_masses << " matched=" << matched << " bonds=" << bonds << " confident=" << (confident ? 1 : 0)
         << " fed=" << (confident ? raw_errors.size() : 0) << " estimate_after=" << std::fixed << std::setprecision(4)
         << ms3.estimate() << " counter=" << (vote_armed ? vote.counter() : -1);
      writeErrors(ev, confident ? raw_errors : std::vector<double> {});
    }
    if (voted)
    {
      EventLine ev;
      ev << "[MASS-OFFSET-VOTE] level=3 tid=" << tracking_id << " result=" << vr.result << " errors=" << vr.errors
         << " scans=" << vr.scans << " runner_up=" << vr.runner_up << " from_ppm=" << std::fixed << std::setprecision(4)
         << used_ppm;
      if (vr.has_centre) ev << " centre_ppm=" << std::fixed << std::setprecision(4) << vr.centre;
      if (vr.result == "jump") ev << " to_ppm=" << std::fixed << std::setprecision(4) << ms3.estimate();
    }
    return confident;
  }

  bool MassOffsets::onMs2(const std::string& tracking_id, double used_ppm, int matched, const std::vector<double>& raw_errors)
  {
    const bool confident = matched >= kMassOffsetConfidentMs2Fragments;
    if (!confident) return false;
    ms2.feed(raw_errors);
    EventLine ev;
    ev << "[MASS-OFFSET] level=2 tid=" << tracking_id << " used_ppm=" << std::fixed << std::setprecision(4) << used_ppm
       << " matched=" << matched << " confident=1 fed=" << raw_errors.size() << " estimate_after=" << std::fixed
       << std::setprecision(4) << ms2.estimate();
    writeErrors(ev, raw_errors);
    return true;
  }

  void MassOffsets::logStart() const
  {
    EventLine() << "[MASS-OFFSET-START] level=2 start_ppm=" << std::fixed << std::setprecision(4) << ms2.start()
                << " sigma_ppm=" << kMassOffsetSigmaMs2Ppm << " vote=off";
    EventLine() << "[MASS-OFFSET-START] level=3 start_ppm=" << std::fixed << std::setprecision(4) << ms3.start()
                << " sigma_ppm=" << kMassOffsetSigmaMs3Ppm << " vote=" << (vote_armed ? "armed" : "off");
  }
} // namespace OpenMS
