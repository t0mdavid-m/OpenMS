// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms1Buffer.h>

#include <algorithm>
#include <cmath>

namespace OpenMS
{
  namespace
  {
    // SPEC 4, and features_v17.py of the package: the same constants, the same order of operations.
    constexpr double kProton = 1.007276;
    constexpr double kIso = 1.00235;
    constexpr double kLambdaPerDa = 0.000623;
    constexpr double kFloor = 0.02;
    constexpr double kPpm = 1e-5;
    constexpr int kZMin = 4;
    constexpr int kZMax = 50;
    constexpr int kMinIso = 3;
    constexpr double kMinCos = 0.7;
    constexpr double kTicScale = 1e9;
    constexpr double kHorizonS = 180.0;
    constexpr double kGapMs = 30000.0;
    constexpr double kTrendMs = 30000.0;

    /// The charges, expected envelope and m/z windows of one precursor mass within one scan range.
    struct Envelope
    {
      std::vector<double> p;            ///< kept isotopes' expected intensities, normalised to the tallest
      std::vector<int> zs;              ///< counting candidates: charges whose whole kept envelope lies in the range
      double pnorm = 0.0;
      std::vector<double> lo, hi;       ///< per (charge, isotope), charge-major: [target(1 - 1e-5), target(1 + 1e-5))

      Envelope(double mass, double range_lo, double range_hi)
      {
        const double lam = kLambdaPerDa * mass;
        const int n = static_cast<int>(lam * 3 + 12);
        std::vector<double> all(static_cast<size_t>(std::max(n, 0)));
        double pmax = 0.0;
        for (int k = 0; k < n; ++k)
        {
          all[static_cast<size_t>(k)] = std::exp(k * std::log(lam) - lam - std::lgamma(k + 1.0));
          pmax = std::max(pmax, all[static_cast<size_t>(k)]);
        }
        std::vector<int> ks;
        for (int k = 0; k < n; ++k)
        {
          if (all[static_cast<size_t>(k)] >= kFloor * pmax)
          {
            ks.push_back(k);
            p.push_back(all[static_cast<size_t>(k)] / pmax);
          }
        }
        if (ks.empty()) return;
        const int kmin = ks.front(), kmax = ks.back();
        for (int z = kZMin; z <= kZMax; ++z)
        {
          if ((mass + kmin * kIso + z * kProton) / z > range_lo && (mass + kmax * kIso + z * kProton) / z < range_hi) zs.push_back(z);
        }
        double ss = 0.0;
        for (double v : p) ss += v * v;
        pnorm = std::sqrt(ss);
        for (int z : zs)
        {
          for (int k : ks)
          {
            const double target = (mass + k * kIso + z * kProton) / z;
            lo.push_back(target * (1 - kPpm));
            hi.push_back(target * (1 + kPpm));
          }
        }
      }
    };

    /// The summed found intensity of the counting charges in one survey (features_v17.profile_value).
    double profileValue(const std::vector<double>& mz, const std::vector<double>& it, const Envelope& env)
    {
      if (env.zs.empty() || mz.empty()) return 0.0;
      const size_t nk = env.p.size();
      double value = 0.0;
      std::vector<double> obs(nk);
      for (size_t zi = 0; zi < env.zs.size(); ++zi)
      {
        int found = 0;
        for (size_t k = 0; k < nk; ++k)
        {
          const size_t w = zi * nk + k;
          // Both ends as a lower-bound search: a centroid exactly at the upper end is outside (SPEC 4).
          const auto a = std::lower_bound(mz.begin(), mz.end(), env.lo[w]);
          const auto b = std::lower_bound(a, mz.end(), env.hi[w]);
          double best = 0.0;
          for (auto i = a; i != b; ++i) best = std::max(best, it[static_cast<size_t>(i - mz.begin())]);
          obs[k] = best;
          if (best > 0) ++found;
        }
        if (found < kMinIso) continue;
        double dot = 0.0, oo = 0.0;
        for (size_t k = 0; k < nk; ++k)
        {
          dot += obs[k] * env.p[k];
          oo += obs[k] * obs[k];
        }
        if (dot / (std::sqrt(oo) * env.pnorm) >= kMinCos)
        {
          for (size_t k = 0; k < nk; ++k) value += obs[k];
        }
      }
      return value;
    }

    double adjusted(const BufferedSurvey& s, const Envelope& env)
    {
      const double v = profileValue(s.mz, s.intensity, env);
      return (v > 0 && s.ion_current > 0) ? kTicScale * v / s.ion_current : 0.0;
    }
  } // namespace

  void Ms1Buffer::add(BufferedSurvey s)
  {
    kind_by_id_[s.tracking_id] = s.kind;
    const uint64_t newest = s.resolve_ts;
    surveys_.push_back(std::move(s));
    while (!surveys_.empty() && surveys_.front().resolve_ts + kHorizonMs < newest) surveys_.pop_front();
  }

  Ms1Buffer::Selection Ms1Buffer::select(const std::string& parent_survey_id, uint64_t t) const
  {
    Selection sel;
    const SurveyKind* ref = nullptr;
    auto it = kind_by_id_.find(parent_survey_id);
    if (!parent_survey_id.empty() && it != kind_by_id_.end())
    {
      ref = &it->second;
      sel.kind_from = "parent survey";
    }
    else
    {
      for (auto r = surveys_.rbegin(); r != surveys_.rend(); ++r)
      {
        if (r->resolve_ts <= t)
        {
          ref = &r->kind;
          sel.kind_from = "last survey";
          break;
        }
      }
    }
    if (ref == nullptr) return sel;   // no survey resolved by t: an empty selection, the range (500, 1800)
    if (!ref->first_mass.empty() && !ref->last_mass.empty())
    {
      sel.lo = std::stod(ref->first_mass);
      sel.hi = std::stod(ref->last_mass);
    }
    for (const BufferedSurvey& s : surveys_)
    {
      if (s.resolve_ts > t) break;
      if (s.resolve_ts + kHorizonMs < t) continue;
      if (s.kind != *ref) continue;
      sel.surveys.push_back(&s);
      sel.centroids += s.centroid_count;
    }
    return sel;
  }

  double adjustedPrecursorValue(const BufferedSurvey& s, double mass, double lo, double hi)
  {
    return adjusted(s, Envelope(mass, lo, hi));
  }

  PeakFeatures peakFeatures(const std::vector<double>& ts, const std::vector<double>& inten, double t)
  {
    std::vector<size_t> seen;
    for (size_t j = 0; j < inten.size(); ++j)
      if (inten[j] > 0) seen.push_back(j);
    if (seen.empty()) return {kHorizonS, kHorizonS, kHorizonS, static_cast<double>(ts.size()), 0.0, 0.0, 0.0, 0.0};
    // The current elution: walking back from the last sighting, it ends before the first pair of consecutive sightings
    // more than 30 s apart with at least one survey between them.
    size_t k = seen.size() - 1;
    while (k > 0 && !(ts[seen[k]] - ts[seen[k - 1]] > kGapMs && seen[k] - seen[k - 1] > 1)) --k;
    const std::vector<size_t> e(seen.begin() + static_cast<std::ptrdiff_t>(k), seen.end());
    size_t apex = e.front();   // the first of equal maxima
    for (size_t j : e)
      if (inten[j] > inten[apex]) apex = j;
    const size_t last = e.back();
    std::vector<size_t> w;
    for (size_t j : e)
      if (ts[j] >= ts[last] - kTrendMs) w.push_back(j);
    double trend = 0.0;
    if (w.size() >= 2)
    {
      // Sequential sums, as the reference's Python sum().
      double xs = 0.0;
      for (size_t j : w) xs += ts[j];
      const double xm = xs / static_cast<double>(w.size());
      double sxx = 0.0;
      for (size_t j : w) sxx += (ts[j] - xm) * (ts[j] - xm);
      if (sxx > 0)
      {
        std::vector<double> ly;
        for (size_t j : w) ly.push_back(std::log10(inten[j]));
        double ys = 0.0;
        for (double v : ly) ys += v;
        const double ym = ys / static_cast<double>(ly.size());
        double sxy = 0.0;
        for (size_t q = 0; q < w.size(); ++q) sxy += (ts[w[q]] - xm) * (ly[q] - ym);
        trend = 1e4 * sxy / sxx;
      }
    }
    return {(t - ts[e.front()]) / 1e3,
            (t - ts[apex]) / 1e3,
            (t - ts[last]) / 1e3,
            static_cast<double>(ts.size() - 1 - last),
            inten[last] / inten[apex],
            std::log10(std::max(inten[apex], 1.0)),
            trend,
            static_cast<double>(e.size())};
  }

  PeakFeatures planPeakFeatures(const Ms1Buffer::Selection& selection, double mass, uint64_t t)
  {
    const Envelope env(mass, selection.lo, selection.hi);
    std::vector<double> ts, vals;
    ts.reserve(selection.surveys.size());
    vals.reserve(selection.surveys.size());
    for (const BufferedSurvey* s : selection.surveys)
    {
      ts.push_back(static_cast<double>(s->resolve_ts));
      vals.push_back(adjusted(*s, env));
    }
    return peakFeatures(ts, vals, static_cast<double>(t));
  }
} // namespace OpenMS
