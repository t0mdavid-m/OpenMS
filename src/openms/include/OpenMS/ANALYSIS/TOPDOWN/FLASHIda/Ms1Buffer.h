// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <array>
#include <cstdint>
#include <deque>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{
  /**
    @brief The kind of an MS1 scan, as scan_commands.tsv logs its command (ADR-0079; package 17, SPEC 3): scan_type,
    faims_cv, faims_enabled, first_mass, last_mass -- the five cells, compared as the strings the log prints.
  */
  struct OPENMS_DLLAPI SurveyKind
  {
    std::string scan_type, faims_cv, faims_enabled, first_mass, last_mass;
    bool operator==(const SurveyKind& o) const
    {
      return scan_type == o.scan_type && faims_cv == o.faims_cv && faims_enabled == o.faims_enabled
             && first_mass == o.first_mass && last_mass == o.last_mass;
    }
    bool operator!=(const SurveyKind& o) const { return !(*this == o); }
  };

  /// One resolved MS1 scan as the W gate keeps it: its centroids (m/z ascending, intensity > 0, 64-bit) and its ion
  /// current (the sum of all its centroid intensities, 64-bit).
  struct OPENMS_DLLAPI BufferedSurvey
  {
    std::string tracking_id;
    uint64_t resolve_ts = 0;   ///< ms of the engine clock, the scan_results.tsv resolve_ts of its row
    SurveyKind kind;
    std::vector<double> mz, intensity;
    double ion_current = 0.0;
    /// The number of centroids of the full scan (buffer_centroids of the log). Live it is mz.size(); a test that feeds
    /// pruned centroids sets the full count.
    uint64_t centroid_count = 0;
  };

  /**
    @brief The MS1 buffer of the W gate (ADR-0079; package 17, SPEC 3): every MS1 the engine resolved in the last 180 s.

    The only state the gate adds. Filled under analysis_mutex_ when an MS1 row is written, read under it when an MS3 is
    planned. A plan at time t selects the surveys with resolve_ts in [t - 180000, t] (both ends included) whose kind is
    the kind of the survey the parent MS2 was picked from, or, when that survey is unknown, of the last survey resolved
    at or before t.
  */
  class OPENMS_DLLAPI Ms1Buffer
  {
  public:
    static constexpr uint64_t kHorizonMs = 180000;

    /// Append @p s (resolve_ts must not decrease) and drop the surveys more than the horizon older than it.
    void add(BufferedSurvey s);

    struct Selection
    {
      std::vector<const BufferedSurvey*> surveys;   ///< in resolve_ts order
      double lo = 500.0, hi = 1800.0;               ///< the scan range of the kind; (500, 1800) when it has none
      std::string kind_from = "none";                ///< "parent survey" | "last survey" | "none"
      uint64_t centroids = 0;                        ///< buffer_centroids: the full-scan centroid counts, summed
    };

    /// The surveys a plan at @p t uses; @p parent_survey_id is the tracking id the parent MS2 was picked from.
    Selection select(const std::string& parent_survey_id, uint64_t t) const;

    size_t size() const { return surveys_.size(); }

  private:
    std::deque<BufferedSurvey> surveys_;
    /// The kind of every survey ever added, by tracking id (a few bytes per survey): the reference looks the parent's
    /// survey up over the whole run, also when it has left the horizon.
    std::unordered_map<std::string, SurveyKind> kind_by_id_;
  };

  /// The eight peak features (inputs 30-37; SPEC 5), 64-bit, in the order pk3_since_start, pk3_since_apex,
  /// pk3_since_last, pk3_missed, pk3_rel_max, pk3_log_max, pk3_trend, pk3_n_seen.
  using PeakFeatures = std::array<double, 8>;

  /// The adjusted value of the precursor of mass @p mass in one survey (SPEC 4): 1e9 x the summed found intensity of the
  /// counting charges / the survey's ion current; 0 when either is 0. @p lo, @p hi: the scan range of the selection.
  OPENMS_DLLAPI double adjustedPrecursorValue(const BufferedSurvey& s, double mass, double lo, double hi);

  /// The eight features from the resolve times (ms, ascending) and adjusted values of the selected surveys (SPEC 5).
  OPENMS_DLLAPI PeakFeatures peakFeatures(const std::vector<double>& resolve_ts, const std::vector<double>& adjusted, double t);

  /// The eight features of a plan: @p selection, the precursor mass @p mass (stage 0 as logged), the plan time @p t.
  OPENMS_DLLAPI PeakFeatures planPeakFeatures(const Ms1Buffer::Selection& selection, double mass, uint64_t t);
} // namespace OpenMS
