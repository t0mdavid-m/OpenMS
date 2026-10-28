// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qvalue.h>
#include <algorithm>
#include <cmath>

namespace OpenMS
{
  /// simple function to generate distribution from input vector values
  Matrix<double> getDistVector(const std::vector<double>& values, Size num_bin, double minv, double maxv)
  {
    Matrix<double> ret(num_bin, 1, .0);
    for (const auto& v : values)
    {
      if (v >= maxv) continue;
      if (v < minv) continue;

      Size bin = Size((v - minv) / (maxv - minv) * num_bin);
      ret.setValue(bin, 0, ret.getValue(bin, 0) + 1);
    }
    return ret;
  }

  std::map<uint, double> Qvalue::updatePeakGroupQvalues(std::vector<DeconvolvedSpectrum>& deconvolved_spectra) // per ms level + precursor update as well.
  {
    std::map<uint, double> noise_weights; // per MS level (the weight of the last level was returned for all levels)
    std::map<uint, std::vector<double>> score_map_target; // target PeakGroupScoring vector per ms level
    std::map<uint, std::vector<double>> score_signal_decoy_map; // signal decoy PeakGroupScoring vector per ms level
    std::map<uint, std::vector<double>> score_noise_decoy_map; // noise decoy PeakGroupScoring vector per ms level
    // qscore -> qvalue per ms level as two arrays, filled with the unique qscores in descending order (a std::map built from them in this
    // insertion order is the former table)
    std::map<uint, std::pair<std::vector<double>, std::vector<double>>> qscore_qvalue_arrays;

    // to calculate qvalues per ms level, store Qscores per ms level
    std::vector<bool> used_feature_indices; // the feature indices seen so far, by index

    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      if (deconvolved_spectrum.empty())
        continue;

      uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      for (auto& pg : deconvolved_spectrum)
      {
        const uint feature_index = pg.getFeatureIndex();
        if (feature_index > 0 && feature_index < used_feature_indices.size() && used_feature_indices[feature_index])
          continue;
        if (feature_index >= used_feature_indices.size()) used_feature_indices.resize((Size)feature_index + 1, false);
        used_feature_indices[feature_index] = true;

        if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target)
        {
          score_map_target[ms_level].push_back(pg.getQscore2D());
        }
        else if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::signal_decoy)
        {
          score_signal_decoy_map[ms_level].push_back(pg.getQscore2D());
        }
        else if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::noise_decoy)
        {
          score_noise_decoy_map[ms_level].push_back(pg.getQscore2D());
        }
      }
    }

    // per ms score, calculate Qvalues
    for (auto& [ms_level, scores_target] : score_map_target)
    {
      auto& scores_signal_decoy = score_signal_decoy_map[ms_level];
      auto& scores_noise_decoy = score_noise_decoy_map[ms_level];

      std::sort(scores_target.begin(), scores_target.end());
      std::sort(scores_signal_decoy.begin(), scores_signal_decoy.end());
      std::sort(scores_noise_decoy.begin(), scores_noise_decoy.end());

      double sum = 0;
      double max_score_for_weight_calculation = .7;
      double min_score_for_weight_calculation = .3;
      double iso_sum = std::accumulate(scores_signal_decoy.begin(), scores_signal_decoy.end(), .0);

      for (int i = scores_signal_decoy.size() - 1; i >= 0; i--)
      {
        sum += scores_signal_decoy[i];
        if (sum > iso_sum * .8 || scores_signal_decoy[i] < .5)
        {
          max_score_for_weight_calculation = std::min(max_score_for_weight_calculation, scores_signal_decoy[i]);
          break;
        }
      }

      Size num_bin = 6;
      // get the score distributions
      auto score_dist_target = getDistVector(scores_target, num_bin, min_score_for_weight_calculation, max_score_for_weight_calculation);
      auto score_dist_noise_decoy = getDistVector(scores_noise_decoy, num_bin, min_score_for_weight_calculation, max_score_for_weight_calculation);
      auto score_dist_signal_decoy = getDistVector(scores_signal_decoy, num_bin, min_score_for_weight_calculation, max_score_for_weight_calculation);

      // noise decoy weight calculation using Least Square
      Matrix<double> left(score_dist_target.rows(), 2, 1);
      for (int r = 0; r < score_dist_target.rows(); r++)
      {
        double v = score_dist_target.getValue(r, 0);
        v -= score_dist_signal_decoy.getValue(r, 0);
        score_dist_target.setValue(r, 0, v);
        left.setValue(r, 0, score_dist_noise_decoy.getValue(r, 0));
      }

      auto calculated_vec = left.completeOrthogonalDecomposition().pseudoInverse() * score_dist_target;
      double noise_weight = calculated_vec.row(0)[0];

      if (calculated_vec.row(1)[0] < 0)
      {
        auto calculated_vec_non_negative = score_dist_noise_decoy.completeOrthogonalDecomposition().pseudoInverse() * score_dist_target;
        noise_weight = calculated_vec_non_negative.row(0)[0];
      }

      // an empty fit window (no noise decoy between the bounds, e.g. the top signal decoy below 0.3) leaves no evidence: weight 1
      if (std::isnan(noise_weight) || score_dist_noise_decoy.sum() == 0) noise_weight = 1.0;
      noise_weight = std::max(noise_weight, 0.01);
      noise_weights[ms_level] = noise_weight;
      std::sort(scores_target.rbegin(), scores_target.rend());
      std::sort(scores_signal_decoy.rbegin(), scores_signal_decoy.rend());
      std::sort(scores_noise_decoy.rbegin(), scores_noise_decoy.rend());

      // now get the qvalues
      auto& [qscore_keys, qvalue_values] = qscore_qvalue_arrays[ms_level];
      qscore_keys.reserve(scores_target.size());
      qvalue_values.reserve(scores_target.size());
      double nom_i = 0, nom_c = 0, nom_n = 0;
      Size j_i = 0, j_n = 0;

      for (Size i = 0; i < scores_target.size(); i++)
      {
        double ts = scores_target[i];
        double di = 0, dc = 0, dn = 0;
        while (i < scores_target.size() - 1 && scores_target[i + 1] == ts)
        {
          i++;
        }

        while (j_n < scores_noise_decoy.size() && scores_noise_decoy[j_n] >= ts)
        {
          dn += noise_weight;
          ++j_n;
        }
        while (j_i < scores_signal_decoy.size() && scores_signal_decoy[j_i] >= ts)
        {
          di++;
          ++j_i;
        }
        nom_n += dn;
        nom_i += di;
        nom_c += dc;
        double tmp_q = (nom_i + nom_c + nom_n) / double(1 + i);
        qscore_keys.push_back(ts);
        qvalue_values.push_back(std::min(1.0, tmp_q));
      }
    }

    // refine qvalues to make them monotonic decreasing
    for (const auto& titem : score_map_target)
    {
      uint ms_level = titem.first;
      auto& [qscore_keys, qvalue_values] = qscore_qvalue_arrays[ms_level];
      bool key_nan = false;
      for (double qs : qscore_keys) key_nan = key_nan || std::isnan(qs);
      // the former table, built by the same insertions (needed only with a NaN among the keys or the qscores looked up below)
      std::map<double, double> map_qvalue;

      double cummin = 1.0;
      if (key_nan)
      {
        for (Size k = 0; k < qscore_keys.size(); k++) map_qvalue[qscore_keys[k]] = qvalue_values[k];
        for (auto&& rit = map_qvalue.begin(); rit != map_qvalue.end(); ++rit)
        {
          cummin = std::min(rit->second, cummin);
          rit->second = cummin;
        }
      }
      else
      {
        // unique keys in descending order: reversed, they are in the ascending order of the map's iteration
        std::reverse(qscore_keys.begin(), qscore_keys.end());
        std::reverse(qvalue_values.begin(), qvalue_values.end());
        for (auto& qvalue : qvalue_values)
        {
          cummin = std::min(qvalue, cummin);
          qvalue = cummin;
        }
      }

      // The q-values below are looked up without inserting, the spectra in parallel: std::map::operator[] gave a Qscore2D without a q-value
      // (of a decoy mass, or of a mass of a feature counted once) 0 and inserted it, which a later lookup of that Qscore2D found as 0 again.
      // That holds for ordered keys; with a NaN among the keys or the Qscore2D looked up, the serial loop with operator[] stays.
      bool has_nan = key_nan;
#pragma omp parallel for schedule(dynamic, 64) reduction(|| : has_nan)
      for (SignedSize s = 0; s < (SignedSize)deconvolved_spectra.size(); s++)
      {
        const auto& deconvolved_spectrum = deconvolved_spectra[s];
        if (deconvolved_spectrum.empty() || deconvolved_spectrum.isDecoy())
          continue;
        if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() == ms_level + 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
        {
          has_nan = has_nan || std::isnan(deconvolved_spectrum.getPrecursorPeakGroup().getQscore2D());
        }
        if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != ms_level)
        {
          continue;
        }
        for (const auto& pg : deconvolved_spectrum)
        {
          has_nan = has_nan || std::isnan(pg.getQscore2D());
        }
      }

      if (! has_nan)
      {
        // binary search in the ascending keys, as std::map::find
        // the q-value of the lowest target score at or above qs (the same targets are counted at both thresholds), the highest target
        // score's above all of them; it was 0 for a score that no target has (a decoy mass, a precursor that is no target mass, a mass of
        // a feature counted once), which let such masses pass every q-value cut
        const auto qvalue_of = [&qscore_keys, &qvalue_values](double qs) {
          if (qscore_keys.empty()) return 1.0;
          const auto it = std::lower_bound(qscore_keys.begin(), qscore_keys.end(), qs);
          return it == qscore_keys.end() ? qvalue_values.back() : qvalue_values[it - qscore_keys.begin()];
        };
#pragma omp parallel for schedule(dynamic, 64)
        for (SignedSize s = 0; s < (SignedSize)deconvolved_spectra.size(); s++)
        {
          auto& deconvolved_spectrum = deconvolved_spectra[s];
          if (deconvolved_spectrum.empty() || deconvolved_spectrum.isDecoy())
            continue;

          // set precursor Qvalue here
          if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() == ms_level + 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
          {
            auto precursor_pg = deconvolved_spectrum.getPrecursorPeakGroup();
            double qs = precursor_pg.getQscore2D();

            precursor_pg.setQvalue(qvalue_of(qs));
            deconvolved_spectrum.setPrecursorPeakGroup(precursor_pg);
          }

          if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != ms_level)
          {
            continue;
          }

          for (auto& pg : deconvolved_spectrum)
          {
            pg.setQvalue(qvalue_of(pg.getQscore2D()));
          }
        }
        continue;
      }

      if (! key_nan)
      {
        // the refined table as a map, for the serial lookups below
        for (Size k = 0; k < qscore_keys.size(); k++) map_qvalue.emplace_hint(map_qvalue.end(), qscore_keys[k], qvalue_values[k]);
      }
      // as qvalue_of above (a NaN Qscore2D: q-value 1)
      const auto map_qvalue_of = [&map_qvalue](double qs) {
        if (map_qvalue.empty() || std::isnan(qs)) return 1.0;
        const auto it = map_qvalue.lower_bound(qs);
        return it == map_qvalue.end() ? map_qvalue.rbegin()->second : it->second;
      };
      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        if (deconvolved_spectrum.empty() || deconvolved_spectrum.isDecoy())
          continue;

        // set precursor Qvalue here
        if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() == ms_level + 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
        {
          auto precursor_pg = deconvolved_spectrum.getPrecursorPeakGroup();
          double qs = precursor_pg.getQscore2D();

          precursor_pg.setQvalue(map_qvalue_of(qs));
          deconvolved_spectrum.setPrecursorPeakGroup(precursor_pg);
        }

        if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != ms_level)
        {
          continue;
        }

        for (auto& pg : deconvolved_spectrum)
        {
          pg.setQvalue(map_qvalue_of(pg.getQscore2D()));
        }
      }
    }
    return noise_weights;
  }
} // namespace OpenMS