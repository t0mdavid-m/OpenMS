// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qvalue.h>

namespace OpenMS
{
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

  double Qvalue::updatePeakGroupQvalues(std::vector<DeconvolvedSpectrum>& deconvolved_spectra) // per ms level + precursor update as well.
  {
    double noise_weight = 1;
    std::map<uint, std::vector<double>> tscore_map; // per ms level
    std::map<uint, std::vector<double>> dscore_iso_decoy_map;
    std::map<uint, std::vector<double>> dscore_noise_decoy_map;
    std::map<uint, std::vector<double>> dscore_charge_decoy_map;
    std::map<uint, std::map<double, double>> qscore_qvalue_map; //

    // to calculate qvalues per ms level, store Qscores per ms level
    std::set<uint> used_feature_indices;

    for (auto& deconvolved_spectrum : deconvolved_spectra)
    {
      if (deconvolved_spectrum.empty() || deconvolved_spectrum.isDecoy())
        continue;

      uint ms_level = deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      for (auto& pg : deconvolved_spectrum)
      {
        if (pg.getFeatureIndex() > 0 && used_feature_indices.find(pg.getFeatureIndex()) != used_feature_indices.end())
          continue;
        used_feature_indices.insert(pg.getFeatureIndex());
        tscore_map[ms_level].push_back(pg.getQscore2D());
      }
    }

    for (auto& decoy_deconvolved_spectrum : deconvolved_spectra)
    {
      if (decoy_deconvolved_spectrum.empty() || !decoy_deconvolved_spectrum.isDecoy())
        continue;

      uint ms_level = decoy_deconvolved_spectrum.getOriginalSpectrum().getMSLevel();
      for (auto& pg : decoy_deconvolved_spectrum)
      {
        if (pg.getFeatureIndex() > 0 && used_feature_indices.find(pg.getFeatureIndex()) != used_feature_indices.end())
          continue;
        used_feature_indices.insert(pg.getFeatureIndex());
        if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::isotope_decoy)
        {
          dscore_iso_decoy_map[ms_level].push_back(pg.getQscore2D());
        }
        else if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::noise_decoy)
        {
          dscore_noise_decoy_map[ms_level].push_back(pg.getQscore2D());
        }
        else if (pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::charge_decoy)
        {
          dscore_charge_decoy_map[ms_level].push_back(pg.getQscore2D());
        }
      }
    }

    for (auto& [ms_level, qscores] : tscore_map)
    {
      auto& dscore_iso = dscore_iso_decoy_map[ms_level];
      auto& dscore_charge = dscore_charge_decoy_map[ms_level];
      auto& dscore_noise = dscore_noise_decoy_map[ms_level];

      std::sort(qscores.begin(), qscores.end());
      std::sort(dscore_iso.begin(), dscore_iso.end());
      std::sort(dscore_noise.begin(), dscore_noise.end());
      std::sort(dscore_charge.begin(), dscore_charge.end());

      double sum = 0;
      double max_score_for_weight_calculation = .5;
      double min_score_for_weight_calculation = .2;
      double iso_sum = std::accumulate(dscore_iso.begin(), dscore_iso.end(), .0);
      double noise_sum = std::accumulate(dscore_noise.begin(), dscore_noise.end(), .0);

      for (int i = dscore_iso.size() - 1; i >= 0; i--)
      {
        sum += dscore_iso[i];
        if (sum > iso_sum * .8 || dscore_iso[i] < .5)
        {
          max_score_for_weight_calculation = std::min(max_score_for_weight_calculation, dscore_iso[i]);
          break;
        }
      }

//      sum = 0;
//      for (double i : dscore_noise)
//      {
//        sum += i;
//        if (sum > noise_sum * .1 || i > .2)
//        {
//          min_score_for_weight_calculation = std::max(min_score_for_weight_calculation, i);
//          break;
//        }
//      }

      Size num_bin = 5;//std::min(Size(5), 1 + dscore_noise.size()/100);
      auto qscore_vec = getDistVector(qscores, num_bin, min_score_for_weight_calculation, max_score_for_weight_calculation);
      auto qscore_charge_vec = getDistVector(dscore_charge, num_bin, min_score_for_weight_calculation, max_score_for_weight_calculation);
      auto qscore_noise_vec = getDistVector(dscore_noise, num_bin, min_score_for_weight_calculation, max_score_for_weight_calculation);
      auto qscore_iso_vec = getDistVector(dscore_iso, num_bin, min_score_for_weight_calculation, max_score_for_weight_calculation);

      Matrix<double> left( qscore_vec.rows(), 2, 1);
      for (int r = 0; r < qscore_vec.rows(); r++)
      {
        double v = qscore_vec.getValue(r, 0);
        v -= qscore_iso_vec.getValue(r, 0) + qscore_charge_vec.getValue(r, 0);
        v = std::max(v, .0);
        qscore_vec.setValue(r, 0, v);
        left.setValue(r, 0, qscore_noise_vec.getValue(r, 0));
      }

      auto calculated_vec = left.completeOrthogonalDecomposition().pseudoInverse() * qscore_vec;
      noise_weight = calculated_vec.row(0)[0];

      if (calculated_vec.row(1)[0] < 0) // &&
      {
        auto calculated_vec_non_negative = qscore_noise_vec.completeOrthogonalDecomposition().pseudoInverse() * qscore_vec;
        noise_weight = calculated_vec_non_negative.row(0)[0];
      }

      if (std::isnan(noise_weight)) noise_weight = 1.0;
      noise_weight = std::max(noise_weight, 0.01);
      //noise_weight = (qscore_noise_vec.completeOrthogonalDecomposition().pseudoInverse() * qscore_vec);

      //std::cout<<qscore_noise_vec.sum() << " vs " << qscore_vec.sum() << "\n" <<  calculated_vec << " vs " << (qscore_noise_vec.completeOrthogonalDecomposition().pseudoInverse() * qscore_vec) << " to " << noise_weight << std::endl;
      //noise_weight = 2;
      /*
      double a = 0, b = 0;
      for (double i : qscores)
      {
        if (i < min_score_for_weight_calculation)
          continue;
        if (i > max_score_for_weight_calculation)
          break;
        b++;
      }
      //std::cout<< min_score_for_weight_calculation << " to " << max_score_for_weight_calculation << " total " << b;
      for (double i : dscore_charge)
      {
        if (i < min_score_for_weight_calculation)
          continue;
        if (i > max_score_for_weight_calculation)
          break;
        b--;
      }
      //std::cout<< " charge removed " << b;
      for (double i : dscore_iso)
      {
        if (i < min_score_for_weight_calculation)
          continue;
        if (i > max_score_for_weight_calculation)
          break;
        b--;
      }
      //std::cout<< " iso removed " << b;
      for (double i : dscore_noise)
      {
        if (i < min_score_for_weight_calculation)
          continue;
        if (i > max_score_for_weight_calculation)
          break;
        a++;
      }
      */
      //std::cout<< " total noise " << a << " ratio " << b/a << std::endl;
      std::sort(qscores.rbegin(), qscores.rend());
      std::sort(dscore_iso.rbegin(), dscore_iso.rend());
      std::sort(dscore_noise.rbegin(), dscore_noise.rend());
      std::sort(dscore_charge.rbegin(), dscore_charge.rend());


      //noise_weight = noise_weight * noise_weight;
      auto& map_qvalue = qscore_qvalue_map[ms_level];
      double nom_i = 0, nom_c = 0, nom_n = 0;
      Size j_i = 0, j_c = 0, j_n = 0;

      for (Size i = 0; i < qscores.size(); i++)
      {
        double ts = qscores[i];
        double di = 0, dc = 0, dn = 0;
        while (i < qscores.size() - 1 && qscores[i + 1] == ts)
        {
          i++;
        }

        while (j_n < dscore_noise.size() && dscore_noise[j_n] >= ts)
        {
          dn += noise_weight;
          ++j_n;
        }
        while (j_i < dscore_iso.size() && dscore_iso[j_i] >= ts)
        {
          di++;
          ++j_i;
        }
        while (j_c < dscore_charge.size() && dscore_charge[j_c] >= ts)
        {
          dc++;
          ++j_c;
        }
        nom_n += dn;
        nom_i += di;
        nom_c += dc;
        double tmp_q = (nom_i + nom_c + nom_n) / double(1 + i);
        map_qvalue[ts] = std::min(1.0, tmp_q);
      }
    }

    for (auto& titem : tscore_map)
    {
      uint ms_level = titem.first;
      auto& map_qvalue = qscore_qvalue_map[ms_level];

      double cummin = 1.0;
      {
        for (auto&& rit = map_qvalue.begin(); rit != map_qvalue.end(); ++rit)
        {
          cummin = std::min(rit->second, cummin);
          rit->second = cummin;
        }
      }

      for (auto& deconvolved_spectrum : deconvolved_spectra)
      {
        if (deconvolved_spectrum.empty() || deconvolved_spectrum.isDecoy())
          continue;

        if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() == ms_level + 1 && !deconvolved_spectrum.getPrecursorPeakGroup().empty())
        {
          auto precursor_pg = deconvolved_spectrum.getPrecursorPeakGroup();
          double qs = precursor_pg.getQscore2D();

          precursor_pg.setQvalue(map_qvalue[qs]);
          deconvolved_spectrum.setPrecursorPeakGroup(precursor_pg);
        }

        if (deconvolved_spectrum.getOriginalSpectrum().getMSLevel() != ms_level)
        {
          continue;
        }

        for (auto& pg : deconvolved_spectrum)
        {
          pg.setQvalue(map_qvalue[pg.getQscore2D()]);
        }
      }
    }
    return noise_weight;
  }
} // namespace OpenMS