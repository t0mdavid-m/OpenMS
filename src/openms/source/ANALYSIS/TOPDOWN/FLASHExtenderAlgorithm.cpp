// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <utility>

namespace OpenMS
{
  inline const Size max_node_cntr = 500;

  FLASHExtenderAlgorithm::FLASHExtenderAlgorithm(): DefaultParamHandler("FLASHExtenderAlgorithm"), ProgressLogger()
  {
    setDefaultParams_();
  }

  FLASHExtenderAlgorithm& FLASHExtenderAlgorithm::operator=(const FLASHExtenderAlgorithm& rhs)
  {
    if (this == &rhs) return *this;

    DefaultParamHandler::operator=(rhs);
    return *this;
  }

  void FLASHExtenderAlgorithm::setDefaultParams_()
  {
    defaults_.setValue("flanking_mass_tol", 700.0, "Flanking mass tolerance in Da.");

    defaults_.setValue("fdr", 1.0, "Protein FDR threshold.");
    defaults_.setMaxFloat("fdr", 1.0);
    defaults_.setMinFloat("fdr", 0.01);

    defaults_.setValue("keep_decoy", "false", "Keep decoy proteins.");
    defaults_.addTag("keep_decoy", "advanced");
    defaults_.setValidStrings("keep_decoy", {"true", "false"});

    defaultsToParam_();
  }

  void FLASHExtenderAlgorithm::updateMembers_()
  {
  }


  int FLASHExtenderAlgorithm::getVertex_(int peak_index, int pro_index, int score, int num_mod, int layer) const
  {
    return (((peak_index * (pro_length_ + 1) + pro_index) * (max_mod_cntr_ + 1) + num_mod) * (max_path_score_ - min_path_score_ + 1)
           + (score - min_path_score_)) * (max_layer_ + 1) + layer;
  }

  int FLASHExtenderAlgorithm::getPeakIndex_(int vertex) const
  {
    return (vertex / (max_layer_ + 1) / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1)) / (pro_length_ + 1);
  }

  int FLASHExtenderAlgorithm::getProIndex_(int vertex) const
  {
    return ((vertex / (max_layer_ + 1) / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1))) % (pro_length_ + 1);
  }

  int FLASHExtenderAlgorithm::getScore_(int vertex) const
  {
    return 0;
  }

  int FLASHExtenderAlgorithm::getModNumber_(int vertex) const
  {
    return 0;
  }

  void FLASHExtenderAlgorithm::constructSubDAG_(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited,
                                           int vertex1, int peak_index2, int pro_index2, int layer, bool allow_truncation)
  {
    if (!visited[vertex1]) return;
    int peak_index1 = getPeakIndex_(vertex1);
    int score1 = getScore_(vertex1);
    int num_mod1 = getModNumber_(vertex1);
    int pro_index1 = getProIndex_(vertex1);

    if (peak_index1 == peak_index2 && pro_index1 == pro_index2)
    {
      for (int score = score1 - 1; score >= min_path_score_; score--)
      {
        int vertex = getVertex_(peak_index1, pro_index1, score, num_mod1, layer);
         // if lower scoring path exist, delete. keep the best score per num_mod
        visited[vertex] = false;
      }
      return;
    }

    if (allow_truncation)
    {
      for (int pro_index = pro_index1 + 1; pro_index <= pro_index2 ; pro_index++)
      {
        int vertex2 = getVertex_(peak_index1, pro_index, score1, num_mod1, layer);
        if (visited[vertex2]) return;
        visited[vertex2] = true;
        constructSubDAG_(dag, visited, vertex2, peak_index2, pro_index2, layer);
      }
      return;
    }

    double delta_mass1 = peak_masses_[peak_index1] - pro_masses_[pro_index1];

    for (int peak_index = peak_index1 + 1; peak_index <= peak_index2 ; peak_index++)
    {
      int score = score1 + peak_scores_[peak_index];
      for (int pro_index = pro_index1 + 1; pro_index <= pro_index2; pro_index++)
      {
        double delta_mass = peak_masses_[peak_index] - pro_masses_[pro_index];
        if (delta_mass - delta_mass1 > max_mod_mass_) break;
        if (delta_mass1 - delta_mass > max_mod_mass_) continue;
        int num_mod = num_mod1;
        if (std::abs(delta_mass - delta_mass1) > tol_ * peak_masses_[peak_index1]) num_mod ++;
        if (num_mod > max_mod_cntr_) continue;
        if (peak_index1 == peak_index2 && pro_index1 == pro_index2) layer = 0; // for sink, go to the layer 0
        int vertex2 = getVertex_(peak_index, pro_index, score, num_mod, layer);
        if (visited[vertex2]) return;
        visited[vertex2] = true;
        constructSubDAG_(dag, visited, vertex2, peak_index2, pro_index2, layer);
      }
    }

  }

} // namespace OpenMS