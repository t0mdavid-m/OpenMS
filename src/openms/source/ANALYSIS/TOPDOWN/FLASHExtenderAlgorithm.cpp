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


  int FLASHExtenderAlgorithm::getVertex_(int node_index, int pro_index, int score, int num_mod) const
  {
    return ((node_index * (pro_length_ + 1) + pro_index) * (max_mod_cntr_ + 1) + num_mod) * (max_path_score_ - min_path_score_ + 1)
           + (score - min_path_score_);
  }

  int FLASHExtenderAlgorithm::getNodeIndex_(int vertex) const
  {
    return (vertex / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1)) / (pro_length_ + 1);
  }

  int FLASHExtenderAlgorithm::getProIndex_(int vertex) const
  {
    return ((vertex / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1))) % (pro_length_ + 1);
  }

  int FLASHExtenderAlgorithm::getScore_(int vertex) const
  {
    return 0;
  }

  int FLASHExtenderAlgorithm::getModNumber_(int vertex) const
  {
    return 0;
  }

  //void FLASHExtenderAlgorithm::constructVertices_()

  void FLASHExtenderAlgorithm::constructDAG_(FLASHHelperClasses::DAG& dag, std::vector<int>& sinks)
  {
    auto visited = boost::dynamic_bitset<>(dag.size());
    sinks.resize(max_mod_cntr_, 0);
    std::vector<int> tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends;

    // TODO get node and protein indices for tag start end positions
    if (tag_node_starts.empty()) return;

    int src = getVertex_(0, 0, 0, 0);
    connectBetweenTags(dag, visited, sinks, src, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends);
  }

  void FLASHExtenderAlgorithm::connectBetweenTags(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited, std::vector<int>& sinks, int vertex, std::vector<int>& tag_node_starts, std::vector<int>& tag_pro_starts, std::vector<int>& tag_node_ends, std::vector<int>& tag_pro_ends)
  {
    int node_index = getNodeIndex_(vertex);
    int pro_index = getProIndex_(vertex);

    if (node_index == tag_node_ends.back() && pro_index == tag_pro_ends.back()) // end of the last tag
    {
      extendBetweenTags(dag, visited, sinks,
                        vertex, -1, -1, true);
      return;
    }

    std::vector<int> next_vertices(max_mod_cntr_, 0);

    auto node_start_itr = std::find(tag_node_starts.begin(), tag_node_starts.end(), node_index);
    auto pro_start_itr = std::find(tag_pro_starts.begin(), tag_pro_starts.end(), pro_index);
    auto node_end_itr = std::find(tag_node_ends.begin(), tag_node_ends.end(), node_index);
    auto pro_end_itr = std::find(tag_pro_ends.begin(), tag_pro_ends.end(), pro_index);

    if (node_start_itr != tag_node_starts.end() && pro_start_itr != tag_pro_starts.end()) // within tag
    {
      int node_end = tag_node_ends[std::distance(tag_node_starts.begin(), node_start_itr)];
      int pro_end = tag_node_ends[std::distance(tag_pro_starts.begin(), pro_start_itr)];
      extendBetweenTags(dag, visited, next_vertices,
                         vertex, node_end, pro_end, true);
    }
    else if (node_end_itr != tag_node_ends.end() && pro_end_itr != tag_pro_ends.end()) // between tag. Never after the last tag.
    {
      int node_start = tag_node_starts[std::distance(tag_node_ends.begin(), node_end_itr) + 1];
      int pro_start = tag_node_ends[std::distance(tag_pro_ends.begin(), pro_end_itr) + 1];
      extendBetweenTags(dag, visited, next_vertices,
                        vertex, node_start, pro_start, false);
    }

    for (auto next_vertex : next_vertices)
    {
      if (next_vertex == 0) continue;
      connectBetweenTags(dag, visited, sinks, next_vertex, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends);
    }
  }

  void FLASHExtenderAlgorithm::extendBetweenTags(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited, std::vector<int>& sinks,
                                           int vertex, int node_index, int pro_index, bool within_tag)
  {
    if (!visited[vertex]) return;
    int node_index1 = getNodeIndex_(vertex);
    int score1 = getScore_(vertex);
    int num_mod1 = getModNumber_(vertex);
    int pro_index1 = getProIndex_(vertex);

    for (int score = score1 + 1; score <= max_path_score_; score++)
    {
      // if the starting point has already taken by a higher scoring path, don't go further.
      int higher_score_vertex = getVertex_(node_index1, pro_index1, score, num_mod1);
      if (visited[higher_score_vertex]) return;
    }

    if (node_index1 == node_index && pro_index1 == pro_index)
    {
      sinks[num_mod1] = vertex;
      return;
    }

    if (vertex == 0)
    {
      for (int pro_i = pro_index1 + 1; pro_i <= pro_index; pro_i++)
      {
        int vertex2 = getVertex_(node_index1, pro_i, score1, num_mod1);
        visited[vertex2] = true;
        extendBetweenTags(dag, visited, sinks, vertex2, node_index, pro_index);
      }
      return;
    }

    double delta_mass1 = node_masses_[node_index1] - pro_masses_[pro_index1];
    double delta_mass2 = node_masses_[node_index] - pro_masses_[pro_index];
    if (std::abs(delta_mass2 - delta_mass1) > max_mod_mass_ * max_mod_cntr_ + 1.0) return;

    for (int node_i = node_index1 + 1; node_i <= node_index; node_i++)
    {
      int score = score1 + node_scores_[node_i];
      for (int pro_i = pro_index1 + 1; pro_i <= pro_index; pro_i++)
      {
        double delta_mass = node_masses_[node_i] - pro_masses_[pro_i];
        if (delta_mass - delta_mass1 > max_mod_mass_) break;
        if (delta_mass1 - delta_mass > max_mod_mass_) continue;
        int num_mod = num_mod1;
        if (std::abs(delta_mass - delta_mass1) > tol_ * node_masses_[node_index1]) num_mod ++;
        if (within_tag && num_mod != num_mod1) continue; // on tag layer. No additional mod allowed.
        if (num_mod > max_mod_cntr_) continue;
        int next_vertex = getVertex_(node_i, pro_i, score, num_mod);
        visited[next_vertex] = true;
        extendBetweenTags(dag, visited, sinks, next_vertex, node_index, pro_index);
      }
    }
  }
} // namespace OpenMS