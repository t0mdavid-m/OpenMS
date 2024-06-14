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

    //registerStringOption_("fragment_type", "<choice>", "full", "For what type of sequence/fragment the mass should be computed\n", false);
    //setValidStrings_("fragment_type", ListUtils::create<String>("full,internal,N-terminal,C-terminal,a-ion,b-ion,c-ion,x-ion,y-ion,z-ion"));

    defaultsToParam_();
  }

  void FLASHExtenderAlgorithm::updateMembers_()
  {
    for (Size i = 0; i < Residue::SizeOfResidueType; i++)
    {
      auto res_type = Residue::ResidueType(i);
      res_type_names_[Residue::getResidueTypeName(res_type)] = res_type;
    }
  }

  int FLASHExtenderAlgorithm::getVertex_(int node_index, int pro_index, int score, int num_mod) const
  {
    return ((node_index * (pro_length_ + 1) + pro_index) * (max_mod_cntr_ + 1) + num_mod) * (max_path_score_ - min_path_score_ + 1)
           + (std::min(max_path_score_, score) - min_path_score_);
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
    return vertex % (max_path_score_ - min_path_score_ + 1) + min_path_score_;
  }

  int FLASHExtenderAlgorithm::getModNumber_(int vertex) const
  {
    return (vertex / (max_path_score_ - min_path_score_ + 1)) % (max_mod_cntr_ + 1);
  }

  void FLASHExtenderAlgorithm::run(const DeconvolvedSpectrum& dspec, const FLASHTaggerAlgorithm& tagger)
  {

    // TODO
    ion_types_str_ = std::vector<String>({"b-ion"});
    tol_ = 5e-6;
    max_mod_cntr_ = 1;

    std::vector<ProteinHit> hits;
    std::vector<double> mzs;
    std::vector<int> scores;
    tagger.getProteinHits(hits);
    auto spec = tagger.getSpectrum();

    std::vector<double> prefix_shifts;
    std::vector<double> suffix_shifts;

    Residue empty;
    for (const auto& ion_str : ion_types_str_)
    {
      auto res_type_ = res_type_names_[ion_str];
      auto letter = Residue::residueTypeToIonLetter(res_type_);
      double shift = empty.getMonoWeight(res_type_) + Residue::getInternalToFull().getMonoWeight();

      if (letter == "a" || letter == "b" || letter == "c")
      {
        prefix_shifts.push_back(shift);
      }
      else
      {
        suffix_shifts.push_back(shift);
      }
    }

    // 0 for suffix 1 for prefix 2 for suffix and prefix if precursor mass is available
    for (int k = 0; k < 2; k++) // TODO make it 3
    {
      MSSpectrum tspec;
      node_scores_.clear();
      node_masses_.clear();
      for (const auto& p : spec)
      {
        if (k == 0)
        {
          for (auto shift : suffix_shifts)
          {
            tspec.emplace_back(p.getMZ() - shift, p.getIntensity());
          }
        }
        else if (k == 1)
        {
          for (auto shift : prefix_shifts)
          {
            tspec.emplace_back(p.getMZ() - shift, p.getIntensity());
          }
        }
        else if (precursor_mass_ > 0)
        {
          for (auto shift : prefix_shifts)
          {
            tspec.emplace_back(p.getMZ() - shift, p.getIntensity());
          }
          for (auto shift : suffix_shifts)
          {
            tspec.emplace_back(precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - (p.getMZ() - shift), p.getIntensity()); // TODO check after all
          }
        }
      }
      if (tspec.empty()) continue;
      tspec.sortByPosition();
      // score, mz makes the nodes. merge with tolerance.
      node_scores_.push_back(0);
      node_masses_.push_back(0);

      for (int i = 0; i < tspec.size();i++)
      {
        double mass = tspec[i].getMZ();
        int score = (int)tspec[i].getIntensity();
        double prev_mass = node_masses_.back();
        int prev_score = node_scores_.back();
        if (mass - prev_mass < tol_ * mass) // they are the same
        {
          if (prev_score < score)
          {
            node_masses_.pop_back();
            node_scores_.pop_back();
          }
          score += prev_score;
        }
        node_masses_.push_back(mass);
        node_scores_.push_back(score);
      }

      max_path_score_ = 100; // TODO what if it hits the ceiling?
      min_path_score_ = 0;

     // #pragma omp parallel for default(none) shared(k, tagger, hits, dspec, prefix_shifts, suffix_shifts) TODO
      for (int i = 0; i < (int)hits.size(); i++)
      {
        auto hit = hits[i];
        // protein masses
        pro_masses_.clear();
        pro_masses_.push_back(0);
        auto seq = hit.getSequence();

        if (k == 0) seq = seq.reverse();
        pro_length_ = seq.length();
        for (const auto& aa : seq)
        {
          pro_masses_.push_back(pro_masses_.back() + AASequence::fromString(aa).getMonoWeight(Residue::Internal));
        }

        std::vector<int> tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends;
        std::vector<FLASHHelperClasses::Tag> tags;
        std::vector<int> sinks;

        tagger.getTagsMatchingTo(hit, tags);

        FLASHHelperClasses::DAG dag((int)(node_scores_.size() * (1 + pro_masses_.size()) * (1 + max_mod_cntr_) * (1 + max_path_score_ - min_path_score_)));

        for (const auto& tag : tags)
        {
          if (k == 0 && tag.getCtermMass() <= 0) continue;
          if (k == 1 && tag.getNtermMass() <= 0) continue;

          std::vector<int> positions;
          std::vector<double> masses;
          tagger.getMatchedPositionsAndFlankingMassDiffs(positions, masses, hit, tag);
          auto tag_masses = tag.getMzs();
          std::sort(tag_masses.begin(), tag_masses.end());

          for (int j = 0; j < positions.size(); j++)
          {
            int pos = positions[j]; // AASequence::fromString(aa).getMonoWeight(Residue::Internal));
            std::vector<double> start_masses, end_masses;
            if (tag.getCtermMass() >= 0) // suffix
            {
              if (precursor_mass_ == 0)
              {
                pos = seq.length() - pos; // invert pos
              }
              for (auto shift : suffix_shifts)
              {
                double start_mass = tag_masses[0] - shift;
                double end_mass = tag_masses.back() - shift;
                if (precursor_mass_ > 0)
                {
                  start_mass = precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - start_mass;
                  end_mass = precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - end_mass;
                }
                start_masses.push_back(start_mass);
                end_masses.push_back(end_mass);
              }
            }
            else if (tag.getNtermMass() >= 0) // prefix
            {
              for (auto shift : prefix_shifts)
              {
                double start_mass = tag_masses[0] - shift;
                double end_mass = tag_masses.back() - shift;
                start_masses.push_back(start_mass);
                end_masses.push_back(end_mass);
              }
            }

            for (int l = 0; l < start_masses.size(); l++)
            {
              double delta_start = start_masses[i] * tol_;
              auto lower_start = std::lower_bound(node_masses_.begin(), node_masses_.end(), start_masses[i] - delta_start);

              double delta_end = end_masses[i] * tol_;
              auto lower_end = std::lower_bound(node_masses_.begin(), node_masses_.end(), end_masses[i] - delta_end);

              int highest_score = -100;
              int highest_score_start = -1;
              int highest_score_end = -1;

              while(lower_start != node_masses_.end())
              {
                if (std::abs(*lower_start - start_masses[i]) < delta_start)
                {
                  int index = std::distance(node_masses_.begin(), lower_start);
                  if (highest_score < node_scores_[index])
                  {
                    highest_score_start = index;
                    highest_score = node_scores_[index];
                  }
                }
                lower_start++;
              }

              highest_score = -100;

              while(lower_end != node_masses_.end())
              {
                if (std::abs(*lower_end - end_masses[i]) < delta_end)
                {
                  int index = std::distance(node_masses_.begin(), lower_end);
                  if (highest_score < node_scores_[index])
                  {
                    highest_score_end = index;
                    highest_score = node_scores_[index];
                  }
                }
                lower_end++;
              }

              if (highest_score_start< 0 || highest_score_end < 0) continue;

              tag_node_starts.push_back(highest_score_start);
              tag_pro_starts.push_back(pos);
              tag_node_ends.push_back(highest_score_end);
              tag_pro_ends.push_back(pos + tag.getLength());
            }
          }
        }

        constructDAG_(dag, sinks, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends);
        std::vector<std::vector<int>> all_paths;
        for (int sink : sinks)
        {
          dag.findAllPaths(0, sink, all_paths, 1);
          std::cout<<"sink " << sink<<std::endl;
        }
        for (auto& path : all_paths)
        {
          for (int v : path)
          {
          //  std::cout<< node_masses_[getNodeIndex_(v)] << " " << pro_masses_[getProIndex_(v)]<<std::endl;
          }
          //std::cout<<std::endl;
        }
      // precursor_mass_

      }// add positive proteoforms all?
    }
  }

  void FLASHExtenderAlgorithm::constructDAG_(FLASHHelperClasses::DAG& dag, std::vector<int>& sinks, std::vector<int>& tag_node_starts, std::vector<int>& tag_pro_starts, std::vector<int>& tag_node_ends, std::vector<int>& tag_pro_ends)
  {
    auto visited = boost::dynamic_bitset<>(dag.size());
    sinks.resize(max_mod_cntr_, 0);

    if (tag_node_starts.empty()) return;

    int src = getVertex_(0, 0, 0, 0);
    visited[src] = true;
    connectBetweenTags(dag, visited, sinks, src, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends);
  }

  void FLASHExtenderAlgorithm::connectBetweenTags(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited, std::vector<int>& sinks, int vertex, std::vector<int>& tag_node_starts, std::vector<int>& tag_pro_starts, std::vector<int>& tag_node_ends, std::vector<int>& tag_pro_ends)
  {
    int node_index = getNodeIndex_(vertex);
    int pro_index = getProIndex_(vertex);

    if (pro_index == tag_pro_ends.back()) // end of the last tag
    {
      extendBetweenTags(dag, visited, sinks,
                        vertex, -1, -1, true);
      return;
    }

    auto node_start_itr = std::find(tag_node_starts.begin(), tag_node_starts.end(), node_index);
    auto pro_start_itr = std::find(tag_pro_starts.begin(), tag_pro_starts.end(), pro_index);
    auto node_end_itr = std::find(tag_node_ends.begin(), tag_node_ends.end(), node_index);
    auto pro_end_itr = std::find(tag_pro_ends.begin(), tag_pro_ends.end(), pro_index);

    if (node_start_itr != tag_node_starts.end() && pro_start_itr != tag_pro_starts.end()) // within tag
    {
      std::vector<int> next_vertices(max_mod_cntr_, 0);
      Size tag_index = std::distance(tag_node_starts.begin(), node_start_itr);
      int node_end = tag_node_ends[tag_index];
      int pro_end = tag_pro_ends[tag_index];
      extendBetweenTags(dag, visited, next_vertices,
                         vertex, node_end, pro_end, true);
      for (auto next_vertex : next_vertices)
      {
        if (next_vertex == 0) continue;
        connectBetweenTags(dag, visited, sinks, next_vertex, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends);
      }
    }
    else if (vertex == 0 || (node_end_itr != tag_node_ends.end() && pro_end_itr != tag_pro_ends.end())) // between tag. Never after the last tag.
    {
      for (int tag_index = 0; tag_index < tag_node_starts.size(); tag_index++) // for all reachable tag starting point, run extension
      {
        std::vector<int> next_vertices(max_mod_cntr_, 0);
        int node_start = tag_node_starts[tag_index];
        int pro_start = tag_pro_ends[tag_index];
        if (node_start <= node_index || pro_start <= pro_index) continue;
        extendBetweenTags(dag, visited, next_vertices, vertex, node_start, pro_start);
        for (auto next_vertex : next_vertices)
        {
          if (next_vertex == 0) continue;
          connectBetweenTags(dag, visited, sinks, next_vertex, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends);
        }
      }
    }
  }

  void FLASHExtenderAlgorithm::extendBetweenTags(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited, std::vector<int>& sinks,
                                           int vertex, int node_index, int pro_index, bool go_diagonal)
  {
    if (!visited[vertex]) return;
    int node_index1 = getNodeIndex_(vertex);
    int pro_index1 = getProIndex_(vertex);
    int score1 = getScore_(vertex);
    int num_mod1 = getModNumber_(vertex);

    if (node_index < 0)
    {
      sinks[num_mod1] = vertex;
      pro_index = std::min(pro_index1 + 5, (int)pro_masses_.size() - 1); // if sink is not specified, stretch up to five amino acids.
    }
    else if (node_index1 > node_index || pro_index1 > pro_index) return;

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
    double delta_mass;

    if (node_index >= 0)
    {
      delta_mass = node_masses_[node_index] - pro_masses_[pro_index];
      if (std::abs(delta_mass - delta_mass1) > max_mod_mass_ * max_mod_cntr_ + 1.0) return;
    }
    for (int node_i = node_index1 + 1; node_i <= (node_index < 0 ? (int)node_scores_.size() - 1: node_index); node_i++)
    {
      int score = score1 + node_scores_[node_i];
      for (int pro_i = pro_index1 + 1; pro_i <= pro_index; pro_i++)
      {
        delta_mass = node_masses_[node_i] - pro_masses_[pro_i];

        if (delta_mass - delta_mass1 > max_mod_mass_) continue;
        if (delta_mass1 - delta_mass > max_mod_mass_) continue;
        int num_mod = num_mod1;
        if (std::abs(delta_mass - delta_mass1) > tol_ * node_masses_[node_index1])
        {
          num_mod ++;
        }//else std::cout<<delta_mass<< " " << delta_mass1 << std::endl; this happens. do later.

        if (go_diagonal && num_mod != num_mod1) continue; // on tag layer. No additional mod allowed.
        if (num_mod > max_mod_cntr_) continue;
        int next_vertex = getVertex_(node_i, pro_i, score, num_mod);
        visited[next_vertex] = true;
        extendBetweenTags(dag, visited, sinks, next_vertex, node_index, pro_index);
      }
    }
  }
} // namespace OpenMS