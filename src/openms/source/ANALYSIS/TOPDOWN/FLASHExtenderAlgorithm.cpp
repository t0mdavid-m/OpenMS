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
  defaults_.setValue("max_mass_shift", 500.0, "Maximum mass shift for blind search.");
  defaults_.setValue("num_blind", 1, "Number of blind modification per protein.");

  defaults_.setValue("fdr", 1.0, "Protein FDR threshold.");
  defaults_.setMaxFloat("fdr", 1.0);
  defaults_.setMinFloat("fdr", 0.01);

  defaults_.setValue("keep_decoy", "false", "Keep decoy proteins.");
  defaults_.addTag("keep_decoy", "advanced");
  defaults_.setValidStrings("keep_decoy", {"true", "false"});

  // registerStringOption_("fragment_type", "<choice>", "full", "For what type of sequence/fragment the mass should be computed\n", false);
  // setValidStrings_("fragment_type", ListUtils::create<String>("full,internal,N-terminal,C-terminal,a-ion,b-ion,c-ion,x-ion,y-ion,z-ion"));

  defaultsToParam_();
}

void FLASHExtenderAlgorithm::updateMembers_()
{
  max_mod_cntr_ = param_.getValue("num_blind");
}

int FLASHExtenderAlgorithm::getVertex_(int node_index, int pro_index, int score, int num_mod) const
{
  return ((node_index * (pro_length_ + 1) + pro_index) * (max_mod_cntr_ + 1) + num_mod) * (max_path_score_ - min_path_score_ + 1)
         + (std::min(max_path_score_, std::max(min_path_score_, score)) - min_path_score_);
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
  setLogType(CMD);

  ion_types_str_ = std::vector<String>({"b-ion", "c-ion", "y-ion", "z-ion"}); // why ion type order matter? fix TODO
  tol_ = 5e-6;

  std::vector<ProteinHit> hits;
  std::vector<double> mzs;
  std::vector<int> scores;
  tagger.getProteinHits(hits);
  startProgress(0, (SignedSize)hits.size() * 3, "running FLASHExtender ...");

  auto spec = tagger.getSpectrum();

  std::vector<double> prefix_shifts;
  std::vector<double> suffix_shifts;

  Residue empty;
  for (const auto& ion_str : ion_types_str_)
  {
    if (ion_str == "a-ion") { prefix_shifts.push_back(Residue::getInternalToAIon().getMonoWeight()); }
    else if (ion_str == "b-ion") { prefix_shifts.push_back(Residue::getInternalToBIon().getMonoWeight()); }
    else if (ion_str == "c-ion") { prefix_shifts.push_back(Residue::getInternalToCIon().getMonoWeight()); }
    else if (ion_str == "x-ion") { suffix_shifts.push_back(Residue::getInternalToXIon().getMonoWeight()); }
    else if (ion_str == "y-ion") { suffix_shifts.push_back(Residue::getInternalToYIon().getMonoWeight()); }
    else if (ion_str == "z-ion") { suffix_shifts.push_back(Residue::getInternalToZIon().getMonoWeight()); }
    else
    {
      continue; // TODO warn
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
        for (const auto& shift : suffix_shifts)
        {
          tspec.emplace_back(p.getMZ() - shift, p.getIntensity());
        }
      }
      else if (k == 1)
      {
        for (const auto& shift : prefix_shifts)
        {
          tspec.emplace_back(p.getMZ() - shift, p.getIntensity());
        }
      }
      else if (precursor_mass_ > 0)
      {
        for (const auto& shift : prefix_shifts)
        {
          tspec.emplace_back(p.getMZ() - shift, p.getIntensity());
        }
        for (const auto& shift : suffix_shifts)
        {
          tspec.emplace_back(precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - (p.getMZ() - shift),
                             p.getIntensity()); // TODO check after all
        }
      }
    }
    if (tspec.empty()) continue;
    tspec.sortByPosition();

    // score, mz makes the nodes. merge with tolerance.
    node_scores_.push_back(0);
    node_masses_.push_back(0);

    max_path_score_ = 500;
    min_path_score_ = -50;

    for (int i = 0; i < tspec.size(); i++)
    {
      double mass = tspec[i].getMZ();
      int score = (int)tspec[i].getIntensity();

      // if (score > 0) max_path_score_ += score;
      // else min_path_score_ += score;

      double prev_mass = node_masses_.back();
      int prev_score = node_scores_.back();
      if (mass - prev_mass < 2 * tol_ * mass) // they are the same
      {
        if (score > 0)
        {
          node_masses_.pop_back();
          node_scores_.pop_back();
        }
        score += prev_score;
      }
      node_masses_.push_back(mass);
      node_scores_.push_back(score);
    }

//#pragma omp parallel for default(none) shared(k, tagger, hits, prefix_shifts, suffix_shifts)
    for (int i = 0; i < (int)hits.size(); i++)
    {
      nextProgress();
      auto hit = hits[i];
      // protein masses
      pro_masses_.clear();
      pro_masses_.push_back(0);
      auto seq = hit.getSequence();

      if (k == 0) seq = seq.reverse();
      for (const auto& aa : seq)
      {
        pro_masses_.push_back(pro_masses_.back() + AASequence::fromString(aa).getMonoWeight(Residue::Internal));
      }
      pro_length_ = pro_masses_.size();

      std::vector<int> tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends;
      std::vector<FLASHHelperClasses::Tag> tags;
      std::set<int> sinks;

      tagger.getTagsMatchingTo(hit, tags);

      // std::cout<<hit.getAccession()<<std::endl;
      // for (const auto& tag : tags) std::cout<<tag.toString()<<std::endl;
      FLASHHelperClasses::DAG dag(
        (int)(node_scores_.size() * (1 + pro_masses_.size()) * (1 + max_mod_cntr_) * (1 + max_path_score_ - min_path_score_)));

      for (const auto& tag : tags)
      {
        if (k == 0 && tag.getCtermMass() <= 0) continue;
        if (k == 1 && tag.getNtermMass() <= 0) continue;
        if (k == 2 && precursor_mass_ <= 0) continue;
        std::vector<int> positions;
        std::vector<double> masses;

        tagger.getMatchedPositionsAndFlankingMassDiffs(positions, masses, hit, tag);
        auto tag_masses = tag.getMzs();
        std::sort(tag_masses.begin(), tag_masses.end());

        for (int j = 0; j < positions.size(); j++)
        {
          int pos = positions[j];
          std::vector<double> start_masses, end_masses;
          if (tag.getCtermMass() >= 0) // suffix
          {
            if (precursor_mass_ == 0)
            {
              pos = pro_length_ - 1 - pos; // invert pos
            }
            for (const auto& shift : suffix_shifts)
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
            double delta_start = start_masses[l] * tol_ * 2;
            auto lower_start = std::lower_bound(node_masses_.begin(), node_masses_.end(), start_masses[l] - delta_start);

            double delta_end = end_masses[l] * tol_ * 2;
            auto lower_end = std::lower_bound(node_masses_.begin(), node_masses_.end(), end_masses[l] - delta_end);

            int highest_score = min_path_score_;
            int highest_score_start = -1;
            int highest_score_end = -1;

            while (lower_start != node_masses_.end())
            {
              if (std::abs(*lower_start - start_masses[l]) < delta_start)
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

            highest_score = -1e5; // TODO

            while (lower_end != node_masses_.end())
            {
              if (std::abs(*lower_end - end_masses[l]) < delta_end)
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

            if (highest_score_start < 0 || highest_score_end < 0) continue;

            tag_node_starts.push_back(highest_score_start);
            tag_node_ends.push_back(highest_score_end);
            if (k == 0) // suffix inverted
            {
              tag_pro_starts.push_back(pos - tag.getLength());
              tag_pro_ends.push_back(pos); // this can be much faster...
            }
            else
            {
              tag_pro_starts.push_back(pos);
              tag_pro_ends.push_back(pos + tag.getLength()); // this can be much faster...
            }
            // std::cout<< " * " << tag_node_starts.back() << " " << node_masses_[tag_node_starts.back()] << " " << node_masses_[tag_node_ends.back()]
            // << " " <<
            //    pro_masses_[tag_pro_starts.back()] << " " << pro_masses_[tag_pro_ends.back()] << std::endl;
          }
        }
      }

      constructDAG_(dag, sinks, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends);

      std::vector<std::vector<int>> all_paths;
      int src = getVertex_(0, 0, 0, 0);
      std::set<int> best_sinks;
      int best_score = min_path_score_;
      for (int sink : sinks)
      {
        if (sink == src || getScore_(sink) < best_score) continue;
        best_score = getScore_(sink);
      }
      for (int sink : sinks)
      {
        if(sink == src || getScore_(sink) < best_score) continue;
        best_sinks.insert(sink);
      }
//#pragma omp critical
      for (const auto& best_sink : best_sinks)
      {
        std::cout<<best_sink<<std::endl;
        all_paths.clear();
        dag.findAllPaths(best_sink, src, all_paths, 3);
        if (all_paths.empty()) continue;
        std::cout << hit.getAccession() << " " << dag.size() << std::endl;
        for (auto& path : all_paths)
        {
          std::cout << path.size() << std::endl;
          for (int v : path)
          {
            std::cout << v << " " << getNodeIndex_(v) << " " << node_masses_[getNodeIndex_(v)] << " " << pro_masses_[getProIndex_(v)] << " " << getModNumber_(v)
                      << " " << getScore_(v) << std::endl;
          }
          std::cout << std::endl;
        }
      }

      // precursor_mass_
    } // add positive proteoforms all?
  }
  endProgress();
}

void FLASHExtenderAlgorithm::constructDAG_(FLASHHelperClasses::DAG& dag,
                                           std::set<int>& sinks,
                                           std::vector<int>& tag_node_starts,
                                           std::vector<int>& tag_pro_starts,
                                           std::vector<int>& tag_node_ends,
                                           std::vector<int>& tag_pro_ends)
{
  auto visited = boost::dynamic_bitset<>(dag.size());

  sinks.clear();

  if (tag_node_starts.empty()) return;

  int src = getVertex_(0, 0, 0, 0);

  visited[src] = true;
  connectBetweenTags(dag, visited, sinks, src, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends);
}

void FLASHExtenderAlgorithm::connectBetweenTags(FLASHHelperClasses::DAG& dag,
                                                boost::dynamic_bitset<>& visited,
                                                std::set<int>& sinks,
                                                int vertex,
                                                std::vector<int>& tag_node_starts,
                                                std::vector<int>& tag_pro_starts,
                                                std::vector<int>& tag_node_ends,
                                                std::vector<int>& tag_pro_ends)
{
  int node_index = getNodeIndex_(vertex);
  int pro_index = getProIndex_(vertex);
  //int num_mod = getModNumber_(vertex);

  if (std::find(tag_pro_ends.begin(), tag_pro_ends.end(), pro_index) != tag_pro_ends.end()) // end of the last tag
  {
    extendBetweenTags(dag, visited, sinks, vertex, -1, -1, 1e5); //
    return;
  }

  int tag_start_index = -1;
  int tag_end_index = -1;

  for (int i = 0; i < tag_node_starts.size(); i++)
  {
    if (tag_node_starts[i] == node_index && tag_pro_starts[i] == pro_index) tag_start_index = i;
    if (tag_node_ends[i] == node_index && tag_pro_ends[i] == pro_index) tag_end_index = i;
  }

  int src = getVertex_(0, 0, 0, 0);
  if (tag_start_index >= 0) // within tag
  {
    std::set<int> next_vertices;
    int node_end = tag_node_ends[tag_start_index];
    int pro_end = tag_pro_ends[tag_start_index];

    extendBetweenTags(dag, visited, next_vertices, vertex, node_end, pro_end, 1e5);

    for (auto next_vertex : next_vertices)
    {
      if (next_vertex == src) continue;
      connectBetweenTags(dag, visited, sinks, next_vertex, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends);
    }
  }
  else if (vertex == src || tag_end_index >= 0) // between tag. Never after the last tag.
  {
    for (int tag_index = 0; tag_index < tag_node_starts.size(); tag_index++) // for all reachable tag starting point, run extension
    {
      std::set<int> next_vertices;
      int node_start = tag_node_starts[tag_index];
      int pro_start = tag_pro_starts[tag_index];
      if (node_start <= node_index || pro_start <= pro_index) continue;
      extendBetweenTags(dag, visited, next_vertices, vertex, node_start, pro_start, 0);

      for (auto next_vertex : next_vertices)
      {
        if (next_vertex == src) continue;
        connectBetweenTags(dag, visited, sinks, next_vertex, tag_node_starts, tag_pro_starts, tag_node_ends, tag_pro_ends);
      }
    }
  }
}

void FLASHExtenderAlgorithm::extendBetweenTags(FLASHHelperClasses::DAG& dag,
                                               boost::dynamic_bitset<>& visited,
                                               std::set<int>& sinks,
                                               int vertex,
                                               int node_index,
                                               int pro_index,
                                               int diagonal_counter)
{
  if (! visited[vertex]) return;
  int node_index1 = getNodeIndex_(vertex);
  int pro_index1 = getProIndex_(vertex);
  int score1 = getScore_(vertex);
  int num_mod1 = getModNumber_(vertex);
  if (num_mod1 == max_mod_cntr_) diagonal_counter = 1e5;

  for (int score = score1 + 1; score <= max_path_score_; score++)
  {
    // if the starting point has already taken by a higher scoring path, don't go further.
    int higher_score_vertex = getVertex_(node_index1, pro_index1, score, num_mod1);
    if (visited[higher_score_vertex]) return;
  }

  if (node_index < 0 && diagonal_counter > 0)
  {
    sinks.insert(vertex);
    pro_index = std::min(pro_index1 + 10, (int)pro_masses_.size() - 1); // if sink is not specified, stretch up to ten amino acids.
  }
  else if (node_index1 > node_index || pro_index1 > pro_index)
    return;

  if (node_index1 == node_index && pro_index1 == pro_index)
  {
    sinks.insert(vertex);
    return;
  }

  if (vertex == getVertex_(0, 0, 0, 0))
  {
    for (int pro_i = pro_index1 + 1; pro_i <= pro_index; pro_i++)
    {
      int vertex2 = getVertex_(node_index1, pro_i, score1, num_mod1);

      if (visited[vertex2] && dag.hasEdge(vertex2, vertex)) continue;
      dag.addEdge(vertex2, vertex, visited);
      if (! visited[vertex2]) continue;
      extendBetweenTags(dag, visited, sinks, vertex2, node_index, pro_index, diagonal_counter);
    }
  }

  double delta_mass1 = node_masses_[node_index1] - pro_masses_[pro_index1];
  double delta_mass;

  if (node_index >= 0)
  {
    delta_mass = node_masses_[node_index] - pro_masses_[pro_index];
    if (std::abs(delta_mass - delta_mass1) > max_mod_mass_ * (max_mod_cntr_ - num_mod1) + 1.0) return;
    if (diagonal_counter > 0)
    {
      if (std::abs(delta_mass - delta_mass1) > 2 * tol_ * node_masses_[node_index]) return;
    }
  }

  for (int node_i = node_index1 + 1; node_i <= (node_index < 0 ? (int)node_scores_.size() - 1 : node_index); node_i++)
  {
    if (node_index < 0 && node_scores_[node_i] < 0) continue;
    int score = score1 + node_scores_[node_i];
    for (int pro_i = pro_index1 + 1; pro_i <= pro_index; pro_i++)
    {
      double t_delta_mass = node_masses_[node_i] - pro_masses_[pro_i];

      if (t_delta_mass - delta_mass1 > max_mod_mass_) continue;
      if (delta_mass1 - t_delta_mass > max_mod_mass_) break;

      int num_mod = num_mod1;
      if (std::abs(t_delta_mass - delta_mass1) > 2 * tol_ * node_masses_[node_i]) { num_mod++; }

      if (diagonal_counter > 0 && num_mod != num_mod1) continue; //
      if (num_mod > max_mod_cntr_) continue;
      if (node_index >= 0 && num_mod == max_mod_cntr_ && std::abs(t_delta_mass - delta_mass) > 2 * tol_ * node_masses_[node_i]) continue;

      int next_vertex = getVertex_(node_i, pro_i, score, num_mod);

      if (visited[next_vertex] && dag.hasEdge(next_vertex, vertex)) continue;
      dag.addEdge(next_vertex, vertex, visited);
      if (! visited[next_vertex]) continue;

      int next_diagonal_counter = diagonal_counter;
      if (diagonal_counter > 0) next_diagonal_counter --;
      else if (num_mod != num_mod1) next_diagonal_counter = 1;
      extendBetweenTags(dag, visited, sinks, next_vertex, node_index, pro_index, next_diagonal_counter);
    }
  }
}
} // namespace OpenMS