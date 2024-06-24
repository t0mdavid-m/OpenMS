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
  defaults_.setValue("max_mass_shift", 10.0, "Maximum mass shift for blind search.");
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
  max_mod_mass_ = param_.getValue("max_mass_shift");
}

Size FLASHExtenderAlgorithm::getVertex_(int node_index, int pro_index, int score, int num_mod, Size pro_length) const
{
  return ((node_index * (pro_length + 1) + pro_index) * (max_mod_cntr_ + 1) + num_mod) * (max_path_score_ - min_path_score_ + 1)
         + (std::min(max_path_score_, std::max(min_path_score_, score)) - min_path_score_);
}

int FLASHExtenderAlgorithm::getNodeIndex_(Size vertex, Size pro_length) const
{
  return (vertex / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1)) / (pro_length + 1);
}

int FLASHExtenderAlgorithm::getProIndex_(Size vertex, Size pro_length) const
{
  return ((vertex / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1))) % (pro_length + 1);
}

int FLASHExtenderAlgorithm::getScore_(Size vertex) const
{
  return vertex % (max_path_score_ - min_path_score_ + 1) + min_path_score_;
}

int FLASHExtenderAlgorithm::getModNumber_(Size vertex) const
{
  return (vertex / (max_path_score_ - min_path_score_ + 1)) % (max_mod_cntr_ + 1);
}

void FLASHExtenderAlgorithm::calcualte_precursor_mass_(const std::vector<MSSpectrum>& node_spectrum_list,
                                                       const std::vector<std::vector<double>>& pro_masses_list,
                                                       const std::map<int, std::vector<Size>>& best_paths)
{
  if (node_spectrum_list.size() < 2 || pro_masses_list.size() < 2) return;
  auto npath = best_paths.find(0);
  if (npath == best_paths.end()) return;

  auto cpath = best_paths.find(1);
  if (cpath == best_paths.end()) return;

  Size npro_index = getProIndex_(npath->second[0], pro_masses_list[0].size());
  Size cpro_index = getProIndex_(cpath->second[0], pro_masses_list[0].size());

  Size nnode_index = getNodeIndex_(npath->second[0], pro_masses_list[0].size());
  Size cnode_index = getNodeIndex_(cpath->second[0], pro_masses_list[0].size());

  precursor_mass_
    = node_spectrum_list[0][nnode_index].getMZ() + node_spectrum_list[1][cnode_index].getMZ() + Residue::getInternalToFull().getMonoWeight();

  if (npro_index + cpro_index >= pro_masses_list[0].size() - 1)
  {
    double overlapping_mass = 0; // TODO consider PTM mass. + how to end the C  term truncation
    while (npro_index >= 0 && npro_index + cpro_index > pro_masses_list[0].size() - 1)
    {
      overlapping_mass += pro_masses_list[0][npro_index];
      npro_index--;
    }
    while (cpro_index >= 0 && npro_index + cpro_index > pro_masses_list[0].size() - 1)
    {
      overlapping_mass += pro_masses_list[1][cpro_index];
      cpro_index--;
    }

    precursor_mass_ -= overlapping_mass;
  }
}

void FLASHExtenderAlgorithm::get_pro_masses_(const ProteinHit& hit, std::vector<double>& pro_masses, int mode)
{
  pro_masses.reserve(hit.getSequence().size() + 1);
  pro_masses.push_back(0);

  auto seq = hit.getSequence();
  pro_masses.reserve(seq.size());

  if (mode == 0) seq = seq.reverse();
  for (const auto& aa : seq)
  {
    pro_masses.push_back(pro_masses.back() + AASequence::fromString(aa).getMonoWeight(Residue::Internal));
  }
}

void FLASHExtenderAlgorithm::define_nodes_(const FLASHTaggerAlgorithm& tagger, MSSpectrum& node_spec, MSSpectrum& tol_spec, double max_mass, int mode)
{
  const auto& spec = tagger.getSpectrum();
  // 0 for suffix 1 for prefix 2 for suffix and prefix if precursor mass is available
  MSSpectrum t_node_spec, t_tol_spec;

  t_node_spec.reserve(spec.size() * ion_types_str_.size() + 1);
  t_node_spec.emplace_back(0, 0);
  t_tol_spec.reserve(spec.size() * ion_types_str_.size() + 1);
  t_tol_spec.emplace_back(0, 0);

  for (const auto& p : spec)
  {
    if (mode == 0)
    {
      for (const auto& shift : suffix_shifts_)
      {
        double mass = p.getMZ() - shift;
        if (mass <= 0 || mass > max_mass + 1) continue;
        t_node_spec.emplace_back(mass, p.getIntensity());
        t_tol_spec.emplace_back(mass, 2 * tol_ * mass);
      }
    }
    else if (mode == 1)
    {
      for (const auto& shift : prefix_shifts_)
      {
        double mass = p.getMZ() - shift;
        if (mass <= 0 || mass > max_mass + 1) continue;
        t_node_spec.emplace_back(mass, p.getIntensity());
        t_tol_spec.emplace_back(mass, 2 * tol_ * mass);
      }
    }
    else if (precursor_mass_ > 0)
    {
      for (const auto& shift : prefix_shifts_)
      {
        double mass = p.getMZ() - shift;
        if (mass <= 0 || mass > precursor_mass_ - Residue::getInternalToFull().getMonoWeight() + 1) continue;
        t_node_spec.emplace_back(mass, p.getIntensity());
        t_tol_spec.emplace_back(mass, 2 * tol_ * mass);
      }
      for (const auto& shift : suffix_shifts_)
      {
        double mass = p.getMZ() - shift;
        if (mass <= 0 || mass > precursor_mass_ - Residue::getInternalToFull().getMonoWeight() + 1) continue;
        t_node_spec.emplace_back(precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - mass, p.getIntensity());
        t_tol_spec.emplace_back(precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - mass, 2 * tol_ * mass);
      }
    }
  }
  if (precursor_mass_ > 0)
  {
    t_node_spec.emplace_back(precursor_mass_ - Residue::getInternalToFull().getMonoWeight(), 1);
    t_tol_spec.emplace_back(precursor_mass_ - Residue::getInternalToFull().getMonoWeight(),
                            2 * tol_ * (precursor_mass_ - Residue::getInternalToFull().getMonoWeight()));
  }

  t_node_spec.sortByPosition();
  t_tol_spec.sortByPosition();
  // Assign the sorted values back to the original vectors
  node_spec.reserve(t_node_spec.size());
  tol_spec.reserve(t_node_spec.size());

  for (int k = 0; k < t_node_spec.size(); k++)
  {
    const auto& p = t_node_spec[k];
    double mass = p.getMZ();
    float score = p.getIntensity();

    // if (score > 0) max_path_score_ += score;
    // else min_path_score_ += score;

    if (mass - node_spec.back().getMZ() < t_tol_spec[k].getIntensity()) // they are the same
    {
      float prev_score = node_spec.back().getIntensity();
      if (node_spec.size() > 1)
      {
        node_spec.pop_back();
        tol_spec.pop_back();
      }
      score = std::max(.0f, score);
      score += std::max(.0f, prev_score);
    }
    node_spec.emplace_back(mass, score);
    tol_spec.emplace_back(mass, t_tol_spec[k].getIntensity());
  }
  node_spec.sortByPosition();
  tol_spec.sortByPosition();
}

void FLASHExtenderAlgorithm::run_(const FLASHTaggerAlgorithm& tagger,
                                  const ProteinHit& hit,
                                  const MSSpectrum& node_spec,
                                  const MSSpectrum& tol_spec,
                                  const std::vector<double>& pro_masses,
                                  std::vector<std::vector<Size>>& all_paths,
                                  int mode) // per hit
{
  std::vector<std::vector<int>> tag_edges(4); // pro start end indices and node start end indices
  std::vector<FLASHHelperClasses::Tag> tags;
  std::set<Size> sinks;

  tagger.getTagsMatchingTo(hit, tags);

  for (auto edge : tag_edges)
    edge.reserve(tags.size() * 2);

  FLASHHelperClasses::DAG dag((1 + node_spec.size()) * (1 + pro_masses.size()) * (1 + max_mod_cntr_) * (1 + max_path_score_ - min_path_score_));

  for (const auto& tag : tags)
  {
    if (mode == 0 && tag.getCtermMass() <= 0) continue;
    if (mode == 1 && tag.getNtermMass() <= 0) continue;
    if (mode == 2 && precursor_mass_ <= 0) continue;
    std::vector<int> positions;
    std::vector<double> masses;

    tagger.getMatchedPositionsAndFlankingMassDiffs(positions, masses, hit, tag);
    auto tag_masses = tag.getMzs();
    std::sort(tag_masses.begin(), tag_masses.end());

    std::vector<double> start_masses, end_masses, start_tols, end_tols;
    if (tag.getCtermMass() >= 0) // suffix
    {
      for (const auto& shift : suffix_shifts_)
      {
        double start_mass = tag_masses[0] - shift;
        double end_mass = tag_masses.back() - shift;

        start_tols.push_back(start_mass * tol_ * 2);
        end_tols.push_back(end_mass * tol_ * 2);

        if (precursor_mass_ > 0)
        {
          start_mass = precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - start_mass;
          end_mass = precursor_mass_ - Residue::getInternalToFull().getMonoWeight() - end_mass;
          start_masses.push_back(end_mass);
          end_masses.push_back(start_mass);
        }
        else
        {
          start_masses.push_back(start_mass);
          end_masses.push_back(end_mass);
        }
      }
    }
    else
    {
      for (const auto& shift : prefix_shifts_)
      {
        double start_mass = tag_masses[0] - shift;
        double end_mass = tag_masses.back() - shift;
        start_tols.push_back(start_mass * tol_ * 2);
        end_tols.push_back(end_mass * tol_ * 2);
        start_masses.push_back(start_mass);
        end_masses.push_back(end_mass);
      }
    }

    for (int l = 0; l < start_masses.size(); l++)
    {
      double delta_start = start_tols[l];
      double delta_end = end_tols[l];

      int highest_score_start = node_spec.findHighestInWindow(start_masses[l], delta_start, delta_start);
      int highest_score_end = node_spec.findHighestInWindow(end_masses[l], delta_end, delta_end);

      if (highest_score_start < 0 || highest_score_end < 0) continue;
      for (int pos : positions)
      {
        if (tag.getCtermMass() >= 0) // suffix
        {
          if (precursor_mass_ == 0)
          {
            pos = pro_masses.size() - 1 - pos; // invert pos
          }
        }

        if (mode == 0) // suffix inverted
        {
          if (pos - tag.getLength() >= 0 && pos < node_spec.size())
          {
            tag_edges[0].emplace_back(pos - tag.getLength());
            tag_edges[1].emplace_back(pos);
            tag_edges[2].emplace_back(highest_score_start);
            tag_edges[3].emplace_back(highest_score_end);
          }
        }
        else
        {
          if (pos >= 0 && pos + tag.getLength() < node_spec.size())
          {
            tag_edges[0].emplace_back(pos);
            tag_edges[1].emplace_back(pos + tag.getLength()); // this can be much faster...
            tag_edges[2].emplace_back(highest_score_start);
            tag_edges[3].emplace_back(highest_score_end);
          }
        }
      }
    }
  }
  constructDAG_(dag, sinks, node_spec, tol_spec, pro_masses, tag_edges, mode);

  Size src = getVertex_(0, 0, 0, 0, pro_masses.size());
  std::set<Size> best_sinks;
  int max_score = 0;
  int min_mod = max_mod_cntr_;
  Size max_len = 0;
  for (Size sink : sinks)
  {
    if (sink == src || getScore_(sink) < max_score) continue; // getModNumber_(sink) == 0 ||
    if (precursor_mass_ > 0 && getNodeIndex_(sink, pro_masses.size()) < node_spec.size() - 1) continue;
    max_score = getScore_(sink);
  }

  for (Size sink : sinks)
  {
    if (sink == src || getScore_(sink) < max_score) continue; // || getModNumber_(sink) > min_mod
    if (precursor_mass_ > 0 && getNodeIndex_(sink, pro_masses.size()) < node_spec.size() - 1) continue;
    min_mod = std::min(min_mod, getModNumber_(sink));
  }

  for (Size sink : sinks)
  {
    if (sink == src || getScore_(sink) < max_score || getModNumber_(sink) > min_mod) continue; //
    if (precursor_mass_ > 0 && getNodeIndex_(sink, pro_masses.size()) < node_spec.size() - 1) continue;
    best_sinks.insert(sink);
  }

  for (const auto& best_sink : best_sinks)
  {
    std::vector<std::vector<Size>> all_paths_per_sink;
    dag.findAllPaths(best_sink, src, all_paths_per_sink, 0);

    for (const auto& path : all_paths_per_sink)
    {
      max_len = std::max(max_len, path.size());
    }
  }

  for (const auto& best_sink : best_sinks)
  {
    std::vector<std::vector<Size>> all_paths_per_sink;
    dag.findAllPaths(best_sink, src, all_paths_per_sink, 0);

    for (const auto& path : all_paths_per_sink)
    {
      if (path.size() < max_len) continue;
      all_paths.push_back(path);
    }
  }

}

void FLASHExtenderAlgorithm::run(const FLASHTaggerAlgorithm& tagger)
{
  setLogType(CMD);

  ion_types_str_ = std::vector<String>({"b-ion", "c-ion", "y-ion", "z-ion"});
  tol_ = 5e-6;
  std::vector<ProteinHit> hits;
  std::vector<double> mzs;
  std::vector<int> scores;
  tagger.getProteinHits(hits);

  startProgress(0, (SignedSize)std::min(5, (int)hits.size()), "running FLASHExtender ...");

  Residue empty;
  for (const auto& ion_str : ion_types_str_)
  {
    if (ion_str == "a-ion") { prefix_shifts_.push_back(Residue::getInternalToAIon().getMonoWeight()); }
    else if (ion_str == "b-ion") { prefix_shifts_.push_back(Residue::getInternalToBIon().getMonoWeight()); }
    else if (ion_str == "c-ion") { prefix_shifts_.push_back(Residue::getInternalToCIon().getMonoWeight()); }
    else if (ion_str == "x-ion") { suffix_shifts_.push_back(Residue::getInternalToXIon().getMonoWeight()); }
    else if (ion_str == "y-ion") { suffix_shifts_.push_back(Residue::getInternalToYIon().getMonoWeight()); }
    else if (ion_str == "z-ion") { suffix_shifts_.push_back(Residue::getInternalToZIon().getMonoWeight()); }
    else
    {
      continue; // TODO warn
    }
  }

  max_path_score_ = 300;
  min_path_score_ = -100;

#pragma omp parallel for default(none) shared(hits, tagger, std::cout)
  for (int i = 0; i < std::min(10, (int)hits.size()); i++)
  {
    nextProgress();
    auto hit = hits[i];
    int total_score = 0;

    std::map<int, std::vector<Size>> best_paths;
    std::vector<MSSpectrum> node_spectrum_list;
    std::vector<std::vector<double>> pro_masses_list;
    for (int mode = 0; mode < 3; mode++)
    {
      if (mode == 2 && precursor_mass_ <= 0)
      {
        calcualte_precursor_mass_(node_spectrum_list, pro_masses_list, best_paths);
        if (precursor_mass_ <= 0) break;
      }

      std::vector<double> pro_masses;
      MSSpectrum node_spec, tol_spec;
      get_pro_masses_(hit, pro_masses, mode);
      define_nodes_(tagger, node_spec, tol_spec, pro_masses.back(), mode);

      if (mode < 2)
      {
        pro_masses_list.push_back(pro_masses);
        node_spectrum_list.push_back(node_spec);
      }

      std::vector<std::vector<Size>> all_paths;

      run_(tagger, hit, node_spec, tol_spec, pro_masses, all_paths, mode);
      // maybe filter out paths
      int score = 0;
      for (const auto& path : all_paths)
      {
        best_paths[mode] = path;
        score = getScore_(path.back());
        for (Size v : path)
        {
          std::cout << getNodeIndex_(v, pro_masses.size()) << " " << node_spec[getNodeIndex_(v, pro_masses.size())].getMZ() << " "
                    << pro_masses[getProIndex_(v, pro_masses.size())] << " " << getModNumber_(v) << " " << getScore_(v) << std::endl;
        }
        std::cout << hit.getAccession() << " " << path.size() << std::endl;
        break;
      }

      if (mode == 2) total_score = score;
      else
        total_score += score;
    }
    hit.setScore(total_score);
  } // add positive proteoforms all?

  endProgress();
}

void FLASHExtenderAlgorithm::constructDAG_(FLASHHelperClasses::DAG& dag,
                                           std::set<Size>& sinks,
                                           const MSSpectrum& node_spec,
                                           const MSSpectrum& tol_spec,
                                           const std::vector<double>& pro_masses,
                                           const std::vector<std::vector<int>>& tag_edges,
                                           int mode)
{
  auto visited = boost::dynamic_bitset<>(dag.size());
  if (tag_edges.empty()) return;
  Size src = getVertex_(0, 0, 0, 0, pro_masses.size());

  visited[src] = true;
  std::set<Size> visited_tag_edges;
  connectBetweenTags_(dag, visited, visited_tag_edges, sinks, src, node_spec, tol_spec, pro_masses, tag_edges, mode);
}

void FLASHExtenderAlgorithm::connectBetweenTags_(FLASHHelperClasses::DAG& dag,
                                                 boost::dynamic_bitset<>& visited,
                                                 std::set<Size>& visited_tag_edges,
                                                 std::set<Size>& sinks,
                                                 Size vertex,
                                                 const MSSpectrum& node_spec,
                                                 const MSSpectrum& tol_spec,
                                                 const std::vector<double>& pro_masses,
                                                 const std::vector<std::vector<int>>& tag_edges,
                                                 int mode)
{
  if (tag_edges.size() < 4) return;
  int node_index = getNodeIndex_(vertex, pro_masses.size());
  int pro_index = getProIndex_(vertex, pro_masses.size());

  int tag_start_index = -1;
  int tag_end_index = -1;

  const auto& tag_pro_starts = tag_edges[0];
  const auto& tag_pro_ends = tag_edges[1];
  const auto& tag_node_starts = tag_edges[2];
  const auto& tag_node_ends = tag_edges[3];

  for (int i = 0; i < tag_pro_starts.size(); i++)
  {
    if (tag_start_index < 0 && tag_node_starts[i] == node_index && tag_pro_starts[i] == pro_index) { tag_start_index = i; }
    if (tag_end_index < 0 && tag_node_ends[i] == node_index && tag_pro_ends[i] == pro_index) { tag_end_index = i; }
  }

  Size src = getVertex_(0, 0, 0, 0, pro_masses.size());

  if (tag_start_index >= 0) // within tag
  {
    int node_end = -1;
    int i = tag_start_index;
    int pro_end = -1;
    while ((i < tag_node_starts.size()) && (tag_node_starts[i] == node_index) && (tag_pro_starts[i] == pro_index))
    {
      if (pro_end < tag_pro_ends[i]) pro_end = tag_pro_ends[i];
      i++;
    }
    i = tag_start_index;
    while ((i < tag_node_starts.size()) && (tag_node_starts[i] == node_index) && (tag_pro_starts[i] == pro_index))
    {
      if (pro_end == tag_pro_ends[i]) { node_end = std::max(tag_node_ends[i], node_end); }
      i++;
    }

    std::set<Size> next_vertices;

    extendBetweenTags_(dag, visited, next_vertices, vertex, node_end, pro_end, 1e5, node_spec, tol_spec, pro_masses, mode);

    for (auto next_vertex : next_vertices)
    {
      if (next_vertex == src) continue;
      connectBetweenTags_(dag, visited, visited_tag_edges, sinks, next_vertex, node_spec, tol_spec, pro_masses, tag_edges, mode);
    }
  }
  else if (vertex == src || tag_end_index >= 0) // between tag.
  {
    std::set<Size> reachable_vertices;

    for (int tag_index = 0; tag_index < tag_node_starts.size(); tag_index++) // for all reachable tag starting point, run extension
    {
      int node_start = tag_node_starts[tag_index];
      int pro_start = tag_pro_starts[tag_index];

      if (node_index > node_start || pro_index > pro_start) continue;
      bool is_visited_start = false;

      int mod_num = getModNumber_(vertex);
      for (int mn = 0; mn <= mod_num; mn++)
        if (visited_tag_edges.find(getVertex_(node_start, pro_start, 0, mn, pro_masses.size())) != visited_tag_edges.end())
        {
          is_visited_start = true;
          break;
        }

      if (is_visited_start) continue;
      visited_tag_edges.insert(getVertex_(node_start, pro_start, 0, mod_num, pro_masses.size()));

      std::set<Size> next_vertices;
      extendBetweenTags_(dag, visited, next_vertices, vertex, node_start, pro_start, 0, node_spec, tol_spec, pro_masses, mode);

      // std::cout <<node_start << " and " << pro_start << ": " << node_mass_map_[mode][node_index] << " to " <<  node_mass_map_[mode][node_start] <<
      // " " << pro_masses_[pro_index] << " to " <<  pro_masses_[pro_start] << " " << next_vertices.size() << std::endl;

      for (auto next_vertex : next_vertices)
      {
        if (next_vertex == src) continue;
        connectBetweenTags_(dag, visited, visited_tag_edges, sinks, next_vertex, node_spec, tol_spec, pro_masses, tag_edges, mode);
        reachable_vertices.insert(next_vertex);
      }
    }
    if (vertex != src && reachable_vertices.empty())
    {
      extendBetweenTags_(dag, visited, sinks, vertex,  -1, -1, 1e5, node_spec, tol_spec, pro_masses, mode);
    }
  }
}

void FLASHExtenderAlgorithm::extendBetweenTags_(FLASHHelperClasses::DAG& dag,
                                                boost::dynamic_bitset<>& visited,
                                                std::set<Size>& sinks,
                                                Size vertex,
                                                int node_index,
                                                int pro_index,
                                                int diagonal_counter,
                                                const MSSpectrum& node_spec,
                                                const MSSpectrum& tol_spec,
                                                const std::vector<double>& pro_masses,
                                                int mode)
{
  if (! visited[vertex]) return;
  int node_index1 = getNodeIndex_(vertex, pro_masses.size());
  int pro_index1 = getProIndex_(vertex, pro_masses.size());
  int score1 = getScore_(vertex);
  int num_mod1 = getModNumber_(vertex);
  if (num_mod1 == max_mod_cntr_) diagonal_counter = 1e5;
  double node_mass1 = node_spec[node_index1].getMZ();

  if (vertex == getVertex_(0, 0, 0, 0, pro_masses.size()))
  {
    for (int pro_i = pro_index1 + 1; pro_i <= pro_index; pro_i++)
    {
      Size vertex2 = getVertex_(node_index1, pro_i, score1, num_mod1, pro_masses.size());
      dag.addEdge(vertex2, vertex, visited);
      if (! visited[vertex2]) continue;
      extendBetweenTags_(dag, visited, sinks, vertex2, node_index, pro_index, diagonal_counter, node_spec, tol_spec, pro_masses, mode);
    }
  }

  double delta_mass1 = node_mass1 - pro_masses[pro_index1];
  double delta_mass;
  double margin;

  if (node_index >= 0)
  {
    double node_mass = node_spec[node_index].getMZ();
    delta_mass = node_mass - pro_masses[pro_index];
    margin = 2 * tol_spec[node_index].getIntensity(); // give some margin
    if (std::abs(delta_mass - delta_mass1) > max_mod_mass_ * (max_mod_cntr_ - num_mod1) + margin) { return; }
    if (diagonal_counter > 0)
    {
      if (std::abs(delta_mass - delta_mass1) > margin) return;
    }
  }

  for (int score = score1 + 1; score <= max_path_score_; score++)
  {
    // if the starting point has already taken by a higher scoring path, don't go further.
    Size higher_score_vertex = getVertex_(node_index1, pro_index1, score, num_mod1, pro_masses.size());
    if (visited[higher_score_vertex]) return;
  }

  if (node_index < 0)
  {
    sinks.insert(vertex);
    pro_index = std::min(pro_index1 + 50, (int)pro_masses.size() - 1); // if sink is not specified, stretch up to ten amino acids.
  }
  else if (node_index1 > node_index || pro_index1 > pro_index)
    return;

  if (node_index1 == node_index && pro_index1 == pro_index)
  {
    sinks.insert(vertex);
    return;
  }

  for (int node_i = node_index1 + 1; node_i <= (node_index < 0 ? (int)node_spec.size() - 1 : node_index); node_i++)
  {
    if (node_index < 0 && node_spec[node_i].getIntensity() < 0) continue;
    int score = score1 + node_spec[node_i].getIntensity();
    double t_node_mass = node_spec[node_i].getMZ();
    double t_margin = tol_spec[node_i].getIntensity();

    for (int pro_i = pro_index1 + (pro_index1 == 0 ? 0 : 1); pro_i <= pro_index; pro_i++)
    {
      double t_delta_mass = t_node_mass - pro_masses[pro_i];

      if (t_delta_mass - delta_mass1 > max_mod_mass_ + t_margin) continue;
      if (delta_mass1 - t_delta_mass > max_mod_mass_ + t_margin) break;

      int num_mod = num_mod1;
      if (std::abs(t_delta_mass - delta_mass1) > t_margin) { num_mod++; }

      if (diagonal_counter > 0 && num_mod != num_mod1) continue; //
      if (num_mod > max_mod_cntr_) continue;
      if (node_index >= 0)
      {
        if (std::abs(t_delta_mass - delta_mass) > max_mod_mass_ * (max_mod_cntr_ - num_mod) + margin) continue;
      }

      Size next_vertex = getVertex_(node_i, pro_i, score, num_mod, pro_masses.size());

      dag.addEdge(next_vertex, vertex, visited);
      if (! visited[next_vertex]) continue;

      int next_diagonal_counter = diagonal_counter;
      if (diagonal_counter > 0) next_diagonal_counter--;
      else if (num_mod != num_mod1)
        next_diagonal_counter = 1;

      extendBetweenTags_(dag, visited, sinks, next_vertex, node_index, pro_index, next_diagonal_counter, node_spec, tol_spec, pro_masses, mode);
    }
  }
}
} // namespace OpenMS