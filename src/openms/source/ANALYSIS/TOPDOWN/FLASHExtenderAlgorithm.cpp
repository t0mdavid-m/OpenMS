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
inline const int max_hit_count = 10;
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
  defaults_.setValue("max_mass_shift", 300.0, "Maximum mass shift for blind search.");
  defaults_.setValue("num_blind", 1, "Number of blind modification per terminal. Per protein, it doubles.");

  defaults_.setValue("ion_type", std::vector<std::string> {"b", "y"}, "Ion types to consider");
  defaults_.setValidStrings("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

  defaultsToParam_();
}

void FLASHExtenderAlgorithm::updateMembers_()
{
  max_mod_cntr_ = param_.getValue("num_blind");
  max_mod_mass_ = param_.getValue("max_mass_shift");
  // fdr_ = param_.getValue("fdr");
  //  keep_decoy_ = param_.getValue("keep_decoy") == "true";
  // std::vector<std::string> ion_strs = param_.getValue("ion_type");
}

Size FLASHExtenderAlgorithm::getVertex_(int node_index, int pro_index, int score, int num_mod, Size pro_length) const
{
  return ((node_index * (pro_length + 1) + pro_index) * (max_mod_cntr_ + 1) + num_mod) * (max_path_score_ - min_path_score_ + 1)
         + (std::min(max_path_score_, std::max(min_path_score_, score)) - min_path_score_);
}

int FLASHExtenderAlgorithm::getNodeIndex_(Size vertex, Size pro_length) const
{
  return ((int)vertex / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1)) / ((int)pro_length + 1);
}

int FLASHExtenderAlgorithm::getProIndex_(Size vertex, Size pro_length) const
{
  return (((int)vertex / (max_path_score_ - min_path_score_ + 1) / (max_mod_cntr_ + 1))) % ((int)pro_length + 1);
}

int FLASHExtenderAlgorithm::getScore_(Size vertex) const
{
  return (int)vertex % (max_path_score_ - min_path_score_ + 1) + min_path_score_;
}

int FLASHExtenderAlgorithm::getModNumber_(Size vertex) const
{
  return ((int)vertex / (max_path_score_ - min_path_score_ + 1)) % (max_mod_cntr_ + 1);
}

// take the hits. Just calculate the mass of truncated protein. Then add modification masses if they are disjoint. If they overlap and the same mass,
// we have a single one. If they are overlapping but different, we add all of them. the max mod count is also adjusted.
double FLASHExtenderAlgorithm::calculate_precursor_mass_(const ProteinHit& hit,
                                                         int protein_start_position,
                                                         int protein_end_position,
                                                         const std::vector<int>& mod_starts,
                                                         const std::vector<int>& mod_ends,
                                                         const std::vector<double>& mod_masses)
{
  double precursor_mass = -1;
  if (protein_start_position < 0 || protein_end_position < 0) return precursor_mass;
  auto seq = hit.getSequence();

  precursor_mass = 0;
  for (int i = protein_start_position; i < protein_end_position; i++)
  {
    if (seq[i] == 'X') continue;
    precursor_mass += AASequence::fromString(seq[i], true).getMonoWeight(Residue::Internal);
  }
  precursor_mass += Residue::getInternalToFull().getMonoWeight();

  if (mod_starts.empty()) { return precursor_mass; }

  std::vector<std::tuple<int, int, double>> ranges;
  ranges.reserve(mod_starts.size());
  for (Size i = 0; i < mod_starts.size(); i++)
  {
    ranges.emplace_back(mod_starts[i], mod_ends[i], mod_masses[i]);
  }

  std::sort(ranges.begin(), ranges.end(), [](const std::tuple<int, int, double>& left, const std::tuple<int, int, double>& right) {
    return std::get<0>(left) != std::get<0>(right) ? std::get<0>(left) < std::get<0>(right) : std::get<1>(left) < std::get<1>(right);
  });

  double total_mod_masses = 0;
  for (Size i = 1; i < ranges.size(); i++)
  {
    if (std::get<0>(ranges[i]) < std::get<1>(ranges[i - 1]))
    {
      // Overlap found
      if (std::abs(std::get<2>(ranges[i]) - std::get<2>(ranges[i - 1])) < Constants::C13C12_MASSDIFF_U) continue;
    }
    total_mod_masses += std::get<2>(ranges[i - 1]);
  }

  total_mod_masses += std::get<2>(ranges.back());

  return precursor_mass + total_mod_masses;
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
    if (aa == 'X') pro_masses.push_back(pro_masses.back()); // repeat the previous mass
    else
      pro_masses.push_back(pro_masses.back() + AASequence::fromString(aa, true).getMonoWeight(Residue::Internal));
  }
}

void FLASHExtenderAlgorithm::define_nodes_(const FLASHTaggerAlgorithm& tagger,
                                           MSSpectrum& node_spec,
                                           MSSpectrum& tol_spec,
                                           double max_mass,
                                           double precursor_mass,
                                           int mode)
{
  const auto& spec = tagger.getSpectrum();
  // 0 for suffix 1 for prefix 2 for suffix and prefix if precursor mass is available
  MSSpectrum t_node_spec, t_tol_spec;

  t_node_spec.reserve(spec.size() * ion_types_str_.size() + 1);
  t_tol_spec.reserve(spec.size() * ion_types_str_.size() + 1);

  for (const auto& p : spec)
  {
    if (mode == 0)
    {
      for (const auto& shift : suffix_shifts_)
      {
        double mass = p.getMZ() - shift;
        if (mass <= 0 || mass > max_mass + 1) continue;
        t_node_spec.emplace_back(mass, p.getIntensity());
        t_tol_spec.emplace_back(mass, 2 * tol_ * p.getMZ());
      }
    }
    else if (mode == 1)
    {
      for (const auto& shift : prefix_shifts_)
      {
        double mass = p.getMZ() - shift;
        if (mass <= 0 || mass > max_mass + 1) continue;
        t_node_spec.emplace_back(mass, p.getIntensity());
        t_tol_spec.emplace_back(mass, 2 * tol_ * p.getMZ());
      }
    }
    else if (mode == 2 && precursor_mass > 0)
    {
      for (const auto& shift : prefix_shifts_)
      {
        double mass = p.getMZ() - shift;
        if (mass <= 0 || mass >= precursor_mass - Residue::getInternalToFull().getMonoWeight()) continue;
        t_node_spec.emplace_back(mass, p.getIntensity());
        t_tol_spec.emplace_back(mass, 2 * tol_ * p.getMZ());
      }
      for (const auto& shift : suffix_shifts_)
      {
        double mass = p.getMZ() - shift;
        if (mass <= 0 || mass >= precursor_mass - Residue::getInternalToFull().getMonoWeight()) continue;
        t_node_spec.emplace_back(precursor_mass - Residue::getInternalToFull().getMonoWeight() - mass, p.getIntensity());
        t_tol_spec.emplace_back(precursor_mass - Residue::getInternalToFull().getMonoWeight() - mass, 2 * tol_ * p.getMZ());
      }
    }
  }
  if (precursor_mass > 0)
  {
    t_node_spec.emplace_back(precursor_mass - Residue::getInternalToFull().getMonoWeight(), 1);
    t_tol_spec.emplace_back(precursor_mass - Residue::getInternalToFull().getMonoWeight(), 2 * tol_ * precursor_mass);
  }

  t_node_spec.sortByPosition();
  t_tol_spec.sortByPosition();
  // Assign the sorted values back to the original vectors
  node_spec.reserve(t_node_spec.size() + 1);
  tol_spec.reserve(t_node_spec.size() + 1);

  node_spec.emplace_back(0, 0);
  tol_spec.emplace_back(0, 0);

  for (Size k = 0; k < t_node_spec.size(); k++)
  {
    const auto& p = t_node_spec[k];
    double mass = p.getMZ();
    float score = p.getIntensity();
    double margin = k == 0 ? .0 : std::max(t_tol_spec[k].getIntensity(), t_tol_spec[k - 1].getIntensity());
    // if (score > 0) max_path_score_ += score;
    // else min_path_score_ += score;

    if (mass - node_spec.back().getMZ() < margin) // they are the same
    {
      float prev_score = node_spec.back().getIntensity();
      margin += mass - node_spec.back().getMZ();
      mass = (mass + node_spec.back().getMZ()) / 2.0;
      if (node_spec.size() > 1)
      {
        node_spec.pop_back();
        tol_spec.pop_back();
      }

      score = 1 + std::max(prev_score, score);
    }

    node_spec.emplace_back(mass, score);
    tol_spec.emplace_back(mass, margin);
  }
  node_spec.sortByPosition();
  tol_spec.sortByPosition();
}

void FLASHExtenderAlgorithm::run_(const ProteinHit& hit,
                                  const MSSpectrum& node_spec,
                                  const MSSpectrum& tol_spec,
                                  const std::vector<double>& pro_masses,
                                  boost::dynamic_bitset<>& visited,
                                  const double precursor_mass,
                                  std::map<int, std::vector<Size>>& all_paths_per_mode,
                                  int mode) // per hit
{
  std::vector<std::vector<int>> tag_edges(4); // pro start end indices and node start end indices
  std::vector<FLASHHelperClasses::Tag> tags;
  std::set<Size> sinks;

  const std::vector<int>& indices = hit.getMetaValue("TagIndices");
  tags.reserve(indices.size());
  int max_tag_length = 3;
  for (auto i : indices)
  {
    max_tag_length = std::max(max_tag_length, (int)tags_[i].getLength());
  }
  for (auto i : indices)
  {
    if (tags_[i].getLength() < std::min(4, max_tag_length - 3)) continue;
    tags.push_back(tags_[i]);
  }

  for (auto& edge : tag_edges)
    edge.reserve(tags.size() * 2);

  FLASHHelperClasses::DAG dag((1 + node_spec.size()) * (1 + pro_masses.size()) * (1 + max_mod_cntr_) * (1 + max_path_score_ - min_path_score_));
  for (const auto& tag : tags)
  {
    // std::cout<<mode<< " "  << tags.size() << " " << tag.getNtermMass() << " "  <<tag.getCtermMass() << " "  <<precursor_mass<< " " <<
    // tag.toString() << " " << hit.getAccession() << std::endl;
    if (mode == 0 && tag.getCtermMass() < 0) continue;
    if (mode == 1 && tag.getNtermMass() < 0) continue;
    if (mode == 2 && precursor_mass <= 0) continue;

    std::vector<int> positions;
    std::vector<double> masses;

    FLASHTaggerAlgorithm::getMatchedPositionsAndFlankingMassDiffs(positions, masses, flanking_mass_tol_, hit, tag);
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

        if (mode == 2)
        {
          start_mass = precursor_mass - Residue::getInternalToFull().getMonoWeight() - start_mass;
          end_mass = precursor_mass - Residue::getInternalToFull().getMonoWeight() - end_mass;
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

    for (Size l = 0; l < start_masses.size(); l++)
    {
      double delta_start = start_tols[l];
      double delta_end = end_tols[l];

      int highest_score_start = node_spec.findHighestInWindow(start_masses[l], delta_start, delta_start);
      int highest_score_end = node_spec.findHighestInWindow(end_masses[l], delta_end, delta_end);

      if (highest_score_start < 0 || highest_score_end < 0 || highest_score_start >= (int)node_spec.size()
          || highest_score_end >= (int)node_spec.size())
        continue;
      for (int pos : positions)
      {
        if (mode == 0) // suffix inverted
        {
          pos = (int)pro_masses.size() - 1 - pos; // invert pos
          if (pos - tag.getLength() >= 0 && pos < (int)pro_masses.size())
          {
            tag_edges[0].emplace_back(pos - tag.getLength());
            tag_edges[1].emplace_back(pos);
            tag_edges[2].emplace_back(highest_score_start);
            tag_edges[3].emplace_back(highest_score_end); // check...
          }
        }
        else
        {
          if (pos >= 0 && pos + tag.getLength() < pro_masses.size())
          {
            tag_edges[0].emplace_back(pos);
            tag_edges[1].emplace_back(pos + tag.getLength()); // this can be much faster...
            tag_edges[2].emplace_back(highest_score_start);
            tag_edges[3].emplace_back(highest_score_end);
          }
        }
      }
    }

    // if (tag_edges.empty() || tag_edges[0].empty()) continue;

    //    if (! tag_edges.empty() && ! tag_edges[0].empty())
    //    {
    //      std::cout << hit.getAccession() << " " << mode << " " << tag.toString() << std::endl;
    //      auto seq = hit.getSequence();
    //      if (mode == 0) seq = seq.reverse();
    //      std::cout << seq.substr(tag_edges[0].back(), tag_edges[1].back() - tag_edges[0].back()) << "\n";
    //
    //
    //      std::cout << node_spec[tag_edges[2].back()].getMZ() << " * " << node_spec[tag_edges[3].back()].getMZ();  std::cout << std::endl;
    //      std::cout << pro_masses[tag_edges[0].back()] << " ** " << pro_masses[tag_edges[1].back()] << " : " << tag_edges[0].back() << "  " <<
    //      tag_edges[1].back(); std::cout << std::endl;
    //
    //      std::cout << pro_masses[tag_edges[1].back()] - pro_masses[tag_edges[0].back()] << " vs. "
    //                << node_spec[tag_edges[3].back()].getMZ() - node_spec[tag_edges[2].back()].getMZ();
    //
    //      std::cout << std::endl;
    //    }
  }
  constructDAG_(dag, sinks, visited, node_spec, tol_spec, pro_masses, tag_edges, mode, ! tags.empty());
  Size src = getVertex_(0, 0, 0, 0, pro_masses.size());
  std::vector<int> max_scores(max_mod_cntr_ + 1, 0);
  for (Size sink : sinks)
  {
    int num_mod = getModNumber_(sink);
    if (sink == src || getScore_(sink) < max_scores[num_mod]) continue; // getModNumber_(sink) == 0 ||
    if (precursor_mass > 0 && getNodeIndex_(sink, pro_masses.size()) < (int)node_spec.size() - 1) continue;
    max_scores[num_mod] = getScore_(sink);
  }

  std::vector<std::vector<std::vector<Size>>> paths(max_mod_cntr_ + 1, std::vector<std::vector<Size>>());
  for (Size sink : sinks)
  {
    int num_mod = getModNumber_(sink);
    if (sink == src || getScore_(sink) < max_scores[num_mod]) continue; //
    if (precursor_mass > 0 && getNodeIndex_(sink, pro_masses.size()) + 1 < (int)node_spec.size()) continue;
    dag.findAllPaths(sink, src, paths[num_mod], 20);
  }

  for (int num_mod = 0; num_mod <= max_mod_cntr_; num_mod++)
  {
    for (const auto& path : paths[num_mod])
    {
      auto iter = all_paths_per_mode.find(num_mod);
      // std::cout<<mode << " "   << num_mod << " "   << getScore_(path[0])<<std::endl;
      if (iter == all_paths_per_mode.end() || iter->second.size() < path.size()) all_paths_per_mode[num_mod] = path;
    }
  }
}

void FLASHExtenderAlgorithm::run(const FLASHTaggerAlgorithm& tagger, double flanking_mass_tol, double ppm, bool multiple_hits_per_spec)
{
  // setLogType(CMD);
  ion_types_str_ = param_.getValue("ion_type").toStringVector();

  Residue empty;
  for (const auto& ion_str : ion_types_str_)
  {
    if (ion_str == "a") { prefix_shifts_.push_back(Residue::getInternalToAIon().getMonoWeight()); }
    else if (ion_str == "b") { prefix_shifts_.push_back(Residue::getInternalToBIon().getMonoWeight()); }
    else if (ion_str == "c") { prefix_shifts_.push_back(Residue::getInternalToCIon().getMonoWeight()); }
    else if (ion_str == "x") { suffix_shifts_.push_back(Residue::getInternalToXIon().getMonoWeight()); }
    else if (ion_str == "y") { suffix_shifts_.push_back(Residue::getInternalToYIon().getMonoWeight()); }
    else if (ion_str == "z") { suffix_shifts_.push_back(Residue::getInternalToZIon().getMonoWeight()); }
    else if (ion_str == "zp1") { suffix_shifts_.push_back(Residue::getInternalToZp1Ion().getMonoWeight()); }
    else if (ion_str == "zp2") { suffix_shifts_.push_back(Residue::getInternalToZp2Ion().getMonoWeight()); }
    else { continue; }
  }

  flanking_mass_tol_ = flanking_mass_tol;
  tol_ = ppm / 1e6;
  tags_.clear();
  proteoform_hits_.clear();

  tagger.getTags(tags_);
  if (tags_.empty()) return;
  std::vector<double> mzs;
  std::vector<int> scores;
  std::vector<ProteinHit> hits;
  tagger.getProteinHits(hits, max_hit_count);
  if (hits.empty()) return;
  proteoform_hits_.reserve(hits.size());
  if (tagger.getSpectrum().metaValueExists("PrecursorMass")) { precursor_mass_ = tagger.getSpectrum().getMetaValue("PrecursorMass"); }
  startProgress(0, (int)hits.size(), "running FLASHExtender ...");

#pragma omp parallel for default(none) shared(hits, tagger, multiple_hits_per_spec, std::cout)
  for (Size i = 0; i < hits.size(); i++)
  {
    nextProgress();
    auto& hit = hits[i];

    int total_score = 0;
    int total_match_cntr = 0;
    std::vector<int> mod_starts, mod_ends;
    std::vector<double> mod_masses;
    int protein_start_position = -1, protein_end_position = -1;
    int max_nterm_index = 0, max_cterm_rindex = 0;

    boost::dynamic_bitset<> visited;
    double precursor_mass = -1;

    std::map<int, std::map<int, std::vector<Size>>> all_path_map; // mode, num_mod, path
    std::map<int, std::vector<Size>> best_path_map;               // mode, best paths
    std::map<int, MSSpectrum> node_spec_map;          // mode, best paths
    std::map<int, std::vector<double>> pro_mass_map;
    for (int mode = 0; mode < 3; mode++) //
    {
      if (mode == 2 && precursor_mass <= 0)
      {
        if (max_nterm_index + max_cterm_rindex >= (int)hit.getSequence().size())
          precursor_mass = calculate_precursor_mass_(hit, protein_start_position, protein_end_position, mod_starts, mod_ends, mod_masses);
        if (precursor_mass <= 0) precursor_mass = precursor_mass_;
        if (precursor_mass <= 0) break;
      }

      auto& pro_masses = pro_mass_map[mode] = std::vector<double>();
      auto& node_spec = node_spec_map[mode] = MSSpectrum();
      MSSpectrum tol_spec;

      get_pro_masses_(hit, pro_masses, mode);
      define_nodes_(tagger, node_spec, tol_spec, pro_masses.back(), precursor_mass, mode);
      if (visited.empty())
        visited = boost::dynamic_bitset<>((1 + tagger.getSpectrum().size() * ion_types_str_.size()) * (1 + pro_masses.size()) * (2 + max_mod_cntr_)
                                          * (1 + max_path_score_ - min_path_score_));
      run_(hit, node_spec, tol_spec, pro_masses, visited, precursor_mass, all_path_map[mode], mode);

      if (mode < 2)
      {
        const auto paths_c = all_path_map.find(0);
        const auto paths_n = all_path_map.find(1);

        std::map<int, int> nscores, cscores; // mod and score

        if (paths_n != all_path_map.end())
        {
          for (const auto& [mod, path] : paths_n->second)
          {
            if (path.empty()) continue;
            nscores[mod] = getScore_(path[0]);
          }
        }
        if (paths_c != all_path_map.end())
        {
          for (const auto& [mod, path] : paths_c->second)
          {
            if (path.empty()) continue;
            cscores[mod] = getScore_(path[0]);
          }
        }
        int max_score = 0;

        if (nscores.empty()) // only c term
        {
          for (const auto& [mod, path] : paths_c->second)
          {
            if (path.empty()) continue;
            if (max_score >= getScore_(path[0])) continue;
            max_score = getScore_(path[0]);
            best_path_map[0] = path;
          }
        }
        else if (cscores.empty()) // only n term
        {
          for (const auto& [mod, path] : paths_n->second)
          {
            if (path.empty()) continue;
            if (max_score >= getScore_(path[0])) continue;
            max_score = getScore_(path[0]);
            best_path_map[1] = path;
          }
        }
        else // both terms
        {
          for (int mc = 0; mc <= max_mod_cntr_; mc++)
          {
            if (cscores.find(mc) == cscores.end()) continue;
            for (int mn = 0; mc + mn <= max_mod_cntr_; mn++)
            {
              if (nscores.find(mn) == nscores.end()) continue;
              int sum_score = nscores[mn] + cscores[mc];

              if (max_score >= sum_score) continue;
              max_score = sum_score;
              best_path_map[0] = paths_c->second[mc];
              best_path_map[1] = paths_n->second[mn];
            }
          }
        }
      }
      else if (mode == 2)
      {
        int max_score = 0;
        const auto paths = all_path_map.find(2);
        for (const auto& [mod, path] : paths->second)
        {
          if (path.empty()) continue;
          if (max_score >= getScore_(path[0])) continue;
          max_score = getScore_(path[0]);
          best_path_map[2] = path;
        }

        mod_starts.clear();
        mod_ends.clear();
        mod_masses.clear();
      }
      if (mode == 0) continue;
      total_score = 0;
      total_match_cntr = 0;

      for (int m = mode == 2 ? 2 : 0; m <= mode; m++)
      {
        if (best_path_map.empty() || best_path_map.find(m) == best_path_map.end() || best_path_map[m].empty()) continue;
        auto& best_path = best_path_map[m];
        auto& t_pro_masses = pro_mass_map[m];
        auto& t_node_spec = node_spec_map[m];

        total_score += getScore_(best_path[0]);
        total_match_cntr += (int)best_path.size() - 1;

        double prev_mass_shift = 0;
        int prev_mod_count = 0;
        int pre_pro_index = 0;
        for (auto iter = best_path.rbegin(); iter != best_path.rend(); iter++)
        {
          auto pro_index = getProIndex_(*iter, t_pro_masses.size());
          auto node_index = getNodeIndex_(*iter, t_pro_masses.size());
          auto mass_shift = t_node_spec[node_index].getMZ() - t_pro_masses[pro_index];
          auto mod_count = getModNumber_(*iter);

//          std::cout << hit.getAccession() << " mode " << m << " input pre " << precursor_mass_ << " cal pre " << precursor_mass << "score"
//                    << getScore_(*iter) << " " << node_index << " " << pro_index << " in " <<  t_node_spec.size() << " " <<
//            t_pro_masses.size()<< " masses " << t_pro_masses[pro_index] << " " << t_node_spec[node_index].getMZ() << " " << t_pro_masses[pro_index]  -  t_node_spec[node_index].getMZ()
//                    << " " << mod_count << std::endl;

          if (node_index == 0)
          {
            if (m > 0) protein_start_position = pro_index;
            if (m == 0) protein_end_position = (int)hit.getSequence().size() - pro_index;
          }

          if (m == 0) max_cterm_rindex = std::max(max_cterm_rindex, pro_index);
          if (m == 1) max_nterm_index = std::max(max_nterm_index, pro_index);

          if (mod_count != prev_mod_count)
          {
            double mod_mass = mass_shift - prev_mass_shift;
            int start = m > 0 ? pre_pro_index : (int)hit.getSequence().size() - 1 - pro_index;
            int end = m > 0 ? pro_index : (int)hit.getSequence().size() - 1 - pre_pro_index;
            mod_starts.push_back(start);
            mod_ends.push_back(end);
            mod_masses.push_back(mod_mass);
          }
          prev_mod_count = mod_count;
          pre_pro_index = pro_index;
          prev_mass_shift = mass_shift;
        }
      }
    }
    if (total_score <= 0) continue;

    // std::cout<<hit.getAccession()<< " " << total_score << " " << hit.getScore() <<std::endl;
    //  if (precursor_mass < 0) continue; // if precursor mass is not specified, skip
    hit.setMetaValue("Modifications", mod_masses);
    hit.setMetaValue("ModificationStarts", mod_starts);
    hit.setMetaValue("ModificationEnds", mod_ends);
    hit.setMetaValue("MatchedAA", total_match_cntr);
    hit.setCoverage((double)total_match_cntr / (double)hit.getSequence().size());
    hit.setScore(total_score);
    hit.setMetaValue("StartPosition", protein_start_position);
    hit.setMetaValue("EndPosition", protein_end_position);
    hit.setMetaValue("Mass", precursor_mass < 0 ? precursor_mass_ : precursor_mass);
    hit.setMetaValue("RT", tagger.getSpectrum().getRT());
    // hit.setMetaValue("Proforma", string)
#pragma omp critical
    {
      bool insert = true;

      if (! multiple_hits_per_spec && ! proteoform_hits_.empty()) // when multiple hits are not allowed
      {
        if (proteoform_hits_.back().getScore() >= hit.getScore()) insert = false;
        else
          proteoform_hits_.pop_back();
      }

      if (insert) { proteoform_hits_.push_back(hit); }
    }
  }
  endProgress();
}

void FLASHExtenderAlgorithm::constructDAG_(FLASHHelperClasses::DAG& dag,
                                           std::set<Size>& sinks,
                                           boost::dynamic_bitset<>& visited,
                                           const MSSpectrum& node_spec,
                                           const MSSpectrum& tol_spec,
                                           const std::vector<double>& pro_masses,
                                           const std::vector<std::vector<int>>& tag_edges,
                                           int mode,
                                           bool use_tags)
{
  // if (tag_edges.empty()) return;
  if (visited.size() < dag.size()) visited.resize(dag.size());
  visited.reset();
  Size src = getVertex_(0, 0, 0, 0, pro_masses.size());
  visited[src] = true;
  std::set<Size> visited_tag_edges;

  connectBetweenTags_(dag, visited, visited_tag_edges, sinks, src, node_spec, tol_spec, pro_masses, tag_edges, mode, use_tags);
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
                                                 int mode,
                                                 bool use_tags)
{
  int node_index = getNodeIndex_(vertex, pro_masses.size());
  int pro_index = getProIndex_(vertex, pro_masses.size());

  int tag_start_index = -1;
  int tag_end_index = -1;

  const auto& tag_pro_starts = tag_edges[0];
  const auto& tag_pro_ends = tag_edges[1];
  const auto& tag_node_starts = tag_edges[2];
  const auto& tag_node_ends = tag_edges[3];
  // use_tags = !tag_pro_starts.empty();
  for (int i = 0; i < (int)tag_pro_starts.size(); i++)
  {
    if (tag_start_index < 0 && tag_node_starts[i] == node_index && tag_pro_starts[i] == pro_index) { tag_start_index = i; }
    if (tag_end_index < 0 && tag_node_ends[i] == node_index && tag_pro_ends[i] == pro_index) { tag_end_index = i; }
  }

  Size src = getVertex_(0, 0, 0, 0, pro_masses.size());

  if (tag_start_index >= 0) // within tag
  {
    int node_end = -1;
    Size i = tag_start_index;
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
      connectBetweenTags_(dag, visited, visited_tag_edges, sinks, next_vertex, node_spec, tol_spec, pro_masses, tag_edges, mode, use_tags);
    }
  }

  if (vertex == src || tag_end_index >= 0) // between tag.
  {
    std::set<Size> reachable_vertices;

    for (Size tag_index = 0; tag_index < tag_node_starts.size(); tag_index++) // for all reachable tag starting point, run extension
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

      for (auto next_vertex : next_vertices)
      {
        if (next_vertex == src) continue;
        connectBetweenTags_(dag, visited, visited_tag_edges, sinks, next_vertex, node_spec, tol_spec, pro_masses, tag_edges, mode, use_tags);
        reachable_vertices.insert(next_vertex);
      }
    }
    if ((vertex != src || ! use_tags) && reachable_vertices.empty())
    {
      extendBetweenTags_(dag, visited, sinks, vertex, -1, -1, mode < 2 && use_tags ? 1e5 : 0, node_spec, tol_spec, pro_masses, mode);
      // if (mode == 2) extendBetweenTags_(dag, visited, sinks, vertex,  node_spec.size() - 1, pro_masses.size() - 1, 1e5, node_spec, tol_spec,
      // pro_masses, mode);
    }
  }
}

void FLASHExtenderAlgorithm::extendBetweenTags_(FLASHHelperClasses::DAG& dag,
                                                boost::dynamic_bitset<>& visited,
                                                std::set<Size>& sinks,
                                                Size vertex,
                                                const int node_index,
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
  auto src = getVertex_(0, 0, 0, 0, pro_masses.size());

  if (node_index < 0)
  {
    pro_index = mode < 2 ? std::min(pro_index1 + 10, (int)pro_masses.size() - 1)
                         : ((int)pro_masses.size() - 1); // if sink is not specified, stretch up to 10 amino acids.
  }

  if (vertex == src)
  {
    for (int pro_i = pro_index1 + 1; pro_i <= pro_index; pro_i++)
    {
      Size vertex2 = getVertex_(0, pro_i, 0, 0, pro_masses.size());
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
    margin = 2 * std::max(tol_spec[node_index].getIntensity(), tol_spec[node_index1].getIntensity()); // give some margin
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
    if (vertex != src) sinks.insert(vertex);
  }

  else if (node_index1 > node_index || pro_index1 > pro_index)
    return;

  if (node_index1 == node_index && pro_index1 == pro_index && vertex != src)
  {
    sinks.insert(vertex);
    return;
  }
  for (int node_i = node_index1 + 1; node_i <= (node_index < 0 ? ((int)node_spec.size() - 1) : node_index); node_i++)
  {
    if (node_index < 0 && node_spec[node_i].getIntensity() < 0 && diagonal_counter > 0) continue;
    int score = score1 + (int)node_spec[node_i].getIntensity();
    double t_node_mass = node_spec[node_i].getMZ();
    double t_margin = std::max(tol_spec[node_i].getIntensity(), tol_spec[node_index1].getIntensity());

    for (int pro_i = pro_index1 + (pro_index1 == 0 ? 0 : 1); pro_i <= pro_index; pro_i++) //
    {
      double t_delta_mass = t_node_mass - pro_masses[pro_i];
      double delta_delta = t_delta_mass - delta_mass1;
      if (delta_delta > max_mod_mass_ + t_margin) continue;
      if (-delta_delta > max_mod_mass_ + t_margin) break;

      int num_mod = num_mod1;
      if (std::abs(delta_delta) > t_margin)
      {
        if (abs(delta_delta) < t_margin + 0.036386) continue;
        num_mod++;
      }

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