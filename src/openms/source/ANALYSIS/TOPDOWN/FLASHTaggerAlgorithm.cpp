// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <utility>

namespace OpenMS
{
inline const Size max_node_cntr = 400;

std::vector<Residue> FLASHTaggerAlgorithm::getAA_(double l, double r, double tol, int iso_offset) const
{
  std::vector<Residue> ret;
  if (l == r) return ret;
  double iso_mass = std::abs(iso_offset * Constants::C13C12_MASSDIFF_U);
  double diff1 = std::abs(std::abs(r - l) - iso_mass);
  double diff2 = std::abs(std::abs(r - l) + iso_mass);
  double abs_tol = std::max(l, r) * tol / 1e6 * 2;
  auto iter = aa_mass_map_.lower_bound(diff1 - abs_tol);

  while (iter != aa_mass_map_.end())
  {
    if (std::abs(diff1 - iter->first) < abs_tol || std::abs(diff2 - iter->first) < abs_tol)
    {
      for (auto& aa : iter->second)
        ret.push_back(aa);
    }
    else if (iter->first - diff2 > abs_tol) { break; }
    iter++;
  }
  return ret;
}

std::vector<std::vector<Residue>> FLASHTaggerAlgorithm::getGap_(double l, double r, double tol, int iso_offset) const
{
  std::vector<std::vector<Residue>> ret;
  if (l == r) return ret;
  double iso_mass = std::abs(iso_offset * Constants::C13C12_MASSDIFF_U);

  double diff1 = std::abs(std::abs(r - l) - iso_mass);
  double diff2 = std::abs(std::abs(r - l) + iso_mass);
  double abs_tol = 2 * std::abs(std::max(l, r) * tol / 1e6);
  auto iter = gap_mass_map_.lower_bound(diff1 - abs_tol);

  while (iter != gap_mass_map_.end())
  {
    if (std::abs(diff1 - iter->first) < abs_tol || std::abs(diff2 - iter->first) < abs_tol)
    {
      for (auto& aa : iter->second)
        ret.push_back(aa);
    }
    else if (iter->first - diff2 > abs_tol) { break; }
    iter++;
  }
  return ret;
}

void FLASHTaggerAlgorithm::updateEdgeMasses_()
{
  aa_mass_map_.clear();
  gap_mass_map_.clear();

  for (const auto& aa : aas_)
  {
    double aa_mass = aa->getMonoWeight(Residue::Internal);
    if (aa_mass_map_.find(aa_mass) == aa_mass_map_.end()) aa_mass_map_[aa_mass] = std::vector<Residue>();
    aa_mass_map_[aa_mass].push_back(*aa);
  }

  if (max_gap_count_ > 0)
  {
    std::map<double, std::vector<std::vector<Residue>>> prev_gap_mass_map_;
    prev_gap_mass_map_[.0] = std::vector<std::vector<Residue>>(1, std::vector<Residue>());
    for (int i = 0; i < max_aa_in_gap_; i++)
    {
      for (const auto& prev : prev_gap_mass_map_)
      {
        for (const auto& current : aa_mass_map_)
        {
          if (gap_mass_map_.find(prev.first + current.first) == gap_mass_map_.end())
            gap_mass_map_[prev.first + current.first] = std::vector<std::vector<Residue>>();

          for (const auto& aa_vec : prev.second)
          {
            for (const auto& aa : current.second)
            {
              auto new_aa_vec(aa_vec);
              new_aa_vec.push_back(aa);
              gap_mass_map_[prev.first + current.first].push_back(new_aa_vec);
            }
          }
        }
      }
      prev_gap_mass_map_ = gap_mass_map_;
    }

    std::map<double, std::vector<std::vector<Residue>>> tmp_gap_mass_map;
    for (auto& e : gap_mass_map_)
    {
      std::vector<std::vector<Residue>> tmp_e;
      for (auto& f : e.second)
      {
        if (f.size() <= 1) continue;
        tmp_e.push_back(f);
      }
      if (tmp_e.empty()) continue;
      tmp_gap_mass_map[e.first] = tmp_e;
    }
    gap_mass_map_ = tmp_gap_mass_map;
  }
}

Size FLASHTaggerAlgorithm::getVertex_(int index, int path_score, int level, int iso_level, int gap_level) const
{
  return (((index * (max_tag_length_ + 1) + level) * (max_iso_in_tag_ + 1) + iso_level) * (max_gap_count_ + 1) + gap_level)
           * (max_path_score_ - min_path_score_ + 1)
         + (path_score - min_path_score_);
}

int FLASHTaggerAlgorithm::getIndex_(Size vertex) const
{
  return (((int)vertex / (max_path_score_ - min_path_score_ + 1)) / (max_gap_count_ + 1)) / (max_iso_in_tag_ + 1) / (max_tag_length_ + 1);
}

void FLASHTaggerAlgorithm::constructDAG_(FLASHHelperClasses::DAG& dag,
                                         const std::vector<double>& mzs,
                                         const std::vector<int>& scores,
                                         int length,
                                         double tol)
{
  // from source to sink, connect but the edge direction is from sink to source.
  edge_aa_map_.clear();
  int start_index = 1; // zeroth = source.
  int end_index = 1;
  boost::dynamic_bitset<> visited(dag.size());
  visited[getVertex_(0, 0, 0, 0, 0)] = true;

  while (end_index < (int)mzs.size())
  {
    auto r = mzs[end_index];

    // first, make edge from r to source.
    Size vertex1 = getVertex_(end_index, scores[end_index], 0, 0, 0);
    Size vertex2 = getVertex_(0, 0, 0, 0, 0);

    dag.addEdge(vertex1, vertex2, visited); // move to DAG?

    // from an edge i, j to class edge.  for each i, j make a unique key. key to an edge.

    while (start_index < end_index && r - mzs[start_index] > max_edge_mass_)
      start_index++;

    for (int g = 0; g < 2; g++) // 0 for only a.a 1 for with gap.
    {
      for (int n = 0; n < 2; n++) // 0 for all a.a 1 for isotope errors. Allow only one isotope errors.
      {
        for (int current_index = start_index; current_index < end_index; current_index++)
        {
          auto l = mzs[current_index];
          int edge_score = scores[end_index];

          // make edge from r to l if they make an a.a. mass.
          std::vector<Residue> aas;
          std::vector<std::vector<Residue>> gaps;
          if (g == 0)
          {
            aas = getAA_(l, r, tol, n);
            if (aas.empty()) continue;
          }
          else
          {
            gaps = getGap_(l, r, tol, n);
            if (gaps.empty()) continue;
          }
          // end_index, current_index to amino acid strings.
          if (edge_aa_map_.find(end_index) == edge_aa_map_.end()) { edge_aa_map_[end_index] = std::map<int, std::vector<String>>(); }
          auto& e = edge_aa_map_[end_index];

          if (e.find(current_index) == e.end()) { e[current_index] = std::vector<String>(); }

          if (g == 0)
          {
            for (auto& aa : aas)
            {
              auto aaStr = n == 0 ? aa.toString() : aa.toString().toLower();
              e[current_index].push_back(aaStr);
            }
          }
          else
          {
            for (auto& gap : gaps)
            {
              std::stringstream ss;
              for (auto& aa : gap)
                ss << aa.toString().toLower();
              e[current_index].emplace_back(ss.str());
            }
          }

          int gap_diff = g == 0 ? 0 : 1;

          for (int d = 0; d + gap_diff <= max_gap_count_; d++)
          {
            for (int iso = 0; iso + n <= max_iso_in_tag_; iso++)
            {
              for (int lvl = 0; lvl < length; lvl++)
              {
                for (int score = min_path_score_; score <= max_path_score_; score++)
                {
                  if (score - edge_score < min_path_score_) continue;
                  if (score - edge_score > max_path_score_) break;

                  vertex1 = getVertex_(end_index, score, lvl + 1, iso + n, d + gap_diff);
                  vertex2 = getVertex_(current_index, score - edge_score, lvl, iso, d);
                  dag.addEdge(vertex1, vertex2, visited);
                }
              }
            }
          }
        }
        if (max_iso_in_tag_ == 0) break;
      }
      if (max_gap_count_ == 0) break;
    }
    //  make edge from sink to r
    if (end_index < (int)mzs.size() - 1)
    {
      for (int d = 0; d <= max_gap_count_; d++)
      {
        for (int iso = 0; iso <= max_iso_in_tag_; iso++)
        {
          for (int score = min_path_score_; score <= max_path_score_; score++)
          {
            vertex1 = getVertex_((int)mzs.size() - 1, score, length, iso, d);
            vertex2 = getVertex_(end_index, score, length, iso, d);
            dag.addEdge(vertex1, vertex2, visited);
          }
        }
      }
    }
    end_index++;
  }
}

FLASHTaggerAlgorithm::FLASHTaggerAlgorithm(): DefaultParamHandler("FLASHTaggerAlgorithm"), ProgressLogger()
{
  setDefaultParams_();
}

FLASHTaggerAlgorithm::FLASHTaggerAlgorithm(const FLASHTaggerAlgorithm& other): DefaultParamHandler(other), ProgressLogger(other)
{
}

FLASHTaggerAlgorithm& FLASHTaggerAlgorithm::operator=(const FLASHTaggerAlgorithm& rhs)
{
  if (this == &rhs) return *this;

  DefaultParamHandler::operator=(rhs);
  return *this;
}

void FLASHTaggerAlgorithm::setDefaultParams_()
{
  defaults_.setValue("max_count", 300,
                     "Maximum number of the tags per length (lengths set by -min_length and -max_length options). The tags with different amino acid "
                     "combinations but with the same masses are counted once. E.g., "
                     "TII, TIL, TLI, TLL are distinct tags even though they have the same mass differences. "
                     "but are counted as one tag. ");
  defaults_.setMinInt("max_count", 0);

  defaults_.setValue(
    "min_length", 4,
    "Minimum length of a tag. Each mass gap contributes to a single length (even if a mass gap is represented by multiple amino acids). ");
  defaults_.setMaxInt("min_length", 9);
  defaults_.setMinInt("min_length", 2);

  defaults_.setValue(
    "max_length", 8,
    "Maximum length of a tag. Each mass gap contributes to a single length (even if a mass gap is represented by multiple amino acids). ");
  defaults_.setMaxInt("max_length", 30);
  defaults_.setMinInt("max_length", 3);

  defaults_.setValue("flanking_mass_tol", 1000000.0, "Flanking mass tolerance (the flanking mass minus protein mass up to the matching amino acid) in Da. This limits the possible terminal modification mass.");

  defaults_.setValue("max_iso_error_count", 0, "Maximum isotope error count per tag.");

  defaults_.setValue("allow_iso_error", "false", "Allow up to one isotope error in each tag.");
  defaults_.setValidStrings("allow_iso_error", {"true", "false"});
  defaults_.addTag("allow_iso_error", "advanced");

  defaults_.setValue("min_matched_aa", 4, "Minimum number of amino acids in matched proteins, covered by tags.");

  defaults_.setValue("allow_gap", "false", "Allow a mass gap (a mass representing multiple consecutive amino acids) in each tag.");
  defaults_.setValidStrings("allow_gap", {"true", "false"});
  defaults_.addTag("allow_gap", "advanced");

  defaults_.setValue("max_aa_in_gap", 2, "Maximum amino acid count in a mass gap.");
  defaults_.setMaxInt("max_aa_in_gap", 3);
  defaults_.setMinInt("max_aa_in_gap", 2);
  defaults_.addTag("max_aa_in_gap", "advanced");

  defaultsToParam_();
}

void FLASHTaggerAlgorithm::updateMembers_()
{
  max_tag_count_ = param_.getValue("max_count");
  min_tag_length_ = param_.getValue("min_length");
  max_tag_length_ = param_.getValue("max_length");
  max_iso_in_tag_ = param_.getValue("allow_iso_error").toString() == "true" ? 1 : 0;
  min_cov_aa_ = (int)param_.getValue("min_matched_aa");
  max_aa_in_gap_ = param_.getValue("max_aa_in_gap");
  max_gap_count_ = param_.getValue("allow_gap").toString() == "true" ? 1 : 0;
  // prsm_fdr_ = param_.getValue("fdr");
  flanking_mass_tol_ = param_.getValue("flanking_mass_tol");
  updateEdgeMasses_();
  max_edge_mass_ = max_iso_in_tag_ * Constants::C13C12_MASSDIFF_U;
  max_edge_mass_ += gap_mass_map_.empty() ? aa_mass_map_.rbegin()->first : std::max(aa_mass_map_.rbegin()->first, gap_mass_map_.rbegin()->first);
}

const MSSpectrum& FLASHTaggerAlgorithm::getSpectrum() const
{
  return spec_;
}

void FLASHTaggerAlgorithm::run(const DeconvolvedSpectrum& deconvolved_spectrum, double ppm)
{
  // setLogType(CMD);

  if (deconvolved_spectrum.empty() || deconvolved_spectrum.isDecoy() || deconvolved_spectrum.getOriginalSpectrum().getMSLevel() == 1) return;

  auto tags = std::vector<FLASHHelperClasses::Tag>();
  tags.reserve(max_tag_count_ * max_tag_length_);

  spec_.setRT(deconvolved_spectrum.getOriginalSpectrum().getRT());
  if (deconvolved_spectrum.getPrecursorPeakGroup().getMonoMass() > 0)
  {
    spec_.setMetaValue("PrecursorMass", deconvolved_spectrum.getPrecursorPeakGroup().getMonoMass());
  }
  getTags_(deconvolved_spectrum, ppm);

  std::sort(tags_.rbegin(), tags_.rend());

  int index = 0;
  for (auto& tag : tags_)
  {
    tag.setIndex(index++);
    tag.setRetentionTime(deconvolved_spectrum.getOriginalSpectrum().getRT());
  }
}

void FLASHTaggerAlgorithm::getTags_(const DeconvolvedSpectrum& dspec, double ppm)
{
  std::vector<double> mzs;
  std::vector<int> scores;
  mzs.reserve(dspec.size());
  scores.reserve(dspec.size());

  for (auto& pg : dspec)
  {
    int score = // (int)round(5 * log10(std::max(1e-1, pg.getQscore() / std::max(1e-1, (1.0 - pg.getQscore())))));
      (int)round(max_score * pg.getQscore());
    // if (score <= -5) continue;
    scores.push_back(score);
    mzs.push_back(pg.getMonoMass());
  }
  getTags_(mzs, scores, dspec.getScanNumber(), ppm);
}

void FLASHTaggerAlgorithm::updateTagSet_(std::set<FLASHHelperClasses::Tag>& tag_set,
                                         std::map<String, std::vector<FLASHHelperClasses::Tag>>& seq_tag,
                                         const std::vector<Size>& path,
                                         const std::vector<double>& mzs,
                                         const std::vector<int>& scores,
                                         int scan,
                                         double ppm)
{
  double flanking_mass = -1;

  std::vector<String> seqs {""};
  std::vector<double> tag_mzs;
  std::vector<int> tag_scores;
  tag_mzs.reserve(path.size() - 1);
  tag_scores.reserve(path.size() - 1);

  for (int j = 1; j < (int)path.size(); j++)
  {
    int i1 = getIndex_(path[j - 1]); // c term side
    int i2 = getIndex_(path[j]);     // n term side

    if (edge_aa_map_.find(i1) != edge_aa_map_.end() && edge_aa_map_[i1].find(i2) != edge_aa_map_[i1].end())
    {
      auto& edge_aa = edge_aa_map_[i1];
      std::vector<String> tmp_seqs;
      tmp_seqs.reserve(seqs.size() * edge_aa[i2].size());
      for (const auto& tmp_seq : seqs)
      {
        for (const auto& seq : edge_aa[i2])
        {
          tmp_seqs.emplace_back(seq + tmp_seq);
        }
      }
      seqs = tmp_seqs;
      tag_mzs.push_back(mzs[i1]);
      tag_scores.push_back(scores[i1]);
    }
    else if (i2 == 0) // nterm
    {
      tag_mzs.push_back(mzs[i1]);
      tag_scores.push_back(scores[i1]);
      flanking_mass = mzs[i1];
    }
  }

  std::vector<double> rev_tag_mzs;
  rev_tag_mzs.reserve(tag_mzs.size());
  std::vector<int> rev_tag_scores;
  rev_tag_scores.reserve(tag_scores.size());

  for (int i = (int)tag_mzs.size() - 1; i >= 0; i--)
  {
    rev_tag_mzs.push_back(tag_mzs[i]);
    rev_tag_scores.push_back(tag_scores[i]);
  }

  for (const auto& seq : seqs)
  {
    auto iter = seq_tag.find(seq);
    bool pass = true;
    if (iter != seq_tag.end()) // remove overlapping tags.
    {
      for (const auto& pt : iter->second)
      {
        if (pt.getNtermMass() < 0) continue;
        if (abs(pt.getNtermMass() - flanking_mass) / std::max(pt.getNtermMass(), flanking_mass) * 1e6 > ppm) continue;
        pass = false;
        break;
      }
    }
    if (pass)
    {
      auto direct_tag = FLASHHelperClasses::Tag(seq, flanking_mass, -1, tag_mzs, tag_scores, scan);
      tag_set.insert(direct_tag);
      seq_tag[seq].push_back(direct_tag);
    }

    pass = true;
    String rev_seq = String(seq).reverse();
    iter = seq_tag.find(rev_seq);
    if (iter != seq_tag.end()) // remove overlapping tags.
    {
      for (const auto& pt : iter->second)
      {
        if (pt.getCtermMass() < 0) continue;
        if (abs(pt.getCtermMass() - flanking_mass) / std::max(pt.getCtermMass(), flanking_mass) * 1e6 > ppm) continue;
        pass = false;
        break;
      }
    }
    if (pass)
    {
      auto reverse_tag = FLASHHelperClasses::Tag(rev_seq, -1, flanking_mass, rev_tag_mzs, rev_tag_scores, scan);
      tag_set.insert(reverse_tag);
      seq_tag[rev_seq].push_back(reverse_tag);
    }
  }
}

void FLASHTaggerAlgorithm::getTags_(const std::vector<double>& mzs, const std::vector<int>& scores, int scan, double ppm)
{
  if (max_tag_count_ == 0) return;

  std::vector<int> _scores;
  std::vector<double> _mzs;
  int threshold;

  if (mzs.size() >= max_node_cntr)
  {
    _scores = scores;
    std::sort(_scores.rbegin(), _scores.rend());
    threshold = _scores[max_node_cntr - 1];
    _scores.clear();

    _mzs.reserve(max_node_cntr + 1);
    _scores.reserve(max_node_cntr + 1);
  }
  else
  {
    _mzs.reserve(mzs.size() + 1);
    _scores.reserve(mzs.size() + 1);
    threshold = *std::min_element(scores.begin(), scores.end());
  }

  _mzs.push_back(.0);
  _scores.push_back(0);
  spec_.reserve(_mzs.size());
  for (int i = 0; i < (int)mzs.size(); i++)
  {
    if (scores[i] < threshold) continue;
    _mzs.push_back(mzs[i]);
    _scores.push_back(scores[i]);
    spec_.emplace_back(mzs[i], scores[i]);
  }

  // filtration of top 500 masses is done

  int max_vertex_score = *std::max_element(_scores.begin(), _scores.end());
  int min_vertex_score = *std::min_element(_scores.begin(), _scores.end());

  max_path_score_ = std::max(max_vertex_score, max_vertex_score) * (max_tag_length_ + 2);
  min_path_score_ = std::max(min_vertex_score, min_vertex_score) * (max_tag_length_ + 2);

  max_path_score_ = std::max(max_path_score_, std::max(max_vertex_score, max_vertex_score) * (min_tag_length_ - 2));
  min_path_score_ = std::min(min_path_score_, std::max(min_vertex_score, min_vertex_score) * (min_tag_length_ - 2));
  min_path_score_ = std::max(0, min_path_score_);

  std::set<FLASHHelperClasses::Tag> tagSet;
  std::map<String, std::vector<FLASHHelperClasses::Tag>> seq_tag;

  for (int length = min_tag_length_; length <= max_tag_length_; length++)
  {
    FLASHHelperClasses::DAG dag(_mzs.size() * (1 + max_tag_length_) * (1 + max_gap_count_) * (1 + max_iso_in_tag_)
                                * (1 + max_path_score_ - min_path_score_));
    constructDAG_(dag, _mzs, _scores, length, ppm);

    std::set<FLASHHelperClasses::Tag> _tagSet;
    for (int score = max_path_score_; score >= min_path_score_ && (int)_tagSet.size() < max_tag_count_; score--)
    {
      std::vector<std::vector<Size>> all_paths;
      all_paths.reserve(max_tag_count_);
      for (int d = 0; d <= max_gap_count_; d++)
      {
        for (int iso = 0; iso <= max_iso_in_tag_; iso++)
        {
          dag.findAllPaths(getVertex_((int)_mzs.size() - 1, score, length, iso, d), getVertex_(0, 0, 0, 0, 0), all_paths, max_tag_count_);
        }
      }
      for (const auto& path : all_paths)
      {
        updateTagSet_(_tagSet, seq_tag, path, _mzs, _scores, scan, ppm);
      }
    }
    tagSet.insert(_tagSet.begin(), _tagSet.end());
  }

  for (int length = min_tag_length_; length <= max_tag_length_; length++)
  {
    // int count = 0;
    for (const auto& tag : tagSet)
    {
      if ((int)tag.getLength() != length) continue;
      tags_.push_back(tag);
      // if (++count == max_tag_count_) break;
    }
  }

  std::sort(tags_.begin(), tags_.end(),
            [](const FLASHHelperClasses::Tag& a, const FLASHHelperClasses::Tag& b) { return a.getScore() > b.getScore(); });
}

Size FLASHTaggerAlgorithm::find_with_X_(const std::string_view& A, const String& B, Size pos) // allow a single X. pos is in A
{
  for (size_t i = pos; i <= A.length() - B.length(); ++i)
  {
    bool match = true;
    int x_cntr = 0;
    for (size_t j = 0; j < B.length(); ++j)
    {
      if (A[i + j] == 'X') x_cntr++;
      if ((A[i + j] != B[j] && A[i + j] != 'X') || x_cntr > 1)
      {
        match = false;
        break;
      }
    }
    if (match) { return i; }
  }
  return String::npos;
}

// Make output struct containing all information about matched entries and tags, coverage, score etc.
void FLASHTaggerAlgorithm::runMatching(const std::vector<FASTAFile::FASTAEntry>& fasta_entry, double max_mod_mass, int tag_length)
{
  std::vector<std::pair<ProteinHit, std::vector<int>>> pairs;
  std::vector<int> start_loc(tags_.size(), 0);
  std::vector<int> end_loc(tags_.size(), 0);
  int scan = 0;
  protein_hits_.clear();
  // for each tag, find the possible start and end locations in the protein sequence. If C term, they are negative values to specify values are from
  // the end of the protein
  for (int i = 0; i < (int)tags_.size(); i++)
  {
    const auto& tag = tags_[i];
    if (tag_length > 0 && tag.getLength() != tag_length) continue;
    scan = tag.getScan();
    auto flanking_mass = std::max(tag.getNtermMass(), tag.getCtermMass());
    start_loc[i] = std::max(0, int(floor(flanking_mass - flanking_mass_tol_) / aa_mass_map_.rbegin()->first));
    end_loc[i] = int(ceil(flanking_mass + flanking_mass_tol_) / aa_mass_map_.begin()->first) + (int)tag.getLength() + 1;
  }

  // int min_hit_tag_score = max_path_score_;
  startProgress(0, (SignedSize)fasta_entry.size(), "Running FLASHTagger: searching database");

  double taget_count = 0;
  double decoy_count = 0;

  for (const auto& fe : fasta_entry)
  {
    if (fe.identifier.hasPrefix("DECOY")) { decoy_count++; }
    else { taget_count++; }
    decoy_factor_ = decoy_count / taget_count;
  }

#pragma omp parallel for default(none) shared(pairs, fasta_entry, taget_count, decoy_count, start_loc, end_loc, scan, tag_length,max_mod_mass, std::cout)
  for (int i = 0; i < (int)fasta_entry.size(); i++)
  {
    const auto& fe = fasta_entry[i];
    nextProgress();
    std::vector<int> matched_tag_indices;
    auto x_pos = fe.sequence.find('X');
    std::set<Size> matched_positions;

    std::map<int, std::map<int, int>> nterm_to_mzscore, cterm_to_mzscore;
    // std::map<Size, double> matched_protein_pos_mass;
    // find range, match allowing X.
    std::set<double> matched_masses;
    for (int j = 0; j < (int)tags_.size(); j++)
    {
      auto& tag = tags_[j];
      if (tag_length > 0 && tag.getLength() != tag_length) continue;
      auto A = tag.getMzs();
      bool isSubset = std::includes(matched_masses.begin(), matched_masses.end(), A.begin(), A.end());
      if (isSubset) { continue; }

      bool isNterm = tag.getNtermMass() > 0;

      int s, t;
      if (isNterm) { s = start_loc[j]; }
      else { s = std::max(0, int(fe.sequence.length()) - 1 - end_loc[j]); }
      t = std::min(end_loc[j] - start_loc[j], int(fe.sequence.length()) - s);
      if (t < (int)tag.getLength()) continue;
      const auto sub_seq = std::string_view(fe.sequence.data() + s, t);

      auto uppercase_tag_seq = tag.getSequence().toUpper();
      std::vector<int> positions;
      Size tpos = 0;
      while (true)
      {
        tpos = sub_seq.find(uppercase_tag_seq, tpos);
        if (tpos == std::string_view::npos) break;
        positions.push_back((int)tpos + s);
        tpos++;
      }

      if (positions.empty() && (int)x_pos >= s && (int)x_pos <= s + t) // only if perfect hits are not found and X exists
      {
        tpos = 0;
        while (true)
        {
          tpos = find_with_X_(sub_seq, uppercase_tag_seq, tpos);
          if (tpos == std::string_view::npos) break;
          positions.push_back((int)tpos + s);
          tpos++;
        }
      }

      bool matched = false;
      for (const auto& pos : positions)
      {
        if (tag.getNtermMass() > 0 && pos >= 0)
        {
          auto nterm = fe.sequence.substr(0, std::min(pos, (int)fe.sequence.length()));
          if (x_pos != String::npos) { nterm.erase(remove(nterm.begin(), nterm.end(), 'X'), nterm.end()); }
          double aamass = nterm.empty() ? 0 : AASequence::fromString(nterm).getMonoWeight(Residue::Internal);
          double flanking_mass = tag.getNtermMass() - aamass;
          if (std::abs(flanking_mass) > flanking_mass_tol_) continue;
          if (max_mod_mass > 0 && flanking_mass > max_mod_mass) continue;

          int nominal_flanking_mass = SpectralDeconvolution::getNominalMass(flanking_mass);
          if (nterm_to_mzscore.find(nominal_flanking_mass) == nterm_to_mzscore.end()) nterm_to_mzscore[nominal_flanking_mass] = std::map<int, int>();
          auto& sub_map = nterm_to_mzscore[nominal_flanking_mass];

          for (int off = 0; off <= (int)tag.getLength(); off++)
          {
            int mz = SpectralDeconvolution::getNominalMass(tag.getMzs()[off]);
            int score = tag.getScore(off);
            if (sub_map.find(mz) == sub_map.end()) sub_map[mz] = 0;
            sub_map[mz] = std::max(sub_map[mz], score);
          }
        }
        else if (tag.getCtermMass() > 0 && pos + tag.getSequence().length() < fe.sequence.length())
        {
          auto cterm = fe.sequence.substr(pos + tag.getSequence().length());
          if (x_pos != String::npos) cterm.erase(remove(cterm.begin(), cterm.end(), 'X'), cterm.end());
          double aamass = cterm.empty() ? 0 : AASequence::fromString(cterm).getMonoWeight(Residue::Internal);
          double flanking_mass = tag.getCtermMass() - aamass;
          if (std::abs(flanking_mass) > flanking_mass_tol_) continue;
          if (max_mod_mass > 0 && flanking_mass > max_mod_mass) continue;

          int nominal_flanking_mass = SpectralDeconvolution::getNominalMass(flanking_mass);
          if (cterm_to_mzscore.find(nominal_flanking_mass) == cterm_to_mzscore.end()) cterm_to_mzscore[nominal_flanking_mass] = std::map<int, int>();
          auto& sub_map = cterm_to_mzscore[nominal_flanking_mass];

          for (int off = 0; off <= (int)tag.getLength(); off++)
          {
            int mz = SpectralDeconvolution::getNominalMass(tag.getMzs()[off]);
            int score = tag.getScore(off);
            if (sub_map.find(mz) == sub_map.end()) sub_map[mz] = 0;
            sub_map[mz] = std::max(sub_map[mz], score);
          }
        }
        else
          continue;

        for (int off = 0; off <= (int)tag.getLength(); off++)
        {
          matched_positions.insert(pos + off);
        }
        matched = true;
      }
      if (matched)
      {
        matched_masses.insert(tag.getMzs().begin(), tag.getMzs().end());
        matched_tag_indices.push_back(j); // tag indices
      }
    }
    if (matched_tag_indices.empty()) continue;

    int match_cntr = 0;
    for (const auto& ps : matched_positions)
    {
      if (fe.sequence[ps] == 'X') continue;
      match_cntr++;
    }

    if (match_cntr < min_cov_aa_ + 1) continue;

    double match_score = 0;
    double max_match_score = 0;
    double total_max_match_score = 0;
    double prev_fm = -1;

    std::set<int> counted;
    for (const auto& [fm, mzscore] : nterm_to_mzscore)
    {
      if (fm - prev_fm > 0)
      {
        match_score = 0;
        counted.clear();
      }

      int score = 0;
      for (const auto& pair : mzscore)
      {
        if (counted.find(pair.first) != counted.end()) continue;
        counted.insert(pair.first);
        score += pair.second;
      }
      //std::cout<< fe.identifier << " " << std::to_string (fm) << " " << score << std::endl;

      match_score += score;
      max_match_score = std::max(max_match_score, match_score);
      prev_fm = fm;
    }
    total_max_match_score = max_match_score;
    max_match_score = 0;
    match_score = 0;
    prev_fm = -1;
    counted.clear();

    for (const auto& [fm, mzscore] : cterm_to_mzscore)
    {
      if (fm - prev_fm > 0) // always...
      {
        match_score = 0;
        counted.clear();
      }

      int score = 0;
      for (const auto& pair : mzscore)
      {
        if (counted.find(pair.first) != counted.end()) continue;
        counted.insert(pair.first);
        score += pair.second;
      }

      match_score += score;
      max_match_score = std::max(max_match_score, match_score);
      prev_fm = fm;
    }
    total_max_match_score = std::max(total_max_match_score, max_match_score);

    //std::cout<<fe.identifier << " score2:  " << total_max_match_score<<std::endl;

    ProteinHit hit(0, 0, fe.identifier, fe.sequence); //
    hit.setDescription(fe.description);
    hit.setMetaValue("Scan", scan);
    hit.setMetaValue("MatchedAA", match_cntr);
    hit.setCoverage(double(match_cntr) / (double)fe.sequence.length());
    hit.setScore(total_max_match_score);
#pragma omp critical
    {
      pairs.emplace_back(hit, matched_tag_indices);
    }
  }

  endProgress();

  protein_hits_.reserve(pairs.size());

  for (auto& [hit, indices] : pairs)
  {
    hit.setMetaValue("TagIndices", indices);
    protein_hits_.push_back(hit);
    // if (!hit.getAccession().hasPrefix("DECOY")) {
    //   std::cout<<hit.getAccession()<< " " << hit.getScore() << std::endl;
    //   for (auto i : indices) std::cout<< tags_[i].toString()<<std::endl;
    // }
  }
  std::sort(protein_hits_.begin(), protein_hits_.end(), [](const ProteinHit& left, const ProteinHit& right) {
    return left.getScore() == right.getScore()
             ? (left.getCoverage() == right.getCoverage() ? (left.getAccession() > right.getAccession()) : (left.getCoverage() > right.getCoverage()))
             : (left.getScore() > right.getScore());
  });
}

void FLASHTaggerAlgorithm::getProteinHits(std::vector<ProteinHit>& hits, int max_target_out) const
{
  hits.reserve(protein_hits_.size());
  int count = 0;
  double prev_score = 0;
  for (const auto& hit : protein_hits_)
  {
    if (hit.getScore() <= 0) break;
    hits.push_back(hit);
    if (hit.getAccession().hasPrefix("DECOY")) continue;
    if (max_target_out > 0 && count++ >= max_target_out)
    {
      if (hit.getScore() < prev_score) break; // keep adding if the score is the same
    }
    prev_score = hit.getScore();
  }
  // std::cout<<hits.size()<< " "  << prev_score << std::endl;
}

void FLASHTaggerAlgorithm::getTags(std::vector<FLASHHelperClasses::Tag>& tags, int tag_length) const
{
  for (const auto& tag : tags_)
  {
    if (tag_length > 0 && tag.getLength() != tag_length) continue;
    tags.push_back(tag);
  }
}

void FLASHTaggerAlgorithm::getMatchedPositionsAndFlankingMassDiffs(std::vector<int>& positions,
                                                                   std::vector<double>& masses,
                                                                   double flanking_mass_tol,
                                                                   const ProteinHit& hit,
                                                                   const FLASHHelperClasses::Tag& tag)
{
  Size pos = 0;
  std::vector<int> indices;
  const auto& seq = hit.getSequence();
  auto tagseq = tag.getSequence().toUpper();
  while (true)
  {
    pos = find_with_X_(seq, tagseq, pos + 1);
    if (pos == String::npos) break;
    double delta_mass = .0;

    auto x_pos = seq.find('X');
    if (tag.getNtermMass() > 0)
    {
      auto nterm = seq.substr(0, pos);
      if (x_pos != String::npos) { nterm.erase(remove(nterm.begin(), nterm.end(), 'X'), nterm.end()); }
      delta_mass = tag.getNtermMass() - (nterm.empty() ? 0 : AASequence::fromString(nterm).getMonoWeight(Residue::Internal));
    }
    else
    {
      auto cterm = seq.substr(pos + tag.getSequence().length());
      if (x_pos != String::npos) cterm.erase(remove(cterm.begin(), cterm.end(), 'X'), cterm.end());
      delta_mass = tag.getCtermMass() - (cterm.empty() ? 0 : AASequence::fromString(cterm).getMonoWeight(Residue::Internal));
    }
    if (std::abs(delta_mass) > flanking_mass_tol) continue;
    masses.push_back(delta_mass);
    positions.push_back((int)pos);
  }
}
} // namespace OpenMS