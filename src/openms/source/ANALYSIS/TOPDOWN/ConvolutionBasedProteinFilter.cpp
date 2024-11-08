// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/ConvolutionBasedProteinFilter.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <queue>
#include <vector>
#include <Eigen/Sparse>

#pragma warning(disable : 4100)

namespace OpenMS
{
ConvolutionBasedProteinFilter::ConvolutionBasedProteinFilter(): DefaultParamHandler("FLASHTaggerAlgorithm"), ProgressLogger()
{
  setDefaultParams_();
}

ConvolutionBasedProteinFilter::ConvolutionBasedProteinFilter(const ConvolutionBasedProteinFilter& other):
    DefaultParamHandler(other),
    ProgressLogger(other)
{
}

ConvolutionBasedProteinFilter& ConvolutionBasedProteinFilter::operator=(const ConvolutionBasedProteinFilter& rhs)
{
  if (this == &rhs) return *this;

  DefaultParamHandler::operator=(rhs);
  return *this;
}

void ConvolutionBasedProteinFilter::setDefaultParams_()
{
  defaults_.setValue(
    "min_length", 4,
    "Minimum length of a tag. Each mass gap contributes to a single length (even if a mass gap is represented by multiple amino acids). ");
  defaults_.setMaxInt("min_length", 9);
  defaults_.setMinInt("min_length", 2);

  defaultsToParam_();
}

void ConvolutionBasedProteinFilter::updateMembers_()
{
  min_tag_length_ = param_.getValue("min_length");
}

const MSSpectrum& ConvolutionBasedProteinFilter::getSpectrum() const
{
  return spec_;
}

void ConvolutionBasedProteinFilter::vectorizeFasta(const std::vector<FASTAFile::FASTAEntry>& fasta_entry,
                                                   std::vector<boost::dynamic_bitset<>>& vectorized_fasta_entry,
                                                   std::vector<std::vector<int>>& vectorized_fasta_entry_indices,
                                                   std::vector<std::map<int, double>>& mass_map,
                                                   std::vector<std::vector<Size>>& bit_protein_indices,
                                                   bool reverse)
{
  vectorized_fasta_entry.reserve(fasta_entry.size());
  vectorized_fasta_entry_indices.reserve(fasta_entry.size());
  mass_map.reserve(fasta_entry.size());
  int max_vec_size = 0;
  for (const auto& i : fasta_entry)
  {
    auto seq = i.sequence;
    seq.erase(remove(seq.begin(), seq.end(), 'X'), seq.end()); // remove all X
    // boost::dynamic_bitset<> vec(1 + round(multi_factor_for_vectorization * AASequence::fromString(seq).getMonoWeight(Residue::Internal)));
    boost::dynamic_bitset<> vec;
    std::vector<int> vec_indices;
    int pro_vec_len = 1 + round(multi_factor_for_vectorization * AASequence::fromString(seq).getMonoWeight(Residue::Internal));
    max_vec_size = std::max(max_vec_size, pro_vec_len);
    vec.resize(pro_vec_len);
    std::map<int, double> masses;
    double nmass = 0;
    vec[0] = true;
    vec_indices.push_back(0);
    masses[reverse ? seq.size() : 0] = .0;
    for (Size j = 0; j < seq.size(); j++)
    {
      Size index = reverse ? (seq.size() - j - 1) : j;
      nmass += AASequence::fromString(seq[index]).getMonoWeight(Residue::Internal);
      masses[index + (reverse ? 0 : 1)] = nmass;
      int vindex = int(round(multi_factor_for_vectorization * nmass));

      vec[vindex] = true;
      vec_indices.push_back(vindex);
      max_vec_size = std::max(max_vec_size, vindex);
      // vec[(l/2 + SpectralDeconvolution::getNominalMass(nmass)) % l] = true;
    }

    mass_map.push_back(masses);
    vectorized_fasta_entry.push_back(vec);
    vectorized_fasta_entry_indices.push_back(vec_indices);
  }

  bit_protein_indices.resize(max_vec_size);
  for (int fasta_index = 0; fasta_index < vectorized_fasta_entry_indices.size(); fasta_index++)
  {
    for (const auto vindex : vectorized_fasta_entry_indices[fasta_index])
    {
      bit_protein_indices[vindex].push_back(fasta_index);
    }
  }
}


void ConvolutionBasedProteinFilter::GetScoreAndMatchCount_(const std::vector<int>& spec_indices,
                                                           const std::vector<int>& pro_indices,
                                                           std::vector<int>& spec_scores,
                                                           int& max_score,
                                                           int& match_cntr) const
{
  // Step 1: Extract indices of set bits from spec_vec and pro_vec
  // Step 2: Prepare a score and match counter map (or array)

  int result_size = 1 + spec_indices.back() + pro_indices.back();
  Eigen::Matrix<short, Eigen::Dynamic, 1> scores = Eigen::Matrix<short, Eigen::Dynamic, 1>::Zero(result_size);
  Eigen::Matrix<short, Eigen::Dynamic, 1> matches = Eigen::Matrix<short, Eigen::Dynamic, 1>::Zero(result_size);

  // Step 3: Convolve the indices using SIMD to speed up score and match count updates
  for (size_t i = 0; i < pro_indices.size(); ++i)
  {
    for (size_t j = 0; j < spec_indices.size(); ++j)
    {
      size_t index_sum = pro_indices[i] + spec_indices[j];

      if (index_sum + 8 <= result_size) // Processing 8 elements at a time (for 16-bit integers)
      {
        // Use Eigen's block operations to vectorize the addition with 16-bit integers
        scores.segment<8>(index_sum) += Eigen::Matrix<short, 8, 1>::Constant(spec_scores[j]);
        matches.segment<8>(index_sum) += Eigen::Matrix<short, 8, 1>::Ones();
      }
      else
      {
        // For leftovers or small cases, handle the score and match count normally
        scores[index_sum] += spec_scores[j];
        matches[index_sum]++;
      }
    }
  }

  // Step 4: Find the maximum score and its corresponding match count
  max_score = 0;
  match_cntr = 0;

  // Simple scalar loop for finding the max score and match count
  for (size_t i = 0; i < scores.size(); ++i)
  {
    if (scores[i] > max_score)
    {
      max_score = scores[i];
      match_cntr = matches[i];
    }
  }

}

// Make output struct containing all information about matched entries and tags, coverage, score etc.
void ConvolutionBasedProteinFilter::runMatching(const DeconvolvedSpectrum& deconvolved_spectrum,
                                                const std::vector<FASTAFile::FASTAEntry>& fasta_entry,
                                               // const std::vector<std::vector<int>>& vectorized_fasta_entry_indices,
                                                //const std::vector<std::vector<int>>& reversed_vectorized_fasta_entry_indices,
                                                const std::vector<std::vector<Size>>& bit_protein_indices
                                               // const std::vector<std::vector<Size>>& reversed_bit_protein_indices,
                                               // int tag_length
                                                )
{
  int scan = deconvolved_spectrum.getScanNumber();
  protein_hits_.clear();

  std::vector<int> spec_indices;
  std::vector<int> spec_scores;

  spec_indices.push_back(0);
  spec_scores.push_back(1);

  auto t_dspec(deconvolved_spectrum);
  t_dspec.sortByQscore();

  for (const auto& pg : t_dspec)
  {
    int mn = multi_factor_for_vectorization > 1 ? round(ConvolutionBasedProteinFilter::multi_factor_for_vectorization * pg.getMonoMass())
      : SpectralDeconvolution::getNominalMass((pg.getMonoMass()));
    spec_indices.push_back(mn);
    spec_scores.push_back(FLASHTaggerAlgorithm::getPeakGroupScore(pg));
    if (spec_indices.size() >= 16) break;
  }
  std::sort(spec_indices.begin(), spec_indices.end());
  startProgress(0, (SignedSize)fasta_entry.size(), "Running Protein filter: searching database");

  int result_size = 1 + spec_indices.back() + bit_protein_indices.size();
  std::vector<short> scores(result_size * (fasta_entry.size() + 1), 0);
  //auto matches = std::map<std::tuple<int, int>, int>();
  std::vector<short> max_scores(fasta_entry.size(), 0);

for (size_t i = 0; i < bit_protein_indices.size(); ++i)
  {
    if (bit_protein_indices[i].empty()) continue;
    const auto& bis = bit_protein_indices[i];
    for (size_t j = 0; j < spec_indices.size(); ++j)
    {
      size_t index_sum = i + spec_indices[j];
      for (auto pi : bis)
      {
        // Ensure the indices are within the dimensions of the matrix
        Size index = index_sum * max_scores.size() + pi;

        scores[index] += spec_scores[j];
        max_scores[pi] = std::max(max_scores[pi], scores[index]);
      }
    }
  }
//  for (size_t i = 0; i < reversed_bit_protein_indices.size(); ++i)
//  {
//    for (size_t pi : reversed_bit_protein_indices[i])
//    {
//      for (size_t j = 0; j < spec_indices.size(); ++j)
//      {
//        size_t index_sum = i + spec_indices[j];
//
//        scores(index_sum, pi) += spec_scores[j];
//        matches(index_sum, pi)++;
//
//      }
//    }
//  }

  std::map<int, int> index_to_score;

//#pragma omp parallel for default(none) shared(fasta_entry, vectorized_fasta_entry_indices, reversed_vectorized_fasta_entry_indices, spec_indices, \
//                                                spec_scores, tag_length, scores, scan, index_to_score, result_size, std::cout)
  for (size_t i = 0; i < fasta_entry.size(); i++)
  {

    index_to_score[i] = max_scores[i];

  }

  std::priority_queue<std::pair<int, int>, std::vector<std::pair<int, int>>, std::greater<>> min_heap;

  // Iterate over the map to populate the min-heap
  for (const auto& [index, score] : index_to_score) {
    min_heap.emplace(score, index);  // Add (score, index) to the heap

    // If heap size exceeds N, remove the smallest element
    if (min_heap.size() > 30) {
      min_heap.pop();
    }
  }

  // Extract the indices from the heap into a vector
  std::vector<int> top_indices;
  while (!min_heap.empty()) {
    top_indices.push_back(min_heap.top().second);  // Get the index
    min_heap.pop();
  }

  // Optional: Sort the indices by score in descending order
  std::reverse(top_indices.begin(), top_indices.end());

  for (const auto i : top_indices)
  {
    const auto& fe = fasta_entry[i];
    std::cout<<fe.identifier << " " << index_to_score[i]<<std::endl;
    int match_cntr = 1;
    ProteinHit hit(0, 0, fe.identifier, fe.sequence); //
    hit.setDescription(fe.description);
    hit.setMetaValue("Scan", scan);
    hit.setMetaValue("MatchedAA", match_cntr);
    hit.setCoverage(double(match_cntr) / (double)fe.sequence.length()); // TODO change.. when truncation happens
    hit.setScore(index_to_score[i]);
    protein_hits_.push_back(hit);
  }
  //
  // #pragma omp parallel for default(none) \
//  shared(fasta_entry, vectorized_fasta_entry_indices, reversed_vectorized_fasta_entry_indices, spec_indices, spec_scores, tag_length, scan,
  //  std::cout) for (int i = 0; i < (int)fasta_entry.size(); i++)
  //  {
  //    const auto& fe = fasta_entry[i];
  //    const auto& pro_vec = vectorized_fasta_entry_indices[i];
  //    const auto& rev_pro_vec = reversed_vectorized_fasta_entry_indices[i];
  //    nextProgress();
  //    int n_max_score, n_match_cntr, c_max_score, c_match_cntr;
  //    GetScoreAndMatchCount_(spec_indices, pro_vec, spec_scores, n_max_score, n_match_cntr);
  //    GetScoreAndMatchCount_(spec_indices, rev_pro_vec, spec_scores, c_max_score, c_match_cntr);
  //
  //    int match_cntr = n_match_cntr + c_match_cntr; // std::max(n_match_cntr, c_match_cntr);
  //    int max_score = n_max_score + c_max_score;    // std::max(n_max_score, c_max_score);
  //
  //    if (match_cntr < 6) continue;
  //
  //    ProteinHit hit(0, 0, fe.identifier, fe.sequence); //
  //    hit.setDescription(fe.description);
  //    hit.setMetaValue("Scan", scan);
  //    hit.setMetaValue("MatchedAA", match_cntr);
  //    hit.setCoverage(double(match_cntr) / (double)fe.sequence.length()); // TODO change.. when truncation happens
  //    hit.setScore(max_score);
  // #pragma omp critical
  //    {
  //      protein_hits_.push_back(hit);
  //      // pairs.emplace_back(hit, matched_tag_indices);
  //    }
  //  }

  endProgress();
  if (! protein_hits_.empty())
  {
    std::sort(protein_hits_.begin(), protein_hits_.end(), [](const ProteinHit& left, const ProteinHit& right) {
      return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getAccession() > right.getAccession())
                                                                                              : (left.getCoverage() > right.getCoverage()))
                                                 : (left.getScore() > right.getScore());
    });
  }
}

void ConvolutionBasedProteinFilter::getProteinHits(std::vector<ProteinHit>& hits, int max_target_out) const
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
} // namespace OpenMS