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
#include <Eigen/Dense>
#include <vector>

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
                                                   bool reverse)
{
  vectorized_fasta_entry.reserve(fasta_entry.size());
  vectorized_fasta_entry_indices.reserve(fasta_entry.size());
  mass_map.reserve(fasta_entry.size());
  for (const auto& i : fasta_entry)
  {
    auto seq = i.sequence;
    seq.erase(remove(seq.begin(), seq.end(), 'X'), seq.end()); // remove all X
    // boost::dynamic_bitset<> vec(1 + round(multi_factor_for_vectorization * AASequence::fromString(seq).getMonoWeight(Residue::Internal)));
    boost::dynamic_bitset<> vec;
    std::vector<int> vec_indices;
    int pro_vec_len = 1 + round(multi_factor_for_vectorization * AASequence::fromString(seq).getMonoWeight(Residue::Internal));

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
      // vec[(l/2 + SpectralDeconvolution::getNominalMass(nmass)) % l] = true;
    }

    mass_map.push_back(masses);
    vectorized_fasta_entry.push_back(vec);
    vectorized_fasta_entry_indices.push_back(vec_indices);
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

      if (index_sum + 8 <= result_size)  // Processing 8 elements at a time (for 16-bit integers)
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

  // int result_size = 1 + spec_indices.back() + pro_indices.back();
  // std::vector<short> scores(result_size, 0);
  // std::vector<short> matches(result_size, 0);

  // // Step 3: Convolve the indices using SIMD to speed up score and match count updates
  // for (size_t i = 0; i < pro_indices.size(); ++i)
  // {
  //   for (size_t j = 0; j < spec_indices.size(); ++j)
  //   {
  //     size_t index_sum = pro_indices[i] + spec_indices[j];

  //     if (index_sum + 8 <= result_size)  // Processing 8 elements at a time (for 16-bit integers)
  //     {
  //       // Use SIMD to process 8 elements at a time
  //       int16x8_t score_val = vdupq_n_s16(spec_scores[j]);  // Broadcast the score value

  //       // Load current scores and matches at the target index_sum
  //       int16x8_t current_scores = vld1q_s16(&scores[index_sum]);
  //       int16x8_t current_matches = vld1q_s16(&matches[index_sum]);

  //       // Increment the scores and matches using SIMD
  //       int16x8_t new_scores = vaddq_s16(current_scores, score_val);
  //       int16x8_t new_matches = vaddq_s16(current_matches, vdupq_n_s16(1));

  //       // Store the updated scores and matches
  //       vst1q_s16(&scores[index_sum], new_scores);
  //       vst1q_s16(&matches[index_sum], new_matches);
  //     }
  //     else
  //     {
  //       // For leftovers or small cases, handle the score and match count normally
  //       scores[index_sum] += spec_scores[j];
  //       matches[index_sum]++;
  //     }
  //   }
  // }

  // Step 4: Find the maximum score and its corresponding match count
  max_score = 0;
  match_cntr = 0;

  // Simple scalar loop for finding the max score and match count
  for (size_t i = 0; i < scores.size(); ++i) {
    if (scores[i] > max_score) {
      max_score = scores[i];
      match_cntr = matches[i];
    }
  }

  //  std::vector<int> scores(spec_vec.size() + pro_vec.size(), 0); // for SIMD-friendly access
  //  std::vector<int> matches(spec_vec.size() + pro_vec.size(), 0); // same for matches
  //
  //  Size pro_vec_index = pro_vec.find_first();
  //  while (pro_vec_index != pro_vec.npos)
  //  {
  //    Size spec_vec_index = spec_vec.find_first();
  //    while (spec_vec_index != spec_vec.npos)
  //    {
  //      Size index_sum = spec_vec_index + pro_vec_index;
  //
  //      // NEON parallel addition for spec_scores
  //      // Assuming both spec_scores and scores are aligned and contiguous
  //      int32x4_t score_val = vdupq_n_s32(spec_scores[spec_vec_index]);  // Broadcast score value
  //      int32x4_t current_score = vld1q_s32(&scores[index_sum]);         // Load current scores
  //      int32x4_t new_score = vaddq_s32(current_score, score_val);       // Vector add
  //      vst1q_s32(&scores[index_sum], new_score);                        // Store result
  //
  //      matches[index_sum]++;  // Simple increment, can be NEON-optimized similarly
  //
  //      spec_vec_index = spec_vec.find_next(spec_vec_index);
  //    }
  //    pro_vec_index = pro_vec.find_next(pro_vec_index);
  //  }
  //
  //  max_score = 0;
  //  match_cntr = 0;
  //  for (Size i = 0; i < scores.size(); ++i)
  //  {
  //    if (scores[i] > max_score)
  //    {
  //      max_score = scores[i];
  //      match_cntr = matches[i];
  //    }
  //  }


  //  std::map<Size, int> scores;
  //  std::map<Size, int> matches;
  //
  //  Size pro_vec_index = pro_vec.find_first();
  //  while (pro_vec_index != pro_vec.npos)
  //  {
  //    Size spec_vec_index = spec_vec.find_first();
  //    while (spec_vec_index != spec_vec.npos)
  //    {
  //      if (scores.find(spec_vec_index + pro_vec_index) == scores.end())
  //      {
  //        scores[spec_vec_index + pro_vec_index] = 0;
  //        matches[spec_vec_index + pro_vec_index] = 0;
  //      }
  ////      scores[spec_vec_index + pro_vec_index] += spec_scores[spec_vec_index];
  ////      matches[spec_vec_index + pro_vec_index] ++;
  //      spec_vec_index = spec_vec.find_next(spec_vec_index);
  //    }
  //    pro_vec_index = pro_vec.find_next(pro_vec_index);
  //  }
  //  max_score = 0;
  //  match_cntr = 0;
  //  for (const auto& [i, s] : scores)
  //  {
  //    if (max_score > s) continue;
  //    max_score = s;
  //    match_cntr = matches[i];
  //  }
}

// Make output struct containing all information about matched entries and tags, coverage, score etc.
void ConvolutionBasedProteinFilter::runMatching(const DeconvolvedSpectrum& deconvolved_spectrum,
                                                const std::vector<FASTAFile::FASTAEntry>& fasta_entry,
                                                const std::vector<std::vector<int>>& vectorized_fasta_entry_indices,
                                                const std::vector<std::vector<int>>& reversed_vectorized_fasta_entry_indices,
                                                int tag_length)
{
  int scan = deconvolved_spectrum.getScanNumber();
  protein_hits_.clear();

  std::vector<int> spec_indices;
  std::vector<int> spec_scores;

  spec_indices.push_back(0);
  spec_scores.push_back(1);
  for (const auto& pg : deconvolved_spectrum)
  {
    int mn = SpectralDeconvolution::getNominalMass(pg.getMonoMass()); // SpectralDeconvolution::getNominalMass(deconvolved_spectrum[deconvolved_spectrum.size()
                                                                      // - 1].getMonoMass() - pg.getMonoMass());
    spec_indices.push_back(mn);
    spec_scores.push_back(FLASHTaggerAlgorithm::getPeakGroupScore(pg));
  }

  startProgress(0, (SignedSize)fasta_entry.size(), "Running Protein filter: searching database");

#pragma omp parallel for default(none) \
  shared(fasta_entry, vectorized_fasta_entry_indices, reversed_vectorized_fasta_entry_indices, spec_indices, spec_scores, tag_length, scan, std::cout)
  for (int i = 0; i < (int)fasta_entry.size(); i++)
  {
    const auto& fe = fasta_entry[i];
    const auto& pro_vec = vectorized_fasta_entry_indices[i];
    const auto& rev_pro_vec = reversed_vectorized_fasta_entry_indices[i];
    nextProgress();
    int n_max_score, n_match_cntr, c_max_score, c_match_cntr;
    GetScoreAndMatchCount_(spec_indices, pro_vec, spec_scores, n_max_score, n_match_cntr);
    GetScoreAndMatchCount_(spec_indices, rev_pro_vec, spec_scores, c_max_score, c_match_cntr);

    int match_cntr = n_match_cntr + c_match_cntr; // std::max(n_match_cntr, c_match_cntr);
    int max_score = n_max_score + c_max_score;    // std::max(n_max_score, c_max_score);

    if (match_cntr < 6) continue;

    ProteinHit hit(0, 0, fe.identifier, fe.sequence); //
    hit.setDescription(fe.description);
    hit.setMetaValue("Scan", scan);
    hit.setMetaValue("MatchedAA", match_cntr);
    hit.setCoverage(double(match_cntr) / (double)fe.sequence.length()); // TODO change.. when truncation happens
    hit.setScore(max_score);
#pragma omp critical
    {
      protein_hits_.push_back(hit);
      // pairs.emplace_back(hit, matched_tag_indices);
    }
  }

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