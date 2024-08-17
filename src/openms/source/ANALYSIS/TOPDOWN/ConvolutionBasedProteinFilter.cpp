// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/ConvolutionBasedProteinFilter.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>

namespace OpenMS
{
ConvolutionBasedProteinFilter::ConvolutionBasedProteinFilter(): DefaultParamHandler("FLASHTaggerAlgorithm"), ProgressLogger()
{
  setDefaultParams_();
}

ConvolutionBasedProteinFilter::ConvolutionBasedProteinFilter(const ConvolutionBasedProteinFilter& other): DefaultParamHandler(other), ProgressLogger(other)
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
                                                   std::vector<std::map<int, double>>& mass_map,
                                                   bool reverse)
{
  vectorized_fasta_entry.reserve(fasta_entry.size());
  mass_map.reserve(fasta_entry.size());
  for (const auto & i : fasta_entry)
  {
    auto seq = i.sequence;
    seq.erase(remove(seq.begin(), seq.end(), 'X'), seq.end()); // remove all X
    boost::dynamic_bitset<> vec(1 + round(multi_factor_for_vectorization * AASequence::fromString(seq).getMonoWeight(Residue::Internal)));
    //int l = 4000;
    //boost::dynamic_bitset<> vec(l);

    std::map<int, double> masses;
    double nmass = 0;
    vec[0] = true;
    masses[reverse? seq.size() : 0] = .0;
    for (Size j = 0; j < seq.size(); j++)
    {
      Size index = reverse? (seq.size() - j - 1): j;
      nmass += AASequence::fromString(seq[index]).getMonoWeight(Residue::Internal);
      masses[index + (reverse? 0 : 1)] = nmass;
      vec[round(multi_factor_for_vectorization * nmass)] = true;
      //vec[(l/2 + SpectralDeconvolution::getNominalMass(nmass)) % l] = true;
    }

    mass_map.push_back(masses);
    vectorized_fasta_entry.push_back(vec);
  }
}

void ConvolutionBasedProteinFilter::GetScoreAndMatchCount_(const boost::dynamic_bitset<>& spec_vec, const boost::dynamic_bitset<>& pro_vec, std::vector<int>& spec_scores, int& max_score, int& match_cntr) const
{
  std::vector<int> scores(spec_vec.size() + pro_vec.size(), 0);
  std::vector<int> matches(spec_vec.size() + pro_vec.size(), 0);

  Size pro_vec_index = pro_vec.find_first();
  while (pro_vec_index != pro_vec.npos)
  {
    Size spec_vec_index = spec_vec.find_first();
    while (spec_vec_index != spec_vec.npos)
    {
      scores[spec_vec_index + pro_vec_index] += spec_scores[spec_vec_index];
      matches[spec_vec_index + pro_vec_index] ++;
      spec_vec_index = spec_vec.find_next(spec_vec_index);
    }
    pro_vec_index = pro_vec.find_next(pro_vec_index);
  }
  max_score = 0;
  match_cntr = 0;
  for (Size i = 0; i < scores.size(); i++)
  {
    if (max_score > scores[i]) continue;
    max_score = scores[i];
    match_cntr = matches[i];
  }
}

// Make output struct containing all information about matched entries and tags, coverage, score etc.
void ConvolutionBasedProteinFilter::runMatching(const DeconvolvedSpectrum& deconvolved_spectrum, const std::vector<FASTAFile::FASTAEntry>& fasta_entry, const std::vector<boost::dynamic_bitset<>>& vectorized_fasta_entry, const std::vector<boost::dynamic_bitset<>>& reversed_vectorized_fasta_entry, double max_mod_mass, int tag_length)
{
  int scan = deconvolved_spectrum.getScanNumber();
  protein_hits_.clear();

  boost::dynamic_bitset<> spec_vec(1 + SpectralDeconvolution::getNominalMass(deconvolved_spectrum[deconvolved_spectrum.size() - 1].getMonoMass()));
  std::vector<int> spec_scores(spec_vec.size(), 0);

  spec_vec[0] = true;
  spec_scores[0] = 1;
  for (const auto& pg : deconvolved_spectrum)
  {
    int mn = SpectralDeconvolution::getNominalMass(pg.getMonoMass());//SpectralDeconvolution::getNominalMass(deconvolved_spectrum[deconvolved_spectrum.size() - 1].getMonoMass() - pg.getMonoMass());
    spec_vec[mn] = true;
    spec_scores[mn] = FLASHTaggerAlgorithm::getPeakGroupScore(pg);
  }

  startProgress(0, (SignedSize)fasta_entry.size(), "Running Protein filter: searching database");

#pragma omp parallel for default(none) shared(fasta_entry, vectorized_fasta_entry, reversed_vectorized_fasta_entry, spec_vec,spec_scores, tag_length, scan,  std::cout)
  for (int i = 0; i < (int)fasta_entry.size(); i++)
  {
    const auto& fe = fasta_entry[i];
    const auto& pro_vec = vectorized_fasta_entry[i];
    const auto& rev_pro_vec = reversed_vectorized_fasta_entry[i];
    nextProgress();
    int n_max_score, n_match_cntr, c_max_score, c_match_cntr;
    GetScoreAndMatchCount_(spec_vec, pro_vec, spec_scores, n_max_score, n_match_cntr);
    GetScoreAndMatchCount_(spec_vec, rev_pro_vec, spec_scores, c_max_score, c_match_cntr);

    int match_cntr = n_match_cntr + c_match_cntr;//std::max(n_match_cntr, c_match_cntr);
    int max_score = n_max_score + c_max_score;//std::max(n_max_score, c_max_score);

    if (match_cntr < 6) continue;

    ProteinHit hit(0, 0, fe.identifier, fe.sequence); //
    hit.setDescription(fe.description);
    hit.setMetaValue("Scan", scan);
    hit.setMetaValue("MatchedAA", match_cntr);
    hit.setCoverage(double(match_cntr) / (double)fe.sequence.length());
    hit.setScore(max_score);
#pragma omp critical
    {
      protein_hits_.push_back(hit);
      //pairs.emplace_back(hit, matched_tag_indices);
    }
  }

  endProgress();
  if (!protein_hits_.empty())
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