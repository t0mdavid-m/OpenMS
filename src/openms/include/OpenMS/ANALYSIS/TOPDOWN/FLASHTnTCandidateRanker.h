// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <array>
#include <vector>

namespace OpenMS
{
  /**
   * @brief The learned candidate cut of FLASHTnT (candidate_ranking "learned"): each candidate protein of a spectrum gets a score from a
   * small model of the spectrum's own evidence for it, and the max_protein_hits best are extended.
   *
   * A candidate's features read the spectrum's deconvolved masses, its tags' shifts (the protein's residue sum before or after a tag minus
   * the tag's flanking mass) and its protein's residue sums; not the precursor mass, not other candidates, not other spectra:
   *   0 score: what the candidate search scored it;
   *   1 best_chance: over the sides with tags, -log10 P(Poisson >= the masses on that side's ion ladder) at the side's offset (the first
   *     0.05 Da window holding the most of its shifts, their median), lambda from the ladder's density at 10 ppm;
   *   2 n_bound_err: the distance (Da, at most 9.99, 9.99 without N-side tags) of the cut before the start (N offset) from a residue
   *     boundary of the protein, allowing a common N-terminal modification (acetyl, formyl, methyl, dimethyl, trimethyl, myristoyl,
   *     palmitoyl, pyro-Glu) and +-1 isotope;
   *   3 c_bound_err: the same for the cut after the end (C offset), allowing amidation;
   *   4 c_intact: the C offset is the intact C-terminus (within 1.1 Da of 0).
   * The model (FLASHTnTCandidateModel.h) is a LightGBM LambdaRank, trained on 11 top-down runs to put the protein of each spectrum's
   * best-scoring PrSM of a wide search (500 candidates) first, wherever that PrSM came from a tag candidate. Left out by dataset, its 10
   * best held that protein in 94.6 % of the spectra (gapped-tag stage 89.2 %), the union of the 100 best by score and by matched tags in
   * 93.8 % (86.8 %).
   */
  class OPENMS_DLLAPI FLASHTnTCandidateRanker
  {
  public:
    /// the spectrum: its deconvolved masses (ascending) and the ion type offsets from the residue sums (prefix: c ions NH3; suffix: z
    /// ions -NH3, z+1 ions -16.018724)
    struct Spectrum
    {
      std::vector<double> masses;
      double prefix_offset = 0, suffix_offset = 0;
    };

    /// a candidate: the residue sums of its protein from each terminus (ascending from 0, fixed modifications included), the shifts of its
    /// tags on each side and what the candidate search scored
    struct Candidate
    {
      const std::vector<double>* prefix = nullptr;
      const std::vector<double>* suffix = nullptr;
      const double* n_shifts = nullptr;
      Size n_count = 0;
      const double* c_shifts = nullptr;
      Size c_count = 0;
      double score = 0;
    };

    static constexpr Size feature_count = 5;

    /// the features of a candidate
    static std::array<double, feature_count> features(const Spectrum& spectrum, const Candidate& candidate);

    /// the model's score of a candidate's features: the higher, the more likely the protein of the spectrum's best PrSM
    static double score(const std::array<double, feature_count>& features);
  };
} // namespace OpenMS
