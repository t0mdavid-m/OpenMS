// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace OpenMS
{
/**
@brief Gapped tags: protein candidates and anchors for spectra in which no sequence tag is found.

A gapped tag is a set of deconvolved masses whose differences are sums of one or more residues. It lands at a residue position of a
protein when windows starting there have masses matching every offset of the tag from its lowest mass (read as a prefix series) or
from its highest mass (read as a suffix series). A landing is returned as protein, first and one-past-last residue, and the tag
masses, which is what the extension needs to anchor a path (FLASHHelperClasses::Tag with an explicit length).

Port of GappedFiltration (T. D. Mueller, https://github.com/t0mdavid-m/GappedFiltration) following its C++ form GappedFiltrationPP
(O. Kohlbacher, https://github.com/okohlbacher/GappedFiltrationPP): every window of 1 to max_gap_mass / 57 residues is binned by mass
into sorted posting lists of start positions, and a bit set with a rank structure marks the masses that are sums of residues.
Differences: the tolerance of an offset is the deconvolution tolerance of each of its two masses (not a fraction of the gap), and
every mass bin the tolerance window touches is looked up; the search is complete rather than a beam over tags: with every mass as
the base (the lowest of a prefix tag, the highest of a suffix tag), each other mass within max_gap_mass votes, through the index, for
the window starts matching its offset, so every start reached by at least min_tag_size - 1 offsets is found, with all the masses that
reach it. Each pair of masses is looked up once, where a beam looks up each tag it keeps. In a crowded spectrum such tags land on most
proteins of a database by chance, so the proteins are ranked by the score FLASHTnT ranks candidates by (FLASHTaggerAlgorithm::
scoreShifts at the landings' shifts) and the landings of the max_proteins best are reported for the extension.
*/
class OPENMS_DLLAPI FLASHGappedTaggerAlgorithm : public DefaultParamHandler
{
public:
  /// a gapped tag placed on a protein
  struct Landing
  {
    Size protein = 0;           ///< index of the protein in the sequences given to buildIndex
    int start = 0;              ///< first residue of the tag in the protein (0-based)
    int end = 0;                ///< one past its last residue
    bool prefix = true;         ///< masses read as a prefix (N-terminal) series; otherwise as a suffix series
    double flank_mass = 0;      ///< residue mass before start (prefix) or after end (suffix)
    std::vector<double> masses; ///< tag masses, ascending
    std::vector<int> scores;    ///< node score of each mass
    double qscore_sum = 0;      ///< summed Qscore of the tag masses
  };

  FLASHGappedTaggerAlgorithm();
  FLASHGappedTaggerAlgorithm(const FLASHGappedTaggerAlgorithm&) = default;
  FLASHGappedTaggerAlgorithm& operator=(const FLASHGappedTaggerAlgorithm&) = default;
  ~FLASHGappedTaggerAlgorithm() override = default;

  /// index the protein sequences; residues other than the 20 standard ones end every window containing them
  void buildIndex(const std::vector<std::string>& protein_seqs);

  /**
    @brief every gapped tag of at least min_tag_size masses (of the top_n) and where it lands, for the max_proteins best proteins. Thread
    safe once the index is built.
    @param dspec deconvolved spectrum
    @param ppm mass tolerance of a deconvolved mass in ppm
    @param landings output, one per protein position and series (the tag with the most masses, then the highest Qscore sum)
  */
  void search(const DeconvolvedSpectrum& dspec, double ppm, std::vector<Landing>& landings) const;

  /// bytes held by the index and the plausibility set
  Size getIndexBytes() const;

protected:
  void updateMembers_() override;

private:
  /// true if a sum of residues lies in [lo, hi] Da
  bool plausible_(double lo, double hi) const;
  std::int64_t key_(double mass) const
  {
    return (std::int64_t)std::llround(mass / bin_width_);
  }

  // parameters
  double max_gap_mass_ = 5000;
  double bin_width_ = 0.1;
  int min_tag_size_ = 5;
  int top_n_ = 400;
  int max_proteins_ = 100;

  // proteins: residues concatenated, prefix sums of residue masses (0 for a non-standard residue), protein start offsets
  std::vector<double> cum_mass_;
  std::vector<std::uint8_t> valid_;
  std::vector<Size> offsets_;

  // mass index: bin -> slot of a non-empty bin; per slot a sorted posting list of window start positions
  std::int64_t key_min_ = 0, key_max_ = 0;
  std::vector<std::uint32_t> slot_;
  std::vector<std::uint64_t> post_offsets_;
  std::vector<std::uint32_t> postings_;

  // residue sums on a 0.1 mDa grid up to max_gap_mass, and the number of set bits before each 64-bit word
  std::vector<std::uint64_t> plausible_bits_;
  std::vector<std::uint32_t> plausible_rank_;
  std::int64_t plausible_units_ = 0;
  static constexpr double plausible_scale_ = 1e4;
};
} // namespace OpenMS
