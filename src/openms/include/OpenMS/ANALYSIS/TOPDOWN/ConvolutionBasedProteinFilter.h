// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <boost/dynamic_bitset.hpp>
#include <iomanip>
#include <iostream>

namespace OpenMS
{
  /**
  @brief
  @ingroup Topdown
  */

  class OPENMS_DLLAPI ConvolutionBasedProteinFilter : public DefaultParamHandler, public ProgressLogger
  {
  public:
    /// constructor
    ConvolutionBasedProteinFilter();

    /// destructor
    ~ConvolutionBasedProteinFilter() override = default;

    /// copy constructor
    ConvolutionBasedProteinFilter(const ConvolutionBasedProteinFilter&);

    /// move constructor
    ConvolutionBasedProteinFilter(ConvolutionBasedProteinFilter&& other) = default;

    /// assignment operator
    ConvolutionBasedProteinFilter& operator=(const ConvolutionBasedProteinFilter& other);

    /// Find sequence tags from @p mzs and @p intensities then store them in @p tags.
    /**
      @brief
      Decoy or MS level 1 spectra are removed by this process.
      Overlapping PeakGroups in merged spectra are also removed.

      @param deconvolved_spectra spectra deconvolved by FLASHDeconv.
      @param ppm The acceptable ppm tolerance for mass
      @param fasta_entry fasta entry to searched against

    */

    void runMatching(const DeconvolvedSpectrum& deconvolved_spectrum, const std::vector<FASTAFile::FASTAEntry>& fasta_entry,
                     //const std::vector<std::vector<int>>& vectorized_fasta_entry_indices,
                     //const std::vector<std::vector<int>>& reversed_vectorized_fasta_entry_indices,
                     const std::vector<std::vector<Size>>& bit_protein_indices
                     //const std::vector<std::vector<Size>>& reversed_bit_protein_indices,
                     //int tag_length = 0
                     );
    const MSSpectrum& getSpectrum() const;
    void getProteinHits(std::vector<ProteinHit>& hits, int max_target_count) const;
    static void vectorizeFasta(const std::vector<FASTAFile::FASTAEntry>& fasta_entry,
                               std::vector<boost::dynamic_bitset<>>& vectorized_fasta_entry,
                               std::vector<std::vector<int>>& vectorized_fasta_entry_indices,
                               std::vector<std::map<int, double>>& mass_map,
                               std::vector<std::vector<Size>>& bit_protein_indices,
                               bool reverse);

    static const int multi_factor_for_vectorization = 1;
  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();

  private:
    void GetScoreAndMatchCount_(const std::vector<int>& spec_indices,
                                const std::vector<int>& pro_indices, std::vector<int>& spec_scores, int& max_score, int& match_cntr) const;

    MSSpectrum spec_;
    std::set<const Residue*> aas_ = ResidueDB::getInstance()->getResidues("Natural20");
    std::map<double, std::vector<Residue>> aa_mass_map_;
    std::vector<ProteinHit> protein_hits_;

    int min_tag_length_ = 0;
    double max_edge_mass_ = 0;
  };
} // namespace OpenMS