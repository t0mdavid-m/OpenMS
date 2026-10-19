// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/config.h>
#include <iomanip>
#include <iostream>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>

namespace OpenMS
{
  class OPENMS_DLLAPI FLASHTnTFile
{
  public:
    /// write header line for tag output file
    static void writeTagHeader(std::fstream& fs);

    /// write header line for PrSM output file
    static void writePrSMHeader(std::fstream& fs);

    /// write header line for Proteoform output file
    static void writeProHeader(std::fstream& fs);

    /// header for the protein-level output, written before proteoform characterization
    static void writeProteinHeader(std::fstream& fs);

    /// the candidate proteins and their protein-level scores
    static void writeProteins(const std::vector<ProteinHit>& hits, std::fstream& fs);

    /// write the features in regular file output
    static void writeTags(const FLASHTnTAlgorithm& tnt, double flanking_mass_tol, std::fstream& fs);

    static void writePrSMs(const std::vector<ProteinHit>& hits, std::fstream& fs);

    /**
      @brief the ProForma 2.0 string of a PrSM or proteoform hit; empty if a terminus is undetermined (ProForma cannot express it)

      A variable (searched) modification is written by its name (two names joined by '+' as two tags), a blind shift by its mass: the
      Unimod names listed for a blind shift (ModID) are only those of matching mass. Overlapping ranges, which ProForma does not
      allow, are written as one range with their summed mass.
    */
    static String getProFormaString(const ProteinHit& hit);

    static void writeProteoforms(const std::vector<ProteinHit>& hits, std::fstream& fs, double pro_fdr);
  };
} // namespace OpenMS