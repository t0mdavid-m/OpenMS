// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace OpenMS
{

  std::vector<std::string> MS3FragmentMatcher::getMS3IonTypes(char precursor_ion_class)
  {
    switch (precursor_ion_class)
    {
      case 'b':
      case 'a':
      case 'c':
        return {"a", "b", "yb", "ya"};
      case 'y':
      case 'x':
      case 'z':
        return {"a", "b", "y"};
      default:
        return {"a", "b", "y"};
    }
  }

  bool MS3FragmentMatcher::isPrefixIonType(const std::string& ion_type)
  {
    return ion_type == "a" || ion_type == "b" || ion_type == "c";
  }

  double MS3FragmentMatcher::getIonShift(const std::string& ion_type)
  {
    if (ion_type == "a") return -CO_MASS_;
    if (ion_type == "b") return 0.0;
    if (ion_type == "y") return WATER_MASS_;
    if (ion_type == "yb") return 0.0;
    if (ion_type == "ya") return -CO_MASS_;
    return 0.0;
  }

  std::vector<MS3FragmentMatcher::TheoreticalMass> MS3FragmentMatcher::computeTheoreticalMasses(
    const std::string& /*subsequence*/,
    const std::vector<std::string>& /*ion_types*/,
    const std::vector<FragmentAnalysis::PTMSite>& /*ptm_sites*/)
  {
    return {}; // stub — implemented in Task 2
  }

  int MS3FragmentMatcher::matchSpectrum(
    const DeconvolvedSpectrum& /*spectrum*/,
    const std::vector<TheoreticalMass>& /*theoretical*/,
    double /*tolerance_ppm*/,
    std::vector<double>* /*ppm_errors*/)
  {
    return 0; // stub — implemented in Task 3
  }

  std::string MS3FragmentMatcher::extractSubsequence(
    const std::string& /*protein_sequence*/,
    const ProteoformContext& /*ctx*/,
    char /*fragment_ion_type*/,
    int /*fragment_ion_index*/)
  {
    return ""; // stub — implemented in Task 4
  }

  std::vector<FragmentAnalysis::PTMSite> MS3FragmentMatcher::rebasePTMSites(
    const std::vector<FragmentAnalysis::PTMSite>& /*ptm_sites*/,
    int /*subseq_start_in_proteoform*/,
    int /*subseq_length*/)
  {
    return {}; // stub — implemented in Task 4
  }

  std::vector<double> MS3FragmentMatcher::calibrateAndScore(
    const std::vector<const DeconvolvedSpectrum*>& /*variant_spectra*/,
    const std::string& /*protein_sequence*/,
    const ProteoformContext& /*ctx*/,
    char /*fragment_ion_type*/,
    int /*fragment_ion_index*/,
    double /*loose_tolerance_ppm*/,
    double /*tight_tolerance_ppm*/)
  {
    return {}; // stub — implemented in Task 4
  }

} // namespace OpenMS
