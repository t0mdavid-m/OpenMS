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
    const std::string& subsequence,
    const std::vector<std::string>& ion_types,
    const std::vector<FragmentAnalysis::PTMSite>& ptm_sites)
  {
    int n = static_cast<int>(subsequence.size());
    if (n < 2) return {};

    // Get residue masses (0-based indexing)
    std::vector<double> res_mass(n, 0.0);
    for (int i = 0; i < n; ++i)
    {
      const Residue* res = ResidueDB::getInstance()->getResidue(subsequence[i]);
      if (res != nullptr)
        res_mass[i] = res->getMonoWeight(Residue::Internal);
    }

    // Separate fixed and ambiguous PTMs (convert 1-based PTMSite to 0-based)
    struct PTM0
    {
      int start;
      int end;
      double mass;
      bool fixed;
    };
    std::vector<PTM0> ptms;
    for (const auto& p : ptm_sites)
    {
      PTM0 pm;
      pm.start = p.start_position - 1;
      pm.end = p.end_position - 1;
      pm.mass = p.mass_shift;
      pm.fixed = (p.start_position == p.end_position);
      ptms.push_back(pm);
    }

    // Pre-compute cumulative fixed PTM contributions
    // fixed_prefix[i] = total fixed PTM mass at positions <= i
    // fixed_suffix[i] = total fixed PTM mass at positions >= i
    std::vector<double> fixed_prefix(n, 0.0);
    std::vector<double> fixed_suffix(n, 0.0);
    for (const auto& pm : ptms)
    {
      if (! pm.fixed) continue;
      if (pm.start < 0 || pm.start >= n) continue;
      for (int i = pm.start; i < n; ++i)
        fixed_prefix[i] += pm.mass;
      for (int i = 0; i <= pm.start; ++i)
        fixed_suffix[i] += pm.mass;
    }

    std::vector<TheoreticalMass> result;

    for (const auto& ion_type : ion_types)
    {
      bool is_prefix = isPrefixIonType(ion_type);
      double shift = getIonShift(ion_type);

      double cumulative = 0.0;
      for (int i = 0; i < n - 1; ++i)
      {
        double base_mass;
        int frag_start_0, frag_end_0; // 0-based inclusive range covered by this ion

        if (is_prefix)
        {
          cumulative += res_mass[i];
          base_mass = cumulative + shift + fixed_prefix[i];
          frag_start_0 = 0;
          frag_end_0 = i;
        }
        else
        {
          int idx = n - 1 - i;
          cumulative += res_mass[idx];
          base_mass = cumulative + shift + fixed_suffix[idx];
          frag_start_0 = idx;
          frag_end_0 = n - 1;
        }

        // Check ambiguous PTMs for this position
        double ambiguous_delta = 0.0;
        bool has_ambiguous = false;
        for (const auto& pm : ptms)
        {
          if (pm.fixed) continue;
          // Is the ambiguous range fully covered by this ion?
          if (pm.start >= frag_start_0 && pm.end <= frag_end_0)
          {
            base_mass += pm.mass; // fully covered, always include
          }
          // Is it partially overlapping?
          else if (pm.end >= frag_start_0 && pm.start <= frag_end_0)
          {
            ambiguous_delta += pm.mass;
            has_ambiguous = true;
          }
          // else: no overlap, skip
        }

        int position = i + 1; // 1-based from the relevant terminus

        if (has_ambiguous)
        {
          TheoreticalMass tm_with;
          tm_with.mass = base_mass + ambiguous_delta;
          tm_with.position = position;
          tm_with.ion_type = ion_type;
          tm_with.includes_ptm = true;
          result.push_back(tm_with);

          TheoreticalMass tm_without;
          tm_without.mass = base_mass;
          tm_without.position = position;
          tm_without.ion_type = ion_type;
          tm_without.includes_ptm = false;
          result.push_back(tm_without);
        }
        else
        {
          TheoreticalMass tm;
          tm.mass = base_mass;
          tm.position = position;
          tm.ion_type = ion_type;
          tm.includes_ptm = false;
          result.push_back(tm);
        }
      }
    }

    return result;
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
