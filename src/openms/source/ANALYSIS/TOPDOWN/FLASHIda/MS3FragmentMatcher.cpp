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
    const DeconvolvedSpectrum& spectrum,
    const std::vector<TheoreticalMass>& theoretical,
    double tolerance_ppm,
    std::vector<double>* ppm_errors)
  {
    if (spectrum.empty() || theoretical.empty()) return 0;

    // Sort theoretical masses for binary search
    std::vector<size_t> sorted_idx(theoretical.size());
    std::iota(sorted_idx.begin(), sorted_idx.end(), 0);
    std::sort(sorted_idx.begin(), sorted_idx.end(),
      [&theoretical](size_t a, size_t b) { return theoretical[a].mass < theoretical[b].mass; });

    std::vector<bool> theo_used(theoretical.size(), false);
    int match_count = 0;

    // For each deconvolved mass, find closest unused theoretical mass within tolerance
    for (Size si = 0; si < spectrum.size(); ++si)
    {
      double obs_mass = spectrum[si].getMonoMass();
      double tol_da = obs_mass * tolerance_ppm * 1e-6;

      int best_theo_idx = -1;
      double best_ppm = tolerance_ppm + 1.0;

      // Binary search for candidates
      for (size_t k = 0; k < sorted_idx.size(); ++k)
      {
        size_t ti = sorted_idx[k];
        if (theo_used[ti]) continue;
        double theo_mass = theoretical[ti].mass;
        if (theo_mass < obs_mass - tol_da) continue;
        if (theo_mass > obs_mass + tol_da) break;

        double ppm_err = std::abs((obs_mass - theo_mass) / theo_mass) * 1e6;
        if (ppm_err < best_ppm)
        {
          best_ppm = ppm_err;
          best_theo_idx = static_cast<int>(ti);
        }
      }

      if (best_theo_idx >= 0)
      {
        theo_used[best_theo_idx] = true;
        ++match_count;
        if (ppm_errors)
        {
          double signed_ppm = (obs_mass - theoretical[best_theo_idx].mass)
                              / theoretical[best_theo_idx].mass * 1e6;
          ppm_errors->push_back(signed_ppm);
        }
      }
    }

    return match_count;
  }

  std::string MS3FragmentMatcher::extractSubsequence(
    const std::string& protein_sequence,
    const ProteoformContext& ctx,
    char fragment_ion_type,
    int fragment_ion_index)
  {
    if (ctx.region_start < 0 || ctx.region_end < 0) return "";
    int proteoform_length = ctx.region_end - ctx.region_start;
    if (fragment_ion_index <= 0 || fragment_ion_index > proteoform_length) return "";

    int subseq_start_0based; // absolute position in protein_sequence
    if (fragment_ion_type == 'b' || fragment_ion_type == 'a' || fragment_ion_type == 'c')
    {
      subseq_start_0based = ctx.region_start;
    }
    else // y, x, z
    {
      subseq_start_0based = ctx.region_end - fragment_ion_index;
    }

    if (subseq_start_0based < 0 ||
        subseq_start_0based + fragment_ion_index > static_cast<int>(protein_sequence.size()))
      return "";

    return protein_sequence.substr(subseq_start_0based, fragment_ion_index);
  }

  std::vector<FragmentAnalysis::PTMSite> MS3FragmentMatcher::rebasePTMSites(
    const std::vector<FragmentAnalysis::PTMSite>& ptm_sites,
    int subseq_start_in_proteoform,
    int subseq_length)
  {
    int subseq_end = subseq_start_in_proteoform + subseq_length - 1; // 1-based inclusive end

    std::vector<FragmentAnalysis::PTMSite> result;
    for (const auto& ptm : ptm_sites)
    {
      // No overlap — skip
      if (ptm.end_position < subseq_start_in_proteoform || ptm.start_position > subseq_end)
        continue;

      FragmentAnalysis::PTMSite rebased;
      rebased.start_position = std::max(ptm.start_position, subseq_start_in_proteoform)
                               - subseq_start_in_proteoform + 1;
      rebased.end_position = std::min(ptm.end_position, subseq_end)
                             - subseq_start_in_proteoform + 1;
      rebased.position = (rebased.start_position + rebased.end_position) / 2;
      rebased.mass_shift = ptm.mass_shift;
      result.push_back(rebased);
    }
    return result;
  }

  std::vector<double> MS3FragmentMatcher::calibrateAndScore(
    const std::vector<const DeconvolvedSpectrum*>& variant_spectra,
    const std::string& protein_sequence,
    const ProteoformContext& ctx,
    char fragment_ion_type,
    int fragment_ion_index,
    double loose_tolerance_ppm,
    double tight_tolerance_ppm)
  {
    std::vector<double> scores(variant_spectra.size(), 0.0);

    // Extract subsequence
    std::string subseq = extractSubsequence(protein_sequence, ctx, fragment_ion_type, fragment_ion_index);
    if (subseq.empty()) return scores;

    // Determine subsequence start in proteoform (1-based)
    int proteoform_length = ctx.region_end - ctx.region_start;
    int subseq_start_1based;
    if (fragment_ion_type == 'b' || fragment_ion_type == 'a' || fragment_ion_type == 'c')
      subseq_start_1based = 1;
    else
      subseq_start_1based = proteoform_length - fragment_ion_index + 1;

    // Rebase PTMs
    auto rebased_ptms = rebasePTMSites(ctx.ptm_sites, subseq_start_1based, fragment_ion_index);

    // Ion types
    auto ion_types = getMS3IonTypes(fragment_ion_type);

    // Compute theoretical masses once
    auto theoretical = computeTheoreticalMasses(subseq, ion_types, rebased_ptms);
    if (theoretical.empty()) return scores;

    // Pass 1: loose matching for calibration
    std::vector<double> all_ppm_errors;
    for (size_t vi = 0; vi < variant_spectra.size(); ++vi)
    {
      if (variant_spectra[vi] == nullptr || variant_spectra[vi]->empty()) continue;
      matchSpectrum(*variant_spectra[vi], theoretical, loose_tolerance_ppm, &all_ppm_errors);
    }

    // Compute median ppm error
    double correction_factor = 1.0;
    if (! all_ppm_errors.empty())
    {
      std::sort(all_ppm_errors.begin(), all_ppm_errors.end());
      double median_ppm;
      size_t mid = all_ppm_errors.size() / 2;
      if (all_ppm_errors.size() % 2 == 0)
        median_ppm = (all_ppm_errors[mid - 1] + all_ppm_errors[mid]) / 2.0;
      else
        median_ppm = all_ppm_errors[mid];

      correction_factor = 1.0 / (1.0 + median_ppm * 1e-6);
    }

    // Pass 2: apply correction, match at tight tolerance
    for (size_t vi = 0; vi < variant_spectra.size(); ++vi)
    {
      if (variant_spectra[vi] == nullptr || variant_spectra[vi]->empty()) continue;

      // Create corrected copy
      DeconvolvedSpectrum corrected(0);
      for (Size pi = 0; pi < variant_spectra[vi]->size(); ++pi)
      {
        PeakGroup pg = (*variant_spectra[vi])[pi];
        pg.setMonoisotopicMass(pg.getMonoMass() * correction_factor);
        corrected.push_back(pg);
      }

      scores[vi] = static_cast<double>(matchSpectrum(corrected, theoretical, tight_tolerance_ppm));
    }

    return scores;
  }

} // namespace OpenMS
