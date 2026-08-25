// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/config.h>

#include <string>
#include <vector>

namespace OpenMS
{

  /**
   * @brief MS3 fragment matching via direct theoretical mass calculation against a precursor subsequence.
   *
   * Scores MS3 exploration CE variants by computing theoretical fragment masses from the
   * MS3 precursor's subsequence (not the full protein), matching deconvolved masses against them,
   * and applying two-pass mass calibration (loose tolerance for calibration, tight for scoring).
   *
   * Handles MS3-specific ion types:
   * - b-precursor: a, b (same direction) + yb, ya (cross-direction, no water)
   * - y-precursor: a, b, y (standard)
   *
   * Supports PTM-aware dual theoretical masses for ambiguous modification regions.
   */
  class OPENMS_DLLAPI MS3FragmentMatcher
  {
  public:
    /// A single theoretical fragment mass entry
    struct TheoreticalMass
    {
      double mass = 0.0;
      int position = 0;          ///< 1-based fragment index from the relevant terminus
      std::string ion_type;      ///< "a", "b", "y", "yb", "ya"
      bool includes_ptm = false; ///< For ambiguous PTMs: true = mass includes PTM shift
      double ambiguous_included = 0.0; ///< Sum of ambiguous/fully-covered PTM mass folded into this ion's mass (for the bare-backbone frame fix)
    };

    /// Detail of a single observed-to-theoretical match
    struct MatchDetail
    {
      double observed_mass = 0.0;    ///< Observed (deconvolved) mass
      double theoretical_mass = 0.0; ///< Matched theoretical mass
      double ppm_error = 0.0;        ///< Signed ppm error: (obs - theo) / theo * 1e6
      int position = 0;              ///< 1-based fragment index from relevant terminus
      std::string ion_type;          ///< "a", "b", "y", "yb", "ya"
      bool includes_ptm = false;     ///< Whether the matched theoretical includes an ambiguous PTM
      double ambiguous_included = 0.0; ///< PTM mass folded into the matched theoretical (copied from TheoreticalMass; used for adjusted/theoretical frame masses)
      /// PeakGroup::getQscore() of the MS3 deconvolved mass this detail was matched from, at the
      /// representative charge. Carries the MS3 leg of identification.tsv's fragment_qscores; the MS2
      /// leg is filled independently in FragmentAnalysis::runTagBasedFragmentMatching_. Unaffected by
      /// the calibration copy: that rescales the monoisotopic mass only and never re-scores the group.
      double qscore = 0.0;
    };

    /// Cached proteoform context from MS2 tag-based matching
    struct ProteoformContext
    {
      int region_start = -1; ///< 0-based start position in protein sequence
      int region_end = -1;   ///< 0-based exclusive end position in protein sequence
      std::vector<FragmentAnalysis::PTMSite> ptm_sites; ///< 1-based positions relative to proteoform
    };

    /// Compile-time constant: loose tolerance for calibration pass
    static constexpr double LOOSE_TOLERANCE_PPM = 500.0;

    // -- Ion type handling --

    /// The ion classes this matcher can project an MS3 sub-fragment through: a, b, c (prefix) and x, y, z
    /// (suffix). Everything else -- notably 'u', the label an exhaustive-mode unassigned target carries
    /// (ADR-0023 decision 5) -- has no frame, and every projection site below must refuse rather than fall
    /// through to its suffix branch.
    ///
    /// The refusal keys on the CLASS and NEVER on the index: ion type and ion index are two independent
    /// fields that travel independently, so an index-only guard leaves a 'u' carrying a plausible index
    /// cutting a real suffix out of the proteoform and matching against it -- confident wrong
    /// identifications instead of none.
    static bool isKnownIonClass(char ion_class);

    /// Select MS3 ion types based on precursor fragment class. Empty for an unknown class (isKnownIonClass).
    static std::vector<std::string> getMS3IonTypes(char precursor_ion_class);

    /// Returns true for prefix ion types (a, b), false for suffix (y, yb, ya)
    static bool isPrefixIonType(const std::string& ion_type);

    /// Returns the ion mass shift in Da
    static double getIonShift(const std::string& ion_type);

    // -- Theoretical mass calculation --

    /// Compute theoretical fragment masses for a subsequence with optional PTM handling
    static std::vector<TheoreticalMass> computeTheoreticalMasses(
      const std::string& subsequence,
      const std::vector<std::string>& ion_types,
      const std::vector<FragmentAnalysis::PTMSite>& ptm_sites = {});

    // -- Matching --

    /// Match deconvolved masses against theoretical, return match count
    static int matchSpectrum(
      const DeconvolvedSpectrum& spectrum,
      const std::vector<TheoreticalMass>& theoretical,
      double tolerance_ppm,
      std::vector<double>* ppm_errors = nullptr,
      std::vector<MatchDetail>* match_details = nullptr);

    // -- Subsequence extraction + PTM rebasing --

    /// Extract the precursor fragment's subsequence from the protein
    static std::string extractSubsequence(
      const std::string& protein_sequence,
      const ProteoformContext& ctx,
      char fragment_ion_type,
      int fragment_ion_index);

    /// Rebase PTM sites from proteoform coordinates to subsequence coordinates
    static std::vector<FragmentAnalysis::PTMSite> rebasePTMSites(
      const std::vector<FragmentAnalysis::PTMSite>& ptm_sites,
      int subseq_start_in_proteoform,
      int subseq_length);

    /// Render the MS3 fragment sub-sequence as ProForma: extractSubsequence + clip/rebase the parent
    /// mods (rebasePTMSites) into the fragment frame + toProForma. Uses the SAME subseq-start logic as
    /// calibrateAndScore. Returns "" when the context isn't populated (ion_type=='\0' | index==0 |
    /// invalid region) — callers leave the proteoform empty, never a parent fallback. Consumed by the
    /// scan_results MS3 rows so the fragment proteoform is present even when identification is deferred.
    static std::string fragmentProForma(
      const std::string& protein_sequence,
      const ProteoformContext& ctx,
      char fragment_ion_type,
      int fragment_ion_index);

    /// The (subseq, wide fragment-frame PTM ranges) that fragmentProForma renders, exposed as the a-priori
    /// wide base the identification.tsv leaf narrower (FragmentAnalysis::narrowFragmentPTMSites) re-brackets
    /// over ONE MS3 scan's matched EQUIVALENT ions. Under Change L the leaf is a per-scan equiv-frame bracket
    /// that MAY exceed this base and MAY merge co-observed mods -- there is NO leaf-subset-of-scan_commands
    /// guarantee. Returns false when the context isn't populated.
    static bool fragmentProFormaSites(
      const std::string& protein_sequence,
      const ProteoformContext& ctx,
      char fragment_ion_type,
      int fragment_ion_index,
      std::string& out_subseq,
      std::vector<FragmentAnalysis::PTMSite>& out_wide_sites);

    // -- Two-pass calibration pipeline --

    /// Score all CE variants via two-pass calibration + matching
    static std::vector<double> calibrateAndScore(
      const std::vector<const DeconvolvedSpectrum*>& variant_spectra,
      const std::string& protein_sequence,
      const ProteoformContext& ctx,
      char fragment_ion_type,
      int fragment_ion_index,
      double loose_tolerance_ppm,
      double tight_tolerance_ppm,
      std::vector<FragmentAnalysis::ProteoformMatch>* detailed_results = nullptr);

    /// Compute the equivalent full-protein ion for an MS3 fragment match.
    static void computeEquivalentIon(
      const std::string& protein_sequence,
      const ProteoformContext& ctx,
      char precursor_ion_type,
      int precursor_ion_index,
      const std::string& ms3_ion_type,
      int ms3_ion_index,
      double ms3_theoretical_mass,
      const std::vector<double>& protein_prefix_masses,
      std::string& equiv_type,
      int& equiv_index,
      double& mass_offset);

    /// Compute prefix mass sums for a protein sequence (fixed PTMs included).
    /// prefix[i] = sum of internal residue masses for positions [0, i).
    static std::vector<double> computeProteinPrefixMasses(
      const std::string& protein_sequence,
      const std::vector<FragmentAnalysis::PTMSite>& ptm_sites = {},
      int proteoform_start = 0);

    /// Sum the ambiguous PTM mass (start!=end; fixed PTMs already live in the prefix masses) that falls
    /// within the EQUIVALENT full-protein ion's own residue coverage (prefix b/a of index k -> absolute
    /// residues [0, k-1]; suffix y of index k -> [P-k, P-1]). This is the mod mass the equivalent ion
    /// itself carries -- correctly DROPPING the mods that belong to the complement when an MS3 suffix ion
    /// maps to a prefix equivalent (and vice versa), while reproducing the sub-frame ambiguous_included
    /// exactly for same-direction maps. @p sub_includes_ptm is the MS3 match's own PTM verdict, used only
    /// for a mod that straddles the equivalent ion's boundary.
    static double coveredAmbiguousInEquivFrame(
      const ProteoformContext& ctx,
      const std::string& equiv_type,
      int equiv_index,
      int protein_length,
      bool sub_includes_ptm,
      bool is_flip);

  private:
    static constexpr double PROTON_MASS_ = 1.007276;
    static constexpr double WATER_MASS_ = 18.010565;
    static constexpr double CO_MASS_ = 27.994915;
  };

} // namespace OpenMS
