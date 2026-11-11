// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/config.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
   * @brief What the annotator exploration metrics read (ADR-0074): one PrSM annotated on one raw spectrum by
   * FLASHAnnotatorAlgorithm, terminal fragments only.
   *
   * Pure static functions; nothing is owned here. The caller owns the FLASHAnnotatorAlgorithm (annotate() caches its
   * isotope table, so one instance per thread).
   *
   * - raw_score: the annotator's score, -log10 of the binomial chance of the matched fragment hypotheses at the
   *   spectrum's own decoy-measured chance rate.
   * - bond_coverage: bonds covered by an accepted terminal envelope over L - 1, a FRACTION (the annotator's coverage / 100).
   * - explained_signal: the raw intensity of the peaks of accepted terminal envelopes, each peak once, over the
   *   intensity of the centroided peak list the annotator matched against.
   *
   * At MS2 raw_score and bond_coverage are the search's own (FLASHTnT annotates all PrSMs of a spectrum in one call, and
   * the chance rate is pooled over them, so a second annotation of the best PrSM alone need not reproduce RawScore); only
   * explained_signal is taken from annotateSearched. At MS3 all three come from annotate on the target's PrSM.
   */
  class OPENMS_DLLAPI AnnotatorMetrics
  {
  public:
    struct OPENMS_DLLAPI Values
    {
      double raw_score = -1;
      double bond_coverage = -1;
      double explained_signal = -1;
    };

    /// the mass a b target lacks against a chain of its own residues: its C-terminus carries no water
    static constexpr double B_TARGET_WATER_LOSS = -18.010565;

    /// The Param FLASHTnT hands its own annotator for the live search: the annotator's defaults, the ion types, the
    /// fixed modifications, internal fragments off; tol = @p tol_ppm when it is positive.
    static Param parameters(const std::vector<std::string>& ion_types, const std::vector<std::string>& fixed_mod, double tol_ppm = -1.0);

    /// The PrSM of a search hit: the entry's sequence cut by the hit's StartPosition / EndPosition (an open terminus runs
    /// to the entry's end and is marked unknown), its Modifications rebased onto the cut. The engine-side counterpart of
    /// FLASHTnT's private toRawPrSM_. @p precursor_mass: the hit's GivenMass, -1 without one. False without a hit or sequence.
    static bool ms2PrSM(const FragmentAnalysis::SearchedPrSM& searched, FLASHAnnotatorAlgorithm::PrSM& out, double& precursor_mass);

    /// An MS3 target as a chain of its own: the target fragment's subsequence (MS3FragmentMatcher::extractSubsequence), the
    /// proteoform's modifications inside it (rebasePTMSites, so a range that leaves the target is clipped to it, as the
    /// MS3 matcher reads it), fragment charges capped at @p charge_cap. A y target is a truncated chain. A b target
    /// carries B_TARGET_WATER_LOSS on its last residue: the annotator adds a modification to every fragment that holds
    /// its whole range, so to every suffix ion and to no prefix ion. False when the ion is not b/y
    /// (MS3FragmentMatcher::isAnnotatableMs3Precursor) or the context gives no subsequence.
    static bool ms3TargetPrSM(const std::string& protein_sequence, const MS3FragmentMatcher::ProteoformContext& ctx, char ion_type,
                              int ion_index, int charge_cap, FLASHAnnotatorAlgorithm::PrSM& out);

    /// @p prsm annotated on @p raw (centroided here if it is a profile spectrum). All -1 when nothing was annotated.
    static Values annotate(FLASHAnnotatorAlgorithm& annotator, const MSSpectrum& raw, FLASHAnnotatorAlgorithm::PrSM prsm);

    /// ms2PrSM, the precursor charge cap FLASHTnT applies (FLASHAnnotatorAlgorithm::precursorChargeCap), then annotate.
    static Values annotateSearched(FLASHAnnotatorAlgorithm& annotator, const MSSpectrum& raw, const FragmentAnalysis::SearchedPrSM& searched);
  };
} // namespace OpenMS
