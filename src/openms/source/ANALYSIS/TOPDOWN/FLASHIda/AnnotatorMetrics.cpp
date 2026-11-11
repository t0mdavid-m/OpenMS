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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/AnnotatorMetrics.h>

#include <algorithm>
#include <numeric>

namespace OpenMS
{
  Param AnnotatorMetrics::parameters(const std::vector<std::string>& ion_types, const std::vector<std::string>& fixed_mod, double tol_ppm)
  {
    Param p = FLASHAnnotatorAlgorithm().getDefaults();
    p.setValue("ion_type", ion_types);
    p.setValue("fixed_mod", fixed_mod);
    p.setValue("internal_ions", "false");
    if (tol_ppm > 0.0) p.setValue("tol", tol_ppm);
    return p;
  }

  bool AnnotatorMetrics::ms2PrSM(const FragmentAnalysis::SearchedPrSM& searched, FLASHAnnotatorAlgorithm::PrSM& out, double& precursor_mass)
  {
    out = FLASHAnnotatorAlgorithm::PrSM();
    precursor_mass = -1.0;
    if (!searched.found || searched.protein_sequence.empty()) return false;
    const ProteinHit& hit = searched.hit;
    const int length = static_cast<int>(searched.protein_sequence.size());
    // 1-based on the hit, <= 0 for an open terminus (FragmentAnalysis::matchFromPrSM reads them the same way)
    const int start = hit.metaValueExists("StartPosition") ? static_cast<int>(hit.getMetaValue("StartPosition")) : -1;
    const int end = hit.metaValueExists("EndPosition") ? static_cast<int>(hit.getMetaValue("EndPosition")) : -1;
    const int first = start > 0 ? std::min(length, start - 1) : 0;
    const int stop = end > 0 ? std::max(first, std::min(length, end)) : length;
    if (stop <= first) return false;
    out.sequence = searched.protein_sequence.substr(static_cast<size_t>(first), static_cast<size_t>(stop - first));
    out.known_n_term = start > 0;
    out.known_c_term = end > 0;
    if (hit.metaValueExists("Modifications") && hit.metaValueExists("ModificationStarts") && hit.metaValueExists("ModificationEnds"))
    {
      const std::vector<double> masses = hit.getMetaValue("Modifications");
      const std::vector<int> starts = hit.getMetaValue("ModificationStarts"), ends = hit.getMetaValue("ModificationEnds");
      for (Size m = 0; m < masses.size() && m < starts.size() && m < ends.size(); m++)
      {
        FLASHAnnotatorAlgorithm::Modification mod;
        mod.mass = masses[m];
        mod.start = starts[m] - first; // 0-based on the entry, as FLASHTnT's toRawPrSM_ reads them
        mod.end = ends[m] - first;
        out.modifications.push_back(mod);
      }
    }
    precursor_mass = hit.metaValueExists("GivenMass") ? static_cast<double>(hit.getMetaValue("GivenMass")) : -1.0;
    return true;
  }

  bool AnnotatorMetrics::ms3TargetPrSM(const std::string& protein_sequence, const MS3FragmentMatcher::ProteoformContext& ctx, char ion_type,
                                       int ion_index, int charge_cap, FLASHAnnotatorAlgorithm::PrSM& out)
  {
    out = FLASHAnnotatorAlgorithm::PrSM();
    if (!MS3FragmentMatcher::isAnnotatableMs3Precursor(ion_type)) return false;
    const std::string subseq = MS3FragmentMatcher::extractSubsequence(protein_sequence, ctx, ion_type, ion_index);
    if (subseq.empty()) return false;
    const int length = static_cast<int>(subseq.size());
    // the target's first residue in the proteoform, 1-based: the frame MS3FragmentMatcher::calibrateAndScore rebases in
    const int proteoform_length = ctx.region_end - ctx.region_start;
    const int start_1based = ion_type == 'b' ? 1 : proteoform_length - ion_index + 1;
    out.sequence = subseq;
    out.precursor_charge = charge_cap;
    for (const FragmentAnalysis::PTMSite& site : MS3FragmentMatcher::rebasePTMSites(ctx.ptm_sites, start_1based, ion_index))
    {
      FLASHAnnotatorAlgorithm::Modification mod;
      mod.mass = site.mass_shift;
      mod.start = site.start_position - 1;
      mod.end = site.end_position - 1;
      out.modifications.push_back(mod);
    }
    if (ion_type == 'b')
    {
      FLASHAnnotatorAlgorithm::Modification water;
      water.mass = B_TARGET_WATER_LOSS;
      water.start = length - 1;
      water.end = length - 1;
      out.modifications.push_back(water);
    }
    return true;
  }

  AnnotatorMetrics::Values AnnotatorMetrics::annotate(FLASHAnnotatorAlgorithm& annotator, const MSSpectrum& raw, FLASHAnnotatorAlgorithm::PrSM prsm)
  {
    Values values;
    std::vector<FLASHAnnotatorAlgorithm::PeakArrays> spectra(1);
    FLASHAnnotatorAlgorithm::toCentroidedPeakArrays(raw, spectra[0]);
    if (spectra[0].mzs.empty() || prsm.sequence.empty()) return values;
    prsm.spectrum_index = 0;
    const std::vector<FLASHAnnotatorAlgorithm::PrSM> prsms(1, prsm);
    std::vector<FLASHAnnotatorAlgorithm::Result> results;
    annotator.annotate(spectra, prsms, results);
    if (results.empty() || !results[0].annotated) return values;
    const FLASHAnnotatorAlgorithm::Result& r = results[0];
    values.raw_score = r.score;
    values.bond_coverage = r.coverage >= 0.0 ? r.coverage / 100.0 : -1.0;
    const double total = std::accumulate(spectra[0].intensities.begin(), spectra[0].intensities.end(), 0.0);
    values.explained_signal = (r.explained_intensity >= 0.0 && total > 0.0) ? r.explained_intensity / total : -1.0;
    return values;
  }

  AnnotatorMetrics::Values AnnotatorMetrics::annotateSearched(FLASHAnnotatorAlgorithm& annotator, const MSSpectrum& raw,
                                                              const FragmentAnalysis::SearchedPrSM& searched)
  {
    FLASHAnnotatorAlgorithm::PrSM prsm;
    double precursor_mass = -1.0;
    if (!ms2PrSM(searched, prsm, precursor_mass)) return Values();
    int reason = 4;
    prsm.precursor_charge = annotator.precursorChargeCap(raw, precursor_mass, reason);
    return annotate(annotator, raw, prsm);
  }
} // namespace OpenMS
