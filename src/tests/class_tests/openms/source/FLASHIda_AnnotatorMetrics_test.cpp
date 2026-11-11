// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

// AnnotatorMetrics (ADR-0074): what the raw_score / bond_coverage / explained_signal exploration metrics read from
// FLASHAnnotatorAlgorithm, and the one field the annotator gained for them (Result::explained_intensity). Every
// spectrum here is synthetic: the predicted isotope envelopes of a known set of fragments, built with the annotator's
// own getIsotopeEnvelope, the way FLASHAnnotatorAlgorithm_test builds its spectra. No engine is driven; the sweep
// behaviour is FLASHIda_exploration_test's.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/AnnotatorMetrics.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
///////////////////////////

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  /// the first 40 residues of horse cytochrome c
  const std::string kChain = "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKT";

  /// a context whose proteoform is the whole of kChain, with the given modifications (1-based, inclusive)
  MS3FragmentMatcher::ProteoformContext wholeChain(const std::vector<FragmentAnalysis::PTMSite>& sites = {})
  {
    MS3FragmentMatcher::ProteoformContext ctx;
    ctx.region_start = 0;
    ctx.region_end = static_cast<int>(kChain.size());
    ctx.protein_sequence = kChain;
    ctx.ptm_sites = sites;
    return ctx;
  }

  FragmentAnalysis::PTMSite site(int start, int end, double mass)
  {
    FragmentAnalysis::PTMSite s;
    s.start_position = start;
    s.end_position = end;
    s.position = (start + end) / 2;
    s.mass_shift = mass;
    return s;
  }

  /// the engine's annotator Param (AnnotatorMetrics::parameters) with these ion types and no fixed modification
  void configure(FLASHAnnotatorAlgorithm& a, const std::vector<std::string>& ion_types)
  {
    a.setParameters(AnnotatorMetrics::parameters(ion_types, std::vector<std::string>()));
  }

  /// A centroided spectrum of the isotope envelopes of @p masses, each at charges 1..@p max_charge.
  MSSpectrum envelopesOf(FLASHAnnotatorAlgorithm& algo, const std::vector<double>& masses, int max_charge)
  {
    MSSpectrum spec;
    std::vector<double> mzs, intensities;
    for (Size i = 0; i < masses.size(); i++)
    {
      for (int z = 1; z <= max_charge; z++)
      {
        algo.getIsotopeEnvelope(masses[i], z, mzs, intensities);
        for (Size j = 0; j < mzs.size(); j++) spec.emplace_back(mzs[j], static_cast<float>(intensities[j] * (1000.0 + 10.0 * i)));
      }
    }
    spec.sortByPosition();
    spec.setType(SpectrumSettings::CENTROID);
    return spec;
  }

  /// the monoisotopic masses of the terminal fragments @p algo predicts for @p prsm
  std::vector<double> fragmentMasses(const FLASHAnnotatorAlgorithm& algo, const FLASHAnnotatorAlgorithm::PrSM& prsm)
  {
    std::vector<FLASHAnnotatorAlgorithm::Fragment> fragments;
    algo.getFragments(prsm, fragments);
    std::vector<double> masses;
    for (const auto& f : fragments) masses.push_back(f.mono_mass);
    return masses;
  }

  double totalIntensity(const MSSpectrum& spec)
  {
    double total = 0.0;
    for (const auto& p : spec) total += p.getIntensity();
    return total;
  }
} // namespace

START_TEST(AnnotatorMetrics, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// The annotator's new field (the owner-approved edit inside the FLASHTnT boundary, ADR-0074 d.10). z, z+1 and z+2 of
// one length put their envelopes on the same peaks. The spectrum holds the z envelopes only; annotated with all three
// types, the envelopes of the other two are accepted on those same peaks, so the envelope intensities summed count
// them more than once, and explained_intensity does not: it can never exceed what the spectrum holds.
START_SECTION(explained_intensity_counts_a_shared_peak_once)
{
  FLASHAnnotatorAlgorithm::PrSM prsm;
  prsm.sequence = kChain;
  prsm.precursor_charge = 2;
  prsm.spectrum_index = 0;

  FLASHAnnotatorAlgorithm z_only;
  configure(z_only, {"z"});
  const MSSpectrum spec = envelopesOf(z_only, fragmentMasses(z_only, prsm), 2);
  const double total = totalIntensity(spec);
  ABORT_IF(spec.empty() || total <= 0.0)

  FLASHAnnotatorAlgorithm all_three;
  configure(all_three, {"z", "zp1", "zp2"});
  std::vector<FLASHAnnotatorAlgorithm::PeakArrays> spectra(1);
  FLASHAnnotatorAlgorithm::toCentroidedPeakArrays(spec, spectra[0]);
  std::vector<FLASHAnnotatorAlgorithm::Result> results;
  all_three.annotate(spectra, std::vector<FLASHAnnotatorAlgorithm::PrSM>(1, prsm), results);
  ABORT_IF(results.size() != 1 || !results[0].annotated)
  const auto& r = results[0];

  double envelope_sum = 0.0;
  int shifted_envelopes = 0;   // envelopes of zp1 / zp2: the ones that sit on the z peaks
  for (const auto& e : r.envelopes)
  {
    envelope_sum += e.intensity;
    if (all_three.getIonTypes()[static_cast<Size>(e.ion_type)] != "z") ++shifted_envelopes;
  }
  TEST_TRUE(shifted_envelopes > 0)                       // non-vacuity: peaks really are shared
  TEST_TRUE(r.explained_intensity > 0.5 * total)         // the spectrum is its z envelopes
  TEST_TRUE(r.explained_intensity <= total * (1.0 + 1e-9))
  TEST_TRUE(envelope_sum > r.explained_intensity * (1.0 + 1e-6))   // the plain sum counts the shared peaks again
}
END_SECTION

// explained_signal = explained_intensity over the whole centroided scan. A spectrum made only of the PrSM's own b/y
// envelopes is almost wholly explained; peaks nobody predicted, added above its m/z range, lower it by exactly their
// share of the new total.
START_SECTION(explained_signal_is_one_on_a_fully_explained_spectrum_and_falls_with_added_peaks)
{
  FLASHAnnotatorAlgorithm::PrSM prsm;
  prsm.sequence = kChain;
  prsm.precursor_charge = 2;

  FLASHAnnotatorAlgorithm by;
  configure(by, {"b", "y"});
  MSSpectrum spec = envelopesOf(by, fragmentMasses(by, prsm), 2);
  const double own = totalIntensity(spec);
  ABORT_IF(spec.empty() || own <= 0.0)

  const AnnotatorMetrics::Values clean = AnnotatorMetrics::annotate(by, spec, prsm);
  TEST_TRUE(clean.explained_signal > 0.8)
  TEST_TRUE(clean.explained_signal <= 1.0 + 1e-9)
  TEST_TRUE(clean.raw_score > 0.0)
  TEST_TRUE(clean.bond_coverage > 0.8)

  const double top = spec.back().getMZ();
  double added = 0.0;
  for (int i = 0; i < 20; ++i)
  {
    spec.emplace_back(top + 50.0 + 7.3 * i, 5000.0f);
    added += 5000.0;
  }
  spec.sortByPosition();
  const AnnotatorMetrics::Values noisy = AnnotatorMetrics::annotate(by, spec, prsm);
  TOLERANCE_RELATIVE(1.001)
  TEST_REAL_SIMILAR(noisy.explained_signal, clean.explained_signal * own / (own + added))
}
END_SECTION

// ADR-0074 d.4. A b target has no C-terminal water, so its suffix ions are the MS3 matcher's "yb": the residues alone.
// The spectrum holds exactly those ions. The target's PrSM covers them; the same chain without the correction
// predicts every suffix ion 18.011 Da too heavy and covers none.
START_SECTION(a_b_target_predicts_its_water_less_suffix_ions)
{
  const int index = 30;
  FLASHAnnotatorAlgorithm::PrSM target;
  TEST_EQUAL(AnnotatorMetrics::ms3TargetPrSM(kChain, wholeChain(), 'b', index, 2, target), true)
  TEST_EQUAL(std::string(target.sequence), kChain.substr(0, index))
  TEST_EQUAL(target.precursor_charge, 2)
  ABORT_IF(target.modifications.size() != 1)
  TEST_REAL_SIMILAR(target.modifications[0].mass, AnnotatorMetrics::B_TARGET_WATER_LOSS)
  TEST_EQUAL(target.modifications[0].start, index - 1)
  TEST_EQUAL(target.modifications[0].end, index - 1)

  std::vector<double> yb;
  for (const auto& t : MS3FragmentMatcher::computeTheoreticalMasses(kChain.substr(0, index), {"yb"})) yb.push_back(t.mass);
  TEST_EQUAL(static_cast<int>(yb.size()), index - 1)

  FLASHAnnotatorAlgorithm by;
  configure(by, {"b", "y"});
  const MSSpectrum spec = envelopesOf(by, yb, 2);
  ABORT_IF(spec.empty())

  const AnnotatorMetrics::Values corrected = AnnotatorMetrics::annotate(by, spec, target);
  TEST_TRUE(corrected.bond_coverage > 0.5)

  FLASHAnnotatorAlgorithm::PrSM plain = target;
  plain.modifications.clear();
  const AnnotatorMetrics::Values uncorrected = AnnotatorMetrics::annotate(by, spec, plain);
  TEST_TRUE(uncorrected.bond_coverage >= 0.0)
  TEST_TRUE(uncorrected.bond_coverage < 0.1)
}
END_SECTION

// A y target is a truncated chain: no correction, and the annotator's b and y masses on it are the MS3 matcher's.
START_SECTION(a_y_target_is_a_plain_chain)
{
  const int index = 25;
  FLASHAnnotatorAlgorithm::PrSM target;
  TEST_EQUAL(AnnotatorMetrics::ms3TargetPrSM(kChain, wholeChain(), 'y', index, 3, target), true)
  const std::string subseq = kChain.substr(kChain.size() - index);
  TEST_EQUAL(std::string(target.sequence), subseq)
  TEST_EQUAL(target.modifications.empty(), true)

  FLASHAnnotatorAlgorithm by;
  configure(by, {"b", "y"});
  std::vector<FLASHAnnotatorAlgorithm::Fragment> fragments;
  by.getFragments(target, fragments);
  const auto theoretical = MS3FragmentMatcher::computeTheoreticalMasses(subseq, {"b", "y"});
  int compared = 0;
  for (const auto& f : fragments)
  {
    const std::string type = by.getIonTypes()[static_cast<Size>(f.ion_type)];
    for (const auto& t : theoretical)
    {
      if (t.ion_type != type || t.position != f.length) continue;
      TEST_REAL_SIMILAR(f.mono_mass, t.mass)
      ++compared;
    }
  }
  TEST_EQUAL(compared, 2 * (index - 1))
}
END_SECTION

// The proteoform's modifications inside the target travel with it, in the target's own coordinates (0-based,
// inclusive); one outside it does not.
START_SECTION(a_target_carries_the_modifications_inside_it)
{
  const auto ctx = wholeChain({site(1, 1, 42.010565), site(14, 20, 615.1659), site(35, 35, 79.966331)});

  FLASHAnnotatorAlgorithm::PrSM b32;
  TEST_EQUAL(AnnotatorMetrics::ms3TargetPrSM(kChain, ctx, 'b', 32, 4, b32), true)
  ABORT_IF(b32.modifications.size() != 3)   // the acetyl, the heme range, the water loss; residue 35 is outside
  TEST_REAL_SIMILAR(b32.modifications[0].mass, 42.010565)
  TEST_EQUAL(b32.modifications[0].start, 0)
  TEST_EQUAL(b32.modifications[0].end, 0)
  TEST_REAL_SIMILAR(b32.modifications[1].mass, 615.1659)
  TEST_EQUAL(b32.modifications[1].start, 13)
  TEST_EQUAL(b32.modifications[1].end, 19)
  TEST_REAL_SIMILAR(b32.modifications[2].mass, AnnotatorMetrics::B_TARGET_WATER_LOSS)
  TEST_EQUAL(b32.modifications[2].start, 31)

  FLASHAnnotatorAlgorithm::PrSM y10;
  TEST_EQUAL(AnnotatorMetrics::ms3TargetPrSM(kChain, ctx, 'y', 10, 2, y10), true)   // residues 31..40
  ABORT_IF(y10.modifications.size() != 1)
  TEST_REAL_SIMILAR(y10.modifications[0].mass, 79.966331)
  TEST_EQUAL(y10.modifications[0].start, 4)   // residue 35 of the proteoform
  TEST_EQUAL(y10.modifications[0].end, 4)
}
END_SECTION

// ADR-0074 d.7: only a b or y target is annotated (ADR-0054), and only against an identified proteoform.
START_SECTION(a_c_or_z_target_and_an_empty_context_give_no_prsm)
{
  FLASHAnnotatorAlgorithm::PrSM out;
  TEST_EQUAL(AnnotatorMetrics::ms3TargetPrSM(kChain, wholeChain(), 'c', 10, 2, out), false)
  TEST_EQUAL(AnnotatorMetrics::ms3TargetPrSM(kChain, wholeChain(), 'z', 10, 2, out), false)
  TEST_EQUAL(AnnotatorMetrics::ms3TargetPrSM(kChain, wholeChain(), 'a', 10, 2, out), false)
  TEST_EQUAL(AnnotatorMetrics::ms3TargetPrSM(kChain, wholeChain(), 'u', 0, 2, out), false)
  TEST_EQUAL(AnnotatorMetrics::ms3TargetPrSM("", MS3FragmentMatcher::ProteoformContext(), 'b', 10, 2, out), false)
  TEST_EQUAL(AnnotatorMetrics::ms3TargetPrSM(kChain, wholeChain(), 'b', 10, 2, out), true)   // the control
}
END_SECTION

// The engine-side counterpart of FLASHTnT's toRawPrSM_: the entry cut by the hit's positions, its modifications
// rebased onto the cut, an open terminus marked unknown.
START_SECTION(ms2_prsm_is_cut_from_the_entry_by_the_hits_positions)
{
  FragmentAnalysis::SearchedPrSM searched;
  searched.found = true;
  searched.protein_sequence = kChain;
  searched.hit.setMetaValue("StartPosition", 6);
  searched.hit.setMetaValue("EndPosition", 30);
  searched.hit.setMetaValue("Modifications", DoubleList {79.966331});
  searched.hit.setMetaValue("ModificationStarts", IntList {10});
  searched.hit.setMetaValue("ModificationEnds", IntList {12});
  searched.hit.setMetaValue("GivenMass", 3000.5);

  FLASHAnnotatorAlgorithm::PrSM prsm;
  double precursor_mass = 0.0;
  TEST_EQUAL(AnnotatorMetrics::ms2PrSM(searched, prsm, precursor_mass), true)
  TEST_EQUAL(std::string(prsm.sequence), kChain.substr(5, 25))
  TEST_EQUAL(prsm.known_n_term, true)
  TEST_EQUAL(prsm.known_c_term, true)
  ABORT_IF(prsm.modifications.size() != 1)
  TEST_EQUAL(prsm.modifications[0].start, 5)
  TEST_EQUAL(prsm.modifications[0].end, 7)
  TEST_REAL_SIMILAR(prsm.modifications[0].mass, 79.966331)
  TEST_REAL_SIMILAR(precursor_mass, 3000.5)

  // an open N-terminus runs to the entry's start and is marked unknown
  searched.hit.setMetaValue("StartPosition", -1);
  TEST_EQUAL(AnnotatorMetrics::ms2PrSM(searched, prsm, precursor_mass), true)
  TEST_EQUAL(std::string(prsm.sequence), kChain.substr(0, 30))
  TEST_EQUAL(prsm.known_n_term, false)
  TEST_EQUAL(prsm.known_c_term, true)
  TEST_EQUAL(prsm.modifications[0].start, 10)

  // nothing identified: no PrSM
  searched.found = false;
  TEST_EQUAL(AnnotatorMetrics::ms2PrSM(searched, prsm, precursor_mass), false)
}
END_SECTION

// bond_coverage is the annotator's coverage as a fraction, and raw_score its score, of the same annotation.
START_SECTION(bond_coverage_is_a_fraction)
{
  FLASHAnnotatorAlgorithm::PrSM prsm;
  prsm.sequence = kChain;
  prsm.precursor_charge = 2;
  prsm.spectrum_index = 0;

  FLASHAnnotatorAlgorithm by;
  configure(by, {"b", "y"});
  const MSSpectrum spec = envelopesOf(by, fragmentMasses(by, prsm), 2);
  ABORT_IF(spec.empty())

  std::vector<FLASHAnnotatorAlgorithm::PeakArrays> spectra(1);
  FLASHAnnotatorAlgorithm::toCentroidedPeakArrays(spec, spectra[0]);
  std::vector<FLASHAnnotatorAlgorithm::Result> results;
  by.annotate(spectra, std::vector<FLASHAnnotatorAlgorithm::PrSM>(1, prsm), results);
  ABORT_IF(results.size() != 1 || !results[0].annotated)

  const AnnotatorMetrics::Values values = AnnotatorMetrics::annotate(by, spec, prsm);
  TEST_TRUE(results[0].coverage > 1.0)   // a percentage
  TEST_REAL_SIMILAR(values.bond_coverage, results[0].coverage / 100.0)
  TEST_TRUE(values.bond_coverage <= 1.0)
  TEST_REAL_SIMILAR(values.raw_score, results[0].score)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
