// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/CONCEPT/Constants.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
  /// A centroided spectrum made of the predicted envelopes of all fragments of @p prsm at charges 1..max_charge.
  /// Every envelope is moved by @p offset isotope spacings and its m/z scaled by 1 + ppm * 1e-6.
  FLASHAnnotatorAlgorithm::PeakArrays syntheticSpectrum(FLASHAnnotatorAlgorithm& algo, const FLASHAnnotatorAlgorithm::PrSM& prsm, int max_charge,
                                                        double offset, double ppm)
  {
    std::vector<FLASHAnnotatorAlgorithm::Fragment> fragments;
    algo.getFragments(prsm, fragments);
    MSSpectrum spec;
    std::vector<double> mzs, intensities;
    for (Size i = 0; i < fragments.size(); i++)
    {
      for (int z = 1; z <= std::min(max_charge, fragments[i].max_charge); z++)
      {
        algo.getIsotopeEnvelope(fragments[i].mono_mass, z, mzs, intensities);
        for (Size j = 0; j < mzs.size(); j++)
        {
          double mz = (mzs[j] + offset * Constants::C13C12_MASSDIFF_U / z) * (1.0 + ppm * 1e-6);
          spec.emplace_back(mz, (float)(intensities[j] * (1000.0 + 10.0 * i)));
        }
      }
    }
    spec.sortByPosition();
    FLASHAnnotatorAlgorithm::PeakArrays peaks;
    FLASHAnnotatorAlgorithm::toPeakArrays(spec, peaks);
    return peaks;
  }
} // namespace

START_TEST(FLASHAnnotatorAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FLASHAnnotatorAlgorithm* ptr = nullptr;
FLASHAnnotatorAlgorithm* null_ptr = nullptr;
START_SECTION(FLASHAnnotatorAlgorithm())
{
  ptr = new FLASHAnnotatorAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FLASHAnnotatorAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION((static double binomialTailScore(int n, int k, double p)))
{
  TOLERANCE_ABSOLUTE(1e-6)
  TEST_REAL_SIMILAR(FLASHAnnotatorAlgorithm::binomialTailScore(10, 0, 0.5), 0.0)
  TEST_REAL_SIMILAR(FLASHAnnotatorAlgorithm::binomialTailScore(1, 1, 0.5), 0.30103)
  TEST_REAL_SIMILAR(FLASHAnnotatorAlgorithm::binomialTailScore(10, 10, 0.5), 3.0103)
  TEST_REAL_SIMILAR(FLASHAnnotatorAlgorithm::binomialTailScore(20, 5, 0.1), 1.364773)
  TEST_REAL_SIMILAR(FLASHAnnotatorAlgorithm::binomialTailScore(200, 30, 0.02), 16.796925)
  TEST_REAL_SIMILAR(FLASHAnnotatorAlgorithm::binomialTailScore(50, 3, 0.2), 0.000559)
  // no underflow far in the tail
  TEST_EQUAL(FLASHAnnotatorAlgorithm::binomialTailScore(2000, 1500, 0.001) > 300, true)
}
END_SECTION

START_SECTION((static double getProteoformMass(const String& sequence, const std::vector<Modification>& modifications)))
{
  TOLERANCE_ABSOLUTE(1e-3)
  std::vector<FLASHAnnotatorAlgorithm::Modification> mods;
  TEST_REAL_SIMILAR(FLASHAnnotatorAlgorithm::getProteoformMass("PEPTIDE", mods), 799.35996)
  mods.push_back({15.9949, 1, 1});
  TEST_REAL_SIMILAR(FLASHAnnotatorAlgorithm::getProteoformMass("PEPTIDE", mods), 815.35486)
  TEST_REAL_SIMILAR(FLASHAnnotatorAlgorithm::getProteoformMass("PEP#TIDE", mods), -1.0)
}
END_SECTION

const String sequence = "AEIYNKDGNKVDLYGKAVGLHYFSKGNGENSYGGNGDMTY"; // 40 residues
FLASHAnnotatorAlgorithm::PrSM prsm;
prsm.sequence = sequence;
prsm.precursor_charge = 3;
prsm.spectrum_index = 0;

START_SECTION((bool getFragments(const PrSM& prsm, std::vector<Fragment>& fragments) const))
{
  FLASHAnnotatorAlgorithm algo;
  std::vector<FLASHAnnotatorAlgorithm::Fragment> fragments;
  TEST_EQUAL(algo.getFragments(prsm, fragments), true)
  // b and y of lengths 1..39 at least 300 Da
  TEST_EQUAL(fragments.size() > 70 && fragments.size() <= 78, true)
  TOLERANCE_ABSOLUTE(1e-3)
  for (const auto& f : fragments)
  {
    TEST_EQUAL(f.max_charge >= 1 && f.max_charge <= 3, true)
    TEST_EQUAL(algo.getIonTypes()[f.ion_type] == "b" ? f.bond == f.length : f.bond == 40 - f.length, true)
    if (algo.getIonTypes()[f.ion_type] == "b" && f.length == 4) { TEST_REAL_SIMILAR(f.mono_mass, 71.03711 + 129.04259 + 113.08406 + 163.06333) }
    if (algo.getIonTypes()[f.ion_type] == "y" && f.length == 3) { TEST_REAL_SIMILAR(f.mono_mass, 131.04049 + 101.04768 + 163.06333 + 18.01056) }
  }

  // a modification on residues 3-4 (0-based): b4 and y36 split it, longer prefixes carry its mass
  FLASHAnnotatorAlgorithm::PrSM modified = prsm;
  modified.modifications.push_back({42.0106, 3, 4});
  std::vector<FLASHAnnotatorAlgorithm::Fragment> mod_fragments;
  algo.getFragments(modified, mod_fragments);
  TEST_EQUAL(mod_fragments.size(), fragments.size() - 2) // b4 and y36 split the range
  for (const auto& f : mod_fragments)
  {
    if (algo.getIonTypes()[f.ion_type] == "b" && f.length == 4) { TEST_EQUAL(true, false) }
    if (algo.getIonTypes()[f.ion_type] == "b" && f.length == 6) { TEST_REAL_SIMILAR(f.mono_mass, 71.03711 + 129.04259 + 113.08406 + 163.06333 + 114.04293 + 128.09496 + 42.0106) }
  }

  // undetermined start: no prefix fragments
  FLASHAnnotatorAlgorithm::PrSM no_start = prsm;
  no_start.known_n_term = false;
  algo.getFragments(no_start, mod_fragments);
  for (const auto& f : mod_fragments)
  {
    TEST_EQUAL(algo.getIonTypes()[f.ion_type], "y")
  }

  FLASHAnnotatorAlgorithm::PrSM unknown = prsm;
  unknown.sequence = "AEIYNKDG#NKVDLYG";
  TEST_EQUAL(algo.getFragments(unknown, mod_fragments), false)
}
END_SECTION

START_SECTION((void getIsotopeEnvelope(double mono_mass, int charge, std::vector<double>& mzs, std::vector<double>& intensities)))
{
  FLASHAnnotatorAlgorithm algo;
  std::vector<double> mzs, intensities;
  algo.getIsotopeEnvelope(5000.0, 5, mzs, intensities);
  TEST_EQUAL(mzs.size() >= 5, true)
  TEST_EQUAL(mzs.size(), intensities.size())
  double top = *std::max_element(intensities.begin(), intensities.end());
  TEST_REAL_SIMILAR(top, 1.0)
  for (Size j = 0; j < intensities.size(); j++)
  {
    TEST_EQUAL(intensities[j] >= 0.1, true)
    if (j > 0) { TEST_REAL_SIMILAR(mzs[j] - mzs[j - 1], Constants::C13C12_MASSDIFF_U / 5) }
  }
}
END_SECTION

START_SECTION((void annotate(const std::vector<PeakArrays>& spectra, const std::vector<PrSM>& prsms, std::vector<Result>& results)))
{
  FLASHAnnotatorAlgorithm algo;
  std::vector<FLASHAnnotatorAlgorithm::Result> results;

  // the algorithm's own envelopes at charges 1-3: every bond is covered
  std::vector<FLASHAnnotatorAlgorithm::PeakArrays> spectra {syntheticSpectrum(algo, prsm, 3, 0.0, 0.0)};
  algo.annotate(spectra, {prsm}, results);
  TEST_EQUAL(results.size(), 1)
  TEST_EQUAL(results[0].annotated, true)
  TEST_EQUAL(results[0].covered_bonds, 39)
  TEST_REAL_SIMILAR(results[0].coverage, 100.0)
  TEST_EQUAL(results[0].matched_fragments > 60, true)
  TEST_EQUAL(results[0].matched_fragments <= results[0].fragment_hypotheses, true)
  TEST_EQUAL(results[0].chance_rate > 0 && results[0].chance_rate < 0.05, true)
  TEST_EQUAL(results[0].score > 20, true)
  TEST_EQUAL(results[0].llr > 0, true)
  TOLERANCE_ABSOLUTE(0.01)
  TEST_REAL_SIMILAR(results[0].calibration_ppm, 0.0)
  for (const auto& e : results[0].envelopes)
  {
    TEST_EQUAL(e.charge >= 1 && e.charge <= 3, true)
    TEST_EQUAL(e.matched_peaks, e.predicted_peaks)
    TEST_EQUAL(e.cosine > 0.99 && e.pearson > 0.99, true)
  }

  // the same envelopes moved off the isotope grid by half an isotope spacing: nothing matches
  spectra[0] = syntheticSpectrum(algo, prsm, 3, 0.5, 0.0);
  algo.annotate(spectra, {prsm}, results);
  TEST_EQUAL(results[0].annotated, true)
  TEST_EQUAL(results[0].matched_fragments, 0)
  TEST_EQUAL(results[0].covered_bonds, 0)
  TEST_EQUAL(results[0].envelopes.size(), 0)
  TEST_REAL_SIMILAR(results[0].score, 0.0)

  // 15 ppm off: lost at 10 ppm without calibration, recovered with it
  spectra[0] = syntheticSpectrum(algo, prsm, 3, 0.0, 15.0);
  Param param = algo.getParameters();
  param.setValue("calibrate", "false");
  algo.setParameters(param);
  algo.annotate(spectra, {prsm}, results);
  TEST_EQUAL(results[0].matched_fragments, 0)
  param.setValue("calibrate", "true");
  algo.setParameters(param);
  algo.annotate(spectra, {prsm}, results);
  TEST_EQUAL(results[0].covered_bonds, 39)
  TEST_REAL_SIMILAR(results[0].calibration_ppm, 15.0)

  // no decoy shifts: no chance rate and no scores
  param.setValue("decoy_shifts", std::vector<double>());
  algo.setParameters(param);
  algo.annotate(spectra, {prsm}, results);
  TEST_EQUAL(results[0].covered_bonds, 39)
  TEST_REAL_SIMILAR(results[0].chance_rate, -1.0)
  TEST_REAL_SIMILAR(results[0].score, -1.0)

  // a PrSM without its spectrum is not annotated
  FLASHAnnotatorAlgorithm::PrSM missing = prsm;
  missing.spectrum_index = -1;
  algo.annotate(spectra, {missing}, results);
  TEST_EQUAL(results[0].annotated, false)
  TEST_EQUAL(results[0].covered_bonds, -1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
