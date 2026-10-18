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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHGappedTaggerAlgorithm.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
  double residueSum(const std::string& seq, Size from, Size to)
  {
    double m = 0;
    for (Size i = from; i < to; ++i) m += ResidueDB::getInstance()->getResidue(seq[i])->getMonoWeight(Residue::Internal);
    return m;
  }

  /// a spectrum of the given masses (Qscore 0.9) plus masses no residue sum connects
  DeconvolvedSpectrum spectrum(const std::vector<double>& masses)
  {
    DeconvolvedSpectrum dspec(1);
    std::vector<double> all = masses;
    for (double m : {1234.567, 2345.678, 3456.789}) all.push_back(m);
    for (double m : all)
    {
      PeakGroup pg(1, 10, true);
      pg.setMonoisotopicMass(m);
      pg.setQscore(0.9);
      dspec.push_back(pg);
    }
    dspec.sort();
    return dspec;
  }
} // namespace

START_TEST(FLASHGappedTaggerAlgorithm, "$Id$")

FLASHGappedTaggerAlgorithm* ptr = nullptr;
FLASHGappedTaggerAlgorithm* null_ptr = nullptr;
START_SECTION(FLASHGappedTaggerAlgorithm())
  ptr = new FLASHGappedTaggerAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
END_SECTION

START_SECTION(~FLASHGappedTaggerAlgorithm())
  delete ptr;
END_SECTION

const std::vector<std::string> proteins = {"MSTNPKPQRKTKRNTNRRPQDVKFPGGGQIVGGVYLLPRRGPRLGVRATRKTSERSQPRGRRQPIPKARRPEGRTWAQPGYPWPLYGNEGCGWAGWLLSPRGSRPSWGPTDPRRRSRNLGKVIDTLTCGFADLMGYIPLVGAPLGGAARALAHGVRVLEDGVNYATGNLPGCSFSIFLLALLSCLTVPASA",
                                           "MKKTAIAIAVALAGFATVAQAAPKDNTWYTGAKLGWSQYHDTGFINNNGPTHENQLGAGAFGGYQVNPYVGFEMGYDWLGRMPYKGSVENGAYKAQGVQLTAKLGYPITDDLDIYTRLGGMVWRADTKSNVYGKNHDTGVSPVFAGGVEYAITPEIATRLEYQWTNNIGDAHTIGTRPDNGMLSLGVSYRFGQGEAAPVVAPAPAPAPE",
                                           "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"};

START_SECTION(void buildIndex(const std::vector<std::string>& protein_seqs))
  FLASHGappedTaggerAlgorithm tagger;
  tagger.buildIndex(proteins);
  TEST_EQUAL(tagger.getIndexBytes() > 0, true)
END_SECTION

START_SECTION(void search(const DeconvolvedSpectrum& dspec, double ppm, std::vector<Landing>& landings) const)
{
  FLASHGappedTaggerAlgorithm tagger;
  tagger.buildIndex(proteins);
  const std::string& seq = proteins[1];

  // b-like masses at residues 30, 32, 33, 36, 38 of protein 1: gaps of 2, 1, 3, 2 residues, 4 ppm off
  const std::vector<Size> cuts = {30, 32, 33, 36, 38};
  std::vector<double> masses;
  for (Size c : cuts) masses.push_back((residueSum(seq, 0, c) + Constants::PROTON_MASS_U) * (1 + 4e-6));
  std::vector<FLASHGappedTaggerAlgorithm::Landing> landings;
  tagger.search(spectrum(masses), 5.0, landings);
  bool found = false;
  for (const auto& l : landings)
  {
    if (l.protein != 1 || l.start != 30 || ! l.prefix) continue;
    found = true;
    TEST_EQUAL(l.end, 38)
    TEST_EQUAL(l.masses.size(), 5)
    TEST_REAL_SIMILAR(l.flank_mass, residueSum(seq, 0, 30))
  }
  TEST_EQUAL(found, true)

  // y-like masses read from the C-terminus: the landing starts at the residue after the highest-mass cut
  masses.clear();
  for (Size c : cuts) masses.push_back(residueSum(seq, c, seq.size()) + EmpiricalFormula("H2O").getMonoWeight());
  tagger.search(spectrum(masses), 5.0, landings);
  found = false;
  for (const auto& l : landings)
  {
    if (l.protein != 1 || l.start != 30 || l.prefix) continue;
    found = true;
    TEST_EQUAL(l.end, 38)
    TEST_REAL_SIMILAR(l.flank_mass, residueSum(seq, 38, seq.size()))
  }
  TEST_EQUAL(found, true)

  // the same masses alternately 0.3 Da too heavy and too light: no gap is a residue sum, nothing lands there
  masses.clear();
  for (Size i = 0; i < cuts.size(); ++i) masses.push_back(residueSum(seq, 0, cuts[i]) + Constants::PROTON_MASS_U + (i % 2 ? 0.3 : -0.3));
  tagger.search(spectrum(masses), 5.0, landings);
  found = false;
  for (const auto& l : landings) found |= l.protein == 1 && l.start == 30;
  TEST_EQUAL(found, false)
}
END_SECTION

END_TEST
