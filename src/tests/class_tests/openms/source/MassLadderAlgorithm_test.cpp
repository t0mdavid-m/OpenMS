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
#include <OpenMS/ANALYSIS/TOPDOWN/MassLadderAlgorithm.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <random>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassLadderAlgorithm, "$Id$")

START_SECTION(([EXTRA] searchSlidingLadders for many spectra equals the search of each spectrum alone))
{
  // a database of 60 proteins, some with a non-standard residue (X) that splits them into two runs of standard residues
  std::mt19937 rng(7);
  const std::string aa = "ACDEFGHIKLMNPQRSTVWY";
  std::vector<std::string> proteins;
  for (int p = 0; p < 60; ++p)
  {
    std::string seq = "M";
    const int length = 80 + (int)(rng() % 220);
    for (int i = 1; i < length; ++i) seq += aa[rng() % aa.size()];
    if (p % 7 == 3) seq[length / 2] = 'X';
    proteins.push_back(seq);
  }
  MassLadderAlgorithm ladders;
  ladders.buildIndex(proteins, {"b", "y"}, true, true);

  // 70 spectra (more than one batch of 64): b and y ladders from an internal start or end of a protein, plus unrelated masses
  const auto* rdb = ResidueDB::getInstance();
  const double b_offset = Residue::getInternalToBIon().getMonoWeight(), y_offset = Residue::getInternalToYIon().getMonoWeight();
  std::vector<std::vector<double>> masses(70);
  std::vector<std::set<Size>> excluded(70);
  for (Size s = 0; s < masses.size(); ++s)
  {
    const auto& seq = proteins[(s * 13) % proteins.size()];
    const Size anchor = 5 + s % 20;
    double sum = 0;
    for (Size x = anchor; x < anchor + 12 && x < seq.size(); ++x)
    {
      if (seq[x] == 'X') break;
      sum += rdb->getResidue(seq[x])->getMonoWeight(Residue::Internal);
      if (x % 2 == 0) masses[s].push_back(sum + b_offset);
    }
    sum = 0;
    for (Size x = seq.size() - 1 - s % 15; x > seq.size() - 14 - s % 15; --x)
    {
      if (seq[x] == 'X') break;
      sum += rdb->getResidue(seq[x])->getMonoWeight(Residue::Internal);
      if (x % 3 != 0) masses[s].push_back(sum + y_offset);
    }
    for (int k = 0; k < 15; ++k) masses[s].push_back(300.0 + (double)(rng() % 900000) / 100.0);
    std::sort(masses[s].begin(), masses[s].end());
    if (s % 3 == 1) excluded[s] = {(s * 13) % proteins.size(), 2, 5}; // also the protein of the ladder
    if (s % 11 == 4) masses[s].clear();                              // a spectrum without masses
  }

  std::vector<const std::vector<double>*> masses_of;
  std::vector<const std::set<Size>*> excluded_of;
  for (Size s = 0; s < masses.size(); ++s)
  {
    masses_of.push_back(&masses[s]);
    excluded_of.push_back(&excluded[s]);
  }
  std::vector<std::vector<MassLadderAlgorithm::Candidate>> batched;
  ladders.searchSlidingLadders(masses_of, 10.0, 20, excluded_of, batched);
  TEST_EQUAL(batched.size(), masses.size())

  Size with_candidates = 0;
  for (Size s = 0; s < masses.size(); ++s)
  {
    std::vector<MassLadderAlgorithm::Candidate> alone;
    ladders.searchSlidingLadders(masses[s], 10.0, 20, excluded[s], alone);
    TEST_EQUAL(batched[s].size(), alone.size())
    if (batched[s].size() != alone.size()) continue;
    with_candidates += ! alone.empty();
    for (Size c = 0; c < alone.size(); ++c)
    {
      TEST_EQUAL(batched[s][c].protein, alone[c].protein)
      TEST_EQUAL(batched[s][c].rank, alone[c].rank)
      TEST_EQUAL(batched[s][c].score == alone[c].score, true) // the same arithmetic, so the same double
      TEST_EQUAL(batched[s][c].ladders.size(), alone[c].ladders.size())
      for (Size l = 0; l < alone[c].ladders.size() && l < batched[s][c].ladders.size(); ++l)
      {
        TEST_EQUAL(batched[s][c].ladders[l].prefix, alone[c].ladders[l].prefix)
        TEST_EQUAL(batched[s][c].ladders[l].start, alone[c].ladders[l].start)
        TEST_EQUAL(batched[s][c].ladders[l].end, alone[c].ladders[l].end)
        TEST_EQUAL(batched[s][c].ladders[l].mass_indices == alone[c].ladders[l].mass_indices, true)
      }
    }
    if (! alone.empty() && excluded[s].empty()) TEST_EQUAL(alone[0].protein, (s * 13) % proteins.size()) // the planted ladder wins
    for (const auto& cand : alone) TEST_EQUAL(excluded[s].count(cand.protein), 0)
  }
  TEST_EQUAL(with_candidates > 30, true)
}
END_SECTION

END_TEST
