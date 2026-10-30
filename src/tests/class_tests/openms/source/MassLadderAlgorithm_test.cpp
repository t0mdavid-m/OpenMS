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
#include <OpenMS/CONCEPT/Constants.h>
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

START_SECTION(([EXTRA] searchPrecursorLadders with prepared window starts equals the sweep of each precursor mass))
{
  // a first protein of 66000 single standard residues between non-standard ones (66000 residues of segments, no window), then 300
  // proteins: the first start of a precursor mass in the first chunk lies more than 65535 residues after the chunk's first residue (the
  // long form of a start). The protein is little work, so the first chunk (up to 128 of about equal work) does not end inside it.
  std::mt19937 rng(11);
  const std::string aa = "ACDEFGHIKLMNPQRSTVWY";
  std::vector<std::string> proteins;
  {
    std::string gap;
    for (int k = 0; k < 66000; ++k) gap += "AX";
    proteins.push_back(gap);
  }
  for (int p = 1; p < 301; ++p)
  {
    std::string seq = "M";
    const int length = 100 + (int)(rng() % 300);
    for (int i = 1; i < length; ++i) seq += aa[rng() % aa.size()];
    if (p % 9 == 4) seq[length / 3] = 'X';
    proteins.push_back(seq);
  }
  MassLadderAlgorithm ladders;
  ladders.buildIndex(proteins, {"b", "y"}, true, true);

  // 300 precursor masses (at least 256, so that they are prepared): windows of 15-60 residues of the proteins, some one isotope off
  const auto* rdb = ResidueDB::getInstance();
  const double h2o = Residue::getInternalToFull().getMonoWeight(), b_offset = Residue::getInternalToBIon().getMonoWeight(),
               y_offset = Residue::getInternalToYIon().getMonoWeight();
  auto internal = [&](const std::string& seq, Size from, Size to) {
    double sum = 0;
    for (Size x = from; x < to; ++x) sum += rdb->getResidue(seq[x])->getMonoWeight(Residue::Internal);
    return sum;
  };
  std::vector<double> precursor_masses;
  while (precursor_masses.size() < 300)
  {
    const Size p = 1 + rng() % (proteins.size() - 1);
    const auto& seq = proteins[p];
    const Size length = 15 + rng() % 46, from = rng() % (seq.size() - length);
    if (seq.substr(from, length).find('X') != std::string::npos) continue;
    precursor_masses.push_back(internal(seq, from, from + length) + h2o + (rng() % 5 == 0 ? Constants::C13C12_MASSDIFF_U : 0.0));
  }
  // a spectrum: b and y masses of a few of those windows, and unrelated masses
  std::vector<double> masses;
  for (int k = 0; k < 6; ++k)
  {
    const auto& seq = proteins[1 + k * 7];
    for (Size x = 3; x < 25; x += 2) masses.push_back(internal(seq, 0, x) + b_offset);
    for (Size x = 4; x < 25; x += 3) masses.push_back(internal(seq, seq.size() - x, seq.size()) + y_offset);
  }
  for (int k = 0; k < 40; ++k) masses.push_back(300.0 + (double)(rng() % 900000) / 100.0);
  std::sort(masses.begin(), masses.end());

  std::vector<MassLadderAlgorithm::PrecursorLadderResult> prepared, swept;
  ladders.preparePrecursorWindows(precursor_masses, 10.0);
  ladders.searchPrecursorLadders(masses, precursor_masses, 10.0, 10.0, 5, prepared);
  ladders.clearPrecursorWindows();
  ladders.searchPrecursorLadders(masses, precursor_masses, 10.0, 10.0, 5, swept);
  TEST_EQUAL(prepared.size(), swept.size())
  Size windows = 0, candidates = 0;
  for (Size n = 0; n < swept.size() && n < prepared.size(); ++n)
  {
    TEST_EQUAL(prepared[n].windows, swept[n].windows)
    TEST_EQUAL(prepared[n].candidates.size(), swept[n].candidates.size())
    windows += swept[n].windows;
    candidates += swept[n].candidates.size();
    for (Size c = 0; c < swept[n].candidates.size() && c < prepared[n].candidates.size(); ++c)
    {
      const auto &a = prepared[n].candidates[c], &b = swept[n].candidates[c];
      TEST_EQUAL(a.protein, b.protein)
      TEST_EQUAL(a.rank, b.rank)
      TEST_EQUAL(a.score == b.score, true) // the same windows in the same order, so the same double
      TEST_EQUAL(a.window_start, b.window_start)
      TEST_EQUAL(a.window_end, b.window_end)
      TEST_EQUAL(a.isotope, b.isotope)
      TEST_EQUAL(a.ladders.size(), b.ladders.size())
      for (Size l = 0; l < a.ladders.size() && l < b.ladders.size(); ++l) TEST_EQUAL(a.ladders[l].mass_indices == b.ladders[l].mass_indices, true)
    }
  }
  TEST_EQUAL(windows >= 300, true) // every precursor mass has at least the window it was taken from
  TEST_EQUAL(candidates > 0, true)
}
END_SECTION

END_TEST
