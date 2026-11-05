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

START_SECTION(([EXTRA] the ladders read the termini with the terminal modifications of the search))
{
  // protein 0 starts with Gln (pyroglutamate), protein 1 with Met-Ser (acetylation after the initiator Met removal); others random
  std::mt19937 rng(11);
  const std::string aa = "ACDEFGHIKLMNPSTVWY"; // no Q, so only protein 0 can carry the pyroglutamate
  std::vector<std::string> proteins;
  for (int p = 0; p < 40; ++p)
  {
    std::string seq = p == 0 ? "Q" : "M";
    if (p == 1) seq += "S";
    while (seq.size() < 120) seq += aa[rng() % aa.size()];
    proteins.push_back(seq);
  }
  const double acetyl = 42.010565, pyro = -17.026549;
  std::vector<MassLadderAlgorithm::TerminalModification> mods(2);
  mods[0].mass = acetyl;
  mods[0].position = 3; // protein N-terminus, any residue
  mods[1].mass = pyro;
  mods[1].position = 3;
  mods[1].residues = "Q";
  MassLadderAlgorithm ladders, plain;
  ladders.buildIndex(proteins, {"b", "y"}, true, true, nullptr, mods, 1);
  plain.buildIndex(proteins, {"b", "y"}, true, true);

  const auto* rdb = ResidueDB::getInstance();
  auto residues = [&](const std::string& seq, Size from, Size to) {
    double sum = 0;
    for (Size x = from; x < to; ++x) sum += rdb->getResidue(seq[x])->getMonoWeight(Residue::Internal);
    return sum;
  };
  const double b_offset = Residue::getInternalToBIon().getMonoWeight(), y_offset = Residue::getInternalToYIon().getMonoWeight();
  const double h2o = Residue::getInternalToFull().getMonoWeight();

  // precursor ladder: protein 1 without its Met, acetylated (residues 1..69), with b ions carrying the acetylation and y ions without
  {
    const auto& seq = proteins[1];
    std::vector<double> masses;
    for (Size x = 4; x < 30; x += 2) masses.push_back(residues(seq, 1, x) + acetyl + b_offset);
    for (Size x = 40; x < 66; x += 3) masses.push_back(residues(seq, x, 70) + y_offset);
    std::sort(masses.begin(), masses.end());
    const double precursor = residues(seq, 1, 70) + acetyl + h2o;
    std::vector<MassLadderAlgorithm::PrecursorLadderResult> with_mods, without;
    ladders.searchPrecursorLadders(masses, {precursor}, 10.0, 10.0, 5, with_mods);
    plain.searchPrecursorLadders(masses, {precursor}, 10.0, 10.0, 5, without);
    TEST_EQUAL(with_mods[0].candidates.empty(), false)
    if (! with_mods[0].candidates.empty())
    {
      const auto& best = with_mods[0].candidates[0];
      TEST_EQUAL(best.protein, 1)
      TEST_EQUAL(best.window_start, 1)
      TEST_EQUAL(best.window_end, 70)
      TEST_EQUAL(best.n_mod, 0)
      TEST_EQUAL(best.c_mod, -1)
    }
    bool plain_found = false; // unmodified, the window of protein 1 does not fit the precursor mass
    for (const auto& c : without[0].candidates) plain_found |= c.protein == 1 && c.window_start == 1 && c.window_end == 70;
    TEST_EQUAL(plain_found, false)
    TEST_EQUAL(with_mods[0].windows > without[0].windows, true) // the acetylated windows are hypotheses too
  }

  // sliding ladder: b ions of protein 0 from its pyroglutamate N-terminus
  {
    const auto& seq = proteins[0];
    std::vector<double> masses;
    for (Size x = 3; x < 40; x += 2) masses.push_back(residues(seq, 0, x) + pyro + b_offset);
    for (int k = 0; k < 10; ++k) masses.push_back(500.0 + (double)(rng() % 500000) / 100.0);
    std::sort(masses.begin(), masses.end());
    std::vector<MassLadderAlgorithm::Candidate> with_mods, without;
    ladders.searchSlidingLadders(masses, 10.0, 5, {}, with_mods);
    plain.searchSlidingLadders(masses, 10.0, 5, {}, without);
    TEST_EQUAL(with_mods.empty(), false)
    if (! with_mods.empty())
    {
      TEST_EQUAL(with_mods[0].protein, 0)
      TEST_EQUAL(with_mods[0].n_mod, 1)
      TEST_EQUAL(with_mods[0].ladders.empty(), false)
      if (! with_mods[0].ladders.empty()) TEST_EQUAL(with_mods[0].ladders[0].mass_indices.size() >= 15, true)
    }
    bool plain_found = false;
    for (const auto& c : without) plain_found |= c.protein == 0;
    TEST_EQUAL(plain_found, false)
  }

  // terminal ladder: the same masses, read from the pyroglutamate N-terminus
  {
    const auto& seq = proteins[0];
    std::vector<double> masses;
    for (Size x = 3; x < 40; x += 2) masses.push_back(residues(seq, 0, x) + pyro + b_offset);
    std::sort(masses.begin(), masses.end());
    std::vector<MassLadderAlgorithm::Candidate> with_mods, without;
    ladders.searchTerminalLadders(masses, 10.0, 5, {}, with_mods);
    plain.searchTerminalLadders(masses, 10.0, 5, {}, without);
    TEST_EQUAL(with_mods.empty(), false)
    if (! with_mods.empty())
    {
      TEST_EQUAL(with_mods[0].protein, 0)
      TEST_EQUAL(with_mods[0].n_mod, 1)
    }
    bool plain_found = false;
    for (const auto& c : without) plain_found |= c.protein == 0;
    TEST_EQUAL(plain_found, false)
  }
}
END_SECTION

END_TEST
