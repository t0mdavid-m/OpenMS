// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>

#include <string>
#include <vector>

using namespace OpenMS;
using PTMSite = FragmentAnalysis::PTMSite;

START_TEST(ProForma, "$Id$")

START_SECTION(toProForma_no_ptms)
{
  std::string result = FragmentAnalysis::toProForma("PEPTIDE", {});
  TEST_STRING_EQUAL(result, "PEPTIDE")
}
END_SECTION

START_SECTION(toProForma_single_localized_ptm)
{
  PTMSite site;
  site.position = 5;
  site.start_position = 5;
  site.end_position = 5;
  site.mass_shift = 79.9663;
  std::string result = FragmentAnalysis::toProForma("PEPTIDE", {site});
  // After residue at 1-based position 5 ('I'): "PEPTI[+79.9663]DE"
  TEST_STRING_EQUAL(result, "PEPTI[+79.9663]DE")
}
END_SECTION

START_SECTION(toProForma_single_ambiguous_ptm)
{
  PTMSite site;
  site.position = 4;
  site.start_position = 3;
  site.end_position = 5;
  site.mass_shift = 19.0523;
  std::string result = FragmentAnalysis::toProForma("PEPTIDE", {site});
  // Region 3-5 = "PTI", wrapped: "PE(PTI)[+19.0523]DE"
  TEST_STRING_EQUAL(result, "PE(PTI)[+19.0523]DE")
}
END_SECTION

START_SECTION(toProForma_negative_mass_shift)
{
  PTMSite site;
  site.position = 5;
  site.start_position = 5;
  site.end_position = 5;
  site.mass_shift = -18.0106;
  std::string result = FragmentAnalysis::toProForma("PEPTIDE", {site});
  TEST_STRING_EQUAL(result, "PEPTI[-18.0106]DE")
}
END_SECTION

START_SECTION(toProForma_ptm_at_first_residue)
{
  PTMSite site;
  site.position = 1;
  site.start_position = 1;
  site.end_position = 1;
  site.mass_shift = 42.0106;
  std::string result = FragmentAnalysis::toProForma("PEPTIDE", {site});
  // After residue at position 1 ('P'): "P[+42.0106]EPTIDE"
  TEST_STRING_EQUAL(result, "P[+42.0106]EPTIDE")
}
END_SECTION

START_SECTION(toProForma_ptm_at_last_residue)
{
  PTMSite site;
  site.position = 7;
  site.start_position = 7;
  site.end_position = 7;
  site.mass_shift = 79.9663;
  std::string result = FragmentAnalysis::toProForma("PEPTIDE", {site});
  // After residue at position 7 ('E', last): "PEPTIDE[+79.9663]"
  TEST_STRING_EQUAL(result, "PEPTIDE[+79.9663]")
}
END_SECTION

START_SECTION(toProForma_multiple_ptms_mixed)
{
  PTMSite localized;
  localized.position = 2;
  localized.start_position = 2;
  localized.end_position = 2;
  localized.mass_shift = 15.9949;

  PTMSite ambiguous;
  ambiguous.position = 6;
  ambiguous.start_position = 5;
  ambiguous.end_position = 7;
  ambiguous.mass_shift = 79.9663;

  std::string result = FragmentAnalysis::toProForma("PEPTIDES", {localized, ambiguous});
  // Right-to-left:
  // 1. Ambiguous 5-7 on "PEPTIDES": "(IDE)" at positions 5-7 -> "PEPT(IDE)[+79.9663]S"
  // 2. Localized at 2: after 'E' -> "PE[+15.9949]PT(IDE)[+79.9663]S"
  TEST_STRING_EQUAL(result, "PE[+15.9949]PT(IDE)[+79.9663]S")
}
END_SECTION

START_SECTION(toProForma_ambiguous_at_start)
{
  PTMSite site;
  site.position = 2;
  site.start_position = 1;
  site.end_position = 3;
  site.mass_shift = 19.0523;
  std::string result = FragmentAnalysis::toProForma("PEPTIDE", {site});
  // Region 1-3 = "PEP", wrapped: "(PEP)[+19.0523]TIDE"
  TEST_STRING_EQUAL(result, "(PEP)[+19.0523]TIDE")
}
END_SECTION

END_TEST
