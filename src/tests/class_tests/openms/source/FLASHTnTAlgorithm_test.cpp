// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <cstdint>
#include <set>
///////////////////////////

using namespace OpenMS;

// The label-blind rules of FLASHTnT that can be checked without a search: the entry order that breaks candidate-cut ties, the protein
// evidence rule (option protein_evidence) and the blind modification count it reads, and the defaults of the options
START_TEST(FLASHTnTAlgorithm, "$Id$")

START_SECTION((static bool proteinEvidenceSupported(Size proteoforms, bool strong)))
{
  // two distinct proteoforms, or one strong proteoform
  TEST_EQUAL(FLASHTnTAlgorithm::proteinEvidenceSupported(1, false), false)
  TEST_EQUAL(FLASHTnTAlgorithm::proteinEvidenceSupported(2, false), true)
  TEST_EQUAL(FLASHTnTAlgorithm::proteinEvidenceSupported(1, true), true)
  TEST_EQUAL(FLASHTnTAlgorithm::proteinEvidenceSupported(0, false), false)
}
END_SECTION

START_SECTION((static Size blindModCount(const ProteinHit& hit)))
{
  ProteinHit h;
  h.setMetaValue("Modifications", DoubleList {15.99, 42.01});
  TEST_EQUAL(FLASHTnTAlgorithm::blindModCount(h), 2) // without the kinds, every modification counts as blind
  h.setMetaValue("ModificationVariable", IntList {1, 0}); // one variable, one blind
  TEST_EQUAL(FLASHTnTAlgorithm::blindModCount(h), 1)
  h.setMetaValue("ModificationVariable", IntList {1, 2}); // one variable, one fixed by the database entry
  TEST_EQUAL(FLASHTnTAlgorithm::blindModCount(h), 0)
  h.setMetaValue("BlindModCount", 1); // the search's own count wins
  TEST_EQUAL(FLASHTnTAlgorithm::blindModCount(h), 1)
}
END_SECTION

START_SECTION((static std::uint64_t entryOrder(std::uint64_t index)))
{
  // deterministic, collision-free on the first entries, and not the database order: about half the consecutive pairs ascend
  std::set<std::uint64_t> seen;
  Size ascending = 0;
  for (std::uint64_t i = 0; i < 1000; i++)
  {
    TEST_EQUAL(FLASHTnTAlgorithm::entryOrder(i), FLASHTnTAlgorithm::entryOrder(i))
    seen.insert(FLASHTnTAlgorithm::entryOrder(i));
    if (i > 0 && FLASHTnTAlgorithm::entryOrder(i) > FLASHTnTAlgorithm::entryOrder(i - 1)) ascending++;
  }
  TEST_EQUAL(seen.size(), 1000)
  TEST_EQUAL(ascending > 400 && ascending < 600, true)
  TEST_EQUAL(FLASHTnTAlgorithm::entryOrder(0) == 0, false)
}
END_SECTION

START_SECTION((FLASHTnTAlgorithm()))
{
  // the defaults of the FDR and mass options, and the values each of them takes
  const Param p = FLASHTnTAlgorithm().getDefaults();
  TEST_STRING_EQUAL(p.getValue("terminal_ladder_class").toString(), "true")
  TEST_EQUAL(p.getValidStrings("terminal_ladder_class").size(), 2)
  TEST_EQUAL(p.hasTag("terminal_ladder_class", "advanced"), true)
  TEST_STRING_EQUAL(p.getValue("protein_evidence").toString(), "proteoforms")
  TEST_EQUAL(p.getValidStrings("protein_evidence").size(), 2)
  TEST_REAL_SIMILAR((double)p.getValue("protein_evidence_min_mass"), 3000.0)
}
END_SECTION

END_TEST
