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
#include <cstdint>
#include <set>
///////////////////////////

using namespace OpenMS;

// The label-blind rules of FLASHTnT that can be checked without a search: the entry order that breaks candidate-cut ties, and the
// defaults of the options
START_TEST(FLASHTnTAlgorithm, "$Id$")

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
}
END_SECTION

END_TEST
