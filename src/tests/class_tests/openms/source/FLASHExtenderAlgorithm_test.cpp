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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
///////////////////////////

using namespace OpenMS;

// the truncation penalty FLASHTnT's precursor re-selection allows for
START_TEST(FLASHExtenderAlgorithm, "$Id$")

START_SECTION((truncation_penalty))
{
  TEST_EQUAL(FLASHExtenderAlgorithm::truncation_penalty, 19)
}
END_SECTION

END_TEST
