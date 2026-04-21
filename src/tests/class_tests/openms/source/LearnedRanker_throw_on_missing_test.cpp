// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Müller $
// $Authors: Tom Müller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/LearnedRanker.h>

#include <stdexcept>

using namespace OpenMS;
using namespace std;

START_TEST(LearnedRanker_throw_on_missing, "LearnedRanker fail-loud on missing model")

START_SECTION([extra] constructor throws on missing model)
{
  bool threw = false;
  try
  {
    LearnedRanker r("/nonexistent/path/to/no_such_model.onnx");
  }
  catch (const std::runtime_error&)
  {
    threw = true;
  }
  TEST_TRUE(threw);
}
END_SECTION

END_TEST
