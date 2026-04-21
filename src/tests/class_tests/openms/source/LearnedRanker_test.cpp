// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Müller $
// $Authors: Tom Müller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/LearnedRanker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <cmath>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(LearnedRanker, "FLASHIda ONNX-backed precursor scorer")

LearnedRanker* ranker = nullptr;

START_SECTION(LearnedRanker(const std::string& model_path))
{
  const string model_path = OPENMS_GET_TEST_DATA_PATH("LearnedRanker_test_policy.onnx");
  ranker = new LearnedRanker(model_path);
  TEST_NOT_EQUAL(ranker, nullptr);
}
END_SECTION

START_SECTION(std::vector<float> score(const std::vector<PeakGroup>& candidates, const LearnedRankerGlobalContext& ctx))
{
  // Build a small synthetic candidate set
  std::vector<PeakGroup> candidates;
  for (int i = 0; i < 5; ++i)
  {
    PeakGroup pg(4, 40, true);
    pg.setMonoisotopicMass(5000.0 + i * 100.0);
    pg.setQscore(0.5 + 0.1 * i);
    pg.setIsotopeCosine(0.9f);
    pg.setChargeScore(0.8f);
    pg.setSNR(10.0f);
    pg.setAvgPPMError(2.0f);
    pg.setRepAbsCharge(10);
    pg.setAbsChargeRange(8, 12);
    pg.setChargeIsotopeCosine(10, 0.95f);
    candidates.push_back(pg);
  }

  LearnedRankerGlobalContext ctx;
  ctx.rt = 10.0f;
  ctx.faims_cv = 0.0f;
  ctx.queue_depth = 3.0f;

  auto scores = ranker->score(candidates, ctx);
  TEST_EQUAL(scores.size(), 5);
  for (float s : scores) { TEST_TRUE(std::isfinite(s)); }
}
END_SECTION

delete ranker;

END_TEST
