// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorSelection.h>

#include "FLASHIda_TestAccess.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

using namespace OpenMS;

// FASTA tag-based targeting adds the masses of a detected protein's proteoforms to the inclusion list
// for rt_window seconds (PrecursorSelection::addDynamicTargets_). These sections pin WHAT an
// inclusion-mode survey then matches against: the reference targets PLUS the tag-expanded ones.
//
// The tag-expanded rows used to be appended to inclusion_targets_, which filterAndRank reads as "a TSV
// inclusion list was loaded" -- so with ida.log reference targets (files.target_logs), the first tag
// detection of a run switched every later survey to the TSV branch, and every ida.log target was gone
// for the rest of the run. Needs no fixture: the reference lists are written by the test itself.

namespace
{
  // Windows temp paths carry backslashes, which are escapes in JSON; every reader here accepts '/'.
  std::string jsonPath(std::string p)
  {
    std::replace(p.begin(), p.end(), '\\', '/');
    return p;
  }

  std::string cfgJson(const std::string& files)
  {
    return std::string(R"({
      "deconvolution": { "tol": [10, 10, 10] },
      "precursor_selection": { "rank_by": "qscore", "max_precursors": 3, "targeting": "inclusion",
                               "strict_inclusion": true, "rt_window": 180 },
      "characterization": { "mode": "off" },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "resolution": 120000 },
        "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29 }
      },
      "files": )") + files + R"(
    })";
  }

  bool has(const std::vector<double>& v, double m)
  {
    return std::any_of(v.begin(), v.end(), [m](double x) { return std::abs(x - m) < 1e-6; });
  }

  // Two ida.log targets, in the engine's own ida.log layout (RT in minutes on the MS1 line).
  void writeIdaLog(const std::string& path)
  {
    std::ofstream o(path);
    o << "MS1 Scan# 3 RT 0.5000 (Access ID !!$) - 1 targets\n"
      << "Mass=3897.6156\tZ=7\tScore=0.34746\tWindow=[557.5548-558.6361]\tPrecursorIntensity=12271.77441\n"
      << "AllMass=3897.6156\n"
      << "MS1 Scan# 9 RT 0.9000 (Access ID !!*) - 1 targets\n"
      << "Mass=10261.3901\tZ=11\tScore=0.81000\tWindow=[933.8000-934.9000]\tPrecursorIntensity=99999.00000\n"
      << "AllMass=10261.3901\n";
  }
}

START_TEST(FLASHIda_DynamicTargets, "$Id$")

START_SECTION(tag_expanded_targets_add_to_ida_log_targets)
{
  std::string log;
  NEW_TMP_FILE(log)
  writeIdaLog(log);
  Config cfg(cfgJson(R"({ "target_logs": [")" + jsonPath(log) + R"("] })"));
  Deconvolution deconv(cfg, cfg.toleranceList());
  PrecursorSelection sel(cfg, deconv);

  const std::vector<double> before = FLASHIdaTestAccess::inclusionTargetMasses(sel, 60.0);
  TEST_EQUAL(before.size(), 2)
  TEST_TRUE(has(before, 3897.6156))
  TEST_TRUE(has(before, 10261.3901))

  // A tag detection at rt 50 s: two proteoform masses, active 50 .. 50 + rt_window s.
  FLASHIdaTestAccess::addDynamicTargets(sel, {12345.6789, 12387.6895}, 50.0, 100);

  const std::vector<double> during = FLASHIdaTestAccess::inclusionTargetMasses(sel, 60.0);
  // ISSUE(bug 2): pre-fix this was 2 -- the two tag-expanded masses ALONE; both ida.log targets gone.
  TEST_EQUAL(during.size(), 4)
  TEST_TRUE(has(during, 3897.6156))
  TEST_TRUE(has(during, 10261.3901))
  TEST_TRUE(has(during, 12345.6789))
  TEST_TRUE(has(during, 12387.6895))
  TEST_TRUE(std::is_sorted(during.begin(), during.end()))   // the selection loop binary-searches it

  // Past the tag-expanded rows' window the reference list is back to itself -- and is still there.
  // Pre-fix: EMPTY, which under strict inclusion means no MS2 for the rest of the run.
  const std::vector<double> after = FLASHIdaTestAccess::inclusionTargetMasses(sel, 50.0 + 180.0 + 1.0);
  TEST_EQUAL(after.size(), 2)
  TEST_TRUE(has(after, 3897.6156))
  TEST_TRUE(has(after, 10261.3901))
}
END_SECTION

START_SECTION(tag_expanded_targets_add_to_a_tsv_inclusion_list)
{
  // The TSV path is the one that always worked; pinned so the fix cannot have moved it.
  std::string tsv;
  NEW_TMP_FILE(tsv)
  {
    std::ofstream o(tsv);
    o << "mass\tcharge\trt_start\trt_end\tpriority\n"
      << "12351.3000\t-1\t0.5\t10.0\t5\n"      // active 30 .. 600 s
      << "8565.7600\t-1\t20.0\t30.0\t1\n";     // active 1200 .. 1800 s
  }
  Config cfg(cfgJson(R"({ "inclusion_list": ")" + jsonPath(tsv) + R"(" })"));
  Deconvolution deconv(cfg, cfg.toleranceList());
  PrecursorSelection sel(cfg, deconv);

  TEST_EQUAL(FLASHIdaTestAccess::inclusionTargetMasses(sel, 60.0).size(), 1)
  FLASHIdaTestAccess::addDynamicTargets(sel, {5000.5}, 50.0, 100);
  const std::vector<double> during = FLASHIdaTestAccess::inclusionTargetMasses(sel, 60.0);
  TEST_EQUAL(during.size(), 2)
  TEST_TRUE(has(during, 12351.3))
  TEST_TRUE(has(during, 5000.5))
  // The row outside its own RT window stays out, tag-expanded rows or not.
  TEST_TRUE(!has(during, 8565.76))
  TEST_EQUAL(FLASHIdaTestAccess::inclusionTargetMasses(sel, 1500.0).size(), 1)
  TEST_TRUE(has(FLASHIdaTestAccess::inclusionTargetMasses(sel, 1500.0), 8565.76))
}
END_SECTION

START_SECTION(tag_expanded_targets_alone_without_a_reference_list)
{
  // No ida.log and no TSV: the tag-expanded rows are the whole list while they are active.
  Config cfg(cfgJson(R"({ "target_logs": [] })"));
  Deconvolution deconv(cfg, cfg.toleranceList());
  PrecursorSelection sel(cfg, deconv);

  TEST_EQUAL(FLASHIdaTestAccess::inclusionTargetMasses(sel, 60.0).size(), 0)
  FLASHIdaTestAccess::addDynamicTargets(sel, {7000.25, 6000.125}, 50.0, 100);
  const std::vector<double> during = FLASHIdaTestAccess::inclusionTargetMasses(sel, 60.0);
  TEST_EQUAL(during.size(), 2)
  TEST_REAL_SIMILAR(during[0], 6000.125)
  TEST_REAL_SIMILAR(during[1], 7000.25)
  TEST_EQUAL(FLASHIdaTestAccess::inclusionTargetMasses(sel, 49.0).size(), 0)    // before rt_start
}
END_SECTION

END_TEST
