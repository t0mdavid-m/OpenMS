// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <cmath>
#include <string>
#include <vector>

using namespace OpenMS;

START_TEST(MS3FragmentMatcher, "$Id$")

START_SECTION(getMS3IonTypes)
{
  // b-precursor: a, b (same direction) + yb, ya (cross-direction)
  auto b_types = MS3FragmentMatcher::getMS3IonTypes('b');
  TEST_EQUAL(b_types.size(), 4)
  TEST_EQUAL(b_types[0], "a")
  TEST_EQUAL(b_types[1], "b")
  TEST_EQUAL(b_types[2], "yb")
  TEST_EQUAL(b_types[3], "ya")

  // a-precursor: same as b-precursor
  auto a_types = MS3FragmentMatcher::getMS3IonTypes('a');
  TEST_EQUAL(a_types.size(), 4)

  // y-precursor: a, b, y (standard)
  auto y_types = MS3FragmentMatcher::getMS3IonTypes('y');
  TEST_EQUAL(y_types.size(), 3)
  TEST_EQUAL(y_types[0], "a")
  TEST_EQUAL(y_types[1], "b")
  TEST_EQUAL(y_types[2], "y")

  // z-precursor: same as y-precursor
  auto z_types = MS3FragmentMatcher::getMS3IonTypes('z');
  TEST_EQUAL(z_types.size(), 3)

  // Unknown: defaults to y-precursor behavior
  auto unk_types = MS3FragmentMatcher::getMS3IonTypes('?');
  TEST_EQUAL(unk_types.size(), 3)
}
END_SECTION

START_SECTION(isPrefixIonType)
{
  TEST_TRUE(MS3FragmentMatcher::isPrefixIonType("a"))
  TEST_TRUE(MS3FragmentMatcher::isPrefixIonType("b"))
  TEST_TRUE(! MS3FragmentMatcher::isPrefixIonType("y"))
  TEST_TRUE(! MS3FragmentMatcher::isPrefixIonType("yb"))
  TEST_TRUE(! MS3FragmentMatcher::isPrefixIonType("ya"))
}
END_SECTION

START_SECTION(getIonShift)
{
  double co = 27.994915;
  double water = 18.010565;
  TEST_REAL_SIMILAR(MS3FragmentMatcher::getIonShift("a"), -co)
  TEST_REAL_SIMILAR(MS3FragmentMatcher::getIonShift("b"), 0.0)
  TEST_REAL_SIMILAR(MS3FragmentMatcher::getIonShift("y"), water)
  TEST_REAL_SIMILAR(MS3FragmentMatcher::getIonShift("yb"), 0.0)
  TEST_REAL_SIMILAR(MS3FragmentMatcher::getIonShift("ya"), -co)
}
END_SECTION

END_TEST
