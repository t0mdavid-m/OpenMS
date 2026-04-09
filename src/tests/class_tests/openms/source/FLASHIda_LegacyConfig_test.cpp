// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include <stdexcept>

using namespace OpenMS;

START_TEST(FLASHIda_LegacyConfig, "$Id$")

START_SECTION(legacy_config_rejected_via_Config)
{
  // Config directly rejects non-JSON input
  bool threw = false;
  try
  {
    Config cfg("10 100 1 10 5 0.5 -1");
  }
  catch (const std::invalid_argument&)
  {
    threw = true;
  }
  TEST_EQUAL(threw, true)
  (void)threw; // MSVC C4189
}
END_SECTION

START_SECTION(empty_config_rejected_via_Config)
{
  bool threw = false;
  try
  {
    Config cfg("");
  }
  catch (const std::invalid_argument&)
  {
    threw = true;
  }
  TEST_EQUAL(threw, true)
  (void)threw; // MSVC C4189
}
END_SECTION

START_SECTION(legacy_config_rejected_via_FLASHIda)
{
  // FLASHIda constructor delegates to Config and propagates the exception
  const char* legacy_input = "10 100 1 10 5 0.5 -1";
  bool threw = false;
  try
  {
    FLASHIda ida(const_cast<char*>(legacy_input));
  }
  catch (const std::invalid_argument&)
  {
    threw = true;
  }
  TEST_EQUAL(threw, true)
  (void)threw; // MSVC C4189
}
END_SECTION

START_SECTION(empty_config_rejected_via_FLASHIda)
{
  const char* empty_input = "";
  bool threw = false;
  try
  {
    FLASHIda ida(const_cast<char*>(empty_input));
  }
  catch (const std::invalid_argument&)
  {
    threw = true;
  }
  TEST_EQUAL(threw, true)
  (void)threw; // MSVC C4189
}
END_SECTION

END_TEST
