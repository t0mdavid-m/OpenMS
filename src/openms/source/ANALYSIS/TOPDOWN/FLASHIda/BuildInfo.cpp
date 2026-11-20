// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/BuildInfo.h>

#include <OpenMS/CONCEPT/VersionInfo.h>

#ifndef FLASHIDA_CI_RUN
#define FLASHIDA_CI_RUN 0
#endif

#define FLASHIDA_STR2(x) #x
#define FLASHIDA_STR(x) FLASHIDA_STR2(x)

namespace OpenMS
{
  std::string FLASHIdaBuildInfo::ciRun()
  {
    const std::string run = FLASHIDA_STR(FLASHIDA_CI_RUN);
    return run == "0" ? std::string("local") : run;
  }

  std::string FLASHIdaBuildInfo::commit()
  {
    return std::string(VersionInfo::getRevision());
  }
} // namespace OpenMS
