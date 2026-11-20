// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <string>

namespace OpenMS
{
  /// Engine build 1: which build of the engine runs. The CI run id comes from a compile definition of BuildInfo.cpp alone
  /// (sources.cmake reads GITHUB_RUN_ID at configure), so a new run recompiles one small file; the commit is OpenMS's own
  /// configured revision (VersionInfo).
  struct OPENMS_DLLAPI FLASHIdaBuildInfo
  {
    /// The CI run that built this engine, "local" outside CI.
    static std::string ciRun();
    /// The OpenMS submodule commit the build was configured from ("" when git tracking was off).
    static std::string commit();
  };
} // namespace OpenMS
