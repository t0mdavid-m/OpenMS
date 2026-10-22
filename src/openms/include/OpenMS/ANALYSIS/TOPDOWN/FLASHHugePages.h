// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#pragma once

#include <cstdlib> // defines __GLIBC__ (via features.h) before the check below

#if defined(__GLIBC__) && defined(__linux__)
  #include <cstring>
  #include <string>
  #include <unistd.h>
#endif

namespace OpenMS
{
  /**
    @brief Called first in main() of FLASHDeconv and FLASHTnT: puts the heap on transparent huge pages (Linux with glibc only).

    glibc backs the heap with transparent huge pages only when asked at start-up (tunable glibc.malloc.hugetlb, read before main), so the
    process executes itself once more with it set. With 4 KiB pages, the per-thread arenas (~20 GB for FLASHDeconv at 192 threads) cost tens
    of seconds of page-fault system time and seconds of page-table teardown at exit; with huge pages FLASHDeconv is ~20 % and FLASHTnT ~23 %
    faster at 192 threads, with identical output. If the re-execution fails, the run continues as before. Skipped:
    - when GLIBC_TUNABLES already names the tunable (set by the user), FLASH_NO_HUGETLB is set, or after the re-execution (FLASH_HUGETLB_EXEC
      marks it, in case glibc drops the variable as in secure-execution mode);
    - with OpenMP thread binding requested (OMP_PROC_BIND, OMP_PLACES, GOMP_CPU_AFFINITY, KMP_AFFINITY): the OpenMP runtime then binds the
      main thread to the first place before main, and the re-executed process would inherit that one-CPU mask, so every thread would run on
      one CPU. With binding, set GLIBC_TUNABLES=glibc.malloc.hugetlb=1 directly;
    - under valgrind, where /proc/self/exe is valgrind's own executable.
  */
  inline void useTransparentHugePages(const char** argv)
  {
#if defined(__GLIBC__) && defined(__linux__)
    if (std::getenv("FLASH_NO_HUGETLB") != nullptr || std::getenv("FLASH_HUGETLB_EXEC") != nullptr) return;
    if (std::getenv("OMP_PROC_BIND") != nullptr || std::getenv("OMP_PLACES") != nullptr || std::getenv("GOMP_CPU_AFFINITY") != nullptr
        || std::getenv("KMP_AFFINITY") != nullptr)
      return;
    const char* preload = std::getenv("LD_PRELOAD");
    if (preload != nullptr && std::strstr(preload, "vgpreload") != nullptr) return;
    const char* tunables = std::getenv("GLIBC_TUNABLES");
    if (tunables != nullptr && std::strstr(tunables, "glibc.malloc.hugetlb") != nullptr) return;
    const std::string tunables_before = tunables != nullptr ? tunables : "";
    const std::string value = tunables != nullptr && *tunables != '\0' ? std::string(tunables) + ":glibc.malloc.hugetlb=1" : std::string("glibc.malloc.hugetlb=1");
    if (setenv("GLIBC_TUNABLES", value.c_str(), 1) != 0 || setenv("FLASH_HUGETLB_EXEC", "1", 1) != 0) return;
    execv("/proc/self/exe", const_cast<char* const*>(argv));
    unsetenv("FLASH_HUGETLB_EXEC"); // the re-execution failed: the environment as it was
    if (tunables != nullptr) { setenv("GLIBC_TUNABLES", tunables_before.c_str(), 1); }
    else { unsetenv("GLIBC_TUNABLES"); }
#else
    (void)argv;
#endif
  }
} // namespace OpenMS
