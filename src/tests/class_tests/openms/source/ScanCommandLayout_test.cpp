// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

// Layout query binary for ScanCommand and IsolationStage structs.
// Outputs sizeof and offsetof values in a parseable format.
// Used by CI to verify C++/C# struct layout agreement.

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <cstddef>
#include <cstdio>

using namespace OpenMS;

int main()
{
  // IsolationStage layout
  std::printf("IsolationStage.sizeof=%zu\n", sizeof(IsolationStage));
  std::printf("IsolationStage.precursor_mz.offset=%zu\n", offsetof(IsolationStage, precursor_mz));
  std::printf("IsolationStage.isolation_width.offset=%zu\n", offsetof(IsolationStage, isolation_width));
  std::printf("IsolationStage.collision_energy.offset=%zu\n", offsetof(IsolationStage, collision_energy));
  std::printf("IsolationStage.reaction_time.offset=%zu\n", offsetof(IsolationStage, reaction_time));
  std::printf("IsolationStage.reagent_max_it.offset=%zu\n", offsetof(IsolationStage, reagent_max_it));
  std::printf("IsolationStage.reagent_agc_target.offset=%zu\n", offsetof(IsolationStage, reagent_agc_target));
  std::printf("IsolationStage.charge_state.offset=%zu\n", offsetof(IsolationStage, charge_state));
  std::printf("IsolationStage.activation_type.offset=%zu\n", offsetof(IsolationStage, activation_type));

  // ScanCommand layout
  std::printf("ScanCommand.sizeof=%zu\n", sizeof(ScanCommand));
  std::printf("ScanCommand.scan_id.offset=%zu\n", offsetof(ScanCommand, scan_id));
  std::printf("ScanCommand.msn_level.offset=%zu\n", offsetof(ScanCommand, msn_level));
  std::printf("ScanCommand.priority.offset=%zu\n", offsetof(ScanCommand, priority));
  std::printf("ScanCommand.is_agc.offset=%zu\n", offsetof(ScanCommand, is_agc));
  std::printf("ScanCommand.num_stages.offset=%zu\n", offsetof(ScanCommand, num_stages));
  std::printf("ScanCommand.orbitrap_resolution.offset=%zu\n", offsetof(ScanCommand, orbitrap_resolution));
  std::printf("ScanCommand.agc_target.offset=%zu\n", offsetof(ScanCommand, agc_target));
  std::printf("ScanCommand.pad1.offset=%zu\n", offsetof(ScanCommand, pad1));
  std::printf("ScanCommand.first_mass.offset=%zu\n", offsetof(ScanCommand, first_mass));
  std::printf("ScanCommand.last_mass.offset=%zu\n", offsetof(ScanCommand, last_mass));
  std::printf("ScanCommand.max_it.offset=%zu\n", offsetof(ScanCommand, max_it));
  std::printf("ScanCommand.analyzer.offset=%zu\n", offsetof(ScanCommand, analyzer));
  std::printf("ScanCommand.scan_description.offset=%zu\n", offsetof(ScanCommand, scan_description));
  std::printf("ScanCommand.stages.offset=%zu\n", offsetof(ScanCommand, stages));

  std::printf("LAYOUT_CHECK_PASSED\n");
  return 0;
}
