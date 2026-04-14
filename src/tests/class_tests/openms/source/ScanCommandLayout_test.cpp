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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
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
  std::printf("ScanCommand.enqueue_timestamp_ms.offset=%zu\n", offsetof(ScanCommand, enqueue_timestamp_ms));
  std::printf("ScanCommand.qscore.offset=%zu\n", offsetof(ScanCommand, qscore));
  std::printf("ScanCommand.mono_mass.offset=%zu\n", offsetof(ScanCommand, mono_mass));
  std::printf("ScanCommand.charge_cos.offset=%zu\n", offsetof(ScanCommand, charge_cos));
  std::printf("ScanCommand.charge_snr.offset=%zu\n", offsetof(ScanCommand, charge_snr));
  std::printf("ScanCommand.iso_cos.offset=%zu\n", offsetof(ScanCommand, iso_cos));
  std::printf("ScanCommand.snr.offset=%zu\n", offsetof(ScanCommand, snr));
  std::printf("ScanCommand.charge_score.offset=%zu\n", offsetof(ScanCommand, charge_score));
  std::printf("ScanCommand.ppm_error.offset=%zu\n", offsetof(ScanCommand, ppm_error));
  std::printf("ScanCommand.precursor_intensity.offset=%zu\n", offsetof(ScanCommand, precursor_intensity));
  std::printf("ScanCommand.peakgroup_intensity.offset=%zu\n", offsetof(ScanCommand, peakgroup_intensity));
  std::printf("ScanCommand.hcd_energy.offset=%zu\n", offsetof(ScanCommand, hcd_energy));
  std::printf("ScanCommand.pad2.offset=%zu\n", offsetof(ScanCommand, pad2));
  std::printf("ScanCommand.faims_cv.offset=%zu\n", offsetof(ScanCommand, faims_cv));
  std::printf("ScanCommand.microscans.offset=%zu\n", offsetof(ScanCommand, microscans));
  std::printf("ScanCommand.pad3.offset=%zu\n", offsetof(ScanCommand, pad3));
  std::printf("ScanCommand.rf_lens.offset=%zu\n", offsetof(ScanCommand, rf_lens));
  std::printf("ScanCommand.source_cid.offset=%zu\n", offsetof(ScanCommand, source_cid));
  std::printf("ScanCommand.source_cid_scaling.offset=%zu\n", offsetof(ScanCommand, source_cid_scaling));
  std::printf("ScanCommand.data_type.offset=%zu\n", offsetof(ScanCommand, data_type));
  std::printf("ScanCommand.scan_rate.offset=%zu\n", offsetof(ScanCommand, scan_rate));
  std::printf("ScanCommand.parent_scan_id.offset=%zu\n", offsetof(ScanCommand, parent_scan_id));
  std::printf("ScanCommand.reserved_.offset=%zu\n", offsetof(ScanCommand, reserved_));

  std::printf("LAYOUT_CHECK_PASSED\n");
  return 0;
}
