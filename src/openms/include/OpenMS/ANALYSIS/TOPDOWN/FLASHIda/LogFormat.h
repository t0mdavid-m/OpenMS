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
  struct ScanCommand;

  /**
    @brief The number formats the five log streams print, in one place (ADR-0079).

    A value the W gate reads "as logged" must be the bytes the logger writes, to the last digit: the model was fitted on
    the log strings (the ms3_esm package, SPEC 6-7). So the logger and the gate call these functions, never two copies of
    the formatting. Every function uses a fresh local stream in the "C" locale (the engine never imbues one), so no
    precision leaks between cells.
  */

  /// @p v as a default-formatted ostream prints it: 6 significant digits (every stage-scored scan_commands column, every
  /// scan_results float).
  OPENMS_DLLAPI std::string formatSig6(double v);
  /// @p v as std::fixed / setprecision(4) prints it (masses, m/z bounds, identification_score).
  OPENMS_DLLAPI std::string formatFixed4(double v);
  /// A logged cell read back as a double; NaN when it does not parse.
  OPENMS_DLLAPI double parseLogged(const std::string& cell);

  /**
    @brief The per-stage and scored cells of one scan_commands.tsv row, as writeScanCommandRow prints them.

    Stages are joined by ';' (stage 0 = the MS2 precursor, the last stage = an MS3's target fragment); co-isolation
    notches by ',' inside their stage's group (charge, m/z and isolation width only). The scored cells carry stage 0 and,
    on an MS3, ';' and stage 1. A stage-less command (MS1, AGC) prints "0" in the per-stage cells and "none" as its
    activation, and "0" as its mono_mass.
  */
  struct OPENMS_DLLAPI ScanCommandCells
  {
    std::string mono_mass, charge, precursor_mz, isolation_width, collision_energy, activation;
    std::string reaction_time, reagent_max_it, reagent_agc_target;
    std::string qscore, charge_cos, charge_snr, iso_cos, snr, charge_score, hcd_energy;
    std::string ppm_error, precursor_intensity, peakgroup_intensity;
    std::string first_mass, last_mass, max_it;
  };

  /// The cells of @p cmd exactly as scan_commands.tsv prints them.
  OPENMS_DLLAPI ScanCommandCells scanCommandCells(const ScanCommand& cmd);

} // namespace OpenMS
