// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/LogFormat.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include <iomanip>
#include <limits>
#include <sstream>

namespace OpenMS
{
  std::string formatSig6(double v)
  {
    std::ostringstream os;
    os << v;
    return os.str();
  }

  std::string formatFixed4(double v)
  {
    std::ostringstream os;
    os << std::fixed << std::setprecision(4) << v;
    return os.str();
  }

  double parseLogged(const std::string& cell)
  {
    try
    {
      return std::stod(cell);
    }
    catch (const std::exception&)
    {
      return std::numeric_limits<double>::quiet_NaN();
    }
  }

  ScanCommandCells scanCommandCells(const ScanCommand& cmd)
  {
    ScanCommandCells c;
    for (int i = 0; i < cmd.num_stages; ++i)
    {
      if (i > 0)
      {
        c.charge += ";"; c.precursor_mz += ";"; c.isolation_width += ";"; c.collision_energy += ";"; c.activation += ";";
        c.reaction_time += ";"; c.reagent_max_it += ";"; c.reagent_agc_target += ";";
      }
      c.charge += std::to_string(cmd.stages[i].charge_state);
      c.precursor_mz += formatSig6(cmd.stages[i].precursor_mz);
      c.isolation_width += formatSig6(cmd.stages[i].isolation_width);
      c.collision_energy += formatSig6(cmd.stages[i].collision_energy);
      c.activation += cmd.stages[i].activation_type;
      // Co-isolation notches, ','-joined INSIDE this stage's ';'-group -- the same two-axis grammar the wire uses, so
      // scan_commands.tsv records what was actually commanded rather than only the anchor window. Only the
      // charge-dependent geometry gains the axis: collision energy and activation are per fragmentation event and all
      // notches of a stage share one. No notches (every mode but multiplexed) emits nothing here.
      auto stage_notches = notchesForStage(cmd, i);
      for (int n = 0; n < stage_notches.second; ++n)
      {
        c.charge += "," + std::to_string(stage_notches.first[n].charge_state);
        c.precursor_mz += "," + formatSig6(stage_notches.first[n].precursor_mz);
        c.isolation_width += "," + formatSig6(stage_notches.first[n].isolation_width);
      }
      c.reaction_time += formatSig6(cmd.stages[i].reaction_time);
      c.reagent_max_it += formatSig6(cmd.stages[i].reagent_max_it);
      c.reagent_agc_target += std::to_string(cmd.stages[i].reagent_agc_target);
    }
    // Stage-less commands (MS1 / AGC) populate no per-stage fields. Explicit placeholders keep the row's column count:
    // a tab-split parser drops a trailing empty field.
    if (cmd.num_stages == 0)
    {
      c.charge = c.precursor_mz = c.isolation_width = c.collision_energy = c.reaction_time = c.reagent_max_it =
          c.reagent_agc_target = "0";
      c.activation = "none";
    }
    // The scored columns are 2-stage 'MS2-precursor;fragment' on an MS3; every other row prints stage 0 only.
    const bool ms3 = cmd.msn_level == 3;
    auto sc = [&](double v0, double v1) { return ms3 ? formatSig6(v0) + ";" + formatSig6(v1) : formatSig6(v0); };
    auto sci = [&](int v0, int v1) { return ms3 ? std::to_string(v0) + ";" + std::to_string(v1) : std::to_string(v0); };
    // mono_mass sets its own format (ADR-0046 decision 10): 4 decimals, both stages. Stage-less rows keep their "0".
    c.mono_mass = cmd.num_stages == 0 ? std::string("0")
                                      : (ms3 ? formatFixed4(cmd.mono_mass) + ";" + formatFixed4(cmd.mono_mass_s1)
                                             : formatFixed4(cmd.mono_mass));
    c.qscore = sc(cmd.qscore, cmd.qscore_s1);
    c.charge_cos = sc(cmd.charge_cos, cmd.charge_cos_s1);
    c.charge_snr = sc(cmd.charge_snr, cmd.charge_snr_s1);
    c.iso_cos = sc(cmd.iso_cos, cmd.iso_cos_s1);
    c.snr = sc(cmd.snr, cmd.snr_s1);
    c.charge_score = sc(cmd.charge_score, cmd.charge_score_s1);
    c.hcd_energy = sci(cmd.hcd_energy, cmd.hcd_energy_s1);
    c.ppm_error = sc(cmd.ppm_error, cmd.ppm_error_s1);
    c.precursor_intensity = sc(cmd.precursor_intensity, cmd.precursor_intensity_s1);
    c.peakgroup_intensity = sc(cmd.peakgroup_intensity, cmd.peakgroup_intensity_s1);
    // Scan-range bounds and max_it: fixed + 4 (ADR-0026 d6, ADR-0072).
    c.first_mass = formatFixed4(cmd.first_mass);
    c.last_mass = formatFixed4(cmd.last_mass);
    c.max_it = formatFixed4(cmd.max_it);
    return c;
  }

} // namespace OpenMS
