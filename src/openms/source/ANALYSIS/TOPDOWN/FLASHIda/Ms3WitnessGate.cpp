// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3WitnessGate.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/LogFormat.h>

#include <algorithm>
#include <cctype>
#include <cmath>

namespace OpenMS
{
  namespace
  {
    const double kNaN = std::numeric_limits<double>::quiet_NaN();

    /// Where a ProForma string's residues start: after a leading <...> (only when it is closed, as the reference's
    /// regular expression ^<[^>]*> removes it).
    size_t afterGlobalMods(const std::string& proforma)
    {
      if (!proforma.empty() && proforma[0] == '<')
      {
        const size_t close = proforma.find('>');
        if (close != std::string::npos) return close + 1;
      }
      return 0;
    }
  } // namespace

  double asLoggedSig6(double v)
  {
    if (std::isnan(v)) return v;
    return parseLogged(formatSig6(v));   // ADR-0079: the logger's own formatter (LogFormat.h), never a second copy
  }

  double asLoggedFixed4(double v)
  {
    if (std::isnan(v)) return v;
    return parseLogged(formatFixed4(v));
  }

  std::string proformaLetters(const std::string& proforma)
  {
    std::string out;
    int depth = 0;
    for (size_t i = afterGlobalMods(proforma); i < proforma.size(); ++i)
    {
      const char c = proforma[i];
      if (c == '[') ++depth;
      else if (c == ']') depth = depth > 0 ? depth - 1 : 0;
      else if (depth == 0 && std::isalpha(static_cast<unsigned char>(c))) out.push_back(c);
    }
    return out;
  }

  int proformaResidueCount(const std::string& proforma)
  {
    // inputs_from_log_v17.residues: re.sub(r"[^A-Z]", "", re.sub(r"\[[^\]]*\]", "", re.sub(r"^<[^>]*>", "", s)))
    int n = 0;
    size_t i = afterGlobalMods(proforma);
    while (i < proforma.size())
    {
      const char c = proforma[i];
      if (c == '[')
      {
        const size_t close = proforma.find(']', i + 1);
        if (close != std::string::npos)
        {
          i = close + 1;
          continue;
        }
      }
      if (c >= 'A' && c <= 'Z') ++n;
      ++i;
    }
    return n;
  }

  bool proformaHasShift(const std::string& proforma) { return proforma.find('[', afterGlobalMods(proforma)) != std::string::npos; }

  double loggedStage(const std::string& cell, int stage)
  {
    std::vector<std::string> groups;
    size_t s = 0;
    while (true)
    {
      const size_t semi = cell.find(';', s);
      groups.push_back(cell.substr(s, semi == std::string::npos ? std::string::npos : semi - s));
      if (semi == std::string::npos) break;
      s = semi + 1;
    }
    const std::string& g = stage >= 0 ? (static_cast<size_t>(stage) < groups.size() ? groups[static_cast<size_t>(stage)] : groups.back())
                                      : groups.back();
    const std::string anchor = g.substr(0, g.find(','));
    if (anchor.empty()) return kNaN;
    return parseLogged(anchor);
  }

  double witnessMs1Plan(const WitnessLastMs1& ms1, double mono_mass)
  {
    static const double kOffsets[] = {0.0, 1.00335, -1.00335};
    double best = 0.0;
    for (const auto& [mass, intensity] : ms1.masses)
      for (double off : kOffsets)
      {
        const double target = mono_mass + off;
        // inputs_from_log_v17: abs(mi - (mass + off)) <= 1e-5 * mass -- the tolerance on the precursor mass itself.
        if (std::fabs(mass - target) <= 1e-5 * mono_mass && intensity > best) best = intensity;
      }
    return best;
  }

  size_t WitnessGate::scoreTarget(const WitnessCandidate& c, int planner_rank)
  {
    ++candidates;
    const size_t first = records.size();
    const std::vector<std::pair<ScanCommand, std::string>> cmds = compose ? compose(c) : std::vector<std::pair<ScanCommand, std::string>>{};
    const WitnessBonds b = bonds ? bonds(c.ion_type, c.ion_index) : WitnessBonds{};
    const char t = c.ion_type.empty() ? '\0' : static_cast<char>(std::tolower(static_cast<unsigned char>(c.ion_type[0])));
    const double frag_n = (t == 'a' || t == 'b' || t == 'c') ? 1.0 : ((t == 'x' || t == 'y' || t == 'z') ? 0.0 : kNaN);

    auto blank = [&](int variant) {
      WitnessRecord r;
      r.candidate = candidates;
      r.variant = variant;
      r.planner_rank = planner_rank;
      r.ion_type = c.ion_type;
      r.ion_index = c.ion_index;
      r.bonds = b;
      r.raw.fill(kNaN);
      r.x.fill(std::numeric_limits<float>::quiet_NaN());
      return r;
    };
    if (cmds.empty())
    {
      WitnessRecord r = blank(0);
      r.failed_input = "command";
      records.push_back(std::move(r));
    }
    for (size_t v = 0; v < cmds.size(); ++v)
    {
      WitnessRecord r = blank(static_cast<int>(v));
      r.cmd = cmds[v].first;
      r.ms3_proteoform = cmds[v].second;
      const ScanCommandCells cells = scanCommandCells(r.cmd);
      auto& raw = r.raw;
      // The candidate's own cells, as its scan_commands.tsv row would print them (inputs_from_log_v17.row_inputs).
      raw[kW17FragN] = frag_n;
      raw[kW17ZMs2] = loggedStage(cells.charge, 0);
      raw[kW17ZMs3] = loggedStage(cells.charge, -1);
      raw[kW17Ce1] = loggedStage(cells.collision_energy, 0);
      raw[kW17Cid3] = loggedStage(cells.collision_energy, -1);
      raw[kW17Intensity] = loggedStage(cells.precursor_intensity, -1);
      raw[kW17L] = r.ms3_proteoform.empty() ? kNaN : static_cast<double>(proformaResidueCount(r.ms3_proteoform));
      raw[kW17IsoCos] = loggedStage(cells.iso_cos, -1);
      raw[kW17Qscore] = loggedStage(cells.qscore, -1);
      raw[kW17Snr] = loggedStage(cells.snr, -1);
      raw[kW17HasShift] = r.ms3_proteoform.empty() ? kNaN : (proformaHasShift(r.ms3_proteoform) ? 1.0 : 0.0);
      raw[kW17MassMs2] = loggedStage(cells.mono_mass, 0);
      raw[kW17MassMs3] = loggedStage(cells.mono_mass, -1);
      raw[kW17PrecIntensity] = loggedStage(cells.precursor_intensity, 0);
      raw[kW17PrecQscore] = loggedStage(cells.qscore, 0);
      raw[kW17PrecIsoCos] = loggedStage(cells.iso_cos, 0);
      raw[kW17PrecSnr] = loggedStage(cells.snr, 0);
      raw[kW17ChargeSnr] = loggedStage(cells.charge_snr, -1);
      raw[kW17ChargeCos] = loggedStage(cells.charge_cos, -1);
      raw[kW17ChargeScore] = loggedStage(cells.charge_score, -1);
      raw[kW17PpmError] = loggedStage(cells.ppm_error, -1);
      // The parent's rows and the last MS1.
      raw[kW17IdScore] = ctx.parent.id_score;
      raw[kW17Ms1Plan] = ctx.x_ms1_plan;
      raw[kW17Ms2MassCount] = ctx.parent.mass_count;
      raw[kW17Ms2TicCoverage] = ctx.parent.tic_coverage;
      raw[kW17Ms2PrecursorShare] = ctx.parent.precursor_share;
      raw[kW17Ms2FragmentCount] = ctx.parent.fragment_count;
      raw[kW17Ms2TagCount] = ctx.parent.tag_count;
      raw[kW17Ms2Ions] = ctx.parent.ms2_ions;

      for (int i = 0; i < kW17RawCount && r.failed_input.empty(); ++i)
        if (std::isnan(raw[static_cast<size_t>(i)])) r.failed_input = w17RawName(i);
      r.x = w17ModelInputs(raw, ctx.pk3);
      if (r.failed_input.empty())
      {
        if (model == nullptr)
          r.failed_input = "model";
        else
        {
          r.w = model->w(r.x);
          if (!std::isfinite(r.w)) r.failed_input = "w";
        }
      }
      records.push_back(std::move(r));
    }
    // One verdict per target: its largest W; none when any of its commands failed (SPEC 8).
    bool failed = false;
    double tw = -1.0;
    for (size_t i = first; i < records.size(); ++i)
    {
      if (!records[i].failed_input.empty()) failed = true;
      else tw = std::max(tw, records[i].w);
    }
    for (size_t i = first; i < records.size(); ++i) records[i].target_w = failed ? kNaN : tw;
    return first;
  }
} // namespace OpenMS
