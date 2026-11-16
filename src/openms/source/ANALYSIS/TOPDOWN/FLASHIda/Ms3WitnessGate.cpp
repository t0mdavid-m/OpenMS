// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3WitnessGate.h>

#include <cctype>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace OpenMS
{
  namespace
  {
    const double kNaN = std::numeric_limits<double>::quiet_NaN();

    double readBack(const std::string& s)
    {
      try
      {
        return std::stod(s);
      }
      catch (const std::exception&)
      {
        return kNaN;
      }
    }

    // A raw value for the [MS3-W] line: as passed to the trees, enough digits to give it back exactly (the inputs are
    // at most 6 significant digits or 4 decimals of a mass); "nan" when missing.
    void put(std::ostringstream& os, double v)
    {
      if (std::isnan(v))
      {
        os << "nan";
        return;
      }
      std::ostringstream one;
      one << std::setprecision(12) << v;
      os << one.str();
    }

    // Which raw slots the line writes once per line, and which per candidate. A slot is a parent's only if every record
    // of a line holds the same value: queue_ahead differs per entry at a re-score, and cid3 and queue_ahead between a
    // plan's Q MS3 and its characterization MS3, which share the line.
    const int kParentSlots[] = {kWZMs2,         kWMassMs2,        kWPrecIntensity,      kWPrecQscore,      kWPrecIsoCos,
                                kWPrecSnr,      kWIdScore,        kWMs1Plan,            kWMs2MassCount,    kWMs2TicCoverage,
                                kWMs2PrecursorShare, kWMs2FragmentCount, kWMs2TagCount, kWMs2Ions,
                                kWPlanAge,      kWRtFrac};
    const int kCandidateSlots[] = {kWRank,      kWQueueAhead, kWCe1,       kWCid3,      kWIntensity,   kWIsoCos,
                                   kWQscore,    kWSnr,        kWChargeSnr, kWChargeCos, kWChargeScore, kWPpmError,
                                   kWMassMs3,   kWHasShift,   kWNMs2Prior, kWNMs3Prior};
  } // namespace

  double asLoggedSig6(double v)
  {
    if (std::isnan(v)) return v;
    std::ostringstream os;   // the default format every stage-scored scan_commands column and scan_results column uses
    os << v;
    return readBack(os.str());
  }

  double asLoggedFixed4(double v)
  {
    if (std::isnan(v)) return v;
    std::ostringstream os;
    os << std::fixed << std::setprecision(4) << v;
    return readBack(os.str());
  }

  std::string proformaLetters(const std::string& proforma)
  {
    size_t i = 0;
    if (!proforma.empty() && proforma[0] == '<')
    {
      const size_t close = proforma.find('>');
      i = (close == std::string::npos) ? proforma.size() : close + 1;
    }
    std::string out;
    int depth = 0;
    for (; i < proforma.size(); ++i)
    {
      const char c = proforma[i];
      if (c == '[') ++depth;
      else if (c == ']') depth = depth > 0 ? depth - 1 : 0;
      else if (depth == 0 && std::isalpha(static_cast<unsigned char>(c))) out.push_back(c);
    }
    return out;
  }

  bool proformaHasShift(const std::string& proforma)
  {
    size_t i = 0;
    if (!proforma.empty() && proforma[0] == '<')
    {
      const size_t close = proforma.find('>');
      i = (close == std::string::npos) ? proforma.size() : close + 1;
    }
    return proforma.find('[', i) != std::string::npos;
  }

  bool witnessPlaced(const std::string& parent_letters, const std::string& fragment_letters, bool n_terminal)
  {
    const size_t L = fragment_letters.size();
    if (L == 0 || L > parent_letters.size()) return false;
    const size_t off = n_terminal ? 0 : parent_letters.size() - L;
    return parent_letters.compare(off, L, fragment_letters) == 0;
  }

  double witnessMs1Plan(const WitnessLastMs1& ms1, double mono_mass)
  {
    static const double kOffsets[] = {0.0, 1.00335, -1.00335};
    double best = 0.0;
    for (const auto& [mass, intensity] : ms1.masses)
      for (double off : kOffsets)
      {
        const double target = mono_mass + off;
        if (target > 0 && std::fabs(mass - target) <= 10e-6 * target && intensity > best) best = intensity;
      }
    return best;
  }

  WitnessRecord WitnessGate::score(const std::string& ion_type, int ion_index, int charge, double ce1,
                                   const FragmentAnalysis::FragmentScores& s1, int rank) const
  {
    WitnessRecord r;
    r.ion_type = ion_type;
    r.ion_index = ion_index;
    r.charge = std::abs(charge);
    r.raw = ctx.parent;
    const char t = ion_type.empty() ? '\0' : static_cast<char>(std::tolower(static_cast<unsigned char>(ion_type[0])));
    r.raw[kWFragN] = (t == 'a' || t == 'b' || t == 'c') ? 1.0 : ((t == 'x' || t == 'y' || t == 'z') ? 0.0 : kNaN);
    r.raw[kWZMs3] = r.charge > 0 ? static_cast<double>(r.charge) : kNaN;
    r.raw[kWCe1] = asLoggedSig6(ce1);
    r.raw[kWCid3] = asLoggedSig6(cid3);
    r.raw[kWIntensity] = asLoggedSig6(s1.precursor_intensity);
    r.raw[kWL] = ion_index > 0 ? static_cast<double>(ion_index) : kNaN;
    r.raw[kWIsoCos] = asLoggedSig6(s1.iso_cos);
    r.raw[kWQscore] = asLoggedSig6(s1.qscore);
    r.raw[kWSnr] = asLoggedSig6(s1.snr);
    r.raw[kWMassMs3] = asLoggedFixed4(s1.mono_mass);
    r.raw[kWChargeSnr] = asLoggedSig6(s1.charge_snr);
    r.raw[kWChargeCos] = asLoggedSig6(s1.charge_cos);
    r.raw[kWChargeScore] = asLoggedSig6(s1.charge_score);
    r.raw[kWPpmError] = asLoggedSig6(s1.ppm_error);
    r.raw[kWQueueAhead] = queue_ahead;
    r.raw[kWRank] = static_cast<double>(rank);

    // has_shift and the placement read the fragment's ProForma string, as scan_commands' ms3_proteoform writes it.
    const std::string proforma = (fragment_proforma && ion_index > 0) ? fragment_proforma(ion_type, ion_index) : std::string();
    r.raw[kWHasShift] = proforma.empty() ? kNaN : (proformaHasShift(proforma) ? 1.0 : 0.0);
    const std::string letters = proformaLetters(proforma);
    r.placed = ctx.has_key && !letters.empty() && !std::isnan(r.raw[kWFragN])
               && witnessPlaced(ctx.parent_letters, letters, r.raw[kWFragN] == 1.0);
    if (r.placed) r.key = ctx.key;
    // The parent is one of the MS2 identifications under its own key (its row is written after its MS3 are planned).
    r.raw[kWNMs2Prior] = r.placed ? static_cast<double>(ctx.prior_ms2 + 1) : 1.0;
    r.raw[kWNMs3Prior] = r.placed ? static_cast<double>(ctx.prior_ms3) : 0.0;
    r.uncovered = uncovered_bonds ? uncovered_bonds(ion_type, ion_index) : -1;

    const int missing = r.raw.firstMissing();
    if (missing >= 0 || model == nullptr)
    {
      r.verdict = "untested";
      r.reason = std::string("missing_input=") + (missing >= 0 ? witnessRawName(missing) : "model");
      return r;
    }
    r.w = model->w(r.raw);
    if (witnessAdmits(r.w, threshold))
      r.verdict = "admitted";
    else
    {
      r.verdict = "refused";
      r.reason = "below_threshold";
    }
    return r;
  }

  std::string witnessLine(const WitnessGate& g, const char* stage, const std::string& extra,
                          const std::vector<WitnessRecord>& records)
  {
    if (records.empty()) return std::string();
    std::ostringstream os;
    os << "[MS3-W] stage=" << stage << extra << " parent=" << g.ctx.parent_id << " gate=" << g.gate_name
       << " model=w_ext sha=" << (g.model != nullptr ? g.model->sha8() : std::string("none")) << " threshold=";
    put(os, g.threshold);
    os << " walk=" << (g.rank_by_w ? "yield" : "intensity");
    for (int slot : kParentSlots)
    {
      os << ' ' << witnessRawName(slot) << '=';
      put(os, records.front().raw[slot]);   // the plan-level slots are the same in every record of the plan
    }
    os << " candidates=";
    for (size_t i = 0; i < records.size(); ++i)
    {
      const WitnessRecord& r = records[i];
      if (i > 0) os << ';';
      os << r.ion_type << r.ion_index << '/' << r.charge << ":w=";
      if (r.w < 0)
        os << "-1";
      else
        os << std::fixed << std::setprecision(8) << r.w << std::defaultfloat;
      os << ',' << r.verdict;
      if (!r.reason.empty()) os << ",reason=" << r.reason;
      os << ",unc=" << r.uncovered;
      for (int slot : kCandidateSlots)
      {
        os << ',' << witnessRawName(slot) << '=';
        put(os, r.raw[slot]);
      }
    }
    return os.str();
  }
} // namespace OpenMS
