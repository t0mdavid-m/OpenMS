// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3YieldModel.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/EventLog.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <limits>
#include <sstream>

namespace OpenMS
{
  namespace
  {
    // Clip to [lo, hi]; a non-finite value (NaN, +-inf) takes @p fallback, one of the two bounds.
    double clipTo(double v, double lo, double hi, double fallback)
    {
      if (!std::isfinite(v)) return fallback;
      return std::min(std::max(v, lo), hi);
    }
  } // anonymous namespace

  YieldPrediction predictBondYield(const YieldCoefficients& c, const YieldInputs& in)
  {
    YieldPrediction out;
    if (in.length <= 0) return out;   // no ion index: the gate does not test it (ADR-0069 decision 7)
    using R = YieldClipRanges;
    const double i_frag = clipTo(in.i_frag, R::kFragIntensityMin, R::kFragIntensityMax, R::kFragIntensityMin);
    const double i_prec = clipTo(in.i_prec, R::kPrecIntensityMin, R::kPrecIntensityMax, R::kPrecIntensityMin);
    const double z_prec = clipTo(in.z_prec, R::kPrecChargeMin, R::kPrecChargeMax, R::kPrecChargeMin);
    // A precursor mass of 0 would divide by zero; it is never commanded, and the ratio then takes its upper bound.
    const double ratio = (in.m_prec > 0.0) ? in.m_frag / in.m_prec : std::numeric_limits<double>::quiet_NaN();
    const double mass_ratio = clipTo(ratio, R::kMassRatioMin, R::kMassRatioMax, R::kMassRatioMax);
    const double iso_cos = clipTo(in.iso_cos, R::kIsoCosMin, R::kIsoCosMax, R::kIsoCosMin);
    const double eta = c.intercept + c.fragment_intensity * std::log10(i_frag) + c.precursor_intensity * std::log10(i_prec)
                       + c.precursor_charge * z_prec + c.mass_ratio * mass_ratio + c.fragment_iso_cos * iso_cos;
    out.tested = true;
    out.p = 1.0 / (1.0 + std::exp(-eta));
    out.e = out.p * static_cast<double>(in.length - 1);
    return out;
  }

  bool yieldAdmits(const YieldGateConfig& gate, const YieldPrediction& pred)
  {
    return !pred.tested || pred.e >= gate.threshold;
  }

  YieldInputs yieldInputsOf(const ScanCommand& ms_ctx, const FragmentAnalysis::FragmentScores& s1, int ion_index)
  {
    // Each field is what buildMS3 (ScanCommandQueue.cpp) copies into the MS3 it builds from ms_ctx and frag_scores:
    // stages[0] = ms_ctx.stages[0], mono_mass / precursor_intensity from ms_ctx, the *_s1 block from frag_scores.
    YieldInputs in;
    in.i_frag = s1.precursor_intensity;
    in.i_prec = ms_ctx.precursor_intensity;
    in.z_prec = (ms_ctx.num_stages > 0) ? static_cast<double>(ms_ctx.stages[0].charge_state) : 0.0;
    in.m_frag = s1.mono_mass;
    in.m_prec = ms_ctx.mono_mass;
    in.iso_cos = s1.iso_cos;
    in.length = ion_index;
    return in;
  }

  YieldInputs yieldInputsOf(const ScanCommand& ms3)
  {
    YieldInputs in;
    in.i_frag = ms3.precursor_intensity_s1;
    in.i_prec = ms3.precursor_intensity;
    in.z_prec = (ms3.num_stages > 0) ? static_cast<double>(ms3.stages[0].charge_state) : 0.0;
    in.m_frag = ms3.mono_mass_s1;
    in.m_prec = ms3.mono_mass;
    in.iso_cos = ms3.iso_cos_s1;
    in.length = ms3IonOf(ms3).second;
    return in;
  }

  std::pair<std::string, int> ms3IonOf(const ScanCommand& cmd)
  {
    std::string ion_type;
    int ion_index = 0;
    if (cmd.msn_level == 3 && std::strlen(cmd.scan_description) > 4)
    {
      std::string desc(cmd.scan_description);
      // rfind: the charge delimiter is the LAST '@'. The 3-char base-94 tracking id can itself contain '@'
      // (0x40, e.g. id "!#@"); find() would lock onto the id's '@' and mis-parse the ion (logging the marker
      // 'R' as ion_type + a mass digit as ion_index). The mass token has no '@', so rfind('@') is exactly the
      // charge delimiter. (Matches the test-side decodeTrailingIon, which also uses rfind.)
      auto at_pos = desc.rfind('@');
      if (at_pos != std::string::npos)
      {
        size_t pos = at_pos + 1;
        while (pos < desc.size() && (std::isdigit(static_cast<unsigned char>(desc[pos])) || desc[pos] == '-'))
          pos++;
        if (pos < desc.size() && std::isalpha(static_cast<unsigned char>(desc[pos])))
        {
          ion_type = std::string(1, desc[pos]);
          pos++;
          if (pos < desc.size())
            ion_index = std::atoi(desc.c_str() + pos);
        }
      }
    }
    return {ion_type, ion_index};
  }

  std::string yieldEventFields(const YieldGateConfig& gate, const YieldInputs& in, const YieldPrediction& pred)
  {
    std::ostringstream os;
    os << std::setprecision(6);
    os << " e=" << pred.e << " p=" << pred.p << " threshold=" << gate.threshold << " i_frag=" << in.i_frag
       << " i_prec=" << in.i_prec << " z_prec=" << in.z_prec << " m_frag=" << in.m_frag << " m_prec=" << in.m_prec
       << " iso_cos=" << in.iso_cos << " length=" << in.length;
    return os.str();
  }

  void logYieldVerdict(const char* gate_name, int precursor_id, const std::string& ion_type, int ion_index, int charge,
                       const char* verdict, const YieldGateConfig& gate, const YieldInputs& in, const YieldPrediction& pred)
  {
    EventLine ev;
    ev << "[MS3-YIELD] gate=" << gate_name << " precursor_id=" << precursor_id << " ion=" << ion_type << ion_index
       << " charge=" << charge << " verdict=" << verdict;
    if (!pred.tested) ev << " reason=no_ion_index";
    ev << yieldEventFields(gate, in, pred);
  }
} // namespace OpenMS
