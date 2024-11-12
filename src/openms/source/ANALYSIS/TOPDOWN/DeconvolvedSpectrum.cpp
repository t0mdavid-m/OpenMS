// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>

namespace OpenMS
{
  DeconvolvedSpectrum::DeconvolvedSpectrum(const int scan_number) : scan_number_(scan_number)
  {
  }

  bool DeconvolvedSpectrum::operator<(const DeconvolvedSpectrum& a) const
  {
    return this->scan_number_ < a.scan_number_;
  }

  bool DeconvolvedSpectrum::operator>(const DeconvolvedSpectrum& a) const
  {
    return this->scan_number_ > a.scan_number_;
  }

  bool DeconvolvedSpectrum::operator==(const DeconvolvedSpectrum& a) const
  {
    return this->scan_number_ == a.scan_number_;
  }

  /*
  std::vector<PeakGroup> DeconvolvedSpectrum::getNonOverlappingPeakGroups() const
  {
    std::map<double, std::set<int>> peak_to_pgs;
    std::vector<PeakGroup> filtered_pg_vec;
    filtered_pg_vec.reserve(size());
    std::vector<PeakGroup> ret;
    ret.reserve(size());

    int i = 0;
    for (const auto& pg : *this)
    {
      for (const auto& p : pg)
      {
        peak_to_pgs[p.mz].insert(i);
      }
      i++;
    }

    std::set<int> max_indices;
    for (const auto& e : peak_to_pgs)
    {
      int max_index = 0;
      double max_qscore = 0;
      for (const auto pg_index : e.second)
      {
        auto pg = peak_groups_[pg_index];
        if (max_qscore > pg.getFeatureQscore())
          continue;

        max_qscore = pg.getQscore2D();
        max_index = pg_index;
      }
      max_indices.insert(max_index);
    }
    for (int mi : max_indices) ret.push_back(peak_groups_[mi]);

    std::sort(ret.begin(), ret.end());
    return ret;
  }*/


  MSSpectrum DeconvolvedSpectrum::toSpectrum(const int to_charge, double tol, bool retain_undeconvolved)
  {
    auto out_spec = MSSpectrum(spec_);
    out_spec.clear(false);

   // if ((spec_.getMSLevel() > min_ms_level && precursor_peak_group_.empty()) || empty())
   // {
    //  return out_spec;
   // }
    double charge_mass_offset = (double)abs(to_charge) * FLASHHelperClasses::getChargeMass(to_charge >= 0);
    std::unordered_set<double> deconvolved_mzs;
    std::stringstream val {};

    val << "tol=" << tol << ";massoffset=" << std::to_string(charge_mass_offset) << ";chargemass=" << std::to_string(FLASHHelperClasses::getChargeMass(peak_groups_[0].isPositive()));
    if (!precursor_peak_group_.empty())
    {
      val << ";precursorscan=" << precursor_scan_number_ << ";precursormass=" << std::to_string(precursor_peak_group_.getMonoMass())
          << ";precursorscore=" << std::to_string(precursor_peak_group_.getQscore2D());
    }
    else
    {
      val << ";precursorscan=0;precursormass=0;precursorscore=0";
    }

    val << ";peaks=";
    out_spec.reserve(size());
    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }

      out_spec.emplace_back(pg.getMonoMass() + charge_mass_offset, pg.getIntensity());
      auto [z1, z2] = pg.getAbsChargeRange();
      int min_iso = -1, max_iso = 0;

      for (auto& p : pg)
      {
        min_iso = min_iso < 0 ? p.isotopeIndex : std::min(min_iso, p.isotopeIndex);
        max_iso = std::max(max_iso, p.isotopeIndex);
      }
      val << z1 << ":" << z2 << "," << min_iso << ":" << max_iso << ";";

      if (retain_undeconvolved)
      {
        for (auto& p : pg)
        {
          deconvolved_mzs.insert(p.mz);
        }
      }
    }

    val << "cos=";
    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      val << pg.getIsotopeCosine() << ",";
    }

    val << ";snr=";
    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      val << pg.getSNR() << ",";
    }

    val << ";qscore=";
    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      val << pg.getQscore2D() << ",";
    }

    val << ";qvalue=";
    for (auto& pg : *this)
    {
      if (pg.empty())
      {
        continue;
      }
      val << pg.getQvalue() << ",";
    }
    out_spec.setMetaValue("DeconvMassInfo", val.str());

    if (retain_undeconvolved)
    {
      for (const auto& p : spec_)
      {
        if (deconvolved_mzs.find(p.getMZ()) != deconvolved_mzs.end()) // if p is deconvolved
        {
          continue;
        }
        out_spec.emplace_back(p.getMZ() + charge_mass_offset - FLASHHelperClasses::getChargeMass(to_charge >= 0), p.getIntensity());
      }
    }
    out_spec.sortByPosition();
    if (!precursor_peak_group_.empty())
    {
      Precursor precursor(precursor_peak_);
      precursor.setCharge(to_charge);
      precursor.setMZ(precursor_peak_group_.getMonoMass() + charge_mass_offset);
      precursor.setIntensity(precursor_peak_group_.getIntensity());

      out_spec.getPrecursors().clear();
      out_spec.getPrecursors().emplace_back(precursor);
    }
    return out_spec;
  }

  const MSSpectrum& DeconvolvedSpectrum::getOriginalSpectrum() const
  {
    return spec_;
  }

  const PeakGroup& DeconvolvedSpectrum::getPrecursorPeakGroup() const
  {
    return precursor_peak_group_;
  }

  int DeconvolvedSpectrum::getPrecursorCharge() const
  {
    return precursor_peak_.getCharge();
  }

  double DeconvolvedSpectrum::getCurrentMaxMass(const double max_mass) const
  {
    if (spec_.getMSLevel() == 1 || precursor_peak_group_.empty())
    {
      return max_mass;
    }
    return std::min(max_mass, precursor_peak_group_.getMonoMass());
  }

  double DeconvolvedSpectrum::getCurrentMinMass(const double min_mass) const
  {
    if (spec_.getMSLevel() == 1)
    {
      return min_mass;
    }
    return 50.0;
  }

  int DeconvolvedSpectrum::getCurrentMaxAbsCharge(const int max_abs_charge) const
  {
    if (spec_.getMSLevel() == 1 || precursor_peak_group_.empty())
    {
      return max_abs_charge;
    }
    return std::min(max_abs_charge, abs(precursor_peak_group_.getRepAbsCharge()));
  }

  const Precursor& DeconvolvedSpectrum::getPrecursor() const
  {
    return precursor_peak_;
  }

  int DeconvolvedSpectrum::getScanNumber() const
  {
    return scan_number_;
  }

  int DeconvolvedSpectrum::getPrecursorScanNumber() const
  {
    return precursor_scan_number_;
  }

  const Precursor::ActivationMethod& DeconvolvedSpectrum::getActivationMethod() const
  {
    return activation_method_;
  }

  void DeconvolvedSpectrum::setPrecursor(const Precursor& precursor)
  {
    precursor_peak_ = precursor;
  }

  void DeconvolvedSpectrum::setPrecursorIntensity(const float i)
  {
    precursor_peak_.setIntensity(i);
  }

  void DeconvolvedSpectrum::setActivationMethod(const Precursor::ActivationMethod& method)
  {
    activation_method_ = method;
  }

  void DeconvolvedSpectrum::setPrecursorPeakGroup(const PeakGroup& pg)
  {
    precursor_peak_group_ = pg;
  }

  void DeconvolvedSpectrum::setOriginalSpectrum(const MSSpectrum& spec)
  {
    auto filter_str = spec.getMetaValue("filter string").toString();
    Size pos = filter_str.find("cv=");

    if (pos != String::npos)
    {
      Size end = filter_str.find(" ", pos);
      if (end == String::npos) end = filter_str.length() - 1;
      cv_ = std::stod(filter_str.substr(pos + 3, end - pos));
    }
    spec_ = spec;
  }

  double DeconvolvedSpectrum::getCV() const
  {
    return cv_;
  }

  void DeconvolvedSpectrum::setPrecursorScanNumber(const int scan_number)
  {
    precursor_scan_number_ = scan_number;
  }

  std::vector<PeakGroup>::const_iterator DeconvolvedSpectrum::begin() const noexcept
  {
    return peak_groups_.begin();
  }

  std::vector<PeakGroup>::const_iterator DeconvolvedSpectrum::end() const noexcept
  {
    return peak_groups_.end();
  }

  std::vector<PeakGroup>::iterator DeconvolvedSpectrum::begin() noexcept
  {
    return peak_groups_.begin();
  }

  std::vector<PeakGroup>::iterator DeconvolvedSpectrum::end() noexcept
  {
    return peak_groups_.end();
  }

  const PeakGroup& DeconvolvedSpectrum::operator[](const Size i) const
  {
    return peak_groups_[i];
  }

  PeakGroup& DeconvolvedSpectrum::operator[](const Size i)
  {
    return peak_groups_[i];
  }

  void DeconvolvedSpectrum::push_back(const PeakGroup& pg)
  {
    peak_groups_.push_back(pg);
  }

  void DeconvolvedSpectrum::emplace_back(const PeakGroup& pg)
  {
    peak_groups_.emplace_back(pg);
  }

  void DeconvolvedSpectrum::pop_back()
  {
    peak_groups_.pop_back();
  }

  PeakGroup& DeconvolvedSpectrum::back()
  {
    return peak_groups_.back();
  }

  Size DeconvolvedSpectrum::size() const noexcept
  {
    return peak_groups_.size();
  }

  void DeconvolvedSpectrum::clear()
  {
    std::vector<PeakGroup>().swap(peak_groups_);
  }

  void DeconvolvedSpectrum::reserve(Size n)
  {
    peak_groups_.reserve(n);
  }

  bool DeconvolvedSpectrum::empty() const
  {
    return peak_groups_.empty();
  }

  bool DeconvolvedSpectrum::isDecoy() const
  {
    if (empty())
      return false;
    if (peak_groups_[0].getTargetDecoyType() != PeakGroup::TargetDecoyType::target)
      return true;
    if (!precursor_peak_group_.empty() && precursor_peak_group_.getTargetDecoyType() != PeakGroup::TargetDecoyType::target)
      return true;
    return false;
  }

  FLASHHelperClasses::IsobaricQuantities DeconvolvedSpectrum::getQuantities() const
  {
    return quantities_;
  }

  void DeconvolvedSpectrum::setQuantities(const FLASHHelperClasses::IsobaricQuantities& quantities)
  {
    quantities_ = quantities;
  }

  void DeconvolvedSpectrum::setPeakGroups(std::vector<PeakGroup>& x)
  {
    std::vector<PeakGroup>().swap(peak_groups_);
    peak_groups_ = x;
  }

  void DeconvolvedSpectrum::sort()
  {
    std::sort(peak_groups_.begin(), peak_groups_.end());
  }

  void DeconvolvedSpectrum::sortByQscore()
  {
    std::sort(peak_groups_.begin(), peak_groups_.end(), [](const PeakGroup& p1, const PeakGroup& p2) { return p1.getQscore2D() > p2.getQscore2D(); });
  }
} // namespace OpenMS