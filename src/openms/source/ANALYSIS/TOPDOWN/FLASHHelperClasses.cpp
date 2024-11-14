// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <utility>

namespace OpenMS
{
FLASHHelperClasses::PrecalculatedAveragine::PrecalculatedAveragine(const double min_mass, const double max_mass, const double delta, CoarseIsotopePatternGenerator& generator,
                                                                           const bool use_RNA_averagine, double decoy_iso_distance, bool is_centroid) :
      mass_interval_(delta),
      min_mass_(min_mass)
  {
    int i = 0;
    while (true)
    {
      double mass = i * mass_interval_;
      i++;
      if (mass < min_mass)
      {
        continue;
      }
      if (mass > max_mass)
      {
        break;
      }

      auto iso = use_RNA_averagine ? generator.estimateFromRNAMonoWeight(mass) : generator.estimateFromPeptideMonoWeight(mass);

      if (decoy_iso_distance > 0)
      {
        auto decoy_iso(iso);
        double max_intensity = iso.getMostAbundant().getIntensity();

        for (Size k = 0; k < iso.size(); k++)
        {
          Size index = (k % 2 == 0 ? k + 1 : k - 1) % iso.size();
          float intensity = iso[index].getIntensity() / max_intensity;
          if (is_centroid)
            intensity = std::pow(intensity, 3);
          else
            intensity = std::pow(intensity, .3333);
          decoy_iso[k].setIntensity(intensity);
          //std::cout<<iso[index].getIntensity() << " " << decoy_iso[k].getIntensity() << "\n";
          decoy_iso[k].setMZ(iso[k].getMZ() * decoy_iso_distance);
        }
        decoy_iso.sortByMass();
        decoy_iso.renormalize();
        iso = decoy_iso;
      }

      const double min_pwr = .9999;
      const Size min_iso_length = 2;
      const int min_left_right_count = 2;
      double total_pwr = .0;
      size_t most_abundant_index_ = 0;
      double most_abundant_int = 0;

      /// sum of squared intensities to see the total power of isotope pattern. The range of isotope pattern is
      /// determined so those within range cover min_pwr of the total power.
      for (Size k = 0; k < iso.size(); k++)
      {
        total_pwr += iso[k].getIntensity() * iso[k].getIntensity();
        if (most_abundant_int >= iso[k].getIntensity())
        {
          continue;
        }
        most_abundant_int = iso[k].getIntensity();
        most_abundant_index_ = k;
      }

      int left_count = 0;
      int right_count = (int)iso.size() - 1;
      int trim_count = 0;

      while (iso.size() - trim_count > min_iso_length && left_count < right_count)
      {
        double lint = iso[left_count].getIntensity();
        double rint = iso[right_count].getIntensity();
        double pwr;
        bool trim_left = true;
        if (lint < rint)
        {
          pwr = lint * lint;
        }
        else
        {
          pwr = rint * rint;
          trim_left = false;
        }
        if (total_pwr - pwr < total_pwr * min_pwr)
        {
          break;
        }
        total_pwr -= pwr;
        trim_count++;
        if (trim_left)
        {
          iso[left_count].setIntensity(0);
          left_count++;
        }
        else
        {
          iso[right_count].setIntensity(0); // for trimming
          right_count--;
        }
      }
      left_count = (int)most_abundant_index_ - left_count;
      right_count = right_count - (int)most_abundant_index_;

      double intensity_sum = 0;
      for (auto& k : iso)
      {
        float ori_int = k.getIntensity();
        k.setIntensity(ori_int / (float)sqrt(total_pwr));
        intensity_sum += k.getIntensity();
      }
      left_count = left_count < min_left_right_count ? min_left_right_count : left_count;
      right_count = right_count < min_left_right_count ? min_left_right_count : right_count;

      apex_index_.push_back(most_abundant_index_);
      right_count_from_apex_.push_back(right_count);
      left_count_from_apex_.push_back(left_count);
      average_mono_mass_difference_.push_back(iso.averageMass() - iso[0].getMZ());
      abundant_mono_mass_difference_.push_back(iso.getMostAbundant().getMZ() - iso[0].getMZ());
      isotopes_.push_back(iso);
      snr_mul_factor_.push_back(intensity_sum * intensity_sum);
    }
  }

  Size FLASHHelperClasses::PrecalculatedAveragine::massToIndex_(const double mass) const
  {
    Size i = (Size)round(std::max(.0, mass - min_mass_) / mass_interval_);
    i = std::min(i, isotopes_.size() - 1);
    return i;
  }

  IsotopeDistribution FLASHHelperClasses::PrecalculatedAveragine::get(const double mass) const
  {
    return isotopes_[massToIndex_(mass)];
  }

  size_t FLASHHelperClasses::PrecalculatedAveragine::getMaxIsotopeIndex() const
  {
    return max_isotope_index_;
  }

  Size FLASHHelperClasses::PrecalculatedAveragine::getLeftCountFromApex(const double mass) const
  {
    return (Size)left_count_from_apex_[massToIndex_(mass)];
  }

  double FLASHHelperClasses::PrecalculatedAveragine::getAverageMassDelta(const double mass) const
  {
    return average_mono_mass_difference_[massToIndex_(mass)];
  }

  double FLASHHelperClasses::PrecalculatedAveragine::getMostAbundantMassDelta(const double mass) const
  {
    return abundant_mono_mass_difference_[massToIndex_(mass)];
  }

  double FLASHHelperClasses::PrecalculatedAveragine::getSNRMultiplicationFactor(const double mass) const
  {
    return snr_mul_factor_[massToIndex_(mass)];
  }

  Size FLASHHelperClasses::PrecalculatedAveragine::getRightCountFromApex(const double mass) const
  {
    return (Size)right_count_from_apex_[massToIndex_(mass)];
  }

  Size FLASHHelperClasses::PrecalculatedAveragine::getApexIndex(const double mass) const
  {
    return apex_index_[massToIndex_(mass)];
  }

  Size FLASHHelperClasses::PrecalculatedAveragine::getLastIndex(const double mass) const
  {
    Size index = massToIndex_(mass);
    return apex_index_[index] + right_count_from_apex_[index];
  }

  void FLASHHelperClasses::PrecalculatedAveragine::setMaxIsotopeIndex(const int index)
  {
    max_isotope_index_ = index;
  }

  FLASHHelperClasses::LogMzPeak::LogMzPeak(const Peak1D& peak, const bool positive) :
      mz(peak.getMZ()), intensity(peak.getIntensity()), logMz(getLogMz(peak.getMZ(), positive)), abs_charge(0), is_positive(positive), isotopeIndex(0)
  {
  }

  double FLASHHelperClasses::LogMzPeak::getUnchargedMass() const
  {
    if (abs_charge == 0)
    {
      return .0;
    }
    if (mass <= 0)
    {
      return (mz - getChargeMass(is_positive)) * (float)abs_charge;
    }
    return mass;
  }

  bool FLASHHelperClasses::LogMzPeak::operator<(const LogMzPeak& a) const
  {
    if (this->logMz == a.logMz)
    {
      return this->intensity < a.intensity;
    }
    return this->logMz < a.logMz;
  }

  bool FLASHHelperClasses::LogMzPeak::operator>(const LogMzPeak& a) const
  {
    if (this->logMz == a.logMz)
    {
      return this->intensity > a.intensity;
    }
    return this->logMz > a.logMz;
  }

  bool FLASHHelperClasses::LogMzPeak::operator==(const LogMzPeak& a) const
  {
    return this->logMz == a.logMz && this->intensity == a.intensity;
  }


  float FLASHHelperClasses::getChargeMass(const bool positive_ioniziation_mode)
  {
    return (float)(positive_ioniziation_mode ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U);
  }

  double FLASHHelperClasses::getLogMz(const double mz, const bool positive)
  {
    return std::log(mz - getChargeMass(positive));
  }

  bool FLASHHelperClasses::IsobaricQuantities::empty() const
  {
    return quantities.empty();
  }
  FLASHHelperClasses::Tag::Tag(String seq, double n_mass, double c_mass, std::vector<double>& mzs, std::vector<int>& scores,  int scan) :
      seq_(std::move(seq)), n_mass_(n_mass), c_mass_(c_mass), mzs_(mzs), scores_(scores), scan_(scan), length_(mzs.size() - 1)
  {
    upper_seq_ = seq_.toUpper();
  }

  const String& FLASHHelperClasses::Tag::getSequence() const
  {
    return seq_;
  }

  const String& FLASHHelperClasses::Tag::getUppercaseSequence() const
  {
    return upper_seq_;
  }

  Size FLASHHelperClasses::Tag::getLength() const
  {
    return length_;
  }


  const std::vector<double>& FLASHHelperClasses::Tag::getMzs() const
  {
    return mzs_;
  }

  double FLASHHelperClasses::Tag::getNtermMass() const
  {
    return n_mass_;
  }

  double FLASHHelperClasses::Tag::getCtermMass() const
  {
    return c_mass_;
  }

  int FLASHHelperClasses::Tag::getScore() const
  {
    return std::accumulate(scores_.begin(), scores_.end(), 0);
  }

  int FLASHHelperClasses::Tag::getScore(int pos) const
  {
    if (pos < 0 || pos >= scores_.size()) return 0;
    return scores_[pos];
  }

  int FLASHHelperClasses::Tag::getScan() const
  {
    return scan_;
  }

  bool FLASHHelperClasses::Tag::operator<(const Tag& a) const
  {
    if (this->length_ == a.length_)
    {
      if (this->seq_ == a.seq_)
        return (this->c_mass_ + this->n_mass_) < (a.c_mass_ + a.n_mass_);
      return this->seq_ > a.seq_;
    }
    return this->length_ < a.length_;
  }

  bool FLASHHelperClasses::Tag::operator>(const Tag& a) const
  {
    if (this->length_ == a.length_)
    {
      if (this->seq_ == a.seq_)
        return (this->c_mass_ + this->n_mass_) > (a.c_mass_ + a.n_mass_);
      return this->seq_ < a.seq_;
    }
    return this->length_ > a.length_;
  }

  bool FLASHHelperClasses::Tag::operator==(const Tag& a) const
  {
    return this->seq_ == a.seq_ && this->n_mass_ == a.n_mass_ && this->c_mass_ == a.c_mass_;
  }


  String FLASHHelperClasses::Tag::toString() const
  {
    String ret;
    if (n_mass_ >= 0)
      ret = '[' + std::to_string(n_mass_) + "]\t";
    ret += seq_;
    if (c_mass_ >= 0)
      ret += "\t[" + std::to_string(c_mass_) + "]";
    ret += "\tscore : " + std::to_string(getScore()) + "\tmzs : ";

    for (auto mz : mzs_)
    {
      ret += std::to_string(mz) + " ";
    }

    return ret;
  }
} // namespace OpenMS
