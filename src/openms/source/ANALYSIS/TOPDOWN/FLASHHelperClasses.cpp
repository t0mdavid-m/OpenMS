// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <exception>
#include <typeinfo>
#include <utility>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
  FLASHHelperClasses::PrecalculatedAveragine::PrecalculatedAveragine(const double min_mass, const double max_mass, const double delta, CoarseIsotopePatternGenerator& generator,
                                                                     const bool use_RNA_averagine, const double decoy_iso_distance) :
      mass_interval_(delta),
      min_mass_(min_mass)
  {
    // the mass grid
    std::vector<double> masses;
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
      masses.push_back(mass);
    }

    const Size bin_count = masses.size();
    apex_index_.resize(bin_count);
    right_count_from_apex_.resize(bin_count);
    left_count_from_apex_.resize(bin_count);
    average_mono_mass_difference_.resize(bin_count);
    abundant_mono_mass_difference_.resize(bin_count);
    isotopes_.resize(bin_count);
    snr_mul_factor_.resize(bin_count);
    std::vector<std::exception_ptr> errors(bin_count);

    // The bins are independent of each other (the left and right counts are made cumulative below), so they are calculated in parallel,
    // each thread with its own copy of the generator.
    const auto calculate_bin = [&](const Size index, CoarseIsotopePatternGenerator& local_generator) {
      const double mass = masses[index];
      auto iso = use_RNA_averagine ? local_generator.estimateFromRNAMonoWeight(mass) : local_generator.estimateFromPeptideMonoWeight(mass);

      if (decoy_iso_distance > 0)
      {
        auto decoy_iso(iso);

        for (Size k = 0; k < iso.size(); k++)
        {
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
      double pwr = 0;
      while (iso.size() - trim_count > min_iso_length && left_count < right_count)
      {
        double lint = iso[left_count].getIntensity();
        double rint = iso[right_count].getIntensity();

        bool trim_left = true;
        if (lint < rint)
        {
          pwr += lint * lint;
        }
        else
        {
          pwr += rint * rint;
          trim_left = false;
        }
        if (total_pwr - pwr < total_pwr * min_pwr)
        {
          break;
        }

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
      total_pwr -= pwr;

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

      apex_index_[index] = most_abundant_index_;
      right_count_from_apex_[index] = right_count;
      left_count_from_apex_[index] = left_count;
      average_mono_mass_difference_[index] = iso.averageMass() - iso[0].getMZ();
      abundant_mono_mass_difference_[index] = iso.getMostAbundant().getMZ() - iso[0].getMZ();
      isotopes_[index] = iso;
      snr_mul_factor_[index] = intensity_sum * intensity_sum;
    };

    // a generator of a derived class would be sliced by a copy: it runs serially on the caller's generator
    if (typeid(generator) != typeid(CoarseIsotopePatternGenerator))
    {
      for (Size index = 0; index < bin_count; ++index)
      {
        calculate_bin(index, generator);
      }
    }
    else
    {
      int worker_count = 1;
#ifdef _OPENMP
      if (! omp_in_parallel()) { worker_count = std::max(1, (int)std::min(bin_count, (Size)omp_get_max_threads())); }
#endif
      std::vector<CoarseIsotopePatternGenerator> generators(worker_count, generator); // copied here: no allocation inside the region
#pragma omp parallel for num_threads(worker_count) schedule(static) if (worker_count > 1)
      for (SignedSize index = 0; index < (SignedSize)bin_count; ++index)
      {
        int worker = 0;
#ifdef _OPENMP
        worker = omp_get_thread_num();
#endif
        try
        {
          calculate_bin((Size)index, generators[worker]);
        }
        catch (...)
        {
          errors[index] = std::current_exception();
        }
      }
    }

    // the first error in mass order; the left and right counts are the maxima up to each mass
    int max_left_count = 0;
    int max_right_count = 0;
    for (Size index = 0; index < bin_count; ++index)
    {
      if (errors[index]) { std::rethrow_exception(errors[index]); }
      max_left_count = std::max(max_left_count, left_count_from_apex_[index]);
      max_right_count = std::max(max_right_count, right_count_from_apex_[index]);
      left_count_from_apex_[index] = max_left_count;
      right_count_from_apex_[index] = max_right_count;
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
  FLASHHelperClasses::Tag::Tag(const String& seq, double n_mass, double c_mass, const std::vector<double>& mzs, const std::vector<int>& scores,  int scan) :
      seq_(std::move(seq)), n_mass_(n_mass), c_mass_(c_mass), mzs_(mzs), scores_(scores), scan_(scan), length_(mzs.size() - 1)
  {
    upper_seq_ = seq_.toUpper();
  }

  FLASHHelperClasses::Tag::Tag(String seq, double n_mass, double c_mass, std::vector<double>& mzs, std::vector<int>& scores, int scan, Size length) :
      seq_(std::move(seq)), n_mass_(n_mass), c_mass_(c_mass), mzs_(mzs), scores_(scores), scan_(scan), length_(length)
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
      {
        if (this->c_mass_ >= 0 && a.c_mass_ >= 0) // c term tag
        {
          return this->c_mass_ < a.c_mass_;
        }
        else if (this->n_mass_ >= 0 && a.n_mass_ >= 0)
        {
          return this->n_mass_ < a.n_mass_;
        }
        else
        {
          return this->n_mass_ < a.n_mass_;
        }
      }
      return this->seq_ < a.seq_;
    }
    return this->length_ < a.length_;
  }

  bool FLASHHelperClasses::Tag::operator>(const Tag& a) const
  {
    if (this->length_ == a.length_)
    {
      if (this->seq_ == a.seq_)
      {
        if (this->c_mass_ >= 0 && a.c_mass_ >= 0) // c term tag
        {
          return this->c_mass_ > a.c_mass_;
        }
        else if (this->n_mass_ >= 0 && a.n_mass_ >= 0)
        {
          return this->n_mass_ > a.n_mass_;
        }
        else
        {
          return this->n_mass_ > a.n_mass_;
        }
      }
      return this->seq_ > a.seq_;
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

  FLASHHelperClasses::MassTag::MassTag(const String& full_seq, const std::vector<int>& positions, int scan)
  {
    if (positions.size() < 2) return;
    seq_ = std::move(full_seq);
    positions_ = positions;
    scan = scan_;
    //is_n_term_ = positions[0] < positions[1];


// String seq_, upper_seq_;
    //       float rt_ = -1;
    //       int scan_;
    //       std::vector<int> positions_;
    //       int index_;
    //       Size length_;

  }

  bool FLASHHelperClasses::MassTag::operator<(const MassTag& a) const
  {
    (void)a;
    return false;
  }

  bool FLASHHelperClasses::MassTag::operator>(const MassTag& a) const
  {
    (void)a;
    return false;
  }

  bool FLASHHelperClasses::MassTag::operator==(const MassTag& a) const
  {
    (void)a;
    return false;
  }

    const String& FLASHHelperClasses::MassTag::getSequence() const
    {
      return seq_;
    }

    Size FLASHHelperClasses::MassTag::getLength() const
    {
      return length_;
    }

    int FLASHHelperClasses::MassTag::getScan() const
    {
      return scan_;
    }

  } // namespace OpenMS