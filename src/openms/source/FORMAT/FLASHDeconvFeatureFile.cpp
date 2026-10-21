// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FLASHDeconvFeatureFile.h>
#include <algorithm>
#include <exception>
#include <iomanip>

namespace OpenMS
{
  /**
    @brief FLASHDeconv Spectrum level output *.tsv, *.msalign (for TopPIC) file formats
     @ingroup FileIO
**/

  void FLASHDeconvFeatureFile::writeHeader(std::fstream& fs, bool report_decoy)
  {
    fs << "FeatureIndex\tFileName\tMSLevel";
    if (report_decoy) fs << "\tIsDecoy";

    fs << "\tMonoisotopicMass\tAverageMass\tMassCount\tStartRetentionTime"
          "\tEndRetentionTime\tRetentionTimeDuration\tApexRetentionTime"
          "\tSumIntensity\tMaxIntensity\tFeatureQuantity\tMinCharge\tMaxCharge\tChargeCount\tIsotopeCosineScore\tQscore2D\tPerChargeIntensity\tPerIsotopeIntensity"
          "\n";
  }

  void FLASHDeconvFeatureFile::writeTopFDFeatureHeader(std::fstream& fs, uint ms_level)
  {
    //  //File_name	Fraction_ID	Spectrum_ID	Scans	MS_one_ID	MS_one_scans	Fraction_feature_ID	Fraction_feature_intensity
    // Fraction_feature_score	Fraction_feature_min_time	Fraction_feature_max_time
    // Fraction_feature_apex_time	Precursor_neutral_monoisotopic_mass	Precursor_monoisotopic_mz	Precursor_average_mz	Precursor_charge	Precursor_intensity
    // (TopPIC 2.x layout; TopPIC 2.x reads the columns by position and crashes on the 1.x layout without Precursor_neutral_monoisotopic_mass)

    if (ms_level == 1)
    {
      fs << "File_name\tFraction_ID\tFeature_ID\tMass\tIntensity\tMin_time\tMax_time\tMin_scan\tMax_scan\tMin_charge\tMax_charge\tApex_time\tApex_scan\tApex_intensity\tRep_charge\tRep_average_mz\tEnvelope_num\tEC_score\n";
    }
    else
    {
      fs << "File_name\tFraction_ID\tSpectrum_ID\tScans\tMS_one_ID\tMS_one_scans\tFraction_feature_ID\tFraction_feature_intensity\tFraction_feature_score\tFraction_feature_min_time\tFraction_feature_max_time\tFraction_feature_apex_time\tPrecursor_neutral_monoisotopic_mass\tPrecursor_monoisotopic_mz\tPrecursor_average_mz\tPrecursor_charge\tPrecursor_intensity\n";
    }
    fs.flush();
  }

  void FLASHDeconvFeatureFile::writeFeatures(const std::vector<FLASHHelperClasses::MassFeature>& mass_features, const String& file_name, std::fstream& fs, bool report_decoy)
  {
    // the row of a feature
    const auto write_row = [&file_name, report_decoy](const FLASHHelperClasses::MassFeature& mass_feature, std::stringstream& ss) {
      const auto& mt = mass_feature.mt;
      double mass = mt.getCentroidMZ() + mass_feature.iso_offset * Constants::ISOTOPE_MASSDIFF_55K_U;
      double avg_mass = mass_feature.avg_mass;
      double sum_intensity = .0;

      for (auto& p : mt)
      {
        sum_intensity += p.getIntensity();
      }

      ss << mass_feature.index << "\t" << file_name << "\t" << mass_feature.ms_level;

      if (report_decoy)
      { ss << "\t" << (mass_feature.is_decoy? 1 : 0);
      }

      ss << "\t" << std::to_string(mass) << "\t" << std::to_string(avg_mass) << "\t" // massdiff
         << mt.getSize() << "\t" << mt.begin()->getRT()/60.0 << "\t" << mt.rbegin()->getRT()/60.0 << "\t" << mt.getTraceLength() << "\t" << mt[mt.findMaxByIntPeak()].getRT()/60.0 << "\t" << sum_intensity << "\t"
         << mt.getMaxIntensity(false) << "\t" << mt.computePeakArea() << "\t" << mass_feature.min_charge << "\t" << mass_feature.max_charge << "\t" << mass_feature.charge_count << "\t"
         << mass_feature.isotope_score << "\t" << std::setprecision (15) << (mass_feature.qscore) << std::setprecision (-1) << "\t";

      for (int i = mass_feature.min_charge; i <= mass_feature.max_charge; i++)
      {
        ss << mass_feature.per_charge_intensity[abs(i)];
        if (i < mass_feature.max_charge)
        { ss << ";";
        }
      }

      ss << "\t";
      int iso_end_index = 0;

      for (Size i = 0; i < mass_feature.per_isotope_intensity.size(); i++)
      {
        if (mass_feature.per_isotope_intensity[i] == 0)
        {
          continue;
        }
        iso_end_index = (int)i;
      }
      for (int i = 0; i <= iso_end_index; i++)
      {
        ss << mass_feature.per_isotope_intensity[i];

        if (i < iso_end_index)
        { ss << ";";
        }
      }
      ss << "\n";
    };

    // The rows are formatted in blocks in parallel and written in order. A block starts with the stream state the rows before it leave:
    // the default precision before the first row, and the precision a row ends with (-1) before any other row.
    const Size block_size = 256;
    const Size block_count = (mass_features.size() + block_size - 1) / block_size;
    std::vector<std::string> blocks(block_count);
    std::vector<std::exception_ptr> errors(block_count);
#pragma omp parallel for schedule(dynamic, 1)
    for (SignedSize b = 0; b < (SignedSize)block_count; b++)
    {
      try
      {
        std::stringstream ss;
        if (b > 0) ss << std::setprecision(-1);
        const Size row_end = std::min(mass_features.size(), (Size)(b + 1) * block_size);
        for (Size row = (Size)b * block_size; row < row_end; row++)
        {
          write_row(mass_features[row], ss);
        }
        blocks[b] = ss.str();
      }
      catch (...)
      {
        errors[b] = std::current_exception();
      }
    }
    // the first error (of the first row that failed), and then nothing written, as before
    for (const auto& error : errors)
    {
      if (error) std::rethrow_exception(error);
    }
    for (const auto& block : blocks)
    {
      fs << block;
    }
  }

  void FLASHDeconvFeatureFile::writeTopFDFeatures(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const std::vector<FLASHHelperClasses::MassFeature>& mass_features,
                                                  const std::map<int, double>& scan_rt_map, const String& file_name, std::fstream& fs, uint ms_level)
  {
    std::stringstream ss;
    static std::map<int, uint> msNscan_to_feature_index;

    if (ms_level == 1)
    {
      uint max_feature_index = 0;
      for (const auto& mass_feature : mass_features)
      {
        if (mass_feature.ms_level != 1 || mass_feature.is_decoy) continue;
        double sum_intensity = .0;

        for (const auto& m : mass_feature.mt)
        {
          sum_intensity += m.getIntensity();
        }

        const auto& apex = mass_feature.mt[mass_feature.mt.findMaxByIntPeak()];
        ss << file_name << "\t0\t" << mass_feature.index << "\t" << std::to_string(mass_feature.mt.getCentroidMZ()) << "\t" << std::to_string(sum_intensity) << "\t" << std::to_string(mass_feature.mt.begin()->getRT()/60.0)
           << "\t" << std::to_string(mass_feature.mt.rbegin()->getRT()/60.0) << "\t" << mass_feature.min_scan_number << "\t" << mass_feature.max_scan_number << "\t"
           << mass_feature.min_charge << "\t" << mass_feature.max_charge << "\t" << std::to_string(apex.getRT()/60.0) << "\t" << mass_feature.scan_number << "\t"
           << std::to_string(apex.getIntensity()) << "\t" << mass_feature.rep_charge << "\t" << mass_feature.rep_mz << "\t0\t0\n";

        max_feature_index = std::max(max_feature_index, mass_feature.index);
      }
      for (auto& dspec : deconvolved_spectra)
      {
        if (dspec.getOriginalSpectrum().getMSLevel() == 1 || dspec.getPrecursorPeakGroup().empty()) continue;
        if (dspec.getPrecursorPeakGroup().getFeatureIndex() != 0) continue;
        auto pg = dspec.getPrecursorPeakGroup();

        double rt = scan_rt_map.at(pg.getScanNumber())/60.0;
        const auto& [z, Z] = pg.getAbsChargeRange();
        pg.setFeatureIndex(++max_feature_index);
        dspec.setPrecursorPeakGroup(pg);

        double rep_mz = 0;
        double p_int = 0;

        for (const auto& p : pg)
        {
          if (p.abs_charge != pg.getRepAbsCharge()) continue;
          if (p.intensity < p_int) continue;
          p_int = p.intensity;
          rep_mz = p.mz;
        }
        ss << file_name << "\t0\t" << max_feature_index << "\t" << std::to_string(pg.getMonoMass()) << "\t" << std::to_string(pg.getIntensity()) << "\t" << std::to_string(rt)
           << "\t" << std::to_string(rt) << "\t" << pg.getScanNumber() << "\t" << pg.getScanNumber() << "\t"
           << z << "\t" << Z << "\t" << std::to_string(rt) << "\t" << pg.getScanNumber() << "\t"
           << std::to_string(pg.getIntensity()) << "\t" << pg.getRepAbsCharge() << "\t" << std::to_string(rep_mz) << "\t0\t0\n";
      }
    }

    if (ms_level == 2)
    {
      // The MS1 target features by their index. A precursor without one got an index beyond them in the MS1 pass above, which may be
      // the index of a decoy feature (report_FDR) or an MSn feature as well, so a feature is looked up by its index among these only.
      std::map<uint, const FLASHHelperClasses::MassFeature*> ms1_features;
      auto mono_mz = [](double mono_mass, int charge) { return charge == 0 ? .0 : mono_mass / charge + Constants::PROTON_MASS_U; };
      for (const auto& mass_feature : mass_features)
      {
        if (mass_feature.ms_level != 1 || mass_feature.is_decoy) continue;
        ms1_features[mass_feature.index] = &mass_feature;
      }
      for (auto& dspec : deconvolved_spectra)
      {
        if (dspec.getOriginalSpectrum().getMSLevel() == 1 || dspec.getPrecursorPeakGroup().empty()) continue;
        if (dspec.getPrecursorPeakGroup().getFeatureIndex() == 0) continue;
        auto pg = dspec.getPrecursorPeakGroup();
        int ms2_scan_number = dspec.getScanNumber();
        ss << file_name << "\t0\t" << ms2_scan_number - 1 << "\t" << ms2_scan_number << "\t" << pg.getScanNumber() - 1 << "\t" << pg.getScanNumber() << "\t" << pg.getFeatureIndex() << "\t";

        const auto fit = ms1_features.find(pg.getFeatureIndex());
        if (fit != ms1_features.end())
        {
          const auto& mass_feature = *fit->second;
          const auto& apex = mass_feature.mt[mass_feature.mt.findMaxByIntPeak()];
          double sum_intensity = .0;

          for (const auto& m : mass_feature.mt)
          {
            sum_intensity += m.getIntensity();
          }
          ss << std::to_string(sum_intensity) << "\t0\t" << std::to_string(mass_feature.mt.begin()->getRT()/60.0)
             << "\t" << std::to_string(mass_feature.mt.rbegin()->getRT()/60.0) << "\t" << std::to_string(apex.getRT()/60.0) << "\t" << std::to_string(mass_feature.mt.getCentroidMZ())
             << "\t" << std::to_string(mono_mz(mass_feature.mt.getCentroidMZ(), mass_feature.rep_charge)) << "\t" << std::to_string(mass_feature.rep_mz) << "\t" << mass_feature.rep_charge << "\t" << std::to_string(dspec.getPrecursor().getIntensity()) << "\n";
        }
        else
        {
          double rt = scan_rt_map.at(pg.getScanNumber()) / 60.0;
          double rep_mz = 0;
          double p_int = 0;

          for (const auto& p : pg)
          {
            if (p.abs_charge != pg.getRepAbsCharge()) continue;
            if (p.intensity < p_int) continue;
            p_int = p.intensity;
            rep_mz = p.mz;
          }

          ss << std::to_string(pg.getIntensity()) << "\t0\t" << std::to_string(rt)
             << "\t" << std::to_string(rt) << "\t" << std::to_string(rt) << "\t" << std::to_string(pg.getMonoMass())
             << "\t" << std::to_string(mono_mz(pg.getMonoMass(), pg.getRepAbsCharge())) << "\t" << std::to_string(rep_mz) << "\t" << pg.getRepAbsCharge() << "\t" << std::to_string(dspec.getPrecursor().getIntensity()) << "\n";
        }
      }
    }
    fs << ss.str();
  }
} // namespace OpenMS