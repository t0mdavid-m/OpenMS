// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/config.h>
#include <iomanip>

namespace OpenMS
{
  /**
    @brief FLASHDeconv Spectrum level output *.tsv, *.msalign (for TopPIC) file formats
     @ingroup FileIO
**/

  class OPENMS_DLLAPI FLASHDeconvSpectrumFile
  {
  public:
    /**
            @brief write the header in the tsv output file (spectrum level)
            @param fs file stream to the output file
            @param ms_level ms level of the spectrum
            @param detail if set true, detailed information of the mass (e.g., peak list for the mass) is written
            @param report_decoy if set true, decoy and qvalue information will be written.
            @param write_abundance if set true (FLASHDeconv -quant), the MS1 header ends with the columns Abundance and AbundanceTraceCount.
       */
    static void writeDeconvolvedMassesHeader(std::fstream& fs,
                                             uint ms_level,
                                             bool detail,
                                             bool report_decoy,
                                             bool write_abundance = false);
    /**
          @brief write the deconvolved masses in the output file (spectrum level)
          @param dspec deconvolved spectrum to write
          @param fs file stream to the output file
          @param file_name the output file name that the deconvolved masses will be written.
          @param avg averagine information to calculate monoisotopic and average mass difference within this function. In PeakGroup (peaks of DeconvolvedSpectrum) only monoisotopic mass is recorded. To write both monoisotopic and average masses, their mass difference should be calculated using this averagine information.
          @param decoy_avg averagine for noise decoy
          @param tol mass tolerance
          @param write_detail if this is set, more detailed information on each mass will be written in the output file.
          @param report_decoy if set true, decoy and qvalue information will be written.
          @param noise_decoy_weight noise decoy weight of the spectrum's MS level: how often each noise decoy mass is written on average (deterministic)
          @param abundances (FLASHDeconv -quant) per mass of @p dspec (by its position): the abundance and the number of raw mass traces it is
          summed over, appended to each MS1 row (nan and 0 for a decoy). Without it, the rows are written as before.
          Default MS1 headers are:
            FileName, ScanNum, TargetDecoyType, RetentionTime, MassCountInSpec, AverageMass, MonoisotopicMass,
            SumIntensity, MinCharge, MaxCharge,
            PeakCount, IsotopeCosine, ChargeScore, MassSNR, ChargeSNR, RepresentativeCharge, RepresentativeMzStart, RepresentativeMzEnd, setQscore, PerChargeIntensity, PerIsotopeIntensity

      Default MS2 headers include MS1 headers plus:
        PrecursorScanNum, PrecursorMz, PrecursorIntensity, PrecursorCharge, PrecursorSNR, PrecursorMonoisotopicMass, PrecursorQscore

      Detailed MS1 and MS2 headers include all corresponding headers above plus:
        PeakMZs, PeakIntensities, PeakCharges, PeakMasses, PeakIsotopeIndices, PeakPPMErrors
    */
    static void writeDeconvolvedMasses(DeconvolvedSpectrum& dspec,
                                       std::fstream& fs,
                                       const String& file_name,
                                       const FLASHHelperClasses::PrecalculatedAveragine& avg,
                                       const FLASHHelperClasses::PrecalculatedAveragine& decoy_avg,
                                       double tol,
                                       bool write_detail,
                                       bool report_decoy, double noise_decoy_weight,
                                       const std::vector<std::pair<double, Size>>* abundances = nullptr);

    /// writeDeconvolvedMasses for every spectrum of @p spectra whose MS level has a file in @p out_files (open in @p streams), in input order and
    /// with the same bytes; the spectra are formatted in parallel (with @p report_decoy, the noise-decoy carries are replayed serially up front; @p noise_decoy_weights by MS level - 1).
    /// @p abundances (FLASHDeconv -quant) is indexed by the position of a spectrum in @p spectra, then of the mass in it
    static void writeDeconvolvedMasses(std::vector<DeconvolvedSpectrum>& spectra, std::vector<std::fstream>& streams, const std::vector<String>& out_files,
                                       const String& file_name, const FLASHHelperClasses::PrecalculatedAveragine& avg,
                                       const FLASHHelperClasses::PrecalculatedAveragine& decoy_avg, const std::vector<double>& tols, bool write_detail,
                                       bool report_decoy, const std::vector<double>& noise_decoy_weights,
                                       const std::vector<std::vector<std::pair<double, Size>>>* abundances = nullptr);

    /**
     *
     * @param map
     * @param deconvolved_spectra
     * @param deconvolved_mzML_file
     * @param annotated_mzML_file
     * @param mzml_charge
     * @param tols
     */
    static void writeMzML(const MSExperiment& map,
                          std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                          const String& deconvolved_mzML_file,
                          const String& annotated_mzML_file,
                          int mzml_charge, DoubleList tols);

    /// The deconvolved map writeMzML stores as the deconvolved mzML (the experimental settings of @p map and one spectrum per non-empty,
    /// non-decoy deconvolved spectrum with masses, in input order), built in memory, e.g. for FLASHTnT on a raw mzML.
    static void toDeconvolvedMap(const MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, int mzml_charge,
                                 const DoubleList& tols, MSExperiment& deconvolved_map);

    /**
     * write isobaric quantification results
     * @param fs
     * @param deconvolved_spectra
     */
    static void writeIsobaricQuantification(std::fstream& fs, std::vector<DeconvolvedSpectrum>& deconvolved_spectra);

    /// write topFD header
    static void writeTopFDHeader(std::fstream& fs, const Param& param);

    /**
      @brief write the deconvolved masses TopFD output (*.msalign)
      @param dspec deconvolved spectrum to write
      @param fs file stream to the output file
      @param filename mzml file name
      @param qval_threshold qvalue threshold to filter out high qvalue precursors.
      @param min_ms_level min ms level of the dataset
    */
    static void writeTopFD(DeconvolvedSpectrum& dspec, std::fstream& fs, const String& filename,
                           const double qval_threshold = 1.0,
                           const uint min_ms_level = 1);

  private:

    /// number of minimum peak count in topFD msalign file
    static const int topFD_min_peak_count_ = 3;
    /// number of maximum peak count in topFD msalign file
    static const int topFD_max_peak_count_ = 500;
  };
}// namespace OpenMS
