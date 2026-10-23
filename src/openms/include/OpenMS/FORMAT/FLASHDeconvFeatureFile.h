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
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/config.h>
#include <iomanip>
#include <map>
#include <iostream>


namespace OpenMS
{
  /**
    @brief FLASHDeconv feature level output *.tsv, *.ms1ft (for Promex), *.feature (for TopPIC) file formats
     @ingroup FileIO
**/

  class OPENMS_DLLAPI FLASHDeconvFeatureFile
  {
  public:
    /// write header line for regular file output; with @p write_abundance, the columns Abundance and AbundanceTraceCount (FLASHDeconv -quant) are appended
    static void writeHeader(std::fstream& fs, bool report_decoy = false, bool write_abundance = false);

    /// write header line for topFD feature file
    static void writeTopFDFeatureHeader(std::fstream& fs, uint ms_level);

    /// write the features in regular file output; with @p abundances (feature index -> abundance and the number of raw mass traces it is
    /// summed over, FLASHDeconv -quant), the two columns are appended to each row (nan and 0 for a feature without an entry)
    static void writeFeatures(const std::vector<FLASHHelperClasses::MassFeature>& mass_features, const String& file_name, std::fstream& fs, bool report_decoy = false,
                              const std::map<uint, std::pair<double, Size>>* abundances = nullptr);

    /**
     * @brief Find mass features and write features in TopFD format files.
     * @param deconvolved_spectra deconvolved spectra - feature indices are updated only for TopFD and TopPIC outputs
     * @param mass_features mass features to be written
     * @param scan_rt_map scan number to retention time map
     * @param file_name input spectrum file name
     * @param fs file stream
     * @param ms_level ms level
     */

    static void writeTopFDFeatures(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const std::vector<FLASHHelperClasses::MassFeature>& mass_features,
                                   const std::map<int, double>& scan_rt_map, const String& file_name, std::fstream& fs, uint ms_level);

  };
} // namespace OpenMS