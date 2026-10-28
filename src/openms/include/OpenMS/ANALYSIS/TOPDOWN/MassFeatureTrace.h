// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <iomanip>
#include <iostream>

namespace OpenMS
{
  /**
  @brief Feature trace in mass dimension for FLASHDeconv
  This class performs mass tracing on the deconvolved masses by SpectralDeconvolution
  In other words, per spectrum deconvolved masses are converted into deconvolved features
  Currently only works for MS1 spectra. (Top-down DIA is not yet used much).
  Every time an MS1 spectrum is deconvolved, the relevant information is stored in this class.
  This class also comes with tsv, TopFD feature file formats.
  @ingroup Topdown
  */

  class OPENMS_DLLAPI MassFeatureTrace : public DefaultParamHandler
  {
  public:
    typedef FLASHHelperClasses::PrecalculatedAveragine PrecalculatedAveragine;
    typedef FLASHHelperClasses::LogMzPeak LogMzPeak;

    /// constructor
    MassFeatureTrace();

    /// destructor
    ~MassFeatureTrace() override = default;

    /// copy constructor
    MassFeatureTrace(const MassFeatureTrace&) = default;

    /// move constructor
    MassFeatureTrace(MassFeatureTrace&& other) = default;

    /// assignment operator
    MassFeatureTrace& operator=(const MassFeatureTrace& fd) = default;
    MassFeatureTrace& operator=(MassFeatureTrace&& fd) = default;

    /**
       @brief Find mass features.
       @param averagine precalculated averagine for cosine calculation
       @param deconvolved_spectra the spectra on which features are found
       @param ms_level ms level to process
       @param is_decoy if set, only process decoy spectra. otherwise only target spectra
       */
    std::vector<FLASHHelperClasses::MassFeature> findFeaturesAndUpdateQscore2D(const PrecalculatedAveragine& averagine, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, int ms_level = 1,
                                                                                     bool is_decoy = false);

    /// The features of one tracing before they are numbered: per feature, the peak groups it labels (index of the spectrum in the traced
    /// vector, index of the peak group in the spectrum) and its 2D Qscore
    struct TracedFeatures
    {
      std::vector<FLASHHelperClasses::MassFeature> features;
      std::vector<std::vector<std::pair<Size, Size>>> peak_groups;
      std::vector<double> qscores;
    };

    /// restart the feature indices at 1 (FLASHDeconvAlgorithm::run calls it per run; the indices continued over runs in one process)
    static void resetFeatureIndex();

    /// The tracing of findFeaturesAndUpdateQscore2D without numbering the features and labelling their peak groups (reads the spectra only)
    TracedFeatures traceFeatures(const PrecalculatedAveragine& averagine, const std::vector<DeconvolvedSpectrum>& deconvolved_spectra, int ms_level = 1,
                                 bool is_decoy = false) const;

    /**
       @brief Numbers traced features with the process-wide feature index, in their order, and labels their peak groups (feature index and
       2D Qscore), which findFeaturesAndUpdateQscore2D does feature by feature.
       @param traced the traced features (moved out)
       @param deconvolved_spectra the traced spectra, or the vector they were taken from when @p spectrum_indices is given
       @param spectrum_indices if not empty, the index in @p deconvolved_spectra of each traced spectrum
       */
    static std::vector<FLASHHelperClasses::MassFeature> numberFeatures(TracedFeatures& traced, std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                                                       const std::vector<Size>& spectrum_indices = {});

  protected:
    void updateMembers_() override;

  private:
    /// peak group information is stored in here for tracing
    std::map<double, std::map<double, PeakGroup>> peak_group_map_; // rt , mono mass, peakgroup
  };
} // namespace OpenMS