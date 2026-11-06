// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <map>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/KERNEL/Peak1D.h>

namespace OpenMS
{
  class PeakGroup;

  /**
@brief   Qvalue : contains functions to calculate Qvalues from deconvolution quality score (Qscore) defined in PeakGroupScoring
@ingroup Topdown
*/

  class OPENMS_DLLAPI Qvalue
  {
  public:
    typedef FLASHHelperClasses::LogMzPeak LogMzPeak;

    /// Calculate and perform a batch update of peak group qvalues using Qscores of target and dummy peak groups in deconvolved spectra, when FDR report is necessary.
    /// A Qscore2D without a target of the same score (a decoy, or a mass of a feature counted once) gets the q-value of the lowest target
    /// score at or above it (the targets counted at that threshold are the same), the highest target score's if it is above all.
    /// @param deconvolved_spectra target and decoy deconvolved spectra
    /// @return the noise decoy weight of each MS level (by MS level) for decoy output
    static std::map<uint, double> updatePeakGroupQvalues(std::vector<DeconvolvedSpectrum>& deconvolved_spectra);

  private:
  };
} // namespace OpenMS