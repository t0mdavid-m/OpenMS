// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Quantification.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

namespace OpenMS
{

  Quantification::Quantification(const Config& config) :
      config_(config)
  {
  }

  bool Quantification::isDifferentiallyAbundant(const double* mzs,
                                                 const double* ints,
                                                 const int length,
                                                 const double rt,
                                                 const int ms_level,
                                                 const char* name,
                                                 double reporter_mz_tol,
                                                 double fold_change_threshold,
                                                 bool only_one_condition)
  {
    // Create spectrum
    MSSpectrum spec;
    for (int i = 0; i < length; i++)
    {
      if (ints[i] > 0) { spec.emplace_back(mzs[i], ints[i]); }
    }
    spec.setMSLevel(ms_level);
    spec.setName(name);
    spec.setRT(rt);

    // Set precursor with HCD activation - neccessary for channel extractor
    OpenMS::Precursor precursor;
    precursor.setActivationMethods({OpenMS::Precursor::ActivationMethod::HCD});
    spec.setPrecursors({precursor});

    // Create experiment
    MSExperiment exp;
    exp.addSpectrum(spec);

    // TODO: Variable channel extractors
    TMTSixPlexQuantitationMethod quant_method;
    IsobaricChannelExtractor channel_extractor(&quant_method);

    // Set parameters
    Param p = channel_extractor.getParameters();
    p.setValue("reporter_mass_shift", reporter_mz_tol);
    channel_extractor.setParameters(p);

    // Extract reporter-ion channels into a ConsensusMap
    ConsensusMap consensus_map_raw;
    channel_extractor.extractChannels(exp, consensus_map_raw);

    // Collect m/z-intensity pairs from the ConsensusFeatures
    std::vector<std::pair<double, double>> mz_int;
    mz_int.reserve(quant_method.getNumberOfChannels());
    for (const auto& cf : consensus_map_raw)
    {
        for (auto& i : cf)
        {
            mz_int.emplace_back(i.getMZ(), i.getIntensity());
        }
    }

    // Something went wrong – bail out early
    if (mz_int.size() != quant_method.getNumberOfChannels())
    {
        std::cout << "FLASHIda - channel extraction failed..." << std::endl;
        return false;
    }

    // Sort by m/z to ensure channel order
    std::sort(mz_int.begin(), mz_int.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    // Extract intensities
    std::vector<double> intensities;
    intensities.reserve(mz_int.size());
    for (const auto& p2 : mz_int) intensities.push_back(p2.second);

    for (auto intensity : intensities) {
      std::cout << intensity << std::endl;
    }

    // TODO: Make configurable
    if (only_one_condition) {
      bool first_sample_present = std::none_of(
        intensities.begin(), intensities.begin()+3, [](double x){ return x < 1e-3; }
      );
      std::cout << first_sample_present << std::endl;
      bool second_sample_present = std::none_of(
        intensities.begin()+3, intensities.end(), [](double x){ return x < 1e-3; }
      );
      std::cout << second_sample_present << std::endl;
      bool first_sample_missing = std::all_of(
        intensities.begin(), intensities.begin()+3, [](double x){ return x < 1e-3; }
      );
      std::cout << first_sample_missing << std::endl;
      bool second_sample_missing = std::all_of(
        intensities.begin()+3, intensities.end(), [](double x){ return x < 1e-3; }
      );
      std::cout << second_sample_missing << std::endl;
      // No signal
      if (first_sample_missing && second_sample_missing) {
        return false;
      }
      // Both are incomplete
      else if (!first_sample_present && !second_sample_present) {
        return false;
      }
      // One signal
      else if (
        (first_sample_missing || second_sample_missing)
        && (first_sample_present || second_sample_present)
      ) {
        return true;
      }
    }

    // Reject spectra with missing / too-low reporter peaks
    if (!only_one_condition && std::any_of(intensities.begin(), intensities.end(),
                    [](double x){ return x < 1e-3; }))
    {
        return false;
    }

    const double sample1_mean = std::accumulate(intensities.begin(),
                                                intensities.begin() + 3, 0.0) / 3.0;

    const double sample2_mean = std::accumulate(intensities.begin() + 3,
                                                intensities.end(), 0.0) / 3.0;

    const double fold_change = sample1_mean / sample2_mean;

    return (fold_change > fold_change_threshold) || ((1/fold_change) > fold_change_threshold);
  }

} // namespace OpenMS
