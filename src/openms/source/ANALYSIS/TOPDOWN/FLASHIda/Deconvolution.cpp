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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/Precursor.h>

#include <cmath>

namespace OpenMS
{

  Deconvolution::Deconvolution(const Config& config)
  {
    Param sd_defaults = SpectralDeconvolution().getDefaults();
    sd_defaults.setValue("min_charge", config.deconvolution().min_charge);
    sd_defaults.setValue("max_charge", config.deconvolution().max_charge);
    sd_defaults.setValue("min_mass", config.deconvolution().min_mass);
    sd_defaults.setValue("max_mass", config.deconvolution().max_mass);
    DoubleList tol_values = {config.level(1).tolerance_ppm, config.level(2).tolerance_ppm};
    sd_defaults.setValue("tol", tol_values);
    fd_.setParameters(sd_defaults);
    fd_.calculateAveragine(false);
  }

  DeconvolvedSpectrum Deconvolution::deconvolveMS1(const double* mzs, const double* ints,
                                                    int length, double rt, double faims_cv)
  {
    auto spec = makeMSSpectrum_(mzs, ints, length, rt, 1, "ms1_spectrum");
    if (faims_cv != 0.0) { spec.setMetaValue("filter string", DataValue("cv=" + std::to_string(faims_cv))); }

    PeakGroup empty;
    fd_.performSpectrumDeconvolution(spec, 0, empty);
    deconvolved_spectrum_ = fd_.getDeconvolvedSpectrum();
    for (PeakGroup pg : deconvolved_spectrum_) {
      std::cout << pg.getMonoMass() << ", ";
    }
    std::cout << std::endl;
    return deconvolved_spectrum_;
  }

  int Deconvolution::deconvolveMSn(const double* mzs, const double* ints, int length,
                                    double rt, double precursor_mass, int precursor_charge)
  {
    // Clear previous state
    ms2_deconvolved_spectrum_.clear();
    ms2_deconv_valid_ = false;
    ms2_deconv_rt_ = rt;

    if (length == 0)
    {
      return 0;
    }

    // Create MSSpectrum from input
    auto spec = makeMSSpectrum_(mzs, ints, length, rt, 2, "ms2_spectrum");

    // Create precursor PeakGroup - only set if precursor_mass > 0 AND precursor_charge != 0
    PeakGroup precursor_pg;
    if (precursor_mass > 0 && precursor_charge != 0)
    {
      int abs_charge = std::abs(precursor_charge);
      bool is_positive = precursor_charge > 0;

      // Calculate precursor m/z from mass and charge
      double charge_mass = FLASHHelperClasses::getChargeMass(is_positive);
      double precursor_mz = (precursor_mass + abs_charge * charge_mass) / abs_charge;

      // Set precursor on the MSSpectrum (required for deconvolution mass range calculation)
      Precursor precursor;
      precursor.setMZ(precursor_mz);
      precursor.setCharge(precursor_charge);
      spec.getPrecursors().push_back(precursor);

      // Construct PeakGroup with proper charge range and polarity
      precursor_pg = PeakGroup(abs_charge, abs_charge, is_positive);
      precursor_pg.push_back(FLASHHelperClasses::LogMzPeak());
      precursor_pg.setMonoisotopicMass(precursor_mass);
      precursor_pg.setRepAbsCharge(abs_charge);
      precursor_pg.setQscore(1.0);  // Known precursor from MS1, high confidence
      precursor_pg.setSNR(1.0);
    }

    // Perform deconvolution (empty precursor_pg if mass <= 0 or charge == 0)
    fd_.performSpectrumDeconvolution(spec, 0, precursor_pg);
    ms2_deconvolved_spectrum_ = fd_.getDeconvolvedSpectrum();

    if (ms2_deconvolved_spectrum_.empty())
    {
      return 0;
    }

    // Sort by qscore (highest first) for getBestMS2Masses
    ms2_deconvolved_spectrum_.sortByQscore();
    ms2_deconv_valid_ = true;

    return static_cast<int>(ms2_deconvolved_spectrum_.size());
  }

  int Deconvolution::deconvolveMS2Py(const std::vector<double>& mzs,
                                      const std::vector<double>& ints,
                                      double rt, double precursor_mass, int precursor_charge)
  {
    if (mzs.empty() || mzs.size() != ints.size())
    {
      return 0;
    }
    return deconvolveMSn(mzs.data(), ints.data(), static_cast<int>(mzs.size()),
                         rt, precursor_mass, precursor_charge);
  }

  int Deconvolution::getMS2PeakGroupCount() const
  {
    return ms2_deconv_valid_ ? static_cast<int>(ms2_deconvolved_spectrum_.size()) : 0;
  }

  void Deconvolution::clearStoredMS2()
  {
    ms2_deconvolved_spectrum_.clear();
    ms2_deconv_valid_ = false;
    ms2_deconv_rt_ = -1.0;
  }

  MSSpectrum Deconvolution::makeMSSpectrum_(const double* mzs, const double* ints,
                                             int length, double rt, int ms_level,
                                             const char* name)
  {
    auto spec = MSSpectrum();
    for (int i = 0; i < length; i++)
    {
      if (ints[i] <= 0) { continue; }
      spec.emplace_back(mzs[i], ints[i]);
    }
    spec.setMSLevel(ms_level);
    spec.setName(name);
    spec.setRT(rt);
    return spec;
  }

} // namespace OpenMS
