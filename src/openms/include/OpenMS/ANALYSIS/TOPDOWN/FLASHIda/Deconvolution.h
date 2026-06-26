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

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <vector>

namespace OpenMS
{

  /**
   * @brief Wraps SpectralDeconvolution and owns MS1/MS2 deconvolution results for FLASHIda.
   *
   * Responsible for:
   * - Building MSSpectrum from raw C arrays
   * - Running SpectralDeconvolution for MS1 and MS2 scans
   * - Storing and exposing the MS2 deconvolution result for downstream use
   *   (fragment analysis, MS3 targeting, exploration)
   */
  class OPENMS_DLLAPI Deconvolution
  {
  public:
    /// Constructor: initialise SpectralDeconvolution with explicit tolerance values
    /// @param tolerance_ppm_values PPM tolerances indexed by MS level (index 0 = MS1, 1 = MS2, ...)
    Deconvolution(const Config& config, const DoubleList& tolerance_ppm_values);

    /**
     * @brief Deconvolve an MS1 spectrum and return the result
     * @param mzs m/z values
     * @param ints intensities
     * @param length number of peaks
     * @param rt retention time in seconds
     * @param faims_cv FAIMS CV value (0.0 if FAIMS not enabled)
     * @return deconvolved spectrum (also accessible via deconvolvedMS1())
     */
    DeconvolvedSpectrum deconvolveMS1(const double* mzs, const double* ints, int length,
                                      double rt, double faims_cv);

    /**
     * @brief Deconvolve an MSn (n>1) spectrum, store result, return peak-group count
     * @param mzs m/z values
     * @param ints intensities
     * @param length number of peaks
     * @param rt retention time in seconds
     * @param precursor_mass precursor monoisotopic mass (<=0 means no precursor)
     * @param precursor_charge precursor charge (0 means no precursor)
     * @return number of peak groups found
     */
    int deconvolveMSn(const double* mzs, const double* ints, int length,
                      double rt, double precursor_mass, int precursor_charge);

    /**
     * @brief Python-friendly overload of deconvolveMS2
     */
    int deconvolveMS2Py(const std::vector<double>& mzs, const std::vector<double>& ints,
                        double rt, double precursor_mass, int precursor_charge);

    /// Access last MS1 deconvolution result, const (only valid after deconvolveMS1())
    const DeconvolvedSpectrum& deconvolvedMS1() const { return deconvolved_spectrum_; }

    /// Access last MS1 deconvolution result, mutable (for in-place sorting during filtering)
    DeconvolvedSpectrum& deconvolvedMS1() { return deconvolved_spectrum_; }

    /// Access stored MS2 deconvolution result (const)
    const DeconvolvedSpectrum& storedMS2() const { return ms2_deconvolved_spectrum_; }

    /// Access stored MS2 deconvolution result (non-const, e.g. for exploration feedResult)
    DeconvolvedSpectrum& storedMS2() { return ms2_deconvolved_spectrum_; }

    /// True if MS2 deconvolution is valid (i.e. deconvolveMS2() succeeded)
    bool hasStoredMS2() const { return ms2_deconv_valid_; }

    /// Number of peak groups in stored MS2 deconvolution (0 if not valid)
    int getMS2PeakGroupCount() const;

    /// Retention time of the stored MS2 deconvolution
    double storedMS2RT() const { return ms2_deconv_rt_; }

    /// Clear stored MS2 deconvolution state
    void clearStoredMS2();

    /**
     * @brief Access the underlying SpectralDeconvolution engine
     *
     * Used by getPeakGroups() / PrecursorSelection to call
     * fd_.setTargetMasses() before triggering a deconvolution.
     */
    SpectralDeconvolution& engine() { return fd_; }

  private:
    /**
     * @brief Build an MSSpectrum from raw C arrays
     * @param mzs m/z values
     * @param ints intensities
     * @param length number of peaks
     * @param rt retention time in seconds
     * @param ms_level MS level (1 or 2)
     * @param name spectrum name (may be nullptr)
     * @return populated MSSpectrum
     */
    static MSSpectrum makeMSSpectrum_(const double* mzs, const double* ints, int length,
                                      double rt, int ms_level, const char* name);

    /**
     * @brief Set the signed charge window on fd_ and re-apply parameters.
     *
     * Mutates the cached sd_param_ ("min_charge"/"max_charge") and pushes it to fd_
     * (which re-runs updateMembers_). Early-returns when the window already matches, so
     * back-to-back calls with the same window do not re-set parameters. The window is signed:
     * a negative sign selects negative-ion mode (SpectralDeconvolution derives polarity from
     * the sign of min_charge before taking absolute values).
     * @param signed_min_charge signed minimum charge (e.g. +1 or -1)
     * @param signed_max_charge signed maximum charge (e.g. +|precursor| or -|precursor|)
     */
    void setChargeWindow_(int signed_min_charge, int signed_max_charge);

    /// SpectralDeconvolution engine
    SpectralDeconvolution fd_;

    /// Cached full parameter set for fd_ (the global MS1 window); restored before each MS1 deconv.
    Param sd_param_;

    /// Global signed charge window (the MS1/default window) captured in the ctor.
    int global_min_charge_ = 0;
    int global_max_charge_ = 0;

    /// Currently-applied signed charge window on fd_ (avoids redundant setParameters calls).
    int cur_min_charge_ = 0;
    int cur_max_charge_ = 0;

    /// Last MS1 deconvolution result
    DeconvolvedSpectrum deconvolved_spectrum_{0};

    /// Stored MS2 deconvolution result
    DeconvolvedSpectrum ms2_deconvolved_spectrum_{0};

    /// Retention time of stored MS2 deconvolution
    double ms2_deconv_rt_ = -1.0;

    /// Whether stored MS2 deconvolution is valid
    bool ms2_deconv_valid_ = false;
  };

} // namespace OpenMS
