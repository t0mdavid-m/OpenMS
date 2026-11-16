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

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/config.h>

namespace OpenMS
{
  /**
   * @brief The charge check of a selected precursor (ADR-0078; CONTEXT.md "Charge check").
   *
   * The charge is read from the isotope spacing of the commanded isolation window on the decision survey. The window's profile
   * is rebuilt from the survey's centroids (one Gaussian per centroid, FWHM = m/z / R(m/z), R(m/z) = R200 * sqrt(200 / (m/z))),
   * autocorrelated, and a comb at the spacing 1.00235 / z, over one lag span shared by every candidate charge, picks the charge.
   * Pure functions: what a measurement changes is selection's decision (PrecursorSelection).
   */
  namespace ChargeCheck
  {
    /// The outcome of a check, as scan_commands.tsv's charge_check_status logs it. measure() returns Confirmed, Changed,
    /// UnverifiedResolution, UnverifiedLowScore, Inconsistent or Empty; selection turns a Changed into Swapped or Remeasured, or
    /// keeps the candidate as UnverifiedNoSpecies.
    enum class Status
    {
      NotChecked,           ///< no check ran on this command (a re-planned one). An empty cell means the check is off
      Confirmed,            ///< the measured charge is the commanded one
      Changed,              ///< a confident measurement 1..max_dz away: selection acts on it
      UnverifiedResolution, ///< isotope spacing / peak width below min_resolvability at the commanded charge
      UnverifiedLowScore,   ///< the comb score is below min_score
      UnverifiedNoSpecies,  ///< Changed, but nothing to swap to and no remeasure where the command is
      Inconsistent,         ///< a confident change of more than max_dz, or one the authored charges refuse
      Empty,                ///< no signal in the window
      Swapped,              ///< acquired as the survey's PeakGroup at the measured charge
      Remeasured            ///< acquired at the measured charge, its mass re-read at it
    };

    /// The authored part of precursor_selection.charge_check; every other parameter is fixed (spec 4.1).
    struct Settings
    {
      double min_score = 0.30;
      int max_dz = 3;
      double min_resolvability = 2.5;
    };

    /// What measure() read off the window. z_meas 0 = nothing measured.
    struct Measurement
    {
      Status status = Status::NotChecked;
      int z_meas = 0;
      double z_cont = 0.0;        ///< 1.00235 / the refined spacing
      double score = 0.0;         ///< the comb score of z_meas, 0..1
      double resolvability = 0.0; ///< (1.00235 / z_cmd) / FWHM at the commanded m/z
    };

    /// The averagine isotope spacing the comb uses (Da): the prototype's, which the replay fixtures pin.
    constexpr double ISO_SPACING = 1.00235;

    /// (1.00235 / @p z) / (m/z / R(m/z)); 0 for a non-positive charge, m/z or resolving power.
    OPENMS_DLLAPI double resolvability(double mz, int z, double r200);

    /// The comb alone on [@p lo, @p hi], over the candidate charges @p z_min..@p z_max: fills z_meas, z_cont and score, with
    /// status Empty when the window holds no signal and NotChecked otherwise. @p survey must be sorted by m/z.
    OPENMS_DLLAPI Measurement comb(const MSSpectrum& survey, double lo, double hi, double r200, int z_min, int z_max);

    /// The check of the window [@p lo, @p hi] commanded at @p mz and charge @p z_cmd, in the decision order of spec 5:
    /// resolution, empty, score, confirmed, changed (1..max_dz), inconsistent.
    OPENMS_DLLAPI Measurement measure(const MSSpectrum& survey, double lo, double hi, double mz, int z_cmd, double r200, int z_min,
                                      int z_max, const Settings& settings);

    /// The log vocabulary: "not_checked", "confirmed", "changed", "unverified_resolution", "unverified_low_score",
    /// "unverified_no_species", "inconsistent", "empty", "swapped", "remeasured".
    OPENMS_DLLAPI const char* toString(Status status);
  } // namespace ChargeCheck
} // namespace OpenMS
