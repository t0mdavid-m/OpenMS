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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>

namespace OpenMS
{

  // ---------------------------------------------------------------------------
  // ProteoformModel free methods
  // ---------------------------------------------------------------------------

  double ProteoformModel::coveragePct() const
  {
    return 0.0;
  }

  std::vector<double> ProteoformModel::combinedMs2FrameMasses() const
  {
    return {};
  }

  // ---------------------------------------------------------------------------
  // ProteoformTracker
  // ---------------------------------------------------------------------------

  ProteoformTracker::ProteoformTracker(const Config& cfg, IdaLogger& logger) :
    config_(cfg), logger_(logger)
  {
  }

  void ProteoformTracker::feedScan(int nominal_mass, uint8_t ms_level, const Ms2Params& params, int scan_id,
                                   const DeconvolvedSpectrum& deconv,
                                   const FragmentAnalysis::ProteoformMatch& match, double id_score)
  {
    ProteoformModel& m = models_[nominal_mass];
    m.nominal_mass = nominal_mass;

    PendingScan ps;
    ps.scan_id = scan_id;
    ps.ms_level = ms_level;
    ps.params = params;
    ps.match = match;
    ps.id_score = id_score;
    // Extract (monoisotopic mass, per-charge intensity) for every deconvolved peak group.
    // Iteration idiom matches IdaLogger / FragmentAnalysis: index into DeconvolvedSpectrum with [i].
    // Intensity: getChargeIntensity(getMaxIntensityAbsCharge()) — the max-charge intensity used
    // elsewhere as the sort key and precursor_intensity field (FragmentAnalysis.cpp:713).
    for (size_t i = 0; i < deconv.size(); ++i)
    {
      const PeakGroup& pg = deconv[i];
      ps.masses_intensities.emplace_back(pg.getMonoMass(),
                                         static_cast<double>(pg.getChargeIntensity(pg.getMaxIntensityAbsCharge())));
    }
    m.pending.push_back(std::move(ps));
  }

  void ProteoformTracker::finalize(int nominal_mass)
  {
    auto it = models_.find(nominal_mass);
    if (it != models_.end())
    {
      it->second.finalized = true;
    }
  }

  std::vector<ScanCommand> ProteoformTracker::planNextScans(int /*nominal_mass*/, ScanCommandQueue& /*queue*/)
  {
    return {};
  }

  const ProteoformModel* ProteoformTracker::model(int nominal_mass) const
  {
    auto it = models_.find(nominal_mass);
    if (it == models_.end())
    {
      return nullptr;
    }
    return &it->second;
  }

  void ProteoformTracker::mapScanOntoModel_(ProteoformModel& /*mdl*/, const PendingScan& /*scan*/)
  {
  }

  void ProteoformTracker::narrowModifications_(ProteoformModel& /*mdl*/)
  {
  }

  void ProteoformTracker::emitRow_(const ProteoformModel& /*mdl*/)
  {
  }

} // namespace OpenMS
