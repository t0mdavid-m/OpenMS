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

#include <cmath>
#include <limits>

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
    if (it == models_.end()) return;
    ProteoformModel& m = it->second;
    if (m.pending.empty()) return;

    // 1) Pick the winner: the staged scan with the highest match.score that actually identified a
    //    proteoform (score >= 0 and a non-empty sequence). Ties keep the first-seen (strictly-greater).
    const PendingScan* win = nullptr;
    for (const PendingScan& ps : m.pending)
    {
      if (ps.match.score < 0 || ps.match.proteoform_sequence.empty()) continue;
      if (win == nullptr || ps.match.score > win->match.score) win = &ps;
    }

    if (win == nullptr)
    {
      // No scan identified anything -> keep the model empty (no header, no row, no plan).
      m.pending.clear();
      m.finalized = true;
      return;
    }

    // 2) Seed the model header from the winner.
    m.proteoform_sequence = win->match.proteoform_sequence;
    m.region_start = win->match.region_start;
    m.region_end = win->match.region_end;
    m.identification_score = win->match.score;
    m.winner_scan_id = win->scan_id;

    // 3) Seed modifications from the winner's PTM sites (full-protein 1-based ranges).
    m.modifications.clear();
    for (const FragmentAnalysis::PTMSite& s : win->match.ptm_sites)
    {
      m.modifications.push_back(ModificationState{s.mass_shift, s.start_position, s.end_position, 0.0, 0.0});
    }

    // 4) Map EVERY staged scan's already-matched fragments onto the model (Strategy A).
    m.fragments.clear();
    for (const PendingScan& ps : m.pending)
    {
      mapScanOntoModel_(m, ps);
    }

    // 5) Narrowing (T6) and row emission (T11) are still stubs.
    narrowModifications_(m);
    ++m.update_index;
    emitRow_(m);

    // 6) Done.
    m.pending.clear();
    m.finalized = true;
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

  void ProteoformTracker::mapScanOntoModel_(ProteoformModel& m, const PendingScan& ps)
  {
    // --- Strategy A ----------------------------------------------------------------------------
    // Pool the scan's ALREADY-matched fragments (ps.match.fragments) onto the winner proteoform.
    // No theoretical masses are recomputed here. The hard part is the coordinate-frame conversion:
    // every MS2/MS3 fragment, regardless of which (possibly truncated) region its own scan matched
    // against, is reduced to ONE backbone-cleavage position expressed in the WINNER region frame so
    // the same cleavage from any scan/level collides into the same MappedFragment.
    //
    // Frame definitions (all derived from the verified source semantics):
    //   region_start (rs/ws): 0-based inclusive start, -1 => full sequence (=> 0).
    //   region_end   (re/we): 0-based EXCLUSIVE end,   -1 => full sequence (=> P).
    //   P = full-protein length = winner proteoform_sequence length.
    //
    //   MS2 ion_index `idx` is 1-based in THAT SCAN's matched region:
    //     prefix b_idx covers region residues 1..idx  -> full-protein prefix cut = rs + idx
    //     suffix y_idx covers the last idx region res. -> full-protein suffix cut = (P - re) + idx
    //   (These match MS3FragmentMatcher::computeEquivalentIon exactly: prefix `start + i`,
    //    suffix `P - end + i`. MS3 equiv_index is therefore ALREADY in the full-protein frame.)
    //
    //   Winner-region frame index (the map key index):
    //     prefix:  winner_idx = full_prefix_cut - ws
    //     suffix:  winner_idx = full_suffix_cut - (P - we)
    //   Identity reduction: scan-region == winner-region == full => rs=ws=0, re=we=P =>
    //     prefix winner_idx = (0 + idx) - 0 = idx; suffix winner_idx = ((P-P) + idx) - 0 = idx. OK.
    // -------------------------------------------------------------------------------------------

    const int P = static_cast<int>(m.proteoform_sequence.size());
    if (P <= 0) return;

    // Winner region resolved to [ws, we) (0-based, exclusive end).
    const int ws = (m.region_start < 0) ? 0 : m.region_start;
    const int we = (m.region_end < 0) ? P : m.region_end;
    const int L = we - ws; // winner-region residue count
    if (L <= 0) return;

    // This scan's matched region resolved to [rs, re) (0-based, exclusive end).
    const int rs = (ps.match.region_start < 0) ? 0 : ps.match.region_start;
    const int re = (ps.match.region_end < 0) ? P : ps.match.region_end;

    const double tol_ppm = config_.level(ps.ms_level).tolerance_ppm;

    for (const FragmentAnalysis::ProteoformMatch::FragmentMatch& fm : ps.match.fragments)
    {
      // 1) Pick the ion identity by level.
      std::string type;
      int idx = 0;
      double obs_mass = 0.0;
      if (ps.ms_level == 3)
      {
        type = fm.equiv_type;       // already full-protein b/y
        idx = fm.equiv_index;       // already full-protein index
        obs_mass = fm.adjusted_mass;
      }
      else
      {
        type = fm.ion_type;         // region-relative
        idx = fm.ion_index;         // region-relative 1-based
        obs_mass = fm.observed_mass;
      }
      if (type.empty() || idx <= 0) continue;

      const char ic = type[0];
      const bool is_prefix = (ic == 'a' || ic == 'b' || ic == 'c');
      const bool is_suffix = (ic == 'x' || ic == 'y' || ic == 'z');
      if (!is_prefix && !is_suffix) continue;

      // 2) Convert to the full-protein cleavage, then to the winner-region frame.
      int winner_idx;
      if (ps.ms_level == 3)
      {
        // MS3 equiv_index is full-protein; map straight into the winner region.
        winner_idx = is_prefix ? (idx - ws) : (idx - (P - we));
      }
      else
      {
        // MS2 idx is relative to this scan's region -> full-protein -> winner region.
        const int full_cut = is_prefix ? (rs + idx) : ((P - re) + idx);
        winner_idx = is_prefix ? (full_cut - ws) : (full_cut - (P - we));
      }

      // DROP fragments whose cleavage falls outside the winner region.
      if (winner_idx < 1 || winner_idx > L) continue;

      // 3) Recover intensity: closest mass within per-level tolerance; otherwise fall back to the
      //    overall closest mass (documented choice -- never drop a matched fragment for lack of an
      //    intensity peak, since the fragment was already validated as a match upstream).
      double intensity = 0.0;
      {
        double best_diff = std::numeric_limits<double>::max();
        bool within_tol = false;
        for (const std::pair<double, double>& mi : ps.masses_intensities)
        {
          const double diff = std::abs(mi.first - obs_mass);
          const double tol_abs = obs_mass * tol_ppm * 1e-6;
          const bool in_tol = (diff <= tol_abs);
          // Prefer the closest in-tolerance entry; fall back to the closest overall.
          if (in_tol && !within_tol)
          {
            within_tol = true;
            best_diff = diff;
            intensity = mi.second;
          }
          else if (in_tol == within_tol && diff < best_diff)
          {
            best_diff = diff;
            intensity = mi.second;
          }
        }
      }

      // 4) Upsert the MappedFragment keyed by (ion_type, winner_frame_index).
      const FragmentKey key{type, winner_idx};
      MappedFragment& mfrag = m.fragments[key];
      mfrag.ion_type = type;
      mfrag.ion_index = winner_idx;
      mfrag.is_prefix = is_prefix;
      // Coverage = the winner-region residue span this ion covers (1-based inclusive, winner frame).
      if (is_prefix)
      {
        mfrag.cover_start = 1;
        mfrag.cover_end = winner_idx;
      }
      else
      {
        mfrag.cover_start = L - winner_idx + 1;
        mfrag.cover_end = L;
      }
      // theoretical_mass is left 0 here; T6 computes theoreticals.

      FragmentObservation obs;
      obs.ms_level = ps.ms_level;
      obs.observed_mass = obs_mass;
      obs.intensity = intensity;
      obs.source_scan_id = ps.scan_id;
      obs.params = ps.params;

      if (ps.ms_level == 3)
      {
        if (!mfrag.best_ms3.has_value() || obs.intensity > mfrag.best_ms3->intensity)
          mfrag.best_ms3 = obs;
        ++mfrag.n_ms3;
      }
      else
      {
        if (!mfrag.best_ms2.has_value() || obs.intensity > mfrag.best_ms2->intensity)
          mfrag.best_ms2 = obs;
        ++mfrag.n_ms2;
      }
    }
  }

  void ProteoformTracker::narrowModifications_(ProteoformModel& /*mdl*/)
  {
  }

  void ProteoformTracker::emitRow_(const ProteoformModel& /*mdl*/)
  {
  }

} // namespace OpenMS
