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
#include <map>
#include <set>
#include <string>
#include <vector>

namespace OpenMS
{

  // Anonymous namespace for ProteoformTracker helper functions
  namespace
  {
    // Tighten the UPPER boundary of an ambiguous range [rs, re] (region 1-based) down to
    // new_upper, justified by intensity I. support_lower / support_upper carry the intensities
    // behind each boundary. A move is a CONFLICT iff it would push the upper past the lower
    // (new_upper < rs, inverting the range): on conflict it only wins if I beats the opposing
    // support (support_lower), in which case it overrides and collapses the range, re-justifying
    // BOTH boundaries with I. A non-conflicting move applies whenever it is strictly inward.
    inline void tightenUpper(int& rs, int& re, double& support_lower, double& support_upper, int new_upper, double I)
    {
      if (new_upper >= re) return;  // not inward -> nothing to do
      if (new_upper >= rs)
      {
        // No conflict: stays at or above the lower boundary.
        re = new_upper;
        support_upper = I;
      }
      else if (I > support_lower)
      {
        // Conflict: would cross the lower boundary; new evidence wins -> collapse & re-justify.
        rs = new_upper;
        re = new_upper;
        support_lower = I;
        support_upper = I;
      }
      // else: conflicting, weaker evidence -> keep the higher-support boundary (no move).
    }

    // Symmetric counterpart: tighten the LOWER boundary up to new_lower; opposing = support_upper.
    inline void tightenLower(int& rs, int& re, double& support_lower, double& support_upper, int new_lower, double I)
    {
      if (new_lower <= rs) return;  // not inward
      if (new_lower <= re)
      {
        // No conflict: stays at or below the upper boundary.
        rs = new_lower;
        support_lower = I;
      }
      else if (I > support_upper)
      {
        rs = new_lower;
        re = new_lower;
        support_lower = I;
        support_upper = I;
      }
    }
  } // anonymous namespace

  // ---------------------------------------------------------------------------
  // ProteoformModel free methods
  // ---------------------------------------------------------------------------

  double ProteoformModel::coveragePct() const
  {
    const int L = (region_start < 0) ? static_cast<int>(proteoform_sequence.size()) : (region_end - region_start);
    if (L <= 0) return 0.0;
    if (fragments.empty()) return 0.0;

    // Collect all [cover_start, cover_end] intervals (1-based inclusive), clamped to [1, L].
    std::vector<std::pair<int, int>> intervals;
    intervals.reserve(fragments.size());
    for (const auto& kv : fragments)
    {
      const MappedFragment& f = kv.second;
      const int cs = std::max(f.cover_start, 1);
      const int ce = std::min(f.cover_end, L);
      if (cs <= ce) intervals.emplace_back(cs, ce);
    }
    if (intervals.empty()) return 0.0;

    // Sort by start, then merge overlapping intervals, sum covered length.
    std::sort(intervals.begin(), intervals.end());
    int covered = 0;
    int cur_start = intervals[0].first;
    int cur_end = intervals[0].second;
    for (size_t i = 1; i < intervals.size(); ++i)
    {
      if (intervals[i].first <= cur_end + 1)
      {
        // Overlapping or adjacent: extend.
        cur_end = std::max(cur_end, intervals[i].second);
      }
      else
      {
        covered += cur_end - cur_start + 1;
        cur_start = intervals[i].first;
        cur_end = intervals[i].second;
      }
    }
    covered += cur_end - cur_start + 1;

    return covered / static_cast<double>(L);
  }

  std::vector<double> ProteoformModel::combinedMs2FrameMasses() const
  {
    std::vector<double> masses;
    masses.reserve(fragments.size() * 2);
    for (const auto& kv : fragments)
    {
      const MappedFragment& f = kv.second;
      if (f.best_ms2.has_value()) masses.push_back(f.best_ms2->observed_mass);
      if (f.best_ms3.has_value()) masses.push_back(f.best_ms3->observed_mass);
    }
    std::sort(masses.begin(), masses.end());
    return masses;
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
                                   const FragmentAnalysis::ProteoformMatch& match, double id_score,
                                   const ScanCommand& ms2_ctx)
  {
    ProteoformModel& m = models_[nominal_mass];
    m.nominal_mass = nominal_mass;

    // Capture the MS2 command context once (used by planNextScans/buildMS3 in the next task).
    if (!m.has_ms2_ctx)
    {
      m.ms2_ctx = ms2_ctx;
      m.has_ms2_ctx = true;
    }

    PendingScan ps;
    ps.scan_id = scan_id;
    ps.ms_level = ms_level;
    ps.params = params;
    ps.match = match;
    ps.id_score = id_score;
    // Extract one PeakRecord per deconvolved peak group.
    // mz idiom: getMzRange(charge) centre — matches buildMS2 (ScanCommandQueue.cpp) and the
    // Exploration.cpp precursor_mz derivation, so the value is ready for direct use in buildMS3.
    // Intensity: getChargeIntensity(getMaxIntensityAbsCharge()) — max-charge intensity used
    // elsewhere as the sort key and precursor_intensity field (FragmentAnalysis.cpp:713).
    for (size_t i = 0; i < deconv.size(); ++i)
    {
      const PeakGroup& pg = deconv[i];
      const int charge = pg.getMaxIntensityAbsCharge();
      auto [mz1, mz2] = pg.getMzRange(charge);
      PeakRecord pr;
      pr.mono_mass = pg.getMonoMass();
      pr.mz = (mz1 + mz2) / 2.0;
      pr.charge = charge;
      pr.intensity = static_cast<double>(pg.getChargeIntensity(charge));
      ps.peaks.push_back(pr);
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

      // 3) Recover intensity (and mz/charge for MS3 targeting): closest mass within per-level
      //    tolerance; fall back to the overall closest mass (never drop a matched fragment for lack
      //    of a peak — the fragment was already validated as a match upstream).
      double intensity = 0.0;
      double matched_mz = 0.0;
      int matched_charge = 0;
      {
        double best_diff = std::numeric_limits<double>::max();
        bool within_tol = false;
        for (const PeakRecord& pr : ps.peaks)
        {
          const double diff = std::abs(pr.mono_mass - obs_mass);
          const double tol_abs = obs_mass * tol_ppm * 1e-6;
          const bool in_tol = (diff <= tol_abs);
          // Prefer the closest in-tolerance entry; fall back to the closest overall.
          if (in_tol && !within_tol)
          {
            within_tol = true;
            best_diff = diff;
            intensity = pr.intensity;
            matched_mz = pr.mz;
            matched_charge = pr.charge;
          }
          else if (in_tol == within_tol && diff < best_diff)
          {
            best_diff = diff;
            intensity = pr.intensity;
            matched_mz = pr.mz;
            matched_charge = pr.charge;
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
      obs.frag_mz = matched_mz;
      obs.frag_charge = matched_charge;

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

  void ProteoformTracker::narrowModifications_(ProteoformModel& m)
  {
    // --- Faithful per-fragment bracketing -------------------------------------------------------
    // For each ambiguous modification (candidate_start < candidate_end), shrink the localization
    // range using bracketing MappedFragments. A fragment brackets the range iff its backbone
    // cleavage falls strictly inside it (it covers some, but not all, of the ambiguous residues).
    // We test "PTM shift INSIDE this fragment's coverage" (base + mass_shift) vs "OUTSIDE" (base)
    // against the fragment's observed mass at the per-level tolerance, then tighten the boundary the
    // fragment constrains. Bidirectional intensity support: each boundary records the intensity that
    // justified it; a later fragment that conflicts (would push a boundary past the opposing one)
    // only wins if its intensity beats the opposing support (conflict -> higher-intensity-wins).
    //
    // FRAME ALIGNMENT (verified against T5 mapScanOntoModel_ + MS3FragmentMatcher::computeEquivalentIon):
    //   m.modifications[i].candidate_start/end : FULL-PROTEIN 1-based residue range (seeded in T5).
    //   m.fragments keys / cover_start / cover_end / ion_index : WINNER-REGION 1-based frame (T5).
    //   Part-A theoreticals are computed on the WINNER-REGION substring, so baseline[type][k-1] is
    //   the region-frame theoretical of ion (type, region-index k) -- the same frame as the keys.
    //
    //   A full-protein 1-based residue position `p` maps to the region 1-based position `p - ws`
    //   (ws = max(region_start, 0)): residue p is at region 0-based (p-1)-ws, i.e. region 1-based
    //   (p-1)-ws+1 = p-ws. Identity reduction: region == full => ws = 0 => region_pos == full_pos.
    //   We tighten in the region frame, then convert the final boundaries back to full-protein
    //   (full_pos = region_pos + ws) before storing, so the model keeps full-protein coordinates
    //   consistent with how T5 seeded them and how the log/ProForma read them.
    // -------------------------------------------------------------------------------------------
    const int P = static_cast<int>(m.proteoform_sequence.size());
    if (P <= 0 || m.modifications.empty() || m.fragments.empty()) return;

    // Winner region resolved to [ws, we) (0-based, exclusive end); L = region residue count.
    const int ws = (m.region_start < 0) ? 0 : m.region_start;
    const int we = (m.region_end < 0) ? P : m.region_end;
    const int L = we - ws;
    if (L <= 0) return;

    // The winner-region substring all theoreticals are computed on.
    const String region_seq =
      (m.region_start < 0) ? String(m.proteoform_sequence) : String(m.proteoform_sequence).substr(ws, L);

    // Distinct ion-type chars actually observed in the mapped fragments (as strings), for Part A.
    std::vector<std::string> ion_types;
    {
      std::set<char> seen;
      for (const auto& kv : m.fragments)
      {
        const MappedFragment& f = kv.second;
        if (f.ion_type.empty()) continue;
        seen.insert(f.ion_type[0]);
      }
      for (char c : seen) ion_types.emplace_back(1, c);
    }
    if (ion_types.empty()) return;

    // Process each ambiguous modification independently.
    for (size_t mi = 0; mi < m.modifications.size(); ++mi)
    {
      ModificationState& mod = m.modifications[mi];
      if (mod.candidate_start >= mod.candidate_end) continue;  // already localized / invalid

      // Baseline = theoreticals of the proteoform with ALL OTHER mods applied, WITHOUT this mod.
      // PTM positions handed to Part A must be in region-1-based coordinates.
      std::vector<FragmentAnalysis::PTMSite> other_ptms;
      for (size_t oi = 0; oi < m.modifications.size(); ++oi)
      {
        if (oi == mi) continue;
        const ModificationState& om = m.modifications[oi];
        const int os_region = om.candidate_start - ws;
        const int oe_region = om.candidate_end - ws;
        FragmentAnalysis::PTMSite site;
        site.start_position = os_region;
        site.end_position = oe_region;
        site.position = (os_region + oe_region) / 2;
        site.mass_shift = om.mass_shift;
        other_ptms.push_back(site);
      }
      std::map<char, std::vector<double>> baseline;
      FragmentAnalysis::computePTMAdjustedFragmentMasses(region_seq, other_ptms, ion_types, baseline);

      // Ambiguous range in the region frame (1-based inclusive).
      int rs_region = mod.candidate_start - ws;
      int re_region = mod.candidate_end - ws;

      for (const auto& kv : m.fragments)
      {
        const MappedFragment& f = kv.second;
        if (f.ion_type.empty() || f.ion_index <= 0) continue;

        // Bracketing test (region frame). f.cover_* are region-1-based inclusive (T5):
        //   prefix b_k covers [1, k]      (k = f.ion_index = f.cover_end): brackets iff rs <= k < re
        //   suffix y_k covers [L-k+1, L]  (c = f.cover_start = L-k+1):     brackets iff rs < c <= re
        if (f.is_prefix)
        {
          if (!(rs_region <= f.cover_end && f.cover_end < re_region)) continue;
        }
        else
        {
          if (!(rs_region < f.cover_start && f.cover_start <= re_region)) continue;
        }

        // Best observation: prefer MS2, else MS3. Skip if neither present.
        const FragmentObservation* obs =
          f.best_ms2.has_value() ? &(*f.best_ms2) : (f.best_ms3.has_value() ? &(*f.best_ms3) : nullptr);
        if (obs == nullptr) continue;

        // Theoretical "without" this mod (region frame, 1-based index -> 0-based vector).
        auto bit = baseline.find(f.ion_type[0]);
        if (bit == baseline.end()) continue;
        const std::vector<double>& masses = bit->second;
        const int vi = f.ion_index - 1;
        if (vi < 0 || vi >= static_cast<int>(masses.size())) continue;
        const double base = masses[vi];
        const double with = base + mod.mass_shift;

        const double tol_abs = obs->observed_mass * config_.level(obs->ms_level).tolerance_ppm * 1e-6;
        const bool matches_with = std::abs(obs->observed_mass - with) <= tol_abs;
        const bool matches_without = std::abs(obs->observed_mass - base) <= tol_abs;

        // Ambiguous fragment (matches both or neither) -> no information, skip.
        if (matches_with == matches_without) continue;

        const double I = obs->intensity;

        if (matches_with)
        {
          // PTM shift IS inside this fragment's coverage.
          //   prefix covering [1, k]: PTM is on residues rs..k -> tighten UPPER to k = f.cover_end.
          //   suffix covering [c, L]: PTM is on residues c..re -> tighten LOWER to c = f.cover_start.
          if (f.is_prefix)
            tightenUpper(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_end, I);
          else
            tightenLower(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_start, I);
        }
        else
        {
          // PTM shift is OUTSIDE this fragment's coverage (complementary tighten).
          //   prefix covering [1, k]: PTM is on residues k+1..re -> tighten LOWER to k+1.
          //   suffix covering [c, L]: PTM is on residues rs..c-1 -> tighten UPPER to c-1.
          if (f.is_prefix)
            tightenLower(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_end + 1, I);
          else
            tightenUpper(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_start - 1, I);
        }
      }

      // Convert the narrowed boundaries back to FULL-PROTEIN before storing.
      mod.candidate_start = rs_region + ws;
      mod.candidate_end = re_region + ws;
    }
  }

  void ProteoformTracker::emitRow_(const ProteoformModel& /*mdl*/)
  {
  }

} // namespace OpenMS
