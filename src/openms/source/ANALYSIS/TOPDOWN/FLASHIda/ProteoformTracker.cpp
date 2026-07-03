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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>  // getNominalMass (read-only) for the pooled nominal_mass column

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
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

    // Does fragment @p f bracket the (region-frame, 1-based inclusive) range [rs, re]? A fragment
    // brackets the range iff its backbone cleavage falls strictly inside it (it covers SOME but not
    // ALL of the ambiguous residues). Identical predicate used by narrowModifications_:
    //   prefix b_k covers [1, k]      (k = f.cover_end):   brackets iff rs <= k < re
    //   suffix y_k covers [L-k+1, L]  (c = f.cover_start):  brackets iff rs < c <= re
    inline bool fragmentBrackets(const MappedFragment& f, int rs, int re)
    {
      if (f.ion_type.empty() || f.ion_index <= 0) return false;
      if (f.is_prefix) return (rs <= f.cover_end && f.cover_end < re);
      return (rs < f.cover_start && f.cover_start <= re);
    }

    // Does fragment @p f fully CONTAIN the (region-frame, 1-based inclusive) range [rs, re]? A
    // fragment contains the range iff its coverage spans all of it: cover_start <= rs && cover_end >= re.
    // (This is the dispatch predicate for MS3 targeting -- the whole ambiguous range must lie inside the
    // fragment's coverage so an MS3 re-feed of that fragment carries internal cleavages able to narrow it.)
    //   prefix b_k covers [1, k]      (k = f.cover_end):   contains iff k >= re   (cover_start == 1 <= rs)
    //   suffix y_k covers [L-k+1, L]  (c = f.cover_start): contains iff c <= rs   (cover_end  == L >= re)
    inline bool fragmentContains(const MappedFragment& f, int rs, int re)
    {
      if (f.ion_type.empty() || f.ion_index <= 0) return false;
      return (f.cover_start <= rs && f.cover_end >= re);
    }
  } // anonymous namespace

  // ---------------------------------------------------------------------------
  // ProteoformModel free methods
  // ---------------------------------------------------------------------------

  double ProteoformModel::coveragePct() const
  {
    // Backbone-cleavage sequence coverage: fraction of the L-1 inter-residue bonds witnessed by >=1
    // fragment (NOT a residue-span union, which saturates to 1.0 from one long prefix + one long suffix).
    // Each mapped fragment reduces to exactly one cleavage in the winner-region frame: a prefix ion
    // covering [1,k] witnesses the bond after residue k (= cover_end); a suffix ion covering [c,L]
    // witnesses the bond before residue c (= cover_start-1).
    const int L = (region_start < 0) ? static_cast<int>(proteoform_sequence.size()) : (region_end - region_start);
    if (L < 2 || fragments.empty()) return 0.0;
    std::set<int> sites;   // distinct backbone-bond indices observed, in [1, L-1]
    for (const auto& kv : fragments)
    {
      const MappedFragment& f = kv.second;
      const int site = f.is_prefix ? f.cover_end : (f.cover_start - 1);
      if (site >= 1 && site <= L - 1) sites.insert(site);   // whole-protein prefix (site==L) / suffix (site==0) excluded
    }
    return static_cast<double>(sites.size()) / static_cast<double>(L - 1);
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

  void ProteoformTracker::feedScan(int precursor_id, uint8_t ms_level, const Ms2Params& params, int scan_id,
                                   const DeconvolvedSpectrum& deconv,
                                   const FragmentAnalysis::ProteoformMatch& match, double id_score,
                                   const ScanCommand& ms2_ctx)
  {
    ProteoformModel& m = models_[precursor_id];
    m.precursor_id = precursor_id;
    m.contributing_scan_ids.insert(scan_id);   // cumulative: every fed scan stays (never dropped on supersede)

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
      // Isolation-window span = the PeakGroup's m/z extent at the selected charge (getMzRange centre-to-
      // edge span). This sizes the CE-sweep sub-window in Exploration.cpp:832. It is NOT the legacy
      // (wend - wstart) isolation width; ScanCommandQueue.cpp:352 floors it at 2.0 Th before emission.
      pr.iso_width = mz2 - mz1;
      pr.stage1_scores = FragmentAnalysis::FragmentScores::fromPeakGroup(pg, charge);
      ps.peaks.push_back(pr);
    }
    m.pending.push_back(std::move(ps));
  }

  void ProteoformTracker::finalize(int precursor_id)
  {
    auto it = models_.find(precursor_id);
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
    emitRow_(m, "MS2", ScanCommandQueue::encode(m.winner_scan_id));

    // 6) Done.
    m.pending.clear();
    m.finalized = true;
  }

  void ProteoformTracker::foldMs3(int precursor_id, const std::string& trigger_ion,
                                  const std::string& trigger_scan_id)
  {
    auto it = models_.find(precursor_id);
    if (it == models_.end()) return;
    ProteoformModel& m = it->second;
    // Fold onto an already-identified MS2 baseline ONLY. No baseline -> nothing to refine; never reset.
    if (!m.finalized || m.proteoform_sequence.empty()) { m.pending.clear(); return; }
    if (m.pending.empty()) return;
    // Additive: map the new MS3 scan(s) onto the EXISTING fragments/winner. Do NOT clear, do NOT re-pick.
    for (const PendingScan& ps : m.pending) mapScanOntoModel_(m, ps);
    narrowModifications_(m);
    ++m.update_index;
    emitRow_(m, trigger_ion, trigger_scan_id);
    m.pending.clear();
  }

  std::vector<Ms3Target> ProteoformTracker::planNextScans(int precursor_id)
  {
    auto it = models_.find(precursor_id);
    if (it == models_.end()) return {};
    ProteoformModel& m = it->second;
    // No identified model / no captured MS2 context -> no plan (ADR-0002).
    if (m.proteoform_sequence.empty() || !m.has_ms2_ctx) return {};

    const CharacterizationObjective obj = config_.characterization().objective;
    const int budget = config_.level(2).max_targets;  // reuse the MS2 max_targets as the MS3 budget
    if (budget <= 0) return {};

    const int P = static_cast<int>(m.proteoform_sequence.size());
    if (P <= 0) return {};
    const int ws = (m.region_start < 0) ? 0 : m.region_start;
    const int we = (m.region_end < 0) ? P : m.region_end;
    const int L = we - ws;  // winner-region residue count (fragment frame)
    if (L <= 0) return {};

    // --- 1) Choose ordered target fragments by objective ----------------------------------------
    // Each chosen fragment MUST carry a best-MS2 observation (the MS3 isolation descriptors come
    // from it). Dedup by FragmentKey. Bounded by budget.
    std::vector<const MappedFragment*> targets;
    std::set<FragmentKey> chosen;

    auto add_target = [&](const MappedFragment* f) -> bool {
      if (f == nullptr || !f->best_ms2.has_value()) return false;
      const FragmentKey key{f->ion_type, f->ion_index};
      if (chosen.count(key)) return false;
      chosen.insert(key);
      targets.push_back(f);
      return static_cast<int>(targets.size()) < budget;  // false once the budget is full
    };

    if (obj == CharacterizationObjective::Ambiguity)
    {
      // For each still-ambiguous modification (widest range first), build the list of best-MS2
      // fragments that fully CONTAIN that range (region frame), strongest intensity first. A
      // containing fragment is the only useful MS3 target: re-fragmenting it produces internal
      // cleavages that can narrow the range. (A merely-bracketing fragment would already have
      // narrowed it, so for a genuinely ambiguous range none exists -> we must use containers.)
      // modifications are full-protein 1-based ranges; convert to the region frame via -ws to match
      // the fragment keys, exactly as narrowModifications_ does.
      std::vector<const ModificationState*> ambiguous;
      for (const ModificationState& mod : m.modifications)
        if (mod.candidate_start < mod.candidate_end) ambiguous.push_back(&mod);
      std::sort(ambiguous.begin(), ambiguous.end(),
                [](const ModificationState* a, const ModificationState* b) {
                  return (a->candidate_end - a->candidate_start) > (b->candidate_end - b->candidate_start);
                });

      // Per-mod candidate container lists, each sorted by best-MS2 intensity descending.
      std::vector<std::vector<const MappedFragment*>> per_mod_candidates;
      per_mod_candidates.reserve(ambiguous.size());
      for (const ModificationState* mod : ambiguous)
      {
        const int rs_region = mod->candidate_start - ws;
        const int re_region = mod->candidate_end - ws;
        std::vector<const MappedFragment*> cands;
        for (const auto& kv : m.fragments)
        {
          const MappedFragment& f = kv.second;
          if (!f.best_ms2.has_value()) continue;
          if (!fragmentContains(f, rs_region, re_region)) continue;
          cands.push_back(&f);
        }
        std::sort(cands.begin(), cands.end(),
                  [](const MappedFragment* a, const MappedFragment* b) {
                    return a->best_ms2->intensity > b->best_ms2->intensity;
                  });
        per_mod_candidates.push_back(std::move(cands));
      }

      // Round-robin: pass 1 gives every mod its strongest container, then continue round-robin with
      // each mod's next-best, until the budget fills. add_target dedups by FragmentKey across mods and
      // returns false once the budget is full. A mod with no container contributes nothing.
      bool budget_full = false;
      for (size_t rank = 0; !budget_full; ++rank)
      {
        bool any_at_rank = false;
        for (const auto& cands : per_mod_candidates)
        {
          if (rank >= cands.size()) continue;
          any_at_rank = true;
          const FragmentKey key{cands[rank]->ion_type, cands[rank]->ion_index};
          if (chosen.count(key)) continue;  // already added for an earlier mod
          if (!add_target(cands[rank])) { budget_full = true; break; }
        }
        if (!any_at_rank) break;  // exhausted every mod's candidate list
      }
    }
    else  // CharacterizationObjective::Coverage
    {
      // Find the largest uncovered residue gaps (complement of the merged coverage intervals in the
      // winner-region frame), then for each gap pick the strongest best-MS2 fragment that fully
      // CONTAINS it (its coverage spans the whole gap), so an MS3 re-feed of that fragment yields the
      // internal cleavages needed to cover the gap. Uses span-interval merging (appropriate for
      // gap finding); coveragePct() uses a distinct backbone-cleavage metric instead.
      std::vector<std::pair<int, int>> intervals;  // covered [cs, ce], region-1-based inclusive, clamped to [1, L]
      for (const auto& kv : m.fragments)
      {
        const MappedFragment& f = kv.second;
        const int cs = std::max(f.cover_start, 1);
        const int ce = std::min(f.cover_end, L);
        if (cs <= ce) intervals.emplace_back(cs, ce);
      }
      std::sort(intervals.begin(), intervals.end());

      // Build the complement (uncovered) gaps over [1, L].
      std::vector<std::pair<int, int>> gaps;
      int cursor = 1;
      for (const auto& iv : intervals)
      {
        if (iv.first > cursor) gaps.emplace_back(cursor, iv.first - 1);
        cursor = std::max(cursor, iv.second + 1);
      }
      if (cursor <= L) gaps.emplace_back(cursor, L);
      // Largest gaps first.
      std::sort(gaps.begin(), gaps.end(),
                [](const std::pair<int, int>& a, const std::pair<int, int>& b) {
                  return (a.second - a.first) > (b.second - b.first);
                });

      for (const auto& gap : gaps)
      {
        const MappedFragment* best = nullptr;
        for (const auto& kv : m.fragments)
        {
          const MappedFragment& f = kv.second;
          if (!f.best_ms2.has_value()) continue;
          if (chosen.count(FragmentKey{f.ion_type, f.ion_index})) continue;
          // The fragment's coverage must fully span the gap to be a useful MS3 re-feed target.
          if (!fragmentContains(f, gap.first, gap.second)) continue;
          if (best == nullptr || f.best_ms2->intensity > best->best_ms2->intensity) best = &f;
        }
        if (best != nullptr && !add_target(best)) break;
      }
    }

    if (targets.empty()) return {};

    // Strongest best-MS2 first within the chosen set (the executor dispatches in this order).
    std::sort(targets.begin(), targets.end(),
              [](const MappedFragment* a, const MappedFragment* b) {
                return a->best_ms2->intensity > b->best_ms2->intensity;
              });

    // --- 2) Emit Ms3Targets, budget-bounded ---------------------------------------------------------
    // OFF (default): one target per fragment from best_ms2 — byte-identical to M1.
    // ON (ms3_all_charges): one target per observed charge state of each fragment, strongest first;
    // stop as soon as the total emitted count reaches budget.
    const bool all_charges = config_.characterization().ms3_all_charges;
    std::vector<Ms3Target> out;
    for (const MappedFragment* f : targets)
    {
      if (static_cast<int>(out.size()) >= budget) break;
      // Build the per-charge observation list to emit for this fragment.
      std::vector<const FragmentObservation*> obs_list;
      if (all_charges && !f->ms2_by_charge.empty())
      {
        for (const auto& kv : f->ms2_by_charge) obs_list.push_back(&kv.second);
        std::sort(obs_list.begin(), obs_list.end(),
                  [](const FragmentObservation* a, const FragmentObservation* b) {
                    return a->intensity > b->intensity;
                  });
      }
      else
      {
        obs_list.push_back(&(*f->best_ms2));  // single best charge (default)
      }
      for (const FragmentObservation* o : obs_list)
      {
        if (static_cast<int>(out.size()) >= budget) break;
        Ms3Target t;
        t.ion_type = f->ion_type;
        t.ion_index = f->ion_index;
        t.frag_mz = o->frag_mz;
        t.frag_charge = o->frag_charge;
        t.frag_mass = o->observed_mass;
        t.iso_width = o->iso_width;
        t.stage0_params = o->params;
        t.stage1_scores = o->stage1_scores;
        out.push_back(std::move(t));
      }
    }
    return out;
  }

  const ProteoformModel* ProteoformTracker::model(int precursor_id) const
  {
    auto it = models_.find(precursor_id);
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
      double matched_iso_width = 0.0;
      FragmentAnalysis::FragmentScores matched_stage1;
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
            matched_iso_width = pr.iso_width;
            matched_stage1 = pr.stage1_scores;
          }
          else if (in_tol == within_tol && diff < best_diff)
          {
            best_diff = diff;
            intensity = pr.intensity;
            matched_mz = pr.mz;
            matched_charge = pr.charge;
            matched_iso_width = pr.iso_width;
            matched_stage1 = pr.stage1_scores;
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
      obs.measured_mass = fm.observed_mass;                                    // raw own-scan-frame mass (MS3 subsequence; MS2 == observed)
      obs.theoretical_mass = (ps.ms_level == 3) ? fm.theoretical_mass : 0.0;   // MS3-only carried theoretical (D1=A)
      obs.intensity = intensity;
      obs.source_scan_id = ps.scan_id;
      obs.params = ps.params;
      obs.frag_mz = matched_mz;
      obs.frag_charge = matched_charge;
      obs.iso_width = matched_iso_width;
      obs.stage1_scores = matched_stage1;
      obs.includes_ptm = fm.includes_ptm;   // carry the MS3 localization verdict (bool about THIS fragment's own coverage; frame-invariant)

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
        // Track the best MS2 observation per charge state for config-gated multi-charge MS3.
        auto cit = mfrag.ms2_by_charge.find(obs.frag_charge);
        if (cit == mfrag.ms2_by_charge.end() || obs.intensity > cit->second.intensity)
          mfrag.ms2_by_charge[obs.frag_charge] = obs;
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

      // Localize the mod from bracketing fragments in TWO passes so MS2 is AUTHORITATIVE and the MS3
      // localization verdict can only refine WITHIN the MS2 result (never override it).
      //
      // Pass A -- MS2 mass-test (authoritative). This is the ORIGINAL logic, restricted to best_ms2
      // observations (whose region-frame observed mass matches the region-frame `baseline` here).
      for (const auto& kv : m.fragments)
      {
        const MappedFragment& f = kv.second;
        if (f.ion_type.empty() || f.ion_index <= 0) continue;
        if (!fragmentBrackets(f, rs_region, re_region)) continue;
        if (!f.best_ms2.has_value()) continue;                 // MS2 only in Pass A
        const FragmentObservation* obs = &(*f.best_ms2);

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
        if (matches_with == matches_without) continue;         // uninformative

        const double I = obs->intensity;
        if (matches_with)
        {
          if (f.is_prefix) tightenUpper(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_end, I);
          else             tightenLower(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_start, I);
        }
        else
        {
          if (f.is_prefix) tightenLower(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_end + 1, I);
          else             tightenUpper(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_start - 1, I);
        }
      }

      // Pass B -- MS3 localization verdict (subordinate). MS3-ONLY fragments (best_ms3 present, no
      // best_ms2 -- MS2 already voted for shared keys in Pass A) that STILL bracket the MS2-narrowed
      // range. Uses the propagated includes_ptm (true = PTM inside this fragment's coverage; false =
      // outside) instead of re-deriving from the folded adjusted_mass. Because Pass A already narrowed
      // [rs,re] and we re-check fragmentBrackets against it (and tighten is inward-only), an MS3
      // fragment can only refine WITHIN the MS2 result -- it can never move an MS2-set boundary.
      for (const auto& kv : m.fragments)
      {
        const MappedFragment& f = kv.second;
        if (f.ion_type.empty() || f.ion_index <= 0) continue;
        if (!f.best_ms3.has_value() || f.best_ms2.has_value()) continue;   // MS3-only
        if (!fragmentBrackets(f, rs_region, re_region)) continue;          // re-check vs the narrowed range
        const FragmentObservation& o = *f.best_ms3;
        const double I = o.intensity;

        if (o.includes_ptm)      // PTM IS inside this fragment's coverage (matched the with-variant)
        {
          if (f.is_prefix) tightenUpper(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_end, I);
          else             tightenLower(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_start, I);
        }
        else                     // PTM is OUTSIDE this fragment's coverage (matched the without-variant)
        {
          if (f.is_prefix) tightenLower(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_end + 1, I);
          else             tightenUpper(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_start - 1, I);
        }
      }

      // Convert the narrowed boundaries back to FULL-PROTEIN before storing.
      mod.candidate_start = rs_region + ws;
      mod.candidate_end = re_region + ws;
    }
  }

  void ProteoformTracker::alignedCombinedLists_(const ProteoformModel& m,
      std::vector<double>& measured, std::vector<double>& adjusted, std::vector<double>& theoretical,
      std::vector<double>& diff_da, std::vector<double>& diff_ppm) const
  {
    std::vector<FragmentKey> keys;
    keys.reserve(m.fragments.size());
    for (const auto& kv : m.fragments) keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());  // FragmentKey = pair<string,int> => (ion_type, ion_index) order
    auto push = [&](const FragmentObservation& o) {
      const double adj = o.observed_mass;      // observed_mass holds the ADJUSTED (MS2-frame) value
      const double th  = o.theoretical_mass;
      measured.push_back(o.measured_mass);
      adjusted.push_back(adj);
      theoretical.push_back(th);
      diff_da.push_back(th > 0.0 ? (adj - th) : 0.0);            // guard: MS2 (th==0) => 0, not adj-0
      diff_ppm.push_back(th > 0.0 ? ((adj - th) / th * 1e6) : 0.0);
    };
    for (const FragmentKey& k : keys)
    {
      const MappedFragment& f = m.fragments.at(k);
      if (f.best_ms2.has_value()) push(*f.best_ms2);
      if (f.best_ms3.has_value()) push(*f.best_ms3);
    }
  }

  void ProteoformTracker::emitRow_(const ProteoformModel& m, const std::string& trigger, const std::string& trigger_scan_id)
  {
    // Guard: no identified model -> nothing to emit.
    if (m.proteoform_sequence.empty()) return;

    IdaLogger::PooledModelDescriptor r;

    // --- precursor_id: the model key (per-MS1-selection identity, logged as a plain decimal) ---
    r.precursor_id = m.precursor_id;

    // --- mono_mass: the MS2 context precursor mass (captured once at feedScan time) ---
    r.mono_mass = m.has_ms2_ctx ? m.ms2_ctx.mono_mass : 0.0;

    // --- nominal_mass: derived from the captured precursor mass (the model is now keyed by
    //     precursor_id, so nominal_mass is no longer the key — recompute the real nominal mass
    //     so the existing column keeps its meaning). getNominalMass is a read-only static. ---
    r.nominal_mass = m.has_ms2_ctx ? SpectralDeconvolution::getNominalMass(m.ms2_ctx.mono_mass) : 0;

    // --- proforma ---
    // m.proteoform_sequence is the FULL protein sequence (FragmentAnalysis sets it from protein_sequence).
    // m.modifications are FULL-PROTEIN 1-based ranges.
    // m.region_start/region_end are 0-based (exclusive end), -1 => full sequence.
    // We render the WINNER-REGION substring with region-relative 1-based PTM positions,
    // consistent with the narrowModifications_ frame convention (region_pos = full_pos - ws).
    {
      const int ws = (m.region_start < 0) ? 0 : m.region_start;
      const int P = static_cast<int>(m.proteoform_sequence.size());
      const int we = (m.region_end < 0) ? P : m.region_end;
      const int L = we - ws;
      const std::string region_seq =
        (ws == 0 && we == P) ? m.proteoform_sequence : m.proteoform_sequence.substr(ws, L);

      std::vector<FragmentAnalysis::PTMSite> ptm_sites;
      ptm_sites.reserve(m.modifications.size());
      for (const ModificationState& mod : m.modifications)
      {
        // Convert full-protein 1-based to region-relative 1-based (no-op when ws==0).
        const int s = mod.candidate_start - ws;
        const int e = mod.candidate_end - ws;
        FragmentAnalysis::PTMSite site;
        site.start_position = s;
        site.end_position = e;
        site.position = (s + e) / 2;
        site.mass_shift = mod.mass_shift;
        ptm_sites.push_back(site);
      }
      r.proforma = FragmentAnalysis::toProForma(region_seq, ptm_sites);
    }

    // --- score / coverage / fragment count ---
    r.score = m.identification_score;
    r.coverage_pct = m.coveragePct();
    r.n_fragments = static_cast<int>(m.fragments.size());

    // --- localized_mods and ambiguous_mods ---
    // Format localized (start==end): "<start>[+<mass>]"
    // Format ambiguous (start<end):  "(<start>-<end>)[+<mass>]"
    for (const ModificationState& mod : m.modifications)
    {
      std::ostringstream ss;
      ss << std::fixed << std::setprecision(4);
      if (mod.candidate_start == mod.candidate_end)
      {
        ss << mod.candidate_start;
        if (mod.mass_shift >= 0)
          ss << "[+" << mod.mass_shift << "]";
        else
          ss << "[" << mod.mass_shift << "]";
        r.localized_mods.push_back(ss.str());
      }
      else
      {
        ss << "(" << mod.candidate_start << "-" << mod.candidate_end << ")";
        if (mod.mass_shift >= 0)
          ss << "[+" << mod.mass_shift << "]";
        else
          ss << "[" << mod.mass_shift << "]";
        r.ambiguous_mods.push_back(ss.str());
      }
    }

    // contributing_scan_ids: the model's CUMULATIVE fed-scan set (monotone; a superseded scan never
    // drops). Was rebuilt here from current best_ms2/best_ms3 sources, which dropped superseded scans.
    r.contributing_scan_ids.assign(m.contributing_scan_ids.begin(), m.contributing_scan_ids.end());

    // --- aligned combined lists: five parallel per-fragment vectors in stable FragmentKey order ---
    std::vector<double> _measured, _adjusted, _theoretical, _diff_da, _diff_ppm;
    alignedCombinedLists_(m, _measured, _adjusted, _theoretical, _diff_da, _diff_ppm);
    r.combined_measured    = std::move(_measured);
    r.combined_masses      = std::move(_adjusted);      // combined_ms2_frame_masses column, now FragmentKey order (not mass-sorted)
    r.combined_theoretical = std::move(_theoretical);
    r.combined_diff_da     = std::move(_diff_da);
    r.combined_diff_ppm    = std::move(_diff_ppm);

    // --- update_index ---
    r.update_index = m.update_index;

    r.trigger = trigger;
    r.trigger_scan_id = trigger_scan_id;

    logger_.writePooledModelRow(r);
  }

} // namespace OpenMS
