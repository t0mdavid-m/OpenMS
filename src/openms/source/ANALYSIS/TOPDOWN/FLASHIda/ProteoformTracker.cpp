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
#include <iostream>
#include <limits>
#include <map>
#include <cstdlib>   // std::abs(int) -- the exhaustive pool's charge floor compares absolute charges
#include <set>
#include <sstream>
#include <string>
#include <vector>
// Identification entry points already live here as statics -- identifyExplorationFragments,
// selectNextLevelTargets, scoreCalibratedVariants -- and Exploration calls into them rather than
// driving FLASHTagger/FLASHExtender itself (Exploration::computeFragmentMatch_ is a marshalling
// wrapper around the first). So the responsibility split is largely in place.
//
// What is NOT in place, recorded here so it is not rediscovered from scratch: an MS3 FragmentCount
// exploration group of N variants identifies each variant's spectrum TWICE. Once per variant on the
// feed path (computeExplorationScore_ -> computeFragmentMatch_, every metric branch calls it), then all
// N again in scoreCalibratedVariants once all_received flips -- which re-matches from scratch by
// design, two-pass calibration plus a tight-tolerance rematch, and deliberately reuses nothing from
// the feed path. MS2 groups pay this once; only MS3 FragmentCount doubles.
//
// Whether the two passes are interchangeable is the open question: they run at different tolerances
// (exploration_tolerance_ppm vs LOOSE_TOLERANCE_PPM then level tolerance) and, since ADR-0002, against
// different contexts (the triggering scan's render context vs the tracker's live winner). So this is a
// caching problem with a correctness precondition, not a code-motion problem, and it is not addressed.
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
    // The label an exhaustive-mode target carries when no mapped fragment claims its mass (ADR-0023
    // decision 5). An IN-ENGINE SENTINEL ONLY: paired with ion_index 0 it takes buildMS3's no-ion
    // branch, so it never reaches the wire and logs as an empty ion_type (D-f). Its whole job is to
    // fail MS3FragmentMatcher::isKnownIonClass, so no projection site can cut a subsequence for it.
    constexpr char kUnassignedIonType[] = "u";

    // The authored spelling of an objective, for the [MS3-PLAN] marker. A switch with no `default:`
    // so the next objective is a compiler warning rather than a silently wrong label.
    inline const char* objectiveName(CharacterizationObjective o)
    {
      switch (o)
      {
        case CharacterizationObjective::Ambiguity: return "ambiguity";
        case CharacterizationObjective::Coverage: return "coverage";
        case CharacterizationObjective::Exhaustive: return "exhaustive";
      }
      return "ambiguity";
    }

    // The exhaustive pool's admission test, factored out because it answers "is this peak a target?"
    // and that question gets asked in more than one place as soon as the escalation ladder (ADR-0022)
    // lands: its stopping condition is "would planExhaustive_ return anything?", which is this test
    // plus the dispatch memory. One definition, so the two can never drift.
    //
    // The dispatch-memory test is deliberately NOT part of this: that is per-mass state each caller
    // checks against the same set, not a property of the peak.
    inline bool exhaustivePoolAdmits(const PeakRecord& pr, double min_target_mass, int charge_floor)
    {
      if (min_target_mass > 0.0 && pr.mono_mass < min_target_mass) return false;
      if (charge_floor > 0 && std::abs(pr.charge) < charge_floor) return false;
      // Unisolatable: no isolation centre, or no charge to isolate it at. Exploration.cpp:950/:957
      // drops both anyway -- but only AFTER planNextScans has returned, so admitting one here would
      // stamp its nominal mass into a monotone dispatch memory and burn that mass for the rest of the
      // Precursor's life without ever acquiring it (ADR-0023 D-d).
      if (pr.mz <= 0.0 || pr.charge == 0) return false;
      return true;
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

  // --- Identification entry points (#46): thin static forwarders; results are byte-identical to the
  //     former direct matcher calls. They take the caller's FragmentAnalysis& because they also run in
  //     const helpers (computeFragmentMatch_) and where no tracker instance exists (initiateNextLevel).

  int ProteoformTracker::identifyExplorationFragments(FragmentAnalysis& frag, const String& protein_sequence, int n,
                                                      double* masses, double* qscores, int* charges,
                                                      double* window_starts, double* window_ends,
                                                      char* ion_types, int* fragment_indices,
                                                      DeconvolvedSpectrum& stored_ms2,
                                                      FragmentAnalysis::ProteoformMatch& result,
                                                      const String& fragmentation_method,
                                                      double tolerance_ppm,
                                                      FragmentAnalysis::FragmentScores* frag_scores)
  {
    return frag.getTopFragmentMatches(protein_sequence, n, masses, qscores, charges,
                                      window_starts, window_ends, ion_types, fragment_indices,
                                      stored_ms2, result, fragmentation_method, tolerance_ppm, frag_scores);
  }

  int ProteoformTracker::selectNextLevelTargets(FragmentAnalysis& frag, SelectionMetric selection,
                                                const String& protein_sequence, int n,
                                                double* masses, double* qscores, int* charges,
                                                double* window_starts, double* window_ends,
                                                char* ion_types, int* fragment_indices,
                                                DeconvolvedSpectrum& stored_ms2,
                                                FragmentAnalysis::ProteoformMatch& result,
                                                const String& fragmentation_method,
                                                double tolerance_ppm,
                                                FragmentAnalysis::FragmentScores* frag_scores)
  {
    switch (selection)
    {
      case SelectionMetric::Intensity:
      case SelectionMetric::QScore:
        return frag.getTopFragmentMatches(protein_sequence, n, masses, qscores, charges,
                                          window_starts, window_ends, ion_types, fragment_indices,
                                          stored_ms2, result, fragmentation_method, tolerance_ppm, frag_scores);
      case SelectionMetric::TerminalFragments:
        return frag.getTerminalFragmentIons(protein_sequence, n, masses, qscores, charges,
                                            window_starts, window_ends, ion_types, fragment_indices,
                                            stored_ms2, result, fragmentation_method, tolerance_ppm, frag_scores);
      case SelectionMetric::AmbiguityResolution:
        return frag.getAmbiguityEnclosingIons(protein_sequence, n, masses, qscores, charges,
                                              window_starts, window_ends, ion_types, fragment_indices,
                                              stored_ms2, result, fragmentation_method, tolerance_ppm, frag_scores);
      default:
        return 0;
    }
  }

  std::vector<double> ProteoformTracker::scoreCalibratedVariants(
      const std::vector<const DeconvolvedSpectrum*>& variant_spectra, const std::string& protein_sequence,
      const MS3FragmentMatcher::ProteoformContext& ctx, char fragment_ion_type, int fragment_ion_index,
      double loose_tolerance_ppm, double tight_tolerance_ppm,
      std::vector<FragmentAnalysis::ProteoformMatch>* detailed_results)
  {
    return MS3FragmentMatcher::calibrateAndScore(variant_spectra, protein_sequence, ctx, fragment_ion_type,
                                                 fragment_ion_index, loose_tolerance_ppm, tight_tolerance_ppm,
                                                 detailed_results);
  }

  void ProteoformTracker::feedScan(int precursor_id, uint8_t ms_level, const Ms2Params& params, int scan_id,
                                   const DeconvolvedSpectrum& deconv,
                                   const FragmentAnalysis::ProteoformMatch& match, double id_score,
                                   const ScanCommand& ms2_ctx)
  {
    ProteoformModel& m = models_[precursor_id];
    m.precursor_id = precursor_id;
    m.contributing_scan_ids.insert(scan_id);

    // Capture the MS2 command context
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
      pr.iso_width = mz2 - mz1;
      pr.stage1_scores = FragmentAnalysis::FragmentScores::fromPeakGroup(pg, charge);
      // The FULL envelope alongside the representative charge, so fragment_charges has real
      // candidates. Same arithmetic as the four lines above, once per present charge. No SNR gate and
      // no ordering here: selectNotches owns both, in one place (ADR-0016).
      auto [z_lo, z_hi] = pg.getAbsChargeRange();
      for (int z = z_lo; z <= z_hi; ++z)
      {
        auto [zmz1, zmz2] = pg.getMzRange(z);
        if (zmz2 <= zmz1) { continue; }   // charge not present in this envelope
        ChargeRecord cr;
        cr.charge = z;
        cr.mz = (zmz1 + zmz2) / 2.0;
        cr.iso_width = zmz2 - zmz1;
        cr.intensity = static_cast<double>(pg.getChargeIntensity(z));
        cr.stage1_scores = FragmentAnalysis::FragmentScores::fromPeakGroup(pg, z);
        pr.by_charge.push_back(cr);
      }
      ps.peaks.push_back(pr);
    }
    m.pending.push_back(std::move(ps));
  }

  void ProteoformTracker::finalizeMS2(int precursor_id)
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

    // Retain the winner scan's RAW deconvolved peaks and the parameters it was acquired under, for
    // characterization.mode == Exhaustive (ADR-0023). Every other objective targets out of
    // `fragments`, a table of theoretical ions, which by construction keeps nothing about the masses
    // that matched no ion -- exactly the masses this mode exists to fragment.
    //
    // Captured HERE, where the winner is chosen, because that is the one point at which the winning
    // scan's peak list is both identified and still in hand. A copy, not a pointer: a pointer into a
    // scan container would tie this mode's correctness to that container's retention policy.
    m.winner_peaks = win->peaks;
    m.winner_params = win->params;

    // 3) Seed modifications from the winner's PTM sites (full-protein 1-based ranges).
    m.modifications.clear();
    for (const FragmentAnalysis::PTMSite& s : win->match.ptm_sites)
    {
      m.modifications.push_back(ModificationState{s.mass_shift, s.start_position, s.end_position, 0.0, 0.0});
    }

    // 4) Map EVERY staged scan's already-matched fragments onto the model
    m.fragments.clear();
    m.rematched_nonwinner_.clear();   // C: rebuilt below by mapNonWinnerMs2_ for empty-own-match contributors
    for (const PendingScan& ps : m.pending)
    {
      mapScanOntoModel_(m, ps);
    }

    // 5) Narrow modifications
    narrowModifications_(m);
    ++m.update_index;
    emitPooledIDRow(m, "MS2", ScanCommandQueue::encode(m.winner_scan_id));

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
    emitPooledIDRow(m, trigger_ion, trigger_scan_id);
    m.pending.clear();
  }

  std::vector<Ms3Target> ProteoformTracker::planNextScans(int precursor_id)
  {
    // "This run fired no MS3" resolves to one of NINE causes below (seven return sites; the last one
    // splits three ways by objective). Every one of them used to return an empty plan silently,
    // leaving the reason to be inferred from empty columns in pooled_identification.tsv. Name it
    // instead (ADR-0020). The zero-target outcome is frequently CORRECT -- ambiguity mode with every
    // modification already localized has nothing an MS3 could narrow -- so this is an explanation,
    // not a warning, and there is deliberately no severity attached to it.
    //
    // stdout by design, and stdout ONLY: the engine has no file stream for decisions, and adding one
    // to ida.log would move all 17 log goldens. Accepted consequence -- an instrument acquisition
    // discards the engine's stdout, so this is visible in offline Flash.exe runs and in CI (which
    // greps regression-stdout.txt) but NOT in an instrument run folder. Every other engine marker
    // ([TRACK-CREATE], [EXPL-WINNER], [EXPL-ABORT]) already has exactly this limitation.
    auto no_plan = [precursor_id](const char* reason) -> std::vector<Ms3Target> {
      std::cout << "[MS3-PLAN] precursor_id=" << precursor_id
                << " targets=0 reason=" << reason << std::endl;
      return {};
    };

    auto it = models_.find(precursor_id);
    if (it == models_.end()) return no_plan("no_model");
    ProteoformModel& m = it->second;
    // No identified model / no captured MS2 context -> no plan (ADR-0002).
    if (m.proteoform_sequence.empty()) return no_plan("unidentified_precursor");
    if (!m.has_ms2_ctx) return no_plan("no_ms2_context");

    const CharacterizationObjective obj = config_.characterization().objective;
    const int budget = config_.level(2).max_targets;  // reuse the MS2 max_targets as the MS3 budget
    if (budget <= 0) return no_plan("zero_budget");

    // The fork sits HERE deliberately: after every identification/context/budget guard, so decision 10
    // ("identification is still required") reports through the SAME reason strings in every mode, and
    // before the objective-keyed pool build below, which exhaustive replaces wholesale (ADR-0023).
    if (obj == CharacterizationObjective::Exhaustive) return planExhaustive_(m, precursor_id, budget);

    const int P = static_cast<int>(m.proteoform_sequence.size());
    if (P <= 0) return no_plan("empty_sequence");
    const int ws = (m.region_start < 0) ? 0 : m.region_start;
    const int we = (m.region_end < 0) ? P : m.region_end;
    const int L = we - ws;  // winner-region residue count (fragment frame)
    if (L <= 0) return no_plan("empty_region");

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

    // Distinguishes the two ambiguity dead-ends in the reason line below: "nothing left to resolve"
    // (every mod localized -- the expected end state of a successful run) versus "something to resolve
    // but no fragment spans it" (a real coverage gap). Both yield zero targets; only the second is a
    // reason to change the method.
    bool had_ambiguous_mod = false;

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
      had_ambiguous_mod = !ambiguous.empty();
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
      // Backbone-cleavage-site coverage (same model as ProteoformModel::coveragePct, NOT a residue-span
      // union). Each observed fragment witnesses exactly ONE bond in MS2 (prefix b_k -> bond cover_end;
      // suffix -> bond cover_start-1). "Uncovered" = the L-1 inter-residue bonds not yet witnessed.
      // MS3-ing a parent re-fragments its span, so it can witness the bonds INTERIOR to that span
      // (prefix [1..k] -> bonds 1..k-1; suffix [c..L] -> bonds c..L-1). Target every observed fragment,
      // strongest best-MS2 first, whose interior adds >=1 still-uncovered bond, marking that interior
      // covered (marginal-skip so we don't re-target the same gap), bounded by budget.
      //
      // NB: a residue-SPAN union (the previous implementation) is unusable here: every fragment's span
      // is a terminal slab [1,k] or [c,L], so the complement is always an INTERNAL gap that no terminal
      // fragment can fully contain (fragmentContains) -> it selected 0 targets for all real data.
      std::set<int> uncovered;
      for (int b = 1; b <= L - 1; ++b) uncovered.insert(b);
      for (const auto& kv : m.fragments)
      {
        const MappedFragment& f = kv.second;
        const int site = f.is_prefix ? f.cover_end : (f.cover_start - 1);
        if (site >= 1 && site <= L - 1) uncovered.erase(site);  // bond already witnessed in MS2
      }

      std::vector<const MappedFragment*> cands;
      for (const auto& kv : m.fragments)
        if (kv.second.best_ms2.has_value()) cands.push_back(&kv.second);
      std::sort(cands.begin(), cands.end(),
                [](const MappedFragment* a, const MappedFragment* b) {
                  return a->best_ms2->intensity > b->best_ms2->intensity;  // strongest first
                });

      for (const MappedFragment* f : cands)
      {
        // Interior bonds an MS3 re-feed of this parent could witness.
        const int lo = f->is_prefix ? 1 : f->cover_start;
        const int hi = f->is_prefix ? (f->cover_end - 1) : (L - 1);
        if (hi < lo) continue;  // no interior bonds (e.g. b1)
        bool adds = false;
        for (int b = lo; b <= hi && !adds; ++b) adds = (uncovered.count(b) > 0);
        if (!adds) continue;  // marginal-skip: contributes no new coverage
        for (int b = lo; b <= hi; ++b) uncovered.erase(b);  // optimistic: MS3 drills its interior
        if (!add_target(f)) break;  // dedup + budget bound (existing lambda; false => budget full)
      }
    }

    if (targets.empty())
    {
      if (obj != CharacterizationObjective::Ambiguity) return no_plan("no_fragment_adds_coverage");
      // The run-B case: MS2 alone localized every modification, so there is no range left for an MS3
      // to narrow and zero targets is the correct answer, not a failure.
      return no_plan(had_ambiguous_mod ? "no_containing_fragment" : "all_mods_localized");
    }

    // Strongest best-MS2 first within the chosen set (the executor dispatches in this order).
    std::sort(targets.begin(), targets.end(),
              [](const MappedFragment* a, const MappedFragment* b) {
                return a->best_ms2->intensity > b->best_ms2->intensity;
              });

    // --- 2) Emit Ms3Targets, budget-bounded ---------------------------------------------------------
    // Single (default): one target per fragment from best_ms2 — byte-identical to the pre-mode build.
    // Separate: one target per observed charge state of the fragment, strongest first — N scans.
    // Multiplexed: ONE target carrying the other SNR-positive charges as notches — 1 scan.
    //
    // THE BUDGET COUNTS FRAGMENTS IN ALL THREE. Separate and Multiplexed both acquire one fragment's
    // whole envelope and both cost one slot for it; they differ only in scan count. Counting emitted
    // targets instead would make max_targets: 3 spend everything on the first fragment that happened to
    // be seen at three charges, characterising ONE cleavage site — which is the pathology the modes
    // exist to avoid, not to reproduce (ADR-0016).
    const ChargeAcquisitionMode charge_mode = config_.characterization().fragment_charges;
    const double snr_threshold = config_.targeting().snr_threshold;
    std::vector<Ms3Target> out;
    int fragments_selected = 0;
    for (const MappedFragment* f : targets)
    {
      if (fragments_selected >= budget) break;
      const size_t emitted_before_fragment = out.size();

      // The acquisition charge set, chosen ONCE for both on-modes: anchor (best_ms2, i.e. the most
      // intense charge) plus the SNR-positive rest, most intense first, capped at a 10-plex.
      // Separate and Multiplexed then acquire exactly this set and differ only in how many scans it
      // takes -- N versus 1. Deriving both from one selectNotches call is what makes that identity
      // structural instead of a comment: they cannot drift into gating or ranking differently.
      std::vector<NotchCandidate> extra_charges;
      if (charge_mode != ChargeAcquisitionMode::Single && !f->ms2_by_charge.empty())
      {
        std::vector<NotchCandidate> cands;
        cands.reserve(f->ms2_by_charge.size());
        for (const auto& kv : f->ms2_by_charge)
          cands.push_back({kv.second.frag_charge, kv.second.frag_mz, kv.second.iso_width,
                           kv.second.stage1_scores.charge_snr, kv.second.intensity});
        extra_charges = selectNotches(cands, f->best_ms2->frag_charge, snr_threshold,
                                      MAX_NOTCHES_PER_STAGE,
                                      "MS3 " + f->ion_type + std::to_string(f->ion_index));
      }

      // Build the per-charge observation list to emit for this fragment.
      std::vector<const FragmentObservation*> obs_list;
      obs_list.push_back(&(*f->best_ms2));   // the anchor, in every mode
      if (charge_mode == ChargeAcquisitionMode::Separate)
      {
        for (const NotchCandidate& n : extra_charges)
        {
          auto zit = f->ms2_by_charge.find(n.charge);
          if (zit != f->ms2_by_charge.end()) obs_list.push_back(&zit->second);
        }
      }

      // Multiplexed: the same set arrives as notches on ONE target, so the fragment stage isolates the
      // whole envelope in a single scan. The fragment stage owns its own MAX_NOTCHES_PER_STAGE block,
      // so this is not shared with the stage-0 precursor set the MS3 inherits from its parent.
      const std::vector<NotchCandidate> frag_notches =
          (charge_mode == ChargeAcquisitionMode::Multiplexed) ? extra_charges
                                                             : std::vector<NotchCandidate>{};

      for (const FragmentObservation* o : obs_list)
      {
        // No budget check on the additional charges of THIS fragment under "separate": the slot was
        // paid for by the per-fragment guard above, and the mode's contract is the whole envelope. In
        // the other two modes obs_list holds exactly one entry, so this is unchanged for them.
        if (charge_mode != ChargeAcquisitionMode::Separate
            && static_cast<int>(out.size()) >= budget) break;
        Ms3Target t;
        t.ion_type = f->ion_type;
        t.ion_index = f->ion_index;
        t.frag_mz = o->frag_mz;
        t.frag_charge = o->frag_charge;
        t.frag_mass = o->observed_mass;
        t.iso_width = o->iso_width;
        t.stage0_params = o->params;
        t.stage1_scores = o->stage1_scores;
        t.notches = frag_notches;  // empty unless fragment_charges == Multiplexed
        out.push_back(std::move(t));
      }
      // One fragment consumed, however many targets it produced.
      if (out.size() > emitted_before_fragment) ++fragments_selected;
    }
    // Counterpart to the no_plan lines: the same marker on the success path, so grepping [MS3-PLAN]
    // yields every MS3 planning decision rather than only the negative ones. fragments and targets
    // differ under separate/multiplexed charge modes -- fragments is what the budget counts (ADR-0016).
    std::cout << "[MS3-PLAN] precursor_id=" << precursor_id
              << " targets=" << out.size()
              << " fragments=" << fragments_selected
              << " budget=" << budget
              << " objective=" << (obj == CharacterizationObjective::Ambiguity ? "ambiguity" : "coverage")
              << std::endl;
    return out;
  }

  std::vector<Ms3Target> ProteoformTracker::planExhaustive_(ProteoformModel& m, int precursor_id, int budget)
  {
    // Same marker and the same zero-target vocabulary as planNextScans' guards, so one grep still
    // yields every MS3 planning decision in every mode.
    auto no_plan = [precursor_id](const char* reason) -> std::vector<Ms3Target> {
      std::cout << "[MS3-PLAN] precursor_id=" << precursor_id
                << " objective=exhaustive targets=0 reason=" << reason << std::endl;
      return {};
    };

    // The winner scan's own peak list, captured at finalize (ADR-0023). Deliberately NOT a lookup into
    // a scan archive: the pool is a property of the model, so the feature does not depend on any other
    // feature's retention policy for its raw material to still exist at plan time.
    if (m.winner_peaks.empty()) return no_plan("no_winner_scan");

    // NOT EMPTY *AND* NOT CAPABLE => refuse. An EMPTY activation counts as CAPABLE, and the asymmetry
    // is deliberate: this mirrors upsertMappedObservation_'s only other capability test, whose comment
    // records the reason at length. "" is not ETD -- it is "no activation recorded", which is what
    // every hand-built C++ fixture and every scan config that omits the key produces. Failing closed on
    // "" would return zero MS3 targets for all of them: a quieter failure than the one being prevented.
    // (ADR-0023 decision 3 reads fail-closed; D-b corrects it.)
    //
    // Spelled out rather than delegated on purpose: the named predicate for this test arrives with the
    // escalation ladder (ADR-0022), which is not in this branch. Collapse these two comparisons into
    // that helper in the same commit that introduces it -- do not add a second definition here, or the
    // two will disagree the first time the supported set changes.
    const std::string& win_act = m.winner_params.activation_type;
    if (!win_act.empty() && win_act != "HCD" && win_act != "CID")
      return no_plan("winner_scan_not_ms3_capable");

    const double min_target_mass = config_.characterization().min_target_mass;
    const int charge_floor = config_.level(2).min_charge;
    const double tol_ppm = config_.level(2).tolerance_ppm;

    // --- 1) The pool: this ONE scan's deconvolved masses, mapped or not ----------------------------
    struct PoolEntry
    {
      const PeakRecord* peak = nullptr;
      int nominal = 0;
    };
    std::vector<PoolEntry> pool;
    pool.reserve(m.winner_peaks.size());
    for (const PeakRecord& pr : m.winner_peaks)
    {
      // ORDER IS LOAD-BEARING: every filter runs HERE, and the dispatch memory is stamped only in the
      // emit loop below. See exhaustivePoolAdmits for what a premature stamp costs (ADR-0023 D-d).
      if (!exhaustivePoolAdmits(pr, min_target_mass, charge_floor)) continue;
      const int nominal = SpectralDeconvolution::getNominalMass(pr.mono_mass);
      if (m.dispatched_nominal_masses.count(nominal) > 0) continue;   // already fragmented this species
      pool.push_back(PoolEntry{&pr, nominal});
    }
    if (pool.empty()) return no_plan("pool_exhausted");

    // Decision 4: intensity, descending, no tiebreak -- the same rule as every other target-ranking
    // site, for ADR-0003's reason (more fragment ion means more MS3 precursor). Mapped and unassigned
    // masses rank on the identical number, so the two halves interleave rather than segregating.
    //
    // stable_sort, not sort: PeakGroup::getChargeIntensity() is 0 for any group that never went through
    // a scoring pass, so an all-ties pool is a real state (every hand-built fixture is one), and an
    // unstable sort would make WHICH masses survive a budget truncation implementation-defined.
    std::stable_sort(pool.begin(), pool.end(),
                     [](const PoolEntry& a, const PoolEntry& b) {
                       return a.peak->intensity > b.peak->intensity;
                     });

    // --- 2) The mapped/unassigned label ------------------------------------------------------------
    // IN-TOLERANCE ONLY, and deliberately unlike mapScanOntoModel_, which falls back to the closest
    // peak OVERALL so that "a matched fragment is never dropped for lack of a peak". That is the right
    // rule for its question (does this KNOWN FRAGMENT have intensity here?) and the wrong rule for this
    // one (does this PEAK deserve a known fragment's name?): replaying the fallback would stamp a real
    // b61 onto an arbitrarily distant peak -- a confident wrong label, which is exactly what decision
    // 5's class guard exists to make impossible. ADR-0023 D-c.
    //
    // Bound against the fragment's own observed MS2-frame mass, i.e. the binding mapScanOntoModel_
    // made, read backwards. NOT against MappedFragment::theoretical_mass: that member is declared and
    // never assigned (upsertMappedObservation_ says so in as many words), so keying on it would label
    // every pool mass 'u' and the mode would never match anything.
    auto labelFor = [&](const PeakRecord& pr) -> std::pair<std::string, int> {
      const std::pair<std::string, int> unassigned{std::string(kUnassignedIonType), 0};
      const double tol_abs = pr.mono_mass * tol_ppm * 1e-6;
      const MappedFragment* best = nullptr;
      double best_diff = 0.0;
      for (const auto& kv : m.fragments)
      {
        const MappedFragment& f = kv.second;
        double mapped_mass = 0.0;
        if (f.best_ms2.has_value()) mapped_mass = f.best_ms2->observed_mass;
        else if (f.best_ms3.has_value()) mapped_mass = f.best_ms3->observed_mass;   // also MS2-frame
        else continue;
        if (mapped_mass <= 0.0) continue;
        const double diff = std::abs(mapped_mass - pr.mono_mass);
        if (diff > tol_abs) continue;
        // Closest wins, FragmentKey breaks an exact tie. `fragments` is an unordered_map, so a
        // first-match-wins rule would make the label depend on hash order -- reproducible within a
        // build and not across them.
        if (best == nullptr || diff < best_diff
            || (diff == best_diff && FragmentKey{f.ion_type, f.ion_index} < FragmentKey{best->ion_type, best->ion_index}))
        {
          best = &f;
          best_diff = diff;
        }
      }
      if (best == nullptr) return unassigned;
      // A label is only usable if the matcher can project through it. mapScanOntoModel_ already admits
      // a/b/c/x/y/z only, so this cannot fire today; it keeps "every target is either a KNOWN class
      // with a positive index, or the 'u'/0 sentinel" a local invariant of this function rather than
      // one inherited from a caller three files away.
      if (best->ion_type.empty() || best->ion_index <= 0
          || !MS3FragmentMatcher::isKnownIonClass(best->ion_type[0]))
        return unassigned;
      return {best->ion_type, best->ion_index};
    };

    // --- 3) Emit, budget-bounded -------------------------------------------------------------------
    // THE BUDGET COUNTS SPECIES, not emitted targets, exactly as it counts fragments on the mapped
    // path: separate and multiplexed acquire one species' whole charge envelope for ONE slot and differ
    // only in scan count (ADR-0016). Counting acquisitions instead is the pathology the modes exist to
    // avoid -- max_targets: 3 spent entirely on the first mass that happened to resolve at 3 charges.
    const ChargeAcquisitionMode charge_mode = config_.characterization().fragment_charges;
    const double snr_threshold = config_.targeting().snr_threshold;

    std::vector<Ms3Target> out;
    int species_selected = 0;
    for (const PoolEntry& e : pool)
    {
      if (species_selected >= budget) break;
      const PeakRecord& pr = *e.peak;
      const std::pair<std::string, int> label = labelFor(pr);

      // The acquisition charge set, chosen ONCE for both on-modes from the peak's own envelope: anchor
      // plus the SNR-positive rest, most intense first, capped at a 10-plex. One selectNotches call for
      // both is what keeps "they differ only in scan count" structural rather than aspirational.
      std::vector<NotchCandidate> extra_charges;
      if (charge_mode != ChargeAcquisitionMode::Single && !pr.by_charge.empty())
      {
        std::vector<NotchCandidate> cands;
        cands.reserve(pr.by_charge.size());
        for (const ChargeRecord& cr : pr.by_charge)
          cands.push_back({cr.charge, cr.mz, cr.iso_width, cr.stage1_scores.charge_snr, cr.intensity});
        extra_charges = selectNotches(cands, pr.charge, snr_threshold, MAX_NOTCHES_PER_STAGE,
                                      "MS3 exhaustive m=" + std::to_string(pr.mono_mass));
      }
      // Multiplexed puts that set on ONE command as notches; separate spends one command per member.
      const std::vector<NotchCandidate> frag_notches =
          (charge_mode == ChargeAcquisitionMode::Multiplexed) ? extra_charges : std::vector<NotchCandidate>{};

      // NOT named `emit`: Qt's qobjectdefs.h defines `emit` as an EMPTY macro, and OpenMS links
      // Qt6Core, so `auto emit = ...` preprocesses to `auto = ...` -- MSVC then reports
      // "no variable declared before '='" on a line that looks perfectly valid.
      auto emitTarget = [&](double mz, int charge, double iso_width, const FragmentAnalysis::FragmentScores& s) {
        Ms3Target t;
        t.ion_type = label.first;
        t.ion_index = label.second;
        t.frag_mz = mz;
        t.frag_charge = charge;
        t.frag_mass = pr.mono_mass;
        t.iso_width = iso_width;
        // stage[0] REPLAYS the MS2 that produced this mass (ADR-0003), and for an unassigned mass there
        // is no per-ion best observation to source it from -- the winner scan's own params ARE the
        // acquisition that produced it, for mapped and unassigned masses alike.
        t.stage0_params = m.winner_params;
        t.stage1_scores = s;
        t.notches = frag_notches;   // empty unless fragment_charges == Multiplexed
        out.push_back(std::move(t));
      };

      emitTarget(pr.mz, pr.charge, pr.iso_width, pr.stage1_scores);
      if (charge_mode == ChargeAcquisitionMode::Separate)
      {
        for (const NotchCandidate& n : extra_charges)
        {
          for (const ChargeRecord& cr : pr.by_charge)
            if (cr.charge == n.charge) { emitTarget(cr.mz, cr.charge, cr.iso_width, cr.stage1_scores); break; }
        }
      }

      // Stamped only now: after every filter, and after the mass has actually become a target. The set
      // is monotone and dispatched-but-never-returned counts as done, so a stamp made any earlier is a
      // species this Precursor can never revisit (ADR-0023 D-d, d7).
      m.dispatched_nominal_masses.insert(e.nominal);
      ++species_selected;
    }

    std::cout << "[MS3-PLAN] precursor_id=" << precursor_id
              << " objective=" << objectiveName(CharacterizationObjective::Exhaustive)
              << " pool=" << pool.size()
              << " targets=" << out.size()
              << " species=" << species_selected
              << " truncated=" << (static_cast<int>(pool.size()) - species_selected)
              << " budget=" << budget
              << std::endl;
    return out;
  }

  const FragmentAnalysis::ProteoformMatch* ProteoformTracker::getRematchedNonWinnerMatch(int precursor_id, int scan_id) const
  {
    auto it = models_.find(precursor_id);
    if (it == models_.end()) return nullptr;
    const auto& map = it->second.rematched_nonwinner_;
    auto mit = map.find(scan_id);
    if (mit == map.end() || mit->second.fragments.empty()) return nullptr;
    return &mit->second;
  }

  const ProteoformModel* ProteoformTracker::getModel(int precursor_id) const
  {
    auto it = models_.find(precursor_id);
    if (it == models_.end())
    {
      return nullptr;
    }
    return &it->second;
  }

  MS3FragmentMatcher::ProteoformContext ProteoformTracker::buildWinnerProteoformContext(int precursor_id) const
  {
    // Default context = region -1/-1, no PTMs. If there is no finalized, non-empty winner, MS3 then
    // matches nothing ("no winner -> nothing to characterize"). This replaces the triggering-scan
    // context at the MS3 build sites so MS3 is scored against the LIVE winner proteoform.
    MS3FragmentMatcher::ProteoformContext ctx;
    const ProteoformModel* m = getModel(precursor_id);
    if (m == nullptr || !m->finalized || m->proteoform_sequence.empty()) return ctx;

    const int P = static_cast<int>(m->proteoform_sequence.size());
    const int ws = (m->region_start < 0) ? 0 : m->region_start;
    const int we = (m->region_end < 0) ? P : m->region_end;
    ctx.region_start = ws;   // 0-based inclusive (frame expected by MS3FragmentMatcher)
    ctx.region_end = we;     // 0-based exclusive

    // EVERY winner mod as a region-1-based PTMSite (same conversion as narrowModifications_ :782-789).
    // Include LOCALIZED mods (candidate_start == candidate_end): a localized point PTM's mass must
    // enter the MS3 theoretical (e.g. localized heme/acetyl) or fragments spanning it stop matching.
    // Ambiguous mods (start < end) drive the matcher's dual with/without.
    for (const ModificationState& mod : m->modifications)
    {
      FragmentAnalysis::PTMSite site;
      site.start_position = mod.candidate_start - ws;   // region 1-based
      site.end_position = mod.candidate_end - ws;
      site.position = (site.start_position + site.end_position) / 2;
      site.mass_shift = mod.mass_shift;
      ctx.ptm_sites.push_back(site);
    }
    return ctx;
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

    // A NON-winner MS2 scan matched against ITS OWN proteoform hypothesis, so trusting its already-
    // matched fragments would pool foreign theoreticals. Instead re-match its RAW deconvolved masses
    // against the WINNER ladder so every pooled MS2 fragment is winner-consistent by construction.
    // The winner MS2 scan (FLASHTnT already did the heavy lifting) and all MS3 scans (winner-matched
    // upstream via buildWinnerProteoformContext) keep the existing already-matched pooling below.
    if (ps.ms_level == 2 && ps.scan_id != m.winner_scan_id)
    {
      mapNonWinnerMs2_(m, ps, ws, L);
      return;
    }

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
      const std::vector<ChargeRecord>* matched_env = nullptr;  ///< the matched group's charge envelope
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
            matched_env = &pr.by_charge;
          }
          else if (in_tol == within_tol && diff < best_diff)
          {
            best_diff = diff;
            intensity = pr.intensity;
            matched_mz = pr.mz;
            matched_charge = pr.charge;
            matched_iso_width = pr.iso_width;
            matched_stage1 = pr.stage1_scores;
            matched_env = &pr.by_charge;
          }
        }
      }

      // 4) Build the observation and upsert it (shared helper: coverage + per-level max-intensity best).
      FragmentObservation obs;
      obs.ms_level = ps.ms_level;
      obs.observed_mass = obs_mass;
      obs.measured_mass = fm.observed_mass;                                    // raw own-scan-frame mass (MS3 subsequence; MS2 == observed)
      obs.theoretical_mass = fm.theoretical_mass;                              // proteoform theoretical: MS3 equiv (calibrateAndScore) or MS2 matcher best_theo
      obs.intensity = intensity;
      obs.source_scan_id = ps.scan_id;
      obs.params = ps.params;
      obs.frag_mz = matched_mz;
      obs.frag_charge = matched_charge;
      obs.iso_width = matched_iso_width;
      obs.stage1_scores = matched_stage1;
      // Deposit the HONEST sub-frame verdict + the complement-flip flag. The mod-aware complement verdict
      // (flip iff a complement map, with the N-terminal-net-loss exception) is applied per-(mod,fragment) in
      // narrowModifications_ Pass B, which alone can see the mod's sign + N-terminal anchor. (A pre-flipped
      // uniform verdict here regressed the N-terminal -89 from M1 to (2-8) while fixing the additive heme.)
      obs.includes_ptm = fm.includes_ptm;
      obs.is_complement_flip = (ps.ms_level == 3) && fm.is_complement_flip;

      upsertMappedObservation_(m, type, is_prefix, winner_idx, L, obs, matched_env);
    }
  }

  void ProteoformTracker::upsertMappedObservation_(ProteoformModel& m, const std::string& ion_type,
      bool is_prefix, int winner_idx, int L, const FragmentObservation& obs,
      const std::vector<ChargeRecord>* envelope)
  {
    const FragmentKey key{ion_type, winner_idx};
    MappedFragment& mfrag = m.fragments[key];
    mfrag.ion_type = ion_type;
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
    // theoretical_mass is left 0 at the key level; alignedCombinedLists_ reports obs.theoretical_mass.

    if (obs.ms_level == 3)
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
      //
      // From the matched PeakGroup's WHOLE envelope when we have it, so fragment_charges sees every
      // charge the fragment was resolved at rather than only the representative one. best_ms2 above
      // is untouched by this: the representative charge is getMaxIntensityAbsCharge(), i.e. already
      // the envelope's most intense, so "single" stays byte-identical.
      if (envelope != nullptr && !envelope->empty())
      {
        for (const ChargeRecord& cr : *envelope)
        {
          FragmentObservation zo = obs;
          zo.frag_charge = cr.charge;
          zo.frag_mz = cr.mz;
          zo.iso_width = cr.iso_width;
          zo.intensity = cr.intensity;
          zo.stage1_scores = cr.stage1_scores;
          auto cit = mfrag.ms2_by_charge.find(cr.charge);
          if (cit == mfrag.ms2_by_charge.end() || zo.intensity > cit->second.intensity)
            mfrag.ms2_by_charge[cr.charge] = zo;
        }
      }
      else   // MS3-sourced observations and any match with no envelope: exactly the former behaviour
      {
        auto cit = mfrag.ms2_by_charge.find(obs.frag_charge);
        if (cit == mfrag.ms2_by_charge.end() || obs.intensity > cit->second.intensity)
          mfrag.ms2_by_charge[obs.frag_charge] = obs;
      }
      ++mfrag.n_ms2;   // one per contributing scan, not per charge
    }
  }

  void ProteoformTracker::mapNonWinnerMs2_(ProteoformModel& m, const PendingScan& ps, int ws, int L)
  {
    // Build the WINNER theoretical ladder ONCE (winner-region frame) from the CURRENT (pre-narrow)
    // winner mods -- the same region_seq + computePTMAdjustedFragmentMasses that narrowModifications_
    // uses (:751-792). Ambiguous mods enter at their range boundary, so a BRACKETING ion carries the
    // "mod OUT" value in `base`; "mod IN" = base + shift, tested per bracketed mod below.
    const String region_seq =
      (m.region_start < 0) ? String(m.proteoform_sequence) : String(m.proteoform_sequence).substr(ws, L);

    // Ion types produced by THIS scan's fragmentation. A CE-sweep variant shares the winner's
    // activation string (only collision energy differs); getIonTypesForFragmentationMethod lower-cases
    // it. We match this scan's peaks against the winner ladder for exactly the ion types it produces.
    const std::vector<std::string> ion_types =
      FragmentAnalysis::getIonTypesForFragmentationMethod(ps.params.activation_type);
    if (ion_types.empty()) return;

    // ALL winner mods (localized AND ambiguous) as region-1-based PTMSites (same conversion as
    // narrowModifications_ :782-789). Localized mods bake into the ladder; ambiguous ones are bracketed.
    std::vector<FragmentAnalysis::PTMSite> winner_ptms;
    for (const ModificationState& mod : m.modifications)
    {
      FragmentAnalysis::PTMSite site;
      site.start_position = mod.candidate_start - ws;   // region 1-based
      site.end_position = mod.candidate_end - ws;
      site.position = (site.start_position + site.end_position) / 2;
      site.mass_shift = mod.mass_shift;
      winner_ptms.push_back(site);
    }

    std::map<char, std::vector<double>> ladder;   // region-frame theoreticals (ambiguous mods at range boundary)
    FragmentAnalysis::computePTMAdjustedFragmentMasses(region_seq, winner_ptms, ion_types, ladder);

    const double tol_ppm = config_.level(2).tolerance_ppm;

    // C: if THIS non-winner scan did NOT self-identify (empty own match), seed a per-scan re-matched
    // identification (winner proteoform + full-protein winner mods, like a self-ID row). Fragments this
    // scan's masses match on the winner ladder are appended below; Exploration emits the row afterwards.
    FragmentAnalysis::ProteoformMatch* rematch = nullptr;
    if (ps.match.fragments.empty())
    {
      FragmentAnalysis::ProteoformMatch& rm = m.rematched_nonwinner_[ps.scan_id];
      rm.proteoform_sequence = m.proteoform_sequence;
      rm.region_start = m.region_start;
      rm.region_end = m.region_end;
      rm.ptm_sites.clear();
      for (const ModificationState& mod : m.modifications)
      {
        FragmentAnalysis::PTMSite s;
        s.start_position = mod.candidate_start;   // FULL-PROTEIN 1-based (matches a self-ID row's frame)
        s.end_position = mod.candidate_end;
        s.position = (s.start_position + s.end_position) / 2;
        s.mass_shift = mod.mass_shift;
        rm.ptm_sites.push_back(s);
      }
      rematch = &rm;
    }

    for (const PeakRecord& pr : ps.peaks)
    {
      const double obs_mass = pr.mono_mass;

      // Assign this raw mass to the unique winner ion it lands on (Q7): keep only ions that match
      // EXACTLY ONE candidate theoretical (base, or base + a subset of the shifts of the ambiguous
      // mods it brackets); an ion matching >=2 candidates (base AND base+shift) is dropped as dodgy;
      // across distinct surviving ions the closest ppm wins.
      std::string best_type;
      int best_idx = 0;
      bool best_is_prefix = false;
      double best_theo = 0.0;
      double best_ppm = std::numeric_limits<double>::max();
      bool have_best = false;

      for (const std::string& ts : ion_types)
      {
        if (ts.empty()) continue;
        const char ic = ts[0];
        const bool is_prefix = (ic == 'a' || ic == 'b' || ic == 'c');
        auto lit = ladder.find(ic);
        if (lit == ladder.end()) continue;
        const std::vector<double>& masses = lit->second;

        for (int k = 1; k <= static_cast<int>(masses.size()); ++k)
        {
          const double base = masses[k - 1];   // ambiguous mods this ion brackets are OUT of `base`
          // Winner-frame coverage of ion (ic, k): prefix b_k covers [1,k]; suffix y_k covers [L-k+1,L].
          const int cover_end = is_prefix ? k : L;
          const int cover_start = is_prefix ? 1 : (L - k + 1);

          // Shifts of the AMBIGUOUS winner mods this ion brackets (region frame; localized never bracket).
          std::vector<double> bracket_shifts;
          for (const ModificationState& mod : m.modifications)
          {
            const int s = mod.candidate_start - ws;
            const int e = mod.candidate_end - ws;
            if (s >= e) continue;   // localized -> baked into `base`, never bracketed
            const bool brackets = is_prefix ? (s <= cover_end && cover_end < e)
                                            : (s < cover_start && cover_start <= e);
            if (brackets) bracket_shifts.push_back(mod.mass_shift);
          }

          // Candidate theoreticals = base + any subset of the bracketed shifts (2^n, n tiny). Deduped
          // by exact value so equal-shift subsets (e.g. two identical PTMs) count as ONE mass
          // hypothesis, while a genuine base-vs-base+shift pair (distinct values) still counts as two.
          const int nb = static_cast<int>(bracket_shifts.size());
          std::vector<double> cand;
          cand.reserve(static_cast<size_t>(1) << nb);
          for (int mask = 0; mask < (1 << nb); ++mask)
          {
            double theo = base;
            for (int b = 0; b < nb; ++b)
              if (mask & (1 << b)) theo += bracket_shifts[b];
            if (theo > 0.0) cand.push_back(theo);
          }
          std::sort(cand.begin(), cand.end());

          int n_match = 0;
          double theo_here = 0.0;
          double ppm_here = std::numeric_limits<double>::max();
          for (size_t i = 0; i < cand.size(); ++i)
          {
            if (i > 0 && !(cand[i] > cand[i - 1])) continue;   // dedup exact-equal (sorted: equal iff not greater)
            const double d = std::abs(obs_mass - cand[i]);
            if (d <= cand[i] * tol_ppm * 1e-6)
            {
              ++n_match;
              const double ppm = d / cand[i] * 1e6;
              if (ppm < ppm_here) { ppm_here = ppm; theo_here = cand[i]; }
            }
          }
          if (n_match != 1) continue;   // 0 -> no match; >=2 -> distinct base AND base+shift -> dodgy, drop

          if (!have_best || ppm_here < best_ppm)   // >=2 distinct ions -> closest-ppm wins
          {
            have_best = true;
            best_ppm = ppm_here;
            best_type = ts;
            best_idx = k;
            best_is_prefix = is_prefix;
            best_theo = theo_here;
          }
        }
      }

      if (!have_best) continue;   // no unique winner-ladder match -> drop this mass

      // Deposit as a best_ms2 observation; theoretical is the winner ladder value by construction.
      FragmentObservation obs;
      obs.ms_level = 2;
      obs.observed_mass = obs_mass;
      obs.measured_mass = obs_mass;                // MS2 fragment: measured == observed (protein frame)
      obs.theoretical_mass = best_theo;
      obs.intensity = pr.intensity;
      obs.source_scan_id = ps.scan_id;
      obs.params = ps.params;
      obs.frag_mz = pr.mz;
      obs.frag_charge = pr.charge;
      obs.iso_width = pr.iso_width;
      obs.stage1_scores = pr.stage1_scores;
      obs.includes_ptm = false;                    // MS2: Pass A re-derives base-vs-shift per mod at narrow time
      // Envelope carried here too: this path and the winner path upsert into the SAME MappedFragment,
      // so a fragment first seen via a non-winner re-match must not arrive charge-flattened.
      upsertMappedObservation_(m, best_type, best_is_prefix, best_idx, L, obs, &pr.by_charge);

      // C: also record this re-match as an identification fragment for the per-scan row (winner-region frame
      // ion index; == proteoform frame for a full-region winner, the only case in current data).
      if (rematch != nullptr)
      {
        FragmentAnalysis::ProteoformMatch::FragmentMatch fmr;
        fmr.ion_type = best_type;
        fmr.ion_index = best_idx;
        fmr.observed_mass = obs_mass;
        fmr.theoretical_mass = best_theo;
        fmr.diff_da = obs_mass - best_theo;
        fmr.diff_ppm = (best_theo > 0.0) ? (fmr.diff_da / best_theo * 1e6) : 0.0;
        rematch->fragments.push_back(fmr);
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

        // Complement-aware, MOD-AWARE localization verdict. o.includes_ptm = the matched sub-ion's honest
        // verdict; on a complement flip the equiv ion (deposited under the prefix key) is the COMPLEMENT, so
        // it carries the mod iff the sub does NOT -- EXCEPT an N-terminal net-loss composite (Met-excision +
        // N-alpha-acetyl = -89, mass_shift<0 anchored at residue 1), which stays localized to the N-terminus.
        bool includes = o.includes_ptm;
        if (o.is_complement_flip)
        {
          const bool nterm_loss = (mod.mass_shift < 0.0) && (mod.candidate_start == 1);
          includes = nterm_loss ? o.includes_ptm : !o.includes_ptm;
        }

        if (includes)            // PTM IS inside this equiv ion's coverage
        {
          if (f.is_prefix) tightenUpper(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_end, I);
          else             tightenLower(rs_region, re_region, mod.support_lower, mod.support_upper, f.cover_start, I);
        }
        else                     // PTM is OUTSIDE this equiv ion's coverage
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
      std::vector<std::string>& ions,
      std::vector<double>& measured, std::vector<double>& adjusted, std::vector<double>& theoretical,
      std::vector<double>& diff_da, std::vector<double>& diff_ppm) const
  {
    std::vector<FragmentKey> keys;
    keys.reserve(m.fragments.size());
    for (const auto& kv : m.fragments) keys.push_back(kv.first);
    std::sort(keys.begin(), keys.end());  // FragmentKey = pair<string,int> => (ion_type, ion_index) order
    auto push = [&](const std::string& label, const FragmentObservation& o) {
      const double adj = o.observed_mass;      // observed_mass holds the ADJUSTED (MS2-frame) value
      const double th  = o.theoretical_mass;
      ions.push_back(label);
      measured.push_back(o.measured_mass);
      adjusted.push_back(adj);
      theoretical.push_back(th);
      diff_da.push_back(th > 0.0 ? (adj - th) : 0.0);            // guard: no matched theoretical (th==0) => 0
      diff_ppm.push_back(th > 0.0 ? ((adj - th) / th * 1e6) : 0.0);
    };
    for (const FragmentKey& k : keys)
    {
      const std::string label = k.first + std::to_string(k.second);   // ion_type + ion_index, e.g. "b22"
      const MappedFragment& f = m.fragments.at(k);
      if (f.best_ms2.has_value()) push(label, *f.best_ms2);
      if (f.best_ms3.has_value()) push(label, *f.best_ms3);
    }
  }

  void ProteoformTracker::emitPooledIDRow(const ProteoformModel& m, const std::string& trigger, const std::string& trigger_scan_id)
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
    // Region-relative PTM positions: the proteoform ProForma above renders the winner-REGION substring with
    // positions numbered from 1 within that region (candidate_start - ws). Emit the localized/ambiguous mod
    // COLUMNS in the SAME frame so they agree with the ProForma string. When ws==0 (every current winner) this
    // is byte-identical; it only differs for a sub-region winner with region_start>0.
    const int mod_ws = (m.region_start < 0) ? 0 : m.region_start;
    for (const ModificationState& mod : m.modifications)
    {
      std::ostringstream ss;
      ss << std::fixed << std::setprecision(4);
      if (mod.candidate_start == mod.candidate_end)
      {
        ss << (mod.candidate_start - mod_ws);
        if (mod.mass_shift >= 0)
          ss << "[+" << mod.mass_shift << "]";
        else
          ss << "[" << mod.mass_shift << "]";
        r.localized_mods.push_back(ss.str());
      }
      else
      {
        ss << "(" << (mod.candidate_start - mod_ws) << "-" << (mod.candidate_end - mod_ws) << ")";
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

    // --- aligned combined lists: parallel per-fragment vectors in stable FragmentKey order ---
    std::vector<std::string> _ions;
    std::vector<double> _measured, _adjusted, _theoretical, _diff_da, _diff_ppm;
    alignedCombinedLists_(m, _ions, _measured, _adjusted, _theoretical, _diff_da, _diff_ppm);
    r.combined_masses               = std::move(_adjusted);      // combined_ms2_frame_masses column (FragmentKey order, not mass-sorted)
    r.combined_ms2_fragment_ions    = std::move(_ions);
    r.combined_measured             = std::move(_measured);
    r.combined_theoretical          = std::move(_theoretical);
    r.combined_diff_da              = std::move(_diff_da);
    r.combined_diff_ppm             = std::move(_diff_ppm);

    // --- update_index ---
    r.update_index = m.update_index;

    r.trigger = trigger;
    r.trigger_scan_id = trigger_scan_id;

    logger_.writePooledModelRow(r);
  }

} // namespace OpenMS
