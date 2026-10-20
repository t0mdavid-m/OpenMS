// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/MassLadderAlgorithm.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <boost/math/special_functions/gamma.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
namespace
{
// the residues of terminal and precursor ladders
const std::string standard_residues = "ACDEFGHIKLMNPQRSTVWYUO";

// BEGIN precursor ladder kernel
// Precursor-anchored ladders (option precursor_ladders): the windows (i, j) of the database proteins whose unmodified mass equals the
// precursor mass, scored by the deconvolved masses on their prefix ladder (read from i) and their suffix ladder (read from j). This block
// uses the standard library, boost::math::gamma_p and MassLadderAlgorithm::Segment (three integers) only, so that it can be compiled on
// its own and checked against its prototype (kernel K4 of final_all/sliding/src/slide.cpp, calibrated by entrapment on Eclip_1480).
namespace precursor_ladder
{
// the residue sums a..b of a protein (residues a..b-1, all standard): a window never crosses a non-standard residue
using Segment = MassLadderAlgorithm::Segment;

// The deconvolved masses of a spectrum read as residue sums of one ion type: a mass m matches the residue sum v when v is within tol ppm
// of max(m, t) of t = m + shift. Masses with t < 50 are left out (no residue is lighter). A bitset of 1.000495 Da bins only pre-filters
// the interval test, so it never changes a result.
class TargetSet
{
public:
  TargetSet(const std::vector<double>& masses, double shift, bool prefix, double tol_ppm): prefix_(prefix)
  {
    cum_weight_.push_back(0.0);
    for (std::size_t i = 0; i < masses.size(); i++)
    {
      const double m = masses[i], t = m + shift;
      if (t < 50) continue;
      const double margin = tol_ppm * 1e-6 * std::max(m, t);
      t_.push_back(t);
      lo_.push_back(t - margin);
      hi_.push_back(t + margin);
      source_.push_back(i);
      cum_weight_.push_back(cum_weight_.back() + 2 * margin);
    }
    max_value_ = hi_.empty() ? 0 : hi_.back();
    const std::size_t bit_count = (std::size_t)(max_value_ * bin_inverse_ + bin_phase_) + 2;
    bits_.assign(bit_count / 64 + 2, 0);
    for (std::size_t i = 0; i < lo_.size(); i++)
    {
      const std::size_t b0 = (std::size_t)(lo_[i] * bin_inverse_ + bin_phase_), b1 = (std::size_t)(hi_[i] * bin_inverse_ + bin_phase_);
      for (std::size_t b = b0; b <= b1; b++) bits_[b >> 6] |= std::uint64_t(1) << (b & 63);
    }
  }

  // a prefix ion type (the ladder is read from the window start), else a suffix ion type (read from the window end)
  bool isPrefix() const
  {
    return prefix_;
  }

  // the largest residue sum a mass can match
  double maxValue() const
  {
    return max_value_;
  }

  // the index (into the masses) of a mass that matches the residue sum v (> 0), or -1; the first by interval end if intervals overlap
  long massIndex(double v) const
  {
    if (! (v <= max_value_)) return -1;
    const std::size_t b = (std::size_t)(v * bin_inverse_ + bin_phase_);
    if (((bits_[b >> 6] >> (b & 63)) & 1) == 0) return -1;
    const std::size_t idx = (std::size_t)(std::lower_bound(hi_.begin(), hi_.end(), v) - hi_.begin());
    return idx < hi_.size() && lo_[idx] <= v ? (long)source_[idx] : -1;
  }

  // the summed tolerance widths of the masses with t <= v: the expected number of chance matches per ladder position, times its mass
  // density, up to v
  double weightUpTo(double v) const
  {
    return cum_weight_[(std::size_t)(std::upper_bound(t_.begin(), t_.end(), v) - t_.begin())];
  }

private:
  static constexpr double bin_inverse_ = 1.0 / 1.000495;
  static constexpr double bin_phase_ = 0.5;
  bool prefix_;
  std::vector<double> t_, lo_, hi_, cum_weight_; // ascending; cum_weight_[i]: the widths of the first i
  std::vector<std::size_t> source_;              // the mass index of each target
  double max_value_ = 0;
  std::vector<std::uint64_t> bits_;
};

// a window: residues i..j-1 (0-based) of a protein
struct Window
{
  double log_p = 0;           // -log10 P(X >= k), X ~ Poisson(lambda)
  double lambda = 0;          // expected chance matches
  int k = 0;                  // matched (position, side) pairs: a ladder position counts once per side, whatever ion types match it
  std::uint32_t protein = 0, i = 0, j = 0;
  int isotope = 0;            // the isotope error (-1, 0, 1) at which the window mass equals the precursor mass
  std::vector<int> counts;    // matched positions per target set (for the tags)
};

// LogP descending, k descending, then protein, i and j ascending
inline bool isBetter(const Window& a, const Window& b)
{
  if (a.log_p != b.log_p) return a.log_p > b.log_p;
  if (a.k != b.k) return a.k > b.k;
  if (a.protein != b.protein) return a.protein < b.protein;
  if (a.i != b.i) return a.i < b.i;
  return a.j < b.j;
}

struct Result
{
  std::size_t hypotheses = 0; // the windows whose mass equals the precursor mass (whatever their count)
  std::vector<Window> top;    // the best top_n windows with at least min_k matched positions, best first
};

// Every window (i, j) inside a segment whose unmodified mass S(j) - S(i) + h2o equals the precursor mass (within precursor_tol_ppm of it,
// also at isotope errors -1 and +1, isotope_mass apart), found by a two-pointer sweep over the residue sums S of each segment
// (sums[protein]). A ladder position x in (i, j) counts once per side: on the prefix side if S(x) - S(i) matches a prefix set, on the
// suffix side if S(j) - S(x) matches a suffix set (any of the run's ion types of that side); k is the number of such (position, side)
// matches, and the whole window (the precursor) is no position. With one ion type per side (b y, the prototype's calibrated setting)
// this is the count of every set. Several types of one side (b y c z) would count a position up to twice per side, with correlated
// chance matches (the b and c ions of one fragment of another proteoform), so k >= min_k could come from one or two positions. A window
// with k >= min_k is scored by the chance probability of k, lambda = (j - i - 1) / V x the tolerance widths of all sets up to
// V = S(j) - S(i); with several types per side this bounds the chance rate of a position on a side from above. Windows of fewer than
// min_length residues are not considered (nor counted as hypotheses).
inline void search(const std::vector<TargetSet>& sets, const std::vector<Segment>& segments, const std::vector<const double*>& sums,
                   double precursor_mass, double precursor_tol_ppm, double h2o, double isotope_mass, int min_k, int min_length, std::size_t top_n,
                   Result& result)
{
  result.hypotheses = 0;
  result.top.clear();
  if (top_n == 0 || precursor_mass <= 0) return;
  const double tol = precursor_tol_ppm * 1e-6 * precursor_mass;
  const double low = precursor_mass - h2o - isotope_mass - tol, high = precursor_mass - h2o + isotope_mass + tol;
  std::vector<int> counts(sets.size(), 0);
  for (const auto& seg : segments)
  {
    const double* S = sums[seg.protein];
    if (S[seg.b] - S[seg.a] < low) continue; // no window of the segment is that heavy (residue sums rise within a segment)
    std::size_t jp = seg.a + 1;
    for (std::size_t i = seg.a; i < seg.b; i++)
    {
      const double si = S[i];
      if (S[seg.b] - si < low) break; // nor from this start or any later one
      if (jp < i + 1) jp = i + 1;
      while (jp <= seg.b && S[jp] - si < low) jp++;
      for (std::size_t j = jp; j <= seg.b && S[j] - si <= high; j++)
      {
        if ((int)(j - i) < min_length) continue; // not a hypothesis either
        const double mass = S[j] - si;
        int isotope = 2; // none
        for (int v = -1; v <= 1; v++)
        {
          if (std::abs(mass + h2o - (precursor_mass + v * isotope_mass)) <= tol) isotope = v;
        }
        if (isotope == 2) continue;
        result.hypotheses++;
        const double sj = S[j];
        std::fill(counts.begin(), counts.end(), 0);
        int k = 0;
        for (std::size_t x = i + 1; x < j; x++)
        {
          bool on_prefix = false, on_suffix = false;
          for (std::size_t s = 0; s < sets.size(); s++)
          {
            const bool prefix = sets[s].isPrefix();
            if (sets[s].massIndex(prefix ? S[x] - si : sj - S[x]) < 0) continue;
            counts[s]++;
            (prefix ? on_prefix : on_suffix) = true;
          }
          k += (int)on_prefix + (int)on_suffix;
        }
        if (k < min_k) continue;
        double weight = 0;
        for (const auto& set : sets) weight += set.weightUpTo(mass);
        Window w;
        w.lambda = mass > 0 ? (double)(j - i - 1) / mass * weight : 0;
        w.log_p = -std::log10(std::max(boost::math::gamma_p((double)k, std::max(w.lambda, 1e-12)), 1e-300));
        w.k = k;
        w.protein = seg.protein;
        w.i = (std::uint32_t)i;
        w.j = (std::uint32_t)j;
        w.isotope = isotope;
        if (result.top.size() >= top_n && ! isBetter(w, result.top.back())) continue;
        w.counts = counts;
        result.top.insert(std::upper_bound(result.top.begin(), result.top.end(), w, isBetter), std::move(w));
        if (result.top.size() > top_n) result.top.pop_back();
      }
    }
  }
}
} // namespace precursor_ladder
// END precursor ladder kernel
} // namespace

void MassLadderAlgorithm::buildIndex(const std::vector<std::string>& protein_seqs, const std::vector<std::string>& ion_types, bool terminal_ladders,
                                     bool precursor_ladders)
{
  ion_residue_shifts_.clear();
  protein_lengths_.clear();
  ladder_pre_.clear();
  ladder_suf_.clear();
  ladder_pre_index_.clear();
  ladder_suf_index_.clear();
  n_ladder_shifts_.clear();
  c_ladder_shifts_.clear();
  segments_.clear();
  window_sums_.clear();
  window_own_sums_.clear();

  // per ion type: whether it is a prefix ion, and the shift that makes its fragment mass a residue sum (mass + shift)
  for (const auto& ion : ion_types)
  {
    double offset = 0;
    bool prefix = true;
    if (ion == "a") offset = Residue::getInternalToAIon().getMonoWeight();
    else if (ion == "b") offset = Residue::getInternalToBIon().getMonoWeight();
    else if (ion == "c") offset = Residue::getInternalToCIon().getMonoWeight();
    else
    {
      prefix = false;
      if (ion == "x") offset = Residue::getInternalToXIon().getMonoWeight();
      else if (ion == "y") offset = Residue::getInternalToYIon().getMonoWeight();
      else if (ion == "z") offset = Residue::getInternalToZIon().getMonoWeight();
      else if (ion == "zp1") offset = Residue::getInternalToZp1Ion().getMonoWeight();
      else if (ion == "zp2") offset = Residue::getInternalToZp2Ion().getMonoWeight();
    }
    ion_residue_shifts_.emplace_back(prefix, -offset);
  }
  protein_lengths_.reserve(protein_seqs.size());
  for (const auto& seq : protein_seqs) protein_lengths_.push_back(seq.size());
  if (! terminal_ladders && ! precursor_ladders) return;

  // Terminal ladders: the residue sums of every protein from its N- and from its C-terminus (up to the first non-standard residue), and
  // all of them sorted by mass. A deconvolved mass plus a shift (of its ion type and, at the N-terminus, of the forms without the
  // initiator Met and with acetylation) that lands on a residue sum reads from the intact terminus; no sequence tag is needed. The sums
  // from the N-terminus are also built for the precursor ladders alone: they are the window sums of a protein of standard residues only.
  const auto* rdb = ResidueDB::getInstance();
  const std::string& standard = standard_residues;
  auto sums = [&](const std::string& seq, bool reverse, std::vector<double>& out) {
    out.assign(1, 0.0);
    for (Size j = 0; j < seq.size(); j++)
    {
      const char aa = seq[reverse ? seq.size() - 1 - j : j];
      if (standard.find(aa) == std::string::npos) break;
      out.push_back(out.back() + rdb->getResidue(aa)->getMonoWeight(Residue::Internal));
    }
  };
  ladder_pre_.resize(protein_seqs.size());
  if (terminal_ladders) ladder_suf_.resize(protein_seqs.size());
  for (Size k = 0; k < protein_seqs.size(); k++)
  {
    sums(protein_seqs[k], false, ladder_pre_[k]);
    if (! terminal_ladders) continue;
    sums(protein_seqs[k], true, ladder_suf_[k]);
    for (Size j = 1; j < ladder_pre_[k].size(); j++) ladder_pre_index_.push_back({ladder_pre_[k][j], (std::uint32_t)k});
    for (Size j = 1; j < ladder_suf_[k].size(); j++) ladder_suf_index_.push_back({ladder_suf_[k][j], (std::uint32_t)k});
  }
  if (terminal_ladders)
  {
    auto by_mass = [](const LadderEntry& a, const LadderEntry& b) { return a.mass < b.mass; };
    std::sort(ladder_pre_index_.begin(), ladder_pre_index_.end(), by_mass);
    std::sort(ladder_suf_index_.begin(), ladder_suf_index_.end(), by_mass);

    const double met = rdb->getResidue('M')->getMonoWeight(Residue::Internal);
    const double acetyl = EmpiricalFormula("C2H2O").getMonoWeight();
    for (const auto& [prefix, shift] : ion_residue_shifts_)
    {
      if (prefix)
      {
        for (double form : {0.0, met, -acetyl, met - acetyl}) n_ladder_shifts_.push_back(form + shift); // form - offset
      }
      else { c_ladder_shifts_.push_back(shift); }
    }
  }

  // Precursor-anchored ladders: the database as runs of standard residues (segments) with their residue sums, for the window search of
  // every spectrum with a precursor mass. A protein of standard residues only uses its sums from the N-terminus (ladder_pre_, no copy); a
  // protein with other residues gets its own sums, in which those add 0 (no window crosses them).
  if (precursor_ladders)
  {
    std::vector<Size> other;
    for (Size k = 0; k < protein_seqs.size(); k++)
    {
      if (protein_seqs[k].find_first_not_of(standard_residues) != std::string::npos) other.push_back(k);
    }
    window_sums_.assign(protein_seqs.size(), nullptr);
    window_own_sums_.resize(other.size()); // not resized below: the pointers into it stay valid
    for (Size n = 0; n < other.size(); n++)
    {
      const auto& seq = protein_seqs[other[n]];
      auto& own = window_own_sums_[n];
      own.reserve(seq.size() + 1);
      own.push_back(0);
      for (const char aa : seq)
      {
        own.push_back(own.back() + (standard_residues.find(aa) == std::string::npos ? 0.0 : rdb->getResidue(aa)->getMonoWeight(Residue::Internal)));
      }
      window_sums_[other[n]] = own.data();
    }
    for (Size k = 0; k < protein_seqs.size(); k++)
    {
      if (window_sums_[k] == nullptr) window_sums_[k] = ladder_pre_[k].data(); // standard residues only: the sums of the whole protein
      const auto& seq = protein_seqs[k];
      long a = -1;
      for (Size x = 0; x <= seq.size(); x++)
      {
        const bool is_standard = x < seq.size() && standard_residues.find(seq[x]) != std::string::npos;
        if (is_standard && a < 0) a = (long)x;
        if (! is_standard && a >= 0)
        {
          segments_.push_back({(std::uint32_t)k, (std::uint32_t)a, (std::uint32_t)x});
          a = -1;
        }
      }
    }
  }
}

// Terminal ladders: per protein and shift, the masses landing on a residue sum from the terminus. A protein is scored as
// FLASHTaggerAlgorithm::runMatchingExact scores a side (the chance probability of the match count, times the number of shifts of the
// side), without the terminus itself: its match would only tell proteins starting with Met from shuffled decoys. The top_n best proteins
// are the candidates, each with the ladder of the best shift of each side.
void MassLadderAlgorithm::searchTerminalLadders(const std::vector<double>& masses, double tol_ppm, Size top_n, const std::set<Size>& excluded,
                                                std::vector<Candidate>& candidates) const
{
  candidates.clear();
  if (masses.empty()) return;
  const Size n_proteins = protein_lengths_.size(), n_n = n_ladder_shifts_.size(), n_shifts = n_n + c_ladder_shifts_.size();
  auto shift_of = [&](Size si) { return si < n_n ? n_ladder_shifts_[si] : c_ladder_shifts_[si - n_n]; };
  std::vector<std::uint16_t> counts(n_proteins * n_shifts, 0);
  std::vector<Size> touched;
  for (Size si = 0; si < n_shifts; si++)
  {
    const auto& index = si < n_n ? ladder_pre_index_ : ladder_suf_index_;
    for (double m : masses)
    {
      const double t = m + shift_of(si);
      const double margin = tol_ppm * 1e-6 * std::max(m, t);
      auto it = std::lower_bound(index.begin(), index.end(), t - margin, [](const LadderEntry& e, double v) { return e.mass < v; });
      for (; it != index.end() && it->mass <= t + margin; ++it)
      {
        const Size c = it->protein * n_shifts + si;
        if (counts[c]++ == 0) touched.push_back(c);
      }
    }
  }
  // per shift: the masses' cumulative tolerance weight, for the expected number of chance matches on a protein
  std::vector<std::vector<double>> cum_weight(n_shifts, std::vector<double>(masses.size() + 1, 0.0));
  for (Size si = 0; si < n_shifts; si++)
  {
    for (Size i = 0; i < masses.size(); i++) cum_weight[si][i + 1] = cum_weight[si][i] + std::max(masses[i], masses[i] + shift_of(si));
  }
  std::map<Size, std::pair<double, double>> side_scores; // protein -> best N side, best C side (-log10 P)
  std::map<Size, std::pair<Size, Size>> best_shifts;      // protein -> its best N and C shift
  for (Size c : touched)
  {
    const int k = counts[c];
    if (k < 2) continue;
    const Size p = c / n_shifts, si = c % n_shifts;
    const auto& sums = si < n_n ? ladder_pre_[p] : ladder_suf_[p];
    const double s = shift_of(si);
    const Size lo = std::lower_bound(masses.begin(), masses.end(), -s) - masses.begin();
    const Size hi = std::upper_bound(masses.begin(), masses.end(), sums.back() - s) - masses.begin();
    const double lambda = hi > lo ? 2.0 * tol_ppm * 1e-6 * (double)(sums.size() - 1) / sums.back() * (cum_weight[si][hi] - cum_weight[si][lo]) : 0;
    const double score = -std::log10(std::max(boost::math::gamma_p((double)k, std::max(lambda, 1e-12)), 1e-300));
    auto& sc = side_scores[p];
    auto& bs = best_shifts.emplace(p, std::make_pair(n_shifts, n_shifts)).first->second;
    if (si < n_n && score > sc.first) { sc.first = score; bs.first = si; }
    else if (si >= n_n && score > sc.second) { sc.second = score; bs.second = si; }
  }
  std::vector<std::pair<double, Size>> ranked; // score, protein
  for (const auto& [p, sc] : side_scores)
  {
    if (excluded.count(p) > 0) continue;
    const double score = std::max(n_n > 0 ? sc.first - std::log10((double)n_n) : 0.0, n_shifts > n_n ? sc.second - std::log10((double)(n_shifts - n_n)) : 0.0);
    if (score > 0) ranked.emplace_back(score, p);
  }
  const Size top = std::min(ranked.size(), top_n);
  std::partial_sort(ranked.begin(), ranked.begin() + top, ranked.end(),
                    [](const auto& a, const auto& b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
  for (Size r = 0; r < top; r++)
  {
    Candidate cand;
    cand.protein = ranked[r].second;
    cand.rank = (int)r + 1;
    cand.score = ranked[r].first;
    const Size p = cand.protein, length = protein_lengths_[p];
    const Size bn = best_shifts[p].first, bc = best_shifts[p].second;
    for (Size si : {bn, bc})
    {
      if (si >= n_shifts) continue;
      const bool prefix = si < n_n;
      const auto& sums = prefix ? ladder_pre_[p] : ladder_suf_[p];
      const double s = shift_of(si);
      Ladder ladder;
      ladder.prefix = prefix;
      int first = -1, last = -1; // residues from the terminus at the lowest and the highest ladder mass
      for (Size mi = 0; mi < masses.size(); mi++)
      {
        const double m = masses[mi], t = m + s, margin = tol_ppm * 1e-6 * std::max(m, t);
        auto it = std::lower_bound(sums.begin() + 1, sums.end(), t - margin);
        if (it == sums.end() || *it > t + margin) continue;
        const int j = (int)(it - sums.begin());
        // the whole sequence is no fragment: that mass is the precursor (y) or the precursor minus water (b)
        if (j <= last || j >= (int)length) continue;
        if (first < 0) first = j;
        last = j;
        ladder.mass_indices.push_back(mi);
      }
      if (ladder.mass_indices.size() < 2) continue;
      ladder.start = prefix ? first : (int)length - last;
      ladder.end = prefix ? last : (int)length - first;
      cand.ladders.push_back(std::move(ladder));
    }
    if (cand.ladders.empty()) continue;
    candidates.push_back(std::move(cand));
  }
}

// Precursor-anchored ladders: the windows of the database (targets and decoys alike) whose unmodified mass is the precursor mass, by the
// chance probability of the masses on their prefix and suffix ladders (precursor_ladder::search). Each of the top_n best windows is a
// candidate with the matched masses of its best prefix and its best suffix ion type as ladders, spanning them.
void MassLadderAlgorithm::searchPrecursorLadders(const std::vector<double>& masses, const std::vector<double>& precursor_masses, double fragment_tol_ppm,
                                                 double precursor_tol_ppm, Size top_n, std::vector<PrecursorLadderResult>& results) const
{
  results.assign(precursor_masses.size(), PrecursorLadderResult());
  if (masses.empty()) return;
  const double h2o = Residue::getInternalToFull().getMonoWeight();
  std::vector<precursor_ladder::TargetSet> sets;
  sets.reserve(ion_residue_shifts_.size());
  for (const auto& [prefix, shift] : ion_residue_shifts_) sets.emplace_back(masses, shift, prefix, fragment_tol_ppm);
  precursor_ladder::Result found;
  for (Size n = 0; n < precursor_masses.size(); n++)
  {
    precursor_ladder::search(sets, segments_, window_sums_, precursor_masses[n], precursor_tol_ppm, h2o, Constants::C13C12_MASSDIFF_U,
                             precursor_ladder_min_k, precursor_ladder_min_length, top_n, found);
    auto& result = results[n];
    result.windows = found.hypotheses;
    const double log_hypotheses = std::log10((double)std::max(found.hypotheses, (std::size_t)1));
    for (Size r = 0; r < found.top.size(); r++)
    {
      const auto& win = found.top[r];
      const double* sums = window_sums_[win.protein];
      // per side, the ion type with the most matched positions (the first in ion_type among equal ones)
      int best_prefix = -1, best_suffix = -1;
      for (int s = 0; s < (int)sets.size(); s++)
      {
        int& best = sets[s].isPrefix() ? best_prefix : best_suffix;
        if (best < 0 || win.counts[s] > win.counts[best]) best = s;
      }
      Candidate cand;
      for (int s : {best_prefix, best_suffix})
      {
        if (s < 0) continue;
        const auto& set = sets[s];
        const bool prefix = set.isPrefix();
        Ladder ladder;
        ladder.prefix = prefix;
        std::vector<int> positions; // ascending masses: prefix positions from the start, suffix positions from the end
        for (Size step = 1; step < win.j - win.i; step++)
        {
          const Size x = prefix ? win.i + step : win.j - step;
          const long mi = set.massIndex(prefix ? sums[x] - sums[win.i] : sums[win.j] - sums[x]);
          if (mi < 0) continue;
          positions.push_back((int)x);
          ladder.mass_indices.push_back((Size)mi);
        }
        if (ladder.mass_indices.size() < 2) continue;
        // the ladder spans the positions of its lowest and its highest mass; its first mass is the lowest (the ladder read from i or j)
        ladder.start = prefix ? positions.front() : positions.back();
        ladder.end = prefix ? positions.back() : positions.front();
        cand.ladders.push_back(std::move(ladder));
      }
      if (cand.ladders.empty()) continue; // only possible with more than one ion type per side
      cand.protein = win.protein;
      cand.rank = (int)r + 1;
      cand.score = win.log_p - log_hypotheses;
      cand.window_start = (int)win.i;
      cand.window_end = (int)win.j;
      cand.isotope = win.isotope;
      result.candidates.push_back(std::move(cand));
    }
  }
}

// Sliding ladders: per protein, every start i of a prefix ladder (masses = S(x) - S(i) on prefix ion types) and every end j of a suffix
// ladder (masses = S(j) - S(x) on suffix ion types) inside the runs of standard residues. A position x counts once per side, whatever ion
// types match it, as in a precursor-ladder window; the anchor itself and the far end of the run (the whole rest of the run, no fragment)
// are no positions. An anchor with k >= precursor_ladder_min_k matched positions is scored by the chance probability of k, lambda =
// positions / V x the tolerance widths of the side's ion types up to V, the residue sum of its farthest position below the largest
// mass. A side of a protein scores the better of its best terminus anchor (the protein start, also without the first residue, or the
// protein end) minus log10 of the terminus anchors tested and its best other anchor minus log10 of the other anchors tested, and the
// protein its better side. The top_n proteins are the candidates, each with the ladder of the best ion type at the best anchor of each
// side (at least two masses).
void MassLadderAlgorithm::searchSlidingLadders(const std::vector<double>& masses, double tol_ppm, Size top_n, const std::set<Size>& excluded,
                                               std::vector<Candidate>& candidates) const
{
  candidates.clear();
  if (masses.empty() || segments_.empty() || top_n == 0) return;
  std::vector<precursor_ladder::TargetSet> sets;
  sets.reserve(ion_residue_shifts_.size());
  for (const auto& [prefix, shift] : ion_residue_shifts_) sets.emplace_back(masses, shift, prefix, tol_ppm);
  double max_value[2] = {0, 0}; // per side (0 prefix, 1 suffix): the largest residue sum a mass can match
  for (const auto& set : sets) max_value[set.isPrefix() ? 0 : 1] = std::max(max_value[set.isPrefix() ? 0 : 1], set.maxValue());

  struct Side
  {
    double log_p = 0;    // the best anchor's -log10 P
    long anchor = -1;    // the best anchor (i of a prefix, j of a suffix ladder)
    std::uint32_t a = 0, b = 0; // the run of standard residues of the best anchor
    Size tested = 0;     // the anchors tested on this side of the protein
  };
  // Per side two classes of anchors, each corrected for its own number of tests: the intact termini (prefix: the protein start, also
  // without the initiator residue; suffix: the protein end) and every other start or end, so that a ladder from an intact terminus does
  // not pay for every start of its protein. Slot side * 2 + class (0 other, 1 terminus).
  // Per segment (its own slots, so segments are searched in parallel when the caller is not parallel itself, e.g. one native spectrum),
  // then per protein: the best anchor over its segments (the first among equal ones) and the anchors tested on all of them.
  std::vector<std::array<Side, 4>> seg_best(segments_.size());
#pragma omp parallel for schedule(dynamic, 256) if (! omp_in_parallel())
  for (long si = 0; si < (long)segments_.size(); si++)
  {
    const auto& seg = segments_[si];
    if (excluded.count(seg.protein) > 0) continue;
    const double* S = window_sums_[seg.protein];
    for (int side = 0; side < 2; side++)
    {
      if (max_value[side] <= 0) continue; // no ion type of this side
      const bool prefix = side == 0;
      // prefix: anchors i = a .. b-2, positions x = i+1 .. b-1; suffix: anchors j = a+2 .. b, positions x = j-1 .. a+1
      for (long anchor = prefix ? (long)seg.a : (long)seg.a + 2; prefix ? anchor + 2 <= (long)seg.b : anchor <= (long)seg.b; anchor++)
      {
        const bool terminus = prefix ? anchor <= 1 : anchor == (long)protein_lengths_[seg.protein];
        auto& b = seg_best[si][side * 2 + (terminus ? 1 : 0)];
        b.tested++;
        int k = 0, positions = 0;
        double v_far = 0;
        for (long x = prefix ? anchor + 1 : anchor - 1; prefix ? x < (long)seg.b : x > (long)seg.a; x += prefix ? 1 : -1)
        {
          const double v = prefix ? S[x] - S[anchor] : S[anchor] - S[x];
          if (v > max_value[side]) break;
          positions++;
          v_far = v;
          for (const auto& set : sets)
          {
            if (set.isPrefix() != prefix || set.massIndex(v) < 0) continue;
            k++;
            break;
          }
        }
        if (k < precursor_ladder_min_k) continue;
        double weight = 0;
        for (const auto& set : sets)
        {
          if (set.isPrefix() == prefix) weight += set.weightUpTo(v_far);
        }
        const double lambda = (double)positions / v_far * weight;
        const double log_p = -std::log10(std::max(boost::math::gamma_p((double)k, std::max(lambda, 1e-12)), 1e-300));
        if (log_p > b.log_p) // the first anchor among equal ones
        {
          b.log_p = log_p;
          b.anchor = anchor;
          b.a = seg.a;
          b.b = seg.b;
        }
      }
    }
  }
  std::vector<std::array<Side, 4>> best(protein_lengths_.size());
  for (Size si = 0; si < segments_.size(); si++)
  {
    for (int side = 0; side < 4; side++)
    {
      auto& b = best[segments_[si].protein][side];
      const auto& sb = seg_best[si][side];
      b.tested += sb.tested;
      if (sb.anchor >= 0 && (b.anchor < 0 || sb.log_p > b.log_p))
      {
        const Size tested = b.tested;
        b = sb;
        b.tested = tested;
      }
    }
  }

  // per protein and side, the class of the better corrected score (its best anchor is the ladder of that side)
  auto corrected = [](const Side& b) { return b.anchor < 0 ? 0.0 : b.log_p - std::log10((double)std::max(b.tested, (Size)1)); };
  std::vector<std::array<int, 2>> chosen(best.size(), {0, 0});
  std::vector<std::pair<double, Size>> ranked; // score, protein
  for (Size p = 0; p < best.size(); p++)
  {
    double score = 0;
    for (int side = 0; side < 2; side++)
    {
      const int c = corrected(best[p][side * 2 + 1]) > corrected(best[p][side * 2]) ? 1 : 0;
      chosen[p][side] = side * 2 + c;
      score = std::max(score, corrected(best[p][side * 2 + c]));
    }
    if (score > 0) ranked.emplace_back(score, p);
  }
  const Size top = std::min(ranked.size(), top_n);
  std::partial_sort(ranked.begin(), ranked.begin() + top, ranked.end(),
                    [](const auto& a, const auto& b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
  for (Size r = 0; r < top; r++)
  {
    const Size p = ranked[r].second;
    const double* S = window_sums_[p];
    Candidate cand;
    for (int side = 0; side < 2; side++)
    {
      const auto& b = best[p][chosen[p][side]];
      if (b.anchor < 0) continue;
      const bool prefix = side == 0;
      const long seg_a = b.a, seg_b = b.b;
      // the ion type of this side with the most matched positions at the anchor (the first in ion_type among equal ones)
      Ladder ladder;
      ladder.prefix = prefix;
      std::vector<int> positions;
      for (const auto& set : sets)
      {
        if (set.isPrefix() != prefix) continue;
        Ladder l;
        l.prefix = prefix;
        std::vector<int> pos;
        for (long x = prefix ? b.anchor + 1 : b.anchor - 1; prefix ? x < seg_b : x > seg_a; x += prefix ? 1 : -1)
        {
          const double v = prefix ? S[x] - S[b.anchor] : S[b.anchor] - S[x];
          if (v > set.maxValue()) break;
          const long mi = set.massIndex(v);
          if (mi < 0) continue;
          pos.push_back((int)x);
          l.mass_indices.push_back((Size)mi);
        }
        if (l.mass_indices.size() > ladder.mass_indices.size())
        {
          ladder = std::move(l);
          positions = std::move(pos);
        }
      }
      if (ladder.mass_indices.size() < 2) continue;
      // as a precursor-ladder window side: the ladder spans the positions of its lowest and its highest mass, its first mass the lowest
      ladder.start = prefix ? positions.front() : positions.back();
      ladder.end = prefix ? positions.back() : positions.front();
      cand.ladders.push_back(std::move(ladder));
    }
    if (cand.ladders.empty()) continue;
    cand.protein = p;
    cand.rank = (int)r + 1;
    cand.score = ranked[r].first;
    candidates.push_back(std::move(cand));
  }
}
} // namespace OpenMS
