// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/MassLadderAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/ProteoformDatabase.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <boost/math/special_functions/gamma.hpp>
#include <algorithm>
#include <array>
#include <bit>
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

  // whether a mass matches the residue sum v: massIndex(v) >= 0
  bool matches(double v) const
  {
    if (! (v <= max_value_)) return false;
    const std::size_t b = (std::size_t)(v * bin_inverse_ + bin_phase_);
    if (((bits_[b >> 6] >> (b & 63)) & 1) == 0) return false;
    const std::size_t idx = (std::size_t)(std::lower_bound(hi_.begin(), hi_.end(), v) - hi_.begin());
    return idx < hi_.size() && lo_[idx] <= v;
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

// a precursor mass: a window of residue sum m fits it when low <= m <= high and m + h2o is within tol of it at an isotope error of -1, 0
// or 1 (isotopeOf)
struct Query
{
  double mass = 0, tol = 0, low = 0, high = 0;
};

inline Query makeQuery(double precursor_mass, double precursor_tol_ppm, double h2o, double isotope_mass)
{
  Query q;
  q.mass = precursor_mass;
  q.tol = precursor_tol_ppm * 1e-6 * precursor_mass;
  q.low = precursor_mass - h2o - isotope_mass - q.tol;
  q.high = precursor_mass - h2o + isotope_mass + q.tol;
  return q;
}

// the isotope error (-1, 0, 1; the last of several) at which a window of residue sum m fits the precursor mass; 2 if none
inline int isotopeOf(double mass, const Query& q, double h2o, double isotope_mass)
{
  int isotope = 2; // none
  for (int v = -1; v <= 1; v++)
  {
    if (std::abs(mass + h2o - (q.mass + v * isotope_mass)) <= q.tol) isotope = v;
  }
  return isotope;
}

// A hypothesis: the window i..j-1 of a protein (S its residue sums, mass = S(j) - S(i)) that fits the precursor mass at the isotope
// error isotope. Counts its matched ladder positions and keeps it among the top_n best windows of result. A window that another entry
// of the protein owns (owners, ProteoformDatabase) is no hypothesis of this one.
inline void scoreWindow(const std::vector<TargetSet>& sets, const double* S, std::uint32_t protein, std::size_t i, std::size_t j, double mass,
                        int isotope, int min_k, std::size_t top_n, std::vector<int>& counts, Result& result,
                        const MassLadderAlgorithm::WindowOwnership* owners)
{
  if (owners != nullptr && ! owners->owns(protein, i, j)) return;
  result.hypotheses++;
  const double si = S[i], sj = S[j];
  std::fill(counts.begin(), counts.end(), 0);
  int k = 0;
  if (sets.size() == 2 && sets[0].isPrefix() != sets[1].isPrefix())
  {
    // one ion type per side (b y): k is the sum of their counts
    const std::size_t p = sets[0].isPrefix() ? 0 : 1;
    const TargetSet &prefix_set = sets[p], &suffix_set = sets[1 - p];
    int on_prefix = 0, on_suffix = 0;
    for (std::size_t x = i + 1; x < j; x++)
    {
      on_prefix += (int)prefix_set.matches(S[x] - si);
      on_suffix += (int)suffix_set.matches(sj - S[x]);
    }
    counts[p] = on_prefix;
    counts[1 - p] = on_suffix;
    k = on_prefix + on_suffix;
  }
  else
  {
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
  }
  if (k < min_k) return;
  double weight = 0;
  for (const auto& set : sets) weight += set.weightUpTo(mass);
  Window w;
  w.lambda = mass > 0 ? (double)(j - i - 1) / mass * weight : 0;
  w.log_p = -std::log10(std::max(boost::math::gamma_p((double)k, std::max(w.lambda, 1e-12)), 1e-300));
  w.k = k;
  w.protein = protein;
  w.i = (std::uint32_t)i;
  w.j = (std::uint32_t)j;
  w.isotope = isotope;
  if (result.top.size() >= top_n && ! isBetter(w, result.top.back())) return;
  w.counts = counts;
  result.top.insert(std::upper_bound(result.top.begin(), result.top.end(), w, isBetter), std::move(w));
  if (result.top.size() > top_n) result.top.pop_back();
}

// Every window (i, j) inside a segment whose unmodified mass S(j) - S(i) + h2o equals the precursor mass (within precursor_tol_ppm of it,
// also at isotope errors -1 and +1, isotope_mass apart), found by a two-pointer sweep over the residue sums S of each segment
// (sums[protein]). A ladder position x in (i, j) counts once per side: on the prefix side if S(x) - S(i) matches a prefix set, on the
// suffix side if S(j) - S(x) matches a suffix set (any of the run's ion types of that side); k is the number of such (position, side)
// matches, and the whole window (the precursor) is no position. With one ion type per side (b y, the prototype's calibrated setting)
// this is the count of every set. Several types of one side (b y c z) would count a position up to twice per side, with correlated
// chance matches (the b and c ions of one fragment of another proteoform), so k >= min_k could come from one or two positions. A window
// with k >= min_k is scored by the chance probability of k, lambda = (j - i - 1) / V x the tolerance widths of all sets up to
// V = S(j) - S(i); with several types per side this bounds the chance rate of a position on a side from above. Windows of fewer than
// min_length residues are not considered (nor counted as hypotheses), nor are those without every fixed modification of their entry
// (i > i_max or j < j_min of the segment).
inline void search(const std::vector<TargetSet>& sets, const std::vector<Segment>& segments, const std::vector<const double*>& sums,
                   double precursor_mass, double precursor_tol_ppm, double h2o, double isotope_mass, int min_k, int min_length, std::size_t top_n,
                   Result& result, const MassLadderAlgorithm::WindowOwnership* owners)
{
  result.hypotheses = 0;
  result.top.clear();
  if (top_n == 0 || precursor_mass <= 0) return;
  const Query q = makeQuery(precursor_mass, precursor_tol_ppm, h2o, isotope_mass);
  const double low = q.low, high = q.high;
  std::vector<int> counts(sets.size(), 0);
  for (const auto& seg : segments)
  {
    const double* S = sums[seg.protein];
    if (S[seg.b] - S[seg.a] < low) continue; // no window of the segment is that heavy (residue sums rise within a segment)
    std::size_t jp = seg.a + 1;
    for (std::size_t i = seg.a; i < seg.b; i++)
    {
      if (i > seg.i_max) break; // a later start leaves a fixed modification out
      const double si = S[i];
      if (S[seg.b] - si < low) break; // nor from this start or any later one
      if (jp < i + 1) jp = i + 1;
      while (jp <= seg.b && S[jp] - si < low) jp++;
      for (std::size_t j = jp; j <= seg.b && S[j] - si <= high; j++)
      {
        if ((int)(j - i) < min_length) continue; // not a hypothesis either
        if (j < seg.j_min) continue;
        const double mass = S[j] - si;
        const int isotope = isotopeOf(mass, q, h2o, isotope_mass);
        if (isotope == 2) continue;
        scoreWindow(sets, S, seg.protein, i, j, mass, isotope, min_k, top_n, counts, result, owners);
      }
    }
  }
}

// The same as search for the windows of the precursor mass listed by findWindows (the hypotheses of the sweep, in its order): those of
// chunk c are windows[begin[c] .. end[c] - 1].
inline void searchListed(const std::vector<TargetSet>& sets, const std::vector<Segment>& segments, const std::vector<const double*>& sums,
                         const MassLadderAlgorithm::PreparedWindow* windows, const std::uint32_t* begin, const std::uint32_t* end, std::size_t chunks,
                         double precursor_mass, double precursor_tol_ppm, double h2o, double isotope_mass, int min_k, std::size_t top_n, Result& result,
                         const MassLadderAlgorithm::WindowOwnership* owners)
{
  result.hypotheses = 0;
  result.top.clear();
  if (top_n == 0 || precursor_mass <= 0) return;
  const Query q = makeQuery(precursor_mass, precursor_tol_ppm, h2o, isotope_mass);
  std::vector<int> counts(sets.size(), 0);
  for (std::size_t c = 0; c < chunks; c++)
  {
    for (std::uint32_t t = begin[c]; t < end[c]; t++)
    {
      const auto& seg = segments[windows[t].segment];
      const double* S = sums[seg.protein];
      const std::size_t i = (std::size_t)seg.a + windows[t].start, j = i + windows[t].length;
      const double mass = S[j] - S[i];
      scoreWindow(sets, S, seg.protein, i, j, mass, isotopeOf(mass, q, h2o, isotope_mass), min_k, top_n, counts, result, owners);
    }
  }
}

// The hypotheses of search (the windows that fit a precursor mass) for many precursor masses (queries, ascending and distinct) in one
// pass over the database instead of one sweep per mass. A bitset of 1/32 Da bins of the residue sums that can fit a query (its three
// isotope intervals, widened far beyond rounding errors) pre-filters the windows of at least min_length residues up to the largest such
// sum; a window in a set bin is tested against the queries whose [low, high] can hold its mass, exactly as search tests it. Pass 1 counts
// the hypotheses of each query per chunk of consecutive segments, pass 2 writes them grouped by query, each group in the order of search
// (segment, i, j) over the chunks. Only the ascending queries whose hypotheses fit max_windows in all are listed (prepared).
// Returns false (nothing listed) if the windows are not listed: a segment of 65536 residues or more, residue sums beyond 10 MDa or
// isotope intervals of more than 2^25 bins in all (precursor tolerances far beyond those of mass spectrometers).
inline bool findWindows(const std::vector<Segment>& segments, const std::vector<const double*>& sums, const std::vector<Query>& queries,
                        double h2o, double isotope_mass, int min_length, std::size_t max_windows, std::size_t& prepared, std::size_t& chunks_out,
                        std::vector<MassLadderAlgorithm::PreparedWindow>& windows, std::vector<std::uint32_t>& offsets)
{
  prepared = 0;
  chunks_out = 0;
  windows.clear();
  offsets.clear();
  const std::size_t n_queries = queries.size();
  if (n_queries == 0 || segments.empty()) return true;
  if (n_queries >= ((std::size_t)1 << 31)) return false;
  max_windows = std::min(max_windows, ((std::size_t)1 << 31) - 1);
  for (const auto& seg : segments)
  {
    if (seg.b - seg.a > 65535) return false;
  }
  // the bins of the residue sums that can fit a query: a sum m that passes isotopeOf at v is within tol of mass + v x isotope_mass - h2o
  // up to rounding errors of a few ulps (about 1e-15 x mass), far below pad
  constexpr double bin_inverse = 32.0;
  auto interval = [&](const Query& q, int v) {
    const double center = q.mass + v * isotope_mass - h2o, pad = 1e-9 * q.mass + 1e-6;
    return std::make_pair(center - q.tol - pad, center + q.tol + pad);
  };
  double max_value = 0; // the largest residue sum that can fit a query
  for (const auto& q : queries)
  {
    for (int v = -1; v <= 1; v++) max_value = std::max(max_value, interval(q, v).second);
  }
  if (! (max_value <= 1e7)) return false;
  auto binRange = [&](const Query& q, int v, std::size_t& b0, std::size_t& b1) {
    const auto [lo, hi] = interval(q, v);
    if (! (hi > 0)) return false;
    b0 = lo > 0 ? (std::size_t)(lo * bin_inverse) : 0;
    b1 = (std::size_t)(hi * bin_inverse);
    return true;
  };
  std::size_t entries = 0; // bins per query, at most
  for (const auto& q : queries)
  {
    for (int v = -1; v <= 1; v++)
    {
      std::size_t b0, b1;
      if (binRange(q, v, b0, b1)) entries += b1 - b0 + 1;
    }
  }
  if (entries > ((std::size_t)1 << 25)) return false;
  const std::size_t bin_count = (std::size_t)(max_value * bin_inverse) + 1;
  std::vector<std::uint64_t> bits(bin_count / 64 + 1, 0);
  for (const auto& q : queries)
  {
    for (int v = -1; v <= 1; v++)
    {
      std::size_t b0, b1;
      if (! binRange(q, v, b0, b1)) continue;
      for (std::size_t b = b0; b <= b1; b++) bits[b >> 6] |= std::uint64_t(1) << (b & 63);
    }
  }
  // per set bin (by its rank among them), the queries that can fit a sum in it (bin_queries[bin_offset[r] .. bin_offset[r + 1] - 1])
  std::vector<std::uint32_t> word_rank(bits.size() + 1, 0); // the set bins in the words before
  for (std::size_t w = 0; w < bits.size(); w++) word_rank[w + 1] = word_rank[w] + (std::uint32_t)std::popcount(bits[w]);
  auto rankOf = [&](std::size_t b) { return word_rank[b >> 6] + (std::uint32_t)std::popcount(bits[b >> 6] & ((std::uint64_t(1) << (b & 63)) - 1)); };
  const std::size_t set_bins = word_rank.back();
  std::vector<std::uint32_t> bin_offset(set_bins + 1, 0), last(set_bins, (std::uint32_t)n_queries), bin_queries;
  for (int fill = 0; fill < 2; fill++)
  {
    if (fill == 1)
    {
      for (std::size_t r = 0; r < set_bins; r++) bin_offset[r + 1] += bin_offset[r];
      bin_queries.resize(bin_offset.back());
      std::fill(last.begin(), last.end(), (std::uint32_t)n_queries);
    }
    std::vector<std::uint32_t> next(bin_offset.begin(), bin_offset.end() - 1);
    for (std::size_t n = 0; n < n_queries; n++)
    {
      for (int v = -1; v <= 1; v++)
      {
        std::size_t b0, b1;
        if (! binRange(queries[n], v, b0, b1)) continue;
        for (std::size_t b = b0; b <= b1; b++)
        {
          const std::uint32_t r = rankOf(b);
          if (last[r] == n) continue; // the isotope intervals of a query overlap
          last[r] = (std::uint32_t)n;
          if (fill == 0) bin_offset[r + 1]++;
          else bin_queries[next[r]++] = (std::uint32_t)n;
        }
      }
    }
  }
  // the windows of the segments s0..s1-1 tested against the queries that can fit them, in the order of search: f(hypothesis of the query,
  // query, segment, i, j); hypotheses are tested without branches (most tests pass, but which ones is unpredictable)
  const std::size_t j_min = (std::size_t)std::max(min_length, 1);
  auto enumerate = [&](std::size_t s0, std::size_t s1, auto&& f) {
    for (std::size_t s = s0; s < s1; s++)
    {
      const auto& seg = segments[s];
      const double* S = sums[seg.protein];
      for (std::size_t i = seg.a; i + j_min <= seg.b && i <= seg.i_max; i++)
      {
        const double si = S[i];
        for (std::size_t j = std::max(i + j_min, (std::size_t)seg.j_min); j <= seg.b; j++)
        {
          const double mass = S[j] - si;
          if (mass > max_value) break; // residue sums rise within a segment
          const std::size_t bin = (std::size_t)(mass * bin_inverse);
          if (((bits[bin >> 6] >> (bin & 63)) & 1) == 0) continue;
          const std::uint32_t r = rankOf(bin);
          for (std::uint32_t t = bin_offset[r]; t < bin_offset[r + 1]; t++)
          {
            const Query& q = queries[bin_queries[t]];
            f((mass >= q.low) & (mass <= q.high) & (isotopeOf(mass, q, h2o, isotope_mass) != 2), bin_queries[t], s, i, j);
          }
        }
      }
    }
  };
  // chunks of consecutive segments of about equal work (windows up to max_value: at most max_value / 57 residues long)
  int threads = 1;
#ifdef _OPENMP
  threads = omp_in_parallel() ? 1 : omp_get_max_threads();
#endif
  const std::size_t chunk_count = std::min(segments.size(), (std::size_t)(16 * std::max(threads, 1)));
  const double max_length = max_value / 57.0 + 1;
  double total_work = 0;
  for (const auto& seg : segments) total_work += (double)(seg.b - seg.a) * std::min((double)(seg.b - seg.a), max_length);
  std::vector<std::size_t> chunk_start {0};
  double work = 0;
  for (std::size_t s = 0; s < segments.size(); s++)
  {
    work += (double)(segments[s].b - segments[s].a) * std::min((double)(segments[s].b - segments[s].a), max_length);
    if (work >= total_work * (double)chunk_start.size() / (double)chunk_count && chunk_start.size() < chunk_count) chunk_start.push_back(s + 1);
  }
  if (chunk_start.back() != segments.size()) chunk_start.push_back(segments.size());
  const std::size_t chunks = chunk_start.size() - 1;
  // pass 1: the hypotheses per chunk and query (counts above max_windows are not needed)
  std::vector<std::uint32_t> table(chunks * n_queries, 0);
#pragma omp parallel for schedule(dynamic, 1) if (threads > 1)
  for (long c = 0; c < (long)chunks; c++)
  {
    std::uint32_t* count = table.data() + (std::size_t)c * n_queries;
    enumerate(chunk_start[c], chunk_start[c + 1], [&](bool fits, std::uint32_t n, std::size_t, std::size_t, std::size_t) {
      count[n] += (std::uint32_t)(fits & (count[n] <= max_windows));
    });
  }
  // the queries listed: from the lowest mass while their hypotheses fit max_windows
  std::size_t total = 0;
  for (; prepared < n_queries; prepared++)
  {
    std::size_t sum = 0;
    for (std::size_t c = 0; c < chunks; c++) sum += table[c * n_queries + prepared];
    if (total + sum > max_windows) break;
    total += sum;
  }
  if (prepared == 0) return true;
  // the windows by chunk, in a chunk by query (so that a chunk writes to its own part only): offsets[n * chunks + c] is the first of
  // query n in chunk c, offsets[prepared * chunks + c] the end of chunk c
  chunks_out = chunks;
  offsets.assign((prepared + 1) * chunks, 0);
  std::size_t position = 0;
  for (std::size_t c = 0; c < chunks; c++)
  {
    for (std::size_t n = 0; n < prepared; n++)
    {
      offsets[n * chunks + c] = (std::uint32_t)position;
      position += table[c * n_queries + n];
      table[c * n_queries + n] = offsets[n * chunks + c];
    }
    offsets[prepared * chunks + c] = (std::uint32_t)position;
  }
  // pass 2: the hypotheses of the listed queries
  windows.resize(total);
#pragma omp parallel for schedule(dynamic, 1) if (threads > 1)
  for (long c = 0; c < (long)chunks; c++)
  {
    std::uint32_t* next = table.data() + (std::size_t)c * n_queries;
    enumerate(chunk_start[c], chunk_start[c + 1], [&](bool fits, std::uint32_t n, std::size_t s, std::size_t i, std::size_t j) {
      if (! fits || n >= prepared) return;
      windows[next[n]++] = {(std::uint32_t)s, (std::uint16_t)(i - segments[s].a), (std::uint16_t)(j - i)};
    });
  }
  return true;
}
} // namespace precursor_ladder
// END precursor ladder kernel
} // namespace

namespace
{
// For a terminal-ladder index sorted by mass (all masses above 0): first[b] is the first entry at or above b Da, for b = 0 .. the largest
// mass + 1, so that the entries at or above a mass v in (0, largest] lie from first[floor(v)] on, and from first[floor(v) + 1] on all are
// above v. Not built (empty) for an index of 2^32 entries or more or a mass beyond 10 MDa.
template<typename Entry> void indexByDalton(const std::vector<Entry>& index, std::vector<std::uint32_t>& first)
{
  first.clear();
  if (index.empty() || index.size() >= (Size)UINT32_MAX || ! (index.front().mass > 0) || ! (index.back().mass <= 1e7)) return;
  const Size buckets = (Size)index.back().mass + 2;
  first.resize(buckets);
  Size i = 0;
  for (Size b = 0; b < buckets; b++)
  {
    while (i < index.size() && index[i].mass < (double)b) i++;
    first[b] = (std::uint32_t)i;
  }
}
} // namespace

void MassLadderAlgorithm::buildIndex(const std::vector<std::string>& protein_seqs, const std::vector<std::string>& ion_types, bool terminal_ladders,
                                     bool precursor_ladders, const ProteoformDatabase* db)
{
  ion_residue_shifts_.clear();
  protein_lengths_.clear();
  entry_protein_.clear();
  entry_fixed_count_.clear();
  entry_offset_.clear();
  ladder_pre_fixed_.clear();
  ladder_suf_fixed_.clear();
  window_fixed_sums_.clear();
  window_owners_ = WindowOwnership();
  ladder_pre_.clear();
  ladder_suf_.clear();
  ladder_pre_index_.clear();
  ladder_suf_index_.clear();
  ladder_pre_first_.clear();
  ladder_suf_first_.clear();
  n_ladder_shifts_.clear();
  c_ladder_shifts_.clear();
  segments_.clear();
  window_sums_.clear();
  window_own_sums_.clear();
  clearPrecursorWindows();

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
  // a database of entries (ProteoformDatabase): their fixed modifications (per entry, (position, mass)), and the protein of each entry
  const std::vector<std::vector<std::pair<int, double>>>* fixed = nullptr;
  if (db != nullptr && ! db->isIdentity() && db->entries.size() == protein_seqs.size())
  {
    entry_protein_ = db->proteinOfEntry();
    entry_fixed_count_.reserve(db->entries.size());
    entry_offset_.reserve(db->entries.size());
    for (const auto& entry : db->entries)
    {
      entry_fixed_count_.push_back(entry.fixed.size());
      entry_offset_.push_back(entry.offset);
    }
    if (db->hasFixedModifications()) fixed = &db->fixedMasses();
  }
  if (! terminal_ladders && ! precursor_ladders) return;

  // Terminal ladders: the residue sums of every protein from its N- and from its C-terminus (up to the first non-standard residue), and
  // all of them sorted by mass. A deconvolved mass plus a shift (of its ion type and, at the N-terminus, of the forms without the
  // initiator Met and with acetylation) that lands on a residue sum reads from the intact terminus; no sequence tag is needed. The sums
  // from the N-terminus are also built for the precursor ladders alone: they are the window sums of a protein of standard residues only.
  const auto* rdb = ResidueDB::getInstance();
  const std::string& standard = standard_residues;
  // the fixed modification masses of an entry by residue (empty: none)
  auto fixedByResidue = [&](Size k, std::vector<double>& extra) {
    extra.clear();
    if (fixed == nullptr || (*fixed)[k].empty()) return;
    extra.assign(protein_seqs[k].size(), 0.0);
    for (const auto& [pos, mass] : (*fixed)[k])
    {
      if (pos >= 0 && pos < (int)extra.size()) extra[pos] += mass;
    }
  };
  auto sums = [&](const std::string& seq, bool reverse, const std::vector<double>& extra, std::vector<double>& out) {
    out.assign(1, 0.0);
    for (Size j = 0; j < seq.size(); j++)
    {
      const Size pos = reverse ? seq.size() - 1 - j : j;
      const char aa = seq[pos];
      if (standard.find(aa) == std::string::npos) break;
      if (extra.empty()) out.push_back(out.back() + rdb->getResidue(aa)->getMonoWeight(Residue::Internal));
      else out.push_back(out.back() + rdb->getResidue(aa)->getMonoWeight(Residue::Internal) + extra[pos]);
    }
  };
  ladder_pre_.resize(protein_seqs.size());
  if (terminal_ladders) ladder_suf_.resize(protein_seqs.size());
  if (fixed != nullptr)
  {
    ladder_pre_fixed_.assign(protein_seqs.size(), 0.0);
    ladder_suf_fixed_.assign(protein_seqs.size(), 0.0);
  }
  // the sums of the proteins in parallel (each its own), then the index entries of each protein at its place: the entries of the proteins
  // one after the other, as appended protein by protein
#pragma omp parallel for schedule(dynamic, 64) if (! omp_in_parallel())
  for (long k = 0; k < (long)protein_seqs.size(); k++)
  {
    std::vector<double> extra;
    fixedByResidue((Size)k, extra);
    sums(protein_seqs[k], false, extra, ladder_pre_[k]);
    if (terminal_ladders) sums(protein_seqs[k], true, extra, ladder_suf_[k]);
    if (extra.empty()) continue;
    // the fixed masses in the sums (up to the first non-standard residue from either end)
    const Size n = extra.size();
    for (Size j = 0; j + 1 < ladder_pre_[k].size(); j++) ladder_pre_fixed_[k] += extra[j];
    if (terminal_ladders)
    {
      for (Size j = 0; j + 1 < ladder_suf_[k].size(); j++) ladder_suf_fixed_[k] += extra[n - 1 - j];
    }
  }
  if (terminal_ladders)
  {
    std::vector<Size> pre_first(protein_seqs.size() + 1, 0), suf_first(protein_seqs.size() + 1, 0);
    for (Size k = 0; k < protein_seqs.size(); k++)
    {
      pre_first[k + 1] = pre_first[k] + ladder_pre_[k].size() - 1;
      suf_first[k + 1] = suf_first[k] + ladder_suf_[k].size() - 1;
    }
    ladder_pre_index_.resize(pre_first.back());
    ladder_suf_index_.resize(suf_first.back());
#pragma omp parallel for schedule(dynamic, 64) if (! omp_in_parallel())
    for (long k = 0; k < (long)protein_seqs.size(); k++)
    {
      for (Size j = 1; j < ladder_pre_[k].size(); j++) ladder_pre_index_[pre_first[k] + j - 1] = {ladder_pre_[k][j], (std::uint32_t)k};
      for (Size j = 1; j < ladder_suf_[k].size(); j++) ladder_suf_index_[suf_first[k] + j - 1] = {ladder_suf_[k][j], (std::uint32_t)k};
    }
    // the two sorts at the same time, each as before
    auto by_mass = [](const LadderEntry& a, const LadderEntry& b) { return a.mass < b.mass; };
#pragma omp parallel sections if (! omp_in_parallel())
    {
#pragma omp section
      std::sort(ladder_pre_index_.begin(), ladder_pre_index_.end(), by_mass);
#pragma omp section
      std::sort(ladder_suf_index_.begin(), ladder_suf_index_.end(), by_mass);
    }

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
    // per Da, the first entry of each index: searchTerminalLadders looks a mass up within its Da
    indexByDalton(ladder_pre_index_, ladder_pre_first_);
    indexByDalton(ladder_suf_index_, ladder_suf_first_);
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
      std::vector<double> extra;
      fixedByResidue(other[n], extra);
      own.reserve(seq.size() + 1);
      own.push_back(0);
      for (Size x = 0; x < seq.size(); x++)
      {
        const char aa = seq[x];
        const double mass = standard_residues.find(aa) == std::string::npos ? 0.0 : rdb->getResidue(aa)->getMonoWeight(Residue::Internal);
        if (extra.empty()) own.push_back(own.back() + mass);
        else own.push_back(own.back() + mass + (mass > 0 ? extra[x] : 0.0));
      }
      window_sums_[other[n]] = own.data();
    }
    if (fixed != nullptr) window_fixed_sums_.resize(protein_seqs.size());
    for (Size k = 0; k < protein_seqs.size(); k++)
    {
      if (window_sums_[k] == nullptr) window_sums_[k] = ladder_pre_[k].data(); // standard residues only: the sums of the whole protein
      const auto& seq = protein_seqs[k];
      // the fixed masses in the window sums up to each residue (a non-standard residue carries none)
      if (fixed != nullptr && ! (*fixed)[k].empty())
      {
        std::vector<double> extra;
        fixedByResidue(k, extra);
        auto& fs = window_fixed_sums_[k];
        fs.assign(1, 0.0);
        for (Size x = 0; x < seq.size(); x++) fs.push_back(fs.back() + (standard_residues.find(seq[x]) == std::string::npos ? 0.0 : extra[x]));
      }
      // a window contains every fixed modification of its entry: no other window is new (the same window of the entry without the
      // modifications outside it is one already)
      std::uint32_t i_max = UINT32_MAX, j_min = 0;
      if (fixed != nullptr && ! (*fixed)[k].empty())
      {
        i_max = (std::uint32_t)(*fixed)[k].front().first;
        j_min = (std::uint32_t)(*fixed)[k].back().first + 1;
      }
      long a = -1;
      for (Size x = 0; x <= seq.size(); x++)
      {
        const bool is_standard = x < seq.size() && standard_residues.find(seq[x]) != std::string::npos;
        if (is_standard && a < 0) a = (long)x;
        if (! is_standard && a >= 0)
        {
          segments_.push_back({(std::uint32_t)k, (std::uint32_t)a, (std::uint32_t)x, i_max, j_min});
          a = -1;
        }
      }
    }
    // windows of several entries of a protein with the same residues and fixed modifications are one hypothesis (WindowOwnership):
    // classes of the entries of a protein with the same fixed modifications at the same protein positions
    if (! entry_protein_.empty())
    {
      std::map<std::pair<Size, std::vector<std::pair<int, double>>>, std::vector<std::uint32_t>> classes;
      for (Size k = 0; k < protein_seqs.size(); k++)
      {
        std::vector<std::pair<int, double>> key;
        if (fixed != nullptr)
        {
          for (const auto& [pos, mass] : (*fixed)[k]) key.emplace_back(db->entries[k].offset + pos, mass);
        }
        classes[{entry_protein_[k], key}].push_back((std::uint32_t)k);
      }
      auto& wo = window_owners_;
      wo.offset.resize(protein_seqs.size());
      wo.length.resize(protein_seqs.size());
      wo.class_first.resize(protein_seqs.size());
      wo.class_size.resize(protein_seqs.size());
      for (Size k = 0; k < protein_seqs.size(); k++)
      {
        wo.offset[k] = db->entries[k].offset;
        wo.length[k] = (int)protein_seqs[k].size();
      }
      bool shared = false;
      for (auto& [key, members] : classes)
      {
        std::stable_sort(members.begin(), members.end(), [&](std::uint32_t l, std::uint32_t r) { return wo.length[l] < wo.length[r]; });
        shared |= members.size() > 1;
        for (std::uint32_t m : members)
        {
          wo.class_first[m] = (std::uint32_t)wo.members.size();
          wo.class_size[m] = (std::uint32_t)members.size();
        }
        wo.members.insert(wo.members.end(), members.begin(), members.end());
      }
      if (! shared) window_owners_ = WindowOwnership();
    }
  }
}

void MassLadderAlgorithm::bestEntryPerProtein_(std::vector<std::pair<double, Size>>& ranked, const std::vector<EntryEvidence>& evidence) const
{
  if (entry_protein_.empty()) return;
  auto simpler = [&](Size r, Size k) {
    const Size e = ranked[r].second, ek = ranked[k].second;
    return entry_fixed_count_[e] < entry_fixed_count_[ek] || (entry_fixed_count_[e] == entry_fixed_count_[ek] && e < ek);
  };
  std::map<Size, std::pair<Size, double>> best; // protein -> position in ranked, best score
  for (Size r = 0; r < ranked.size(); r++)
  {
    const auto [it, inserted] = best.emplace(entry_protein_[ranked[r].second], std::make_pair(r, ranked[r].first));
    if (inserted) continue;
    const Size b = it->second.first;
    if (evidence[r].score_sum > evidence[b].score_sum || (evidence[r].score_sum == evidence[b].score_sum && simpler(r, b))) it->second.first = r;
    it->second.second = std::max(it->second.second, ranked[r].first);
  }
  std::map<Size, Size> kept; // protein -> position in ranked
  for (const auto& [p, b] : best) kept[p] = b.first;
  for (Size r = 0; r < ranked.size(); r++)
  {
    const Size e = ranked[r].second, p = entry_protein_[e], b = best[p].first, eb = ranked[b].second;
    Size& k = kept[p];
    if (entry_offset_[e] == entry_offset_[eb] && protein_lengths_[e] == protein_lengths_[eb] && evidence[r].n == evidence[b].n
        && evidence[r].c == evidence[b].c && simpler(r, k))
      k = r;
  }
  std::vector<std::pair<double, Size>> out;
  out.reserve(kept.size());
  for (Size r = 0; r < ranked.size(); r++)
  {
    const Size p = entry_protein_[ranked[r].second];
    if (kept[p] == r) out.emplace_back(best[p].second, ranked[r].second);
  }
  ranked.swap(out);
}

void MassLadderAlgorithm::addEntriesOfTopProteins_(std::vector<std::pair<double, Size>>& ranked, Size& top,
                                                   const std::vector<std::pair<double, Size>>& all) const
{
  if (entry_protein_.empty()) return;
  std::map<Size, std::vector<std::pair<double, Size>>> others; // protein -> its other entries with evidence
  for (Size r = 0; r < top; r++) others[entry_protein_[ranked[r].second]];
  for (const auto& a : all)
  {
    auto it = others.find(entry_protein_[a.second]);
    if (it != others.end()) it->second.push_back(a);
  }
  std::vector<std::pair<double, Size>> out;
  for (Size r = 0; r < top; r++)
  {
    out.push_back(ranked[r]);
    for (const auto& a : others[entry_protein_[ranked[r].second]])
      if (a.second != ranked[r].second) out.push_back(a);
  }
  top = out.size();
  ranked.swap(out);
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
    const auto& first = si < n_n ? ladder_pre_first_ : ladder_suf_first_;
    auto below = [](const LadderEntry& e, double v) { return e.mass < v; };
    for (double m : masses)
    {
      const double t = m + shift_of(si);
      const double margin = tol_ppm * 1e-6 * std::max(m, t);
      const double v = t - margin;
      // the first entry at or above v (std::lower_bound of the whole index): within the Da of v (indexByDalton); none above the largest mass
      auto it = index.end();
      if (first.empty() || ! (v > 0)) { it = std::lower_bound(index.begin(), index.end(), v, below); }
      else if (v <= index.back().mass)
      {
        const Size b = (Size)v;
        it = std::lower_bound(index.begin() + first[b], index.begin() + first[b + 1], v, below);
      }
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
    // the residues of an entry with fixed modifications spread without them (FLASHTaggerAlgorithm::scoreShifts)
    const double spread = ladder_pre_fixed_.empty() ? sums.back() : sums.back() - (si < n_n ? ladder_pre_fixed_[p] : ladder_suf_fixed_[p]);
    const double lambda = hi > lo ? 2.0 * tol_ppm * 1e-6 * (double)(sums.size() - 1) / spread * (cum_weight[si][hi] - cum_weight[si][lo]) : 0;
    const double score = -std::log10(std::max(boost::math::gamma_p((double)k, std::max(lambda, 1e-12)), 1e-300));
    auto& sc = side_scores[p];
    auto& bs = best_shifts.emplace(p, std::make_pair(n_shifts, n_shifts)).first->second;
    if (si < n_n && score > sc.first) { sc.first = score; bs.first = si; }
    else if (si >= n_n && score > sc.second) { sc.second = score; bs.second = si; }
  }
  std::vector<std::pair<double, Size>> ranked; // score, protein
  std::vector<EntryEvidence> ranked_evidence; // with several entries per protein
  for (const auto& [p, sc] : side_scores)
  {
    if (excluded.count(p) > 0) continue;
    const double n_score = n_n > 0 ? sc.first - std::log10((double)n_n) : 0.0;
    const double c_score = n_shifts > n_n ? sc.second - std::log10((double)(n_shifts - n_n)) : 0.0;
    const double score = std::max(n_score, c_score);
    if (! (score > 0)) continue;
    ranked.emplace_back(score, p);
    if (entry_protein_.empty()) continue;
    const auto& bs = best_shifts.at(p);
    ranked_evidence.push_back({std::max(0.0, n_score) + std::max(0.0, c_score), bs.first < n_shifts ? (int)counts[p * n_shifts + bs.first] : 0,
                               bs.second < n_shifts ? (int)counts[p * n_shifts + bs.second] : 0});
  }
  const std::vector<std::pair<double, Size>> with_evidence = ranked; // every entry with evidence (several per protein: PEFF)
  bestEntryPerProtein_(ranked, ranked_evidence);
  Size top = std::min(ranked.size(), top_n);
  std::partial_sort(ranked.begin(), ranked.begin() + top, ranked.end(),
                    [](const auto& a, const auto& b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
  // the top proteins are the cut; each of their entries with evidence is a candidate (annotated residues are residues of the entry)
  addEntriesOfTopProteins_(ranked, top, with_evidence);
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
  const WindowOwnership* owners = window_owners_.empty() ? nullptr : &window_owners_;
  for (Size n = 0; n < precursor_masses.size(); n++)
  {
    // the windows of a prepared mass are listed (preparePrecursorWindows); any other mass is swept
    const auto prepared = std::lower_bound(prepared_masses_.begin(), prepared_masses_.end(), precursor_masses[n]);
    if (prepared != prepared_masses_.end() && *prepared == precursor_masses[n] && precursor_tol_ppm == prepared_tol_ppm_)
    {
      const Size q = (Size)(prepared - prepared_masses_.begin());
      precursor_ladder::searchListed(sets, segments_, window_sums_, prepared_windows_.data(), prepared_offsets_.data() + q * prepared_chunks_,
                                     prepared_offsets_.data() + (q + 1) * prepared_chunks_, prepared_chunks_, precursor_masses[n], precursor_tol_ppm,
                                     h2o, Constants::C13C12_MASSDIFF_U, precursor_ladder_min_k, top_n, found, owners);
    }
    else
    {
      precursor_ladder::search(sets, segments_, window_sums_, precursor_masses[n], precursor_tol_ppm, h2o, Constants::C13C12_MASSDIFF_U,
                               precursor_ladder_min_k, precursor_ladder_min_length, top_n, found, owners);
    }
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

// Precursor ladders of many spectra: the windows of all their precursor masses in one pass over the database (precursor_ladder::findWindows)
// instead of one sweep per mass, which costs more than the rest of the search. Masses up to 1 MDa at precursor tolerances up to 1000 ppm;
// at most 2^27 windows (1 GB) are listed, those of the lowest masses first. The other masses are swept.
void MassLadderAlgorithm::preparePrecursorWindows(const std::vector<double>& precursor_masses, double precursor_tol_ppm)
{
  clearPrecursorWindows();
  if (segments_.empty() || ! (precursor_tol_ppm >= 0 && precursor_tol_ppm <= 1000)) return;
  std::vector<double> masses;
  for (double m : precursor_masses)
  {
    if (m > 0 && m <= 1e6) masses.push_back(m);
  }
  std::sort(masses.begin(), masses.end());
  masses.erase(std::unique(masses.begin(), masses.end()), masses.end());
  if (masses.size() < 256) return; // the pass costs about as much as 150 sweeps (E. coli database), so a few masses are swept
  const double h2o = Residue::getInternalToFull().getMonoWeight();
  std::vector<precursor_ladder::Query> queries;
  queries.reserve(masses.size());
  for (double m : masses) queries.push_back(precursor_ladder::makeQuery(m, precursor_tol_ppm, h2o, Constants::C13C12_MASSDIFF_U));
  Size prepared = 0;
  if (! precursor_ladder::findWindows(segments_, window_sums_, queries, h2o, Constants::C13C12_MASSDIFF_U, precursor_ladder_min_length,
                                      (Size)1 << 27, prepared, prepared_chunks_, prepared_windows_, prepared_offsets_)
      || prepared == 0)
  {
    clearPrecursorWindows();
    return;
  }
  masses.resize(prepared);
  prepared_masses_ = std::move(masses);
  prepared_tol_ppm_ = precursor_tol_ppm;
}

void MassLadderAlgorithm::clearPrecursorWindows()
{
  prepared_tol_ppm_ = -1;
  std::vector<double>().swap(prepared_masses_);
  prepared_chunks_ = 0;
  std::vector<PreparedWindow>().swap(prepared_windows_);
  std::vector<std::uint32_t>().swap(prepared_offsets_);
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
    int k = 0;           // the best anchor's matched positions
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
    // an entry with fixed modifications: its window sums without them spread its positions (FLASHTaggerAlgorithm::scoreShifts)
    const double* F = window_fixed_sums_.empty() || window_fixed_sums_[seg.protein].empty() ? nullptr : window_fixed_sums_[seg.protein].data();
    // the second residue of an entry counts as its terminus (initiator Met removed) only for a form that starts at residue 1
    const bool met_start = entry_offset_.empty() || entry_offset_[seg.protein] == 0;
    for (int side = 0; side < 2; side++)
    {
      if (max_value[side] <= 0) continue; // no ion type of this side
      const bool prefix = side == 0;
      // prefix: anchors i = a .. b-2, positions x = i+1 .. b-1; suffix: anchors j = a+2 .. b, positions x = j-1 .. a+1
      for (long anchor = prefix ? (long)seg.a : (long)seg.a + 2; prefix ? anchor + 2 <= (long)seg.b : anchor <= (long)seg.b; anchor++)
      {
        const bool terminus = prefix ? (anchor == 0 || (anchor == 1 && met_start)) : anchor == (long)protein_lengths_[seg.protein];
        auto& b = seg_best[si][side * 2 + (terminus ? 1 : 0)];
        b.tested++;
        int k = 0, positions = 0;
        double v_far = 0, f_far = 0;
        for (long x = prefix ? anchor + 1 : anchor - 1; prefix ? x < (long)seg.b : x > (long)seg.a; x += prefix ? 1 : -1)
        {
          const double v = prefix ? S[x] - S[anchor] : S[anchor] - S[x];
          if (v > max_value[side]) break;
          positions++;
          v_far = v;
          if (F != nullptr) f_far = prefix ? F[x] - F[anchor] : F[anchor] - F[x];
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
        const double lambda = (double)positions / (F != nullptr ? v_far - f_far : v_far) * weight;
        const double log_p = -std::log10(std::max(boost::math::gamma_p((double)k, std::max(lambda, 1e-12)), 1e-300));
        if (log_p > b.log_p) // the first anchor among equal ones
        {
          b.log_p = log_p;
          b.anchor = anchor;
          b.a = seg.a;
          b.b = seg.b;
          b.k = k;
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
  std::vector<EntryEvidence> ranked_evidence;  // per side: its ladder's score and matched positions
  for (Size p = 0; p < best.size(); p++)
  {
    double score = 0;
    EntryEvidence evidence;
    for (int side = 0; side < 2; side++)
    {
      const int c = corrected(best[p][side * 2 + 1]) > corrected(best[p][side * 2]) ? 1 : 0;
      chosen[p][side] = side * 2 + c;
      score = std::max(score, corrected(best[p][side * 2 + c]));
      if (best[p][side * 2 + c].anchor < 0) continue;
      evidence.score_sum += std::max(0.0, corrected(best[p][side * 2 + c]));
      (side == 0 ? evidence.n : evidence.c) = best[p][side * 2 + c].k;
    }
    if (score > 0)
    {
      ranked.emplace_back(score, p);
      ranked_evidence.push_back(evidence);
    }
  }
  const std::vector<std::pair<double, Size>> with_evidence = ranked; // every entry with evidence (several per protein: PEFF)
  bestEntryPerProtein_(ranked, ranked_evidence);
  Size top = std::min(ranked.size(), top_n);
  std::partial_sort(ranked.begin(), ranked.begin() + top, ranked.end(),
                    [](const auto& a, const auto& b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
  // the top proteins are the cut; each of their entries with evidence is a candidate (annotated residues are residues of the entry)
  addEntriesOfTopProteins_(ranked, top, with_evidence);
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
