// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHGappedTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <algorithm>
#include <bit>
#include <cmath>
#include <numeric>
#include <unordered_map>

namespace OpenMS
{
namespace
{
const char standard_residues[] = "ACDEFGHIKLMNPQRSTVWY";
/// glycine, the lightest residue: no gap is lighter
constexpr double min_residue_mass = 57.02146;
/// residue masses are rounded to the 0.1 mDa grid of the plausibility set; their sum over up to max_gap_mass / 57 residues can be off by
/// this much (35 x 0.05 mDa), so a plausibility window is widened by it
constexpr double plausible_margin = 0.002;

/// dst |= src << shift (bits), both of the same length
void orShifted(std::vector<std::uint64_t>& dst, const std::vector<std::uint64_t>& src, std::int64_t shift)
{
  const Size n = dst.size();
  const Size word_shift = (Size)(shift >> 6);
  const int bit_shift = (int)(shift & 63);
  if (word_shift >= n) return;
  for (Size i = n; i-- > word_shift;)
  {
    const Size j = i - word_shift;
    std::uint64_t v = src[j] << bit_shift;
    if (bit_shift != 0 && j > 0) v |= src[j - 1] >> (64 - bit_shift);
    dst[i] |= v;
  }
}

/// values of sorted @p a found in sorted @p b (exponential then binary search of each value of the shorter list in the longer one),
/// appended to @p out in ascending order
void intersectSorted(const std::uint32_t* a, Size na, const std::uint32_t* b, Size nb, std::vector<std::uint32_t>& out)
{
  if (na == 0 || nb == 0) return;
  if (na > nb)
  {
    std::swap(a, b);
    std::swap(na, nb);
  }
  Size lo = 0;
  for (Size i = 0; i < na && lo < nb; ++i)
  {
    const std::uint32_t v = a[i];
    if (b[lo] > v) continue;
    Size step = 1;
    while (lo + step < nb && b[lo + step] < v)
    {
      lo += step;
      step <<= 1;
    }
    const Size hi = std::min(lo + step, nb - 1);
    lo = (Size)(std::lower_bound(b + lo, b + hi + 1, v) - b);
    if (lo < nb && b[lo] == v)
    {
      out.push_back(v);
      ++lo;
    }
  }
}
} // namespace

FLASHGappedTaggerAlgorithm::FLASHGappedTaggerAlgorithm() : DefaultParamHandler("FLASHGappedTaggerAlgorithm")
{
  defaults_.setValue("max_gap_mass", 2000.0,
                     "Largest residue mass a gapped tag spans, from its lowest to its highest mass (Da); the protein windows are indexed up "
                     "to this mass.");
  defaults_.setMinFloat("max_gap_mass", 100.0);
  defaults_.setValue("min_tag_size", 4, "Smallest number of masses in a gapped tag (with 3, a tag lands on a decoy about as often as on the right protein).");
  defaults_.setMinInt("min_tag_size", 3);
  defaults_.setValue("max_tag_size", 7, "Largest number of masses in a gapped tag.");
  defaults_.setMinInt("max_tag_size", 3);
  defaults_.setValue("beam_width", 100, "Number of gapped tags of each size proposed per spectrum, highest summed Qscore first.");
  defaults_.setMinInt("beam_width", 1);
  defaults_.setValue("top_n", 100, "Number of masses per spectrum (highest Qscore) used for gapped tags.");
  defaults_.setMinInt("top_n", 3);
  defaults_.setValue("stop_at_first_size", "false",
                     "Use only the largest tag size that lands anywhere (as in GappedFiltration); false: every size, and the protein score ranks "
                     "the candidates.");
  defaults_.setValidStrings("stop_at_first_size", {"true", "false"});
  defaults_.setValue("bin_width", 0.1, "Mass bin width of the protein window index (Da). A tolerance window is looked up in every bin it touches.");
  defaults_.setMinFloat("bin_width", 0.001);
  defaults_.addTag("bin_width", "advanced");
  defaultsToParam_();
}

void FLASHGappedTaggerAlgorithm::updateMembers_()
{
  max_gap_mass_ = param_.getValue("max_gap_mass");
  min_tag_size_ = param_.getValue("min_tag_size");
  max_tag_size_ = std::max(min_tag_size_, (int)param_.getValue("max_tag_size"));
  beam_width_ = param_.getValue("beam_width");
  top_n_ = param_.getValue("top_n");
  stop_at_first_size_ = param_.getValue("stop_at_first_size").toString() == "true";
  bin_width_ = param_.getValue("bin_width");
}

Size FLASHGappedTaggerAlgorithm::getIndexBytes() const
{
  return postings_.size() * sizeof(std::uint32_t) + post_offsets_.size() * sizeof(std::uint64_t) + slot_.size() * sizeof(std::uint32_t)
         + cum_mass_.size() * sizeof(double) + valid_.size() + plausible_bits_.size() * sizeof(std::uint64_t)
         + plausible_rank_.size() * sizeof(std::uint32_t);
}

void FLASHGappedTaggerAlgorithm::buildIndex(const std::vector<std::string>& protein_seqs)
{
  double aa_mass[128] = {0};
  for (const char* c = standard_residues; *c; ++c)
  {
    aa_mass[(int)*c] = ResidueDB::getInstance()->getResidue(*c)->getMonoWeight(Residue::Internal);
  }

  // residues concatenated
  Size total = 0;
  for (const auto& s : protein_seqs) total += s.size();
  cum_mass_.assign(total + 1, 0);
  valid_.assign(total, 0);
  offsets_.assign(1, 0);
  offsets_.reserve(protein_seqs.size() + 1);
  Size pos = 0;
  for (const auto& s : protein_seqs)
  {
    for (char ch : s)
    {
      const int c = std::toupper((unsigned char)ch);
      const double m = c < 128 ? aa_mass[c] : 0;
      valid_[pos] = m > 0;
      cum_mass_[pos + 1] = cum_mass_[pos] + m;
      ++pos;
    }
    offsets_.push_back(pos);
  }

  // every window of 1 .. max_length residues, binned by mass; a non-standard residue ends the windows containing it, and windows
  // never cross proteins. Start positions are walked in the outer loop, so each bin receives its positions in ascending order.
  key_min_ = 0;
  key_max_ = key_(max_gap_mass_);
  const int max_length = (int)(max_gap_mass_ / min_residue_mass);
  std::vector<std::uint64_t> counts((Size)key_max_ + 1, 0), cursor;
  // first pass counts the windows per bin, the second fills the posting lists
  auto enumerate = [&](bool fill) {
    for (Size p = 0; p + 1 < offsets_.size(); ++p)
    {
      const Size begin = offsets_[p], end = offsets_[p + 1];
      for (Size start = begin; start < end; ++start)
      {
        if (! valid_[start]) continue;
        const Size last = std::min(end, start + (Size)max_length);
        for (Size stop = start + 1; stop <= last; ++stop)
        {
          if (! valid_[stop - 1]) break;
          const std::int64_t key = key_(cum_mass_[stop] - cum_mass_[start]);
          if (key > key_max_) break;
          if (fill) { postings_[cursor[slot_[(Size)key]]++] = (std::uint32_t)start; }
          else { ++counts[(Size)key]; }
        }
      }
    }
  };
  enumerate(false);
  slot_.assign((Size)key_max_ + 1, std::numeric_limits<std::uint32_t>::max());
  post_offsets_.assign(1, 0);
  for (Size key = 0; key < counts.size(); ++key)
  {
    if (counts[key] == 0) continue;
    slot_[key] = (std::uint32_t)(post_offsets_.size() - 1);
    post_offsets_.push_back(post_offsets_.back() + counts[key]);
  }
  postings_.assign(post_offsets_.back(), 0);
  cursor.assign(post_offsets_.begin(), post_offsets_.end() - 1);
  enumerate(true);

  // sums of 1 .. max_length residues on a 0.1 mDa grid (subset sum by shifted bit sets), with the set bits counted before each word
  plausible_units_ = (std::int64_t)std::llround(max_gap_mass_ * plausible_scale_);
  const Size words = (Size)((plausible_units_ + 1 + 63) / 64);
  std::vector<std::int64_t> units;
  for (const char* c = standard_residues; *c; ++c)
  {
    units.push_back((std::int64_t)std::llround(aa_mass[(int)*c] * plausible_scale_));
  }
  std::sort(units.begin(), units.end());
  units.erase(std::unique(units.begin(), units.end()), units.end());
  std::vector<std::uint64_t> current(words, 0), next(words, 0);
  plausible_bits_.assign(words, 0);
  current[0] = 1;
  const Size tail = (Size)((plausible_units_ + 1) & 63);
  for (int len = 1; len <= max_length; ++len)
  {
    std::fill(next.begin(), next.end(), 0);
    for (auto u : units) orShifted(next, current, u);
    if (tail) next.back() &= (1ULL << tail) - 1;
    current.swap(next);
    bool any = false;
    for (Size i = 0; i < words; ++i)
    {
      plausible_bits_[i] |= current[i];
      any |= current[i] != 0;
    }
    if (! any) break;
  }
  plausible_rank_.assign(words + 1, 0);
  for (Size i = 0; i < words; ++i)
  {
    plausible_rank_[i + 1] = plausible_rank_[i] + (std::uint32_t)std::popcount(plausible_bits_[i]);
  }
}

bool FLASHGappedTaggerAlgorithm::plausible_(double lo, double hi) const
{
  std::int64_t lo_u = (std::int64_t)std::ceil(lo * plausible_scale_), hi_u = (std::int64_t)std::floor(hi * plausible_scale_);
  lo_u = std::max<std::int64_t>(lo_u, 0);
  hi_u = std::min(hi_u, plausible_units_);
  if (hi_u < lo_u) return false;
  // set bits before index i
  auto rank = [this](std::int64_t i) -> std::uint32_t {
    const Size w = (Size)(i >> 6);
    if (w >= plausible_bits_.size()) return plausible_rank_.back();
    const int b = (int)(i & 63);
    const std::uint64_t mask = b == 0 ? 0 : (~0ULL >> (64 - b));
    return plausible_rank_[w] + (std::uint32_t)std::popcount(plausible_bits_[w] & mask);
  };
  return rank(hi_u + 1) > rank(lo_u);
}

Size FLASHGappedTaggerAlgorithm::countRange_(std::int64_t key_lo, std::int64_t key_hi) const
{
  Size n = 0;
  for (std::int64_t key = std::max(key_lo, key_min_); key <= std::min(key_hi, key_max_); ++key)
  {
    const auto s = slot_[(Size)key];
    if (s != std::numeric_limits<std::uint32_t>::max()) n += post_offsets_[s + 1] - post_offsets_[s];
  }
  return n;
}

void FLASHGappedTaggerAlgorithm::collect_(std::int64_t key_lo, std::int64_t key_hi, std::vector<std::uint32_t>& out) const
{
  out.clear();
  for (std::int64_t key = std::max(key_lo, key_min_); key <= std::min(key_hi, key_max_); ++key)
  {
    const auto s = slot_[(Size)key];
    if (s == std::numeric_limits<std::uint32_t>::max()) continue;
    const Size before = out.size();
    out.insert(out.end(), postings_.begin() + (std::ptrdiff_t)post_offsets_[s], postings_.begin() + (std::ptrdiff_t)post_offsets_[s + 1]);
    if (before) std::inplace_merge(out.begin(), out.begin() + (std::ptrdiff_t)before, out.end());
  }
}

void FLASHGappedTaggerAlgorithm::search(const DeconvolvedSpectrum& dspec, double ppm, std::vector<Landing>& landings) const
{
  landings.clear();
  if (postings_.empty()) return;
  const double tol = ppm * 1e-6;

  // masses with the highest Qscore, ascending in mass
  struct Node
  {
    double mass;
    double qscore;
    int score;
  };
  std::vector<Node> nodes;
  nodes.reserve(dspec.size());
  for (const auto& pg : dspec)
  {
    nodes.push_back({pg.getMonoMass(), pg.getQscore(), FLASHTaggerAlgorithm::getNodeScore(pg)});
  }
  if ((int)nodes.size() > top_n_)
  {
    std::nth_element(nodes.begin(), nodes.begin() + top_n_, nodes.end(),
                     [](const Node& a, const Node& b) { return a.qscore != b.qscore ? a.qscore > b.qscore : a.mass < b.mass; });
    nodes.resize(top_n_);
  }
  std::sort(nodes.begin(), nodes.end(), [](const Node& a, const Node& b) { return a.mass < b.mass; });
  const int M = (int)nodes.size();
  if (M < min_tag_size_) return;

  // mass graph: an edge i -> j when the difference is a sum of residues within the tolerance of both masses
  std::vector<std::vector<int>> fwd(M);
  for (int i = 0; i < M; ++i)
  {
    for (int j = i + 1; j < M; ++j)
    {
      const double gap = nodes[j].mass - nodes[i].mass;
      const double d = tol * (nodes[i].mass + nodes[j].mass);
      if (gap + d < min_residue_mass) continue;
      if (gap - d > max_gap_mass_) break;
      if (! plausible_(gap - d - plausible_margin, gap + d + plausible_margin)) continue;
      fwd[i].push_back(j);
    }
  }
  // longest path (in edges) from each mass: a tag of n masses can start at i only if it is at least n - 1
  std::vector<int> fdepth(M, 0);
  for (int i = M - 1; i >= 0; --i)
  {
    int best = -1;
    for (int j : fwd[i]) best = std::max(best, fdepth[j]);
    fdepth[i] = best + 1;
  }

  // landing per (global start position, series): the largest tag, then the highest Qscore sum
  std::unordered_map<std::uint64_t, Landing> best_landing;
  std::vector<std::uint32_t> result, partial, operand;
  std::vector<int> rows, cand_rows, order;
  std::vector<double> score, cand_score;
  std::vector<std::pair<std::int64_t, std::int64_t>> ranges;
  std::vector<Size> sizes;

  // keep the `width` best rows (score, then row order) of a row-major block, preserving row order
  auto keepTop = [&](std::vector<int>& block, std::vector<double>& sc, int cols) {
    const int n_rows = (int)sc.size();
    if (n_rows <= beam_width_) return;
    order.resize(n_rows);
    std::iota(order.begin(), order.end(), 0);
    std::nth_element(order.begin(), order.begin() + beam_width_, order.end(),
                     [&](int a, int b) { return sc[a] != sc[b] ? sc[a] > sc[b] : a < b; });
    order.resize(beam_width_);
    std::sort(order.begin(), order.end());
    for (int r = 0; r < beam_width_; ++r)
    {
      if (order[r] == r) continue;
      std::copy(block.begin() + (std::ptrdiff_t)order[r] * cols, block.begin() + (std::ptrdiff_t)(order[r] + 1) * cols,
                block.begin() + (std::ptrdiff_t)r * cols);
      sc[r] = sc[order[r]];
    }
    block.resize((Size)beam_width_ * cols);
    sc.resize(beam_width_);
  };

  // positions (global, ascending) where windows start with a mass in every key range
  auto intersectRanges = [&]() {
    result.clear();
    sizes.resize(ranges.size());
    for (Size k = 0; k < ranges.size(); ++k)
    {
      sizes[k] = countRange_(ranges[k].first, ranges[k].second);
      if (sizes[k] == 0) return;
    }
    std::vector<Size> idx(ranges.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](Size a, Size b) { return sizes[a] < sizes[b]; });
    collect_(ranges[idx[0]].first, ranges[idx[0]].second, result);
    for (Size k = 1; k < idx.size() && ! result.empty(); ++k)
    {
      partial.clear();
      const auto [lo, hi] = ranges[idx[k]];
      for (std::int64_t key = std::max(lo, key_min_); key <= std::min(hi, key_max_); ++key)
      {
        const auto s = slot_[(Size)key];
        if (s == std::numeric_limits<std::uint32_t>::max()) continue;
        // a start position has at most one window per bin of a tolerance window (lengths differ by a residue), so the parts are disjoint
        intersectSorted(result.data(), result.size(), postings_.data() + post_offsets_[s], (Size)(post_offsets_[s + 1] - post_offsets_[s]),
                        partial);
      }
      std::sort(partial.begin(), partial.end());
      result.swap(partial);
    }
  };

  const int max_length = (int)(max_gap_mass_ / min_residue_mass);
  for (int n = std::min(max_tag_size_, M); n >= min_tag_size_; --n)
  {
    // beam over paths of n masses, highest summed Qscore first; a path spans at most max_gap_mass
    rows.clear();
    score.clear();
    for (int i = 0; i < M; ++i)
    {
      if (fdepth[i] < n - 1) continue;
      rows.push_back(i);
      score.push_back(nodes[i].qscore);
    }
    if (rows.empty()) continue;
    keepTop(rows, score, 1);
    for (int k = 2; k <= n && ! score.empty(); ++k)
    {
      const int cols = k - 1;
      cand_rows.clear();
      cand_score.clear();
      for (Size r = 0; r < score.size(); ++r)
      {
        const int first = rows[r * cols], last = rows[r * cols + cols - 1];
        for (int child : fwd[last])
        {
          const double span = nodes[child].mass - nodes[first].mass;
          if (span - tol * (nodes[child].mass + nodes[first].mass) > max_gap_mass_) break;
          if (fdepth[child] < n - k) continue;
          cand_rows.insert(cand_rows.end(), rows.begin() + (std::ptrdiff_t)r * cols, rows.begin() + (std::ptrdiff_t)(r + 1) * cols);
          cand_rows.push_back(child);
          cand_score.push_back(score[r] + nodes[child].qscore);
        }
      }
      rows.swap(cand_rows);
      score.swap(cand_score);
      keepTop(rows, score, k);
    }
    if (score.empty()) continue;

    bool landed = false;
    for (Size r = 0; r < score.size(); ++r)
    {
      const int* path = rows.data() + r * n;
      const double base = nodes[path[0]].mass, top = nodes[path[n - 1]].mass;
      for (int series = 0; series < 2; ++series)
      {
        const bool prefix = series == 0;
        // offsets of the masses from the lowest (prefix) or the highest (suffix); every offset is a window starting at one position
        ranges.clear();
        for (int j = 1; j < n; ++j)
        {
          const double m = prefix ? nodes[path[j]].mass : nodes[path[n - 1 - j]].mass;
          const double off = prefix ? m - base : top - m;
          const double d = tol * (m + (prefix ? base : top));
          ranges.emplace_back(key_(off - d), key_(off + d));
        }
        intersectRanges();
        if (result.empty()) continue;
        const auto [span_lo, span_hi] = ranges.back();
        for (std::uint32_t p : result)
        {
          const Size protein = (Size)(std::upper_bound(offsets_.begin(), offsets_.end(), (Size)p) - offsets_.begin()) - 1;
          const Size protein_end = offsets_[protein + 1];
          // the residue length of the span offset
          int length = 0;
          for (int len = 1; len <= max_length && p + len <= protein_end && valid_[p + len - 1]; ++len)
          {
            const std::int64_t key = key_(cum_mass_[p + len] - cum_mass_[p]);
            if (key > span_hi) break;
            if (key >= span_lo)
            {
              length = len;
              break;
            }
          }
          if (length == 0) continue;
          const std::uint64_t id = ((std::uint64_t)p << 1) | (prefix ? 1 : 0);
          auto it = best_landing.find(id);
          double qsum = score[r];
          if (it != best_landing.end() && ((int)it->second.masses.size() > n || it->second.qscore_sum >= qsum)) continue;
          Landing l;
          l.protein = protein;
          l.start = (int)(p - offsets_[protein]);
          l.end = l.start + length;
          l.prefix = prefix;
          l.flank_mass = prefix ? cum_mass_[p] - cum_mass_[offsets_[protein]] : cum_mass_[protein_end] - cum_mass_[p + length];
          for (int j = 0; j < n; ++j)
          {
            l.masses.push_back(nodes[path[j]].mass);
            l.scores.push_back(nodes[path[j]].score);
          }
          l.qscore_sum = qsum;
          best_landing[id] = std::move(l);
          landed = true;
        }
      }
    }
    if (landed && stop_at_first_size_) break;
  }

  landings.reserve(best_landing.size());
  for (auto& [id, l] : best_landing) landings.push_back(std::move(l));
  std::sort(landings.begin(), landings.end(), [](const Landing& a, const Landing& b) {
    return a.protein != b.protein ? a.protein < b.protein : (a.start != b.start ? a.start < b.start : a.prefix > b.prefix);
  });
}
} // namespace OpenMS
