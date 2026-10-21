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
#include <unordered_map>
#ifdef _OPENMP
  #include <omp.h>
#endif

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

/// dst |= src << shift (bits) in the words [begin, end) of dst, both of the same length and distinct
void orShifted(std::vector<std::uint64_t>& dst, const std::vector<std::uint64_t>& src, std::int64_t shift, Size begin, Size end)
{
  const Size word_shift = (Size)(shift >> 6);
  const int bit_shift = (int)(shift & 63);
  for (Size i = std::max(begin, word_shift); i < end; ++i)
  {
    const Size j = i - word_shift;
    std::uint64_t v = src[j] << bit_shift;
    if (bit_shift != 0 && j > 0) v |= src[j - 1] >> (64 - bit_shift);
    dst[i] |= v;
  }
}
} // namespace

FLASHGappedTaggerAlgorithm::FLASHGappedTaggerAlgorithm() : DefaultParamHandler("FLASHGappedTaggerAlgorithm")
{
  defaults_.setValue("max_gap_mass", 5000.0,
                     "Largest residue mass a gapped tag spans, from its lowest to its highest mass (Da); the protein windows are indexed up "
                     "to this mass. Fragments of native top-down spectra are often 8-20 residues apart, so a tag of min_tag_size masses "
                     "needs several kDa.");
  defaults_.setMinFloat("max_gap_mass", 100.0);
  defaults_.setValue("min_tag_size", 5,
                     "Smallest number of masses in a gapped tag; every tag of at least this many masses is found, with all the masses that "
                     "land with it. Gaps of many residues are sums of residues almost always, so a short tag lands on a decoy about as often "
                     "as on the right protein (with 3, nearly always; with 4, often in crowded spectra); chance landings also take the "
                     "candidate places of other sources.");
  defaults_.setMinInt("min_tag_size", 3);
  defaults_.setMaxInt("min_tag_size", 100);
  defaults_.setValue("top_n", 400, "Number of masses per spectrum (highest Qscore) used for gapped tags.");
  defaults_.setMinInt("top_n", 3);
  defaults_.setValue("max_proteins", 100,
                     "Number of proteins per spectrum whose gapped-tag landings are reported. Every tag is found, and the proteins they land "
                     "on are ranked by the score FLASHTnT ranks candidates by (the masses of the spectrum against the protein at the "
                     "landings' shifts); in a crowded spectrum, tags of min_tag_size masses land on most proteins of a database by chance.");
  defaults_.setMinInt("max_proteins", 1);
  defaults_.setValue("bin_width", 0.1, "Mass bin width of the protein window index (Da). A tolerance window is looked up in every bin it touches.");
  defaults_.setMinFloat("bin_width", 0.001);
  defaults_.addTag("bin_width", "advanced");
  defaultsToParam_();
}

void FLASHGappedTaggerAlgorithm::updateMembers_()
{
  max_gap_mass_ = param_.getValue("max_gap_mass");
  min_tag_size_ = param_.getValue("min_tag_size");
  top_n_ = param_.getValue("top_n");
  max_proteins_ = param_.getValue("max_proteins");
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
  const Size n_bins = (Size)key_max_ + 1, n_proteins = offsets_.size() - 1;
  int threads = 1;
#ifdef _OPENMP
  if (! omp_in_parallel()) threads = std::max(1, omp_get_max_threads());
#endif
  // the proteins in chunks of about equal residue counts, counted and filled in parallel: a bin receives the positions of the chunks
  // in chunk order, so they stay ascending
  const int n_chunks = threads;
  std::vector<Size> chunk_begin(n_chunks + 1, n_proteins);
  for (int c = 0; c < n_chunks; ++c)
  {
    const auto first = std::lower_bound(offsets_.begin(), offsets_.begin() + (std::ptrdiff_t)n_proteins, total * c / n_chunks);
    chunk_begin[c] = (Size)(first - offsets_.begin());
  }
  // per chunk and bin the number of windows, then where the next window goes in the posting list of the bin
  std::vector<std::vector<std::uint64_t>> counts(n_chunks, std::vector<std::uint64_t>(n_bins, 0));
  // first pass counts the windows per bin, the second fills the posting lists
  auto enumerate = [&](int chunk, bool fill) {
    auto& cursor = counts[chunk];
    for (Size p = chunk_begin[chunk]; p < chunk_begin[chunk + 1]; ++p)
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
          if (fill) { postings_[cursor[(Size)key]++] = (std::uint32_t)start; }
          else { ++cursor[(Size)key]; }
        }
      }
    }
  };
#pragma omp parallel for schedule(dynamic, 1) if (threads > 1)
  for (int c = 0; c < n_chunks; ++c) enumerate(c, false);
  slot_.assign(n_bins, std::numeric_limits<std::uint32_t>::max());
  post_offsets_.assign(1, 0);
  for (Size key = 0; key < n_bins; ++key)
  {
    std::uint64_t cursor = post_offsets_.back();
    for (auto& chunk_counts : counts)
    {
      const std::uint64_t count = chunk_counts[key];
      chunk_counts[key] = cursor;
      cursor += count;
    }
    if (cursor == post_offsets_.back()) continue;
    slot_[key] = (std::uint32_t)(post_offsets_.size() - 1);
    post_offsets_.push_back(cursor);
  }
  postings_.assign(post_offsets_.back(), 0);
#pragma omp parallel for schedule(dynamic, 1) if (threads > 1)
  for (int c = 0; c < n_chunks; ++c) enumerate(c, true);

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
  // the words in blocks, in parallel (bits ORed, in any order)
  const Size block_words = 4096;
  const int n_blocks = (int)((words + block_words - 1) / block_words);
  for (int len = 1; len <= max_length; ++len)
  {
    bool any = false;
#pragma omp parallel for schedule(static) reduction(|| : any) if (threads > 1)
    for (int b = 0; b < n_blocks; ++b)
    {
      const Size begin = (Size)b * block_words, end = std::min(words, begin + block_words);
      std::fill(next.begin() + (std::ptrdiff_t)begin, next.begin() + (std::ptrdiff_t)end, 0);
      for (auto u : units) orShifted(next, current, u, begin, end);
      if (tail && end == words) next.back() &= (1ULL << tail) - 1;
      for (Size i = begin; i < end; ++i)
      {
        plausible_bits_[i] |= next[i];
        any = any || next[i] != 0;
      }
    }
    current.swap(next);
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
  const int need = min_tag_size_ - 1; // masses of a tag besides its base
  const int max_length = (int)(max_gap_mass_ / min_residue_mass);

  // the offsets of a base to the other masses within max_gap_mass, ascending, as the index keys of their tolerance windows, if they are
  // sums of residues (no window has any other mass). The base is the lowest mass of a prefix tag, the highest of a suffix tag.
  struct Offset
  {
    int node;
    std::int64_t key_lo, key_hi;
  };
  std::vector<Offset> offs;
  auto offsetsOf = [&](int b, bool prefix) {
    offs.clear();
    const int step = prefix ? 1 : -1;
    for (int j = b + step; j >= 0 && j < M; j += step)
    {
      const double off = prefix ? nodes[j].mass - nodes[b].mass : nodes[b].mass - nodes[j].mass;
      const double d = tol * (nodes[j].mass + nodes[b].mass);
      if (off - d > max_gap_mass_) break;
      if (off + d < min_residue_mass || ! plausible_(off - d - plausible_margin, off + d + plausible_margin)) continue;
      offs.push_back({j, std::max(key_(off - d), key_min_), std::min(key_(off + d), key_max_)});
    }
  };
  // the masses whose offsets the windows starting at p match, one per window (the highest Qscore: two offsets in the range of one
  // window would otherwise both count), and the window length of each
  std::vector<int> aligned, lengths;
  auto align = [&](std::uint32_t p, Size protein_end) {
    aligned.clear();
    lengths.clear();
    int len = 1;
    for (const auto& o : offs)
    {
      // the offsets ascend, so the window matching the next one is not shorter
      std::int64_t key = 0;
      for (; len <= max_length && p + len <= protein_end && valid_[p + len - 1]; ++len)
      {
        key = key_(cum_mass_[p + len] - cum_mass_[p]);
        if (key >= o.key_lo) break;
      }
      if (len > max_length || p + len > protein_end || ! valid_[p + len - 1]) break;
      if (key > o.key_hi) continue;
      if (! lengths.empty() && lengths.back() == len)
      {
        if (nodes[o.node].qscore > nodes[aligned.back()].qscore) aligned.back() = o.node;
        continue;
      }
      aligned.push_back(o.node);
      lengths.push_back(len);
    }
  };
  auto proteinOf = [&](std::uint32_t p) {
    return (Size)(std::upper_bound(offsets_.begin(), offsets_.end(), (Size)p) - offsets_.begin()) - 1;
  };

  // Every tag: each offset of a base votes for the window starts with a window in its key range (a range narrower than a residue holds
  // a start at most once, as the windows of a start differ by a residue). The votes are counted for a block of starts (positions in the
  // concatenated proteins) at a time, small enough for the first-level cache: the posting list of each bin ascends, so all lists of a
  // base are consumed block by block, each from where it stopped. Per start in the block, the offsets that reach it (saturating; zero
  // again after the block).
  constexpr std::uint32_t block_size = 1u << 15;
  const std::uint32_t n_positions = (std::uint32_t)(cum_mass_.size() - 1);
  std::vector<std::uint8_t> votes(block_size, 0);
  std::vector<std::uint32_t> voted(block_size + 1), reached; // + 1: the branchless append writes before it counts
  struct Stream
  {
    const std::uint32_t* q;
    const std::uint32_t* end;
  };
  std::vector<Stream> streams;
  // a landing: window start, base, its lowest mass and series; its other masses are recovered for the reported landings only
  struct Found
  {
    std::uint32_t p;
    int base, low;
    bool prefix;
    Size protein;
    double shift; ///< as FLASHTnT turns a landing into a shift: residue mass before (prefix) or after (suffix) the tag minus its lowest mass
  };
  std::vector<Found> found;
  for (int series = 0; series < 2; ++series)
  {
    const bool prefix = series == 0;
    for (int b = 0; b < M; ++b)
    {
      offsetsOf(b, prefix);
      if ((int)offs.size() < need) continue;
      streams.clear();
      for (const auto& o : offs)
      {
        for (std::int64_t key = o.key_lo; key <= o.key_hi; ++key)
        {
          const auto s = slot_[(Size)key];
          if (s == std::numeric_limits<std::uint32_t>::max()) continue;
          streams.push_back({postings_.data() + post_offsets_[s], postings_.data() + post_offsets_[s + 1]});
        }
      }
      reached.clear();
      for (std::uint32_t block = 0; block < n_positions && ! streams.empty(); block += block_size)
      {
        const std::uint32_t block_end = std::min(n_positions, block + block_size);
        Size n_voted = 0;
        for (Size i = 0; i < streams.size();)
        {
          auto& st = streams[i];
          for (; st.q != st.end && *st.q < block_end; ++st.q)
          {
            const std::uint32_t at = *st.q - block;
            std::uint8_t& v = votes[at];
            voted[n_voted] = at; // kept if this is the first vote of the start (without a branch)
            n_voted += v == 0;
            v += v != std::numeric_limits<std::uint8_t>::max();
          }
          if (st.q == st.end)
          {
            st = streams.back(); // exhausted
            streams.pop_back();
          }
          else ++i;
        }
        for (Size i = 0; i < n_voted; ++i)
        {
          if (votes[voted[i]] >= need) reached.push_back(block + voted[i]);
          votes[voted[i]] = 0;
        }
      }
      for (std::uint32_t p : reached)
      {
        const Size protein = proteinOf(p);
        const Size protein_end = offsets_[protein + 1];
        align(p, protein_end);
        if ((int)aligned.size() < need) continue;
        const int low = prefix ? b : aligned.back();
        const double flank_mass = prefix ? cum_mass_[p] - cum_mass_[offsets_[protein]] : cum_mass_[protein_end] - cum_mass_[p + lengths.back()];
        found.push_back({p, b, low, prefix, protein, flank_mass - nodes[low].mass});
      }
    }
  }

  // Per protein, the score FLASHTnT ranks candidates by: the masses of the spectrum against the residue sums of a side at the shifts of
  // its landings (FLASHTaggerAlgorithm::scoreShifts; shifts within the tolerance count once), the better side. Only the landings of the
  // max_proteins best proteins are reported.
  std::vector<double> score_masses {0.0}; // mass 0 stands for the terminus, as in FLASHTnT
  score_masses.reserve(dspec.size() + 1);
  for (const auto& pg : dspec) score_masses.push_back(pg.getMonoMass());
  std::sort(score_masses.begin() + 1, score_masses.end());
  std::sort(found.begin(), found.end(), [](const Found& a, const Found& c) {
    return a.protein != c.protein ? a.protein < c.protein : (a.prefix != c.prefix ? a.prefix > c.prefix : a.shift < c.shift);
  });
  std::vector<std::pair<double, Size>> protein_scores; // score, first landing of the protein in found
  std::vector<double> residue_sums, shifts, flanks;
  for (Size f0 = 0; f0 < found.size();)
  {
    Size f1 = f0;
    while (f1 < found.size() && found[f1].protein == found[f0].protein) ++f1;
    const Size begin = offsets_[found[f0].protein], end = offsets_[found[f0].protein + 1];
    double score = 0;
    for (Size s0 = f0; s0 < f1;)
    {
      Size s1 = s0;
      while (s1 < f1 && found[s1].prefix == found[s0].prefix) ++s1;
      const bool prefix = found[s0].prefix;
      shifts.clear();
      flanks.clear();
      for (Size k = s0; k < s1; ++k)
      {
        const double flank = nodes[found[k].low].mass;
        if (! shifts.empty() && found[k].shift - shifts.back() <= tol * (flank + flanks.back())) continue;
        shifts.push_back(found[k].shift);
        flanks.push_back(flank);
      }
      residue_sums.clear(); // ascending, 0 included
      for (Size q = begin; q <= end; ++q)
      {
        residue_sums.push_back(prefix ? cum_mass_[q] - cum_mass_[begin] : cum_mass_[end] - cum_mass_[end - (q - begin)]);
      }
      score = std::max(score, FLASHTaggerAlgorithm::scoreShifts(score_masses, residue_sums, shifts.data(), flanks.data(), shifts.size(), ppm));
      s0 = s1;
    }
    protein_scores.emplace_back(score, f0);
    f0 = f1;
  }
  if ((int)protein_scores.size() > max_proteins_)
  {
    std::nth_element(protein_scores.begin(), protein_scores.begin() + max_proteins_, protein_scores.end(),
                     [](const auto& a, const auto& c) { return a.first != c.first ? a.first > c.first : a.second < c.second; });
    protein_scores.resize(max_proteins_);
  }

  // the reported landings with their masses: per window start and series, the tag with the most masses, then the highest Qscore sum
  std::unordered_map<std::uint64_t, Landing> best_landing;
  for (const auto& [score, f0] : protein_scores)
  {
    for (Size k = f0; k < found.size() && found[k].protein == found[f0].protein; ++k)
    {
      const auto& f = found[k];
      const Size protein_end = offsets_[f.protein + 1];
      offsetsOf(f.base, f.prefix);
      align(f.p, protein_end);
      double qsum = nodes[f.base].qscore;
      for (int j : aligned) qsum += nodes[j].qscore;
      const int n = (int)aligned.size() + 1;
      const std::uint64_t id = ((std::uint64_t)f.p << 1) | (f.prefix ? 1 : 0);
      auto it = best_landing.find(id);
      if (it != best_landing.end()
          && ((int)it->second.masses.size() > n || ((int)it->second.masses.size() == n && it->second.qscore_sum >= qsum)))
        continue;
      Landing l;
      l.protein = f.protein;
      l.start = (int)(f.p - offsets_[f.protein]);
      l.end = l.start + lengths.back();
      l.prefix = f.prefix;
      l.flank_mass = f.prefix ? cum_mass_[f.p] - cum_mass_[offsets_[f.protein]] : cum_mass_[protein_end] - cum_mass_[f.p + lengths.back()];
      // ascending: the base first (prefix) or last (suffix, whose offsets ascend as the masses descend)
      l.masses.reserve(n);
      l.scores.reserve(n);
      if (f.prefix)
      {
        l.masses.push_back(nodes[f.base].mass);
        l.scores.push_back(nodes[f.base].score);
      }
      for (int i = 0; i < (int)aligned.size(); ++i)
      {
        const int j = aligned[f.prefix ? i : (int)aligned.size() - 1 - i];
        l.masses.push_back(nodes[j].mass);
        l.scores.push_back(nodes[j].score);
      }
      if (! f.prefix)
      {
        l.masses.push_back(nodes[f.base].mass);
        l.scores.push_back(nodes[f.base].score);
      }
      l.qscore_sum = qsum;
      if (it != best_landing.end()) { it->second = std::move(l); }
      else { best_landing.emplace(id, std::move(l)); }
    }
  }

  landings.reserve(best_landing.size());
  for (auto& [id, l] : best_landing) landings.push_back(std::move(l));
  std::sort(landings.begin(), landings.end(), [](const Landing& a, const Landing& b) {
    return a.protein != b.protein ? a.protein < b.protein : (a.start != b.start ? a.start < b.start : a.prefix > b.prefix);
  });
}
} // namespace OpenMS
