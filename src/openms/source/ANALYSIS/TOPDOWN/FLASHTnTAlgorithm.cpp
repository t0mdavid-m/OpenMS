// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHGappedTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <boost/math/special_functions/gamma.hpp>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>
#include <unordered_map>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
// One matched tag is enough to propose a protein: a correct protein often has a single tag in a spectrum,
// and the protein-level score and the extension decide between the candidates.
inline const int min_tag_count = 1;

namespace
{
// FDR classes of the PrSMs (meta value FDRClass): determined mass and termini by number of blind modifications (classes 0 to
// fdr_class_max_blind, the last one with that many or more), and undetermined by whether both termini are known. The "; N blind" form is
// kept for parsers of the label. The PrSMs of precursor-ladder candidates (all determined) form a class of their own by their source,
// so that their decoys are not hidden in the large, nearly decoy-free determined classes.
constexpr int fdr_class_max_blind = 2;
constexpr int fdr_class_both_termini = fdr_class_max_blind + 1;
constexpr int fdr_class_open_terminus = fdr_class_both_termini + 1;
constexpr int fdr_class_precursor_ladder = fdr_class_open_terminus + 1;
constexpr int fdr_class_count = fdr_class_precursor_ladder + 1;
const char* const fdr_class_names[] = {"determined; 0 blind", "determined; 1 blind", "determined; 2 blind", "undetermined, both termini",
                                       "undetermined, open terminus", "precursor ladder"};
static_assert(std::size(fdr_class_names) == (Size)fdr_class_count, "one label per FDR class");
// protein level: the null of a class is exponential above this quantile of its decoy scores
constexpr double null_tail_quantile = 0.8;
// protein level: with fewer decoys the null is the empirical survival function only, without the tail
constexpr Size null_min_tail_decoys = 5;
// protein level: a class with fewer decoys uses the null of the decoys of all classes pooled
constexpr Size null_min_class_decoys = 20;
// protein level: floor of a null p-value and of a local FDR before its log
constexpr double fdr_p_floor = 1e-300;
// precursor ladders: a window needs this many matched ladder positions to be a candidate (the prototype's calibrated setting)
constexpr int precursor_ladder_min_k = 3;
// variable modification penalties: learned from a pilot extension of every var_mod_pilot_stride-th spectrum in input order (a fixed
// subset that samples the whole run, the same for every thread count), so that only that share of the spectra is extended twice
constexpr int var_mod_pilot_stride = 4;
// variable modification penalties: with fewer decoy PrSMs, chance usage is not estimated and the penalties are not learned
constexpr Size var_mod_min_decoy_prsms = 50;

// mass and both termini known (unknown: -1)
bool isDetermined(const ProteinHit& hit)
{
  return ! ((double)hit.getMetaValue("Mass") < 0 || (int)hit.getMetaValue("StartPosition") < 0 || (int)hit.getMetaValue("EndPosition") < 0);
}

// an item of a target-decoy ranking
template<typename Key>
struct RankedItem
{
  Key key;       // equal keys (exact) are tied
  bool decoy;
  double weight; // of a decoy
};

// q-values of a target-decoy ranking (best first). After each tie group, FDR = (D + 1) / (D + T + 1) with D the decoy weight and T the
// target count so far: with the +1 pseudo-count, a ranking without a decoy above an item does not give q = 0 there (Eclip_1480
// entrapment: such classes passed false proteins at several times the nominal rate). q is the running minimum of the FDR from the bottom.
template<typename Key>
std::vector<double> targetDecoyQvalues(const std::vector<RankedItem<Key>>& items)
{
  std::vector<double> qvalues(items.size(), 1.0);
  double decoy_count = 0, target_count = 0;
  for (Size i = 0; i < items.size();)
  {
    Size j = i;
    for (; j < items.size() && items[j].key == items[i].key; j++)
    {
      if (items[j].decoy) decoy_count += items[j].weight;
      else target_count++;
    }
    std::fill(qvalues.begin() + i, qvalues.begin() + j, std::min(1.0, (decoy_count + 1.0) / (decoy_count + target_count + 1.0)));
    i = j;
  }
  double cummin = 1.0;
  for (Size k = items.size(); k-- > 0;)
  {
    cummin = std::min(cummin, qvalues[k]);
    qvalues[k] = cummin;
  }
  return qvalues;
}

// Null distribution of the PrSM score from a sample of decoy scores: the empirical survival function with +1 and, from
// null_min_tail_decoys on, an exponential tail above the null_tail_quantile quantile (mean excess as its scale), which goes on decreasing
// beyond the best decoy.
class DecoyScoreNull
{
public:
  explicit DecoyScoreNull(std::vector<double> scores): scores_(std::move(scores))
  {
    std::sort(scores_.begin(), scores_.end());
    const Size n = scores_.size();
    if (n < null_min_tail_decoys) return;
    with_tail_ = true;
    threshold_ = scores_[(Size)(null_tail_quantile * (double)(n - 1))];
    // mean excess over the threshold; Neumaier-compensated sum, bit-identical to the reference implementation of this FDR
    double excess = 0, compensation = 0;
    Size excess_count = 0;
    for (double y : scores_)
    {
      if (y <= threshold_) continue;
      const double e = y - threshold_, t = excess + e;
      compensation += std::abs(excess) >= std::abs(e) ? (excess - t) + e : (e - t) + excess;
      excess = t;
      excess_count++;
    }
    scale_ = excess_count > 0 ? (excess + compensation) / (double)excess_count : 1.0;
    survival_at_threshold_ = (countAtLeast_(threshold_) + 1.0) / ((double)n + 1.0);
  }

  // the chance of a decoy score >= x, floored at fdr_p_floor
  double pValue(double x) const
  {
    const double survival = with_tail_ && x > threshold_ ? survival_at_threshold_ * std::exp(-(x - threshold_) / scale_)
                                                         : (countAtLeast_(x) + 1.0) / ((double)scores_.size() + 1.0);
    return std::max(survival, fdr_p_floor);
  }

private:
  // number of decoy scores >= x
  double countAtLeast_(double x) const
  {
    return (double)(scores_.end() - std::lower_bound(scores_.begin(), scores_.end(), x));
  }

  std::vector<double> scores_; // ascending
  bool with_tail_ = false;
  double threshold_ = 0, scale_ = 1.0, survival_at_threshold_ = 1.0;
};

// Density of the targets' null p-values: the slope of the least concave majorant of their ECDF through (0,0) and (1,1) (Grenander
// estimate), non-increasing in p. A segment is left-open and right-closed: a p at a vertex takes the slope to its left.
class TargetPvalueDensity
{
public:
  explicit TargetPvalueDensity(std::vector<double> pvalues)
  {
    std::sort(pvalues.begin(), pvalues.end());
    const double n = (double)pvalues.size();
    std::vector<std::pair<double, double>> hull;
    auto add = [&hull](double x, double y) {
      while (hull.size() >= 2)
      {
        const auto [x1, y1] = hull[hull.size() - 2];
        const auto [x2, y2] = hull.back();
        if ((x2 - x1) * (y - y1) - (y2 - y1) * (x - x1) >= 0) hull.pop_back(); // also drops collinear points
        else break;
      }
      hull.emplace_back(x, y);
    };
    add(0.0, 0.0);
    for (Size k = 0; k < pvalues.size(); k++) add(pvalues[k], ((double)k + 1.0) / n); // ties stay separate points
    add(1.0, 1.0);
    for (Size j = 0; j < hull.size(); j++)
    {
      xs_.push_back(hull[j].first);
      if (j + 1 == hull.size()) break;
      const double dx = hull[j + 1].first - hull[j].first;
      slopes_.push_back(dx > 0 ? (hull[j + 1].second - hull[j].second) / dx : std::numeric_limits<double>::infinity());
    }
  }

  double density(double p) const
  {
    // the hull never drops (0,0) and ends with (1,1), so there is at least one segment
    const SignedSize j = (SignedSize)(std::lower_bound(xs_.begin(), xs_.end(), p) - xs_.begin()) - 1;
    return slopes_[(Size)std::clamp(j, (SignedSize)0, (SignedSize)slopes_.size() - 1)];
  }

private:
  std::vector<double> xs_;     // hull vertices
  std::vector<double> slopes_; // of the segments between them
};

// BEGIN precursor ladder kernel
// Precursor-anchored ladders (option precursor_ladders): the windows (i, j) of the database proteins whose unmodified mass equals the
// precursor mass, scored by the deconvolved masses on their prefix ladder (read from i) and their suffix ladder (read from j). This block
// uses the standard library and boost::math::gamma_p only, so that it can be compiled on its own and checked against its prototype
// (kernel K4 of final_all/sliding/src/slide.cpp, calibrated by entrapment on Eclip_1480).
namespace precursor_ladder
{
// the residue sums a..b of a protein (residues a..b-1, all standard): a window never crosses a non-standard residue
struct Segment
{
  std::uint32_t protein, a, b;
};

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
// V = S(j) - S(i); with several types per side this bounds the chance rate of a position on a side from above.
inline void search(const std::vector<TargetSet>& sets, const std::vector<Segment>& segments, const std::vector<const double*>& sums,
                   double precursor_mass, double precursor_tol_ppm, double h2o, double isotope_mass, int min_k, std::size_t top_n, Result& result)
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

FLASHTnTAlgorithm::FLASHTnTAlgorithm(): DefaultParamHandler("FLASHTnTAlgorithm"), ProgressLogger()
{
  setDefaultParams_();
}

FLASHTnTAlgorithm::FLASHTnTAlgorithm(const FLASHTnTAlgorithm& other): DefaultParamHandler(other), ProgressLogger(other)
{
}

FLASHTnTAlgorithm& FLASHTnTAlgorithm::operator=(const FLASHTnTAlgorithm& rhs)
{
  if (this == &rhs) return *this;

  DefaultParamHandler::operator=(rhs);
  return *this;
}

void FLASHTnTAlgorithm::setDefaultParams_()
{
  defaults_.setValue("pro_fdr", 1.0,
                     "Specifies the proteoform-level FDR; proteins are reported at the same level. The reported PrSMs are those of the reported "
                     "proteoforms (PrSM q-values only rank the PrSMs within their FDR class; they are no filter). Also the PrSM q-value up to "
                     "which a determined PrSM is preferred as proteoform representative and counts as passing in the protein ranking (0.01 "
                     "when 1).");
  defaults_.setMinFloat("pro_fdr", 0.0);

  defaults_.setValue("only_single_hit", "false", "Allows only a single hit per spectrum.");
  defaults_.setValidStrings("only_single_hit", {"true", "false"});

  defaults_.setValue("discard_underdetermined", "false",
                     "Discards underdetermined proteoform IDs (e.g., those without exact precursor masses or start/end positions).");
  defaults_.setValidStrings("discard_underdetermined", {"true", "false"});

  defaults_.setValue("keep_decoy", "false", "Retains decoy hits in the results.");
  defaults_.setValidStrings("keep_decoy", {"true", "false"});

  defaults_.setValue("protein_level_only", "false",
                     "Stops after the protein-level search, skipping proteoform characterization. The candidate proteins "
                     "and their scores are still reported, so this can be used to shortlist proteins cheaply.");
  defaults_.setValidStrings("protein_level_only", {"true", "false"});

  defaults_.setValue("max_protein_hits", 50,
                     "Number of candidate proteins kept per spectrum by the protein-level search. When proteoform search "
                     "runs, only these are characterized.");
  defaults_.setMinInt("max_protein_hits", 1);

  defaults_.setValue("protein_tol", 0.0,
                     "Mass tolerance (ppm) of the protein-level search: a deconvolved mass, moved by a tag-derived shift, matches a "
                     "protein residue sum within it. 0: the MS2 deconvolution tolerance recorded in the input file.");
  defaults_.setMinFloat("protein_tol", 0.0);

  defaults_.setValue("min_protein_hits", 5,
                     "Number of best-scoring candidate proteins per spectrum that are always characterized. The other "
                     "candidates are characterized only if their protein, or its target-decoy partner, passes protein_fdr. Used only when the database "
                     "contains decoys (DECOY prefix); otherwise all candidates kept by max_protein_hits are characterized.");
  defaults_.setMinInt("min_protein_hits", 1);

  defaults_.setValue("protein_fdr", 0.1,
                     "Run-level protein FDR, from each protein's best tag-based score over all spectra, that a candidate "
                     "beyond min_protein_hits (or its target-decoy partner) must pass to be characterized. Needs decoys in the "
                     "database. 1 characterizes all candidates kept by max_protein_hits.");
  defaults_.setMinFloat("protein_fdr", 0.0);
  defaults_.setMaxFloat("protein_fdr", 1.0);

  defaults_.setValue("ion_type", std::vector<std::string> {"b", "y"}, "Specifies ion types to consider.");
  defaults_.setValidStrings("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

  auto tparam = FLASHTaggerAlgorithm().getDefaults();
  tparam.remove("ion_type");
  defaults_.insert("tag:", tparam);
  auto eparam = FLASHExtenderAlgorithm().getDefaults();
  eparam.remove("ion_type");

  defaults_.insert("ex:", eparam);

  defaults_.setValue("gapped_tags", "fallback",
                     "Candidate proteins from gapped tags (masses whose differences are sums of one or more residues): off; fallback (spectra in "
                     "which no sequence tag proposes a protein); both (every spectrum, proteins not proposed by sequence tags); only (gapped "
                     "tags instead of sequence tags). Where a gapped tag lands on a protein anchors the extension as a sequence tag does.");
  defaults_.setValidStrings("gapped_tags", {"off", "fallback", "both", "only"});
  defaults_.setValue("terminal_ladders", "true",
                     "Candidate proteins from terminal ladders: deconvolved masses that are residue sums from an intact protein terminus (the N-terminus "
                     "also without the initiator Met and with acetylation), found without a sequence tag. In every spectrum the min_protein_hits best "
                     "proteins by the chance probability of their ladder are characterized, anchored on the ladder.");
  defaults_.setValidStrings("terminal_ladders", {"true", "false"});
  defaults_.setValue("precursor_ladders", "true",
                     "Candidate proteins from precursor-anchored ladders: in each MSn spectrum with a precursor mass, every stretch of a protein (between "
                     "non-standard residues) whose unmodified mass equals the precursor mass (within the MS1 deconvolution tolerance of the input, or "
                     "protein_tol's tolerance if the input has no MS1 spectra; isotope errors -1 to +1) is scored by the chance probability of the "
                     "deconvolved masses on its prefix and suffix fragment ladders (a position counts once per side, whatever ion types match it). The "
                     "stretches are unmodified: N-terminal acetylation and other modifications are not proposed (initiator Met removal is a stretch "
                     "from residue 2). The min_protein_hits best stretches are characterized from tags on their ladders. Only their PrSMs of "
                     "determined mass and termini that are not the same proteoform as a determined PrSM of another source in the spectrum are kept; "
                     "they form their own FDR class (\"precursor ladder\"). Not used with protein_level_only. Without any precursor mass in the input, "
                     "results are unchanged; otherwise the PrSMs of this source can also change the others' results (the pooled protein-level null, "
                     "the best PrSM of a spectrum, proteoform clusters).");
  defaults_.setValidStrings("precursor_ladders", {"true", "false"});
  defaults_.setValue("decoy_precursor", "false",
                     "Decoy precursors: spectra whose PrSM took its mass from MS1 are searched again with the decoy precursor mass of their isolation "
                     "window instead of the precursor mass (FLASHDeconv -FD:report_FDR true writes it); a PrSM reporting that mass competes with the "
                     "spectrum's own and counts as a decoy in the FDR. Controls wrong MS1 precursor masses; PrSMs with masses from fragment pairs "
                     "or without a mass are not affected.");
  defaults_.setValidStrings("decoy_precursor", {"true", "false"});
  defaults_.setValue("learn_var_mod_penalty", "true",
                     "Raises the penalty of a variable modification that target PrSMs use per opportunity (fitting residue or terminus) no more often than decoy PrSMs "
                     "(chance), up to that of a blind shift. The usage is counted in a pilot pass over every 4th spectrum (all spectra if the pilot "
                     "has fewer than 50 decoy PrSMs); with a raised penalty, all spectra are then extended once with the learned penalties. "
                     "Needs decoys and at least 50 decoy PrSMs.");
  defaults_.setValidStrings("learn_var_mod_penalty", {"true", "false"});
  defaults_.addTag("learn_var_mod_penalty", "advanced");
  defaults_.insert("gapped:", FLASHGappedTaggerAlgorithm().getDefaults());
  defaults_.setValue("raw_score_weight", 8.0,
                     "With the raw mzML (-in_raw), every PrSM is annotated on its raw spectrum (FLASHAnnotator: isotope envelopes of the proteoform's "
                     "fragments at every charge) and ranked for the FDR by its extension score plus this weight times its raw score (-log10 of the "
                     "chance of the matched fragments). 0: the raw score is reported but not used in the ranking.");
  defaults_.setMinFloat("raw_score_weight", 0.0);
  auto aparam = FLASHAnnotatorAlgorithm().getDefaults();
  aparam.remove("ion_type");
  defaults_.insert("raw:", aparam);
  defaultsToParam_();
}

void FLASHTnTAlgorithm::updateMembers_()
{
  tagger_param_ = param_.copy("tag:", true);
  tagger_param_.setValue("ion_type", param_.getValue("ion_type"));
  extender_param_ = param_.copy("ex:", true);
  extender_param_.setValue("ion_type", param_.getValue("ion_type"));
  gapped_param_ = param_.copy("gapped:", true);
  annotator_param_ = param_.copy("raw:", true);
  annotator_param_.setValue("ion_type", param_.getValue("ion_type"));
  raw_score_weight_ = param_.getValue("raw_score_weight");
  gapped_mode_ = param_.getValue("gapped_tags").toString();
  terminal_ladders_ = param_.getValue("terminal_ladders").toString() == "true";
  precursor_ladders_ = param_.getValue("precursor_ladders").toString() == "true";
  learn_var_mod_penalty_ = param_.getValue("learn_var_mod_penalty").toString() == "true";
  decoy_precursor_ = param_.getValue("decoy_precursor").toString() == "true";
  pro_fdr_ = param_.getValue("pro_fdr");
  keep_decoy_ = param_.getValue("keep_decoy").toString() == "true";
  keep_underdetermined_ = param_.getValue("discard_underdetermined").toString() == "false";
  multiple_hits_per_spec_ = param_.getValue("only_single_hit").toString() == "false";
  protein_level_only_ = param_.getValue("protein_level_only").toString() == "true";
  max_protein_hits_ = param_.getValue("max_protein_hits");
  protein_tol_ = param_.getValue("protein_tol");
  min_protein_hits_ = param_.getValue("min_protein_hits");
  protein_fdr_ = param_.getValue("protein_fdr");
}

std::set<Size> FLASHTnTAlgorithm::getProteinsPassingFDR_(const std::vector<FASTAFile::FASTAEntry>& fasta_entry) const
{
  std::map<Size, double> best_score; // fasta index to the protein's best tag-based score over the run
  for (const auto& hit : protein_hits_)
  {
    Size fi = (int)hit.getMetaValue("FastaIndex");
    auto it = best_score.find(fi);
    if (it == best_score.end() || it->second < hit.getScore()) best_score[fi] = hit.getScore();
  }

  std::vector<std::pair<double, Size>> ranked;
  ranked.reserve(best_score.size());
  for (const auto& [fi, score] : best_score)
  {
    ranked.emplace_back(score, fi);
  }
  std::sort(ranked.begin(), ranked.end(), [](const auto& a, const auto& b) { return a.first > b.first; });

  // target-decoy FDR at each distinct score (ties together), then the running minimum from the bottom
  std::vector<double> fdr(ranked.size(), 1.0);
  double taget_count = 0, decoy_count = 0;
  for (Size i = 0; i < ranked.size();)
  {
    Size j = i;
    for (; j < ranked.size() && ranked[j].first == ranked[i].first; j++)
    {
      if (fasta_entry[ranked[j].second].identifier.hasPrefix("DECOY")) { decoy_count++; }
      else { taget_count++; }
    }
    double f = taget_count > 0 ? decoy_count / decoy_factor_ / taget_count : 1.0;
    std::fill(fdr.begin() + i, fdr.begin() + j, f);
    i = j;
  }

  std::set<Size> passing; // decoys included, so that the decoy hits extended stay comparable to the targets
  double qvalue = 1.0;
  for (Size k = ranked.size(); k-- > 0;)
  {
    qvalue = std::min(qvalue, fdr[k]);
    if (qvalue <= protein_fdr_) passing.insert(ranked[k].second);
  }

  // Add the target-decoy partner of every passing protein. A protein present in the sample passes through its true matches in
  // other spectra, which a decoy never has; extended in a spectrum that is not its own, it can give a false match that no decoy
  // mirrors (Eclip_1480: beyond the top candidates 237 target vs 0 decoy matches). With its decoy partner extended under the same
  // condition, that kind of false match is counted by the target-decoy competition too.
  std::unordered_map<std::string, Size> index_of;
  index_of.reserve(fasta_entry.size());
  for (Size i = 0; i < fasta_entry.size(); i++)
  {
    index_of.emplace(fasta_entry[i].identifier, i);
  }
  std::set<Size> partners;
  for (Size i : passing)
  {
    const String& id = fasta_entry[i].identifier;
    std::vector<std::string> names;
    if (id.hasPrefix("DECOY_")) { names.push_back(id.substr(6)); }
    else if (id.hasPrefix("DECOY")) { names.push_back(id.substr(5)); }
    else { names = {"DECOY_" + id, "DECOY" + id}; }
    for (const auto& name : names)
    {
      auto it = index_of.find(name);
      if (it == index_of.end()) continue;
      partners.insert(it->second);
      break;
    }
  }
  passing.insert(partners.begin(), partners.end());
  return passing;
}

Size FLASHTnTAlgorithm::blindModCount_(const ProteinHit& hit)
{
  if (hit.metaValueExists("BlindModCount")) return (Size)(int)hit.getMetaValue("BlindModCount");
  return hit.getMetaValue("Modifications").toDoubleList().size();
}

bool FLASHTnTAlgorithm::areConsistent_(const ProteinHit& a, const ProteinHit& b, double tol)
{
  // Two PrSMs of one protein are counted as one proteoform when
  // - both masses are determined: the masses agree (within twice the tolerance), or differ by one isotope (1.002371 Da, within twice
  //   the tolerance), and the ranges overlap. The monoisotopic mass of a precursor is off by one isotope often enough that the same
  //   proteoform was counted twice: on Eclip_1480, 1450 determined proteoforms became 1191 with the one-isotope rule (proteins 299
  //   either way). A deamidation (+0.984 Da) cannot be told from a one-isotope error at this tolerance and is merged with it. Two
  //   isotopes apart is not merged: a reduced disulfide (+2.016 Da) is another proteoform, and the two cannot be told apart here.
  //   Where the shifts sit does not matter: a blind shift has no exact position, and the same mass placed elsewhere (+1.0 Da at
  //   residue 26, 29-30, 31 or 33 of RL32 on Eclip_1480) is no evidence of another proteoform. Requiring the same positions and
  //   names counted one proteoform up to four times (one in five proteoform rows of the lysate runs were such repeats);
  // - otherwise (an undetermined mass): every terminus known in both is the same. The modifications of an undetermined PrSM are
  //   not confirmed by a precursor mass, and an open terminus may be any position, so such a PrSM is a separate proteoform only if
  //   no other PrSM of the protein fits its known termini. Counting it apart whenever its blind shifts or its open end differed
  //   counted a determined proteoform again (232 of 511 undetermined proteoform rows of Eclip_1480 shared their known termini with a
  //   determined proteoform).
  if (a == b) return true;
  const double mass1 = a.getMetaValue("Mass"), mass2 = b.getMetaValue("Mass");
  const int sp1 = a.getMetaValue("StartPosition"), ep1 = a.getMetaValue("EndPosition");
  const int sp2 = b.getMetaValue("StartPosition"), ep2 = b.getMetaValue("EndPosition");
  if (mass1 > 0 && mass2 > 0)
  {
    const bool overlap = (sp1 < 0 || ep2 < 0 || sp1 <= ep2) && (sp2 < 0 || ep1 < 0 || sp2 <= ep1);
    const double diff = std::abs(mass1 - mass2), max_diff = std::max(mass1, mass2) * tol * 2e-6;
    return overlap && (diff <= max_diff || std::abs(diff - Constants::ISOTOPE_MASSDIFF_55K_U) <= max_diff);
  }
  return (sp1 < 0 || sp2 < 0 || sp1 == sp2) && (ep1 < 0 || ep2 < 0 || ep1 == ep2);
}

std::vector<std::vector<Size>> FLASHTnTAlgorithm::markRepresentativeProteoformHits_(double tol, const std::vector<int>& hit_class,
                                                                                   const std::vector<double>& qvalues, double rep_q, bool by_qvalue)
{
  std::map<String, std::vector<int>> proteoform_map;
  const int hit_count = (int)proteoform_hits_.size();
  std::vector<bool> passing_determined(hit_count, false), with_mass(hit_count, false);

  for (int i = 0; i < hit_count; i++)
  {
    auto& hit = proteoform_hits_[i];
    hit.setMetaValue("Index", i);
    with_mass[i] = (double)hit.getMetaValue("Mass") > 0;
    if (by_qvalue)
    {
      // A PrSM in no FDR class (undetermined, with discard_underdetermined) has no q-value and is in no proteoform. Sorted as q = 0, it
      // became the representative of determined PrSMs, and their proteoform was lost when the output filter dropped it.
      if (hit_class[i] < 0) continue;
      passing_determined[i] = qvalues[i] <= rep_q && isDetermined(hit);
    }
    proteoform_map[hit.getAccession()].push_back(i);
  }

  std::vector<std::vector<Size>> clusters;
  for (auto& [acc, hit_idx] : proteoform_map)
  {
    // The representative of a proteoform is its row and is tested for consistency with the other PrSMs. It is a determined PrSM passing
    // rep_q if there is one (rather than an undetermined PrSM of another class with a lower q-value), else the PrSM with the lowest
    // q-value; then one with a mass, then the best scoring; then the order of proteoform_hits_. So a representative passes the PrSM
    // filter if one of its PrSMs does, and a proteoform is not lost while one of its PrSMs is reported.
    std::stable_sort(hit_idx.begin(), hit_idx.end(), [&](int l, int r) -> bool {
      if (passing_determined[l] != passing_determined[r]) return passing_determined[l];
      if (by_qvalue && qvalues[l] != qvalues[r]) return qvalues[l] < qvalues[r];
      if (with_mass[l] != with_mass[r]) return with_mass[l];
      return proteoform_hits_[l].getScore() > proteoform_hits_[r].getScore();
    });
    for (Size i = 0; i < hit_idx.size(); i++)
    {
      auto& hit = proteoform_hits_[hit_idx[i]];
      if (hit.metaValueExists("Representative")) continue;
      std::vector<int> prsm_indices {(int)hit.getMetaValue("Index")};
      for (Size j = i + 1; j < hit_idx.size(); j++)
      {
        auto& hit2 = proteoform_hits_[hit_idx[j]];
        if (hit2.metaValueExists("Representative")) continue;
        // a decoy precursor PrSM carries the target accession: merged into a target proteoform, the decoy would no longer be counted
        if (hit.metaValueExists("DecoyPrecursor") != hit2.metaValueExists("DecoyPrecursor")) continue;
        if (areConsistent_(hit, hit2, tol))
        {
          prsm_indices.push_back(hit2.getMetaValue("Index"));
          hit2.setMetaValue("Representative", 0);
        }
      }
      hit.setMetaValue("Representative", 1);
      hit.setMetaValue("PrSMIndices", prsm_indices);
      clusters.emplace_back(prsm_indices.begin(), prsm_indices.end());
    }
  }
  // in the order of the representatives in proteoform_hits_
  std::sort(clusters.begin(), clusters.end(), [](const std::vector<Size>& l, const std::vector<Size>& r) { return l[0] < r[0]; });
  return clusters;
}

void FLASHTnTAlgorithm::setProteinQvalues_(const std::vector<int>& hit_class, const std::vector<double>& qvalues, const std::vector<bool>& decoys,
                                           const std::vector<double>& decoy_weights, const std::vector<std::vector<Size>>& clusters, double rep_q)
{
  // A protein is ranked by the evidence of all its proteoforms: the sum over them of the log of the lowest local FDR among the proteoform's
  // PrSMs, then by its best score among its passing PrSMs (rep_q; among all of them if none passes); proteins are then target-decoy counted.
  // The sum is a ranking statistic, not a probability, and is not reported. The local FDR of a PrSM comes from a two-group model of its
  // class: the null p-value of its score under the decoy scores of the class (DecoyScoreNull; the decoys of all classes pooled for a class
  // with fewer than null_min_class_decoys, as the determined classes with 0 to 3 decoys each in the E. coli lysate runs), the density of
  // the targets' p-values (TargetPvalueDensity) and the share of false targets pi0 = (decoy weight + 1) / targets. A decoy precursor
  // carries the target accession: it is a decoy of the precursor, not of the protein, and ranks no protein.
  const Size hit_count = proteoform_hits_.size();

  // null of the PrSM score per class
  std::vector<std::vector<double>> class_decoy_scores(fdr_class_count);
  std::vector<double> class_decoy_weight(fdr_class_count, 0.0), pooled_decoy_scores;
  for (Size i = 0; i < hit_count; i++)
  {
    if (hit_class[i] < 0 || ! decoys[i]) continue;
    class_decoy_scores[hit_class[i]].push_back(proteoform_hits_[i].getScore());
    class_decoy_weight[hit_class[i]] += decoy_weights[i]; // a decoy precursor of weight 0 still enters the null's sample
    pooled_decoy_scores.push_back(proteoform_hits_[i].getScore());
  }
  std::vector<DecoyScoreNull> nulls; // nulls[0]: the pooled one
  std::vector<Size> class_null(fdr_class_count, 0);
  nulls.reserve(fdr_class_count + 1);
  nulls.emplace_back(std::move(pooled_decoy_scores));
  for (int c = 0; c < fdr_class_count; c++)
  {
    if (class_decoy_scores[c].size() < null_min_class_decoys) continue;
    class_null[c] = nulls.size();
    nulls.emplace_back(std::move(class_decoy_scores[c]));
  }
  std::vector<double> pvalues(hit_count, 1.0);
  std::vector<std::vector<double>> class_target_pvalues(fdr_class_count);
  for (Size i = 0; i < hit_count; i++)
  {
    if (hit_class[i] < 0) continue;
    pvalues[i] = nulls[class_null[hit_class[i]]].pValue(proteoform_hits_[i].getScore());
    if (! decoys[i]) class_target_pvalues[hit_class[i]].push_back(pvalues[i]);
  }

  // two-group local FDR of every PrSM, targets and decoys, in its class
  std::vector<double> pi0;
  std::vector<TargetPvalueDensity> densities;
  densities.reserve(fdr_class_count);
  for (int c = 0; c < fdr_class_count; c++)
  {
    pi0.push_back(std::min(1.0, (class_decoy_weight[c] + 1.0) / std::max(1.0, (double)class_target_pvalues[c].size())));
    densities.emplace_back(std::move(class_target_pvalues[c]));
  }
  std::vector<double> local_fdr(hit_count, 1.0);
  for (Size i = 0; i < hit_count; i++)
  {
    if (hit_class[i] < 0) continue;
    const double density = densities[hit_class[i]].density(pvalues[i]); // +inf gives 0
    local_fdr[i] = density <= 0 ? 1.0 : std::min(1.0, pi0[hit_class[i]] / density);
  }

  // protein key: (evidence, -best score), smaller is better
  std::map<String, std::vector<double>> evidence_terms; // accession -> log local FDR of each of its proteoforms
  for (const auto& members : clusters)
  {
    const auto& rep = proteoform_hits_[members[0]];
    if (rep.metaValueExists("DecoyPrecursor")) continue;
    double min_local_fdr = local_fdr[members[0]];
    for (Size i : members) min_local_fdr = std::min(min_local_fdr, local_fdr[i]);
    evidence_terms[rep.getAccession()].push_back(std::log(std::max(fdr_p_floor, min_local_fdr)));
  }
  std::map<String, double> best_score, best_passing_score;
  for (Size i = 0; i < hit_count; i++)
  {
    const auto& hit = proteoform_hits_[i];
    if (hit_class[i] < 0 || hit.metaValueExists("DecoyPrecursor")) continue;
    const auto it = best_score.find(hit.getAccession());
    if (it == best_score.end() || it->second < hit.getScore()) best_score[hit.getAccession()] = hit.getScore();
    if (qvalues[i] > rep_q) continue;
    const auto itp = best_passing_score.find(hit.getAccession());
    if (itp == best_passing_score.end() || itp->second < hit.getScore()) best_passing_score[hit.getAccession()] = hit.getScore();
  }
  std::vector<std::pair<std::pair<double, double>, String>> ranked;
  ranked.reserve(best_score.size());
  for (const auto& [acc, best] : best_score)
  {
    double evidence = 0;
    if (const auto ite = evidence_terms.find(acc); ite != evidence_terms.end())
    {
      auto& terms = ite->second;
      std::sort(terms.begin(), terms.end()); // summed from the most negative on, so that the sum does not depend on the proteoform order
      for (double t : terms) evidence += t;
    }
    const auto it = best_passing_score.find(acc);
    ranked.emplace_back(std::make_pair(evidence, -(it == best_passing_score.end() ? best : it->second)), acc);
  }
  std::sort(ranked.begin(), ranked.end());

  // protein q-values; an accession with decoy precursor PrSMs only is not ranked and gets 1
  std::vector<RankedItem<std::pair<double, double>>> items;
  items.reserve(ranked.size());
  const double protein_decoy_weight = 1.0 / std::max(decoy_factor_, 1e-9);
  for (const auto& [key, acc] : ranked) items.push_back({key, acc.hasPrefix("DECOY"), protein_decoy_weight});
  const auto protein_qvalues = targetDecoyQvalues(items);
  std::map<String, double> protein_q;
  for (Size k = 0; k < ranked.size(); k++) protein_q[ranked[k].second] = protein_qvalues[k];
  for (auto& hit : proteoform_hits_)
  {
    const auto it = protein_q.find(hit.getAccession());
    hit.setMetaValue("ProteinQvalue", it == protein_q.end() ? 1.0 : it->second);
  }
}

void FLASHTnTAlgorithm::vectorizeProteinSequence_(const std::vector<std::string>& cleaned_protein_seqs,
                                                     std::vector<std::unordered_set<int>>& vec_pro,
                                                     std::vector<std::unordered_set<int>>& rev_vec_pro)
{
  vec_pro.reserve(cleaned_protein_seqs.size());
  rev_vec_pro.reserve(cleaned_protein_seqs.size());
  for (const auto& seq : cleaned_protein_seqs)
  { //ignore X
    std::unordered_set<int> vec, rev_vec;
    double mass = 0;
    vec.insert(0);
    rev_vec.insert(0);

    for (Size j = 0; j < seq.size(); j++)
    {
      Size index = j;
      mass += ResidueDB::getInstance()->getResidue(seq[index])->getMonoWeight(Residue::Internal);
      int vindex = int(round(mass));
      vec.insert(vindex);
    }
    mass = 0;
    for (Size j = 0; j < seq.size(); j++)
    {
      Size index = seq.size() - j - 1;
      mass += ResidueDB::getInstance()->getResidue(seq[index])->getMonoWeight(Residue::Internal);
      int vindex = int(round(mass));
      rev_vec.insert(vindex);
    }
    vec_pro.push_back(vec);
    rev_vec_pro.push_back(rev_vec);
  }
}

void FLASHTnTAlgorithm::run(const MSExperiment& map, const std::vector<FASTAFile::FASTAEntry>& original_fasta_entry)
{
  setLogType(CMD);

  int max_mod_cntr = extender_param_.getValue("max_mod_count");
  double max_mod_mass = max_mod_cntr * (double)extender_param_.getValue("max_mod_mass") + 1.0;
  std::map<double, std::vector<ResidueModification>> blind_mod_map;
  const auto inst = ModificationsDB::getInstance();             // give this from outside ...
  std::map<String, std::vector<Size>> tag_to_protein_indices;   // tag to protein index in fasta
  std::map<String, std::vector<Size>> tag_to_protein_positions; // tag to protein position
  std::map<String, std::vector<double>> tag_to_n_flanking_masses, tag_to_c_flanking_masses; // tag to protein flanking masses

  std::map<int, std::vector<Size>> scan_to_tag_indices;         // scan number to tag indices
  std::vector<std::string> tag_seqs, protein_seqs, cleaned_protein_seqs;
  std::vector<DeconvolvedSpectrum> dspecs;

  std::vector<FASTAFile::FASTAEntry> fasta_entry;
  fasta_entry.reserve(original_fasta_entry.size());
  for (const auto& fe : original_fasta_entry)
  {
    auto seq = fe.sequence;
    std::replace(seq.begin(), seq.end(), 'I', 'L');
    fasta_entry.emplace_back(fe.identifier, fe.description, seq);
  }

  {
    protein_seqs.reserve(fasta_entry.size());
    cleaned_protein_seqs.reserve(fasta_entry.size());
    double taget_count = 0;
    double decoy_count = 0;
    for (const auto& fe : fasta_entry)
    {
      if (fe.identifier.hasPrefix("DECOY")) { decoy_count++; }
      else { taget_count++; }
      protein_seqs.push_back(fe.sequence);
      String cleaned_seq;
      std::remove_copy(fe.sequence.begin(), fe.sequence.end(),
                       std::back_inserter(cleaned_seq), 'X');
      cleaned_protein_seqs.push_back(cleaned_seq);
    }

    decoy_factor_ = decoy_count / taget_count;
  }

  std::vector<String> mod_strs;
  inst->getAllSearchModifications(mod_strs);
  for (int i = 0; i < mod_strs.size(); i++)
  {
    const auto mod = *inst->getModification(mod_strs[i]);
    if (std::abs(mod.getDiffMonoMass()) > max_mod_mass) continue;
    blind_mod_map[mod.getDiffMonoMass()].push_back(mod);
  }

  double precursor_tol = -1, tol, ms2_tol = -1;

  // Only MSn spectra are searched: the masses of an MS1 spectrum are whole proteoforms, not fragments, and its PrSMs (2208 of 10495 on
  // Eclip_1480, 4 of them passing 1%) only fed the decoy counts of the FDR. A file without any MSn spectrum is searched as it is.
  bool has_msn = false;
  for (const auto& spec : map)
  {
    if (spec.getMSLevel() > 1)
    {
      has_msn = true;
      break;
    }
  }

  startProgress(0, (SignedSize)map.size(), "Finding sequence tags ...");

  for (int index = 0; index < map.size(); index++)
  {
    auto spec = map[index];
    nextProgress();

    if (spec.size() < 5) continue;

    int scan = FLASHDeconvAlgorithm::getScanNumber(map, index);

    DeconvolvedSpectrum dspec(scan);
    dspec.setOriginalSpectrum(spec);
    String deconv_meta_str = spec.getMetaValue("DeconvMassInfo").toString();

    int tol_loc_s = deconv_meta_str.find("tol=") + 4;
    int tol_loc_e = deconv_meta_str.find(";", tol_loc_s);

    tol = stod(deconv_meta_str.substr(tol_loc_s, tol_loc_e - tol_loc_s));
    if (spec.getMSLevel() > 1 && ms2_tol < 0) ms2_tol = tol;
    if (spec.getMSLevel() == 1)
    {
      if (precursor_tol < 0) precursor_tol = tol;
      if (has_msn) continue;
    }

    int q_loc_s = deconv_meta_str.find("qscore=") + 7;
    int q_loc_e = deconv_meta_str.find(";", q_loc_s);
    auto q_str = deconv_meta_str.substr(q_loc_s, q_loc_e - q_loc_s);
    Size pos = 0;
    std::vector<double> qscores;
    while (true)
    {
      Size pos_t = q_str.find(",", pos);
      if (pos_t == String::npos) break;
      auto token = q_str.substr(pos, pos_t - pos);
      qscores.push_back(stod(token));
      pos = pos_t + 1;
    }

    int s_loc_s = deconv_meta_str.find("snr=") + 4;
    int s_loc_e = deconv_meta_str.find(";", s_loc_s);
    auto s_str = deconv_meta_str.substr(s_loc_s, s_loc_e - s_loc_s);
    pos = 0;
    std::vector<float> snrs;
    while (true)
    {
      Size pos_t = s_str.find(",", pos);
      if (pos_t == String::npos) break;
      auto token = s_str.substr(pos, pos_t - pos);
      snrs.push_back(stof(token));
      pos = pos_t + 1;
    }

    int s_loc_pre_s = deconv_meta_str.find("precursorscan=") + 14;
    int s_loc_pre_e = deconv_meta_str.find(";", s_loc_pre_s);
    int precursor_scan = stoi(deconv_meta_str.substr(s_loc_pre_s, s_loc_pre_e - s_loc_pre_s));

    double precursor_mass = 0;
    if (deconv_meta_str.hasSubstring("precursormass="))
    {
      int s_loc_prem_s = deconv_meta_str.find("precursormass=") + 14;
      int s_loc_prem_e = deconv_meta_str.find(";", s_loc_prem_s);
      precursor_mass = stod(deconv_meta_str.substr(s_loc_prem_s, s_loc_prem_e - s_loc_prem_s));
    }

    if (deconv_meta_str.hasSubstring("decoyprecursormass="))
    {
      int s_loc_dm_s = deconv_meta_str.find("decoyprecursormass=") + 19;
      int s_loc_dm_e = deconv_meta_str.find(";", s_loc_dm_s);
      dspec.setDecoyPrecursorMass(stod(deconv_meta_str.substr(s_loc_dm_s, s_loc_dm_e == (int)String::npos ? String::npos : s_loc_dm_e - s_loc_dm_s)));
      if (deconv_meta_str.hasSubstring("decoyprecursorweight="))
      {
        int s_loc_dw_s = deconv_meta_str.find("decoyprecursorweight=") + 21;
        int s_loc_dw_e = deconv_meta_str.find(";", s_loc_dw_s);
        dspec.setDecoyPrecursorWeight(stod(deconv_meta_str.substr(s_loc_dw_s, s_loc_dw_e == (int)String::npos ? String::npos : s_loc_dw_e - s_loc_dw_s)));
      }
    }

    // a precursor mass may be given (e.g. by precursor_charge/precursor_mz options) even if no precursor scan is present
    if (precursor_scan > 0 || precursor_mass > 0)
    {
      PeakGroup pg;
      pg.setMonoisotopicMass(precursor_mass);
      if (deconv_meta_str.hasSubstring("precursorQscore="))
      {
        int s_loc_preq_s = deconv_meta_str.find("precursorQscore=") + 16;
        int s_loc_preq_e = deconv_meta_str.find(";", s_loc_preq_s);
        double precursor_qscore = stod(deconv_meta_str.substr(s_loc_preq_s, s_loc_preq_e - s_loc_preq_s));
        pg.setQscore2D(precursor_qscore);
      }
      if (deconv_meta_str.hasSubstring("precursorSNR="))
      {
        int s_loc_pres_s = deconv_meta_str.find("precursorSNR=") + 13;
        int s_loc_pres_e = deconv_meta_str.find(";", s_loc_pres_s);
        double precursor_snr = stod(deconv_meta_str.substr(s_loc_pres_s, s_loc_pres_e - s_loc_pres_s));
        pg.setSNR(precursor_snr);
      }

      dspec.setPrecursorPeakGroup(pg);
    }
    for (Size i = 0; i < spec.size(); i++)
    {
      PeakGroup peak;
      peak.setQscore(qscores[i]);
      peak.setSNR(snrs[i]);
      peak.setMonoisotopicMass(spec[i].getMZ());
      peak.setScanNumber(scan);
      dspec.push_back(peak);
    }

    if (dspec.size() > max_node_cntr_)
    {
      dspec.sortByQscore();
      std::vector<PeakGroup> filtered_peaks;
      filtered_peaks.reserve(max_node_cntr_);
      for (const auto& pg : dspec)
      {
        if (filtered_peaks.size() >= max_node_cntr_) break;
        filtered_peaks.push_back(pg);
      }
      dspec.setPeakGroups(filtered_peaks);
    }

    dspec.sort();
    dspecs.push_back(dspec);
    FLASHTaggerAlgorithm tagger;
    // Run tagger
    tagger.setParameters(tagger_param_);
    tagger.run(dspec, tol);
    tagger.fillTags(tags_);
  }

  endProgress();

  for (Size i = 0; i < tags_.size(); i++)
  {
    const auto& seq = tags_[i].getUppercaseSequence();
    scan_to_tag_indices[tags_[i].getScan()].push_back(i);
    tags_[i].setIndex((int)i);
    if (tag_to_protein_indices.find(seq) != tag_to_protein_indices.end()) continue;
    tag_to_protein_indices[seq] = std::vector<Size>();
    tag_to_protein_positions[seq] = std::vector<Size>();
    tag_to_n_flanking_masses[seq] = std::vector<double>();
    tag_to_c_flanking_masses[seq] = std::vector<double>();

    tag_seqs.push_back(seq);
  } // now we have all the tags

  // tag_seqs protein_seqs

  ACTrie ac_trie(0, 0); // no ambiguities, no mismatches
  ac_trie.addNeedlesAndCompress(tag_seqs);
  OpenMS::ACTrieState ac_state;

  startProgress(0, (SignedSize)protein_seqs.size(), "Searching tags against database ...");

  for (int i = 0; i < protein_seqs.size(); i++) // // Natural19WithoutI
  {
    nextProgress();
    const String& pseq = protein_seqs[i];
    const String& cleanseq = cleaned_protein_seqs[i];
    ac_state.setQuery(pseq);
    ac_trie.getAllHits(ac_state);

    for (const auto& h : ac_state.hits)
    {
      String seq = tag_seqs[h.needle_index];
      tag_to_protein_indices[seq].push_back(i);
      tag_to_protein_positions[seq].push_back(h.query_pos);

      double nmass = h.query_pos <= 0 ? .0 : AASequence::fromString(cleanseq.substr(0, h.query_pos)).getMonoWeight(Residue::Internal);
      double cmass = cleanseq.length() - seq.length() <= h.query_pos ? .0 :
         AASequence::fromString(cleanseq.substr(h.query_pos + seq.length(), cleanseq.length() - h.query_pos - seq.length())).getMonoWeight(Residue::Internal);

      tag_to_n_flanking_masses[seq].push_back(nmass);
      tag_to_c_flanking_masses[seq].push_back(cmass);
    }
  }

  endProgress();
  startProgress(0, (SignedSize)dspecs.size(), "Running candidate protein filtration ...");

  std::vector<std::unordered_set<int>> vec_pro, rev_vec_pro;
  vectorizeProteinSequence_(cleaned_protein_seqs, vec_pro, rev_vec_pro);

  // exact residue sums from the N- and the C-terminus per protein, for the protein-level score (FLASHTaggerAlgorithm::runMatchingExact)
  std::vector<std::vector<double>> pro_prefix(cleaned_protein_seqs.size()), pro_suffix(cleaned_protein_seqs.size());
  for (Size k = 0; k < cleaned_protein_seqs.size(); k++)
  {
    const auto& seq = cleaned_protein_seqs[k];
    auto& pre = pro_prefix[k];
    auto& suf = pro_suffix[k];
    pre.reserve(seq.size() + 1);
    suf.reserve(seq.size() + 1);
    pre.push_back(0);
    suf.push_back(0);
    for (Size j = 0; j < seq.size(); j++) pre.push_back(pre.back() + ResidueDB::getInstance()->getResidue(seq[j])->getMonoWeight(Residue::Internal));
    for (Size j = 0; j < seq.size(); j++) suf.push_back(suf.back() + ResidueDB::getInstance()->getResidue(seq[seq.size() - 1 - j])->getMonoWeight(Residue::Internal));
  }
  const double protein_tol = protein_tol_ > 0 ? protein_tol_ : (ms2_tol > 0 ? ms2_tol : tol);

  // The protein-level search runs on every spectrum before any extension, because the run-level protein FDR
  // that decides which candidates are extended needs all of it. Each spectrum's candidates stay in
  // protein_hits_, in score order, at [first, second).
  std::vector<std::pair<Size, Size>> candidate_ranges(dspecs.size(), {0, 0});
  // Spectra are independent here, so they run in parallel (the shared maps are only read); each spectrum's candidates are collected
  // and appended to protein_hits_ in spectrum order afterwards, so the result does not depend on the thread schedule.
  std::vector<std::vector<ProteinHit>> spec_hits(dspecs.size());
  // candidates from terminal ladders, kept apart from the others (see the extension); their ranges in ladder_candidates
  std::vector<std::vector<ProteinHit>> spec_ladder_hits(dspecs.size());
  std::vector<ProteinHit> ladder_candidates;
  std::vector<std::pair<Size, Size>> ladder_ranges(dspecs.size(), {0, 0});
  // Meta values set by name look the name up in the global registry under a critical section; register the names once and set
  // them by index inside the parallel loop.
  auto& registry = MetaInfoInterface::metaRegistry();
  const UInt mi_scan = registry.registerName("Scan"), mi_matched_aa = registry.registerName("MatchedAA"),
             mi_tag_indices = registry.registerName("TagIndices"), mi_tag_positions = registry.registerName("TagPositions"),
             mi_fasta_index = registry.registerName("FastaIndex"), mi_protein_rank = registry.registerName("ProteinRank"),
             mi_gapped = registry.registerName("GappedTag"), mi_terminal = registry.registerName("TerminalLadder"),
             mi_precursor_ladder = registry.registerName("PrecursorLadder"),
             mi_pl_start = registry.registerName("PrecursorLadderStart"), mi_pl_end = registry.registerName("PrecursorLadderEnd"),
             mi_pl_isotope = registry.registerName("PrecursorLadderIsotope");
  // gapped tags (only for spectra in which no sequence tag proposes a protein); their tags are kept per spectrum until appended to tags_
  FLASHGappedTaggerAlgorithm gapped_tagger;
  if (gapped_mode_ != "off")
  {
    gapped_tagger.setParameters(gapped_param_);
    gapped_tagger.buildIndex(protein_seqs);
    OPENMS_LOG_INFO << "Gapped tag index: " << gapped_tagger.getIndexBytes() / (1024 * 1024) << " MB" << std::endl;
  }
  std::vector<std::vector<FLASHHelperClasses::Tag>> spec_gapped_tags(dspecs.size());

  // per ion type of ion_type: whether it is a prefix ion, and the shift that makes its fragment mass a residue sum (mass + shift)
  std::vector<std::pair<bool, double>> ion_residue_shifts;
  for (const auto& ion : param_.getValue("ion_type").toStringVector())
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
    ion_residue_shifts.emplace_back(prefix, -offset);
  }
  const std::string standard_residues = "ACDEFGHIKLMNPQRSTVWYUO"; // the residues of terminal and precursor ladders

  // Terminal ladders: the residue sums of every protein from its N- and from its C-terminus (up to the first non-standard residue), and
  // all of them sorted by mass. A deconvolved mass plus a shift (of its ion type and, at the N-terminus, of the forms without the
  // initiator Met and with acetylation) that lands on a residue sum reads from the intact terminus; no sequence tag is needed.
  struct LadderEntry
  {
    double mass;
    std::uint32_t protein;
  };
  std::vector<std::vector<double>> ladder_pre, ladder_suf;
  std::vector<LadderEntry> ladder_pre_index, ladder_suf_index;
  std::vector<double> n_ladder_shifts, c_ladder_shifts; // observed mass + shift = residue sum from the terminus
  if (terminal_ladders_)
  {
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
    ladder_pre.resize(protein_seqs.size());
    ladder_suf.resize(protein_seqs.size());
    for (Size k = 0; k < protein_seqs.size(); k++)
    {
      sums(protein_seqs[k], false, ladder_pre[k]);
      sums(protein_seqs[k], true, ladder_suf[k]);
      for (Size j = 1; j < ladder_pre[k].size(); j++) ladder_pre_index.push_back({ladder_pre[k][j], (std::uint32_t)k});
      for (Size j = 1; j < ladder_suf[k].size(); j++) ladder_suf_index.push_back({ladder_suf[k][j], (std::uint32_t)k});
    }
    auto by_mass = [](const LadderEntry& a, const LadderEntry& b) { return a.mass < b.mass; };
    std::sort(ladder_pre_index.begin(), ladder_pre_index.end(), by_mass);
    std::sort(ladder_suf_index.begin(), ladder_suf_index.end(), by_mass);

    const double met = rdb->getResidue('M')->getMonoWeight(Residue::Internal);
    const double acetyl = EmpiricalFormula("C2H2O").getMonoWeight();
    for (const auto& [prefix, shift] : ion_residue_shifts)
    {
      if (prefix)
      {
        for (double form : {0.0, met, -acetyl, met - acetyl}) n_ladder_shifts.push_back(form + shift); // form - offset
      }
      else { c_ladder_shifts.push_back(shift); }
    }
  }

  // Precursor-anchored ladders: the database as runs of standard residues (segments) with their residue sums, for the window search of
  // every spectrum with a precursor mass. A protein of standard residues only uses its pro_prefix sums (no copy); a protein with other
  // residues gets its own sums, in which those add 0 (no window crosses them). Built only if a spectrum has a precursor mass.
  bool any_precursor_mass = false;
  for (const auto& dspec : dspecs) any_precursor_mass |= dspec.getPrecursorPeakGroup().getMonoMass() > 0;
  const bool use_precursor_ladders = precursor_ladders_ && ! protein_level_only_ && any_precursor_mass;
  std::vector<precursor_ladder::Segment> pl_segments;
  std::vector<const double*> pl_sums;
  std::vector<std::vector<double>> pl_own_sums;
  // the precursor tolerance of the windows: the MS1 tolerance of the input, else (no MS1 spectra) the fragment tolerance
  const double pl_precursor_tol = precursor_tol > 0 ? precursor_tol : protein_tol;
  const double h2o = Residue::getInternalToFull().getMonoWeight();
  if (use_precursor_ladders)
  {
    const auto* rdb = ResidueDB::getInstance();
    std::vector<Size> other;
    for (Size k = 0; k < protein_seqs.size(); k++)
    {
      if (protein_seqs[k].find_first_not_of(standard_residues) != std::string::npos) other.push_back(k);
    }
    pl_sums.assign(protein_seqs.size(), nullptr);
    pl_own_sums.resize(other.size()); // not resized below: the pointers into it stay valid
    for (Size n = 0; n < other.size(); n++)
    {
      const auto& seq = protein_seqs[other[n]];
      auto& sums = pl_own_sums[n];
      sums.reserve(seq.size() + 1);
      sums.push_back(0);
      for (const char aa : seq)
      {
        sums.push_back(sums.back() + (standard_residues.find(aa) == std::string::npos ? 0.0 : rdb->getResidue(aa)->getMonoWeight(Residue::Internal)));
      }
      pl_sums[other[n]] = sums.data();
    }
    for (Size k = 0; k < protein_seqs.size(); k++)
    {
      if (pl_sums[k] == nullptr) pl_sums[k] = pro_prefix[k].data(); // cleaned_protein_seqs[k] == protein_seqs[k] (no X)
      const auto& seq = protein_seqs[k];
      long a = -1;
      for (Size x = 0; x <= seq.size(); x++)
      {
        const bool standard = x < seq.size() && standard_residues.find(seq[x]) != std::string::npos;
        if (standard && a < 0) a = (long)x;
        if (! standard && a >= 0)
        {
          pl_segments.push_back({(std::uint32_t)k, (std::uint32_t)a, (std::uint32_t)x});
          a = -1;
        }
      }
    }
  }
  // per spectrum: the windows that fit its precursor mass (hypotheses), and the candidates; the same at the decoy precursor mass
  // (decoy_precursor), so that the decoy copy of a spectrum has this source as its target does
  std::vector<std::vector<ProteinHit>> spec_pl_hits(dspecs.size()), spec_decoy_pl_hits(dspecs.size());
  std::vector<ProteinHit> pl_candidates, decoy_pl_candidates;
  std::vector<std::pair<Size, Size>> pl_ranges(dspecs.size(), {0, 0}), decoy_pl_ranges(dspecs.size(), {0, 0});
  std::vector<Size> spec_pl_windows(dspecs.size(), 0);
  // Parallel over spectra only when there are at least as many spectra with tags as threads; otherwise (e.g. a targeted native run with
  // one or a few spectra) the threads go to the scoring of each spectrum's candidates (runMatchingExact).
  int spectra_with_tags = 0;
  for (const auto& dspec : dspecs)
  {
    spectra_with_tags += gapped_mode_ != "off" || terminal_ladders_ || (use_precursor_ladders && dspec.getPrecursorPeakGroup().getMonoMass() > 0)
                         || scan_to_tag_indices.count(dspec.getScanNumber()) > 0;
  }
#ifdef _OPENMP
  const bool parallel_over_spectra = spectra_with_tags >= omp_get_max_threads();
#else
  const bool parallel_over_spectra = false;
#endif
#pragma omp parallel for schedule(dynamic) if (parallel_over_spectra)
  for (int d = 0; d < (int)dspecs.size(); d++)
  {
    const auto& dspec = dspecs[d];
    nextProgress();
    const auto tag_it = scan_to_tag_indices.find(dspec.getScanNumber());
    std::vector<ProteinHit> hits;
    std::vector<FLASHTaggerAlgorithm::ProteinShifts> hit_shifts; // in the order of hits
    if (gapped_mode_ != "only" && tag_it != scan_to_tag_indices.end())
    {
      const auto& tag_indices = tag_it->second;
      std::map<Size, std::set<Size>> pi_to_pos;         // protein index to covered positions
      std::map<Size, std::set<Size>> pi_to_tag_indices; // protein index to tag indices
      // per protein: shifts (protein residue sum up to the tag minus the tag flanking mass) and the tag flanking masses
      std::map<Size, std::vector<double>> pi_to_n_shifts, pi_to_n_flanks, pi_to_c_shifts, pi_to_c_flanks;
      for (const auto& tag_index : tag_indices)
      {
        const auto& tag = tags_[tag_index];
        const auto& seq = tag.getUppercaseSequence();

        const auto& pis = tag_to_protein_indices.at(seq);
        const auto& pps = tag_to_protein_positions.at(seq);

        for (Size i = 0; i < pis.size(); i++)
        {
          const auto& pseq = cleaned_protein_seqs[pis[i]];

          if (tag.getNtermMass() >= 0 && pps[i] > 0 && pps[i] <= pseq.length() + seq.length())
          {
            double nmass = tag_to_n_flanking_masses.at(seq)[i];
            if (tag.getNtermMass() > nmass + max_mod_mass + 200) continue;
            pi_to_n_shifts[pis[i]].push_back(nmass - tag.getNtermMass());
            pi_to_n_flanks[pis[i]].push_back(tag.getNtermMass());
          }
          if (tag.getCtermMass() >= 0 && pps[i] + seq.length() >= 0 && pseq.length() > pps[i] + seq.length())
          {
            double cmass = tag_to_c_flanking_masses.at(seq)[i];
            if (tag.getCtermMass() > cmass + max_mod_mass + 200) continue;
            pi_to_c_shifts[pis[i]].push_back(cmass - tag.getCtermMass());
            pi_to_c_flanks[pis[i]].push_back(tag.getCtermMass());
          }

          pi_to_tag_indices[pis[i]].insert(tag_index);
          pi_to_pos[pis[i]].insert(pps[i]);
        }
      }

      for (const auto& [pi, pp] : pi_to_pos)
      {
        const auto& fe = fasta_entry[pi];
        Size covered_pos = pp.size();
        if (covered_pos < min_tag_count) continue;

        // check flanking masses ...
        double coverage = (double)covered_pos / (double)fe.sequence.length();
        std::vector<int> positions(pp.begin(), pp.end());

        //if (positions.size() < min_tag_count) continue;
        ProteinHit hit(coverage, 0, fe.identifier, fe.sequence); //
        hit.setDescription(fe.description);
        hit.setMetaValue(mi_scan, dspec.getScanNumber());
        hit.setMetaValue(mi_matched_aa, covered_pos);
        hit.setCoverage(coverage);
        FLASHTaggerAlgorithm::ProteinShifts sh;
        sh.fasta_index = pi;
        sh.n_shifts = std::move(pi_to_n_shifts[pi]);
        sh.n_flanks = std::move(pi_to_n_flanks[pi]);
        sh.c_shifts = std::move(pi_to_c_shifts[pi]);
        sh.c_flanks = std::move(pi_to_c_flanks[pi]);
        hit_shifts.push_back(std::move(sh));

        hit.setMetaValue(mi_tag_indices, std::vector<int>(pi_to_tag_indices[pi].begin(), pi_to_tag_indices[pi].end()));
        hit.setMetaValue(mi_tag_positions, positions);
        hit.setMetaValue(mi_fasta_index, pi);
        hits.push_back(hit);
      }

    }

    // gapped tags: proteins proposed where gapped tags land. Candidates are scored and cut first; only the survivors get their
    // sequence and their landings as tags (a landing becomes a tag spanning the landed residues, with the masses at the gap
    // boundaries), so the protein score (shifts) and the extension (anchors) use it as they use a sequence tag.
    std::vector<FLASHGappedTaggerAlgorithm::Landing> landings;
    std::vector<std::vector<Size>> gapped_groups; // per gapped candidate, its landings
    if (gapped_mode_ == "only" || gapped_mode_ == "both" || (gapped_mode_ == "fallback" && hits.empty()))
    {
      std::set<Size> proposed; // proteins already proposed by sequence tags
      for (const auto& h : hits) proposed.insert((Size)(int)h.getMetaValue(mi_fasta_index));
      gapped_tagger.search(dspec, ms2_tol > 0 ? ms2_tol : tol, landings);
      for (Size l0 = 0; l0 < landings.size();)
      {
        Size l1 = l0;
        while (l1 < landings.size() && landings[l1].protein == landings[l0].protein) ++l1; // landings are sorted by protein
        const Size pi = landings[l0].protein;
        FLASHTaggerAlgorithm::ProteinShifts sh;
        sh.fasta_index = pi;
        std::vector<Size> group;
        std::vector<int> positions;
        for (Size l = l0; l < l1 && proposed.count(pi) == 0; ++l)
        {
          const auto& la = landings[l];
          // the lowest mass borders the residues before (prefix) or after (suffix) the landed span
          const double flank = la.masses[0];
          if (flank > la.flank_mass + max_mod_mass + 200) continue;
          if (la.prefix)
          {
            sh.n_shifts.push_back(la.flank_mass - flank);
            sh.n_flanks.push_back(flank);
          }
          else
          {
            sh.c_shifts.push_back(la.flank_mass - flank);
            sh.c_flanks.push_back(flank);
          }
          group.push_back(l);
          positions.push_back(la.start);
        }
        l0 = l1;
        if (group.empty()) continue;
        std::sort(positions.begin(), positions.end());
        positions.erase(std::unique(positions.begin(), positions.end()), positions.end());
        const auto& fe = fasta_entry[pi];
        const double coverage = (double)positions.size() / (double)fe.sequence.length();
        ProteinHit hit(coverage, 0, fe.identifier, ""); // the sequence is set if the candidate survives the cut
        hit.setDescription(fe.description);
        hit.setMetaValue(mi_scan, dspec.getScanNumber());
        hit.setMetaValue(mi_matched_aa, (int)positions.size());
        hit.setCoverage(coverage);
        hit.setMetaValue(mi_tag_indices, std::vector<int>(group.size(), 0)); // as many as its tags: the cut below counts them
        hit.setMetaValue(mi_tag_positions, positions);
        hit.setMetaValue(mi_fasta_index, pi);
        hit.setMetaValue(mi_gapped, (int)gapped_groups.size());
        gapped_groups.push_back(std::move(group));
        hits.push_back(hit);
        hit_shifts.push_back(std::move(sh));
      }
    }

    // exact masses at the deconvolution tolerance, scored by the chance probability of the match count: 1 Da integer bins and summed
    // node scores let chance matches (more of them on longer proteins) rank noise proteins above the correct one
    FLASHTaggerAlgorithm::runMatchingExact(hits, hit_shifts, dspec, pro_prefix, pro_suffix, protein_tol);

    std::sort(hits.begin(), hits.end(), [](const ProteinHit& left, const ProteinHit& right) {
      return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getDescription() > right.getDescription())
                                                                                              : (left.getCoverage() > right.getCoverage()))
                                                 : (left.getScore() > right.getScore());
    });

    if ((int)hits.size() > max_protein_hits_)
    {
      // The score above uses the masses of the whole spectrum against the protein; a correct protein with few of them can end up below
      // the cut. The number of distinct matched tags is sequence-specific evidence, so keep the best hits by that criterion as well and
      // let the extension decide between them.
      std::vector<bool> keep(hits.size(), false);
      for (int i = 0; i < max_protein_hits_; i++)
      {
        keep[i] = true; // hits are already sorted by score
      }

      std::vector<std::pair<int, Size>> by_tag_count;
      by_tag_count.reserve(hits.size());
      for (Size i = 0; i < hits.size(); i++)
      {
        by_tag_count.emplace_back(-(int)hits[i].getMetaValue(mi_tag_indices).toIntList().size(), i);
      }
      std::partial_sort(by_tag_count.begin(), by_tag_count.begin() + max_protein_hits_, by_tag_count.end());
      for (int i = 0; i < max_protein_hits_; i++)
      {
        keep[by_tag_count[i].second] = true;
      }

      std::vector<ProteinHit> kept;
      kept.reserve(2 * max_protein_hits_);
      for (Size i = 0; i < hits.size(); i++)
      {
        if (keep[i]) kept.push_back(hits[i]);
      }
      hits.swap(kept); // still in score order
    }

    // the surviving gapped candidates: sequence and tags (indices local to this spectrum until appended to tags_ below)
    auto& gtags = spec_gapped_tags[d];
    for (auto& hit : hits)
    {
      if (! hit.metaValueExists(mi_gapped)) continue;
      const Size pi = (Size)(int)hit.getMetaValue(mi_fasta_index);
      hit.setSequence(fasta_entry[pi].sequence);
      std::vector<int> indices;
      for (Size l : gapped_groups[(int)hit.getMetaValue(mi_gapped)])
      {
        auto& la = landings[l];
        const double flank = la.masses[0];
        gtags.emplace_back(String(protein_seqs[pi].substr(la.start, la.end - la.start)), la.prefix ? flank : -1.0, la.prefix ? -1.0 : flank,
                           la.masses, la.scores, dspec.getScanNumber(), (Size)(la.end - la.start));
        gtags.back().setRetentionTime((float)dspec.getOriginalSpectrum().getRT());
        indices.push_back((int)gtags.size() - 1);
      }
      hit.setMetaValue(mi_tag_indices, indices);
      hit.setMetaValue(mi_gapped, 1);
    }

    // Terminal ladders: per protein and shift, the masses landing on a residue sum from the terminus. A protein is scored as
    // runMatchingExact scores a side (the chance probability of the match count, times the number of shifts of the side), without the
    // terminus itself: its match would only tell proteins starting with Met from shuffled decoys. The min_protein_hits best proteins are
    // characterized; the ladder of the best shift of each side becomes a tag spanning it (as a gapped tag does), so the extension is
    // anchored there. They are kept apart from the other candidates, also when a sequence tag proposed the same protein.
    if (terminal_ladders_ && ! dspec.empty())
    {
      const Size n_proteins = protein_seqs.size(), n_n = n_ladder_shifts.size(), n_shifts = n_n + c_ladder_shifts.size();
      std::vector<double> masses; // ascending
      masses.reserve(dspec.size());
      for (const auto& pg : dspec) masses.push_back(pg.getMonoMass());
      auto shift_of = [&](Size si) { return si < n_n ? n_ladder_shifts[si] : c_ladder_shifts[si - n_n]; };
      std::vector<std::uint16_t> counts(n_proteins * n_shifts, 0);
      std::vector<Size> touched;
      for (Size si = 0; si < n_shifts; si++)
      {
        const auto& index = si < n_n ? ladder_pre_index : ladder_suf_index;
        for (double m : masses)
        {
          const double t = m + shift_of(si);
          const double margin = protein_tol * 1e-6 * std::max(m, t);
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
        const auto& sums = si < n_n ? ladder_pre[p] : ladder_suf[p];
        const double s = shift_of(si);
        const Size lo = std::lower_bound(masses.begin(), masses.end(), -s) - masses.begin();
        const Size hi = std::upper_bound(masses.begin(), masses.end(), sums.back() - s) - masses.begin();
        const double lambda = hi > lo ? 2.0 * protein_tol * 1e-6 * (double)(sums.size() - 1) / sums.back() * (cum_weight[si][hi] - cum_weight[si][lo]) : 0;
        const double score = -std::log10(std::max(boost::math::gamma_p((double)k, std::max(lambda, 1e-12)), 1e-300));
        auto& sc = side_scores[p];
        auto& bs = best_shifts.emplace(p, std::make_pair(n_shifts, n_shifts)).first->second;
        if (si < n_n && score > sc.first) { sc.first = score; bs.first = si; }
        else if (si >= n_n && score > sc.second) { sc.second = score; bs.second = si; }
      }
      std::vector<std::pair<double, Size>> ranked; // score, protein
      for (const auto& [p, sc] : side_scores)
      {
        const double score = std::max(n_n > 0 ? sc.first - std::log10((double)n_n) : 0.0, n_shifts > n_n ? sc.second - std::log10((double)(n_shifts - n_n)) : 0.0);
        if (score > 0) ranked.emplace_back(score, p);
      }
      const Size top = std::min(ranked.size(), (Size)min_protein_hits_);
      std::partial_sort(ranked.begin(), ranked.begin() + top, ranked.end(),
                        [](const auto& a, const auto& b) { return a.first != b.first ? a.first > b.first : a.second < b.second; });
      auto& gtags = spec_gapped_tags[d];
      for (Size r = 0; r < top; r++)
      {
        const double score = ranked[r].first;
        const Size p = ranked[r].second;
        const auto& pseq = protein_seqs[p];
        const Size bn = best_shifts[p].first, bc = best_shifts[p].second;
        std::vector<int> indices, positions; // tag indices: -(local index + 1), made global with the gapped tags below
        int matched = 0;
        for (Size si : {bn, bc})
        {
          if (si >= n_shifts) continue;
          const bool prefix = si < n_n;
          const auto& sums = prefix ? ladder_pre[p] : ladder_suf[p];
          const double s = shift_of(si);
          std::vector<double> tag_masses;
          std::vector<int> tag_scores;
          int first = -1, last = -1; // residues from the terminus at the lowest and the highest ladder mass
          for (const auto& pg : dspec)
          {
            const double m = pg.getMonoMass(), t = m + s, margin = protein_tol * 1e-6 * std::max(m, t);
            auto it = std::lower_bound(sums.begin() + 1, sums.end(), t - margin);
            if (it == sums.end() || *it > t + margin) continue;
            const int j = (int)(it - sums.begin());
            // the whole sequence is no fragment: that mass is the precursor (y) or the precursor minus water (b)
            if (j <= last || j >= (int)pseq.size()) continue;
            if (first < 0) first = j;
            last = j;
            tag_masses.push_back(m);
            tag_scores.push_back(FLASHTaggerAlgorithm::getNodeScore(pg));
          }
          if (tag_masses.size() < 2) continue;
          matched += (int)tag_masses.size();
          const int start = prefix ? first : (int)pseq.size() - last, end = prefix ? last : (int)pseq.size() - first;
          gtags.emplace_back(String(pseq.substr(start, end - start)), prefix ? tag_masses[0] : -1.0, prefix ? -1.0 : tag_masses[0], tag_masses,
                             tag_scores, dspec.getScanNumber(), (Size)(end - start));
          gtags.back().setRetentionTime((float)dspec.getOriginalSpectrum().getRT());
          indices.push_back(-(int)gtags.size());
          positions.push_back(start);
        }
        if (indices.empty()) continue;
        const auto& fe = fasta_entry[p];
        const double coverage = (double)matched / (double)fe.sequence.length();
        ProteinHit hit(score, 0, fe.identifier, fe.sequence);
        hit.setDescription(fe.description);
        hit.setMetaValue(mi_scan, dspec.getScanNumber());
        hit.setMetaValue(mi_matched_aa, matched);
        hit.setCoverage(coverage);
        hit.setMetaValue(mi_tag_indices, indices);
        std::sort(positions.begin(), positions.end());
        hit.setMetaValue(mi_tag_positions, positions);
        hit.setMetaValue(mi_fasta_index, p);
        hit.setMetaValue(mi_protein_rank, (int)r + 1);
        hit.setMetaValue(mi_terminal, 1);
        spec_ladder_hits[d].push_back(std::move(hit));
      }
    }

    // Precursor-anchored ladders: the windows of the database (targets and decoys alike) whose unmodified mass is the precursor mass, by
    // the chance probability of the masses on their prefix and suffix ladders (precursor_ladder::search). The min_protein_hits best
    // windows are characterized; the matched masses of the best prefix and the best suffix ion type become tags spanning them (as
    // terminal-ladder tags), from which the extension can start and end the proteoform at the window. FLASHExtender places the tags by
    // their sequence and reaches an internal start or end only within its truncation range around the protein's best mass shift
    // (getBestMassShift_), not at the window it is given: the window is stored on the candidate (PrecursorLadderStart and End, 1-based as
    // StartPosition and EndPosition, and PrecursorLadderIsotope), and the log line after the extension counts the PrSMs at its ends. Kept
    // apart from the other candidates; searched before the check for tag candidates, so spectra without any are searched too.
    if (use_precursor_ladders && ! dspec.empty() && dspec.getPrecursorPeakGroup().getMonoMass() > 0)
    {
      std::vector<double> masses; // ascending
      masses.reserve(dspec.size());
      for (const auto& pg : dspec) masses.push_back(pg.getMonoMass());
      std::vector<precursor_ladder::TargetSet> sets;
      sets.reserve(ion_residue_shifts.size());
      for (const auto& [prefix, shift] : ion_residue_shifts) sets.emplace_back(masses, shift, prefix, protein_tol);
      auto& gtags = spec_gapped_tags[d];
      precursor_ladder::Result found;
      auto candidatesAt = [&](double precursor_mass, std::vector<ProteinHit>& candidates) {
        precursor_ladder::search(sets, pl_segments, pl_sums, precursor_mass, pl_precursor_tol, h2o, Constants::C13C12_MASSDIFF_U,
                                 precursor_ladder_min_k, (Size)min_protein_hits_, found);
        const double log_hypotheses = std::log10((double)std::max(found.hypotheses, (std::size_t)1));
        for (Size r = 0; r < found.top.size(); r++)
        {
          const auto& win = found.top[r];
          const double* sums = pl_sums[win.protein];
          const auto& pseq = protein_seqs[win.protein];
          // per side, the ion type with the most matched positions (the first in ion_type among equal ones)
          int best_prefix = -1, best_suffix = -1;
          for (int s = 0; s < (int)sets.size(); s++)
          {
            int& best = sets[s].isPrefix() ? best_prefix : best_suffix;
            if (best < 0 || win.counts[s] > win.counts[best]) best = s;
          }
          std::vector<int> indices, positions; // tag indices: -(local index + 1), made global with the gapped tags below
          int matched = 0;
          for (int s : {best_prefix, best_suffix})
          {
            if (s < 0) continue;
            const auto& set = sets[s];
            const bool prefix = set.isPrefix();
            std::vector<double> tag_masses; // ascending: prefix positions from the start, suffix positions from the end
            std::vector<int> tag_scores, tag_positions;
            for (Size step = 1; step < win.j - win.i; step++)
            {
              const Size x = prefix ? win.i + step : win.j - step;
              const long mi = set.massIndex(prefix ? sums[x] - sums[win.i] : sums[win.j] - sums[x]);
              if (mi < 0) continue;
              tag_positions.push_back((int)x);
              tag_masses.push_back(masses[mi]);
              tag_scores.push_back(FLASHTaggerAlgorithm::getNodeScore(dspec[mi]));
            }
            if (tag_masses.size() < 2) continue;
            matched += (int)tag_masses.size();
            // the tag spans the positions of its lowest and its highest mass; its flanking mass is the lowest (the ladder read from i or j)
            const int start = prefix ? tag_positions.front() : tag_positions.back(), end = prefix ? tag_positions.back() : tag_positions.front();
            gtags.emplace_back(String(pseq.substr(start, end - start)), prefix ? tag_masses[0] : -1.0, prefix ? -1.0 : tag_masses[0], tag_masses,
                               tag_scores, dspec.getScanNumber(), (Size)(end - start));
            gtags.back().setRetentionTime((float)dspec.getOriginalSpectrum().getRT());
            indices.push_back(-(int)gtags.size());
            positions.push_back(start);
          }
          if (indices.empty()) continue; // only possible with more than one ion type per side
          const auto& fe = fasta_entry[win.protein];
          ProteinHit hit(win.log_p - log_hypotheses, 0, fe.identifier, fe.sequence);
          hit.setDescription(fe.description);
          hit.setMetaValue(mi_scan, dspec.getScanNumber());
          hit.setMetaValue(mi_matched_aa, matched);
          hit.setCoverage((double)matched / (double)fe.sequence.length());
          hit.setMetaValue(mi_tag_indices, indices);
          std::sort(positions.begin(), positions.end());
          hit.setMetaValue(mi_tag_positions, positions);
          hit.setMetaValue(mi_fasta_index, (Size)win.protein);
          hit.setMetaValue(mi_protein_rank, (int)r + 1);
          hit.setMetaValue(mi_precursor_ladder, 1);
          hit.setMetaValue(mi_pl_start, (int)win.i + 1); // residues i..j-1 (0-based) are i+1..j (1-based, inclusive)
          hit.setMetaValue(mi_pl_end, (int)win.j);
          hit.setMetaValue(mi_pl_isotope, win.isotope);
          candidates.push_back(std::move(hit));
        }
      };
      candidatesAt(dspec.getPrecursorPeakGroup().getMonoMass(), spec_pl_hits[d]);
      spec_pl_windows[d] = found.hypotheses;
      if (decoy_precursor_ && dspec.getDecoyPrecursorMass() > 0) candidatesAt(dspec.getDecoyPrecursorMass(), spec_decoy_pl_hits[d]);
    }

    if (hits.empty()) continue;
    for (Size i = 0; i < hits.size(); i++)
    {
      hits[i].setMetaValue(mi_protein_rank, (int)i + 1);
    }
    spec_hits[d] = std::move(hits);
  }
  for (Size d = 0; d < dspecs.size(); d++)
  {
    // gapped tags of this spectrum join tags_; the tag indices of its candidates become global
    const int tag_base = (int)tags_.size();
    for (auto& t : spec_gapped_tags[d])
    {
      t.setIndex((int)tags_.size());
      tags_.push_back(std::move(t));
    }
    candidate_ranges[d].first = protein_hits_.size();
    for (auto& hit : spec_hits[d])
    {
      if (hit.metaValueExists(mi_gapped))
      {
        auto indices = hit.getMetaValue(mi_tag_indices).toIntList();
        for (auto& i : indices) i += tag_base;
        hit.setMetaValue(mi_tag_indices, indices);
      }
      protein_hits_.push_back(std::move(hit));
    }
    candidate_ranges[d].second = protein_hits_.size();
    // ladder candidates (terminal and precursor ladders): their tag indices -(local index + 1) become global
    auto appendLadderCandidates = [&](std::vector<ProteinHit>& from, std::vector<ProteinHit>& to, std::pair<Size, Size>& range) {
      range.first = to.size();
      for (auto& hit : from)
      {
        auto indices = hit.getMetaValue(mi_tag_indices).toIntList();
        for (auto& i : indices) i = tag_base - i - 1;
        hit.setMetaValue(mi_tag_indices, indices);
        to.push_back(std::move(hit));
      }
      range.second = to.size();
    };
    appendLadderCandidates(spec_ladder_hits[d], ladder_candidates, ladder_ranges[d]);
    appendLadderCandidates(spec_pl_hits[d], pl_candidates, pl_ranges[d]);
    appendLadderCandidates(spec_decoy_pl_hits[d], decoy_pl_candidates, decoy_pl_ranges[d]);
  }
  std::vector<std::vector<ProteinHit>>().swap(spec_hits);
  std::vector<std::vector<ProteinHit>>().swap(spec_ladder_hits);
  std::vector<std::vector<ProteinHit>>().swap(spec_pl_hits);
  std::vector<std::vector<ProteinHit>>().swap(spec_decoy_pl_hits);

  endProgress();

  if (use_precursor_ladders)
  {
    std::vector<Size> windows; // per searched spectrum
    for (Size d = 0; d < dspecs.size(); d++)
    {
      if (dspecs[d].getPrecursorPeakGroup().getMonoMass() > 0 && ! dspecs[d].empty()) windows.push_back(spec_pl_windows[d]);
    }
    Size total = 0;
    for (Size w : windows) total += w;
    std::nth_element(windows.begin(), windows.begin() + windows.size() / 2, windows.end());
    String at_decoy_masses;
    if (decoy_precursor_) at_decoy_masses = ", " + String(decoy_pl_candidates.size()) + " at the decoy precursor masses";
    OPENMS_LOG_INFO << "Precursor ladders: " << pl_segments.size() << " segments; " << windows.size() << " spectra with a precursor mass, " << total
                    << " windows (median " << (windows.empty() ? 0 : windows[windows.size() / 2]) << " per spectrum), " << pl_candidates.size()
                    << " candidates (top " << min_protein_hits_ << ", at least " << precursor_ladder_min_k << " matched positions)" << at_decoy_masses
                    << std::endl;
  }

  if (! protein_level_only_)
  {
    // Beyond the min_protein_hits best candidates of a spectrum (a fair target-decoy competition: they are chosen by this spectrum
    // alone), extend only those whose protein, or its target-decoy partner, passes the run-level protein FDR. That needs decoys;
    // without them every candidate is extended, as before.
    const bool use_protein_fdr = decoy_factor_ > 0 && protein_fdr_ < 1.0;
    std::set<Size> passing_proteins;
    if (use_protein_fdr) { passing_proteins = getProteinsPassingFDR_(fasta_entry); }

    // spectra in parallel, each with its own extender; proteoforms are appended in spectrum order afterwards
    std::vector<std::vector<ProteinHit>> spec_proteoforms(dspecs.size());
    // as above: over spectra only with enough spectra to extend, otherwise over each spectrum's candidates inside the extender
    auto hasCandidates = [&](int d) {
      return candidate_ranges[d].first != candidate_ranges[d].second || ladder_ranges[d].first != ladder_ranges[d].second
             || pl_ranges[d].first != pl_ranges[d].second;
    };
    auto overSpectra = [&](const std::vector<int>& spectra) {
      [[maybe_unused]] const auto with_candidates = std::count_if(spectra.begin(), spectra.end(), hasCandidates);
#ifdef _OPENMP
      return with_candidates >= omp_get_max_threads();
#else
      return false;
#endif
    };
    std::vector<int> all_spectra(dspecs.size());
    std::iota(all_spectra.begin(), all_spectra.end(), 0);
    const bool extend_over_spectra = overSpectra(all_spectra);
    const bool keep_all_candidates = (bool)raw_loader_;
    // precursor-ladder PrSMs of the extension that is kept: kept (determined; of them at both window ends, at the window start only, at the
    // window end only), dropped as the same proteoform as a determined PrSM of another source, dropped as undetermined
    std::atomic<Size> pl_determined {0}, pl_at_both_ends {0}, pl_at_start {0}, pl_at_end {0}, pl_duplicates {0}, pl_undetermined {0};
    // The extension of one spectrum's candidates with the given extension parameters. With decoy_mass, the spectrum is the decoy copy at
    // the decoy precursor mass (decoy_precursor), and its precursor-ladder candidates are those of the windows that fit that mass.
    auto extendSpectrumInto = [&](int d, const Param& eparam, const DeconvolvedSpectrum& dspec, std::vector<ProteinHit>& out, bool decoy_mass) {
      out.clear();
      const auto [first, second] = candidate_ranges[d];
      const auto [pl_first, pl_second] = decoy_mass ? decoy_pl_ranges[d] : pl_ranges[d];
      if (first == second && ladder_ranges[d].first == ladder_ranges[d].second && pl_first == pl_second) return;

      std::vector<ProteinHit> hits;
      // the min_protein_hits best candidates of each origin (sequence tags, gapped tags) are always extended, so gapped candidates
      // add chances without taking the places of sequence-tag candidates
      int rank_tag = 0, rank_gapped = 0;
      for (Size i = first; i < second; i++)
      {
        const bool gapped = protein_hits_[i].metaValueExists(mi_gapped);
        bool top = (gapped ? rank_gapped++ : rank_tag++) < min_protein_hits_;
        if (use_protein_fdr && ! top && passing_proteins.count((Size)(int)protein_hits_[i].getMetaValue(mi_fasta_index)) == 0) continue;
        hits.push_back(protein_hits_[i]);
      }

      std::vector<int> spec_vec;
      spec_vec.reserve(dspec.size() + 1);
      spec_vec.push_back(0);
      for (const auto& pg : dspec)
      {
        spec_vec.push_back(int(round(pg.getMonoMass())));
      }

      // Terminal ladder candidates are extended on their own, with their ladder tags only. A ladder of a few masses reads as an
      // undetermined proteoform on any protein: mixed with the others, such PrSMs (2226 decoys vs 8 without them on Eclip_1480) raised the
      // score threshold of the undetermined proteoforms from 2 to 33, and ladder tags added to the sequence tag candidates still raised it
      // to 22. So a ladder PrSM of determined mass and termini competes with the others (the better PrSM is kept as one extension would
      // keep it: the best score, the earlier candidate among equal ones), and an undetermined one is reported only for a spectrum without
      // any other PrSM, in its own FDR class.
      std::vector<ProteinHit> ladder_hits(ladder_candidates.begin() + ladder_ranges[d].first, ladder_candidates.begin() + ladder_ranges[d].second);
      auto extend = [&](std::vector<ProteinHit>& candidates, std::vector<ProteinHit>& proteoforms) {
        proteoforms.clear();
        if (candidates.empty()) return;
        FLASHExtenderAlgorithm extender;
        extender.setParameters(eparam);
        extender.setCandidateBlindModificationMap(blind_mod_map);
        // with raw-spectrum scoring, every candidate's PrSM is kept here and the best one per spectrum is chosen by the combined score in
        // annotateOnRawSpectra_: the extension score alone (deconvolved masses) picked another protein in 1 of 15 spectra of Eclip_1480
        // where TopPIC's protein had the better raw evidence (offline: 5252 -> 5451 PrSMs, 1435 -> 1500 proteoforms at 1%)
        extender.run(candidates, dspec, spec_vec, vec_pro, rev_vec_pro, tags_, tol, multiple_hits_per_spec_ || keep_all_candidates);
        extender.fillProteoforms(proteoforms);
      };
      extend(hits, out);
      auto is_determined = [](const ProteinHit& h) {
        return (double)h.getMetaValue("Mass") > 0 && (int)h.getMetaValue("StartPosition") > 0 && (int)h.getMetaValue("EndPosition") > 0;
      };
      if (! ladder_hits.empty())
      {
        std::vector<ProteinHit> ladder_out;
        extend(ladder_hits, ladder_out);
        std::vector<ProteinHit> determined, undetermined;
        for (auto& h : ladder_out)
        {
          h.setMetaValue(mi_terminal, 1);
          (is_determined(h) ? determined : undetermined).push_back(h);
        }
        const bool fill = out.empty();
        if (multiple_hits_per_spec_ || keep_all_candidates) { out.insert(out.end(), determined.begin(), determined.end()); }
        else if (! determined.empty() && (out.empty() || ! is_determined(out[0]) || determined[0].getScore() > out[0].getScore()))
        {
          // a ladder PrSM that explains the precursor mass also replaces an undetermined one of a higher score
          out.swap(determined);
        }
        if (fill && out.empty()) out.swap(undetermined);
        else if (fill && (multiple_hits_per_spec_ || keep_all_candidates)) { out.insert(out.end(), undetermined.begin(), undetermined.end()); }
      }
      // Precursor-ladder candidates: extended on their own, as the terminal-ladder ones. Only their PrSMs of determined mass and termini
      // are kept (the same rule for targets and decoys): the source proposes a window that explains the precursor mass, and an
      // undetermined PrSM from it means the extension could not; it would only add chance PrSMs to a class meant to be determined.
      // A kept PrSM that is the same proteoform (same accession, areConsistent_) as a determined PrSM another source found in this spectrum
      // is dropped too: the proteoform stays with that source and its class, so the precursor-ladder class holds what this source adds,
      // and correct duplicates neither pad the class meant to show this source's own chance rate nor count twice where several PrSMs per
      // spectrum are kept. The rule looks at no label, so it treats targets and decoys alike.
      // With one PrSM per spectrum, a precursor-ladder PrSM replaces the spectrum's PrSM only with a higher score. A terminal-ladder PrSM
      // also replaces an undetermined one of a higher score, because an intact-terminus ladder is strong evidence; a window with a few
      // matched positions among the windows that fit a precursor mass (about 4e4 per spectrum on a human-size database) is not, and such
      // a rule would replace most undetermined
      // PrSMs of the spectra with a precursor mass by chance windows (harness: rank-1 windows in 92 of 190 E. coli spectra, 26 of them
      // decoys). Being last, they lose the ties (here, and in the best PrSM per spectrum by the combined score: determined, then coverage,
      // then the earlier PrSM), so a PrSM is in their class only if their source gave the better one.
      if (pl_first != pl_second)
      {
        const auto& pl_candidates_now = decoy_mass ? decoy_pl_candidates : pl_candidates;
        std::vector<ProteinHit> pl_hits(pl_candidates_now.begin() + pl_first, pl_candidates_now.begin() + pl_second), pl_out, pl_kept;
        extend(pl_hits, pl_out);
        Size undetermined = 0, duplicates = 0, at_both_ends = 0, at_start = 0, at_end = 0;
        for (auto& h : pl_out)
        {
          h.setMetaValue(mi_precursor_ladder, 1);
          if (! is_determined(h))
          {
            undetermined++;
            continue;
          }
          const bool duplicate = std::any_of(out.begin(), out.end(), [&](const ProteinHit& o) {
            return is_determined(o) && o.getAccession() == h.getAccession() && areConsistent_(o, h, pl_precursor_tol);
          });
          if (duplicate)
          {
            duplicates++;
            continue;
          }
          const bool start_ok = (int)h.getMetaValue("StartPosition") == (int)h.getMetaValue(mi_pl_start);
          const bool end_ok = (int)h.getMetaValue("EndPosition") == (int)h.getMetaValue(mi_pl_end);
          at_both_ends += start_ok && end_ok;
          at_start += start_ok && ! end_ok;
          at_end += end_ok && ! start_ok;
          pl_kept.push_back(std::move(h));
        }
        pl_determined += pl_kept.size();
        pl_at_both_ends += at_both_ends;
        pl_at_start += at_start;
        pl_at_end += at_end;
        pl_duplicates += duplicates;
        pl_undetermined += undetermined;
        if (multiple_hits_per_spec_ || keep_all_candidates) { out.insert(out.end(), pl_kept.begin(), pl_kept.end()); }
        else if (! pl_kept.empty() && (out.empty() || pl_kept[0].getScore() > out[0].getScore())) { out.swap(pl_kept); }
      }
    };
    // Extends the given spectra with the given parameters. Each spectrum is extended on its own and writes only its own slot, so the result
    // depends neither on the thread count and schedule nor on which other spectra are extended in the same call.
    auto extendSpectra = [&](const std::vector<int>& spectra, const Param& eparam) {
      const bool over_spectra = overSpectra(spectra);
#pragma omp parallel for schedule(dynamic) if (over_spectra)
      for (int k = 0; k < (int)spectra.size(); k++)
      {
        nextProgress();
        const int d = spectra[k];
        extendSpectrumInto(d, eparam, dspecs[d], spec_proteoforms[d], false);
      }
    };

    // Variable modification penalties from this run: a decoy (shuffled) protein can use a modification only by chance, so a modification
    // that target PrSMs use per opportunity (a fitting residue in the proteoform, or the proteoform terminus) no more often than decoy PrSMs
    // explains nothing and gets a penalty up to that of a blind shift. The usage is counted in a pilot pass: every var_mod_pilot_stride-th
    // spectrum in input order, extended with the first-pass penalties. If the pilot has fewer than var_mod_min_decoy_prsms decoy PrSMs, the
    // other spectra are extended with the first-pass penalties too and the usage is counted over all spectra (a small run gains nothing
    // from a pilot). Then the main pass makes the result exactly one extension of every spectrum at the penalties that are used:
    // - a raised penalty: every spectrum is extended with the learned penalties, the pilot spectra again. Re-extending only the spectra
    //   whose kept PrSMs used a raised modification is not the same: a candidate whose best path used one can have been dropped afterwards
    //   (e.g. the terminal- and precursor-ladder rules drop undetermined PrSMs and the same proteoform as another source's) while, at the
    //   raised penalty, another path of it is kept (Eclip_1480: a better unmodified 1% PrSM missed in 25 spectra, CISY lost, a false
    //   RS11_BACTN reported).
    // - no raised penalty (or none learned): only the spectra not extended yet are extended, with the first-pass penalties. A spectrum's
    //   extension depends only on the spectrum, its candidates and the parameters, none of which changes between the passes, so the pilot
    //   spectra keep exactly the PrSMs a single pass gives them.
    // Without learning, the main pass is the only one.
    Param learned_extender_param; // empty unless penalties were raised
    FLASHExtenderAlgorithm probe;
    probe.setParameters(extender_param_);
    const auto var_names = probe.getVariableModificationNames();
    const auto var_residues = probe.getVariableModificationResidues();
    const auto var_positions = probe.getVariableModificationPositions();
    std::vector<int> not_extended = all_spectra; // the spectra without PrSMs at the first-pass penalties yet
    if (learn_var_mod_penalty_ && ! var_names.empty() && decoy_factor_ > 0)
    {
      std::vector<int> pilot;
      not_extended.clear();
      for (int d : all_spectra)
      {
        (d % var_mod_pilot_stride == 0 ? pilot : not_extended).push_back(d);
      }
      startProgress(0, (SignedSize)pilot.size(), "Running extension algorithm (pilot pass for the variable modification penalties) ...");
      extendSpectra(pilot, extender_param_);
      endProgress();

      const Size n_mods = var_names.size();
      std::map<String, Size> name_to_mod;
      for (Size i = n_mods; i-- > 0;) name_to_mod[var_names[i]] = i;
      const UInt mi_variable_modifications = registry.registerName("VariableModifications");
      // the modifications a PrSM used (every alternative named for a shift counts)
      auto usedMods = [&](const ProteinHit& h, std::vector<char>& used) {
        used.assign(n_mods, 0);
        if (! h.metaValueExists(mi_variable_modifications)) return;
        for (const auto& entry : h.getMetaValue(mi_variable_modifications).toStringList())
        {
          std::vector<String> alternatives;
          String(entry).split(',', alternatives);
          for (const auto& alt : alternatives)
          {
            std::vector<String> names;
            alt.split('+', names);
            for (const auto& nm : names)
            {
              auto it = name_to_mod.find(nm);
              if (it != name_to_mod.end()) used[it->second] = 1;
            }
          }
        }
      };
      // uses and opportunities per modification, targets and decoys, over the PrSMs of the given spectra (in spectrum order, so the sums
      // are the same in every run). Precursor-ladder PrSMs take no part: their source proposes only windows whose unmodified mass is the
      // precursor mass, so they are unmodified by selection, not by chance, and their opportunities without uses would move the
      // target/decoy enrichment by the share of targets and decoys among them
      std::vector<double> use_t, use_d, opp_t, opp_d;
      Size n_d = 0;
      auto countUsage = [&](const std::vector<int>& spectra) {
        use_t.assign(n_mods, 0);
        use_d.assign(n_mods, 0);
        opp_t.assign(n_mods, 0);
        opp_d.assign(n_mods, 0);
        n_d = 0;
        std::vector<char> used;
        for (int d : spectra)
        {
          for (const auto& h : spec_proteoforms[d])
          {
            if (h.metaValueExists(mi_precursor_ladder)) continue;
            const bool decoy = h.getAccession().hasPrefix("DECOY");
            const String& seq = h.getSequence();
            int s = h.metaValueExists("StartPosition") ? (int)h.getMetaValue("StartPosition") - 1 : 0;
            int e = h.metaValueExists("EndPosition") ? (int)h.getMetaValue("EndPosition") : (int)seq.size();
            if (s < 0) s = 0;
            if (e <= s || e > (int)seq.size()) e = (int)seq.size();
            n_d += decoy;
            usedMods(h, used);
            for (Size m = 0; m < n_mods; m++)
            {
              double opportunities = 0;
              if (var_positions[m] == 1) opportunities = var_residues[m].empty() || var_residues[m].has(seq[s]);
              else if (var_positions[m] == 2) opportunities = var_residues[m].empty() || var_residues[m].has(seq[e - 1]);
              else
              {
                for (int p = s; p < e; p++) opportunities += var_residues[m].empty() || var_residues[m].has(seq[p]);
              }
              (decoy ? opp_d : opp_t)[m] += opportunities;
              if (used[m]) (decoy ? use_d : use_t)[m] += 1;
            }
          }
        }
      };
      countUsage(pilot);
      const std::vector<int>* learned_from = &pilot;
      if (n_d < var_mod_min_decoy_prsms)
      {
        OPENMS_LOG_INFO << "Variable modification penalties: " << n_d << " decoy PrSMs in the pilot pass (fewer than " << var_mod_min_decoy_prsms
                        << "), so the other spectra are extended with the first-pass penalties too and all spectra are used" << std::endl;
        startProgress(0, (SignedSize)not_extended.size(), "Running extension algorithm (the other spectra for the variable modification penalties) ...");
        extendSpectra(not_extended, extender_param_);
        endProgress();
        not_extended.clear();
        countUsage(all_spectra);
        learned_from = &all_spectra;
      }
      // enough decoy PrSMs to estimate chance usage
      if (n_d >= var_mod_min_decoy_prsms)
      {
        OPENMS_LOG_INFO << "Variable modification penalties from the PrSMs of " << learned_from->size() << " of " << dspecs.size() << " spectra ("
                        << n_d << " decoy PrSMs)" << std::endl;
        std::vector<int> penalties(n_mods, 2);
        bool raised = false;
        for (Size m = 0; m < n_mods; m++)
        {
          // usage per opportunity, target over decoy, with one pseudo-use on each side; 8x enrichment or more keeps the default
          // penalty 2, none gives 19 (a blind shift)
          if (opp_t[m] <= 0 || opp_d[m] <= 0) continue;
          const double enrichment = ((use_t[m] + 1) / opp_t[m]) / ((use_d[m] + 1) / opp_d[m]);
          penalties[m] = std::max(2, std::min(19, (int)std::lround(19.0 - 17.0 / 3.0 * std::log2(std::max(1.0, enrichment)))));
          raised |= penalties[m] > 2;
          OPENMS_LOG_INFO << "Variable modification " << var_names[m] << ": used by " << use_t[m] << " target and " << use_d[m]
                          << " decoy PrSMs, enrichment per opportunity " << enrichment << ", penalty " << penalties[m] << std::endl;
        }
        if (raised)
        {
          learned_extender_param = extender_param_;
          learned_extender_param.setValue("var_mod_penalty", penalties);
        }
      }
      else
      {
        OPENMS_LOG_INFO << "Variable modification penalties not learned: " << n_d << " decoy PrSMs in all spectra (fewer than "
                        << var_mod_min_decoy_prsms << ")" << std::endl;
      }
    }

    // Main pass: with learned penalties every spectrum, otherwise the spectra not extended yet with the first-pass penalties. The
    // precursor-ladder counts describe the extension that is kept: they restart before a pass that replaces the pilot's PrSMs, and otherwise
    // the counts of the passes add up to those of a single pass.
    const bool penalties_learned = ! learned_extender_param.empty();
    if (penalties_learned)
    {
      OPENMS_LOG_INFO << "Extending all " << dspecs.size() << " spectra with the learned penalties" << std::endl;
      pl_determined = 0;
      pl_at_both_ends = 0;
      pl_at_start = 0;
      pl_at_end = 0;
      pl_duplicates = 0;
      pl_undetermined = 0;
    }
    const std::vector<int>& main_spectra = penalties_learned ? all_spectra : not_extended;
    startProgress(0, (SignedSize)main_spectra.size(), "Running extension algorithm ...");
    extendSpectra(main_spectra, penalties_learned ? learned_extender_param : extender_param_);
    if (use_precursor_ladders)
    {
      OPENMS_LOG_INFO << "Precursor-ladder PrSMs: " << pl_determined.load() << " of determined mass and termini kept (" << pl_at_both_ends.load()
                      << " at both window ends, " << pl_at_start.load() << " at the window start only, " << pl_at_end.load()
                      << " at the window end only), " << pl_duplicates.load() << " dropped as the same proteoform as a determined PrSM of another source, "
                      << pl_undetermined.load() << " undetermined dropped" << std::endl;
    }
    // Decoy precursors: a spectrum whose PrSM took its mass from MS1 is searched again as a decoy spectrum whose precursor mass is the
    // decoy precursor mass (noise or signal decoy) of its isolation window, as FLASHDeconv chose it with report_FDR. A PrSM of that copy
    // that reports the decoy mass (not a mass from fragment pairs) is a decoy of the precursor: it competes with the spectrum's MS1-based
    // PrSMs (the higher score stays) and counts as a decoy in the FDR. PrSMs with masses from fragment pairs or without a mass (-1) are
    // not affected.
    if (decoy_precursor_)
    {
      const double iso = Constants::C13C12_MASSDIFF_U;
      auto ms1Based = [](const ProteinHit& h) {
        return (double)h.getMetaValue("Mass") > 0 && (int)h.getMetaValue("ProteoformMassByFragmentMass") == 0;
      };
      const Param& eparam_now = learned_extender_param.empty() ? extender_param_ : learned_extender_param;
      std::vector<std::vector<ProteinHit>> decoy_proteoforms(dspecs.size());
      int copies = 0, without_decoy = 0;
      std::vector<int> todo;
      for (int d = 0; d < (int)dspecs.size(); d++)
      {
        if (dspecs[d].getPrecursorPeakGroup().getMonoMass() <= 0) continue;
        if (std::none_of(spec_proteoforms[d].begin(), spec_proteoforms[d].end(), ms1Based)) continue;
        if (dspecs[d].getDecoyPrecursorMass() <= 0)
        {
          without_decoy++;
          continue;
        }
        todo.push_back(d);
      }
      copies = (int)todo.size();
      if (copies == 0)
      {
        OPENMS_LOG_WARN << "Decoy precursors: the input has no decoy precursor masses; run FLASHDeconv with -FD:report_FDR true" << std::endl;
      }
#pragma omp parallel for schedule(dynamic) if (extend_over_spectra)
      for (int k = 0; k < (int)todo.size(); k++)
      {
        const int d = todo[k];
        DeconvolvedSpectrum decoy_spec = dspecs[d];
        auto pg = decoy_spec.getPrecursorPeakGroup();
        const double decoy_mass = decoy_spec.getDecoyPrecursorMass();
        pg.setMonoisotopicMass(decoy_mass);
        decoy_spec.setPrecursorPeakGroup(pg);
        std::vector<ProteinHit> out;
        extendSpectrumInto(d, eparam_now, decoy_spec, out, true);
        for (auto& h : out)
        {
          if (! ms1Based(h)) continue;
          if (std::abs((double)h.getMetaValue("Mass") - decoy_mass) > 3 * iso + 0.1) continue; // the moved mass, up to isotope errors
          h.setMetaValue("DecoyPrecursor", decoy_spec.getDecoyPrecursorWeight()); // its weight as a decoy (FLASHDeconv's noise decoy weight)
          decoy_proteoforms[d].push_back(std::move(h));
        }
      }
      // per spectrum, the decoy precursor competes with the spectrum's PrSMs that took their mass from MS1: the higher score stays
      int decoy_wins = 0;
      for (int d : todo)
      {
        if (decoy_proteoforms[d].empty()) continue;
        double best_decoy = -1e10, best_own = -1e10;
        for (const auto& h : decoy_proteoforms[d]) best_decoy = std::max(best_decoy, h.getScore());
        for (const auto& h : spec_proteoforms[d])
        {
          if (ms1Based(h)) best_own = std::max(best_own, h.getScore());
        }
        if (best_decoy <= best_own) continue;
        decoy_wins++;
        std::vector<ProteinHit> kept;
        for (auto& h : spec_proteoforms[d])
        {
          if (! ms1Based(h)) kept.push_back(std::move(h));
        }
        for (auto& h : decoy_proteoforms[d])
        {
          if (h.getScore() == best_decoy) kept.push_back(std::move(h));
        }
        spec_proteoforms[d].swap(kept);
      }
      OPENMS_LOG_INFO << "Decoy precursors: " << copies << " spectra with an MS1-based PrSM and an MS1 decoy mass (" << without_decoy
                      << " without one), decoy precursor wins in " << decoy_wins << std::endl;
    }

    for (auto& v : spec_proteoforms)
    {
      for (auto& h : v) proteoform_hits_.push_back(std::move(h));
    }
    endProgress();
  }

  // raw-spectrum evidence of every PrSM, targets and decoys alike, before the FDR
  if (raw_loader_ && ! proteoform_hits_.empty()) annotateOnRawSpectra_();

  // the order below and every ranking of the FDR need finite scores: a NaN breaks the sorts and never closes a tie group, and an infinite
  // decoy score leaves the null of the protein level undefined
  for (const auto& hit : proteoform_hits_)
  {
    if (std::isfinite(hit.getScore())) continue;
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "PrSM score must be finite (scan " + hit.getMetaValue("Scan").toString() + ")", std::to_string(hit.getScore()));
  }
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), [](const ProteinHit& left, const ProteinHit& right) {
    return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getMetaValue("Scan") > right.getMetaValue("Scan"))
                                                                                            : (left.getCoverage() > right.getCoverage()))
                                               : (left.getScore() > right.getScore());
  });

  // FDR by target-decoy counting with a +1 pseudo-count (targetDecoyQvalues) at three levels:
  // - PrSMs, per class (fdr_class_names): determined mass and termini with 0, 1, or 2 or more blind modifications (variable modifications
  //   do not count), and undetermined with both termini known (the mass is then open only inside) or with an open terminus. Their decoy
  //   rates differ widely (Eclip_1480 with decoys kept: 2 decoys for 3770 determined targets, 97 for 1070 undetermined with both
  //   termini, 1030 for 2228 with an open terminus). The PrSMs of precursor-ladder candidates (meta value PrecursorLadder; only
  //   determined ones that no other source found in the spectrum are kept) are a class of their own, assigned by that source before the
  //   rules above: a source that searches every window fitting the precursor mass has its own chance rate, which the nearly decoy-free
  //   determined classes would hide. The class does not split by blind modifications: a window fits the unmodified mass, so such a PrSM
  //   with blind shifts left its window (the log line of the class counts them, targets and decoys, to decide whether it should be split).
  //   At the protein level the class uses the pooled null while it has fewer than null_min_class_decoys decoys, and its decoys join the
  //   pooled null (so the null of the sparse determined classes changes with the option; not yet measured). Each PrSM is in exactly one
  //   class, so no PrSM passes twice. With discard_underdetermined the undetermined PrSMs are in no class and take part in nothing below,
  //   the proteoforms included.
  // - proteoforms: all representatives in one ranking by their PrSM q-value, which carries the class, then by score;
  // - proteins: setProteinQvalues_.
  // Without decoys in the database and without discard_underdetermined no q-value is set.
  const bool use_fdr = decoy_factor_ > 0 || ! keep_underdetermined_;
  auto is_decoy = [](const ProteinHit& hit) { return hit.getAccession().hasPrefix("DECOY") || hit.metaValueExists("DecoyPrecursor"); };
  // a decoy precursor counts with its weight (1 for a signal decoy, FLASHDeconv's noise decoy weight for a noise decoy), a protein decoy
  // for 1 / decoy_factor_
  auto decoy_weight = [this](const ProteinHit& hit) {
    return hit.metaValueExists("DecoyPrecursor") && ! hit.getAccession().hasPrefix("DECOY") ? (double)hit.getMetaValue("DecoyPrecursor") : 1.0 / decoy_factor_;
  };
  // PrSM q-value of a passing PrSM for the representatives and the protein ranking: the reported level, or 1% without a filter
  // (pro_fdr 1, as in keep_decoy runs for offline evaluation), so that such a run has the proteoforms of a 1% run
  const double rep_q = pro_fdr_ < 1 ? pro_fdr_ : 0.01;
  const Size hit_count = proteoform_hits_.size();
  std::vector<int> hit_class(hit_count, -1); // index into fdr_class_names; -1: in no class
  std::vector<double> qvalues(hit_count, 1.0), decoy_weights(hit_count, 0.0);
  std::vector<bool> decoys(hit_count, false);
  if (use_fdr)
  {
    std::vector<std::vector<Size>> class_hits(fdr_class_count); // in the order above
    for (Size i = 0; i < hit_count; i++)
    {
      const auto& hit = proteoform_hits_[i];
      decoys[i] = is_decoy(hit);
      if (decoys[i]) decoy_weights[i] = decoy_weight(hit);
      int c;
      if (hit.metaValueExists(mi_precursor_ladder)) c = fdr_class_precursor_ladder; // by source; all of them are determined
      else if (isDetermined(hit)) c = (int)std::min(blindModCount_(hit), (Size)fdr_class_max_blind);
      else if (! keep_underdetermined_) continue;
      else
      {
        const bool both_termini = (int)hit.getMetaValue("StartPosition") > 0 && (int)hit.getMetaValue("EndPosition") > 0;
        c = both_termini ? fdr_class_both_termini : fdr_class_open_terminus;
      }
      hit_class[i] = c;
      class_hits[c].push_back(i);
    }
    // per class: targets and decoys (count and weight), and the protein-level null it uses (its own from null_min_class_decoys decoys on,
    // else the pooled one); for the precursor-ladder class also by blind modification count
    for (int c = 0; c < fdr_class_count; c++)
    {
      Size targets = 0, decoy_count = 0, blind_targets[3] = {0, 0, 0}, blind_decoys[3] = {0, 0, 0};
      double weight = 0;
      for (Size i : class_hits[c])
      {
        const Size b = std::min(blindModCount_(proteoform_hits_[i]), (Size)2);
        if (decoys[i])
        {
          decoy_count++;
          weight += decoy_weights[i];
          blind_decoys[b]++;
        }
        else
        {
          targets++;
          blind_targets[b]++;
        }
      }
      String line = String("FDR class \"") + fdr_class_names[c] + "\": " + String(targets) + " target and " + String(decoy_count)
                    + " decoy PrSMs (decoy weight " + String(weight) + "); protein-level null: "
                    + (decoy_count >= null_min_class_decoys ? "own" : "pooled");
      if (c == fdr_class_precursor_ladder && ! class_hits[c].empty())
      {
        line += "; target/decoy PrSMs by blind modifications: 0: " + String(blind_targets[0]) + "/" + String(blind_decoys[0]) + ", 1: "
                + String(blind_targets[1]) + "/" + String(blind_decoys[1]) + ", 2 or more: " + String(blind_targets[2]) + "/" + String(blind_decoys[2]);
      }
      OPENMS_LOG_INFO << line << std::endl;
    }
    for (int c = 0; c < fdr_class_count; c++)
    {
      std::vector<RankedItem<double>> ranked;
      ranked.reserve(class_hits[c].size());
      for (Size i : class_hits[c]) ranked.push_back({proteoform_hits_[i].getScore(), decoys[i], decoy_weights[i]});
      const auto class_qvalues = targetDecoyQvalues(ranked);
      for (Size k = 0; k < class_hits[c].size(); k++)
      {
        const Size i = class_hits[c][k];
        qvalues[i] = class_qvalues[k];
        proteoform_hits_[i].setMetaValue("qvalue", qvalues[i]);
        proteoform_hits_[i].setMetaValue("FDRClass", fdr_class_names[c]);
      }
    }
  }
  // the proteoforms: positions of their PrSMs in proteoform_hits_ (Index), the representative first; representatives in the order above
  const auto clusters = markRepresentativeProteoformHits_(precursor_tol, hit_class, qvalues, rep_q, use_fdr);

  if (use_fdr)
  {
    // Proteoform q-values: all representatives (targets and decoys of every class) in one ranking. Every PrSM of a proteoform carries
    // the proteoform's q-value, which can be lower than its own PrSM q-value.
    {
      std::vector<Size> order(clusters.size());
      std::iota(order.begin(), order.end(), 0);
      std::stable_sort(order.begin(), order.end(), [&](Size l, Size r) {
        const Size rep_l = clusters[l][0], rep_r = clusters[r][0];
        if (qvalues[rep_l] != qvalues[rep_r]) return qvalues[rep_l] < qvalues[rep_r];
        return proteoform_hits_[rep_l].getScore() > proteoform_hits_[rep_r].getScore();
      });
      std::vector<RankedItem<std::pair<double, double>>> ranked;
      ranked.reserve(order.size());
      for (Size k : order)
      {
        const Size rep = clusters[k][0];
        ranked.push_back({{qvalues[rep], -proteoform_hits_[rep].getScore()}, decoys[rep], decoy_weights[rep]});
      }
      const auto proteoform_qvalues = targetDecoyQvalues(ranked);
      for (Size k = 0; k < order.size(); k++)
      {
        for (Size i : clusters[order[k]]) proteoform_hits_[i].setMetaValue("proqvalue", proteoform_qvalues[k]);
      }
    }
    setProteinQvalues_(hit_class, qvalues, decoys, decoy_weights, clusters, rep_q);
    std::vector<ProteinHit> filtered_proteoform_hits;
    filtered_proteoform_hits.reserve(proteoform_hits_.size());
    for (Size i = 0; i < hit_count; i++)
    {
      if (hit_class[i] < 0) continue; // in no class, without a q-value (undetermined, with keep_underdetermined_ false)
      if (! keep_decoy_ && decoys[i]) continue;
      // no PrSM-level filter: a PrSM is reported when its proteoform passes (every member carries its cluster's proteoform q-value)
      if (pro_fdr_ < 1 && (double)proteoform_hits_[i].getMetaValue("proqvalue", 1.0) > pro_fdr_) continue;
      filtered_proteoform_hits.push_back(proteoform_hits_[i]);
    }
    proteoform_hits_.swap(filtered_proteoform_hits);
  }
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(),
            [](const ProteinHit& left, const ProteinHit& right) { return left.getMetaValue("RT") < right.getMetaValue("RT"); });


  //  define prsm index and tag to prsm indices.
  int prsm_index = 0;
  std::map<int, int> prsm_index_map;

  for (auto& hit : proteoform_hits_)
  {
    int fasta_index = hit.getMetaValue("FastaIndex");
    hit.setSequence(original_fasta_entry[fasta_index].sequence);

    prsm_index_map[hit.getMetaValue("Index")] = prsm_index; // from previous to new index
    hit.setMetaValue("Index", prsm_index);
    for (int tag_index : (std::vector<int>)hit.getMetaValue("TagIndices").toIntList())
    {
      matching_hits_indices_[tag_index].push_back(prsm_index);
    }
    prsm_index++;
  }

  // redefine prsm indices for representative proteoforms
  for (auto& hit : proteoform_hits_)
  {
    if ((int)hit.getMetaValue("Representative") == 0) continue;
    auto prsm_indices = hit.getMetaValue("PrSMIndices").toIntList();
    std::vector<int> new_prsm_indices;
    for (int i : prsm_indices)
    {
      // PrSMs of the proteoform that did not pass the PrSM-level FDR are not reported (operator[] made them PrSM 0 of the output)
      const auto it = prsm_index_map.find(i);
      if (it != prsm_index_map.end()) new_prsm_indices.push_back(it->second);
    }
    std::sort(new_prsm_indices.begin(), new_prsm_indices.end());
    hit.setMetaValue("PrSMIndices", new_prsm_indices);
  }
  // the index of the reported proteoform (as FLASHTnTFile::writeProteoforms numbers them) on every PrSM of it; -1: not reported
  {
    const double pro_fdr = param_.getValue("pro_fdr");
    std::map<int, int> prsm_to_proteoform;
    int proteoform_index = 0;
    for (const auto& hit : proteoform_hits_)
    {
      if ((int)hit.getMetaValue("Representative") == 0) continue;
      if (hit.metaValueExists("proqvalue") && (double)hit.getMetaValue("proqvalue") > pro_fdr) continue;
      for (int i : hit.getMetaValue("PrSMIndices").toIntList()) prsm_to_proteoform[i] = proteoform_index;
      prsm_to_proteoform[(int)hit.getMetaValue("Index")] = proteoform_index;
      proteoform_index++;
    }
    for (auto& hit : proteoform_hits_)
    {
      const auto it = prsm_to_proteoform.find((int)hit.getMetaValue("Index"));
      hit.setMetaValue("ProteoformIndex", it == prsm_to_proteoform.end() ? -1 : it->second);
    }
  }
  // TODO
  // per scan, get the scan number and the precursor mass - maybe given or maybe deconvolved.
  // get all protoeforms - mass and score.
  // get all combinations with a certain tolerance.
  // record the best ranking combination(s) along with summed score (can score be non overlapping?) and mass difference
}

void FLASHTnTAlgorithm::annotateOnRawSpectra_()
{
  // The extension score counts deconvolved masses on the path; the raw score (FLASHAnnotatorAlgorithm) is the chance of the fragment isotope
  // envelopes found in the raw spectrum, at every charge, with the chance rate measured on mass-shifted decoy fragments of the same spectrum.
  // Decoy PrSMs are annotated the same way, so the FDR of the combined score stays target-decoy controlled. Eclip_1480: at 1%, 4846 PrSMs
  // by the extension score, 5189 by the raw score alone, 5241 by extension + 4 x raw, 5247 at 8 x raw (5240 at 16, 5012 at 1); the
  // three lysates: 4278 proteoforms / 991 proteins at 4, 4293 / 1002 at 8.
  std::set<int> scans;
  for (const auto& hit : proteoform_hits_) scans.insert((int)hit.getMetaValue("Scan"));
  std::map<int, MSSpectrum> raw;
  raw_loader_(scans, raw);
  std::vector<int> scan_list;
  std::vector<MSSpectrum*> raw_list;
  std::map<int, int> slot;
  for (auto& [scan, spec] : raw)
  {
    slot[scan] = (int)scan_list.size();
    scan_list.push_back(scan);
    raw_list.push_back(&spec);
  }
  std::vector<FLASHAnnotatorAlgorithm::PeakArrays> spectra(scan_list.size());
  int profile_count = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : profile_count)
  for (int i = 0; i < (int)scan_list.size(); i++)
  {
    profile_count += FLASHAnnotatorAlgorithm::toCentroidedPeakArrays(*raw_list[i], spectra[i]) ? 1 : 0;
    raw_list[i]->clear(false); // the precursor information stays
    raw_list[i]->shrink_to_fit();
  }

  FLASHAnnotatorAlgorithm annotator;
  annotator.setParameters(annotator_param_);
  std::vector<FLASHAnnotatorAlgorithm::PrSM> prsms(proteoform_hits_.size());
  std::vector<double> precursor_masses(proteoform_hits_.size(), -1.0);
  double max_mass = 0;
  Size without_spectrum = 0;
  for (Size i = 0; i < proteoform_hits_.size(); i++)
  {
    const auto& hit = proteoform_hits_[i];
    auto& prsm = prsms[i];
    const String& seq = hit.getSequence();
    const int start = hit.getMetaValue("StartPosition"), end = hit.getMetaValue("EndPosition");
    const int start_in_seq = start <= 0 ? 0 : std::min((int)seq.size(), start - 1);
    const int end_in_seq = end <= 0 ? (int)seq.size() : std::max(start_in_seq, std::min((int)seq.size(), end));
    prsm.sequence = seq.substr(start_in_seq, end_in_seq - start_in_seq);
    prsm.known_n_term = start > 0;
    prsm.known_c_term = end > 0;
    if (hit.metaValueExists("Modifications"))
    {
      const std::vector<double> masses = hit.getMetaValue("Modifications");
      const std::vector<int> starts = hit.getMetaValue("ModificationStarts"), ends = hit.getMetaValue("ModificationEnds");
      for (Size m = 0; m < masses.size() && m < starts.size() && m < ends.size(); m++)
      {
        FLASHAnnotatorAlgorithm::Modification mod;
        mod.mass = masses[m];
        mod.start = starts[m] - start_in_seq;
        mod.end = ends[m] - start_in_seq;
        prsm.modifications.push_back(mod);
      }
    }
    precursor_masses[i] = hit.metaValueExists("GivenMass") ? (double)hit.getMetaValue("GivenMass") : -1.0;
    max_mass = std::max({max_mass, precursor_masses[i], FLASHAnnotatorAlgorithm::getProteoformMass(prsm.sequence, prsm.modifications)});
    const auto it = slot.find((int)hit.getMetaValue("Scan"));
    if (it == slot.end()) without_spectrum++;
    else prsm.spectrum_index = it->second;
  }
  // the isotope table once, beyond the heaviest PrSM; then the precursor charge caps
  std::vector<double> iso_mzs, iso_intensities;
  annotator.getIsotopeEnvelope(max_mass + 200.0, 1, iso_mzs, iso_intensities);
  for (Size i = 0; i < prsms.size(); i++)
  {
    if (prsms[i].spectrum_index < 0) continue;
    int reason = 4;
    prsms[i].precursor_charge = annotator.precursorChargeCap(*raw_list[prsms[i].spectrum_index], precursor_masses[i], reason);
  }
  std::vector<FLASHAnnotatorAlgorithm::Result> results;
  annotator.annotate(spectra, prsms, results);

  Size annotated = 0;
  for (Size i = 0; i < proteoform_hits_.size(); i++)
  {
    auto& hit = proteoform_hits_[i];
    const auto& res = results[i];
    annotated += res.annotated ? 1 : 0;
    hit.setMetaValue("ExtensionScore", hit.getScore());
    hit.setMetaValue("RawScore", res.score);
    hit.setMetaValue("RawMatchedFragments", res.matched_fragments);
    hit.setMetaValue("RawFragmentHypotheses", res.fragment_hypotheses);
    hit.setMetaValue("RawCoverage", res.coverage);
    hit.setScore(hit.getScore() + raw_score_weight_ * std::max(0.0, res.score));
  }
  OPENMS_LOG_INFO << "Raw spectra: " << spectra.size() << " loaded (" << profile_count << " centroided); PrSMs annotated: " << annotated << " of "
                  << proteoform_hits_.size() << " (" << without_spectrum << " without their scan in the raw file)" << std::endl;

  // one PrSM per spectrum: the best by the combined score (a determined mass and termini break a tie, then the coverage); the
  // extension kept every candidate's PrSM for this choice
  if (! multiple_hits_per_spec_)
  {
    std::map<int, Size> best;
    auto is_determined = [](const ProteinHit& h) {
      return (double)h.getMetaValue("Mass") > 0 && (int)h.getMetaValue("StartPosition") > 0 && (int)h.getMetaValue("EndPosition") > 0;
    };
    auto better = [&](const ProteinHit& a, const ProteinHit& b) {
      if (a.getScore() != b.getScore()) return a.getScore() > b.getScore();
      if (is_determined(a) != is_determined(b)) return is_determined(a);
      return a.getCoverage() > b.getCoverage();
    };
    for (Size i = 0; i < proteoform_hits_.size(); i++)
    {
      const int scan = proteoform_hits_[i].getMetaValue("Scan");
      auto it = best.find(scan);
      if (it == best.end()) best[scan] = i;
      else if (better(proteoform_hits_[i], proteoform_hits_[it->second])) it->second = i;
    }
    std::vector<ProteinHit> kept;
    kept.reserve(best.size());
    for (Size i = 0; i < proteoform_hits_.size(); i++)
    {
      if (best[(int)proteoform_hits_[i].getMetaValue("Scan")] == i) kept.push_back(std::move(proteoform_hits_[i]));
    }
    OPENMS_LOG_INFO << "PrSMs: " << proteoform_hits_.size() << " of all candidates, " << kept.size() << " best per spectrum by the combined score" << std::endl;
    proteoform_hits_.swap(kept);
  }
}

void FLASHTnTAlgorithm::getProteoformHitsMatchedBy(const FLASHHelperClasses::Tag& tag, std::vector<ProteinHit>& hits) const
{
  int index = tag.getIndex();

  if (index < 0 || matching_hits_indices_.find(index) == matching_hits_indices_.end()) return;

  for (auto i : matching_hits_indices_.at(index))
  {
    hits.push_back(proteoform_hits_[i]);
  }
}

void FLASHTnTAlgorithm::getTags(std::vector<FLASHHelperClasses::Tag>& tags) const
{
  for (const auto& tag : tags_)
  {
    tags.push_back(tag);
  }
}

} // namespace OpenMS