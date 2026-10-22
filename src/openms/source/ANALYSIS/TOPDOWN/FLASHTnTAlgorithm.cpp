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
#include <OpenMS/ANALYSIS/TOPDOWN/MassLadderAlgorithm.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <exception>
#include <iterator>
#include <limits>
#include <numeric>
#include <functional>
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
// so that their decoys are not hidden in the large, nearly decoy-free determined classes. So do the PrSMs of each residual round of the
// chimeric search (option chimeric_rounds; meta value ChimericRound r > 0): a search of what is left of a spectrum after its best PrSM
// has its own chance rate.
constexpr int fdr_class_max_blind = 2;
constexpr int fdr_class_both_termini = fdr_class_max_blind + 1;
constexpr int fdr_class_open_terminus = fdr_class_both_termini + 1;
constexpr int fdr_class_precursor_ladder = fdr_class_open_terminus + 1;
// the residual rounds of the chimeric search: round r is class fdr_class_residual + r - 1
constexpr int max_chimeric_rounds = 2;
constexpr int fdr_class_residual = fdr_class_precursor_ladder + 1;
constexpr int fdr_class_count = fdr_class_residual + max_chimeric_rounds;
const char* const fdr_class_names[] = {"determined; 0 blind", "determined; 1 blind", "determined; 2 blind", "undetermined, both termini",
                                       "undetermined, open terminus", "precursor ladder", "residual 1", "residual 2"};
static_assert(std::size(fdr_class_names) == (Size)fdr_class_count, "one label per FDR class");
// chimeric search: with a mass matched by the best PrSM, the masses within the spectrum's tolerance of it moved by up to this many isotope
// spacings are removed as well (FLASHDeconv's monoisotopic errors of the same fragment). On the 20 native spectra of the prototype, the
// best PrSM's fragments matched 162 and 211 masses at +1 and -1 isotope, 80 and 85 at +2 and -2, against about 6 per offset by chance
// (fragments moved by +-9, +-13, +-21 Da; final_all/chimeric/isotope_offsets.tsv); the prototype removed them up to +-2.
constexpr int chimeric_isotope_offsets = 2;
// protein level: the null of a class is exponential above this quantile of its decoy scores
constexpr double null_tail_quantile = 0.8;
// protein level: with fewer decoys the null is the empirical survival function only, without the tail
constexpr Size null_min_tail_decoys = 5;
// protein level: a class with fewer decoys uses the null of the decoys of all classes pooled
constexpr Size null_min_class_decoys = 20;
// protein level: floor of a null p-value and of a local FDR before its log
constexpr double fdr_p_floor = 1e-300;
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

// The better of two PrSMs of a spectrum, as the best PrSM per spectrum is chosen after raw-spectrum scoring: the higher score, then a
// determined mass and termini, then the higher coverage (an exact tie keeps the earlier one).
bool isBetterPrSM(const ProteinHit& a, const ProteinHit& b)
{
  auto is_determined = [](const ProteinHit& h) {
    return (double)h.getMetaValue("Mass") > 0 && (int)h.getMetaValue("StartPosition") > 0 && (int)h.getMetaValue("EndPosition") > 0;
  };
  if (a.getScore() != b.getScore()) return a.getScore() > b.getScore();
  if (is_determined(a) != is_determined(b)) return is_determined(a);
  return a.getCoverage() > b.getCoverage();
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
  defaults_.setValue("search_cascade", "false",
                     "Search the candidate sources in stages: sequence tags in every spectrum; gapped tags (if gapped_tags is not off) only in the "
                     "spectra not identified so far (no target PrSM whose proteoform passes pro_fdr, 1% without a filter); terminal and precursor "
                     "ladders (if on) only in the spectra still not identified. An identified spectrum is searched by no later stage. The PrSM q-values are computed per stage "
                     "(and FDR class), so the PrSMs of a later stage compete only with that stage's decoys. The variable modification penalties are "
                     "learned in the first stage. Not used with native, with protein_level_only or without decoys in the database; the chimeric "
                     "rounds search with all sources.");
  defaults_.setValidStrings("search_cascade", {"true", "false"});
  defaults_.setValue("one_proteoform_per_feature", "true",
                     "Co-isolation: when PrSMs of different proteins report the same mass taken from MS1 (within twice the precursor tolerance, "
                     "one isotope apart allowed) within 2 min, the mass is one MS1 feature of one proteoform, and only the protein that explains it "
                     "best keeps it (by its PrSM with the fewest blind modifications, then the higher score, then more PrSMs at that mass); the "
                     "other proteins' PrSMs are reported without a mass. Masses from fragment pairs and decoy precursors are not affected.");
  defaults_.setValidStrings("one_proteoform_per_feature", {"true", "false"});
  defaults_.setValue("native", "false",
                     "Native and targeted top-down (one or a few spectra, often several proteins per spectrum, no MS1 survey): the candidate "
                     "proteins of each source (sequence tags; sliding ladders, i.e. masses read as residue sums from any start or to any end of a "
                     "protein, instead of terminal ladders; precursor ladders) are characterized when their candidate score "
                     "(-log10 of the chance probability of their matches) exceeds that of the best decoy candidate of the same source in the same "
                     "spectrum by native_margin (at most max_protein_hits per source), and every such candidate's PrSM is reported, also without a "
                     "determined mass, one per protein (the better one). Gapped tags are not used, only_single_hit is ignored, and at least one "
                     "residual round of chimeric_rounds is searched (co-isolated proteins). Needs decoys in the database (DECOY prefix); without "
                     "them, the candidates are chosen as without the option.");
  defaults_.setValidStrings("native", {"true", "false"});
  defaults_.setValue("native_margin", 2.0,
                     "With native: the margin (-log10 units) by which a candidate's score must exceed that of the best decoy candidate of its source "
                     "in its spectrum.");
  defaults_.setMinFloat("native_margin", 0.0);
  defaults_.setValue("decoy_precursor", "false",
                     "Decoy precursors: spectra whose PrSM took its mass from MS1 are searched again with the decoy precursor mass of their isolation "
                     "window instead of the precursor mass (FLASHDeconv -FD:report_FDR true writes it and its score). When that search gives a PrSM "
                     "reporting the decoy mass, the decoy precursor wins if its precursor score (FLASHDeconv's 2D Qscore) is higher than the "
                     "precursor's; the winning decoy PrSM replaces the spectrum's MS1-based PrSMs and counts as a decoy in the FDR, at the "
                     "proteoform level only against proteoforms whose mass came from MS1. Controls wrong MS1 precursor masses; PrSMs with masses "
                     "from fragment pairs or without a mass are not affected.");
  defaults_.setValidStrings("decoy_precursor", {"true", "false"});
  defaults_.setValue("chimeric_rounds", 0,
                     "Residual search of chimeric spectra: rounds after the first search (0: off). In each round, every spectrum with a PrSM (target or "
                     "decoy alike) loses the deconvolved masses its best PrSM matched (the masses on the path, also those one or two isotope "
                     "spacings away, and its proteoform mass), the 400-mass cap is applied again to what is left, and the residual is searched "
                     "again: tags, candidate proteins (sequence and gapped tags, terminal ladders; no precursor ladders) without the proteins of "
                     "the spectrum's earlier best PrSMs and their target-decoy partners, and the extension, without the precursor mass (it belongs "
                     "to the first proteoform). Every candidate's PrSM is kept as with -in_raw (scored on the raw spectrum with it), then, with "
                     "only_single_hit, the best one. Residual PrSMs are added to the scan (meta value ChimericRound r) and form their own FDR class "
                     "per round (\"residual 1\", \"residual 2\"). The evidence in a single spectrum is reported, not filtered: ResidualNullScore, the "
                     "best extension score of the other decoy PrSMs of the same residual search, ResidualScoreRatio, the extension score over it, "
                     "and MatchedMasses, the masses on the path. Without the option, results are unchanged.");
  defaults_.setMinInt("chimeric_rounds", 0);
  defaults_.setMaxInt("chimeric_rounds", max_chimeric_rounds);
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
  search_cascade_ = param_.getValue("search_cascade").toString() == "true";
  learn_var_mod_penalty_ = param_.getValue("learn_var_mod_penalty").toString() == "true";
  decoy_precursor_ = param_.getValue("decoy_precursor").toString() == "true";
  one_proteoform_per_feature_ = param_.getValue("one_proteoform_per_feature").toString() == "true";
  chimeric_rounds_ = param_.getValue("chimeric_rounds");
  pro_fdr_ = param_.getValue("pro_fdr");
  keep_decoy_ = param_.getValue("keep_decoy").toString() == "true";
  keep_underdetermined_ = param_.getValue("discard_underdetermined").toString() == "false";
  multiple_hits_per_spec_ = param_.getValue("only_single_hit").toString() == "false";
  protein_level_only_ = param_.getValue("protein_level_only").toString() == "true";
  max_protein_hits_ = param_.getValue("max_protein_hits");
  protein_tol_ = param_.getValue("protein_tol");
  min_protein_hits_ = param_.getValue("min_protein_hits");
  protein_fdr_ = param_.getValue("protein_fdr");
  native_ = param_.getValue("native").toString() == "true";
  native_margin_ = param_.getValue("native_margin");
  if (native_)
  {
    gapped_mode_ = "off";
    multiple_hits_per_spec_ = true;
    // and at least one residual round (chimeric_rounds): on 25 native and targeted test spectra it added 4 expected proteins (three
    // antibody heavy chains, a co-isolated OMPF) without another protein called against the null databases (29 -> 33 of 43; other
    // proteins 8 -> 8; null-database calls 1 -> 1), at 4.2 instead of 3.3 s per spectrum
    chimeric_rounds_ = std::max(chimeric_rounds_, 1);
  }
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
      clusters.emplace_back(prsm_indices.begin(), prsm_indices.end());
      // The row that describes the proteoform: the PrSM with the fewest blind shifts among its determined PrSMs passing rep_q (then in the
      // order above). A blind shift equal to the residues a moved terminus left out explains the same mass worse (Eclip_1480 SODM: 2-206
      // unmodified and 6-206 +464.215, the mass of residues 2-5); which one represented the proteoform followed q-values differing in the
      // fourth digit. The cluster and its q-value stay with the representative chosen above (clusters keeps it first).
      Size row = 0;
      if (passing_determined[hit_idx[i]])
      {
        for (Size k = 1; k < prsm_indices.size(); k++)
        {
          const int m = prsm_indices[k];
          if (passing_determined[m] && blindModCount_(proteoform_hits_[m]) < blindModCount_(proteoform_hits_[prsm_indices[row]])) row = k;
        }
      }
      std::swap(prsm_indices[0], prsm_indices[row]);
      for (Size k = 0; k < prsm_indices.size(); k++) proteoform_hits_[prsm_indices[k]].setMetaValue("Representative", k == 0 ? 1 : 0);
      proteoform_hits_[prsm_indices[0]].setMetaValue("PrSMIndices", prsm_indices);
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
  //
  // A determined PrSM that is not the whole protein (an internal terminus: start > 2 or end before the protein's last residue), and every
  // precursor-ladder PrSM, is scored here with the model of the class "undetermined, both termini" (its null, pi0 and target density; that
  // class itself is unchanged), not with that of its determined class. Its mass agreement is closed by where the free terminus is put
  // (many windows of a long protein fit any mass), so it is no more evidence than termini without a mass; the determined classes have
  // 0 to 2 decoys, so their pi0 is about 1 / targets and their null the pooled one, and such a PrSM got a local FDR of about 0.002 at any
  // score. Over five E. coli lysate runs, 63 of 65 chance determined PrSMs (decoys and B. theta entrapment) were truncated; of the
  // determined PrSMs scoring 85 or less, whole-protein ones were 533 E. coli, 0 B. theta and 0 decoys. On Eclip_1480 against E. coli +
  // B. theta, the false proteins at 1% went from 5 to 4 (all four homologs of E. coli proteins, which shuffled decoys cannot model) and
  // the E. coli-only search from 299 to 292 proteins; on the held-out search against E. coli + human, from 2 false proteins to 0 (287 ->
  // 289 E. coli proteins). The PrSM and proteoform q-values are unchanged.
  const Size hit_count = proteoform_hits_.size();
  const int protein_class_count = fdr_class_count + 1, c3_model = fdr_class_count;
  std::vector<int> protein_class(hit_count, -1);
  for (Size i = 0; i < hit_count; i++)
  {
    const int c = hit_class[i];
    if (c < 0) continue;
    const auto& hit = proteoform_hits_[i];
    const int start = hit.getMetaValue("StartPosition"), end = hit.getMetaValue("EndPosition");
    const bool whole_protein = isDetermined(hit) && start >= 1 && start <= 2 && end == (int)hit.getSequence().size();
    protein_class[i] = (c == fdr_class_precursor_ladder || (c <= fdr_class_max_blind && ! whole_protein)) ? c3_model : c;
  }

  // null of the PrSM score per class
  std::vector<std::vector<double>> class_decoy_scores(protein_class_count);
  std::vector<double> class_decoy_weight(protein_class_count, 0.0), pooled_decoy_scores;
  for (Size i = 0; i < hit_count; i++)
  {
    if (protein_class[i] < 0 || ! decoys[i]) continue;
    class_decoy_scores[protein_class[i]].push_back(proteoform_hits_[i].getScore());
    class_decoy_weight[protein_class[i]] += decoy_weights[i]; // a decoy precursor of weight 0 still enters the null's sample
    pooled_decoy_scores.push_back(proteoform_hits_[i].getScore());
  }
  std::vector<DecoyScoreNull> nulls; // nulls[0]: the pooled one
  std::vector<Size> class_null(protein_class_count, 0);
  nulls.reserve(protein_class_count + 1);
  nulls.emplace_back(std::move(pooled_decoy_scores));
  for (int c = 0; c < protein_class_count; c++)
  {
    if (class_decoy_scores[c].size() < null_min_class_decoys) continue;
    class_null[c] = nulls.size();
    nulls.emplace_back(std::move(class_decoy_scores[c]));
  }
  std::vector<double> pvalues(hit_count, 1.0);
  std::vector<std::vector<double>> class_target_pvalues(protein_class_count);
  for (Size i = 0; i < hit_count; i++)
  {
    if (protein_class[i] < 0) continue;
    pvalues[i] = nulls[class_null[protein_class[i]]].pValue(proteoform_hits_[i].getScore());
    if (! decoys[i]) class_target_pvalues[protein_class[i]].push_back(pvalues[i]);
  }

  // two-group local FDR of every PrSM, targets and decoys, in its class
  std::vector<double> pi0;
  std::vector<TargetPvalueDensity> densities;
  densities.reserve(protein_class_count);
  for (int c = 0; c < protein_class_count; c++)
  {
    pi0.push_back(std::min(1.0, (class_decoy_weight[c] + 1.0) / std::max(1.0, (double)class_target_pvalues[c].size())));
    densities.emplace_back(std::move(class_target_pvalues[c]));
  }
  // the PrSMs scored with the model of "undetermined, both termini"
  for (Size i = 0; i < hit_count; i++)
  {
    if (protein_class[i] == c3_model) pvalues[i] = nulls[class_null[fdr_class_both_termini]].pValue(proteoform_hits_[i].getScore());
  }
  std::vector<double> local_fdr(hit_count, 1.0);
  for (Size i = 0; i < hit_count; i++)
  {
    if (protein_class[i] < 0) continue;
    const int model = protein_class[i] == c3_model ? fdr_class_both_termini : protein_class[i];
    const double density = densities[model].density(pvalues[i]); // +inf gives 0
    local_fdr[i] = density <= 0 ? 1.0 : std::min(1.0, pi0[model] / density);
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
  // proteins in parallel; each set is a copy of one built by the same insertions, as push_back copied it, so that it iterates in the same
  // order (FLASHExtenderAlgorithm::getBestMassShift_ breaks ties by that order)
  const Size first = vec_pro.size(), rev_first = rev_vec_pro.size();
  vec_pro.resize(first + cleaned_protein_seqs.size());
  rev_vec_pro.resize(rev_first + cleaned_protein_seqs.size());
#pragma omp parallel for schedule(dynamic, 64)
  for (SignedSize k = 0; k < (SignedSize)cleaned_protein_seqs.size(); k++)
  { //ignore X
    const auto& seq = cleaned_protein_seqs[k];
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
    vec_pro[first + k] = std::unordered_set<int>(vec);
    rev_vec_pro[rev_first + k] = std::unordered_set<int>(rev_vec);
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
  // per spectrum in dspecs: the tolerance of its deconvolution and, for the chimeric search only, all its masses (ascending) if the cap
  // max_node_cntr_ left some out (else empty: dspecs has them all)
  std::vector<double> spectrum_tols;
  std::vector<std::vector<PeakGroup>> uncapped_masses;

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

  // Spectra in parallel, each into its own slot, collected in input order below, so that dspecs, the tolerances and tags_ are those of one
  // pass in input order. An exception (e.g. a malformed DeconvMassInfo) is rethrown for the first spectrum in input order that raised one.
  struct TaggedSpectrum
  {
    bool with_masses = false; // at least 5 masses: ms_level and tol are read
    UInt ms_level = 0;
    double tol = 0;
    bool searched = false; // dspec, all_masses and tags are set
    DeconvolvedSpectrum dspec;
    std::vector<PeakGroup> all_masses;
    std::vector<FLASHHelperClasses::Tag> tags;
    std::exception_ptr error;
  };
  std::vector<TaggedSpectrum> tagged(map.size());
  auto tagSpectrum = [&](int index, TaggedSpectrum& out) {
    const auto& spec = map[index];
    nextProgress();

    if (spec.size() < 5) return;

    int scan = FLASHDeconvAlgorithm::getScanNumber(map, index);

    DeconvolvedSpectrum dspec(scan);
    dspec.setOriginalSpectrum(spec);
    String deconv_meta_str = spec.getMetaValue("DeconvMassInfo").toString();

    int tol_loc_s = deconv_meta_str.find("tol=") + 4;
    int tol_loc_e = deconv_meta_str.find(";", tol_loc_s);

    const double spec_tol = stod(deconv_meta_str.substr(tol_loc_s, tol_loc_e - tol_loc_s));
    out.with_masses = true;
    out.ms_level = spec.getMSLevel();
    out.tol = spec_tol;
    if (spec.getMSLevel() == 1 && has_msn) return;

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
      if (deconv_meta_str.hasSubstring("decoyprecursorQscore1D="))
      {
        int s_loc_dq1_s = deconv_meta_str.find("decoyprecursorQscore1D=") + 23;
        int s_loc_dq1_e = deconv_meta_str.find(";", s_loc_dq1_s);
        dspec.setDecoyPrecursorQscore1D(stod(deconv_meta_str.substr(s_loc_dq1_s, s_loc_dq1_e == (int)String::npos ? String::npos : s_loc_dq1_e - s_loc_dq1_s)));
      }
      if (deconv_meta_str.hasSubstring("decoyprecursorQscore="))
      {
        int s_loc_dq_s = deconv_meta_str.find("decoyprecursorQscore=") + 21;
        int s_loc_dq_e = deconv_meta_str.find(";", s_loc_dq_s);
        dspec.setDecoyPrecursorQscore(stod(deconv_meta_str.substr(s_loc_dq_s, s_loc_dq_e == (int)String::npos ? String::npos : s_loc_dq_e - s_loc_dq_s)));
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
      if (deconv_meta_str.hasSubstring("precursorQscore1D="))
      {
        int s_loc_pq1_s = deconv_meta_str.find("precursorQscore1D=") + 18; // the precursor's comes before the decoy precursor's
        int s_loc_pq1_e = deconv_meta_str.find(";", s_loc_pq1_s);
        pg.setQscore(stod(deconv_meta_str.substr(s_loc_pq1_s, s_loc_pq1_e - s_loc_pq1_s)));
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

    std::vector<PeakGroup> all_masses; // the chimeric search applies the cap again to what is left of all masses
    if (chimeric_rounds_ > 0 && dspec.size() > max_node_cntr_) all_masses.assign(dspec.begin(), dspec.end());
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
    if (chimeric_rounds_ > 0)
    {
      std::sort(all_masses.begin(), all_masses.end()); // as dspec.sort()
      out.all_masses = std::move(all_masses);
    }
    FLASHTaggerAlgorithm tagger;
    // Run tagger
    tagger.setParameters(tagger_param_);
    tagger.run(dspec, spec_tol);
    tagger.fillTags(out.tags);
    out.dspec = std::move(dspec);
    out.searched = true;
  };
#pragma omp parallel for schedule(dynamic)
  for (int index = 0; index < (int)map.size(); index++)
  {
    try
    {
      tagSpectrum(index, tagged[index]);
    }
    catch (...)
    {
      tagged[index].error = std::current_exception();
    }
  }
  for (auto& out : tagged)
  {
    if (out.error) std::rethrow_exception(out.error);
    if (! out.with_masses) continue;
    tol = out.tol;
    if (out.ms_level > 1 && ms2_tol < 0) ms2_tol = tol;
    if (out.ms_level == 1 && precursor_tol < 0) precursor_tol = tol;
    if (! out.searched) continue;
    dspecs.push_back(std::move(out.dspec));
    spectrum_tols.push_back(tol);
    if (chimeric_rounds_ > 0) uncapped_masses.push_back(std::move(out.all_masses));
    for (auto& tag : out.tags) tags_.push_back(std::move(tag));
    out = TaggedSpectrum();
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

  // Internal residue masses by one-letter code, for the flanking masses of the tag hits. AASequence::fromString() of a string of these
  // letters is one residue per letter (ResidueDB::getResidue(char)), and AASequence::getMonoWeight(Residue::Internal) adds their
  // masses to 0, from the first residue to the last; the same additions in the same order give the same doubles. A protein with any
  // other character (e.g. X, which cleaned_protein_seqs drops, so that its positions are no longer those of the hits) keeps the
  // AASequence computation.
  std::array<double, 256> aa_internal_mass {};
  std::array<bool, 256> aa_plain {};
  {
    const ResidueDB* rdb = ResidueDB::getInstance();
    const Residue* unknown = rdb->getResidue(String("X"));
    for (int c = 'A'; c <= 'Z'; c++)
    {
      if (c == 'X') continue;
      const Residue* r = rdb->getResidue((unsigned char)c);
      if (r == nullptr || r == unknown) continue;
      aa_plain[c] = true;
      aa_internal_mass[c] = r->getMonoWeight(Residue::Internal);
    }
  }

  // The database positions and flanking masses of the tag sequences @p seqs (their entries in tag_to_protein_indices and the other maps
  // exist, empty): an Aho-Corasick trie of them over the proteins. Also used for the new tag sequences of the chimeric search.
  // Proteins run in parallel in blocks (each thread with its own trie state); the hits of a block are appended to the maps in protein
  // order, so every vector is the same as that of one pass over the proteins.
  auto searchTagSequences = [&](const std::vector<std::string>& seqs) {
    ACTrie ac_trie(0, 0); // no ambiguities, no mismatches
    ac_trie.addNeedlesAndCompress(seqs);

    // the vectors of every tag sequence (std::map elements do not move)
    std::vector<std::vector<Size>*> indices_of(seqs.size()), positions_of(seqs.size());
    std::vector<std::vector<double>*> n_masses_of(seqs.size()), c_masses_of(seqs.size());
    for (Size j = 0; j < seqs.size(); j++)
    {
      const String seq = seqs[j];
      indices_of[j] = &tag_to_protein_indices[seq];
      positions_of[j] = &tag_to_protein_positions[seq];
      n_masses_of[j] = &tag_to_n_flanking_masses[seq];
      c_masses_of[j] = &tag_to_c_flanking_masses[seq];
    }

    struct TagHit
    {
      OpenMS::Hit::T needle_index;
      OpenMS::Hit::T query_pos;
      bool summed; // false: the flanking masses are computed with AASequence when the hit is appended (serially, in protein order)
      double nmass;
      double cmass;
    };

    startProgress(0, (SignedSize)protein_seqs.size(), "Searching tags against database ...");

    const Size block_size = 4096;
    std::vector<std::vector<TagHit>> block_hits(block_size);
    for (Size block_start = 0; block_start < protein_seqs.size(); block_start += block_size)
    {
      const Size block_end = std::min(protein_seqs.size(), block_start + block_size);
#pragma omp parallel
      {
        OpenMS::ACTrieState ac_state;
        std::vector<double> prefix; // residue sums from the N-terminus
#pragma omp for schedule(dynamic, 16)
        for (SignedSize i = (SignedSize)block_start; i < (SignedSize)block_end; i++) // // Natural19WithoutI
        {
          auto& hits = block_hits[i - block_start];
          hits.clear();
          const std::string& pseq = protein_seqs[i];
          ac_state.setQuery(pseq);
          ac_trie.getAllHits(ac_state);
          if (ac_state.hits.empty()) continue;
          hits.reserve(ac_state.hits.size());

          bool plain = true;
          for (const char c : pseq)
          {
            if (! aa_plain[(unsigned char)c])
            {
              plain = false;
              break;
            }
          }
          const Size n = pseq.size();
          if (plain)
          {
            prefix.resize(n + 1);
            prefix[0] = .0;
            for (Size k = 0; k < n; k++) prefix[k + 1] = prefix[k] + aa_internal_mass[(unsigned char)pseq[k]];
          }
          for (const auto& h : ac_state.hits)
          {
            const Size end = (Size)h.query_pos + seqs[h.needle_index].length();
            if (! plain || end > n)
            {
              hits.push_back({h.needle_index, h.query_pos, false, .0, .0});
              continue;
            }
            const double nmass = h.query_pos <= 0 ? .0 : prefix[h.query_pos];
            double cmass = .0;
            for (Size k = end; k < n; k++) cmass += aa_internal_mass[(unsigned char)pseq[k]];
            hits.push_back({h.needle_index, h.query_pos, true, nmass, cmass});
          }
        }
      }

      for (Size i = block_start; i < block_end; i++)
      {
        nextProgress();
        for (const auto& h : block_hits[i - block_start])
        {
          double nmass = h.nmass, cmass = h.cmass;
          if (! h.summed)
          {
            const String& cleanseq = cleaned_protein_seqs[i];
            const String& seq = seqs[h.needle_index];
            nmass = h.query_pos <= 0 ? .0 : AASequence::fromString(cleanseq.substr(0, h.query_pos)).getMonoWeight(Residue::Internal);
            cmass = cleanseq.length() - seq.length() <= h.query_pos ? .0 :
              AASequence::fromString(cleanseq.substr(h.query_pos + seq.length(), cleanseq.length() - h.query_pos - seq.length())).getMonoWeight(Residue::Internal);
          }
          indices_of[h.needle_index]->push_back(i);
          positions_of[h.needle_index]->push_back(h.query_pos);
          n_masses_of[h.needle_index]->push_back(nmass);
          c_masses_of[h.needle_index]->push_back(cmass);
        }
      }
    }

    endProgress();
  };
  searchTagSequences(tag_seqs);
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
             mi_sliding = registry.registerName("SlidingLadder"),
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

  // Mass ladders (MassLadderAlgorithm): terminal ladders from the residue sums of every protein from its termini, and precursor-anchored
  // ladders from the runs of standard residues of the database, for the window search of every spectrum with a precursor mass. The
  // precursor ladders are built only if a spectrum has a precursor mass.
  // Native mode (native): the candidates of a source in a spectrum are kept when their score exceeds that of the best decoy candidate of
  // the same source in the spectrum by native_margin_ (decoys thereby never), best first, at most max_protein_hits_. A spectrum is its own
  // null here: native runs have one or a few spectra, too few for an FDR over spectra. Native mode reads ladders from any start or end of a
  // protein (sliding ladders, MassLadderAlgorithm::searchSlidingLadders) instead of the intact termini only: the outer membrane proteins
  // of the targeted test spectra lack their signal peptides (OMPA from residue 22). Gapped-tag candidates were never an expected protein
  // of the native test spectra, so they are not used. Margin 2 (a hundred times the chance probability of the best decoy) was set
  // beforehand. On 25 native and targeted test spectra (OMPA, OMPF, TSX, DLDH and GapA mixtures, DLDH, AQPZ, NISTmAb 1-6, SIGMAmAb,
  // Infliximab, carbonic anhydrase, VAMP2, semiSWEET; 43 expected proteins) against E. coli + the expected proteins, native mode found as
  // many expected proteins as without it (29) and 8 other proteins instead of 27 (6 instead of 8 of them chains of another antibody or
  // porins next to OMPF); against databases without the expected proteins, it called 1 protein instead of 33.
  const bool native_rule = native_ && decoy_factor_ > 0;
  if (native_ && ! native_rule)
  {
    OPENMS_LOG_WARN << "native: the database has no decoys (DECOY prefix); the candidates are chosen as without the option." << std::endl;
  }
  bool any_precursor_mass = false;
  for (const auto& dspec : dspecs) any_precursor_mass |= dspec.getPrecursorPeakGroup().getMonoMass() > 0;
  const bool use_precursor_ladders = precursor_ladders_ && ! protein_level_only_ && any_precursor_mass;
  MassLadderAlgorithm mass_ladders;
  // the sliding ladders of native mode read the runs of standard residues of the precursor ladders
  mass_ladders.buildIndex(protein_seqs, param_.getValue("ion_type").toStringVector(), terminal_ladders_, use_precursor_ladders || native_rule);
  // the precursor tolerance of the windows: the MS1 tolerance of the input, else (no MS1 spectra) the fragment tolerance
  const double pl_precursor_tol = precursor_tol > 0 ? precursor_tol : protein_tol;
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
  const Size native_ladder_top = 1000; // ladder candidates ranked per spectrum and source for the decoy margin
  const std::vector<String> native_sources {"sequence tags", "sliding ladders", "precursor ladders"};
  std::vector<std::atomic<Size>> native_kept(native_sources.size());
  for (auto& k : native_kept) k = 0;
  // @p truncated: the ranking stopped at native_ladder_top, so a decoy missing from it scored at most the last candidate's score
  auto nativeCut = [&](auto& cands, auto score_of, auto protein_of, bool truncated, int source, String& note) {
    auto is_decoy = [&](const auto& c) { return fasta_entry[protein_of(c)].identifier.hasPrefix("DECOY"); };
    std::stable_sort(cands.begin(), cands.end(), [&](const auto& a, const auto& b) { return score_of(a) > score_of(b); });
    double best_decoy = -1;
    for (const auto& c : cands)
    {
      if (is_decoy(c)) { best_decoy = std::max(best_decoy, (double)score_of(c)); }
    }
    if (best_decoy < 0) best_decoy = truncated && ! cands.empty() ? (double)score_of(cands.back()) : 0.0;
    const Size n = cands.size();
    std::remove_reference_t<decltype(cands)> kept;
    String names;
    for (auto& c : cands)
    {
      if ((int)kept.size() >= max_protein_hits_ || score_of(c) < best_decoy + native_margin_) break;
      if (is_decoy(c)) continue;
      names += String(kept.empty() ? " " : ", ") + fasta_entry[protein_of(c)].identifier + " " + String::number(score_of(c), 1);
      kept.push_back(std::move(c));
    }
    native_kept[source] += kept.size();
    note += "; " + native_sources[source] + ": " + String(n) + " candidates, best decoy " + String::number(best_decoy, 1) + ", kept " +
            String(kept.size()) + names;
    cands.swap(kept);
  };
  // A mass-ladder candidate as a protein hit; its ladders become tags (@p gtags) spanning them, with the lowest mass as flanking mass. The
  // tag indices of the hit are -(local index + 1), made global with the gapped tags when appended to tags_.
  auto ladderCandidateHit = [&](const DeconvolvedSpectrum& dspec, const MassLadderAlgorithm::Candidate& cand,
                                std::vector<FLASHHelperClasses::Tag>& gtags) {
    const auto& pseq = protein_seqs[cand.protein];
    std::vector<int> indices, positions;
    int matched = 0;
    for (const auto& ladder : cand.ladders)
    {
      std::vector<double> tag_masses;
      std::vector<int> tag_scores;
      for (Size mi : ladder.mass_indices)
      {
        tag_masses.push_back(dspec[mi].getMonoMass());
        tag_scores.push_back(FLASHTaggerAlgorithm::getNodeScore(dspec[mi]));
      }
      matched += (int)tag_masses.size();
      gtags.emplace_back(String(pseq.substr(ladder.start, ladder.end - ladder.start)), ladder.prefix ? tag_masses[0] : -1.0,
                         ladder.prefix ? -1.0 : tag_masses[0], tag_masses, tag_scores, dspec.getScanNumber(), (Size)(ladder.end - ladder.start));
      gtags.back().setRetentionTime((float)dspec.getOriginalSpectrum().getRT());
      indices.push_back(-(int)gtags.size());
      positions.push_back(ladder.start);
    }
    const auto& fe = fasta_entry[cand.protein];
    ProteinHit hit(cand.score, 0, fe.identifier, fe.sequence);
    hit.setDescription(fe.description);
    hit.setMetaValue(mi_scan, dspec.getScanNumber());
    hit.setMetaValue(mi_matched_aa, matched);
    hit.setCoverage((double)matched / (double)fe.sequence.length());
    hit.setMetaValue(mi_tag_indices, indices);
    std::sort(positions.begin(), positions.end());
    hit.setMetaValue(mi_tag_positions, positions);
    hit.setMetaValue(mi_fasta_index, cand.protein);
    hit.setMetaValue(mi_protein_rank, cand.rank);
    return hit;
  };
  // The candidate proteins of one spectrum: from its sequence tags (@p tag_indices_of_scan, indices into tags_; none: nullptr), gapped tags,
  // terminal ladders and, with @p with_precursor_ladders, precursor ladders. The proteins in @p excluded (fasta indices) are not proposed
  // by any source (the chimeric search: the proteins of the spectrum's earlier best PrSMs). Writes only its outputs: the sequence- and
  // gapped-tag candidates (@p hits_out, ranked), the terminal- and precursor-ladder candidates, the gapped and ladder tags (@p gtags, local
  // indices until appended to tags_), the number of precursor-ladder windows and, with native, what the decoy margin kept (@p note).
  // @p stage: 0 every source as configured; a stage of the search cascade only its sources: 1 sequence tags, 2 gapped tags (in every
  // spectrum searched, as with gapped_tags only), 3 terminal and precursor ladders.
  auto searchCandidates = [&](const DeconvolvedSpectrum& dspec, const std::vector<Size>* tag_indices_of_scan, const std::set<Size>& excluded,
                              bool with_precursor_ladders, int stage, std::vector<ProteinHit>& hits_out, std::vector<ProteinHit>& ladder_out,
                              std::vector<ProteinHit>& pl_out, std::vector<ProteinHit>& decoy_pl_out,
                              std::vector<FLASHHelperClasses::Tag>& gtags, Size& pl_windows, String& note) {
    const bool sequence_stage = stage == 0 || stage == 1, ladder_stage = stage == 0 || stage == 3;
    with_precursor_ladders &= ladder_stage;
    // The sequence-tag candidates (proteins hit by a tag of the spectrum; thousands per spectrum) stay records until the cut below, and only
    // the kept ones become ProteinHits. Per candidate, ascending by protein: its covered positions and its tags (distinct, ascending), and
    // its shifts (protein residue sum up to the tag minus the tag flanking mass) with the tag flanking masses, in the order of the tags and
    // of their hits (as maps per protein collected them).
    struct SeqCandidate
    {
      Size protein;
      Size covered_pos;                  // distinct covered positions, in cand_positions from positions_first
      Size positions_first, tags_first;  // tags: in cand_tags
      Size tag_count;
      Size n_first, n_count, c_first, c_count; // shifts: in shift_values and flank_values
    };
    std::vector<SeqCandidate> seq_candidates;
    std::vector<Size> cand_positions, cand_tags;
    std::vector<double> shift_values, flank_values;
    std::vector<ProteinHit> hits; // the gapped-tag candidates, which come after the sequence-tag ones
    std::vector<FLASHTaggerAlgorithm::ProteinShifts> hit_shifts; // in the order of hits
    if (sequence_stage && gapped_mode_ != "only" && tag_indices_of_scan != nullptr)
    {
      const auto& tag_indices = *tag_indices_of_scan;
      // the hits of the tags in tag order, each with the slot of its protein (protein order of first hit)
      struct ShiftHit
      {
        Size slot;
        double shift, flank;
      };
      struct PositionHit
      {
        Size slot, tag, position;
      };
      std::vector<ShiftHit> n_hits, c_hits;
      std::vector<PositionHit> position_hits;
      std::vector<int> slot_of(cleaned_protein_seqs.size(), -1);
      std::vector<Size> slot_protein;
      for (const auto& tag_index : tag_indices)
      {
        const auto& tag = tags_[tag_index];
        const auto& seq = tag.getUppercaseSequence();

        const auto& pis = tag_to_protein_indices.at(seq);
        const auto& pps = tag_to_protein_positions.at(seq);
        const auto& n_flanking_masses = tag_to_n_flanking_masses.at(seq);
        const auto& c_flanking_masses = tag_to_c_flanking_masses.at(seq);

        for (Size i = 0; i < pis.size(); i++)
        {
          if (excluded.count(pis[i]) > 0) continue;
          const auto& pseq = cleaned_protein_seqs[pis[i]];
          int& slot = slot_of[pis[i]];
          if (slot < 0)
          {
            slot = (int)slot_protein.size();
            slot_protein.push_back(pis[i]);
          }

          if (tag.getNtermMass() >= 0 && pps[i] > 0 && pps[i] <= pseq.length() + seq.length())
          {
            double nmass = n_flanking_masses[i];
            if (tag.getNtermMass() > nmass + max_mod_mass + 200) continue;
            n_hits.push_back({(Size)slot, nmass - tag.getNtermMass(), tag.getNtermMass()});
          }
          if (tag.getCtermMass() >= 0 && pps[i] + seq.length() >= 0 && pseq.length() > pps[i] + seq.length())
          {
            double cmass = c_flanking_masses[i];
            if (tag.getCtermMass() > cmass + max_mod_mass + 200) continue;
            c_hits.push_back({(Size)slot, cmass - tag.getCtermMass(), tag.getCtermMass()});
          }
          position_hits.push_back({(Size)slot, tag_index, pps[i]});
        }
      }

      // the hits of each slot, in their order: positions [first[s], first[s + 1]) of the returned order
      const Size slot_count = slot_protein.size();
      auto bySlot = [slot_count](const auto& slot_hits, std::vector<Size>& first) {
        first.assign(slot_count + 1, 0);
        for (const auto& h : slot_hits) first[h.slot + 1]++;
        for (Size s = 0; s < slot_count; s++) first[s + 1] += first[s];
        std::vector<Size> next(first.begin(), first.end() - 1), order(slot_hits.size());
        for (Size k = 0; k < slot_hits.size(); k++) order[next[slot_hits[k].slot]++] = k;
        return order;
      };
      std::vector<Size> n_first, c_first, p_first;
      const auto n_order = bySlot(n_hits, n_first), c_order = bySlot(c_hits, c_first), p_order = bySlot(position_hits, p_first);

      std::vector<Size> covered_slots; // with a covered position, ascending by protein
      for (Size s = 0; s < slot_count; s++)
      {
        if (p_first[s + 1] > p_first[s]) covered_slots.push_back(s);
      }
      std::sort(covered_slots.begin(), covered_slots.end(), [&](Size a, Size b) { return slot_protein[a] < slot_protein[b]; });
      seq_candidates.reserve(covered_slots.size());
      for (Size s : covered_slots)
      {
        SeqCandidate cand;
        cand.protein = slot_protein[s];
        cand.positions_first = cand_positions.size();
        cand.tags_first = cand_tags.size();
        for (Size k = p_first[s]; k < p_first[s + 1]; k++)
        {
          cand_positions.push_back(position_hits[p_order[k]].position);
          cand_tags.push_back(position_hits[p_order[k]].tag);
        }
        std::sort(cand_positions.begin() + cand.positions_first, cand_positions.end());
        cand_positions.erase(std::unique(cand_positions.begin() + cand.positions_first, cand_positions.end()), cand_positions.end());
        std::sort(cand_tags.begin() + cand.tags_first, cand_tags.end());
        cand_tags.erase(std::unique(cand_tags.begin() + cand.tags_first, cand_tags.end()), cand_tags.end());
        cand.covered_pos = cand_positions.size() - cand.positions_first;
        cand.tag_count = cand_tags.size() - cand.tags_first;
        if (cand.covered_pos < min_tag_count) continue;
        cand.n_first = shift_values.size();
        for (Size k = n_first[s]; k < n_first[s + 1]; k++)
        {
          shift_values.push_back(n_hits[n_order[k]].shift);
          flank_values.push_back(n_hits[n_order[k]].flank);
        }
        cand.n_count = shift_values.size() - cand.n_first;
        cand.c_first = shift_values.size();
        for (Size k = c_first[s]; k < c_first[s + 1]; k++)
        {
          shift_values.push_back(c_hits[c_order[k]].shift);
          flank_values.push_back(c_hits[c_order[k]].flank);
        }
        cand.c_count = shift_values.size() - cand.c_first;
        seq_candidates.push_back(cand);
      }
    }

    // gapped tags: proteins proposed where gapped tags land. Candidates are scored and cut first; only the survivors get their
    // sequence and their landings as tags (a landing becomes a tag spanning the landed residues, with the masses at the gap
    // boundaries), so the protein score (shifts) and the extension (anchors) use it as they use a sequence tag.
    std::vector<FLASHGappedTaggerAlgorithm::Landing> landings;
    std::vector<std::vector<Size>> gapped_groups; // per gapped candidate, its landings
    if (stage == 2 || (stage == 0 && (gapped_mode_ == "only" || gapped_mode_ == "both" || (gapped_mode_ == "fallback" && seq_candidates.empty()))))
    {
      std::set<Size> proposed; // proteins already proposed by sequence tags
      for (const auto& c : seq_candidates) proposed.insert(c.protein);
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
        for (Size l = l0; l < l1 && proposed.count(pi) == 0 && excluded.count(pi) == 0; ++l)
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
    // node scores let chance matches (more of them on longer proteins) rank noise proteins above the correct one. The candidates in the
    // order they had as hits: the sequence-tag ones, then the gapped-tag ones.
    const Size seq_count = seq_candidates.size(), candidate_count = seq_count + hits.size();
    std::vector<double> score_masses {0.0}; // mass 0 stands for the terminus (FLASHTaggerAlgorithm::runMatchingExact)
    score_masses.reserve(dspec.size() + 1);
    for (const auto& pg : dspec) score_masses.push_back(pg.getMonoMass());
    std::vector<double> scores(candidate_count, 0.0);
    // a loop that is already parallel over spectra runs this serially
#pragma omp parallel for schedule(dynamic) if (! omp_in_parallel())
    for (int i = 0; i < (int)candidate_count; i++)
    {
      if ((Size)i < seq_count)
      {
        const auto& c = seq_candidates[i];
        scores[i] = std::max(FLASHTaggerAlgorithm::scoreShifts(score_masses, pro_prefix[c.protein], shift_values.data() + c.n_first,
                                                               flank_values.data() + c.n_first, c.n_count, protein_tol),
                             FLASHTaggerAlgorithm::scoreShifts(score_masses, pro_suffix[c.protein], shift_values.data() + c.c_first,
                                                               flank_values.data() + c.c_first, c.c_count, protein_tol));
        continue;
      }
      const auto& sh = hit_shifts[i - seq_count];
      scores[i] = std::max(FLASHTaggerAlgorithm::scoreShifts(score_masses, pro_prefix[sh.fasta_index], sh.n_shifts.data(), sh.n_flanks.data(),
                                                             sh.n_shifts.size(), protein_tol),
                           FLASHTaggerAlgorithm::scoreShifts(score_masses, pro_suffix[sh.fasta_index], sh.c_shifts.data(), sh.c_flanks.data(),
                                                             sh.c_shifts.size(), protein_tol));
    }

    // Ranked by score, coverage and description with std::sort. The keys are those the hits had, so every comparison and with it the
    // order, ties included, is that of sorting the hits (std::sort's moves depend on the comparisons only).
    struct RankedCandidate
    {
      double score, coverage;
      const String* description;
      Size index;
    };
    std::vector<RankedCandidate> ranked(candidate_count);
    for (Size i = 0; i < candidate_count; i++)
    {
      if (i < seq_count)
      {
        const auto& fe = fasta_entry[seq_candidates[i].protein];
        ranked[i] = {scores[i], (double)seq_candidates[i].covered_pos / (double)fe.sequence.length(), &fe.description, i};
      }
      else { ranked[i] = {scores[i], hits[i - seq_count].getCoverage(), &fasta_entry[hit_shifts[i - seq_count].fasta_index].description, i}; }
    }
    std::sort(ranked.begin(), ranked.end(), [](const RankedCandidate& left, const RankedCandidate& right) {
      return left.score == right.score ? (left.coverage == right.coverage ? (*left.description > *right.description) : (left.coverage > right.coverage))
                                       : (left.score > right.score);
    });

    // the candidate as a hit with its score: a sequence-tag one is built as the maps per protein built it, a gapped-tag one is moved
    auto candidateHit = [&](Size i) {
      if (i >= seq_count)
      {
        ProteinHit hit = std::move(hits[i - seq_count]);
        hit.setScore(scores[i]);
        return hit;
      }
      const auto& c = seq_candidates[i];
      const auto& fe = fasta_entry[c.protein];
      Size covered_pos = c.covered_pos;
      double coverage = (double)covered_pos / (double)fe.sequence.length();
      std::vector<int> positions(cand_positions.begin() + c.positions_first, cand_positions.begin() + (c.positions_first + covered_pos));
      ProteinHit hit(coverage, 0, fe.identifier, fe.sequence);
      hit.setDescription(fe.description);
      hit.setMetaValue(mi_scan, dspec.getScanNumber());
      hit.setMetaValue(mi_matched_aa, covered_pos);
      hit.setCoverage(coverage);
      hit.setMetaValue(mi_tag_indices, std::vector<int>(cand_tags.begin() + c.tags_first, cand_tags.begin() + (c.tags_first + c.tag_count)));
      hit.setMetaValue(mi_tag_positions, positions);
      hit.setMetaValue(mi_fasta_index, c.protein);
      hit.setScore(scores[i]);
      return hit;
    };
    std::vector<ProteinHit> ranked_hits;
    if (native_rule)
    {
      ranked_hits.reserve(candidate_count);
      for (const auto& r : ranked) ranked_hits.push_back(candidateHit(r.index));
      nativeCut(ranked_hits, [](const ProteinHit& h) { return h.getScore(); }, [&](const ProteinHit& h) { return (Size)(int)h.getMetaValue(mi_fasta_index); },
                false, 0, note);
    }
    else
    {
      // The score above uses the masses of the whole spectrum against the protein; a correct protein with few of them can end up below
      // the cut. The number of distinct matched tags is sequence-specific evidence, so keep the best hits by that criterion as well and
      // let the extension decide between them.
      std::vector<bool> keep(candidate_count, (int)candidate_count <= max_protein_hits_);
      if ((int)candidate_count > max_protein_hits_)
      {
        for (int i = 0; i < max_protein_hits_; i++)
        {
          keep[i] = true; // ranked by score
        }

        std::vector<std::pair<int, Size>> by_tag_count;
        by_tag_count.reserve(candidate_count);
        for (Size i = 0; i < candidate_count; i++)
        {
          const Size c = ranked[i].index;
          const Size tag_count = c < seq_count ? seq_candidates[c].tag_count : hits[c - seq_count].getMetaValue(mi_tag_indices).toIntList().size();
          by_tag_count.emplace_back(-(int)tag_count, i);
        }
        std::partial_sort(by_tag_count.begin(), by_tag_count.begin() + max_protein_hits_, by_tag_count.end());
        for (int i = 0; i < max_protein_hits_; i++)
        {
          keep[by_tag_count[i].second] = true;
        }
      }
      for (Size i = 0; i < candidate_count; i++)
      {
        if (keep[i]) ranked_hits.push_back(candidateHit(ranked[i].index)); // in score order
      }
    }
    hits = std::move(ranked_hits);

    // the surviving gapped candidates: sequence and tags (indices local to this spectrum until appended to tags_ below)
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

    // the deconvolved masses (ascending) for the mass ladders
    std::vector<double> masses;
    if (ladder_stage && (terminal_ladders_ || with_precursor_ladders || native_rule) && ! dspec.empty())
    {
      masses.reserve(dspec.size());
      for (const auto& pg : dspec) masses.push_back(pg.getMonoMass());
    }

    // Terminal ladders (MassLadderAlgorithm::searchTerminalLadders): proteins by the chance probability of the masses landing on a residue
    // sum from a terminus. The min_protein_hits best proteins are characterized; the ladder of the best shift of each side becomes a tag
    // spanning it (as a gapped tag does), so the extension is anchored there. They are kept apart from the other candidates, also when a
    // sequence tag proposed the same protein.
    // Native mode: sliding ladders (from any start or to any end of a protein) instead, kept by the decoy margin.
    if (ladder_stage && (terminal_ladders_ || native_rule) && ! dspec.empty())
    {
      std::vector<MassLadderAlgorithm::Candidate> found;
      if (native_rule)
      {
        mass_ladders.searchSlidingLadders(masses, protein_tol, native_ladder_top, excluded, found);
        nativeCut(found, [](const MassLadderAlgorithm::Candidate& c) { return c.score; }, [](const MassLadderAlgorithm::Candidate& c) { return c.protein; },
                  found.size() >= native_ladder_top, 1, note);
      }
      else { mass_ladders.searchTerminalLadders(masses, protein_tol, (Size)min_protein_hits_, excluded, found); }
      for (const auto& cand : found)
      {
        ProteinHit hit = ladderCandidateHit(dspec, cand, gtags);
        hit.setMetaValue(native_rule ? mi_sliding : mi_terminal, 1);
        ladder_out.push_back(std::move(hit));
      }
    }

    // Precursor-anchored ladders: the windows of the database (targets and decoys alike) whose unmodified mass is the precursor mass, by
    // the chance probability of the masses on their prefix and suffix ladders (MassLadderAlgorithm::searchPrecursorLadders). The
    // min_protein_hits best windows are characterized; the matched masses of the best prefix and the best suffix ion type become tags
    // spanning them (as terminal-ladder tags), from which the extension can start and end the proteoform at the window. FLASHExtender places
    // the tags by their sequence and reaches an internal start or end only within its truncation range around the protein's best mass shift
    // (getBestMassShift_), not at the window it is given: the window is stored on the candidate (PrecursorLadderStart and End, 1-based as
    // StartPosition and EndPosition, and PrecursorLadderIsotope), and the log line after the extension counts the PrSMs at its ends. Kept
    // apart from the other candidates; searched before the check for tag candidates, so spectra without any are searched too. The same at
    // the decoy precursor mass (decoy_precursor).
    if (with_precursor_ladders && ! dspec.empty() && dspec.getPrecursorPeakGroup().getMonoMass() > 0)
    {
      std::vector<double> precursor_masses {dspec.getPrecursorPeakGroup().getMonoMass()};
      if (decoy_precursor_ && dspec.getDecoyPrecursorMass() > 0) precursor_masses.push_back(dspec.getDecoyPrecursorMass());
      std::vector<MassLadderAlgorithm::PrecursorLadderResult> found;
      mass_ladders.searchPrecursorLadders(masses, precursor_masses, protein_tol, pl_precursor_tol, native_rule ? native_ladder_top : (Size)min_protein_hits_, found);
      pl_windows = found[0].windows;
      for (auto& result : found)
      {
        if (native_rule)
        {
          nativeCut(result.candidates, [](const MassLadderAlgorithm::Candidate& c) { return c.score; },
                    [](const MassLadderAlgorithm::Candidate& c) { return c.protein; }, result.candidates.size() >= native_ladder_top, 2, note);
        }
      }
      for (Size n = 0; n < found.size(); n++)
      {
        for (const auto& cand : found[n].candidates)
        {
          ProteinHit hit = ladderCandidateHit(dspec, cand, gtags);
          hit.setMetaValue(mi_precursor_ladder, 1);
          hit.setMetaValue(mi_pl_start, cand.window_start + 1); // residues i..j-1 (0-based) are i+1..j (1-based, inclusive)
          hit.setMetaValue(mi_pl_end, cand.window_end);
          hit.setMetaValue(mi_pl_isotope, cand.isotope);
          (n == 0 ? pl_out : decoy_pl_out).push_back(std::move(hit));
        }
      }
    }

    if (hits.empty()) return;
    for (Size i = 0; i < hits.size(); i++)
    {
      hits[i].setMetaValue(mi_protein_rank, (int)i + 1);
    }
    hits_out = std::move(hits);
  };
  const std::set<Size> no_exclusion;
  // gapped and ladder tags of a spectrum join tags_; returns the index of the first, from which the tag indices of its candidates become global
  auto appendSpectrumTags = [&](std::vector<FLASHHelperClasses::Tag>& gtags) {
    const int tag_base = (int)tags_.size();
    for (auto& t : gtags)
    {
      t.setIndex((int)tags_.size());
      tags_.push_back(std::move(t));
    }
    return tag_base;
  };
  // sequence- and gapped-tag candidates: the local tag indices of gapped candidates become global
  auto appendTagCandidates = [&](std::vector<ProteinHit>& from, std::vector<ProteinHit>& to, std::pair<Size, Size>& range, int tag_base) {
    range.first = to.size();
    for (auto& hit : from)
    {
      if (hit.metaValueExists(mi_gapped))
      {
        auto indices = hit.getMetaValue(mi_tag_indices).toIntList();
        for (auto& i : indices) i += tag_base;
        hit.setMetaValue(mi_tag_indices, indices);
      }
      to.push_back(std::move(hit));
    }
    range.second = to.size();
  };
  // ladder candidates (terminal and precursor ladders): their tag indices -(local index + 1) become global
  auto appendLadderCandidates = [&](std::vector<ProteinHit>& from, std::vector<ProteinHit>& to, std::pair<Size, Size>& range, int tag_base) {
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

  // Search cascade (search_cascade): the sources in stages, each searching only the spectra that no stage before identified (a target
  // PrSM passing the PrSM-level FDR of its stage, see identifiedInStage below): 1 sequence tags, 2 gapped tags, 3 terminal and precursor
  // ladders; a stage whose sources are off is left out. Without the cascade, one stage (0) searches every spectrum with every source.
  // The stages are decided by the target-decoy FDR, so the cascade needs decoys.
  const bool cascade = search_cascade_ && ! native_rule && ! protein_level_only_ && decoy_factor_ > 0;
  if (search_cascade_ && ! cascade)
  {
    OPENMS_LOG_WARN << "search_cascade: not used with native, with protein_level_only or without decoys in the database" << std::endl;
  }
  std::vector<int> stages {0};
  if (cascade)
  {
    stages = {1};
    if (gapped_mode_ != "off") stages.push_back(2);
    if (terminal_ladders_ || use_precursor_ladders) stages.push_back(3);
  }
  // The candidates of @p spectra (indices into dspecs) from the sources of @p stage, appended to protein_hits_ and the ladder candidates
  // in spectrum order; the ranges of the other spectra are left as they are.
  auto searchStage = [&](int stage, const std::vector<int>& spectra) {
    const bool with_precursor_ladders = use_precursor_ladders && (stage == 0 || stage == 3);
    std::vector<String> native_notes(dspecs.size());
    // the precursor-ladder windows of all precursor masses searched below, in one pass over the database (the same windows as a sweep per mass)
    if (with_precursor_ladders)
    {
      std::vector<double> all_precursor_masses;
      for (int d : spectra)
      {
        const auto& dspec = dspecs[d];
        if (dspec.empty() || dspec.getPrecursorPeakGroup().getMonoMass() <= 0) continue;
        all_precursor_masses.push_back(dspec.getPrecursorPeakGroup().getMonoMass());
        if (decoy_precursor_ && dspec.getDecoyPrecursorMass() > 0) all_precursor_masses.push_back(dspec.getDecoyPrecursorMass());
      }
      mass_ladders.preparePrecursorWindows(all_precursor_masses, pl_precursor_tol);
    }
    const bool over_spectra = parallel_over_spectra;
#pragma omp parallel for schedule(dynamic) if (over_spectra)
    for (int k = 0; k < (int)spectra.size(); k++)
    {
      nextProgress();
      const int d = spectra[k];
      const auto tag_it = scan_to_tag_indices.find(dspecs[d].getScanNumber());
      searchCandidates(dspecs[d], tag_it == scan_to_tag_indices.end() ? nullptr : &tag_it->second, no_exclusion, with_precursor_ladders, stage,
                       spec_hits[d], spec_ladder_hits[d], spec_pl_hits[d], spec_decoy_pl_hits[d], spec_gapped_tags[d], spec_pl_windows[d],
                       native_notes[d]);
    }
    if (with_precursor_ladders) mass_ladders.clearPrecursorWindows();
    if (native_rule)
    {
      OPENMS_LOG_INFO << "Native candidates (score above the best decoy of the source in the spectrum + " << native_margin_ << "):";
      for (Size k = 0; k < native_sources.size(); k++)
      {
        OPENMS_LOG_INFO << (k ? ", " : " ") << native_sources[k] << " " << (Size)native_kept[k];
      }
      OPENMS_LOG_INFO << " (" << dspecs.size() << " spectra)" << std::endl;
      // the evidence per spectrum, for the few spectra of a native run
      for (Size d = 0; d < dspecs.size() && dspecs.size() <= 20; d++)
      {
        if (! native_notes[d].empty()) OPENMS_LOG_INFO << "  scan " << dspecs[d].getScanNumber() << native_notes[d] << std::endl;
      }
    }
    const Size pl_before = pl_candidates.size(), decoy_pl_before = decoy_pl_candidates.size();
    for (int d : spectra)
    {
      const int tag_base = appendSpectrumTags(spec_gapped_tags[d]);
      appendTagCandidates(spec_hits[d], protein_hits_, candidate_ranges[d], tag_base);
      appendLadderCandidates(spec_ladder_hits[d], ladder_candidates, ladder_ranges[d], tag_base);
      appendLadderCandidates(spec_pl_hits[d], pl_candidates, pl_ranges[d], tag_base);
      appendLadderCandidates(spec_decoy_pl_hits[d], decoy_pl_candidates, decoy_pl_ranges[d], tag_base);
      std::vector<FLASHHelperClasses::Tag>().swap(spec_gapped_tags[d]);
      std::vector<ProteinHit>().swap(spec_hits[d]);
      std::vector<ProteinHit>().swap(spec_ladder_hits[d]);
      std::vector<ProteinHit>().swap(spec_pl_hits[d]);
      std::vector<ProteinHit>().swap(spec_decoy_pl_hits[d]);
    }

    if (with_precursor_ladders)
    {
      std::vector<Size> windows; // per searched spectrum
      for (int d : spectra)
      {
        if (dspecs[d].getPrecursorPeakGroup().getMonoMass() > 0 && ! dspecs[d].empty()) windows.push_back(spec_pl_windows[d]);
      }
      Size total = 0;
      for (Size w : windows) total += w;
      std::nth_element(windows.begin(), windows.begin() + windows.size() / 2, windows.end());
      String at_decoy_masses;
      if (decoy_precursor_) at_decoy_masses = ", " + String(decoy_pl_candidates.size() - decoy_pl_before) + " at the decoy precursor masses";
      OPENMS_LOG_INFO << "Precursor ladders: " << mass_ladders.getSegmentCount() << " segments; " << windows.size() << " spectra with a precursor mass, " << total
                      << " windows (median " << (windows.empty() ? 0 : windows[windows.size() / 2]) << " per spectrum), " << pl_candidates.size() - pl_before
                      << " candidates (top " << min_protein_hits_ << ", at least " << MassLadderAlgorithm::precursor_ladder_min_k << " matched positions)" << at_decoy_masses
                      << std::endl;
    }
  };
  std::vector<int> all_spectra(dspecs.size());
  std::iota(all_spectra.begin(), all_spectra.end(), 0);
  searchStage(stages[0], all_spectra);
  endProgress();

  // the FDR below: decoys, their weights, the PrSM level of a passing PrSM and the FDR class of a PrSM; defined here for the search cascade,
  // which decides from each stage's PrSM q-values which spectra the next stage searches
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
  const UInt mi_chimeric_round = registry.registerName("ChimericRound");
  // the FDR class of a PrSM (index into fdr_class_names); -1: in no class
  auto fdrClassOf = [&](const ProteinHit& hit) {
    const int chimeric_round = hit.metaValueExists(mi_chimeric_round) ? (int)hit.getMetaValue(mi_chimeric_round) : 0;
    if (chimeric_round > 0)
    {
      // by source: a residual round of the chimeric search; undetermined ones only with keep_underdetermined_, as below
      if (! isDetermined(hit) && ! keep_underdetermined_) return -1;
      return fdr_class_residual + std::min(chimeric_round, max_chimeric_rounds) - 1;
    }
    if (hit.metaValueExists(mi_precursor_ladder)) return (int)fdr_class_precursor_ladder; // by source; all of them are determined
    if (isDetermined(hit)) return (int)std::min(blindModCount_(hit), (Size)fdr_class_max_blind);
    if (! keep_underdetermined_) return -1;
    const bool both_termini = (int)hit.getMetaValue("StartPosition") > 0 && (int)hit.getMetaValue("EndPosition") > 0;
    return (int)(both_termini ? fdr_class_both_termini : fdr_class_open_terminus);
  };

  // the order of the PrSMs for the FDR: best score first, then the higher coverage, then the higher scan number
  auto byScore = [](const ProteinHit& left, const ProteinHit& right) {
    return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getMetaValue("Scan") > right.getMetaValue("Scan"))
                                                                                            : (left.getCoverage() > right.getCoverage()))
                                               : (left.getScore() > right.getScore());
  };
  // The PrSM and proteoform q-values of proteoform_hits_ (sorted by byScore), as described at the FDR below: per PrSM its FDR class
  // (hit_class), PrSM q-value, decoy flag and decoy weight; the proteoforms (markRepresentativeProteoformHits_) are returned and every
  // PrSM gets its proteoform's q-value (proqvalue). With @p log, the classes are logged. Used by the FDR below and, on a copy, by the
  // search cascade between its stages.
  auto proteoformQvalues = [&](std::vector<int>& hit_class, std::vector<double>& qvalues, std::vector<bool>& decoys,
                               std::vector<double>& decoy_weights, bool log) {
    const Size hit_count = proteoform_hits_.size();
    hit_class.assign(hit_count, -1);
    qvalues.assign(hit_count, 1.0);
    decoy_weights.assign(hit_count, 0.0);
    decoys.assign(hit_count, false);
    if (use_fdr)
    {
      std::vector<std::vector<Size>> class_hits(fdr_class_count); // in the order above
      for (Size i = 0; i < hit_count; i++)
      {
        const auto& hit = proteoform_hits_[i];
        decoys[i] = is_decoy(hit);
        if (decoys[i]) decoy_weights[i] = decoy_weight(hit);
        const int c = fdrClassOf(hit);
        if (c < 0) continue;
        hit_class[i] = c;
        class_hits[c].push_back(i);
      }
      // per class: targets and decoys (count and weight), and the protein-level null it uses (its own from null_min_class_decoys decoys on,
      // else the pooled one); for the precursor-ladder class also by blind modification count, for a residual class by determined mass and
      // termini. The residual classes are logged only with the chimeric search.
      for (int c = 0; log && c < fdr_class_residual + chimeric_rounds_; c++)
      {
        Size targets = 0, decoy_count = 0, blind_targets[3] = {0, 0, 0}, blind_decoys[3] = {0, 0, 0}, determined_targets = 0, determined_decoys = 0;
        double weight = 0;
        for (Size i : class_hits[c])
        {
          const Size b = std::min(blindModCount_(proteoform_hits_[i]), (Size)2);
          const bool determined = isDetermined(proteoform_hits_[i]);
          if (decoys[i])
          {
            decoy_count++;
            weight += decoy_weights[i];
            blind_decoys[b]++;
            determined_decoys += determined;
          }
          else
          {
            targets++;
            blind_targets[b]++;
            determined_targets += determined;
          }
        }
        String line = String("FDR class \"") + fdr_class_names[c] + "\": " + String(targets) + " target and " + String(decoy_count)
                      + " decoy PrSMs (decoy weight " + String(weight) + "); protein-level null: "
                      + (c == fdr_class_precursor_ladder ? String("that of \"undetermined, both termini\" (setProteinQvalues_)")
                         : String(decoy_count >= null_min_class_decoys ? "own" : "pooled")
                             + (c <= fdr_class_max_blind ? "; those not the whole protein: that of \"undetermined, both termini\"" : ""));
        if (c == fdr_class_precursor_ladder && ! class_hits[c].empty())
        {
          line += "; target/decoy PrSMs by blind modifications: 0: " + String(blind_targets[0]) + "/" + String(blind_decoys[0]) + ", 1: "
                  + String(blind_targets[1]) + "/" + String(blind_decoys[1]) + ", 2 or more: " + String(blind_targets[2]) + "/" + String(blind_decoys[2]);
        }
        if (c >= fdr_class_residual && ! class_hits[c].empty())
        {
          line += "; target/decoy PrSMs determined: " + String(determined_targets) + "/" + String(determined_decoys) + ", undetermined: "
                  + String(targets - determined_targets) + "/" + String(decoy_count - determined_decoys);
        }
        OPENMS_LOG_INFO << line << std::endl;
      }
      const UInt mi_search_stage = registry.registerName("SearchStage");
      std::map<int, std::pair<Size, Size>> stage_counts; // search cascade: target and decoy PrSMs in a class per stage
      for (int c = 0; c < fdr_class_count; c++)
      {
        // the class's PrSMs by stage (0: no search cascade, or a residual round), each in the order above
        std::map<int, std::vector<Size>> by_stage;
        for (Size i : class_hits[c])
        {
          const int stage = proteoform_hits_[i].metaValueExists(mi_search_stage) ? (int)proteoform_hits_[i].getMetaValue(mi_search_stage) : 0;
          by_stage[stage].push_back(i);
          if (stage > 0) (decoys[i] ? stage_counts[stage].second : stage_counts[stage].first)++;
        }
        for (const auto& [stage, hits] : by_stage)
        {
          std::vector<RankedItem<double>> ranked;
          ranked.reserve(hits.size());
          for (Size i : hits) ranked.push_back({proteoform_hits_[i].getScore(), decoys[i], decoy_weights[i]});
          const auto class_qvalues = targetDecoyQvalues(ranked);
          for (Size k = 0; k < hits.size(); k++)
          {
            const Size i = hits[k];
            qvalues[i] = class_qvalues[k];
            proteoform_hits_[i].setMetaValue("qvalue", qvalues[i]);
            proteoform_hits_[i].setMetaValue("FDRClass", fdr_class_names[c]);
          }
        }
      }
      for (const auto& [stage, counts] : stage_counts)
      {
        if (! log) break;
        OPENMS_LOG_INFO << "Search cascade stage " << stage << ": " << counts.first << " target and " << counts.second
                        << " decoy PrSMs in FDR classes (q-values per class within the stage)" << std::endl;
      }
    }
    // the proteoforms: positions of their PrSMs in proteoform_hits_ (Index), the representative first; representatives in the order above
    auto clusters = markRepresentativeProteoformHits_(precursor_tol, hit_class, qvalues, rep_q, use_fdr);
    if (use_fdr)
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
        // Decoy precursors model wrong MS1 masses, which only a proteoform whose mass came from MS1 can have: they count against those, and a
        // proteoform without a mass or with a mass from fragment pairs takes its q-value from the same ranking without them (identical when
        // there are no decoy precursors).
        std::vector<RankedItem<std::pair<double, double>>> ranked_no_dp;
        std::vector<Size> no_dp_pos(order.size(), (Size)-1);
        for (Size k = 0; k < order.size(); k++)
        {
          const Size rep = clusters[order[k]][0];
          if (decoys[rep] && proteoform_hits_[rep].metaValueExists("DecoyPrecursor") && ! proteoform_hits_[rep].getAccession().hasPrefix("DECOY")) continue;
          no_dp_pos[k] = ranked_no_dp.size();
          ranked_no_dp.push_back(ranked[k]);
        }
        const auto proteoform_qvalues_no_dp = ranked_no_dp.size() < ranked.size() ? targetDecoyQvalues(ranked_no_dp) : proteoform_qvalues;
        for (Size k = 0; k < order.size(); k++)
        {
          const Size rep = clusters[order[k]][0];
          const bool ms1_mass = (double)proteoform_hits_[rep].getMetaValue("Mass", -1.0) > 0
                                && (int)proteoform_hits_[rep].getMetaValue("ProteoformMassByFragmentMass", 0) == 0;
          const double pq = (ms1_mass || no_dp_pos[k] == (Size)-1) ? proteoform_qvalues[k] : proteoform_qvalues_no_dp[ranked_no_dp.size() < ranked.size() ? no_dp_pos[k] : k];
          for (Size i : clusters[order[k]]) proteoform_hits_[i].setMetaValue("proqvalue", pq);
        }
      }
    return clusters;
  };

  if (! protein_level_only_)
  {
    // Beyond the min_protein_hits best candidates of a spectrum (a fair target-decoy competition: they are chosen by this spectrum
    // alone), extend only those whose protein, or its target-decoy partner, passes the run-level protein FDR. That needs decoys;
    // without them every candidate is extended, as before.
    // In the search cascade, it is recomputed before each stage's extension from the tag candidates of the stages so far.
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
    const bool keep_all_candidates = (bool)raw_loader_;
    // precursor-ladder PrSMs of the extension that is kept: kept (determined; of them at both window ends, at the window start only, at the
    // window end only), dropped as the same proteoform as a determined PrSM of another source, dropped as undetermined
    std::atomic<Size> pl_determined {0}, pl_at_both_ends {0}, pl_at_start {0}, pl_at_end {0}, pl_duplicates {0}, pl_undetermined {0};
    // The extension of one spectrum's candidates with the given extension parameters: the sequence- and gapped-tag candidates
    // tag_candidates[tag_range], the terminal-ladder candidates ladder_candidates_now[ladder_range] and the precursor-ladder candidates
    // pl_candidates_now[pl_range]. With keep_all, every candidate's PrSM is kept (multiple_hits_per_spec_, raw-spectrum scoring, the
    // residual rounds of the chimeric search); otherwise the extension keeps one.
    auto extendCandidatesInto = [&](const Param& eparam, const DeconvolvedSpectrum& dspec, const std::vector<ProteinHit>& tag_candidates,
                                    std::pair<Size, Size> tag_range, const std::vector<ProteinHit>& ladder_candidates_now, std::pair<Size, Size> ladder_range,
                                    const std::vector<ProteinHit>& pl_candidates_now, std::pair<Size, Size> pl_range, std::vector<ProteinHit>& out,
                                    bool keep_all) {
      out.clear();
      const auto [first, second] = tag_range;
      const auto [pl_first, pl_second] = pl_range;
      if (first == second && ladder_range.first == ladder_range.second && pl_first == pl_second) return;

      std::vector<ProteinHit> hits;
      // the min_protein_hits best candidates of each origin (sequence tags, gapped tags) are always extended, so gapped candidates
      // add chances without taking the places of sequence-tag candidates
      int rank_tag = 0, rank_gapped = 0;
      for (Size i = first; i < second; i++)
      {
        const bool gapped = tag_candidates[i].metaValueExists(mi_gapped);
        // native: the decoy margin chose the candidates
        bool top = (gapped ? rank_gapped++ : rank_tag++) < min_protein_hits_ || native_rule;
        if (use_protein_fdr && ! top && passing_proteins.count((Size)(int)tag_candidates[i].getMetaValue(mi_fasta_index)) == 0) continue;
        hits.push_back(tag_candidates[i]);
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
      std::vector<ProteinHit> ladder_hits(ladder_candidates_now.begin() + ladder_range.first, ladder_candidates_now.begin() + ladder_range.second);
      auto extend = [&](std::vector<ProteinHit>& candidates, std::vector<ProteinHit>& proteoforms) {
        proteoforms.clear();
        if (candidates.empty()) return;
        FLASHExtenderAlgorithm extender;
        extender.setParameters(eparam);
        extender.setCandidateBlindModificationMap(blind_mod_map);
        // the chimeric search removes the masses of a spectrum's best PrSM before the next round
        extender.setReportMatchedMasses(chimeric_rounds_ > 0);
        // with raw-spectrum scoring, every candidate's PrSM is kept here and the best one per spectrum is chosen by the combined score in
        // annotateOnRawSpectra_: the extension score alone (deconvolved masses) picked another protein in 1 of 15 spectra of Eclip_1480
        // where TopPIC's protein had the better raw evidence (offline: 5252 -> 5451 PrSMs, 1435 -> 1500 proteoforms at 1%)
        extender.run(candidates, dspec, spec_vec, vec_pro, rev_vec_pro, tags_, tol, keep_all);
        extender.fillProteoforms(proteoforms);
      };
      extend(hits, out);
      // native: every candidate's PrSM is reported, one per protein (the better one)
      auto mergeNative = [&out](std::vector<ProteinHit>& from) {
        for (auto& h : from)
        {
          auto it = std::find_if(out.begin(), out.end(), [&h](const ProteinHit& o) { return o.getAccession() == h.getAccession(); });
          if (it == out.end()) { out.push_back(std::move(h)); }
          else if (h.getScore() > it->getScore()) { *it = std::move(h); }
        }
      };
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
          h.setMetaValue(native_rule ? mi_sliding : mi_terminal, 1);
          (is_determined(h) ? determined : undetermined).push_back(h);
        }
        const bool fill = out.empty();
        if (native_rule)
        {
          mergeNative(determined);
          mergeNative(undetermined);
        }
        else if (keep_all) { out.insert(out.end(), determined.begin(), determined.end()); }
        else if (! determined.empty() && (out.empty() || ! is_determined(out[0]) || determined[0].getScore() > out[0].getScore()))
        {
          // a ladder PrSM that explains the precursor mass also replaces an undetermined one of a higher score
          out.swap(determined);
        }
        if (native_rule) {}
        else if (fill && out.empty()) out.swap(undetermined);
        else if (fill && keep_all) { out.insert(out.end(), undetermined.begin(), undetermined.end()); }
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
        if (native_rule) { mergeNative(pl_kept); }
        else if (keep_all) { out.insert(out.end(), pl_kept.begin(), pl_kept.end()); }
        else if (! pl_kept.empty() && (out.empty() || pl_kept[0].getScore() > out[0].getScore())) { out.swap(pl_kept); }
      }
    };
    // The extension of spectrum d's candidates. With decoy_mass, the spectrum is the decoy copy at the decoy precursor mass
    // (decoy_precursor), and its precursor-ladder candidates are those of the windows that fit that mass.
    auto extendSpectrumInto = [&](int d, const Param& eparam, const DeconvolvedSpectrum& dspec, std::vector<ProteinHit>& out, bool decoy_mass) {
      extendCandidatesInto(eparam, dspec, protein_hits_, candidate_ranges[d], ladder_candidates, ladder_ranges[d],
                           decoy_mass ? decoy_pl_candidates : pl_candidates, decoy_mass ? decoy_pl_ranges[d] : pl_ranges[d], out,
                           multiple_hits_per_spec_ || keep_all_candidates);
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

    // The PrSMs of the spectra of one stage (@p spectra; the first stage: every spectrum), appended to proteoform_hits_ in spectrum order.
    // Main pass of the first stage: with learned penalties every spectrum, otherwise the spectra not extended yet with the first-pass
    // penalties; a later stage of the search cascade extends its spectra with the penalties in use. The precursor-ladder counts describe
    // the extension that is kept: they restart before a pass that replaces the pilot's PrSMs, and otherwise the counts of the passes add up
    // to those of a single pass.
    const bool penalties_learned = ! learned_extender_param.empty();
    const Param& eparam_in_use = penalties_learned ? learned_extender_param : extender_param_;
    const UInt mi_stage = registry.registerName("SearchStage");
    auto extendStage = [&](int stage, const std::vector<int>& spectra) {
      const bool first_stage = stage == stages[0];
      if (! first_stage)
      {
        startProgress(0, (SignedSize)spectra.size(), "Running extension algorithm (search cascade stage " + String(stage) + ") ...");
        extendSpectra(spectra, eparam_in_use);
      }
      else
      {
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
        extendSpectra(main_spectra, eparam_in_use);
      }
      const bool extend_over_spectra = overSpectra(spectra);
      if (use_precursor_ladders && (stage == 0 || stage == 3))
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
        const Param& eparam_now = eparam_in_use;
        std::vector<std::vector<ProteinHit>> decoy_proteoforms(dspecs.size());
        int copies = 0, without_decoy = 0;
        std::vector<int> todo;
        for (int d : spectra)
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
        if (copies == 0 && first_stage)
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
        // Per spectrum, the decoy precursor competes with the spectrum's PrSMs that took their mass from MS1, on the precursor: the decoy
        // spectrum has to give a PrSM that reports the decoy mass (the try), and that decoy wins when its precursor score (FLASHDeconv's 2D
        // Qscore, the score its own FDR compares target and decoy masses with) is higher than the precursor's. The fragments cannot decide
        // the mass: a related or wrong mass explains nearly the same ladders, only the MS1 signal tells them apart. Inputs from FLASHDeconv
        // versions that do not write the decoy precursor's score fall back to the fragment score (the higher PrSM score wins).
        int decoy_wins = 0, decoy_tries = 0, decoy_by_fragments = 0;
        for (int d : todo)
        {
          if (decoy_proteoforms[d].empty()) continue;
          decoy_tries++;
          double best_decoy = -1e10, best_own = -1e10;
          for (const auto& h : decoy_proteoforms[d]) best_decoy = std::max(best_decoy, h.getScore());
          for (const auto& h : spec_proteoforms[d])
          {
            if (ms1Based(h)) best_own = std::max(best_own, h.getScore());
          }
          // on the Qscores in the survey scan alone (1D) when the input has them: a decoy is rarely traced over several scans, so the 2D
          // Qscore favours the precursor (Eclip_1480: of 137 decoy spectra giving a PrSM with the decoy mass, 12 decoys beat the precursor on
          // the 1D Qscore, 6 on the 2D one)
          const double decoy_qscore1d = dspecs[d].getDecoyPrecursorQscore1D(), decoy_qscore = dspecs[d].getDecoyPrecursorQscore();
          if (decoy_qscore1d >= 0)
          {
            if (decoy_qscore1d <= dspecs[d].getPrecursorPeakGroup().getQscore()) continue;
          }
          else if (decoy_qscore >= 0)
          {
            if (decoy_qscore <= dspecs[d].getPrecursorPeakGroup().getQscore2D()) continue;
          }
          else
          {
            decoy_by_fragments++;
            if (best_decoy <= best_own) continue;
          }
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
                        << " without one); the decoy spectrum gave a PrSM with the decoy mass in " << decoy_tries << ", the decoy precursor wins in "
                        << decoy_wins << (decoy_by_fragments > 0 ? " (" + String(decoy_by_fragments) + " decided on the fragments: no decoy precursor score in the input)" : String(""))
                        << std::endl;
      }

      std::vector<ProteinHit> stage_hits;
      for (int d : spectra)
      {
        for (auto& h : spec_proteoforms[d]) stage_hits.push_back(std::move(h));
        std::vector<ProteinHit>().swap(spec_proteoforms[d]);
      }
      endProgress();

      // raw-spectrum evidence of every PrSM, targets and decoys alike, before the FDR
      if (raw_loader_ && ! stage_hits.empty())
      {
        std::set<int> scans;
        for (const auto& hit : stage_hits) scans.insert((int)hit.getMetaValue("Scan"));
        loadRawSpectra_(scans);
        annotateOnRawSpectra_(stage_hits);
        // one PrSM per spectrum: the best by the combined score; the extension kept every candidate's PrSM for this choice
        if (! multiple_hits_per_spec_) keepBestPerSpectrum_(stage_hits, stage == 0 ? String("") : " (search cascade stage " + String(stage) + ")", "combined score");
      }
      for (auto& h : stage_hits)
      {
        if (stage > 0) h.setMetaValue(mi_stage, stage);
        proteoform_hits_.push_back(std::move(h));
      }
    };

    // The spectra identified after a stage of the search cascade: with a target PrSM, of any stage so far, whose proteoform q-value is at
    // most rep_q, the level at which proteoforms are reported. The q-values are those the FDR below gives the PrSMs so far (per class and
    // stage, the proteoforms over all of them), computed on a copy; the one proteoform per MS1 feature rule is applied at the end only.
    // An identified spectrum is searched by no later stage.
    std::map<int, int> spectrum_of_scan_number;
    for (int d = 0; d < (int)dspecs.size(); d++) spectrum_of_scan_number[dspecs[d].getScanNumber()] = d;
    auto identifiedSoFar = [&]() {
      std::vector<ProteinHit> kept = proteoform_hits_;
      std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), byScore);
      std::vector<int> hit_class;
      std::vector<double> qvalues, decoy_weights;
      std::vector<bool> decoys;
      proteoformQvalues(hit_class, qvalues, decoys, decoy_weights, false);
      std::set<int> identified;
      for (Size i = 0; i < proteoform_hits_.size(); i++)
      {
        const auto& hit = proteoform_hits_[i];
        if (decoys[i] || hit_class[i] < 0 || (double)hit.getMetaValue("proqvalue", 1.0) > rep_q) continue;
        const auto it = spectrum_of_scan_number.find((int)hit.getMetaValue(mi_scan));
        if (it != spectrum_of_scan_number.end()) identified.insert(it->second);
      }
      proteoform_hits_.swap(kept);
      return identified;
    };

    std::vector<int> stage_spectra = all_spectra;
    for (Size s = 0; s < stages.size(); s++)
    {
      const int stage = stages[s];
      if (s > 0)
      {
        startProgress(0, (SignedSize)stage_spectra.size(), "Running candidate protein filtration (search cascade stage " + String(stage) + ") ...");
        searchStage(stage, stage_spectra);
        endProgress();
        if (use_protein_fdr) passing_proteins = getProteinsPassingFDR_(fasta_entry);
      }
      const Size from = proteoform_hits_.size();
      extendStage(stage, stage_spectra);
      if (stage == 0) continue;
      if (s + 1 == stages.size())
      {
        OPENMS_LOG_INFO << "Search cascade stage " << stage << " (" << (stage == 1 ? "sequence tags" : stage == 2 ? "gapped tags" : "terminal and precursor ladders")
                        << "): " << stage_spectra.size() << " spectra searched, " << proteoform_hits_.size() - from << " PrSMs" << std::endl;
        continue;
      }
      const auto identified = identifiedSoFar();
      std::vector<int> rest;
      for (int d : stage_spectra)
      {
        if (identified.count(d) == 0) rest.push_back(d);
      }
      OPENMS_LOG_INFO << "Search cascade stage " << stage << " (" << (stage == 1 ? "sequence tags" : "gapped tags") << "): " << stage_spectra.size()
                      << " spectra searched, " << proteoform_hits_.size() - from << " PrSMs, " << stage_spectra.size() - rest.size()
                      << " of the spectra identified (a target PrSM at proteoform q-value <= " << rep_q << "); the next stage searches the other "
                      << rest.size() << std::endl;
      stage_spectra.swap(rest);
    }
    // the chimeric rounds score their PrSMs on the raw spectra of every stage
    if (cascade && raw_loader_ && chimeric_rounds_ > 0 && ! proteoform_hits_.empty())
    {
      std::set<int> scans;
      for (const auto& hit : proteoform_hits_) scans.insert((int)hit.getMetaValue("Scan"));
      loadRawSpectra_(scans);
    }

    // Chimeric spectra (chimeric_rounds): a spectrum of two or more co-isolated proteoforms is searched again without the masses its best
    // PrSM explains. In the full spectrum, a second proteoform competes with chance PrSMs that take many of their masses from the first one
    // (native prototype, final_all/chimeric: 55% of the masses of chance PrSMs were the first proteoform's; 6% in the residual), so its
    // per-spectrum null is lower in the residual. Every spectrum with a PrSM is searched again, whether its best PrSM is a target or a
    // decoy, so that the residual PrSMs of targets and decoys come from the same procedure; they form an FDR class per round. Per round:
    // - removed: the masses on the path of the spectrum's best PrSM of the last round (the extension's MatchedMassIndices; with
    //   only_single_hit its only PrSM, else the best by score, then a determined mass and termini, then coverage), each also up to
    //   chimeric_isotope_offsets isotope spacings away, and its proteoform mass, all within the spectrum's tolerance. They are removed from
    //   all masses of the spectrum, and the cap max_node_cntr_ is applied again to the rest, as to the input.
    // - searched: tags of the residual (a tag found before in the scan keeps its index; new tag sequences are looked up in the database),
    //   the candidates of searchCandidates without precursor ladders and without the proteins of the spectrum's earlier best PrSMs and
    //   their target-decoy partners, and the extension without a precursor mass (the MS1 mass is the first proteoform's; a mass can still
    //   come from complementary fragment pairs). Every candidate's PrSM is kept, scored on the raw spectrum with -in_raw, and with
    //   only_single_hit the best one per spectrum is reported (by the combined score with -in_raw, as in the first search).
    // - reported per residual PrSM: ChimericRound, RemovedMasses (removed from the spectrum before this round), ResidualNullScore (the
    //   best extension score of the other decoy PrSMs of the same residual search; -1: none) and ResidualScoreRatio (extension score /
    //   ResidualNullScore; -1: no null). They are not a filter: the FDR is the target-decoy competition of the round's class, and the
    //   native prototype's rule (ratio >= 2.5, >= 8 matched masses) was calibrated with a stronger null (also three shuffled-only
    //   databases), so it is applied by the user, e.g. to a single native spectrum, where no class has enough PrSMs for a q-value.
    if (chimeric_rounds_ > 0 && ! proteoform_hits_.empty())
    {
      const Param& eparam_now = learned_extender_param.empty() ? extender_param_ : learned_extender_param;
      const UInt mi_round = registry.registerName("ChimericRound"), mi_matched_masses = registry.registerName("MatchedMassIndices"),
                 mi_removed = registry.registerName("RemovedMasses"), mi_null = registry.registerName("ResidualNullScore"),
                 mi_ratio = registry.registerName("ResidualScoreRatio"), mi_mass = registry.registerName("Mass");
      for (auto& hit : proteoform_hits_) hit.setMetaValue(mi_round, 0);
      std::map<int, int> spectrum_of_scan;
      for (int d = 0; d < (int)dspecs.size(); d++) spectrum_of_scan[dspecs[d].getScanNumber()] = d;
      // the target-decoy partner of a protein (fasta index; -1: none), named as in getProteinsPassingFDR_
      std::unordered_map<std::string, Size> index_of;
      index_of.reserve(fasta_entry.size());
      for (Size i = 0; i < fasta_entry.size(); i++) index_of.emplace(fasta_entry[i].identifier, i);
      auto partnerOf = [&](Size fi) -> long {
        const String& id = fasta_entry[fi].identifier;
        std::vector<std::string> names;
        if (id.hasPrefix("DECOY_")) { names.push_back(id.substr(6)); }
        else if (id.hasPrefix("DECOY")) { names.push_back(id.substr(5)); }
        else { names = {"DECOY_" + id, "DECOY" + id}; }
        for (const auto& name : names)
        {
          const auto it = index_of.find(name);
          if (it != index_of.end()) return (long)it->second;
        }
        return -1;
      };
      const double isotope_spacing = Constants::ISOTOPE_MASSDIFF_55K_U; // FLASHDeconv's isotope spacing
      // per spectrum: all masses (ascending) and the searched masses (after the cap) of its last residual, the proteins excluded and the
      // number of masses removed so far
      std::vector<std::vector<PeakGroup>> last_all(dspecs.size());
      std::vector<DeconvolvedSpectrum> last_searched(dspecs.size());
      std::vector<std::set<Size>> excluded(dspecs.size());
      std::vector<Size> removed_count(dspecs.size(), 0);
      // the sequence tags found so far in each scan; a residual tag equal to one of them keeps its index
      std::map<int, std::vector<Size>> tags_of_scan = scan_to_tag_indices;
      const std::vector<ProteinHit> no_candidates;
      Size last_first = 0, last_end = proteoform_hits_.size(); // the PrSMs of the last round in proteoform_hits_
      for (int round = 1; round <= chimeric_rounds_; round++)
      {
        // the best PrSM of each spectrum in the last round, targets and decoys alike
        std::vector<long> best(dspecs.size(), -1);
        for (Size i = last_first; i < last_end; i++)
        {
          const auto it = spectrum_of_scan.find((int)proteoform_hits_[i].getMetaValue(mi_scan));
          if (it == spectrum_of_scan.end()) continue;
          long& b = best[it->second];
          if (b < 0 || isBetterPrSM(proteoform_hits_[i], proteoform_hits_[(Size)b])) b = (long)i;
        }
        // the residual spectra
        std::vector<int> searched; // ascending
        std::vector<DeconvolvedSpectrum> residual(dspecs.size());
        std::vector<std::vector<PeakGroup>> residual_all(dspecs.size());
        Size removed_total = 0;
        for (int d = 0; d < (int)dspecs.size(); d++)
        {
          if (best[d] < 0) continue;
          const auto& best_hit = proteoform_hits_[(Size)best[d]];
          const DeconvolvedSpectrum& last = round == 1 ? dspecs[d] : last_searched[d];
          std::vector<PeakGroup> all; // the last round's spectrum with the masses its cap left out
          if (round > 1) all = std::move(last_all[d]);
          else if (! uncapped_masses[d].empty()) all = std::move(uncapped_masses[d]);
          else all.assign(last.begin(), last.end());
          std::vector<double> centers; // the masses on the path of the best PrSM, and its proteoform mass
          if (best_hit.metaValueExists(mi_matched_masses))
          {
            for (int i : best_hit.getMetaValue(mi_matched_masses).toIntList())
            {
              if (i >= 0 && i < (int)last.size()) centers.push_back(last[i].getMonoMass());
            }
          }
          const double proteoform_mass = best_hit.getMetaValue(mi_mass);
          if (proteoform_mass > 0) centers.push_back(proteoform_mass);
          std::vector<double> removed_masses;
          removed_masses.reserve(centers.size() * (2 * chimeric_isotope_offsets + 1));
          for (double c : centers)
          {
            for (int k = -chimeric_isotope_offsets; k <= chimeric_isotope_offsets; k++) removed_masses.push_back(c + k * isotope_spacing);
          }
          std::sort(removed_masses.begin(), removed_masses.end());
          std::vector<PeakGroup> kept;
          kept.reserve(all.size());
          for (const auto& pg : all)
          {
            const double m = pg.getMonoMass(), margin = spectrum_tols[d] * 1e-6 * m;
            const auto it = std::lower_bound(removed_masses.begin(), removed_masses.end(), m - margin);
            if (it != removed_masses.end() && *it <= m + margin) continue;
            kept.push_back(pg);
          }
          removed_count[d] += all.size() - kept.size();
          removed_total += all.size() - kept.size();
          // the protein of the best PrSM and its target-decoy partner are not proposed again in this spectrum; in native mode, where every
          // kept candidate's PrSM is reported, none of the proteins reported in the spectrum (one PrSM per protein)
          std::vector<Size> reported {(Size)(int)best_hit.getMetaValue(mi_fasta_index)};
          if (native_)
          {
            for (Size i = last_first; i < last_end; i++)
            {
              const auto& h = proteoform_hits_[i];
              const auto it = spectrum_of_scan.find((int)h.getMetaValue(mi_scan));
              if (it != spectrum_of_scan.end() && it->second == d && h.metaValueExists(mi_fasta_index)) reported.push_back((Size)(int)h.getMetaValue(mi_fasta_index));
            }
          }
          for (Size fasta_index : reported)
          {
            excluded[d].insert(fasta_index);
            if (const long partner = partnerOf(fasta_index); partner >= 0) excluded[d].insert((Size)partner);
          }
          // the residual: no precursor mass, and of the original spectrum only what the search reads (retention time, MS level)
          DeconvolvedSpectrum res(dspecs[d].getScanNumber());
          MSSpectrum original;
          original.setRT(dspecs[d].getOriginalSpectrum().getRT());
          original.setMSLevel(dspecs[d].getOriginalSpectrum().getMSLevel());
          res.setOriginalSpectrum(original);
          res.setPeakGroups(kept);
          if (res.size() > max_node_cntr_) // the cap of the input
          {
            res.sortByQscore();
            std::vector<PeakGroup> filtered_peaks;
            filtered_peaks.reserve(max_node_cntr_);
            for (const auto& pg : res)
            {
              if (filtered_peaks.size() >= max_node_cntr_) break;
              filtered_peaks.push_back(pg);
            }
            res.setPeakGroups(filtered_peaks);
          }
          res.sort();
          residual[d] = std::move(res);
          residual_all[d] = std::move(kept);
          searched.push_back(d);
        }

        // tags of the residual spectra, in spectrum order; new tag sequences are looked up in the database
        std::vector<std::vector<Size>> residual_tags(dspecs.size());
        std::vector<std::string> new_seqs;
        Size tag_count = 0, new_tag_count = 0;
        for (int d : searched)
        {
          FLASHTaggerAlgorithm tagger;
          tagger.setParameters(tagger_param_);
          tagger.run(residual[d], spectrum_tols[d]);
          std::vector<FLASHHelperClasses::Tag> found;
          tagger.fillTags(found);
          tag_count += found.size();
          auto& known = tags_of_scan[residual[d].getScanNumber()];
          for (auto& tag : found)
          {
            Size index = tags_.size();
            for (Size ti : known)
            {
              if (tags_[ti] == tag && tags_[ti].getMzs() == tag.getMzs())
              {
                index = ti;
                break;
              }
            }
            if (index == tags_.size())
            {
              const String seq = tag.getUppercaseSequence();
              if (tag_to_protein_indices.find(seq) == tag_to_protein_indices.end())
              {
                tag_to_protein_indices[seq] = std::vector<Size>();
                tag_to_protein_positions[seq] = std::vector<Size>();
                tag_to_n_flanking_masses[seq] = std::vector<double>();
                tag_to_c_flanking_masses[seq] = std::vector<double>();
                new_seqs.push_back(seq);
              }
              tag.setIndex((int)index);
              tags_.push_back(tag);
              known.push_back(index);
              new_tag_count++;
            }
            residual_tags[d].push_back(index);
          }
        }
        if (! new_seqs.empty()) searchTagSequences(new_seqs);

        // candidates of the residual spectra (spectra in parallel, each writing its own slots), appended in spectrum order
        std::vector<std::vector<ProteinHit>> r_hits(dspecs.size()), r_ladder_hits(dspecs.size()), r_pl_hits(dspecs.size()), r_decoy_pl_hits(dspecs.size());
        std::vector<std::vector<FLASHHelperClasses::Tag>> r_gapped_tags(dspecs.size());
        std::vector<Size> r_windows(dspecs.size(), 0);
#ifdef _OPENMP
        const bool residual_over_spectra = (int)searched.size() >= omp_get_max_threads();
#else
        const bool residual_over_spectra = false;
#endif
#pragma omp parallel for schedule(dynamic) if (residual_over_spectra)
        for (int k = 0; k < (int)searched.size(); k++)
        {
          const int d = searched[k];
          String note; // the native note of the first search is the one logged
          searchCandidates(residual[d], &residual_tags[d], excluded[d], false, 0, r_hits[d], r_ladder_hits[d], r_pl_hits[d], r_decoy_pl_hits[d],
                           r_gapped_tags[d], r_windows[d], note);
        }
        std::vector<ProteinHit> r_candidates, r_ladder_candidates;
        std::vector<std::pair<Size, Size>> r_ranges(dspecs.size(), {0, 0}), r_ladder_ranges(dspecs.size(), {0, 0});
        int with_candidates = 0;
        for (int d : searched)
        {
          const int tag_base = appendSpectrumTags(r_gapped_tags[d]);
          appendTagCandidates(r_hits[d], r_candidates, r_ranges[d], tag_base);
          appendLadderCandidates(r_ladder_hits[d], r_ladder_candidates, r_ladder_ranges[d], tag_base);
          with_candidates += r_ranges[d].first != r_ranges[d].second || r_ladder_ranges[d].first != r_ladder_ranges[d].second;
        }

        // the extension of the residual spectra, every candidate's PrSM kept
        std::vector<std::vector<ProteinHit>> r_prsms(dspecs.size());
#ifdef _OPENMP
        const bool extend_residual_over_spectra = with_candidates >= omp_get_max_threads();
#else
        const bool extend_residual_over_spectra = false;
#endif
#pragma omp parallel for schedule(dynamic) if (extend_residual_over_spectra)
        for (int k = 0; k < (int)searched.size(); k++)
        {
          const int d = searched[k];
          extendCandidatesInto(eparam_now, residual[d], r_candidates, r_ranges[d], r_ladder_candidates, r_ladder_ranges[d], no_candidates,
                               std::pair<Size, Size>(0, 0), r_prsms[d], true);
        }

        // the null of each PrSM: the best extension score of the other decoy PrSMs of the same residual search
        std::vector<ProteinHit> round_hits;
        for (int d : searched)
        {
          auto& prsms = r_prsms[d];
          double best_decoy = -1, second_decoy = -1;
          long best_decoy_at = -1;
          for (Size i = 0; i < prsms.size(); i++)
          {
            if (! prsms[i].getAccession().hasPrefix("DECOY")) continue;
            const double score = prsms[i].getScore();
            if (score > best_decoy)
            {
              second_decoy = best_decoy;
              best_decoy = score;
              best_decoy_at = (long)i;
            }
            else if (score > second_decoy) { second_decoy = score; }
          }
          for (Size i = 0; i < prsms.size(); i++)
          {
            auto& h = prsms[i];
            h.setMetaValue(mi_round, round);
            h.setMetaValue(mi_removed, (int)removed_count[d]);
            h.setMetaValue(mi_null, (long)i == best_decoy_at ? second_decoy : best_decoy);
            round_hits.push_back(std::move(h));
          }
        }
        const Size all_prsms = round_hits.size();
        if (raw_loader_ && ! round_hits.empty()) annotateOnRawSpectra_(round_hits);
        if (! multiple_hits_per_spec_) keepBestPerSpectrum_(round_hits, " (chimeric round " + String(round) + ")", raw_loader_ ? "combined score" : "score");
        Size targets = 0, decoys = 0;
        for (auto& h : round_hits)
        {
          const double extension_score = h.metaValueExists("ExtensionScore") ? (double)h.getMetaValue("ExtensionScore") : h.getScore();
          const double null_score = h.getMetaValue(mi_null);
          h.setMetaValue(mi_ratio, null_score > 0 ? extension_score / null_score : -1.0);
          (h.getAccession().hasPrefix("DECOY") ? decoys : targets)++;
        }
        OPENMS_LOG_INFO << "Chimeric round " << round << ": " << searched.size() << " spectra searched without the masses of their best PrSM ("
                        << removed_total << " masses removed), " << tag_count << " tags (" << new_tag_count << " new, " << new_seqs.size()
                        << " new sequences), " << r_candidates.size() << " tag and " << r_ladder_candidates.size() << " ladder candidates, "
                        << all_prsms << " PrSMs of all candidates; reported: " << targets << " target and " << decoys << " decoy PrSMs" << std::endl;

        last_first = proteoform_hits_.size();
        for (auto& h : round_hits) proteoform_hits_.push_back(std::move(h));
        last_end = proteoform_hits_.size();
        for (int d : searched)
        {
          last_searched[d] = std::move(residual[d]);
          last_all[d] = std::move(residual_all[d]);
        }
      }
    }
    releaseRawSpectra_();
  }

  // the order below and every ranking of the FDR need finite scores: a NaN breaks the sorts and never closes a tie group, and an infinite
  // decoy score leaves the null of the protein level undefined
  for (const auto& hit : proteoform_hits_)
  {
    if (std::isfinite(hit.getScore())) continue;
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "PrSM score must be finite (scan " + hit.getMetaValue("Scan").toString() + ")", std::to_string(hit.getScore()));
  }
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), byScore);

  // One MS1 feature, one proteoform: a co-isolated precursor is one real mass of one proteoform. When PrSMs of different proteins report
  // the same MS1-based mass in the same elution (within 2 min), FLASHExtender fitted at least one of them to it by truncation or blind
  // shifts (Eclip_1480: 55 PrSMs at 1% reported a mass of another protein's 1% proteoform nearby and never their own, against 8 with the
  // retention times moved by 20 min; e.g. TATA 41-89 +16.02 and MALE 328-373 -90.1 +473.1 both at 5549.62 in scans 2325 and 2326). Per
  // group of such PrSMs, each protein is ranked by its PrSM with the fewest blind modifications, then by its best score, then by its
  // number of PrSMs at the mass; the first keeps the mass, the others' PrSMs lose it (Mass -1: undetermined, their own FDR class).
  if (one_proteoform_per_feature_ && precursor_tol > 0)
  {
    const double iso = Constants::ISOTOPE_MASSDIFF_55K_U, rt_window = 120.0;
    std::vector<Size> idx;
    for (Size i = 0; i < proteoform_hits_.size(); i++)
    {
      const auto& h = proteoform_hits_[i];
      if ((double)h.getMetaValue("Mass") <= 0 || (int)h.getMetaValue("ProteoformMassByFragmentMass", 0) != 0 || h.metaValueExists("DecoyPrecursor")) continue;
      idx.push_back(i);
    }
    std::sort(idx.begin(), idx.end(), [&](Size l, Size r) { return (double)proteoform_hits_[l].getMetaValue("Mass") < (double)proteoform_hits_[r].getMetaValue("Mass"); });
    // conflicts: pairs of PrSMs of different proteins at the same mass (or one isotope apart) and elution; union-find over them
    std::vector<Size> parent(idx.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::function<Size(Size)> root = [&](Size a) { return parent[a] == a ? a : parent[a] = root(parent[a]); };
    std::vector<bool> in_conflict(idx.size(), false);
    for (Size a = 0; a < idx.size(); a++)
    {
      const auto& ha = proteoform_hits_[idx[a]];
      const double ma = ha.getMetaValue("Mass"), ra = ha.getMetaValue("RT"), max_diff = ma * precursor_tol * 2e-6;
      for (Size b = a + 1; b < idx.size(); b++)
      {
        const auto& hb = proteoform_hits_[idx[b]];
        const double mb = hb.getMetaValue("Mass");
        if (mb - ma > iso + max_diff) break;
        if (std::abs((double)hb.getMetaValue("RT") - ra) > rt_window) continue;
        const double diff = mb - ma;
        if (std::abs(diff) > max_diff && std::abs(diff - iso) > max_diff) continue;
        // the same protein is the same proteoform here; a target and its decoy are different proteins
        if (ha.getAccession() == hb.getAccession()) { parent[root(a)] = root(b); continue; }
        in_conflict[a] = in_conflict[b] = true;
        parent[root(a)] = root(b);
      }
    }
    std::map<Size, std::vector<Size>> groups;
    for (Size a = 0; a < idx.size(); a++) groups[root(a)].push_back(a);
    Size demoted = 0, groups_with_conflict = 0;
    for (const auto& [g, members] : groups)
    {
      if (std::none_of(members.begin(), members.end(), [&](Size a) { return (bool)in_conflict[a]; })) continue;
      groups_with_conflict++;
      // per protein: fewest blind modifications, best score, PrSM count
      std::map<String, std::tuple<Size, double, Size>> by_protein;
      for (Size a : members)
      {
        const auto& h = proteoform_hits_[idx[a]];
        auto it = by_protein.find(h.getAccession());
        const Size blind = blindModCount_(h);
        if (it == by_protein.end()) { by_protein[h.getAccession()] = {blind, h.getScore(), 1}; continue; }
        auto& [b0, s0, n0] = it->second;
        b0 = std::min(b0, blind);
        s0 = std::max(s0, h.getScore());
        n0++;
      }
      String owner;
      std::tuple<Size, double, Size> best {};
      bool first = true;
      for (const auto& [acc, key] : by_protein)
      {
        const auto& [b, s, n] = key;
        const auto& [bb, bs, bn] = best;
        if (first || b < bb || (b == bb && (s > bs || (s == bs && n > bn))))
        {
          owner = acc;
          best = key;
          first = false;
        }
      }
      for (Size a : members)
      {
        auto& h = proteoform_hits_[idx[a]];
        if (h.getAccession() == owner) continue;
        h.setMetaValue("Mass", -1.0);
        h.setMetaValue("MassLostToFeatureOwner", owner);
        demoted++;
      }
    }
    OPENMS_LOG_INFO << "One proteoform per MS1 feature: " << groups_with_conflict << " masses claimed by PrSMs of several proteins; " << demoted
                    << " PrSMs of the other proteins lost the mass" << std::endl;
  }

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
  //   pooled null (so the null of the sparse determined classes changes with the option; not yet measured). The PrSMs of each residual
  //   round of the chimeric search (meta value ChimericRound r > 0) are class "residual r", also by source and before the rules above;
  //   determined and undetermined ones share it (the log line of the class counts both), and its decoys join the pooled null like those
  //   of the precursor-ladder class. Each PrSM is in exactly one class, so no PrSM passes twice. With discard_underdetermined the undetermined PrSMs are in no class and take part in nothing below,
  //   the proteoforms included.
  // - proteoforms: all representatives in one ranking by their PrSM q-value, which carries the class, then by score;
  // - proteins: setProteinQvalues_.
  // Without decoys in the database and without discard_underdetermined no q-value is set.
  //   With search_cascade, the PrSM q-values are computed per class and stage (meta value SearchStage), so that the PrSMs of a later
  //   stage, from the spectra the stages before did not identify, compete only with that stage's decoys.
  std::vector<int> hit_class; // index into fdr_class_names; -1: in no class
  std::vector<double> qvalues, decoy_weights;
  std::vector<bool> decoys;
  // the proteoforms: positions of their PrSMs in proteoform_hits_ (Index), the representative first; representatives in the order above
  const auto clusters = proteoformQvalues(hit_class, qvalues, decoys, decoy_weights, true);
  const Size hit_count = proteoform_hits_.size();

  if (use_fdr)
  {
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

void FLASHTnTAlgorithm::loadRawSpectra_(const std::set<int>& scans)
{
  releaseRawSpectra_();
  std::map<int, MSSpectrum> raw;
  raw_loader_(scans, raw);
  raw_spectra_.reserve(raw.size());
  for (auto& [scan, spec] : raw)
  {
    raw_slot_[scan] = (int)raw_spectra_.size();
    raw_spectra_.push_back(std::move(spec));
  }
  raw_peaks_.resize(raw_spectra_.size());
  int profile_count = 0;
#pragma omp parallel for schedule(dynamic) reduction(+ : profile_count)
  for (int i = 0; i < (int)raw_spectra_.size(); i++)
  {
    profile_count += FLASHAnnotatorAlgorithm::toCentroidedPeakArrays(raw_spectra_[i], raw_peaks_[i]) ? 1 : 0;
    raw_spectra_[i].clear(false); // the precursor information stays
    raw_spectra_[i].shrink_to_fit();
  }
  raw_profile_count_ = profile_count;
}

void FLASHTnTAlgorithm::releaseRawSpectra_()
{
  std::vector<MSSpectrum>().swap(raw_spectra_);
  std::vector<FLASHAnnotatorAlgorithm::PeakArrays>().swap(raw_peaks_);
  raw_slot_.clear();
  raw_profile_count_ = 0;
}

void FLASHTnTAlgorithm::annotateOnRawSpectra_(std::vector<ProteinHit>& hits) const
{
  // The extension score counts deconvolved masses on the path; the raw score (FLASHAnnotatorAlgorithm) is the chance of the fragment isotope
  // envelopes found in the raw spectrum, at every charge, with the chance rate measured on mass-shifted decoy fragments of the same spectrum.
  // Decoy PrSMs are annotated the same way, so the FDR of the combined score stays target-decoy controlled. Eclip_1480: at 1%, 4846 PrSMs
  // by the extension score, 5189 by the raw score alone, 5241 by extension + 4 x raw, 5247 at 8 x raw (5240 at 16, 5012 at 1); the
  // three lysates: 4278 proteoforms / 991 proteins at 4, 4293 / 1002 at 8. The raw spectra are those of loadRawSpectra_.
  FLASHAnnotatorAlgorithm annotator;
  // only the terminal results are used below (RawScore, its fragment counts, RawCoverage): internal fragments are not annotated
  Param aparam = annotator_param_;
  aparam.setValue("internal_ions", "false");
  annotator.setParameters(aparam);
  std::vector<FLASHAnnotatorAlgorithm::PrSM> prsms(hits.size());
  std::vector<double> precursor_masses(hits.size(), -1.0);
  double max_mass = 0;
  Size without_spectrum = 0;
  for (Size i = 0; i < hits.size(); i++)
  {
    const auto& hit = hits[i];
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
    const auto it = raw_slot_.find((int)hit.getMetaValue("Scan"));
    if (it == raw_slot_.end()) without_spectrum++;
    else prsm.spectrum_index = it->second;
  }
  // the isotope table once, beyond the heaviest PrSM; then the precursor charge caps
  std::vector<double> iso_mzs, iso_intensities;
  annotator.getIsotopeEnvelope(max_mass + 200.0, 1, iso_mzs, iso_intensities);
  for (Size i = 0; i < prsms.size(); i++)
  {
    if (prsms[i].spectrum_index < 0) continue;
    int reason = 4;
    prsms[i].precursor_charge = annotator.precursorChargeCap(raw_spectra_[prsms[i].spectrum_index], precursor_masses[i], reason);
  }
  std::vector<FLASHAnnotatorAlgorithm::Result> results;
  annotator.annotate(raw_peaks_, prsms, results);

  Size annotated = 0;
  for (Size i = 0; i < hits.size(); i++)
  {
    auto& hit = hits[i];
    const auto& res = results[i];
    annotated += res.annotated ? 1 : 0;
    hit.setMetaValue("ExtensionScore", hit.getScore());
    hit.setMetaValue("RawScore", res.score);
    hit.setMetaValue("RawMatchedFragments", res.matched_fragments);
    hit.setMetaValue("RawFragmentHypotheses", res.fragment_hypotheses);
    hit.setMetaValue("RawCoverage", res.coverage);
    hit.setScore(hit.getScore() + raw_score_weight_ * std::max(0.0, res.score));
  }
  OPENMS_LOG_INFO << "Raw spectra: " << raw_peaks_.size() << " loaded (" << raw_profile_count_ << " centroided); PrSMs annotated: " << annotated
                  << " of " << hits.size() << " (" << without_spectrum << " without their scan in the raw file)" << std::endl;
}

void FLASHTnTAlgorithm::keepBestPerSpectrum_(std::vector<ProteinHit>& hits, const String& label, const String& score_name)
{
  std::map<int, Size> best;
  for (Size i = 0; i < hits.size(); i++)
  {
    const int scan = hits[i].getMetaValue("Scan");
    auto it = best.find(scan);
    if (it == best.end()) best[scan] = i;
    else if (isBetterPrSM(hits[i], hits[it->second])) it->second = i;
  }
  // A PrSM without a mass is replaced by the best PrSM of the same spectrum and protein that has a mass and the same known termini: it is
  // the same proteoform, and the precursor mass confirms it. Its higher score is no evidence against it: the terminal paths are scored
  // apart, so the b and y ions of one bond count twice there and once on the precursor-mass path, which is also penalised for an internal
  // end (Eclip_1480: RL3 180-209 at the right MS1 mass scored 116, the same range without a mass 169, and the latter was reported).
  {
    std::map<int, Size> best_with_mass;
    for (Size i = 0; i < hits.size(); i++)
    {
      const auto& h = hits[i];
      const int scan = h.getMetaValue("Scan");
      const auto& b = hits[best[scan]];
      if ((double)b.getMetaValue("Mass") > 0 || (double)h.getMetaValue("Mass") <= 0 || ! isDetermined(h)) continue;
      if (h.getAccession() != b.getAccession() || h.metaValueExists("DecoyPrecursor") != b.metaValueExists("DecoyPrecursor")) continue;
      const int bs = b.getMetaValue("StartPosition"), be = b.getMetaValue("EndPosition");
      if ((bs >= 0 && bs != (int)h.getMetaValue("StartPosition")) || (be >= 0 && be != (int)h.getMetaValue("EndPosition"))) continue;
      auto it = best_with_mass.find(scan);
      if (it == best_with_mass.end() || isBetterPrSM(h, hits[it->second])) best_with_mass[scan] = i;
    }
    for (const auto& [scan, i] : best_with_mass) best[scan] = i;
    if (! best_with_mass.empty())
      OPENMS_LOG_INFO << "PrSMs" << label << ": " << best_with_mass.size() << " spectra report the PrSM with a mass instead of the same proteoform without one"
                      << std::endl;
  }
  std::vector<ProteinHit> kept;
  kept.reserve(best.size());
  for (Size i = 0; i < hits.size(); i++)
  {
    if (best[(int)hits[i].getMetaValue("Scan")] == i) kept.push_back(std::move(hits[i]));
  }
  OPENMS_LOG_INFO << "PrSMs" << label << ": " << hits.size() << " of all candidates, " << kept.size() << " best per spectrum by the " << score_name
                  << std::endl;
  hits.swap(kept);
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