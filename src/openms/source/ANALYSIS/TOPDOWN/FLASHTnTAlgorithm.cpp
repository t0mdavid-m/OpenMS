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
#include <cmath>
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
// kept for parsers of the label.
constexpr int fdr_class_max_blind = 2;
constexpr int fdr_class_both_termini = fdr_class_max_blind + 1;
constexpr int fdr_class_open_terminus = fdr_class_both_termini + 1;
constexpr int fdr_class_count = fdr_class_open_terminus + 1;
const char* const fdr_class_names[] = {"determined; 0 blind", "determined; 1 blind", "determined; 2 blind", "undetermined, both termini",
                                       "undetermined, open terminus"};
static_assert(std::size(fdr_class_names) == (Size)fdr_class_count, "one label per FDR class");
// protein level: the null of a class is exponential above this quantile of its decoy scores
constexpr double null_tail_quantile = 0.8;
// protein level: with fewer decoys the null is the empirical survival function only, without the tail
constexpr Size null_min_tail_decoys = 5;
// protein level: a class with fewer decoys uses the null of the decoys of all classes pooled
constexpr Size null_min_class_decoys = 20;
// protein level: floor of a null p-value and of a local FDR before its log
constexpr double fdr_p_floor = 1e-300;

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
  defaults_.setValue("decoy_precursor", "false",
                     "Decoy precursors: spectra whose PrSM took its mass from MS1 are searched again with the decoy precursor mass of their isolation "
                     "window instead of the precursor mass (FLASHDeconv -FD:report_FDR true writes it); a PrSM reporting that mass competes with the "
                     "spectrum's own and counts as a decoy in the FDR. Controls wrong MS1 precursor masses; PrSMs with masses from fragment pairs "
                     "or without a mass are not affected.");
  defaults_.setValidStrings("decoy_precursor", {"true", "false"});
  defaults_.setValue("learn_var_mod_penalty", "true",
                     "Raises the penalty of a variable modification that target PrSMs use per opportunity (fitting residue or terminus) no more often than decoy PrSMs "
                     "(chance), up to that of a blind shift, and extends the spectra that used it again. Needs decoys and at least 50 decoy PrSMs.");
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
  // - both masses are determined: the masses agree (within twice the tolerance) and the ranges overlap, as TopPIC clusters PrSMs.
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
    return overlap && std::abs(mass1 - mass2) <= std::max(mass1, mass2) * tol * 2e-6;
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
             mi_gapped = registry.registerName("GappedTag"), mi_terminal = registry.registerName("TerminalLadder");
  // gapped tags (only for spectra in which no sequence tag proposes a protein); their tags are kept per spectrum until appended to tags_
  FLASHGappedTaggerAlgorithm gapped_tagger;
  if (gapped_mode_ != "off")
  {
    gapped_tagger.setParameters(gapped_param_);
    gapped_tagger.buildIndex(protein_seqs);
    OPENMS_LOG_INFO << "Gapped tag index: " << gapped_tagger.getIndexBytes() / (1024 * 1024) << " MB" << std::endl;
  }
  std::vector<std::vector<FLASHHelperClasses::Tag>> spec_gapped_tags(dspecs.size());

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
    const std::string standard = "ACDEFGHIKLMNPQRSTVWYUO";
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
      if (prefix)
      {
        for (double form : {0.0, met, -acetyl, met - acetyl}) n_ladder_shifts.push_back(form - offset);
      }
      else { c_ladder_shifts.push_back(-offset); }
    }
  }
  // Parallel over spectra only when there are at least as many spectra with tags as threads; otherwise (e.g. a targeted native run with
  // one or a few spectra) the threads go to the scoring of each spectrum's candidates (runMatchingExact).
  int spectra_with_tags = 0;
  for (const auto& dspec : dspecs) spectra_with_tags += gapped_mode_ != "off" || terminal_ladders_ || scan_to_tag_indices.count(dspec.getScanNumber()) > 0;
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
    ladder_ranges[d].first = ladder_candidates.size();
    for (auto& hit : spec_ladder_hits[d])
    {
      auto indices = hit.getMetaValue(mi_tag_indices).toIntList(); // -(local index + 1)
      for (auto& i : indices) i = tag_base - i - 1;
      hit.setMetaValue(mi_tag_indices, indices);
      ladder_candidates.push_back(std::move(hit));
    }
    ladder_ranges[d].second = ladder_candidates.size();
  }
  std::vector<std::vector<ProteinHit>>().swap(spec_hits);
  std::vector<std::vector<ProteinHit>>().swap(spec_ladder_hits);

  endProgress();

  if (! protein_level_only_)
  {
    // Beyond the min_protein_hits best candidates of a spectrum (a fair target-decoy competition: they are chosen by this spectrum
    // alone), extend only those whose protein, or its target-decoy partner, passes the run-level protein FDR. That needs decoys;
    // without them every candidate is extended, as before.
    const bool use_protein_fdr = decoy_factor_ > 0 && protein_fdr_ < 1.0;
    std::set<Size> passing_proteins;
    if (use_protein_fdr) { passing_proteins = getProteinsPassingFDR_(fasta_entry); }

    startProgress(0, (SignedSize)dspecs.size(), "Running extension algorithm ...");
    // spectra in parallel, each with its own extender; proteoforms are appended in spectrum order afterwards
    std::vector<std::vector<ProteinHit>> spec_proteoforms(dspecs.size());
    // as above: over spectra only with enough spectra to extend, otherwise over each spectrum's candidates inside the extender
    int spectra_with_candidates = 0;
    for (Size d = 0; d < dspecs.size(); d++) spectra_with_candidates += candidate_ranges[d].first != candidate_ranges[d].second || ladder_ranges[d].first != ladder_ranges[d].second;
#ifdef _OPENMP
    const bool extend_over_spectra = spectra_with_candidates >= omp_get_max_threads();
#else
    const bool extend_over_spectra = false;
#endif
    const bool keep_all_candidates = (bool)raw_loader_;
    // the extension of one spectrum's candidates with the given extension parameters
    auto extendSpectrumInto = [&](int d, const Param& eparam, const DeconvolvedSpectrum& dspec, std::vector<ProteinHit>& out) {
      out.clear();
      const auto [first, second] = candidate_ranges[d];
      if (first == second && ladder_ranges[d].first == ladder_ranges[d].second) return;

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
      if (ladder_hits.empty()) return;
      std::vector<ProteinHit> ladder_out;
      extend(ladder_hits, ladder_out);
      std::vector<ProteinHit> determined, undetermined;
      for (auto& h : ladder_out)
      {
        h.setMetaValue(mi_terminal, 1);
        const bool det = (double)h.getMetaValue("Mass") > 0 && (int)h.getMetaValue("StartPosition") > 0 && (int)h.getMetaValue("EndPosition") > 0;
        (det ? determined : undetermined).push_back(h);
      }
      const bool fill = out.empty();
      auto is_determined = [](const ProteinHit& h) {
        return (double)h.getMetaValue("Mass") > 0 && (int)h.getMetaValue("StartPosition") > 0 && (int)h.getMetaValue("EndPosition") > 0;
      };
      if (multiple_hits_per_spec_ || keep_all_candidates) { out.insert(out.end(), determined.begin(), determined.end()); }
      else if (! determined.empty() && (out.empty() || ! is_determined(out[0]) || determined[0].getScore() > out[0].getScore()))
      {
        // a ladder PrSM that explains the precursor mass also replaces an undetermined one of a higher score
        out.swap(determined);
      }
      if (fill && out.empty()) out.swap(undetermined);
      else if (fill && (multiple_hits_per_spec_ || keep_all_candidates)) { out.insert(out.end(), undetermined.begin(), undetermined.end()); }
    };
    auto extendSpectrum = [&](int d, const Param& eparam) { extendSpectrumInto(d, eparam, dspecs[d], spec_proteoforms[d]); };
#pragma omp parallel for schedule(dynamic) if (extend_over_spectra)
    for (int d = 0; d < (int)dspecs.size(); d++)
    {
      nextProgress();
      extendSpectrum(d, extender_param_);
    }

    // Variable modification penalties from this run: a decoy (shuffled) protein can use a modification only by chance, so a modification
    // that target PrSMs use per opportunity (a fitting residue in the proteoform, or the proteoform terminus) no more often than decoy PrSMs
    // explains nothing and gets a penalty up to that of a blind shift. Penalties only go up, so only the spectra in which some PrSM used a
    // raised modification are extended again: the paths of the others keep their scores and the result equals a second full extension.
    Param learned_extender_param; // empty unless penalties were raised
    FLASHExtenderAlgorithm probe;
    probe.setParameters(extender_param_);
    const auto var_names = probe.getVariableModificationNames();
    const auto var_residues = probe.getVariableModificationResidues();
    const auto var_positions = probe.getVariableModificationPositions();
    if (learn_var_mod_penalty_ && ! var_names.empty() && decoy_factor_ > 0)
    {
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
      // uses and opportunities per modification, targets and decoys
      std::vector<double> use_t(n_mods, 0), use_d(n_mods, 0), opp_t(n_mods, 0), opp_d(n_mods, 0);
      double n_d = 0;
      std::vector<char> used;
      for (const auto& v : spec_proteoforms)
      {
        for (const auto& h : v)
        {
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
      // enough decoy PrSMs to estimate chance usage
      if (n_d >= 50)
      {
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
          Param eparam = extender_param_;
          eparam.setValue("var_mod_penalty", penalties);
          learned_extender_param = eparam;
          std::vector<int> again;
          for (int d = 0; d < (int)dspecs.size(); d++)
          {
            bool hit_raised = false;
            for (const auto& h : spec_proteoforms[d])
            {
              usedMods(h, used);
              for (Size m = 0; m < n_mods && ! hit_raised; m++) hit_raised = used[m] && penalties[m] > 2;
              if (hit_raised) break;
            }
            if (hit_raised) again.push_back(d);
          }
          OPENMS_LOG_INFO << "Extending " << again.size() << " spectra again with the learned penalties" << std::endl;
#pragma omp parallel for schedule(dynamic) if (extend_over_spectra)
          for (int k = 0; k < (int)again.size(); k++)
          {
            extendSpectrum(again[k], eparam);
          }
        }
      }
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
        extendSpectrumInto(d, eparam_now, decoy_spec, out);
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
  //   termini, 1030 for 2228 with an open terminus). Each PrSM is in exactly one class, so no PrSM passes twice. With
  //   discard_underdetermined the undetermined PrSMs are in no class and take part in nothing below, the proteoforms included.
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
      if (isDetermined(hit)) c = (int)std::min(blindModCount_(hit), (Size)fdr_class_max_blind);
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
    // else the pooled one)
    for (int c = 0; c < fdr_class_count; c++)
    {
      Size targets = 0, decoy_count = 0;
      double weight = 0;
      for (Size i : class_hits[c])
      {
        if (decoys[i])
        {
          decoy_count++;
          weight += decoy_weights[i];
        }
        else { targets++; }
      }
      String line = String("FDR class \"") + fdr_class_names[c] + "\": " + String(targets) + " target and " + String(decoy_count)
                    + " decoy PrSMs (decoy weight " + String(weight) + "); protein-level null: "
                    + (decoy_count >= null_min_class_decoys ? "own" : "pooled");
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