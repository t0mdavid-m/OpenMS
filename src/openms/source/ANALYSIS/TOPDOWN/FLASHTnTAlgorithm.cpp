// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/AhoCorasickAmbiguous.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHGappedTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <unordered_map>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{
// One matched tag is enough to propose a protein: a correct protein often has a single tag in a spectrum,
// and the protein-level score and the extension decide between the candidates.
inline const int min_tag_count = 1;
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
  defaults_.setValue("prsm_fdr", 1.0, "Specifies the PrSM-level FDR.");
  defaults_.setMinFloat("prsm_fdr", 0.0);

  defaults_.setValue("pro_fdr", 1.0, "Specifies the proteoform-level FDR.");
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

  defaults_.setValue("fixed_mod", std::vector<std::string> {""}, "Specifies the fixed modifications.");
  defaults_.setValue("var_mod", std::vector<std::string> {""}, "Specifies the fixed modifications.");

  defaults_.setValue("max_var_mod_count", 3, "Maximum number of variable modifications per proteoform.");
  defaults_.setValue("max_blind_mod_count", 2, "Maximum number of blind modifications per proteoform.");
  defaults_.setValue("max_mod_mass", 500.0, "Maximum absolute mass of blind modifications.");

  auto tparam = FLASHTaggerAlgorithm().getDefaults();
  tparam.remove("ion_type");
  tparam.remove("fixed_mod");
  defaults_.insert("tag:", tparam);

  auto eparam = FLASHExtenderAlgorithm().getDefaults();
  eparam.remove("ion_type");
  eparam.remove("max_var_mod_count");
  eparam.remove("max_blind_mod_count");
  eparam.remove("max_mod_mass");
  eparam.remove("fixed_mod");
  eparam.remove("var_mod");

  defaults_.insert("ex:", eparam);

  defaults_.setValue("gapped_tags", "fallback",
                     "Candidate proteins from gapped tags (masses whose differences are sums of one or more residues): off; fallback (spectra in "
                     "which no sequence tag proposes a protein); both (every spectrum, proteins not proposed by sequence tags); only (gapped "
                     "tags instead of sequence tags). Where a gapped tag lands on a protein anchors the extension as a sequence tag does.");
  defaults_.setValidStrings("gapped_tags", {"off", "fallback", "both", "only"});
  defaults_.setValue("decoy_precursor", "false",
                     "Decoy precursor masses: spectra whose PrSM took its mass from MS1 are searched again with the precursor mass moved by 5 Da "
                     "to ex:max_mod_mass; a PrSM reporting the moved mass competes with the spectrum's own and counts as a decoy in the FDR. Controls "
                     "wrong MS1 precursor masses; PrSMs with masses from fragment pairs are not affected.");
  defaults_.setValidStrings("decoy_precursor", {"true", "false"});
  defaults_.setValue("learn_var_mod_penalty", "true",
                     "Raises the penalty of a variable modification that target PrSMs use per opportunity (fitting residue or terminus) no more often than decoy PrSMs "
                     "(chance), up to that of a blind shift, and extends the spectra that used it again. Needs decoys and at least 50 decoy PrSMs.");
  defaults_.setValidStrings("learn_var_mod_penalty", {"true", "false"});
  defaults_.addTag("learn_var_mod_penalty", "advanced");
  defaults_.insert("gapped:", FLASHGappedTaggerAlgorithm().getDefaults());
  defaultsToParam_();
}

void FLASHTnTAlgorithm::updateMembers_()
{
  tagger_param_ = param_.copy("tag:", true);
  tagger_param_.setValue("ion_type", param_.getValue("ion_type"));
  tagger_param_.setValue("fixed_mod", param_.getValue("fixed_mod"));

  extender_param_ = param_.copy("ex:", true);
  extender_param_.setValue("ion_type", param_.getValue("ion_type"));
  extender_param_.setValue("max_var_mod_count", param_.getValue("max_var_mod_count"));
  extender_param_.setValue("max_blind_mod_count", param_.getValue("max_blind_mod_count"));
  extender_param_.setValue("max_mod_mass", param_.getValue("max_mod_mass"));
  extender_param_.setValue("fixed_mod", param_.getValue("fixed_mod"));
  extender_param_.setValue("var_mod", param_.getValue("var_mod"));

  gapped_param_ = param_.copy("gapped:", true);
  gapped_mode_ = param_.getValue("gapped_tags").toString();
  learn_var_mod_penalty_ = param_.getValue("learn_var_mod_penalty").toString() == "true";
  decoy_precursor_ = param_.getValue("decoy_precursor").toString() == "true";
  prsm_fdr_ = param_.getValue("prsm_fdr");
  keep_decoy_ = param_.getValue("keep_decoy").toString() == "true";
  keep_underdetermined_ = param_.getValue("discard_underdetermined").toString() == "false";
  multiple_hits_per_spec_ = param_.getValue("only_single_hit").toString() == "false";

  max_total_mod_mass_ = (double)param_.getValue("max_blind_mod_count") * (double)extender_param_.getValue("max_mod_mass") + 1.1; /// TODO

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
  if (a == b) return true;
  double mass1 = a.getMetaValue("Mass");
  double mass2 = b.getMetaValue("Mass");
  if (mass1 * mass2 < 0) return false;

  int sp1 = a.getMetaValue("StartPosition");
  int ep1 = a.getMetaValue("EndPosition");

  int sp2 = b.getMetaValue("StartPosition");
  int ep2 = b.getMetaValue("EndPosition");

  std::vector<String> mod_accs1, mod_accs2;

  if (a.metaValueExists("Modifications"))
  {
    mod_accs1 = a.getMetaValue("ModificationACCs");
    std::sort(mod_accs1.begin(), mod_accs1.end());
  }

  if (b.metaValueExists("Modifications"))
  {
    mod_accs2 = b.getMetaValue("ModificationACCs");
    std::sort(mod_accs2.begin(), mod_accs2.end());
  }

  if (mod_accs1.empty() && mod_accs2.empty()) return sp1 == sp2 && ep1 == ep2; // unmodified ones

  // at least one is modified
  bool mass_matched = std::abs(mass1 - mass2) < std::max(mass1, mass2) * tol / 1e6 * 2;
  bool mod_matched = mod_accs1.size() == mod_accs2.size();

  for (int j = 0; mod_matched && j < mod_accs1.size(); j++)
  {
    if (!mod_accs1[j].empty() && !mod_accs2[j].empty())
    {
      std::set<std::string> set1, set2;
      std::stringstream ss1(mod_accs1[j]), ss2(mod_accs2[j]);
      std::string item;

      while (std::getline(ss1, item, ',')) {
        set1.insert(item);
      }
      while (std::getline(ss2, item, ',')) {
        set2.insert(item);
      }
      std::set<std::string> result;
      std::set_intersection(
        set1.begin(), set1.end(),
        set2.begin(), set2.end(),
        std::inserter(result, result.begin())
      );
      if (result.empty())
      {
        mod_matched = false;
      }
    }
    else
    {
      auto mod_masses1 = a.getMetaValue("Modifications").toDoubleList();
      auto mod_masses2 = b.getMetaValue("Modifications").toDoubleList();
      for (int k = 0; k < mod_masses1.size(); k++)
      {
        if (std::abs(mod_masses1[k] - mod_masses2[k]) > std::max(mass1, mass2) * tol / 1e6 * 2) // mass difference margin by precursor masses
        {
          mod_matched = false;
          break;
        }
      }
    }
  }

  bool mod_loc_matched = mod_matched;

  if (mod_matched && mod_accs1.size() == mod_accs2.size())
  {
    std::vector<int> mod_starts1 = a.getMetaValue("ModificationStarts");
    std::vector<int> mod_starts2 = b.getMetaValue("ModificationStarts");
    std::vector<int> mod_ends1 = a.getMetaValue("ModificationEnds");
    std::vector<int> mod_ends2 = b.getMetaValue("ModificationEnds");

    for (Size i = 0; i < mod_starts1.size(); i++)
    {
      if (mod_starts1[i] > mod_ends2[i] || mod_starts2[i] > mod_ends1[i]) { mod_loc_matched = false; break; }
    }
  }

  if (mass1 > 0 && mass_matched)
  {
    //false if mods are different or mods are the same but mod locs are different
    if (!mod_matched || !mod_loc_matched) return false;
    return true;
  }

  // masses are underdetermined
  if (sp1 == sp2 && ep1 == ep2 && mod_matched && mod_loc_matched)  return true;
  // true if sp ep are the same and mod accs are the same and locs are consistent.
  return false;
}

void FLASHTnTAlgorithm::markRepresentativeProteoformHits_(double tol)
{
  std::map<String, std::vector<int>> proteoform_map;

  for (int i = 0; i < proteoform_hits_.size(); i++)
  {
    auto& hit = proteoform_hits_[i];
    hit.setMetaValue("Index", i);
    const auto& acc = hit.getAccession();
    proteoform_map[acc].push_back(i);
  }

  for (auto& [acc, hit_idx] : proteoform_map)
  {
    for (int i = 0; i < hit_idx.size(); i++)
    {
      std::vector<int> prsm_indices;
      auto& hit = proteoform_hits_[hit_idx[i]];
      if (hit.metaValueExists("Representative")) continue;
      prsm_indices.push_back(hit.getMetaValue("Index"));
      for (int j = i + 1; j < hit_idx.size(); j++)
      {
        auto& hit2 = proteoform_hits_[hit_idx[j]];
        if (areConsistent_(hit, hit2, tol))
        {
          prsm_indices.push_back(hit2.getMetaValue("Index"));
          hit2.setMetaValue("Representative", 0);
        }
      }
      hit.setMetaValue("Representative", 1);
      hit.setMetaValue("PrSMIndices", prsm_indices);
    }
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
  int max_blind_mod_count = param_.getValue("max_blind_mod_count");

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

  double precursor_tol = -1, tol, ms2_tol = -1;

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
    if (spec.getMSLevel() == 1 && precursor_tol < 0)
    {
      precursor_tol = tol;
      // continue;
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

    // a precursor mass may be given (e.g. by precursor_charge/precursor_mz options) even if no precursor scan is present
    if (precursor_scan > 0 || precursor_mass > 0)
    {
      PeakGroup pg;
      pg.setMonoisotopicMass(precursor_mass);
      if (deconv_meta_str.hasSubstring("precursorscore="))
      {
        int s_loc_preq_s = deconv_meta_str.find("precursorscore=") + 15;
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
      dspec.sortByQscore2D();
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
  // Meta values set by name look the name up in the global registry under a critical section; register the names once and set
  // them by index inside the parallel loop.
  auto& registry = MetaInfoInterface::metaRegistry();
  const UInt mi_scan = registry.registerName("Scan"), mi_matched_aa = registry.registerName("MatchedAA"),
             mi_tag_indices = registry.registerName("TagIndices"), mi_tag_positions = registry.registerName("TagPositions"),
             mi_fasta_index = registry.registerName("FastaIndex"), mi_protein_rank = registry.registerName("ProteinRank"),
             mi_gapped = registry.registerName("GappedTag");
  // gapped tags (only for spectra in which no sequence tag proposes a protein); their tags are kept per spectrum until appended to tags_
  FLASHGappedTaggerAlgorithm gapped_tagger;
  if (gapped_mode_ != "off")
  {
    gapped_tagger.setParameters(gapped_param_);
    gapped_tagger.buildIndex(protein_seqs);
    OPENMS_LOG_INFO << "Gapped tag index: " << gapped_tagger.getIndexBytes() / (1024 * 1024) << " MB" << std::endl;
  }
  std::vector<std::vector<FLASHHelperClasses::Tag>> spec_gapped_tags(dspecs.size());
  // Parallel over spectra only when there are at least as many spectra with tags as threads; otherwise (e.g. a targeted native run with
  // one or a few spectra) the threads go to the scoring of each spectrum's candidates (runMatchingExact).
  int spectra_with_tags = 0;
  for (const auto& dspec : dspecs) spectra_with_tags += gapped_mode_ != "off" || scan_to_tag_indices.count(dspec.getScanNumber()) > 0;
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
            if (tag.getNtermMass() > nmass + max_total_mod_mass_ + 200) continue;
            pi_to_n_shifts[pis[i]].push_back(nmass - tag.getNtermMass());
            pi_to_n_flanks[pis[i]].push_back(tag.getNtermMass());
          }
          if (tag.getCtermMass() >= 0 && pps[i] + seq.length() >= 0 && pseq.length() > pps[i] + seq.length())
          {
            double cmass = tag_to_c_flanking_masses.at(seq)[i];
            if (tag.getCtermMass() > cmass + max_total_mod_mass_ + 200) continue;
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
          if (flank > la.flank_mass + max_total_mod_mass_ + 200) continue;
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

    if (hits.empty()) continue;

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

    for (Size i = 0; i < hits.size(); i++)
    {
      hits[i].setMetaValue(mi_protein_rank, (int)i + 1);
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
  }
  std::vector<std::vector<ProteinHit>>().swap(spec_hits);

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
    for (const auto& r : candidate_ranges) spectra_with_candidates += r.first != r.second;
#ifdef _OPENMP
    const bool extend_over_spectra = spectra_with_candidates >= omp_get_max_threads();
#else
    const bool extend_over_spectra = false;
#endif
    // the extension of one spectrum's candidates with the given extension parameters
    auto extendSpectrumInto = [&](int d, const Param& eparam, const DeconvolvedSpectrum& dspec, std::vector<ProteinHit>& out) {
      out.clear();
      const auto [first, second] = candidate_ranges[d];
      if (first == second) return;

      std::vector<ProteinHit> hits;
      for (Size i = first; i < second; i++)
      {
        bool top = (int)(i - first) < min_protein_hits_;
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

      FLASHExtenderAlgorithm extender;
      extender.setParameters(eparam);

      extender.run(hits, dspec, spec_vec, vec_pro, rev_vec_pro, tags_, tol, multiple_hits_per_spec_);

      extender.fillProteoforms(out);
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
    // Decoy precursors: a spectrum whose PrSM took its mass from MS1 is extended again with the precursor mass moved by 5 to max_mod_mass
    // Da (a fixed pseudo-random shift per scan). A PrSM of that copy that reports the moved MS1 mass (not a mass from fragment pairs) is a
    // decoy of the precursor; it competes with the spectrum's own PrSMs (the best one stays) and counts as a decoy in the FDR. If the MS1
    // mass is right, the copy needs a blind shift or a truncation and loses; if it is wrong, both are alike and the decoy wins half the time.
    if (decoy_precursor_)
    {
      const double max_shift = std::max(10.0, (double)extender_param_.getValue("max_mod_mass"));
      const double iso = Constants::C13C12_MASSDIFF_U;
      auto ms1Based = [](const ProteinHit& h) {
        return (double)h.getMetaValue("Mass") > 0 && (int)h.getMetaValue("ProteoformMassByFragmentMass") == 0;
      };
      const Param& eparam_now = learned_extender_param.empty() ? extender_param_ : learned_extender_param;
      std::vector<std::vector<ProteinHit>> decoy_proteoforms(dspecs.size());
      int copies = 0;
      std::vector<int> todo;
      for (int d = 0; d < (int)dspecs.size(); d++)
      {
        if (dspecs[d].getPrecursorPeakGroup().getMonoMass() <= 0) continue;
        if (std::none_of(spec_proteoforms[d].begin(), spec_proteoforms[d].end(), ms1Based)) continue;
        todo.push_back(d);
      }
      copies = (int)todo.size();
#pragma omp parallel for schedule(dynamic) if (extend_over_spectra)
      for (int k = 0; k < (int)todo.size(); k++)
      {
        const int d = todo[k];
        DeconvolvedSpectrum decoy_spec = dspecs[d];
        auto pg = decoy_spec.getPrecursorPeakGroup();
        // pseudo-random shift in [5, max_shift] Da, either sign, fixed by the scan number
        std::uint64_t x = (std::uint64_t)decoy_spec.getScanNumber() * 0x9E3779B97F4A7C15ULL;
        x ^= x >> 29;
        x *= 0xBF58476D1CE4E5B9ULL;
        x ^= x >> 32;
        const double u = (double)(x % 1000003ULL) / 1000003.0;
        const double shift = (5.0 + u * (max_shift - 5.0)) * ((x >> 40) & 1 ? 1.0 : -1.0);
        const double decoy_mass = pg.getMonoMass() + shift;
        pg.setMonoisotopicMass(decoy_mass);
        decoy_spec.setPrecursorPeakGroup(pg);
        std::vector<ProteinHit> out;
        extendSpectrumInto(d, eparam_now, decoy_spec, out);
        for (auto& h : out)
        {
          if (! ms1Based(h)) continue;
          if (std::abs((double)h.getMetaValue("Mass") - decoy_mass) > 3 * iso + 0.1) continue; // the moved mass, up to isotope errors
          h.setMetaValue("DecoyPrecursor", 1);
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
      OPENMS_LOG_INFO << "Decoy precursors: " << copies << " spectra with an MS1-based PrSM, decoy precursor wins in " << decoy_wins << std::endl;
    }

    for (auto& v : spec_proteoforms)
    {
      for (auto& h : v) proteoform_hits_.push_back(std::move(h));
    }
    endProgress();
  }

  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), [](const ProteinHit& left, const ProteinHit& right) {
    return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getMetaValue("Scan") > right.getMetaValue("Scan"))
                                                                                            : (left.getCoverage() > right.getCoverage()))
                                               : (left.getScore() > right.getScore());
  });

  markRepresentativeProteoformHits_(precursor_tol);

  if (decoy_factor_ > 0 || ! keep_underdetermined_)
  {
    std::vector<ProteinHit> filtered_proteoform_hits;
    filtered_proteoform_hits.reserve(proteoform_hits_.size());

    for (int k = 0; k < (keep_underdetermined_ ? 2 : 1); k++)
    {
      for (Size mod = 0; mod <= max_blind_mod_count; mod++)
      {
        double taget_count = 0;
        double decoy_count = 0;

        double taget_count_pro = 0;
        double decoy_count_pro = 0;

        std::map<double, double> map_qvalue;
        std::map<double, double> map_qvalue_pro;

        for (auto& hit : proteoform_hits_)
        {
          // FDR per number of blind modifications; variable modifications do not count
          if (blindModCount_(hit) != mod) continue;
          if (k == 0
              && ((double)hit.getMetaValue("Mass") < 0 || (int)hit.getMetaValue("StartPosition") < 0 || (int)hit.getMetaValue("EndPosition") < 0))
          {
            continue;
          }
          else if (k == 1
                   && ((double)hit.getMetaValue("Mass") > 0 && (int)hit.getMetaValue("StartPosition") > 0
                       && (int)hit.getMetaValue("EndPosition") > 0))
            continue;

          const bool decoy_precursor = hit.metaValueExists("DecoyPrecursor");
          bool is_decoy = hit.getAccession().hasPrefix("DECOY") || decoy_precursor;
          bool is_rep = (int)hit.getMetaValue("Representative") == 1;
          if (is_decoy)
          {
            // one decoy precursor per spectrum: it stands for one target, a protein decoy for 1 / decoy_factor_
            const double w = decoy_precursor && ! hit.getAccession().hasPrefix("DECOY") ? 1.0 : 1.0 / decoy_factor_;
            decoy_count += w;
            if (is_rep) decoy_count_pro += w;
          }
          else
          {
            taget_count++;
            if (is_rep) taget_count_pro++;
          }

          double tmp_qvalue = decoy_count / (decoy_count + taget_count);
          map_qvalue[hit.getScore()] = std::min(1.0, tmp_qvalue);

          if (! is_rep) continue;
          double tmp_qvalue_pro = decoy_count_pro / (decoy_count_pro + taget_count_pro);
          map_qvalue_pro[hit.getScore()] = std::min(1.0, tmp_qvalue_pro);
        }

        double cummin = 1.0;
        for (auto&& rit = map_qvalue.begin(); rit != map_qvalue.end(); ++rit)
        {
          cummin = std::min(rit->second, cummin);
          rit->second = cummin;
        }

        cummin = 1.0;
        for (auto&& rit = map_qvalue_pro.begin(); rit != map_qvalue_pro.end(); ++rit)
        {
          cummin = std::min(rit->second, cummin);
          rit->second = cummin;
        }

        for (auto& hit : proteoform_hits_)
        {
          // FDR per number of blind modifications; variable modifications do not count
          if (blindModCount_(hit) != mod) continue;
          if (k == 0
              && ((double)hit.getMetaValue("Mass") < 0 || (int)hit.getMetaValue("StartPosition") < 0 || (int)hit.getMetaValue("EndPosition") < 0))
            continue;
          if (k == 1
              && ((double)hit.getMetaValue("Mass") > 0 && (int)hit.getMetaValue("StartPosition") > 0 && (int)hit.getMetaValue("EndPosition") > 0))
            continue;

          bool is_decoy = hit.getAccession().hasPrefix("DECOY") || hit.metaValueExists("DecoyPrecursor");
          double qvalue = map_qvalue[hit.getScore()];
          hit.setMetaValue("qvalue", qvalue);
          if (! keep_decoy_ && is_decoy) continue;

          if (prsm_fdr_ < 1 && qvalue > prsm_fdr_) continue;

          auto iter = map_qvalue_pro.lower_bound(hit.getScore());
          if (iter != map_qvalue_pro.end())
          {
            double qvalue_pro = iter->second;
            hit.setMetaValue("proqvalue", qvalue_pro);
          }
          else
            hit.setMetaValue("proqvalue", 1.0);
          // if (pro_fdr_ < 1 && qvalue > pro_fdr_) continue;

          filtered_proteoform_hits.push_back(hit);
        }
      }
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
      new_prsm_indices.push_back(prsm_index_map[i]);
    }
    std::sort(new_prsm_indices.begin(), new_prsm_indices.end());
    hit.setMetaValue("PrSMIndices", new_prsm_indices);
  }
  // TODO
  // per scan, get the scan number and the precursor mass - maybe given or maybe deconvolved.
  // get all protoeforms - mass and score.
  // get all combinations with a certain tolerance.
  // record the best ranking combination(s) along with summed score (can score be non overlapping?) and mass difference
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