// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/ConvolutionBasedProteinFilter.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

namespace OpenMS
{
inline const int max_hit_count = 10;
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
  defaults_.setValue("prsm_fdr", 1.0, "PrSM level FDR");
  defaults_.setMinFloat("prsm_fdr", 0.0);

  defaults_.setValue("pro_fdr", 1.0, "Proteoform level FDR");
  defaults_.setMinFloat("pro_fdr", 0.0);

  defaults_.setValue("only_single_hit", "false", "Allow only a single hit per spectrum");
  defaults_.setValidStrings("only_single_hit", {"true", "false"});

  defaults_.setValue("discard_underdetermined", "false",
                     "To discard underdetermined proteoform hits (i.e., proteoforms without total mass or start end positions");
  defaults_.setValidStrings("discard_underdetermined", {"true", "false"});

  defaults_.setValue("keep_decoy", "false", "To keep decoy hits");
  defaults_.setValidStrings("keep_decoy", {"true", "false"});

  defaults_.setValue("ion_type", std::vector<std::string> {"b", "y"}, "Ion types to consider. Write from the most to the least dominant ion types");
  defaults_.setValidStrings("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

  auto tparam = FLASHTaggerAlgorithm().getDefaults();
  tparam.remove("ion_type");
  defaults_.insert("tag:", tparam);
  auto eparam = FLASHExtenderAlgorithm().getDefaults();
  eparam.remove("ion_type");

  defaults_.insert("ex:", eparam);
  defaultsToParam_();
}

void FLASHTnTAlgorithm::updateMembers_()
{
  tagger_param_ = param_.copy("tag:", true);
  tagger_param_.setValue("ion_type", param_.getValue("ion_type"));
  extender_param_ = param_.copy("ex:", true);
  extender_param_.setValue("ion_type", param_.getValue("ion_type"));
  prsm_fdr_ = param_.getValue("prsm_fdr");
  pro_fdr_ = param_.getValue("pro_fdr");
  keep_decoy_ = param_.getValue("keep_decoy").toString() == "true";
  keep_underdetermined_ = param_.getValue("discard_underdetermined").toString() == "false";
  multiple_hits_per_spec_ = param_.getValue("only_single_hit").toString() == "false";
}

bool FLASHTnTAlgorithm::areConsistent_(const ProteinHit& a, const ProteinHit& b, double tol) const
{
  double mass1 = a.getMetaValue("Mass");
  double mass2 = b.getMetaValue("Mass");
  if (mass1 * mass2 < 0) return false;
  if (mass1 > 0 && mass2 > 0 && std::abs(mass1 - mass2) > std::max(mass1, mass2) * tol / 1e6 * 2) return false;

  //double rt1 = a.getMetaValue("RT");
  //double rt2 = b.getMetaValue("RT");
  //if (std::abs(rt1 - rt2) < 30.0) return true; // if rts are within 30 sec, true

  if (a.metaValueExists("mod_masses") && b.metaValueExists("mod_masses"))
  {
    std::vector<double> mod_masses1 = a.getMetaValue("mod_masses");
    std::vector<double> mod_masses2 = b.getMetaValue("mod_masses");
    if (mod_masses1.size() != mod_masses2.size()) return false;
    else
    {
      std::vector<int> mod_starts1 = a.getMetaValue("ModificationStarts");
      std::vector<int> mod_starts2 = b.getMetaValue("ModificationStarts");
      std::vector<int> mod_ends1 = a.getMetaValue("ModificationEnds");
      std::vector<int> mod_ends2 = b.getMetaValue("ModificationEnds");

      for (Size i = 0; i < mod_masses1.size(); i++)
      {
        if (std::abs(mod_masses1[i] - mod_masses2[i]) > std::max(mod_masses1[i], mod_masses2[i]) * tol / 1e6 * 2)
        {
          return false;
        }
        else if (mod_starts1[i] > mod_ends2[i] || mod_starts2[i] > mod_ends1[i])
        {
          return false;
        }
      }
    }
    return true;
  }
  else if (! a.metaValueExists("mod_masses") && ! a.metaValueExists("mod_masses"))
    return true;
  return false;
}


void FLASHTnTAlgorithm::markRepresentativeProteoformHits_(double tol)
{
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(),
            [](const ProteinHit& left, const ProteinHit& right) { return left.getScore() > right.getScore(); });
  std::map<String, std::vector<ProteinHit>> proteoform_map;
  for (auto& hit : proteoform_hits_)
  {
    const auto& acc = hit.getAccession();

    if (proteoform_map.find(acc) != proteoform_map.end())
    {
      bool skip = false;
      for (const auto& hit2 : proteoform_map[acc])
      {
        if (areConsistent_(hit, hit2, tol))
        {
          skip = true;
          break;
        }
      }
      if (skip) continue;
    }
    hit.setMetaValue("Representative", true);

    proteoform_map[acc].push_back(hit);
    // std::cout<<acc << " " <<  (proteoform_map.find("sp|P02359|RS7_ECOLI") != proteoform_map.end()) << " " << tmp.size() << std::endl;
  }
}


void FLASHTnTAlgorithm::run(const MSExperiment& map, const std::vector<FASTAFile::FASTAEntry>& fasta_entry, double flanking_mass_tol)
{
  setLogType(CMD);
  startProgress(0, (SignedSize)map.size(), "Running FLASHTnT ...");
  //int max_tag_length = tagger_param_.getValue("max_length");
  //int min_tag_length = tagger_param_.getValue("min_length");
  int max_mod_cntr = extender_param_.getValue("max_mod_count");
  double max_mod_mass = max_mod_cntr * (double)extender_param_.getValue("max_mod_mass") + 1.0;
  std::map<double, std::vector<ResidueModification>> mod_map;
  const auto inst = ModificationsDB::getInstance(); // give this from outside ...
  std::vector<String> mod_strs;
  inst->getAllSearchModifications(mod_strs);
  for (int i = 0; i < mod_strs.size(); i++)
  {
    const auto mod = *inst->getModification(mod_strs[i]);
    if (std::abs(mod.getDiffMonoMass()) > max_mod_mass) continue;
    mod_map[mod.getDiffMonoMass()].push_back(mod);
  }
  double precursor_tol = -1;
  std::vector<boost::dynamic_bitset<>> vectorized_fasta_entry, rev_vectorized_fasta_entry;
  std::vector<std::vector<int>> vectorized_fasta_entry_indices, rev_vectorized_fasta_entry_indices;
  std::vector<std::map<int, double>> mass_map, rev_mass_map;
  std::vector<std::vector<Size>> bit_protein_indices, rev_bit_protein_indices;
  ConvolutionBasedProteinFilter::vectorizeFasta(fasta_entry, vectorized_fasta_entry, vectorized_fasta_entry_indices, mass_map, bit_protein_indices, false);
  ConvolutionBasedProteinFilter::vectorizeFasta(fasta_entry, rev_vectorized_fasta_entry, rev_vectorized_fasta_entry_indices, rev_mass_map, rev_bit_protein_indices, true);

  std::vector<std::map<int,std::set<Size>>> fasta_index, rev_fasta_index;

  for (int index = 0; index < map.size(); index++)
  {
    auto spec = map[index];
    nextProgress();
    int scan = FLASHDeconvAlgorithm::getScanNumber(map, index);

    //if (scan > 1300) continue; // TODO
    if (spec.getMSLevel() == 1 && precursor_tol > 0) { continue; }

    DeconvolvedSpectrum dspec(scan);
    dspec.setOriginalSpectrum(spec);
    String deconv_meta_str = spec.getMetaValue("DeconvMassInfo").toString();
    int tol_loc_s = deconv_meta_str.find("tol=") + 4;
    int tol_loc_e = deconv_meta_str.find(";", tol_loc_s);

    double tol = stod(deconv_meta_str.substr(tol_loc_s, tol_loc_e - tol_loc_s));
    if (spec.getMSLevel() == 1 && precursor_tol < 0)
    {
      precursor_tol = tol;
      continue;
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

    int s_loc_pre_s = deconv_meta_str.find("precursorscan=") + 14;
    int s_loc_pre_e = deconv_meta_str.find(";", s_loc_pre_s);
    int precursor_scan = stoi(deconv_meta_str.substr(s_loc_pre_s, s_loc_pre_e - s_loc_pre_s));

    if (precursor_scan > 0)
    {
      int s_loc_prem_s = deconv_meta_str.find("precursormass=") + 14;
      int s_loc_prem_e = deconv_meta_str.find(";", s_loc_prem_s);
      double precursor_mass = stod(deconv_meta_str.substr(s_loc_prem_s, s_loc_prem_e - s_loc_prem_s));
      PeakGroup pg;
      pg.setMonoisotopicMass(precursor_mass);
      dspec.setPrecursorPeakGroup(pg);
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

    for (int i = 0; i < spec.size(); i++)
    {
      PeakGroup peak;
      peak.setQscore(qscores[i]);
      peak.setSNR(snrs[i]);
      peak.setMonoisotopicMass(spec[i].getMZ());
      peak.setScanNumber(scan);
      dspec.push_back(peak);
    }
    dspec.sort();

    if (dspec.size() < 5) continue;

    FLASHTaggerAlgorithm tagger;
    // Run tagger
    tagger.setParameters(tagger_param_);
    tagger.run(dspec, tol);

    FLASHExtenderAlgorithm extender;
    extender.setParameters(extender_param_);
    extender.setModificationMap(mod_map);
    tagger.getTags(tags_);
    std::vector<ProteinHit> hits;
    std::vector<FLASHHelperClasses::Tag> tags;
    bool hit_by_tag = false;
    bool proteoform_found = false;

    tagger.runMatching(fasta_entry, dspec,  vectorized_fasta_entry, rev_vectorized_fasta_entry, mass_map, rev_mass_map, max_mod_mass);
    tagger.getTags(tags);
    tagger.getProteinHits(hits, max_hit_count);
    hit_by_tag |= !hits.empty();
    extender.run(hits, tags, dspec,
                 tagger.getSpectrum(), flanking_mass_tol, tol, multiple_hits_per_spec_);
    extender.getProteoforms(proteoform_hits_);
    proteoform_found = extender.hasProteoforms();

//    if (false && !hit_by_tag && !proteoform_found) // TODO
//    {
//      ConvolutionBasedProteinFilter filter;
//      hits.clear();
//      tags.clear();
//
//      filter.runMatching(dspec, fasta_entry, vectorized_fasta_entry_indices, rev_vectorized_fasta_entry_indices,
//                         bit_protein_indices, rev_bit_protein_indices,min_tag_length);
//      filter.getProteinHits(hits, max_hit_count);
//      extender.run(hits, tags, dspec, tagger.getSpectrum(), flanking_mass_tol, tol, multiple_hits_per_spec_);
//      extender.getProteoforms(proteoform_hits_);
//    }

    decoy_factor_ = tagger.getDecoyFactor();
    // if (! out_tag_file.empty()) { FLASHTnTFile::writeTags(tagger, extender, out_tagger_stream); }
  }
  endProgress();
  // redefine tag index and proteoform to tag indics
  std::map<int, int> scan_index_offset;
  int offset = 0;
  for (auto& tag : tags_)
  {
    int scan = tag.getScan();
    if (scan_index_offset.find(scan) == scan_index_offset.end()) scan_index_offset[scan] = offset;
    offset++;
  }

  for (auto& tag : tags_)
  {
    int scan = tag.getScan();
    int index = tag.getIndex();
    tag.setIndex(index + scan_index_offset[scan]);
  }

  for (auto& hit : proteoform_hits_)
  {
    std::vector<int> tag_indices;
    if (hit.metaValueExists("TagIndices")) tag_indices = hit.getMetaValue("TagIndices").toIntList();
    int scan = hit.getMetaValue("Scan");
    for (int& tag_index : tag_indices)
    {
      tag_index += scan_index_offset[scan];
    }
    hit.setMetaValue("TagIndices", tag_indices);
  }

  markRepresentativeProteoformHits_(precursor_tol);
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), [](const ProteinHit& left, const ProteinHit& right) {
    return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getMetaValue("Scan") > right.getMetaValue("Scan"))
                                                                                            : (left.getCoverage() > right.getCoverage()))
                                               : (left.getScore() > right.getScore());
  });

  if (decoy_factor_ > 0 || !keep_underdetermined_)
  {
    std::vector<ProteinHit> filtered_proteoform_hits;
    filtered_proteoform_hits.reserve(proteoform_hits_.size());

    for (int k = 0; k < (keep_underdetermined_ ? 2 : 1); k++)
    {
      for (Size mod = 0; mod <= max_mod_cntr; mod++)
      {
        double taget_count = 0;
        double decoy_count = 0;

        double taget_count_pro = 0;
        double decoy_count_pro = 0;

        std::map<double, double> map_qvalue;
        std::map<double, double> map_qvalue_pro;

        for (auto& hit : proteoform_hits_)
        {
          std::vector<double> mod_masses = hit.getMetaValue("Modifications");
          if (mod_masses.size() != mod) continue;
          if (k == 0
              && ((double)hit.getMetaValue("Mass") < 0 || (int)hit.getMetaValue("StartPosition") < 0 || (int)hit.getMetaValue("EndPosition") < 0))
          {
            continue;
          }
          else if (k == 1
                   && ((double)hit.getMetaValue("Mass") > 0 && (int)hit.getMetaValue("StartPosition") > 0
                       && (int)hit.getMetaValue("EndPosition") > 0))
            continue;

          bool is_decoy = hit.getAccession().hasPrefix("DECOY");
          bool is_rep = hit.metaValueExists("Representative");
          if (is_decoy)
          {
            decoy_count += 1.0 / decoy_factor_;
            if (is_rep) decoy_count_pro += 1.0 / decoy_factor_;
          }
          else
          {
            taget_count++;
            if (is_rep) taget_count_pro ++;
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
          std::vector<double> mod_masses = hit.getMetaValue("Modifications");
          if (mod_masses.size() != mod) continue;
          if (k == 0
              && ((double)hit.getMetaValue("Mass") < 0 || (int)hit.getMetaValue("StartPosition") < 0 || (int)hit.getMetaValue("EndPosition") < 0))
            continue;
          else if (k == 1
                   && ((double)hit.getMetaValue("Mass") > 0 && (int)hit.getMetaValue("StartPosition") > 0
                       && (int)hit.getMetaValue("EndPosition") > 0))
            continue;

          bool is_decoy = hit.getAccession().hasPrefix("DECOY");
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

  //  define proteoform index and tag to proteoform indices.
  int proteoform_index = 0;
  for (auto& hit : proteoform_hits_)
  {
    hit.setMetaValue("Index", proteoform_index);
    for (int tag_index : (std::vector<int>)hit.getMetaValue("TagIndices").toIntList())
    {
      matching_hits_indices_[tag_index].push_back(proteoform_index);
    }
    proteoform_index++;
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