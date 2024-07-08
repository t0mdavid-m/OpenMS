// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

namespace OpenMS
{
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

  defaults_.setValue("allow_multi_hits", "false", "Allow multiple hits per spectrum");
  defaults_.setValidStrings("allow_multi_hits", {"true", "false"});

  defaults_.setValue("keep_underdetermined", "false", "To keep underdetermined proteoform hits (i.e., proteoforms without total mass or start end positions");
  defaults_.setValidStrings("keep_underdetermined", {"true", "false"});

  defaults_.setValue("keep_decoy", "false", "To keep decoy hits");
  defaults_.setValidStrings("keep_decoy", {"true", "false"});

  defaults_.insert("tag:", FLASHTaggerAlgorithm().getDefaults());
  defaults_.insert("ex:", FLASHExtenderAlgorithm().getDefaults());
  defaultsToParam_();
}

void FLASHTnTAlgorithm::updateMembers_()
{
  tagger_param_ = param_.copy("tag:", true);
  extender_param_ = param_.copy("ex:", true);
  prsm_fdr_ = param_.getValue("prsm_fdr");
  pro_fdr_ = param_.getValue("pro_fdr");
  keep_decoy_ = param_.getValue("keep_decoy").toString() == "true";
  keep_underdetermined_ = param_.getValue("keep_underdetermined").toString() == "true";
  multiple_hits_per_spec_ = param_.getValue("allow_multi_hits").toString() == "true";
}

bool FLASHTnTAlgorithm::areConsistent_(const ProteinHit& a, const ProteinHit& b, double tol) const
{
  //   hit.setMetaValue("ModificationIDs", mod_ids);
  //    hit.setMetaValue("ModificationACCs", mod_accs);
  //    hit.setMetaValue("Modifications", mod_masses);
  //    hit.setMetaValue("ModificationStarts", mod_starts);
  //    hit.setMetaValue("ModificationEnds", mod_ends);
  //    hit.setMetaValue("MatchedAA", total_match_cntr);
  //    hit.setCoverage((double)total_match_cntr / (double)hit.getSequence().size());
  //    hit.setScore(total_score);
  //    hit.setMetaValue("StartPosition", protein_start_position);
  //    hit.setMetaValue("EndPosition", protein_end_position);
  //    hit.setMetaValue("Mass", precursor_mass);
  //    hit.setMetaValue("RT", tagger.getSpectrum().getRT());

  double mass1 = a.getMetaValue("Mass");
  double mass2 = b.getMetaValue("Mass");
  if (std::abs(mass1 - mass2) > std::max(mass1, mass2) * tol / 1e6 * 2) return false;

  double rt1 = a.getMetaValue("RT");
  double rt2 = b.getMetaValue("RT");
  if (std::abs(rt1 - rt2) < 30.0) return true; // if rts are within 30 sec, true

  if (a.metaValueExists("mod_masses") && b.metaValueExists("mod_masses"))
  {
    std::vector<double> mod_masses1 = a.getMetaValue("mod_masses");
    std::vector<double> mod_masses2 = b.getMetaValue("mod_masses");
    if (mod_masses1.size() != mod_masses2.size()) return false;
    else return true;
  }
  else if(!a.metaValueExists("mod_masses") && !a.metaValueExists("mod_masses")) return true;
  // TODO check the modification locations etc.
  return false;
}


void FLASHTnTAlgorithm::markRepresentativeProteoformHits_(double tol)
{
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), [](const ProteinHit& left, const ProteinHit& right) {
    return left.getScore() > right.getScore(); });
  std::map<String, std::vector<ProteinHit>> proteoform_map;
  for (auto& hit : proteoform_hits_)
  {
    String acc = hit.getAccession();

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
    //std::cout<<acc << " " <<  (proteoform_map.find("sp|P02359|RS7_ECOLI") != proteoform_map.end()) << " " << tmp.size() << std::endl;
  }
}

void FLASHTnTAlgorithm::run(const MSExperiment& map, const std::vector<FASTAFile::FASTAEntry>& fasta_entry, double flanking_mass_tol)
{
  setLogType(CMD);
  startProgress(0, (SignedSize)map.size(), "Running FLASHTnT ...");
  int max_tag_length = tagger_param_.getValue("max_length");
  int min_tag_length = tagger_param_.getValue("min_length");
  int max_mod_cntr = extender_param_.getValue("max_mod_count");
  double max_diff_mass = (double)extender_param_.getValue("max_mod_mass") + 1.0;
  std::map<double, std::vector<ResidueModification>> mod_map;
  const auto inst = ModificationsDB::getInstance();  // give this from outside ...
  std::vector<String> mod_strs;
  inst->getAllSearchModifications(mod_strs);
  for (int i = 0; i < mod_strs.size(); i++)
  {
    const auto mod = *inst->getModification(mod_strs[i]);
    if (std::abs(mod.getDiffMonoMass()) > max_diff_mass) continue;
    mod_map[mod.getDiffMonoMass()].push_back(mod);
  }
  double precursor_tol = -1;
    // collect statistics for information
  for (int index = 0; index < map.size(); index++)
  {
    auto spec = map[index];
    nextProgress();
    int scan = FLASHDeconvAlgorithm::getScanNumber(map, index);

    //if (scan < 2248 || scan > 2253) continue;
    if (spec.getMSLevel() == 1 && precursor_tol > 0)
    {
      continue;
    }

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

    FLASHTaggerAlgorithm tagger;
    // Run tagger
    tagger.setParameters(tagger_param_);
    tagger.run(dspec, tol);
    FLASHExtenderAlgorithm extender;
    extender.setParameters(extender_param_);
    extender.setModificationMap(mod_map);
    tagger.getTags(tags_);

    if (multiple_hits_per_spec_)
    {
      tagger.runMatching(fasta_entry);
      extender.run(tagger, flanking_mass_tol, tol, multiple_hits_per_spec_);
      extender.getProteoforms(proteoform_hits_);
    }
    else
    {
      for (int tag_length = max_tag_length; tag_length >= min_tag_length; tag_length--)
      {
        tagger.runMatching(fasta_entry, tag_length);
        extender.run(tagger, flanking_mass_tol, tol, multiple_hits_per_spec_);
        extender.getProteoforms(proteoform_hits_);
        if (extender.hasProteoforms()) break;
      }
    }

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
    std::vector<int> tag_indices = hit.getMetaValue("TagIndices");
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

  if (decoy_factor_ > 0)
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
          if (k == 0 && ((int)hit.getMetaValue("StartPosition") < 0 || (int)hit.getMetaValue("EndPosition") < 0)) continue;
          else if (k == 1 && ((int)hit.getMetaValue("StartPosition") > 0 && (int)hit.getMetaValue("EndPosition") > 0)) continue;

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
            if (is_rep) taget_count_pro += 1.0 / decoy_factor_;
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
          if (k == 0 && ((int)hit.getMetaValue("StartPosition") < 0 || (int)hit.getMetaValue("EndPosition") < 0)) continue;
          else if (k == 1 && ((int)hit.getMetaValue("StartPosition") > 0 && (int)hit.getMetaValue("EndPosition") > 0)) continue;

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

  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), [](const ProteinHit& left, const ProteinHit& right) {
    return left.getMetaValue("RT") < right.getMetaValue("RT");
  });

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