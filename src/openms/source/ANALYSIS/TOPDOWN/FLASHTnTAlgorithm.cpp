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
  defaults_.setValue("fdr", 1.0, "Protein level FDR");
  defaults_.setMinFloat("fdr", 0.0);

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
  fdr_ = param_.getValue("fdr");
  keep_decoy_ = param_.getValue("keep_decoy").toString() == "true";
}

void FLASHTnTAlgorithm::run(const MSExperiment& map, const std::vector<FASTAFile::FASTAEntry>& fasta_entry, double flanking_mass_tol)
{
  setLogType(CMD);
  startProgress(0, (SignedSize)map.size(), "Running FLASHTnT ...");
  // collect statistics for information
  for (int index = 0; index < map.size(); index++)
  {
    auto spec = map[index];
    nextProgress();
    if (spec.getMSLevel() == 1) continue;
    int scan = FLASHDeconvAlgorithm::getScanNumber(map, index); // TODO precursor

    //if (scan != 8621) continue; //
    DeconvolvedSpectrum dspec(scan);
    dspec.setOriginalSpectrum(spec);
    String deconv_meta_str = spec.getMetaValue("DeconvMassInfo").toString();
    int tol_loc_s = deconv_meta_str.find("tol=") + 4;
    int tol_loc_e = deconv_meta_str.find(";", tol_loc_s);
    double tol = stod(deconv_meta_str.substr(tol_loc_s, tol_loc_e - tol_loc_s));

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
    tagger.run(dspec, tol, fasta_entry);
    tagger.getTags(tags_);
    decoy_factor_ = tagger.getDecoyFactor();
    FLASHExtenderAlgorithm extender;

    extender.setParameters(extender_param_);
    extender.run(tagger, flanking_mass_tol, tol);
    extender.getProteoforms(proteoform_hits_);
    // if (! out_tag_file.empty()) { FLASHTaggerFile::writeTags(tagger, extender, out_tagger_stream); }
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

  //std::cout<<tags_.size()<<" " << proteoform_hits_.size()<<std::endl;
  std::sort(proteoform_hits_.begin(), proteoform_hits_.end(), [](const ProteinHit& left, const ProteinHit& right) {
    return left.getScore() == right.getScore() ? (left.getCoverage() == right.getCoverage() ? (left.getMetaValue("Scan") > right.getMetaValue("Scan"))
                                                                                            : (left.getCoverage() > right.getCoverage()))
                                               : (left.getScore() > right.getScore());
  });

  if (decoy_factor_ > 0)
  {
    double taget_count = 0;
    double decoy_count = 0;

    std::vector<ProteinHit> filtered_proteoform_hits;
    filtered_proteoform_hits.reserve(proteoform_hits_.size());
    std::map<double, double> map_qvalue;

    for (auto& hit : proteoform_hits_)
    {
      bool is_decoy = hit.getAccession().hasPrefix("DECOY");
      if (is_decoy) decoy_count += 1.0 / decoy_factor_;
      else
        taget_count++;

      double tmp_qvalue = decoy_count / (decoy_count + taget_count);
      map_qvalue[hit.getScore()] = std::min(1.0, tmp_qvalue);
    }

    double cummin = 1.0;
    for (auto&& rit = map_qvalue.begin(); rit != map_qvalue.end(); ++rit)
    {
      cummin = std::min(rit->second, cummin);
      rit->second = cummin;
    }

    for (auto& hit : proteoform_hits_)
    {
      bool is_decoy = hit.getAccession().hasPrefix("DECOY");
      double qvalue = map_qvalue[hit.getScore()];
      hit.setMetaValue("qvalue", qvalue);
      if (! keep_decoy_ && is_decoy) continue;
      if (fdr_ < 1 && qvalue > fdr_) continue;

      filtered_proteoform_hits.push_back(hit);
    }

    proteoform_hits_.swap(filtered_proteoform_hits);
  }
  //std::cout<<tags_.size()<<" " << proteoform_hits_.size()<<std::endl; //
  // define proteoform index and tag to proteoform indices.
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