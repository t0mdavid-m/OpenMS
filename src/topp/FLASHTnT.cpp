// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FLASHTaggerFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <QFileInfo>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
// We do not want this class to show up in the docu:

class TOPPFLASHTnT : public TOPPBase
{
public:
  TOPPFLASHTnT(): TOPPBase("FLASHTnT", "FLASHTnT to generate de novo sequence tags from TDP spectrum and match them against proteome DB for proteoform identification.", false)
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (deconv.mzML from FLASHDeconv mzML output)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("fasta", "<file>", "", "Input proteome database file (fasta)");
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));

    registerOutputFile_("out_protein", "<file>", "", "Default output protein level tsv file containing matched proteins");
    setValidFormats_("out_protein", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_tag", "<file>", "", "Default output tag level tsv file containing matched tags");
    setValidFormats_("out_tag", ListUtils::create<String>("tsv"));

    registerDoubleOption_("fdr", "value", 1.0, "Proteoform FDR threshold.", false, false);
    setMinFloat_("fdr", 0.01);
    setMaxFloat_("fdr", 1.0);

    registerFlag_("keep_decoy", "Keep decoy proteoforms.", false);

    registerSubsection_("TAG", "FLASHTagger algorithm parameters");
    registerSubsection_("EX", "FLASHExtender algorithm parameters");
  }

  Param getSubsectionDefaults_(const String& prefix) const override
  {
    if (prefix == "TAG")
    {
      auto tagger_param = FLASHTaggerAlgorithm().getDefaults();
      return tagger_param;
    }else if (prefix == "EX")
    {
      auto tagger_param = FLASHExtenderAlgorithm().getDefaults();
      return tagger_param;
    }
    else { throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", prefix); }
  }

  // the main_ function is called after all parameters are read
  ExitCodes main_(int, const char**) override
  {
    OPENMS_LOG_INFO << "Initializing ... " << endl;
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------

    String in_file = getStringOption_("in");
    String in_fasta = getStringOption_("fasta");

    String out_tag_file = getStringOption_("out_tag");
    String out_protein_file = getStringOption_("out_protein");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    MSExperiment map;
    MzMLFile mzml;

    // reading mzMLs with m/z and rt criteria.

    mzml.setLogType(log_type_);
    mzml.load(in_file, map);

    auto tagger_param = getParam_().copy("TAG:", true);
    auto extender_param = getParam_().copy("EX:", true);
    OPENMS_LOG_INFO << "Finding sequence tags from deconvolved MS2 spectra ..." << endl;

    FASTAFile fasta_file;
    std::vector<FASTAFile::FASTAEntry> fasta_entry;
    fasta_file.load(in_fasta, fasta_entry);
    std::vector<ProteinHit> proteoform_hits;
    double decoy_factor = 0;
    double fdr = getDoubleOption_("fdr");
    bool keep_decoy = getFlag_("keep_decoy");
    double tol(10);
    // Run FLASHDeconvAlgorithm here!
    OPENMS_LOG_INFO << "Processing : " << in_file << endl;

    fstream out_tagger_stream;
    fstream out_protein_stream;

    if (! out_tag_file.empty())
    {
      out_tagger_stream = fstream(out_tag_file, fstream::out);
      FLASHTaggerFile::writeTagHeader(out_tagger_stream);
    }

    if (! out_protein_file.empty())
    {
      out_protein_stream = fstream(out_protein_file, fstream::out);
      FLASHTaggerFile::writeProteinHeader(out_protein_stream);
    }
    // collect statistics for information
    for (int index = 0; index < map.size(); index++)
    {
      auto spec = map[index];
      if (spec.getMSLevel() != 2) continue;
      int scan = FLASHDeconvAlgorithm::getScanNumber(map, index);
      DeconvolvedSpectrum dspec(scan);
      dspec.setOriginalSpectrum(spec);
      String deconv_meta_str = spec.getMetaValue("DeconvMassInfo").toString();
      int tol_loc_s = deconv_meta_str.find("tol=") + 4;
      int tol_loc_e = deconv_meta_str.find(";", tol_loc_s);
      tol = stod(deconv_meta_str.substr(tol_loc_s, tol_loc_e - tol_loc_s));

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
      tagger.setParameters(tagger_param);
      tagger.run(dspec, tol, fasta_entry);
      decoy_factor = tagger.getDecoyFactor();
      FLASHExtenderAlgorithm extender;

      extender.setParameters(extender_param);
      extender.run(tagger, tol);
      extender.getProteoforms(proteoform_hits);
      if (! out_tag_file.empty())
      {
        FLASHTaggerFile::writeTags(tagger, extender, out_tagger_stream);
      }
    }

    std::sort(proteoform_hits.begin(), proteoform_hits.end(),
              [](const ProteinHit& left, const ProteinHit& right) {
                return left.getScore() > right.getScore();
              });

    double taget_count = 0;
    double decoy_count = 0;
    std::vector<ProteinHit> filtered_proteoform_hits;
    filtered_proteoform_hits.reserve(proteoform_hits.size());

    for (auto& hit : proteoform_hits)
    {
      bool is_decoy = hit.getAccession().hasPrefix("DECOY");
      if (is_decoy) decoy_count++;
      else taget_count++;

      double qvalue = decoy_factor != 0 ? (decoy_count / (decoy_count + taget_count)) : -1.0; // TODO make the minima thing
      hit.setMetaValue("qvalue", qvalue);
      if (!keep_decoy && is_decoy) continue;
      if (fdr < 1 && qvalue > fdr) continue;
      filtered_proteoform_hits.push_back(hit);
    }

    if (! out_protein_file.empty())
    {
      FLASHTaggerFile::writeProteins(filtered_proteoform_hits, out_protein_stream);
      out_protein_stream.close();
    }

    if (! out_tag_file.empty())
    {
      out_tagger_stream.close();
    }

    OPENMS_LOG_INFO << "FLASHTnT run complete. Now writing the results in output files ..." << endl;

    return EXECUTION_OK;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  TOPPFLASHTnT tool;
  return tool.main(argc, argv);
}
