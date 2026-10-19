// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FLASHTnTFile.h>
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
  TOPPFLASHTnT():
      TOPPBase("FLASHTnT",
               "FLASHTnT to generate de novo sequence tags from TDP spectrum and match them against proteome DB for proteoform identification.",
               false)
  {
  }

protected:
  // this function will be used to register the tool parameters
  // it gets automatically called on tool execution
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input file (deconv.mzML from FLASHDeconv mzML output).", true);
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerInputFile_("fasta", "<file>", "", "Input proteome database file (fasta)", true);
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));

    registerInputFile_("in_raw", "<file>", "",
                       "Raw mzML file (profile or centroid) FLASHDeconv was run on. With it, every PrSM (targets and decoys) is annotated on its raw "
                       "spectrum (fragment isotope envelopes at every charge, as FLASHAnnotator) before the FDR, and ranked by its extension score "
                       "plus -raw_score_weight times its raw score. The raw annotation columns are written in the PrSM and proteoform files.",
                       false);
    setValidFormats_("in_raw", ListUtils::create<String>("mzML"));

    registerOutputFile_("out_pro", "<file>", "", "Output proteoform-level tsv file containing proteoform IDs. ", true);
    setValidFormats_("out_pro", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_prsm", "<file>", "", "Output PrSM-level tsv file containing PrSMs.");
    setValidFormats_("out_prsm", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_tag", "<file>", "", "Output tag-level tsv file containing matched tags.");
    setValidFormats_("out_tag", ListUtils::create<String>("tsv"));

    registerOutputFile_("out_protein", "<file>", "", "Output protein-level tsv file containing the candidate proteins and "
                        "their scores, as ranked before proteoform characterization.", false);
    setValidFormats_("out_protein", ListUtils::create<String>("tsv"));

    // Register proteoform-level FDR parameter (no PrSM-level filter: the reported PrSMs are those of the reported proteoforms)
    registerDoubleOption_("pro_fdr", "Specifies the proteoform-level FDR.", 1.0,
                          "Specifies the proteoform-level FDR; proteins are reported at the same level. The reported PrSMs are those of the "
                          "reported proteoforms (there is no PrSM-level filter).",
                          false);
    setMinFloat_("pro_fdr", 0.0);

    // Register single-hit-only option
    registerFlag_("only_single_hit", "Allows only a single hit per spectrum.", "Allows only a single hit per spectrum.");
    //setValidStrings_("only_single_hit", {"true", "false"});

    // Register underdetermined proteoform discard option
    registerFlag_("discard_underdetermined",
                          "Discards underdetermined proteoform IDs (e.g., those without exact precursor masses or start/end positions).");
    //setValidStrings_("discard_underdetermined", {"true", "false"});

    // Register protein-level-only option
    registerFlag_("protein_level_only",
                  "Stops after the protein-level search, skipping proteoform characterization. Use with -out_protein to "
                  "shortlist proteins cheaply.");

    registerIntOption_("max_protein_hits", "<number>", 50,
                       "Number of candidate proteins kept per spectrum by the protein-level search. When proteoform "
                       "search runs, only these are characterized.", false);
    setMinInt_("max_protein_hits", 1);

    registerDoubleOption_("protein_tol", "<ppm>", 0.0,
                          "Mass tolerance (ppm) of the protein-level search: a deconvolved mass, moved by a tag-derived shift, "
                          "matches a protein residue sum within it. 0: the MS2 deconvolution tolerance recorded in the input file.",
                          false);
    setMinFloat_("protein_tol", 0.0);

    registerIntOption_("min_protein_hits", "<number>", 5,
                       "Number of best-scoring candidate proteins per spectrum that are always characterized. The other "
                       "candidates are characterized only if their protein, or its target-decoy partner, passes -protein_fdr. Used only when the database "
                       "contains decoys (DECOY prefix).", false);
    setMinInt_("min_protein_hits", 1);

    registerDoubleOption_("protein_fdr", "<value>", 0.1,
                          "Run-level protein FDR, from each protein's best tag-based score, that a candidate beyond "
                          "-min_protein_hits (or its target-decoy partner) must pass to be characterized. Needs decoys in the database; without them all "
                          "candidates kept by -max_protein_hits are characterized. 1 characterizes them all.", false);
    setMinFloat_("protein_fdr", 0.0);
    setMaxFloat_("protein_fdr", 1.0);

    // Register decoy retention option
    registerFlag_("keep_decoy", "Retains decoy hits in the results.");
    //setValidStrings_("keep_decoy", {"true", "false"});

    // Register ion type parameter
    registerStringList_("ion_type", "Specifies ion types to consider.", ListUtils::toStringList(std::vector<std::string>{"b", "y"}), "Specifies ion types to consider.", false);
    setValidStrings_("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

    registerStringOption_("gapped_tags", "<mode>", "fallback",
                          "Candidate proteins from gapped tags (masses whose differences are sums of one or more residues): off; fallback "
                          "(spectra in which no sequence tag proposes a protein); both (every spectrum, proteins not proposed by sequence tags); "
                          "only (gapped tags instead of sequence tags). Where a gapped tag lands on a protein anchors the extension as a sequence "
                          "tag does.", false);
    setValidStrings_("gapped_tags", {"off", "fallback", "both", "only"});
    registerStringOption_("terminal_ladders", "<bool>", "true",
                          "Candidate proteins from terminal ladders: deconvolved masses that are residue sums from an intact protein terminus (the "
                          "N-terminus also without the initiator Met and with acetylation), found without a sequence tag. In every spectrum the "
                          "min_protein_hits best proteins by the chance probability of their ladder are characterized, anchored on the ladder.", false);
    setValidStrings_("terminal_ladders", {"true", "false"});
    registerStringOption_("precursor_ladders", "<bool>", "true",
                          "Candidate proteins from precursor-anchored ladders: in each MSn spectrum with a precursor mass, every stretch of a protein "
                          "(between non-standard residues) whose unmodified mass equals the precursor mass (within the MS1 deconvolution tolerance of "
                          "the input, or -protein_tol's tolerance if the input has no MS1 spectra; isotope errors -1 to +1) is scored by the chance "
                          "probability of the deconvolved masses on its prefix and suffix fragment ladders (a position counts once per side, whatever "
                          "ion types match it). The stretches are unmodified: N-terminal acetylation and other modifications are not proposed "
                          "(initiator Met removal is a stretch from residue 2). The -min_protein_hits best stretches are characterized from tags on "
                          "their ladders. Only their PrSMs of determined mass and termini that are not the same proteoform as a determined PrSM of "
                          "another source in the spectrum are kept; they form their own FDR class (\"precursor ladder\"). Not used with "
                          "-protein_level_only. Without any precursor mass in the input, results are unchanged; otherwise the PrSMs of this source can "
                          "also change the others' results (the pooled protein-level null, the best PrSM of a spectrum, proteoform clusters).", false);
    setValidStrings_("precursor_ladders", {"true", "false"});
    registerDoubleOption_("raw_score_weight", "<weight>", 8.0,
                          "With -in_raw: every PrSM is ranked for the FDR by its extension score plus this weight times its raw score (-log10 of the "
                          "chance of the fragment isotope envelopes matched in the raw spectrum). 0: the raw score is reported but not used in the ranking.",
                          false);
    setMinFloat_("raw_score_weight", 0.0);
    registerStringOption_("learn_var_mod_penalty", "<bool>", "true",
                          "Raises the penalty of a variable modification that target PrSMs use per opportunity (fitting residue or terminus) no more often than decoy PrSMs "
                          "(chance), up to that of a blind shift. The usage is counted in a pilot pass over every 4th spectrum (all spectra if the "
                          "pilot has fewer than 50 decoy PrSMs); with a raised penalty, all spectra are then extended once with the learned penalties. "
                          "Needs decoys and at least 50 decoy PrSMs.", false, true);
    setValidStrings_("learn_var_mod_penalty", {"true", "false"});
    registerFlag_("decoy_precursor",
                  "Decoy precursors: spectra whose PrSM took its mass from MS1 are searched again with the decoy precursor mass of their isolation "
                  "window instead of the precursor mass (FLASHDeconv -FD:report_FDR true writes it); a PrSM reporting that mass competes with the "
                  "spectrum's own and counts as a decoy in the FDR.");

    registerSubsection_("tag", "Tag algorithm parameters");
    registerSubsection_("ex", "Extension algorithm parameters");
    registerSubsection_("gapped", "Gapped tag parameters (with -gapped_tags)");
    registerSubsection_("raw", "Raw spectrum annotation parameters (with -in_raw)");

  }

  Param getSubsectionDefaults_(const String& prefix) const override
  {
    if (prefix == "tag")
    {
      auto tnt_param = FLASHTnTAlgorithm().getDefaults();
      return tnt_param.copy("tag:", true);
    }
    else if (prefix == "ex") {
      auto tnt_param = FLASHTnTAlgorithm().getDefaults();
      return tnt_param.copy("ex:", true);
    }
    else if (prefix == "gapped")
    {
      auto tnt_param = FLASHTnTAlgorithm().getDefaults();
      return tnt_param.copy("gapped:", true);
    }
    else if (prefix == "raw")
    {
      auto tnt_param = FLASHTnTAlgorithm().getDefaults();
      return tnt_param.copy("raw:", true);
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
    String in_raw = getStringOption_("in_raw");

    String out_tag_file = getStringOption_("out_tag");
    String out_prsm_file = getStringOption_("out_prsm");
    String out_pro_file = getStringOption_("out_pro");
    String out_protein_file = getStringOption_("out_protein");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    MSExperiment map;
    MzMLFile mzml;

    // reading mzMLs with m/z and rt criteria.

    mzml.setLogType(log_type_);
    mzml.load(in_file, map);

    auto tnt_param = getParam_();
    tnt_param.remove("in");
    tnt_param.remove("fasta");
    tnt_param.remove("in_raw");
    tnt_param.remove("out_tag");
    tnt_param.remove("out_prsm");
    tnt_param.remove("out_pro");
    tnt_param.remove("out_protein");
    tnt_param.remove("log");
    tnt_param.remove("debug");
    tnt_param.remove("threads");
    tnt_param.remove("no_progress");
    tnt_param.remove("force");
    tnt_param.remove("test");
    double max_mod_mass = tnt_param.getValue("ex:max_mod_mass");
    int max_mod_count = tnt_param.getValue("ex:max_mod_count");
    double pro_fdr = tnt_param.getValue("pro_fdr");
    OPENMS_LOG_INFO << "Finding sequence tags from deconvolved spectra ..." << endl;

    FASTAFile fasta_file;
    std::vector<FASTAFile::FASTAEntry> fasta_entry;
    fasta_file.load(in_fasta, fasta_entry);
    std::vector<ProteinHit> proteoform_hits;
    // Run here!
    OPENMS_LOG_INFO << "Processing : " << in_file << endl;

    fstream out_tagger_stream;
    fstream out_prsm_stream;
    fstream out_pro_stream;
    fstream out_protein_stream;

    if (! out_tag_file.empty())
    {
      out_tagger_stream = fstream(out_tag_file, fstream::out);
      FLASHTnTFile::writeTagHeader(out_tagger_stream);
    }

    if (! out_prsm_file.empty())
    {
      out_prsm_stream = fstream(out_prsm_file, fstream::out);
      FLASHTnTFile::writePrSMHeader(out_prsm_stream);
    }

    if (! out_pro_file.empty())
    {
      out_pro_stream = fstream(out_pro_file, fstream::out);
      FLASHTnTFile::writeProHeader(out_pro_stream);
    }

    if (! out_protein_file.empty())
    {
      out_protein_stream = fstream(out_protein_file, fstream::out);
      FLASHTnTFile::writeProteinHeader(out_protein_stream);
    }

    FLASHTnTAlgorithm tnt;
    tnt.setParameters(tnt_param);
    if (! in_raw.empty())
    {
      tnt.setRawSpectrumLoader([&in_raw](const std::set<int>& scans, std::map<int, MSSpectrum>& spectra) {
        FLASHAnnotatorAlgorithm::loadSpectra(in_raw, scans, spectra);
      });
    }
    tnt.run(map, fasta_entry);
    tnt.getProteoforms(proteoform_hits);

    OPENMS_LOG_INFO << "FLASHTnT run complete. Now writing the results in output files ..." << endl;
    if (! out_tag_file.empty())
    {
      FLASHTnTFile::writeTags(tnt, max_mod_count * max_mod_mass + 1, out_tagger_stream);
      out_tagger_stream.close();
    }
    if (! out_prsm_file.empty())
    {
      FLASHTnTFile::writePrSMs(proteoform_hits, out_prsm_stream);
      out_prsm_stream.close();
    }
    if (! out_pro_file.empty())
    {
      FLASHTnTFile::writeProteoforms(proteoform_hits, out_pro_stream, pro_fdr);
      out_pro_stream.close();
    }
    if (! out_protein_file.empty())
    {
      std::vector<ProteinHit> protein_hits;
      tnt.getProteinHits(protein_hits);
      FLASHTnTFile::writeProteins(protein_hits, out_protein_stream);
      out_protein_stream.close();
    }

    return EXECUTION_OK;
  }
};

// the actual main function needed to create an executable
int main(int argc, const char** argv)
{
  TOPPFLASHTnT tool;
  return tool.main(argc, argv);
}
