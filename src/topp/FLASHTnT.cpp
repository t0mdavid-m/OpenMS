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
#include <OpenMS/ANALYSIS/TOPDOWN/ProteoformDatabase.h>
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

    registerInputFile_("fasta", "<file>", "",
                       "Input proteome database file (FASTA or PEFF). A PEFF file (PSI Extended FASTA, e.g. from UniPEFF; first line '# PEFF' or "
                       "the extension .peff) gives every annotated proteoform of a protein as a database entry: its processed forms (mature chains, "
                       "without the initiator Met, signal, transit peptide or propeptide) with up to peff_max_mods annotated modifications fixed "
                       "on them (a residue with any of its modifications on distinct sites of it, or all disulfide bonds of a form, counts as one); variable and blind "
                       "modifications are searched on top. Every target protein without a decoy (DECOY prefix) in the file gets one: its residues "
                       "shuffled within the processed forms, their first and last residues kept, the modifications moved with their residues.",
                       true);
    setValidFormats_("fasta", ListUtils::create<String>("fasta"));

    registerIntOption_("peff_max_mods", "<number>", 1,
                       "PEFF input: at most this many annotated modification units fixed on an entry (0: the processed forms only).", false);
    setMinInt_("peff_max_mods", 0);
    registerStringList_("peff_mod_masses", "<accession=mass>",
                        ListUtils::create<String>("MOD:01819=100.016044,MOD:01893=86.000394,MOD:01614=329.05252,MOD:01956=15.994915,"
                                                  "MOD:00657=14.999666,MOD:00017=skip,MOD:00223=skip,MOD:01976=skip,MOD:01780=skip"),
                        "PEFF input: masses of modification accessions that replace those of the modification database (PSI-MOD, Unimod), or "
                        "'skip' to leave the modification out. The default fills the gaps of the bundled PSI-MOD for E. coli (UniProt ptmlist "
                        "masses). A modification without a mass is left out (counted in the log).",
                        false);

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

    auto tnt_param = FLASHTnTAlgorithm().getDefaults();
    tnt_param.removeAll("tag:");
    tnt_param.removeAll("ex:");

    registerFullParam_(tnt_param);

    //registerSubsection_("tnt", "TnT algorithm parameters");
    registerSubsection_("tag", "Tag algorithm parameters");
    registerSubsection_("ex", "Extension algorithm parameters");
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
    tnt_param.remove("peff_max_mods");
    tnt_param.remove("peff_mod_masses");
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

    double pro_fdr = tnt_param.getValue("pro_fdr");
    // the chimeric columns are written only with the residual search, so that other outputs keep their format
    const bool chimeric = (int)tnt_param.getValue("chimeric_rounds") > 0 || tnt_param.getValue("native").toString() == "true"; // -native searches a residual round
    OPENMS_LOG_INFO << "Finding sequence tags from deconvolved spectra ..." << endl;

    std::vector<FASTAFile::FASTAEntry> fasta_entry;
    // a PEFF file: its proteoform entries (ProteoformDatabase); a FASTA file: its entries are the proteins
    ProteoformDatabase database;
    const bool peff = ProteoformDatabase::isPEFF(in_fasta);
    if (peff)
    {
      database.loadPEFF(in_fasta, getIntOption_("peff_max_mods"), getStringList_("peff_mod_masses"), fasta_entry);
      const auto& st = database.getLoadStatistics();
      Size targets = 0, with_mods = 0;
      for (const auto& p : database.proteins) targets += ! p.decoy;
      for (const auto& e : database.entries) with_mods += ! e.fixed.empty();
      OPENMS_LOG_INFO << "PEFF database: " << st.input_targets << " target and " << st.input_decoys << " decoy proteins in the file, "
                      << st.generated_decoys << " decoys generated; " << targets << " target and " << database.proteins.size() - targets
                      << " decoy proteins, " << database.entries.size() << " entries (" << with_mods << " with fixed modifications; "
                      << st.forms << " processed forms of the targets, peff_max_mods " << getIntOption_("peff_max_mods") << "); "
                      << st.modification_sites << " modification sites and " << st.disulfide_bonds << " disulfide bonds of the targets used"
                      << endl;
      if (st.unpaired_half_cystines > 0)
      {
        OPENMS_LOG_WARN << "PEFF database: " << st.unpaired_half_cystines
                        << " half cystines without a disulfide bond pair are left out (write the PEFF with UniPEFF -annotation_identifiers)" << endl;
      }
      if (st.length_mismatches > 0)
      {
        OPENMS_LOG_WARN << "PEFF database: " << st.length_mismatches << " entries whose \\Length differs from their sequence" << endl;
      }
      if (st.empty_sequences > 0) { OPENMS_LOG_WARN << "PEFF database: " << st.empty_sequences << " entries without residues left out" << endl; }
      if (st.unpaired_decoys > 0)
      {
        OPENMS_LOG_WARN << "PEFF database: " << st.unpaired_decoys << " decoys of the file without a target of their name or accession" << endl;
      }
      if (st.unmirrored_decoys > 0)
      {
        OPENMS_LOG_WARN << "PEFF database: " << st.unmirrored_decoys
                        << " decoys of the file whose entries do not mirror their target's (their number or masses); without decoys in the "
                           "file, mirrored ones are generated"
                        << endl;
      }
      auto siteList = [](const std::map<String, Size>& counts, Size& sites) {
        String list;
        sites = 0;
        for (const auto& [acc, n] : counts)
        {
          sites += n;
          list += (list.empty() ? "" : ", ") + acc + " (" + String(n) + ")";
        }
        return list;
      };
      if (! st.unresolved.empty())
      {
        Size sites = 0;
        const String list = siteList(st.unresolved, sites);
        OPENMS_LOG_INFO << "PEFF database: " << sites << " modification sites without a mass left out (peff_mod_masses): " << list << endl;
      }
      if (! st.wrong_residue.empty())
      {
        Size sites = 0;
        const String list = siteList(st.wrong_residue, sites);
        OPENMS_LOG_WARN << "PEFF database: " << sites << " modification sites on a residue that is not the modification's origin left out: " << list
                        << endl;
      }
    }
    else
    {
      FASTAFile fasta_file;
      fasta_file.load(in_fasta, fasta_entry);
    }
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
      FLASHTnTFile::writePrSMHeader(out_prsm_stream, chimeric);
    }

    if (! out_pro_file.empty())
    {
      out_pro_stream = fstream(out_pro_file, fstream::out);
      FLASHTnTFile::writeProHeader(out_pro_stream, chimeric);
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
    tnt.run(map, fasta_entry, peff ? &database : nullptr);
    tnt.getProteoforms(proteoform_hits);

    OPENMS_LOG_INFO << "FLASHTnT run complete. Now writing the results in output files ..." << endl;
    if (! out_tag_file.empty())
    {
      /// TODO flanking_mass_tol should be calculated within FLASHTagger
      FLASHTnTFile::writeTags(tnt, tnt.getMaxTotalModificationMass(), out_tagger_stream);
      out_tagger_stream.close();
    }
    if (! out_prsm_file.empty())
    {
      FLASHTnTFile::writePrSMs(proteoform_hits, out_prsm_stream, chimeric);
      out_prsm_stream.close();
    }
    if (! out_pro_file.empty())
    {
      FLASHTnTFile::writeProteoforms(proteoform_hits, out_pro_stream, pro_fdr, chimeric);
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
