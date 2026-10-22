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
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FLASHTnTFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/UniPEFFConverter.h>
#include <OpenMS/SYSTEM/File.h>
#include <QFileInfo>
#include <cctype>
#include <fstream>
#include <set>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
// We do not want this class to show up in the docu:

namespace
{
  /// the format of a protein database file
  enum class DatabaseFormat
  {
    FASTA,
    PEFF,
    UNIPROT_XML
  };

  /// the format of @p file by its content (gzip: a UniProtKB .xml.gz download; '# PEFF'; an XML declaration or <uniprot>), else by its
  /// name (.peff; .xml); FASTA otherwise
  DatabaseFormat databaseFormat_(const String& file)
  {
    std::ifstream in(file, std::ios::binary);
    char magic[2] = {0, 0};
    in.read(magic, 2);
    if (in.gcount() == 2 && (unsigned char)magic[0] == 0x1f && (unsigned char)magic[1] == 0x8b) return DatabaseFormat::UNIPROT_XML;
    in.clear();
    in.seekg(0);
    String line;
    while (std::getline(in, line))
    {
      line.trim();
      if (! line.empty()) break;
    }
    if (line.hasPrefix("# PEFF")) return DatabaseFormat::PEFF;
    if (line.hasPrefix("<?xml") || line.hasPrefix("<uniprot")) return DatabaseFormat::UNIPROT_XML;
    const String name = String(file).toLower();
    if (name.hasSuffix(".peff")) return DatabaseFormat::PEFF;
    if (name.hasSuffix(".xml") || name.hasSuffix(".xml.gz")) return DatabaseFormat::UNIPROT_XML;
    return DatabaseFormat::FASTA;
  }
} // namespace

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

    registerInputFile_("db", "<file>", "",
                       "Protein database. Three formats, recognised from the file content (else from the file name):\n"
                       "- FASTA (e.g. .fasta, .fa): one entry per protein. The file must contain the decoys (accessions with the prefix DECOY_).\n"
                       "- PEFF, PSI Extended FASTA (first line '# PEFF', or .peff; e.g. written by UniPEFF -annotation_identifiers): every annotated "
                       "proteoform is a database entry, i.e. each processed form of a protein (full length, mature chains, without the initiator Met, "
                       "signal, transit peptide or propeptide) with up to -peff_max_mods annotated modification units fixed on it (\\ModResPsi, "
                       "\\ModResUnimod; disulfide bonds from \\DisulfideBond pairs). Variable and blind modifications are searched on top. A target "
                       "without a decoy in the file gets one (its residues shuffled within the processed forms, the ends of each form kept, the fixed "
                       "modifications moved with their residues); decoys in the file are used as given. Sequence variants are not used. Every "
                       "entry that the tags, gapped tags or mass ladders support is extended, and the protein keeps the PrSM of a determined mass, "
                       "then the fewest blind modifications, the best score and the most annotated modifications. With -max_blind_mod_count 0 only "
                       "the annotated (and variable) modifications are used: a stricter search (E. coli lysates: about half the proteoforms, almost "
                       "no entrapment hits).\n"
                       "- UniProtKB XML (.xml or .xml.gz, e.g. a UniProt proteome download in XML format): converted on the fly to PEFF as "
                       "'UniPEFF -annotation_identifiers' writes it, then read as PEFF.\n"
                       "A UniProtKB XML file gives the same results as the PEFF converted from it.",
                       true);

    registerStringList_("fixed_mod", "<name>", std::vector<String>(),
                        "Fixed modifications, by name with the residue (e.g. 'Carbamidomethyl (C)' for alkylated cysteines, 'TMT6plex (K)') or by Unimod "
                        "accession number with the residue (e.g. '4 (C)'; 'UNIMOD:4 (C)' is the same). "
                        "The residue's mass is changed for the whole search: every occurrence in targets and decoys, in sequence tags (a tag reads "
                        "the modified residue as its letter), mass ladders, the extension, the raw-spectrum scoring and the reported masses. Only "
                        "modifications of a residue anywhere, one per residue (not terminal ones). Written as ProForma global modifications, "
                        "e.g. <[+57.021464|UNIMOD:4]@C>. With a PEFF database an annotated modification of such a residue adds to it.",
                        false);

    registerIntOption_("peff_max_mods", "<number>", 1,
                       "PEFF or UniProtKB XML database: at most this many annotated modification units fixed on an entry (0: the processed forms "
                       "only). A unit is a residue with one or more of its annotated modifications on distinct sites of it, or all disulfide bonds of "
                       "a form.", false);
    setMinInt_("peff_max_mods", 0);
    registerStringList_("peff_mod_masses", "<accession=mass>",
                        ListUtils::create<String>("MOD:01819=100.016044,MOD:01893=86.000394,MOD:01614=329.05252,MOD:01956=15.994915,"
                                                  "MOD:00657=14.999666,MOD:00017=skip,MOD:00223=skip,MOD:01976=skip,MOD:01780=skip"),
                        "PEFF or UniProtKB XML database: masses of modification accessions that replace those of the modification database (PSI-MOD, Unimod), or "
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

    registerOutputFile_("out_tag", "<file>", "", "Output tag-level tsv file containing matched tags.", false);
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
    const String in_db = getStringOption_("db");

    // Fixed modifications change the residue masses of the search, the in-process amino acid set (ResidueDB), before anything reads them:
    // every mass of the search -- tags, ladders, the extension, the raw annotation, the reported masses -- is that residue's.
    String fixed_proforma;
    {
      std::set<char> fixed_residues;
      for (const String& name : getStringList_("fixed_mod"))
      {
        if (name.empty()) continue;
        const ResidueModification* mod = nullptr;
        try
        {
          mod = FLASHHelperClasses::resolveModification(name); // a name or a Unimod accession, with its residue
        }
        catch (const Exception::BaseException& e)
        {
          OPENMS_LOG_ERROR << "-fixed_mod " << name << ": " << e.what() << endl;
          return ILLEGAL_PARAMETERS;
        }
        const char origin = mod->getOrigin();
        if (mod->getTermSpecificity() != ResidueModification::ANYWHERE || ! std::isupper((unsigned char)origin) || origin == 'X')
        {
          OPENMS_LOG_ERROR << "-fixed_mod " << name << ": only a modification of a residue anywhere can be fixed (not a terminal one)." << endl;
          return ILLEGAL_PARAMETERS;
        }
        if (! fixed_residues.insert(origin).second)
        {
          OPENMS_LOG_ERROR << "-fixed_mod: two fixed modifications of residue " << origin << "." << endl;
          return ILLEGAL_PARAMETERS;
        }
        Residue* residue = const_cast<Residue*>(ResidueDB::getInstance()->getResidue((unsigned char)origin));
        const double before = residue->getMonoWeight(Residue::Internal);
        residue->setMonoWeight(residue->getMonoWeight(Residue::Full) + mod->getDiffMonoMass());
        String accession = mod->getUniModAccession();
        accession.substitute("UniMod:", "UNIMOD:");
        fixed_proforma += "<[" + ProteoformDatabase::proformaTag(accession, mod->getDiffMonoMass()) + "]@" + String(origin) + ">";
        OPENMS_LOG_INFO << "Fixed modification " << mod->getFullId() << ": residue " << origin << " " << before << " -> "
                        << residue->getMonoWeight(Residue::Internal) << " Da" << endl;
      }
    }
    FLASHTnTFile::setGlobalFixedModifications(fixed_proforma);
    // variable modifications: checked here, so that a name or an accession the database does not know stops the run with a clear message
    for (const auto& name : getParam_().getValue("var_mod").toStringVector())
    {
      if (String(name).trim().empty()) continue;
      try
      {
        FLASHHelperClasses::resolveModification(name);
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "-var_mod " << name << ": " << e.what() << endl;
        return ILLEGAL_PARAMETERS;
      }
    }
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
    tnt_param.remove("db");
    tnt_param.remove("fixed_mod");
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
    // PEFF (or UniProtKB XML, converted to PEFF): its proteoform entries (ProteoformDatabase); FASTA: its entries are the proteins
    ProteoformDatabase database;
    const DatabaseFormat db_format = databaseFormat_(in_db);
    String peff_file = in_db;
    if (db_format == DatabaseFormat::UNIPROT_XML)
    {
      peff_file = File::getTemporaryFile();
      UniPEFFConverter::Options conversion;
      conversion.annotation_identifiers = true; // the disulfide bond pairs
      const Size converted = UniPEFFConverter::convert(in_db, peff_file, conversion);
      OPENMS_LOG_INFO << "UniProtKB XML database " << in_db << ": " << converted << " entries converted to PEFF" << endl;
    }
    const bool peff = db_format != DatabaseFormat::FASTA;
    if (peff)
    {
      database.loadPEFF(peff_file, getIntOption_("peff_max_mods"), getStringList_("peff_mod_masses"), fasta_entry);
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
      fasta_file.load(in_db, fasta_entry);
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
