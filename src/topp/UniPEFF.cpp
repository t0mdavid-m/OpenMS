// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/UniPEFFConverter.h>

using namespace OpenMS;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_UniPEFF UniPEFF

@brief Converts a UniProtKB XML protein database to PEFF 1.0 with per-entry
modification, processing, variant and disulfide-bond annotations.

This is an OpenMS-native port of the standalone C# UniPEFF tool
(David L. Tabb, UMC Groningen). For each entry it emits a PEFF descriptor
line with \\PName \\GName \\NcbiTaxId \\TaxName \\Length \\SV \\EV \\PE \\ID
\\AltAC \\ModResPsi \\ModResUnimod \\ModRes \\VariantSimple \\VariantComplex
\\Processed, plus \\DisulfideBond when @c -annotation_identifiers is set.

UniProt's <code>\<feature type="disulfide bond"\></code> entries are translated
into half-cystine modifications (PSI-MOD:00798) which are then merged into
the main modified-residue list in position-sorted order; in Option B
(@c -AnnotationIdentifiers) the resulting half-cystines are referenced by
annotation id from <code>\\DisulfideBond=(id1,id2)</code>.

Modification accession lookup uses UniProt's <em>ptmlist.txt</em> (a snapshot
is bundled under @c share/OpenMS/CHEMISTRY/UniProt_ptmlist.txt); override
with @c -ptmlist. Canonical OBO names come from <em>PSI-MOD.obo</em> (bundled)
and an optional <em>unimod.obo</em>; without them, names fall back to the
UniProt @c ptmlist.txt ID and a warning is printed.

Both plain @c .xml and @c .xml.gz UniProt inputs are accepted (gzip is
auto-detected by the underlying parser).

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_UniPEFF.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_UniPEFF.html
*/

/// @cond TOPPCLASSES

// ──────────────────────────────────────────────────────────────────
// TOPP wrapper
// ──────────────────────────────────────────────────────────────────

class TOPPUniPEFF :
  public TOPPBase
{
public:
  TOPPUniPEFF() :
    TOPPBase("UniPEFF", "Convert a UniProtKB XML protein database to PEFF 1.0 with rich annotations.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input UniProtKB XML file (plain or gzip; gzip is auto-detected).");
    setValidFormats_("in", {"xml"});  // .xml.gz is accepted: FileHandler strips .gz and re-checks

    registerOutputFile_("out", "<file>", "", "Output PEFF 1.0 file.");
    // no setValidFormats_: this branch's FileTypes has no PEFF type

    registerInputFile_("ptmlist", "<file>", "", "UniProt ptmlist.txt; defaults to the bundled snapshot.", false);
    setValidFormats_("ptmlist", {"txt"});

    registerInputFile_("psimod_obo", "<file>", "", "PSI-MOD.obo for canonical modification names; defaults to the bundled OpenMS PSI-MOD.obo.", false);
    setValidFormats_("psimod_obo", {"obo"});

    registerInputFile_("unimod_obo", "<file>", "", "Optional unimod.obo for canonical Unimod names; if absent, names fall back to the UniProt ptmlist IDs.", false);
    setValidFormats_("unimod_obo", {"obo"});

    registerStringOption_("prefix", "<string>", "", "Force a single PEFF prefix for every entry (e.g. 'sp'); if empty, sp/tr is derived from the UniProt dataset.", false);
    registerStringOption_("dbversion", "<string>", "unknown", "Value for the mandatory '# DbVersion=' PEFF header line.", false);

    registerFlag_("annotation_identifiers", "Emit PEFF Option B: assign a sequential id: prefix to every annotation tuple and emit \\DisulfideBond connectivity.");
    registerFlag_("omit_molecular_processing", "Skip the \\Processed annotations (initiator methionine, signal/transit peptide, propeptide, chain).");
    registerFlag_("omit_amino_acid_modifications", "Skip \\ModResPsi / \\ModResUnimod / \\ModRes and \\DisulfideBond; ptmlist is not read.");
    registerFlag_("omit_sequence_variations", "Skip \\VariantSimple and \\VariantComplex annotations.");
  }

  ExitCodes main_(int, const char**) override
  {
    UniPEFFConverter::Options options;
    options.ptmlist = getStringOption_("ptmlist");
    options.psimod_obo = getStringOption_("psimod_obo");
    options.unimod_obo = getStringOption_("unimod_obo");
    options.prefix = getStringOption_("prefix");
    options.dbversion = getStringOption_("dbversion");
    options.annotation_identifiers = getFlag_("annotation_identifiers");
    options.omit_molecular_processing = getFlag_("omit_molecular_processing");
    options.omit_amino_acid_modifications = getFlag_("omit_amino_acid_modifications");
    options.omit_sequence_variations = getFlag_("omit_sequence_variations");
    // the conversion is the library's UniPEFFConverter (FLASHTnT reads UniProtKB XML with it, too)
    try
    {
      UniPEFFConverter::convert(getStringOption_("in"), getStringOption_("out"), options);
    }
    catch (const Exception::FileNotFound&)
    {
      writeLogError_("UniPEFF: cannot locate bundled CHEMISTRY/UniProt_ptmlist.txt and -ptmlist not given.");
      return INPUT_FILE_NOT_FOUND;
    }
    catch (const Exception::UnableToCreateFile& e)
    {
      writeLogError_(e.what());
      return CANNOT_WRITE_OUTPUT_FILE;
    }
    catch (const Exception::ParseError& e)
    {
      writeLogError_(e.what());
      return PARSE_ERROR;
    }
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPUniPEFF tool;
  return tool.main(argc, argv);
}

/// @endcond
