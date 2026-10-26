// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
    @brief Converts a UniProtKB XML protein database (plain or gzip) to PEFF 1.0 with modification, processing, variant and
    disulfide-bond annotations: the conversion of the UniPEFF tool (an OpenMS port of David L. Tabb's C# UniPEFF), as a library
    function so that other tools (FLASHTnT) can read UniProtKB XML directly.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI UniPEFFConverter
  {
  public:
    /// the options of the UniPEFF tool
    struct Options
    {
      String ptmlist;                             ///< UniProt ptmlist.txt; empty: the bundled snapshot
      String psimod_obo;                          ///< PSI-MOD.obo for canonical names; empty: the bundled one
      String unimod_obo;                          ///< optional unimod.obo for canonical Unimod names
      String prefix;                              ///< one PEFF prefix for every entry; empty: sp/tr from the UniProt dataset
      String dbversion = "unknown";               ///< value of the '# DbVersion=' header line
      bool annotation_identifiers = false;        ///< PEFF Option B: annotation ids and \DisulfideBond connectivity
      bool omit_molecular_processing = false;     ///< no \Processed annotations
      bool omit_amino_acid_modifications = false; ///< no \ModResPsi / \ModResUnimod / \ModRes / \DisulfideBond
      bool omit_sequence_variations = false;      ///< no \VariantSimple / \VariantComplex
    };

    /**
      @brief Writes the PEFF of @p in_xml to @p out_peff (a plain, uncompressed file name).

      @return the number of entries written
      @exception Exception::FileNotFound the bundled ptmlist cannot be found
      @exception Exception::UnableToCreateFile the output or a temporary spool file cannot be written, or @p out_peff is compressed
      @exception Exception::ParseError the input holds no writable entry
    */
    static Size convert(const String& in_xml, const String& out_peff, const Options& options);
  };
} // namespace OpenMS
