// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

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
      std::string ptmlist;                             ///< UniProt ptmlist.txt; empty: the bundled snapshot
      std::string psimod_obo;                          ///< PSI-MOD.obo for canonical names; empty: the bundled one
      std::string unimod_obo;                          ///< optional unimod.obo for canonical Unimod names
      std::string prefix;                              ///< one PEFF prefix for every entry; empty: sp/tr from the UniProt dataset
      std::string dbversion = "unknown";               ///< value of the '# DbVersion=' header line
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
    static Size convert(const std::string& in_xml, const std::string& out_peff, const Options& options);
  };
} // namespace OpenMS
