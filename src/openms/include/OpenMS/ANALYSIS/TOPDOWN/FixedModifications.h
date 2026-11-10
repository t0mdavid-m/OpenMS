// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
   * @brief The fixed modifications of a search: one residue-mass table that every reader of a residue mass shares (ADR-0067).
   *
   * A fixed modification sits on every occurrence of one residue, anywhere in the sequence, for the whole search: it changes what the
   * residue weighs, so a tag reads the modified residue as its letter, every ladder, extension and annotation adds the modified mass, and a
   * PrSM never lists it among its modifications. It replaces the FLASHTnT tool's former process-wide ResidueDB overwrite.
   *
   * mass() of a residue without a fixed modification is read live from ResidueDB, as every caller read it before, so an empty table
   * gives the same doubles. A fixed residue's mass is computed as the overwrite did: a copy of the Residue with
   * setMonoWeight(full weight + modification mass).
   */
  class OPENMS_DLLAPI FixedModifications
  {
  public:
    /// one fixed modification
    struct Entry
    {
      char residue = 0;  ///< the one-letter code it sits on
      double mass = 0;   ///< its monoisotopic mass difference
      String accession;  ///< "UNIMOD:4" (empty if it has none)
      String name;       ///< as given, e.g. "Carbamidomethyl (C)"
    };

    /// none
    FixedModifications() = default;

    /// an empty table (shared; nothing cached in it)
    static const FixedModifications& none();

    /**
     * @brief The table of @p names: each a modification name or Unimod accession with its residue ("Carbamidomethyl (C)", "4 (C)",
     * "UNIMOD:4 (C)"); empty strings are skipped.
     * @throw Exception::InvalidValue an unknown name, a terminal modification (only a residue anywhere can be fixed), or two fixed
     * modifications of one residue
     */
    static FixedModifications fromNames(const std::vector<std::string>& names);

    /**
     * @brief The table of the global fixed modifications that prefix a ProForma string, e.g. "<[+57.021464|UNIMOD:4]@C>PEPTCDE"; each
     * group is resolved by its accession ("UNIMOD:4 (C)"). A string without a prefix gives an empty table.
     * @throw Exception::InvalidValue a group without a Unimod accession, or one fromNames refuses
     */
    static FixedModifications fromProFormaPrefix(const String& proforma);

    /// no fixed modification
    bool empty() const;

    /// the internal monoisotopic mass of residue @p c with its fixed modification (ResidueDB's, read live, if it has none)
    double mass(char c) const;

    /// the fixed modification mass on residue @p c (0 if none)
    double delta(char c) const;

    /// @p seq's internal monoisotopic mass plus the fixed modification masses of its residues
    double internalMass(const AASequence& seq) const;

    /// the fixed modifications, in the order given
    const std::vector<Entry>& entries() const;

    /// their names as given, for a fixed_mod parameter
    std::vector<std::string> names() const;

    /// their ProForma global modifications, e.g. "<[+57.021464|UNIMOD:4]@C>" (empty if none)
    String proformaPrefix() const;

    /// the same residues with the same masses
    bool operator==(const FixedModifications& other) const;
    bool operator!=(const FixedModifications& other) const;

  private:
    /// the residue tables by one-letter code
    struct Tables
    {
      std::array<bool, 256> fixed {};
      std::array<double, 256> delta {};
      std::array<double, 256> modified {};
    };
    /// Shared by copies and immutable once built; null when nothing is fixed. Held by pointer so the type stays small: every
    /// Config and every algorithm instance holds one, and a test's main() holds many Configs on the stack.
    std::shared_ptr<const Tables> tables_;
    std::vector<Entry> entries_;
  };
} // namespace OpenMS
