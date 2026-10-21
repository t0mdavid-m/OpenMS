// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

namespace OpenMS
{
/**
@brief The protein database of FLASHTnT: proteins and their database entries.

From a FASTA file (fromFASTA) every protein is one entry of its own letters: the identity database, with which FLASHTnT searches as it
always did. From a PEFF file (PSI Extended FASTA, as UniPEFF writes it from UniProtKB XML; loadPEFF) every annotated proteoform of a
protein is an entry: the letters of a processed form (\\Processed: the full length, every mature chain, the full length without an
initiator Met, signal, transit peptide or propeptide at either end) with up to max_mods annotated modification units fixed on it
(\\ModResPsi, \\ModResUnimod; disulfide bonds from the \\DisulfideBond pairs of UniPEFF -annotation_identifiers). A unit is a residue with
a nonempty subset of its modifications that sit on distinct sites of it (side chain, alpha-amino, alpha-carboxyl group: the PSI-MOD
TermSpec; alternatives on one site, e.g. N6-acetyl and N6-succinyl, are never combined, the lipoprotein pair N-palmitoyl and
S-diacylglycerol is), or all disulfide bonds inside the form. A modification on a residue where a processed form starts is used only in
forms that start there (N-terminal, kind 2), one on a residue where a form ends (and none starts) only in forms that end there
(C-terminal, kind 3); every other one is a side-chain modification (kind 1). An annotation whose residue is not its modification's origin
(PSI-MOD Origin; ModificationsDB for Unimod) is left out.

All entries of a protein carry its identifier (sp|ACC|ID), so that FLASHTnT counts proteins and proteoforms by protein. A decoy (DECOY
prefix) of the file is the partner of the target of its identifier (the rebuilt sp|ACC|ID, or the accession). A target protein without
a decoy in the file gets one, written right after it: one permutation of its residues that keeps the first and the last residue of
every processed form (and every non-standard residue) and shuffles the others within the segments between the form boundaries
(Fisher-Yates, std::mt19937_64 seeded by FNV-1a 64 of the accession). Its entries mirror the target's: the same forms, the fixed
modifications moved with their residues.
@ingroup Topdown
*/
class OPENMS_DLLAPI ProteoformDatabase
{
public:
  /// a modification fixed on an entry
  struct FixedModification
  {
    int position = 0;  ///< residue in the entry, 0-based
    double mass = 0;   ///< monoisotopic mass shift
    String accession;  ///< MOD:nnnnn or UNIMOD:n
    String name;
    int kind = 1;      ///< 1 side chain, 2 N-terminal, 3 C-terminal (the forms it is used in)
    int sites = 1;     ///< the sites of the residue it occupies, as bits: 1 side chain, 2 alpha-amino, 4 alpha-carboxyl group
    int bond = -1;     ///< a half cystine of a disulfide bond: the residue of the other one, in the entry (0-based); else -1
  };

  /// a database entry: the letters of a processed form of a protein (FASTA entry of the same index) with fixed modifications
  struct ProteoformEntry
  {
    Size protein = 0;  ///< index into proteins
    int offset = 0;    ///< residues of the protein before the form
    int length = 0;    ///< residues of the form
    String id;         ///< identifier and form range, e.g. "sp|P69776|LPP_ECOLI 21-78"
    String proforma;   ///< the form with its fixed modifications (empty: not written in ProForma, e.g. with X)
    std::vector<FixedModification> fixed; ///< ascending by position
  };

  /// a protein: its full-length sequence and its entries
  struct DatabaseProtein
  {
    String identifier, description, sequence;
    bool decoy = false;
    long partner = -1;        ///< the target-decoy partner (index into proteins), -1: none
    std::vector<Size> entries; ///< indices into entries, ascending
    std::vector<std::pair<int, int>> forms; ///< its processed forms, 1-based first and last residue (a FASTA protein: its length)
  };

  /// counts of loadPEFF for the log
  struct LoadStatistics
  {
    Size input_targets = 0, input_decoys = 0, generated_decoys = 0;
    Size forms = 0;                       ///< processed forms of the targets
    Size modification_sites = 0;          ///< usable modification sites of the targets
    Size disulfide_bonds = 0;             ///< bonds paired by \\DisulfideBond of the targets
    Size unpaired_half_cystines = 0;      ///< half cystines without a bond partner (skipped)
    Size length_mismatches = 0;           ///< entries whose \\Length differs from their sequence
    Size empty_sequences = 0;             ///< entries without residues (left out)
    Size unpaired_decoys = 0;             ///< decoys in the file without a target (or a second decoy of one)
    Size unmirrored_decoys = 0;           ///< decoys in the file whose entries differ from their target's (count or masses)
    std::map<String, Size> unresolved;    ///< accession -> sites skipped (no mass, or skipped by mod_masses)
    std::map<String, Size> wrong_residue; ///< accession -> sites skipped (the residue is not the modification's origin)
  };

  ProteoformDatabase() = default;

  /// the identity database: one entry per FASTA entry, its own protein, no fixed modifications; the target-decoy partner by name
  /// (DECOY_ or DECOY prefix; the first entry of a name)
  static ProteoformDatabase fromFASTA(const std::vector<FASTAFile::FASTAEntry>& fasta_entries);

  /// whether @p file is a PEFF file: its first line starts with "# PEFF", or its name ends with .peff
  static bool isPEFF(const String& file);

  /**
    @brief read a PEFF file and enumerate its entries
    @param file the PEFF file (FASTA entries without PEFF keys are plain proteins of one entry)
    @param max_mods at most this many modification units per entry (0: processed forms only)
    @param mod_masses masses by accession that replace (or, "skip", leave out) those of ModificationsDB, e.g. "MOD:01819=100.016044"
    @param fasta_entries output: per entry, the protein's identifier and description and the form's letters
  */
  void loadPEFF(const String& file, int max_mods, const std::vector<String>& mod_masses, std::vector<FASTAFile::FASTAEntry>& fasta_entries);

  /// the entries of a parsed PEFF (or FASTA) entry list, as loadPEFF makes them from its file; for tests and callers that read the file themselves
  void build(const std::vector<FASTAFile::FASTAEntry>& input, int max_mods, const std::vector<String>& mod_masses,
             std::vector<FASTAFile::FASTAEntry>& fasta_entries);

  /// one entry per protein, at offset 0 over its whole length, without fixed modifications (a FASTA search)
  bool isIdentity() const
  {
    return identity_;
  }

  /// whether any entry has a fixed modification
  bool hasFixedModifications() const
  {
    return with_fixed_;
  }

  /// per entry: (position in the entry, summed mass of its fixed modifications), ascending by position
  const std::vector<std::vector<std::pair<int, double>>>& fixedMasses() const
  {
    return fixed_masses_;
  }

  /// per entry: (position in the entry, what its fixed modifications occupy as bits: 1 side chain, 2 N-terminus, 4 C-terminus), ascending.
  /// A kind-2 (3) modification occupies the N- (C-) terminus, and also the side chain unless it is an N- (C-) terminal one by its
  /// PSI-MOD TermSpec (e.g. S-diacylglycerol-Cys at a mature start, N6-succinyl-Lys at a form end); a kind-1 one the side chain.
  const std::vector<std::vector<std::pair<int, std::uint8_t>>>& occupancy() const
  {
    return occupancy_;
  }

  /// the protein of each entry
  const std::vector<Size>& proteinOfEntry() const
  {
    return protein_of_;
  }

  const LoadStatistics& getLoadStatistics() const
  {
    return statistics_;
  }

  /// FNV-1a 64 of the bytes of @p s (the decoy seed)
  static std::uint64_t fnv1a64(const String& s);

  /// the ProForma tag of a fixed modification: its mass with its accession as an alternative, e.g. "+42.010565|MOD:00064". The mass comes
  /// first: PSI-MOD's own DiffMono can differ from the mass searched (MOD:00083 trimethyl-K: 43.054775 for the charged form, 42.04695 for a
  /// neutral mass; or none in the bundled PSI-MOD), so a reader outside OpenMS gets the mass that was searched.
  static String proformaTag(const String& accession, double mass);

  std::vector<ProteoformEntry> entries;
  std::vector<DatabaseProtein> proteins;

private:
  /// fill the per-entry tables and the identity flag from entries and proteins
  void finish_();

  bool identity_ = true;
  bool with_fixed_ = false;
  std::vector<std::vector<std::pair<int, double>>> fixed_masses_;
  std::vector<std::vector<std::pair<int, std::uint8_t>>> occupancy_;
  std::vector<Size> protein_of_;
  LoadStatistics statistics_;
};
} // namespace OpenMS
