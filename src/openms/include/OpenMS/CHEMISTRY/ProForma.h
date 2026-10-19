// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <vector>

namespace OpenMS
{
/**
  @brief A proteoform in ProForma 2.0 notation (HUPO-PSI, https://github.com/HUPO-PSI/ProForma): the part used to report and to
  read top-down proteoform-spectrum matches.

  One peptidoform whose modifications are given by name, accession, formula, glycan composition or mass, and sit on a residue,
  on a range of residues, on a terminus, at an unknown position, or are labile.

  Read and written:
  - residues: the 20 standard amino acids, U and O
  - modification tags: a name looked up in Unimod, then PSI-MOD ([Oxidation], [L-methionine sulfoxide]); a name with its source
    ([U:Oxidation], [M:L-methionine sulfoxide]); an accession ([UNIMOD:35], [MOD:00719]); a formula ([Formula:C2H2O]); a glycan
    composition ([Glycan:Hex1HexNAc2]); a mass ([+15.9949], [Obs:+15.995]). Alternatives are separated by | and the first one that
    gives a mass is used ([Phospho|+79.966]); INFO: parts are ignored. A localisation group or score suffix (#g1, (0.9)) is dropped:
    the modification stays where its tag is written.
  - on a residue: PEM[Oxidation]K, several on one residue: S[Phospho][+1.0]; on a range: P(EMK)[+15.995]; on a terminus:
    [Acetyl]-PEPTIDE, PEPTIDE-[Amidated]; at an unknown position: [Phospho]?PEPTIDE, [Phospho]^2?PEPTIDE; labile:
    {Glycan:Hex}PEPTIDE; global fixed: <[Carbamidomethyl]@C>PEPTCDE; charge: PEPTIDE/2
  - written in this order: global, unknown position, labile, N-terminal, sequence, C-terminal, charge (the order of ProForma 2.0);
    read in any order before the sequence

  Not supported (Exception::ParseError): chimeras (+), cross-links and branches (#XL, #BRANCH), isotope labels (<13C>), adducts as
  charge (/[Na:z+1]), the ambiguous residues B, Z, J and X, and overlapping ranges.
*/
class OPENMS_DLLAPI ProForma
{
public:
  /// where a modification sits
  enum class Site
  {
    RESIDUE, ///< on residues start..end (a range if end > start)
    N_TERM,
    C_TERM,
    UNKNOWN, ///< at an unknown position of the sequence
    LABILE   ///< lost on fragmentation: in the precursor mass only
  };

  /// a modification
  struct Modification
  {
    Site site = Site::RESIDUE;
    int start = 0;       ///< first residue, 0-based (Site::RESIDUE)
    int end = 0;         ///< last residue, 0-based (Site::RESIDUE)
    double mass = 0;     ///< monoisotopic mass shift in Da
    String tag;          ///< the tag content as written (name, accession, formula, ...); empty: the mass is written
    bool global = false; ///< from a global fixed modification <[tag]@residues> (written there, not on the residue)
  };

  ProForma() = default;

  /// a sequence of one-letter residues without modifications; Exception::ParseError on residues other than the standard ones, U and O
  explicit ProForma(const String& sequence);

  /// parse a ProForma string; Exception::ParseError if it is invalid or uses an unsupported feature
  static ProForma fromString(const String& proforma);

  /// the ProForma string
  String toString() const;

  /**
    @brief add a modification on residues @p start .. @p end (0-based; a range if end > start)

    With an empty @p tag the modification is written as its mass. Exception::InvalidValue if the residues are outside the sequence or
    the range overlaps another range (ProForma has no overlapping ranges); several modifications on one residue are allowed.
  */
  void addModification(int start, int end, double mass, const String& tag = "");

  /// add a modification at a site other than residues (N_TERM, C_TERM, UNKNOWN, LABILE)
  void addModification(Site site, double mass, const String& tag = "");

  /// the sequence (one-letter residues, no modifications)
  const String& getSequence() const
  {
    return sequence_;
  }

  /// all modifications, global fixed ones expanded onto their residues
  const std::vector<Modification>& getModifications() const
  {
    return mods_;
  }

  /// the charge (/z), 0 if none
  int getCharge() const
  {
    return charge_;
  }

  void setCharge(int charge)
  {
    charge_ = charge;
  }

  /// monoisotopic neutral mass: residues, water and all modifications (labile and unknown-position ones included)
  double getMonoMass() const;

  /// the mass of a tag content (see the class description); @p residue (one-letter code or 0) narrows a name lookup.
  /// Exception::ParseError if no alternative gives a mass.
  static double getTagMass(const String& tag, char residue = 0);

  /// a mass as written in a tag: sign, at most six decimals, no trailing zeros (+15.994915, -18.010565, +1.0); six decimals keep a written
  /// and read-back PrSM identical for annotation (four moved the mass calibration of FLASHAnnotator in 32 of 4245 PrSMs)
  static String formatMass(double mass);

private:
  String sequence_;
  std::vector<Modification> mods_;
  int charge_ = 0;
};
} // namespace OpenMS
