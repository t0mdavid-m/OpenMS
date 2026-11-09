// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

// FixedModifications (ADR-0067): the one residue-mass table of a search. Pure: nothing here writes ResidueDB.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FixedModifications.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Exception.h>

using namespace OpenMS;

START_TEST(FixedModifications, "$Id$")

// An empty table is ResidueDB, read live: the same doubles every reader read before.
START_SECTION(none_is_residuedb)
{
  const ResidueDB* rdb = ResidueDB::getInstance();
  const Residue* unknown = rdb->getResidue(String("X"));
  for (int c = 'A'; c <= 'Z'; ++c)
  {
    if (c == 'X') continue;
    const Residue* r = rdb->getResidue((unsigned char)c);
    if (r == nullptr || r == unknown) continue;
    TEST_EQUAL(FixedModifications::none().mass((char)c) == r->getMonoWeight(Residue::Internal), true) // bitwise
    TEST_EQUAL(FixedModifications::none().delta((char)c), 0.0)
  }
  TEST_EQUAL(FixedModifications::none().empty(), true)
  TEST_STRING_EQUAL(FixedModifications::none().proformaPrefix(), "")
}
END_SECTION

START_SECTION(fromNames_resolves_a_name_and_an_accession)
{
  for (const std::string s : {std::string("Carbamidomethyl (C)"), std::string("4 (C)"), std::string("UNIMOD:4 (C)")})
  {
    const FixedModifications t = FixedModifications::fromNames({s});
    TEST_EQUAL(t.empty(), false)
    TEST_REAL_SIMILAR(t.delta('C'), 57.021464)
    TEST_EQUAL(t.mass('A') == FixedModifications::none().mass('A'), true) // other residues untouched
    TEST_EQUAL(t.entries().size(), 1)
  }
  // empty strings are skipped (the Tagger's and Extender's fixed_mod default is {""})
  TEST_EQUAL(FixedModifications::fromNames({"", "Carbamidomethyl (C)"}) == FixedModifications::fromNames({"Carbamidomethyl (C)"}), true)
  TEST_EQUAL(FixedModifications::fromNames({""}).empty(), true)
  TEST_EQUAL(FixedModifications::fromNames({}).empty(), true)
}
END_SECTION

// mass() of a fixed residue is the mass the FLASHTnT tool's ResidueDB overwrite gave it, to the bit.
START_SECTION(mass_is_the_tool_overwrite_value)
{
  const FixedModifications t = FixedModifications::fromNames({"Carbamidomethyl (C)"});
  Residue copy(*ResidueDB::getInstance()->getResidue((unsigned char)'C'));
  copy.setMonoWeight(copy.getMonoWeight(Residue::Full) + t.delta('C'));
  TEST_EQUAL(t.mass('C') == copy.getMonoWeight(Residue::Internal), true)
  TEST_REAL_SIMILAR(t.mass('C'), ResidueDB::getInstance()->getResidue((unsigned char)'C')->getMonoWeight(Residue::Internal) + 57.021464)
  // ResidueDB itself is untouched
  TEST_EQUAL(FixedModifications::none().mass('C') == t.mass('C'), false)
}
END_SECTION

START_SECTION(fromNames_refuses_what_the_tool_refused)
{
  TEST_EXCEPTION(Exception::InvalidValue, FixedModifications::fromNames({"Acetyl (Protein N-term)"}))              // terminal
  TEST_EXCEPTION(Exception::InvalidValue, FixedModifications::fromNames({"Carbamidomethyl (C)", "iodoTMT6plex (C)"})) // two on C
  TEST_EXCEPTION(Exception::InvalidValue, FixedModifications::fromNames({"NoSuchModification (C)"}))              // unknown
}
END_SECTION

START_SECTION(internalMass_adds_the_deltas)
{
  const FixedModifications t = FixedModifications::fromNames({"iodoTMT6plex (C)"});
  const AASequence s = AASequence::fromString("PEPCTCDE");
  TEST_REAL_SIMILAR(t.internalMass(s), s.getMonoWeight(Residue::Internal) + 2 * 329.226595)
  TEST_EQUAL(FixedModifications::none().internalMass(s) == s.getMonoWeight(Residue::Internal), true) // bitwise when empty
}
END_SECTION

START_SECTION(the_proforma_prefix_round_trips)
{
  const FixedModifications t = FixedModifications::fromNames({"Carbamidomethyl (C)", "Oxidation (M)"});
  TEST_EQUAL(t.proformaPrefix().hasPrefix("<["), true)
  TEST_EQUAL(FixedModifications::fromProFormaPrefix(t.proformaPrefix() + "PEPTCMDE") == t, true)
  TEST_EQUAL(FixedModifications::fromProFormaPrefix("PEPTCDE").empty(), true)
  TEST_EXCEPTION(Exception::InvalidValue, FixedModifications::fromProFormaPrefix("<[+57.021464]@C>PEPTCDE")) // no accession
}
END_SECTION

END_TEST
