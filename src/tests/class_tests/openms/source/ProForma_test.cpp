// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ProForma.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ProForma, "$Id$")

START_SECTION(static ProForma fromString(const String& proforma))
{
  // round trips of the supported notation
  for (const String s : {"PEPTIDE", "PEM[Oxidation]K", "S[Phospho][+1.0]EK", "P(EMK)[+15.9949]R", "[Acetyl]-PEPTIDE", "PEPTIDE-[Amidated]",
                         "[Phospho]?PEPTIDE", "{Glycan:Hex}PEPTIDE", "<[Carbamidomethyl]@C>PEPTCDE", "PEPTIDE/2", "EM[MOD:00719]EVEES[MOD:00046]PEK",
                         "EM[UNIMOD:35]K", "M[Formula:O]K", "(PEP)[-18.0106]TIDE"})
  {
    TEST_EQUAL(ProForma::fromString(s).toString(), s)
  }

  // masses: named, accession, PSI-MOD, formula, glycan, alternatives, observed
  TEST_REAL_SIMILAR(ProForma::fromString("PEM[Oxidation]K").getModifications()[0].mass, 15.994915)
  TEST_REAL_SIMILAR(ProForma::fromString("PEM[U:Oxidation]K").getModifications()[0].mass, 15.994915)
  TEST_REAL_SIMILAR(ProForma::fromString("PEM[UNIMOD:35]K").getModifications()[0].mass, 15.994915)
  TEST_REAL_SIMILAR(ProForma::fromString("PEM[MOD:00719]K").getModifications()[0].mass, 15.994915)
  TEST_REAL_SIMILAR(ProForma::fromString("PEM[Formula:O]K").getModifications()[0].mass, 15.994915)
  TEST_REAL_SIMILAR(ProForma::fromString("PEN[Glycan:Hex1HexNAc1]K").getModifications()[0].mass, 365.132197)
  TEST_REAL_SIMILAR(ProForma::fromString("PES[INFO:x|Phospho|+79.97]K").getModifications()[0].mass, 79.966331)
  TEST_REAL_SIMILAR(ProForma::fromString("PES[Obs:+79.97]K").getModifications()[0].mass, 79.97)
  TEST_REAL_SIMILAR(ProForma::fromString("PEK[Cation:Ca[II]]").getModifications()[0].mass, 37.946941)

  // sites and positions
  TEST_EXCEPTION(Exception::ParseError, ProForma::fromString("PEPTIDE#"))
  ProForma pf = ProForma::fromString("<[Carbamidomethyl]@C>[Acetyl]-PC(EMK)[+15.9949]S[Phospho]C-[Amidated]/3");
  TEST_EQUAL(pf.getSequence(), "PCEMKSC")
  TEST_EQUAL(pf.getCharge(), 3)
  TEST_EQUAL(pf.getModifications().size(), 6)
  int globals = 0, ranges = 0, n_term = 0, c_term = 0;
  for (const auto& m : pf.getModifications())
  {
    globals += m.global;
    if (m.site == ProForma::Site::RESIDUE && m.end > m.start)
    {
      ranges++;
      TEST_EQUAL(m.start, 2)
      TEST_EQUAL(m.end, 4)
    }
    n_term += m.site == ProForma::Site::N_TERM;
    c_term += m.site == ProForma::Site::C_TERM;
  }
  TEST_EQUAL(globals, 2)
  TEST_EQUAL(ranges, 1)
  TEST_EQUAL(n_term, 1)
  TEST_EQUAL(c_term, 1)
  TEST_EQUAL(pf.toString(), "<[Carbamidomethyl]@C>[Acetyl]-PC(EMK)[+15.9949]S[Phospho]C-[Amidated]/3")

  // unknown position with a count, a localisation group (the tag with the name carries the mass)
  pf = ProForma::fromString("[Phospho]^2?PEPTIDE");
  TEST_EQUAL(pf.getModifications().size(), 2)
  pf = ProForma::fromString("PES[Phospho#g1(0.8)]T[#g1(0.2)]K");
  TEST_EQUAL(pf.getModifications().size(), 1)
  TEST_EQUAL(pf.getModifications()[0].start, 2)

  // invalid or unsupported
  TEST_EXCEPTION(Exception::ParseError, ProForma::fromString("A[PhosphoCDEFGHIK"))
  TEST_EXCEPTION(Exception::ParseError, ProForma::fromString("A[NoSuchModification]C"))
  TEST_EXCEPTION(Exception::ParseError, ProForma::fromString("PEPTIDE+PEPTIDE"))
  TEST_EXCEPTION(Exception::ParseError, ProForma::fromString("PEPTBDE"))
  TEST_EXCEPTION(Exception::ParseError, ProForma::fromString("PEP(TI(D)[+1])[+1]E"))
  TEST_EXCEPTION(Exception::ParseError, ProForma::fromString("PEPTIDE/[Na:z+1]"))
  TEST_EXCEPTION(Exception::ParseError, ProForma::fromString("<13C>PEPTIDE"))
  TEST_EXCEPTION(Exception::ParseError, ProForma::fromString("PEK[#XL1]"))
}
END_SECTION

START_SECTION(String toString() const)
{
  ProForma pf("PEPTIDEK");
  pf.addModification(1, 1, 15.994915);
  pf.addModification(3, 5, -18.010565, "");
  pf.addModification(3, 5, 1.0);
  pf.addModification(ProForma::Site::N_TERM, 42.010565, "Acetyl");
  TEST_EQUAL(pf.toString(), "[Acetyl]-PE[+15.994915]P(TID)[-18.010565][+1.0]EK")
  TEST_EXCEPTION(Exception::InvalidValue, pf.addModification(4, 6, 1.0))
  TEST_EXCEPTION(Exception::InvalidValue, pf.addModification(7, 8, 1.0))
  pf.addModification(4, 4, 2.0); // a residue inside a range
  TEST_EQUAL(pf.toString(), "[Acetyl]-PE[+15.994915]P(TI[+2.0]D)[-18.010565][+1.0]EK")
}
END_SECTION

START_SECTION(double getMonoMass() const)
{
  const ProForma pf = ProForma::fromString("{Glycan:Hex}[Acetyl]-PEM[Oxidation]K");
  const double expected = AASequence::fromString("PEMK").getMonoWeight() + 162.052824 + 42.010565 + 15.994915;
  TEST_REAL_SIMILAR(pf.getMonoMass(), expected)
}
END_SECTION

START_SECTION(static String formatMass(double mass))
{
  TEST_EQUAL(ProForma::formatMass(15.994915), "+15.994915")
  TEST_EQUAL(ProForma::formatMass(-18.010565), "-18.010565")
  TEST_EQUAL(ProForma::formatMass(80.9787), "+80.9787")
  TEST_EQUAL(ProForma::formatMass(1.0), "+1.0")
}
END_SECTION

END_TEST
