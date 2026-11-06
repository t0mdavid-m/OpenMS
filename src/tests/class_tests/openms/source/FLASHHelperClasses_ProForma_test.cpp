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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProForma.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

using PFM = FLASHHelperClasses::ProFormaModification;

namespace
{
  /// the modifications of a ProForma string as FLASHAnnotator reads them
  std::vector<PFM> read(const std::string& proforma, std::string& sequence)
  {
    return FLASHHelperClasses::getProFormaModifications(ProForma::parseIon(proforma), sequence);
  }

  std::vector<PFM> read(const std::string& proforma)
  {
    std::string sequence;
    return read(proforma, sequence);
  }

  /// read and written again
  std::string roundTrip(const std::string& proforma)
  {
    std::string sequence;
    const auto mods = read(proforma, sequence);
    return ProForma::toString(FLASHHelperClasses::toPeptidoform(sequence, mods));
  }

  PFM residueMod(int start, int end, double mass, const std::string& tag = "")
  {
    PFM m;
    m.start = start;
    m.end = end;
    m.mass = mass;
    m.tag = tag;
    return m;
  }
}

START_TEST(FLASHHelperClasses_ProForma, "$Id$")

START_SECTION(static std::vector<ProFormaModification> getProFormaModifications(const ProForma::PeptidoformIon& ion, std::string& sequence))
{
  // round trips of the supported notation (a glycan is written with its counts, the charge is not written)
  for (const std::string s : {"PEPTIDE", "PEM[Oxidation]K", "S[Phospho][+1.0]EK", "P(EMK)[+15.9949]R", "[Acetyl]-PEPTIDE", "PEPTIDE-[Amidated]",
                         "[Phospho]?PEPTIDE", "{Glycan:Hex1}PEPTIDE", "<[Carbamidomethyl]@C>PEPTCDE", "EM[MOD:00719]EVEES[MOD:00046]PEK",
                         "EM[UNIMOD:35]K", "M[Formula:O]K", "(PEP)[-18.0106]TIDE"})
  {
    TEST_EQUAL(roundTrip(s), s)
  }
  TEST_EQUAL(roundTrip("PEPTIDE/2"), "PEPTIDE")
  TEST_EQUAL(std::get<int>(ProForma::parseIon("PEPTIDE/2").charge.value()), 2)

  // masses: named, accession, PSI-MOD, formula, glycan, alternatives, observed
  TEST_REAL_SIMILAR(read("PEM[Oxidation]K")[0].mass, 15.994915)
  TEST_REAL_SIMILAR(read("PEM[U:Oxidation]K")[0].mass, 15.994915)
  TEST_REAL_SIMILAR(read("PEM[UNIMOD:35]K")[0].mass, 15.994915)
  TEST_REAL_SIMILAR(read("PEM[MOD:00719]K")[0].mass, 15.994915)
  TEST_REAL_SIMILAR(read("PEM[Formula:O]K")[0].mass, 15.994915)
  TEST_REAL_SIMILAR(read("PEN[Glycan:Hex1HexNAc1]K")[0].mass, 365.132197)
  TEST_REAL_SIMILAR(read("PES[INFO:x|Phospho|+79.97]K")[0].mass, 79.966331)
  TEST_REAL_SIMILAR(read("PES[Obs:+79.97]K")[0].mass, 79.97)
  TEST_REAL_SIMILAR(read("PEK[Cation:Ca[II]]")[0].mass, 37.946941)

  // sites and positions
  TEST_EXCEPTION(Exception::ParseError, read("PEPTIDE#"))
  const std::string full = "<[Carbamidomethyl]@C>[Acetyl]-PC(EMK)[+15.9949]S[Phospho]C-[Amidated]/3";
  std::string sequence;
  const auto mods = read(full, sequence);
  TEST_EQUAL(sequence, "PCEMKSC")
  TEST_EQUAL(std::get<int>(ProForma::parseIon(full).charge.value()), 3)
  TEST_EQUAL(mods.size(), 6)
  int globals = 0, ranges = 0, n_term = 0, c_term = 0;
  for (const auto& m : mods)
  {
    globals += m.global;
    if (m.site == PFM::Site::RESIDUE && m.end > m.start)
    {
      ranges++;
      TEST_EQUAL(m.start, 2)
      TEST_EQUAL(m.end, 4)
    }
    n_term += m.site == PFM::Site::N_TERM;
    c_term += m.site == PFM::Site::C_TERM;
  }
  TEST_EQUAL(globals, 2)
  TEST_EQUAL(ranges, 1)
  TEST_EQUAL(n_term, 1)
  TEST_EQUAL(c_term, 1)
  TEST_EQUAL(ProForma::toString(FLASHHelperClasses::toPeptidoform(sequence, mods)), "<[Carbamidomethyl]@C>[Acetyl]-PC(EMK)[+15.9949]S[Phospho]C-[Amidated]")

  // unknown position with a count, a localisation group (the tag with the name carries the mass)
  TEST_EQUAL(read("[Phospho]^2?PEPTIDE").size(), 2)
  const auto grouped = read("PES[Phospho#g1(0.8)]T[#g1(0.2)]K");
  TEST_EQUAL(grouped.size(), 1)
  TEST_EQUAL(grouped[0].start, 2)

  // invalid or unsupported
  TEST_EXCEPTION(Exception::ParseError, read("A[PhosphoCDEFGHIK"))
  TEST_EXCEPTION(Exception::ParseError, read("A[NoSuchModification]C"))
  TEST_EXCEPTION(Exception::ParseError, read("PEPTIDE+PEPTIDE"))
  TEST_EXCEPTION(Exception::ParseError, read("PEPTBDE"))
  TEST_EXCEPTION(Exception::ParseError, read("PEP(TI(D)[+1])[+1]E"))
  TEST_EXCEPTION(Exception::ParseError, read("PEPTIDE/[Na:z+1]"))
  TEST_EXCEPTION(Exception::ParseError, read("<13C>PEPTIDE"))
  TEST_EXCEPTION(Exception::ParseError, read("PEK[#XL1]"))
}
END_SECTION

START_SECTION(static ProForma::PeptidoformIon parseProForma(const std::string& text))
{
  // FLASHTnT's fixed modification without an accession: the "-" alternative is dropped, the mass stays on its residue
  std::string sequence;
  const auto mods = FLASHHelperClasses::getProFormaModifications(FLASHHelperClasses::parseProForma("PEPC[+57.021464|-]TIDEK[+42.010565|-|INFO:x]"), sequence);
  TEST_EQUAL(sequence, "PEPCTIDEK")
  TEST_EQUAL(mods.size(), 2)
  TEST_REAL_SIMILAR(mods[0].mass, 57.021464)
  TEST_EQUAL(mods[0].start, 3)
  TEST_REAL_SIMILAR(mods[1].mass, 42.010565)
  TEST_EQUAL(mods[1].start, 8)
  TEST_EXCEPTION(Exception::ParseError, ProForma::parseIon("PEPC[+57.021464|-]TIDE"))
}
END_SECTION

START_SECTION(static ProForma::Peptidoform toPeptidoform(const std::string& sequence, const std::vector<ProFormaModification>& modifications))
{
  std::vector<PFM> mods = {residueMod(1, 1, 15.994915), residueMod(3, 5, -18.010565, ""), residueMod(3, 5, 1.0)};
  PFM acetyl;
  acetyl.site = PFM::Site::N_TERM;
  acetyl.mass = 42.010565;
  acetyl.tag = "Acetyl";
  mods.push_back(acetyl);
  TEST_EQUAL(ProForma::toString(FLASHHelperClasses::toPeptidoform("PEPTIDEK", mods)), "[Acetyl]-PE[+15.994915]P(TID)[-18.010565][+1.0]EK")
  auto overlapping = mods;
  overlapping.push_back(residueMod(4, 6, 1.0));
  TEST_EXCEPTION(Exception::InvalidValue, FLASHHelperClasses::toPeptidoform("PEPTIDEK", overlapping))
  auto outside = mods;
  outside.push_back(residueMod(7, 8, 1.0));
  TEST_EXCEPTION(Exception::InvalidValue, FLASHHelperClasses::toPeptidoform("PEPTIDEK", outside))
  mods.push_back(residueMod(4, 4, 2.0)); // a residue inside a range
  TEST_EQUAL(ProForma::toString(FLASHHelperClasses::toPeptidoform("PEPTIDEK", mods)), "[Acetyl]-PE[+15.994915]P(TI[+2.0]D)[-18.010565][+1.0]EK")
  TEST_EXCEPTION(Exception::ParseError, FLASHHelperClasses::toPeptidoform("PEPTXDE", {}))

  // the tags FLASHTnT writes: a mass with an accession as an alternative, also one OpenMS::ProForma does not read
  const std::vector<PFM> with_accession = {residueMod(1, 1, 42.010565, "+42.010565|MOD:00064")};
  TEST_EQUAL(ProForma::toString(FLASHHelperClasses::toPeptidoform("AKK", with_accession)), "AK[+42.010565|MOD:00064]K")
  const std::vector<PFM> without_accession = {residueMod(1, 1, 42.010565, "+42.010565|-")};
  TEST_EQUAL(ProForma::toString(FLASHHelperClasses::toPeptidoform("AKK", without_accession)), "AK[+42.010565|-]K")
}
END_SECTION

START_SECTION(static double getProFormaTagMass(const std::string& tag, char residue = 0))
{
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("Oxidation", 'M'), 15.994915)
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("UNIMOD:35"), 15.994915)
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("UniMod:35"), 15.994915)
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("MOD:00719"), 15.994915)
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("Formula:O"), 15.994915)
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("Glycan:Hex1HexNAc1"), 365.132197)
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("INFO:x|Phospho|+79.97", 'S'), 79.966331)
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("Obs:+79.97"), 79.97)
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("Cation:Ca[II]", 'K'), 37.946941)
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("+14.01565|MOD:00080", 'Q'), 14.01565)
  TEST_REAL_SIMILAR(FLASHHelperClasses::getProFormaTagMass("2-methylbutyryl", 'K'), 84.057515)
  TEST_EXCEPTION(Exception::ParseError, FLASHHelperClasses::getProFormaTagMass("U:35"))
  TEST_EXCEPTION(Exception::ParseError, FLASHHelperClasses::getProFormaTagMass("NoSuchModification"))
}
END_SECTION

START_SECTION(neutral mass: residues and water plus all modifications)
{
  std::string sequence;
  const auto mods = read("{Glycan:Hex}[Acetyl]-PEM[Oxidation]K", sequence);
  double mass = AASequence::fromString(sequence).getMonoWeight();
  for (const auto& m : mods) mass += m.mass;
  const double expected = AASequence::fromString("PEMK").getMonoWeight() + 162.052824 + 42.010565 + 15.994915;
  TEST_REAL_SIMILAR(mass, expected)
}
END_SECTION

START_SECTION(static std::string formatProFormaMass(double mass))
{
  TEST_EQUAL(FLASHHelperClasses::formatProFormaMass(15.994915), "+15.994915")
  TEST_EQUAL(FLASHHelperClasses::formatProFormaMass(-18.010565), "-18.010565")
  TEST_EQUAL(FLASHHelperClasses::formatProFormaMass(80.9787), "+80.9787")
  TEST_EQUAL(FLASHHelperClasses::formatProFormaMass(1.0), "+1.0")
}
END_SECTION

END_TEST
