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
#include <OpenMS/ANALYSIS/TOPDOWN/ProteoformDatabase.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassLadderAlgorithm.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/FORMAT/FLASHTnTFile.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
// the default of FLASHTnT -peff_mod_masses
const std::vector<String> default_mod_masses {"MOD:01819=100.016044", "MOD:01893=86.000394", "MOD:01614=329.05252", "MOD:01956=15.994915",
                                              "MOD:00657=14.999666",  "MOD:00017=skip",       "MOD:00223=skip",      "MOD:01976=skip",
                                              "MOD:01780=skip"};

const String lpp_seq = "MKATKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNMATKYRK";
const FASTAFile::FASTAEntry lpp("sp:P69776",
                                "\\PName=(Major outer membrane lipoprotein Lpp) \\GName=lpp \\NcbiTaxId=83333 \\TaxName=Escherichia coli (strain K12) "
                                "\\Length=78 \\SV=1 \\EV=157 \\PE=1 \\ID=LPP_ECOLI \\AltAC=(P02937)(Q53272) \\ModResPsi=(0:1|MOD:00482|N-formyl-L-methionine "
                                "(Met))(1:21|MOD:00069|N-palmitoyl-L-cysteine)(2:21|MOD:00116|S-diacylglycerol-L-cysteine)(3:78|MOD:00223|N6-mureinyl-L-lysine) "
                                "\\VariantSimple=(4:29|E) \\Processed=(5:1|20|PEFF:0001021|signal peptide)(6:21|78|PEFF:0001020|mature protein)",
                                lpp_seq);

const String dsba_seq = "MKKIWLALAGLVLAFSASAAQYEDGKQYTTLEKPVAGAPQVLEFFSFFCPHCYQFEEVLHISDNVKKKLPEGVKMTKYHVNFMGGDLGKDLTQAWAVAMALGVEDKVTVPLFEGVQKTQT"
                        "IRSASDIRDVFINAGIKGEEYDAAWNSFVVKSLVAQQEKAAADVQLRGVPAMFVNGKYQLNPQGMDTSNMDVFVQQYADTVKYLSEKK";
const String dsba_keys = "\\PName=(Thiol:disulfide interchange protein DsbA) \\GName=dsbA \\NcbiTaxId=83333 \\TaxName=Escherichia coli (strain K12) "
                         "\\Length=208 \\ID=DSBA_ECOLI \\ModResPsi=(0:49|MOD:00798|half cystine)(1:52|MOD:00798|half cystine) "
                         "\\Processed=(2:1|19|PEFF:0001021|signal peptide)(3:20|208|PEFF:0001020|mature protein)";

const String if3_seq = "MKGGKRVQTARPNRINGEIRAQEVRLTGLEGEQLGIVSLREALEKAEEAGVDLVEISPNAEPPVCRIMDYGKFLYEKSKSSKEQKKKQKVIQVKEIKFRPGTDEGDYQVKLRSLIRFLEE"
                       "GDKAKITLRFRGREMAHQQIGMEVLNRVKDDLQELAVVESFPTKIEGRQMIMVLAPKKKQ";
const FASTAFile::FASTAEntry if3("sp:P0A707",
                                "\\PName=(Translation initiation factor IF-3) \\ID=IF3_ECOLI \\ModResPsi=(0:1|MOD:00073|N-methyl-L-methionine) "
                                "\\Processed=(1:1|180|PEFF:0001020|mature protein)(2:1|1|PEFF:0001035|initiator methionine)(3:2|180|PEFF:0001020|mature "
                                "protein)(4:7|180|PEFF:0001020|mature protein)",
                                if3_seq);

double residueSum(const String& seq)
{
  double sum = 0;
  for (char c : seq) sum += ResidueDB::getInstance()->getResidue(c)->getMonoWeight(Residue::Internal);
  return sum;
}

// the mass of an entry: its letters and its fixed modifications
double entryMass(const ProteoformDatabase& db, const std::vector<FASTAFile::FASTAEntry>& fasta, Size e)
{
  double mass = residueSum(fasta[e].sequence);
  for (const auto& fm : db.entries[e].fixed) mass += fm.mass;
  return mass;
}

String sorted(String s)
{
  std::sort(s.begin(), s.end());
  return s;
}
} // namespace

START_TEST(ProteoformDatabase, "$Id$")

START_SECTION(static std::uint64_t fnv1a64(const String& s))
{
  TEST_EQUAL(ProteoformDatabase::fnv1a64(""), 14695981039346656037ULL)
  TEST_EQUAL(ProteoformDatabase::fnv1a64("a"), 12638187200555641996ULL)
}
END_SECTION

START_SECTION(static ProteoformDatabase fromFASTA(const std::vector<FASTAFile::FASTAEntry>& fasta_entries))
{
  std::vector<FASTAFile::FASTAEntry> fasta {{"P1", "one", "PEPTIDEK"}, {"DECOY_P1", "one", "KEDITPEP"}, {"P2", "two", "AAAK"}};
  const auto db = ProteoformDatabase::fromFASTA(fasta);
  TEST_EQUAL(db.isIdentity(), true)
  TEST_EQUAL(db.hasFixedModifications(), false)
  TEST_EQUAL(db.proteins.size(), 3)
  TEST_EQUAL(db.entries.size(), 3)
  TEST_EQUAL(db.proteins[0].partner, 1)
  TEST_EQUAL(db.proteins[1].partner, 0)
  TEST_EQUAL(db.proteins[2].partner, -1)
  TEST_EQUAL(db.proteins[1].decoy, true)
  TEST_EQUAL(db.proteinOfEntry()[2], 2)
}
END_SECTION

START_SECTION(void build(const std::vector<FASTAFile::FASTAEntry>& input, int max_mods, const std::vector<String>& mod_masses, std::vector<FASTAFile::FASTAEntry>& fasta_entries))
{
  // keys, escaping and the FASTA-like identifier and description
  {
    std::vector<FASTAFile::FASTAEntry> input {{"sp:P00001",
                                               "\\PName=(Kitchen-sink protein (test) \\| alpha\\\\beta :-\\))(Second name) \\GName=tst \\NcbiTaxId=83333 "
                                               "\\TaxName=Escherichia coli (strain K12) \\PE=1 \\SV=2 \\ID=TEST_ECOLI \\ModResPsi=(?|MOD:00046|O-phospho-L-serine) "
                                               "\\ModRes=(3||Blocked amino end (Met))",
                                               "MSSKPEPTIDE"}};
    ProteoformDatabase db;
    std::vector<FASTAFile::FASTAEntry> fasta;
    db.build(input, 1, default_mod_masses, fasta);
    TEST_EQUAL(db.proteins.size(), 2)
    TEST_EQUAL(db.proteins[0].identifier, "sp|P00001|TEST_ECOLI")
    TEST_EQUAL(db.proteins[0].description, "Kitchen-sink protein (test) | alpha\\beta :-) OS=Escherichia coli (strain K12) OX=83333 GN=tst PE=1 SV=2")
    TEST_EQUAL(db.proteins[1].identifier, "DECOY_sp|P00001|TEST_ECOLI")
    // '?' positions and generic ModRes are left out: one entry, no fixed modification
    TEST_EQUAL(db.entries.size(), 2)
    TEST_EQUAL(db.entries[0].fixed.size(), 0)
    TEST_EQUAL(fasta[0].identifier, "sp|P00001|TEST_ECOLI")
    TEST_EQUAL(fasta[0].sequence, "MSSKPEPTIDE")
    TEST_EQUAL(db.isIdentity(), true) // one entry per protein over its length, nothing fixed
  }

  // Lpp: fMet only on the full length (N-terminal rule), the lipid pair one unit on the mature form, murein-K skipped (mod_masses)
  {
    ProteoformDatabase db;
    std::vector<FASTAFile::FASTAEntry> fasta;
    db.build({lpp}, 1, default_mod_masses, fasta);
    TEST_EQUAL(db.proteins.size(), 2)
    TEST_EQUAL(db.entries.size(), 12) // [1,78] {} {fMet}; [21,78] {} {pal} {dag} {pal dag}; and the decoy's
    TEST_EQUAL(db.isIdentity(), false)
    TEST_EQUAL(db.hasFixedModifications(), true)
    TEST_EQUAL(db.getLoadStatistics().unresolved.at("MOD:00223"), 1)
    TEST_EQUAL(db.entries[0].offset, 0)
    TEST_EQUAL(db.entries[0].fixed.size(), 0)
    TEST_EQUAL(db.entries[1].fixed.size(), 1)
    TEST_EQUAL(db.entries[1].fixed[0].accession, "MOD:00482")
    TEST_EQUAL(db.entries[1].fixed[0].kind, 2)
    TEST_EQUAL(db.entries[1].fixed[0].position, 0)
    TEST_EQUAL(db.entries[2].offset, 20)
    TEST_EQUAL(db.entries[2].length, 58)
    TEST_EQUAL(fasta[2].sequence, lpp_seq.substr(20))
    TEST_EQUAL(db.entries[2].id, "sp|P69776|LPP_ECOLI 21-78")
    TEST_EQUAL(db.entries[5].fixed.size(), 2)
    TEST_EQUAL(db.fixedMasses()[5].size(), 1) // the pair on one residue
    TEST_REAL_SIMILAR(db.fixedMasses()[5][0].second, 814.74)
    // the pair sits on two sites of C21 (PSI-MOD TermSpec): N-palmitoyl on the amino group, S-diacylglycerol on the side chain
    TEST_EQUAL(db.entries[5].fixed[0].sites, 2)
    TEST_EQUAL(db.entries[5].fixed[1].sites, 1)
    // occupied: the N-terminus (kind 2) and, by S-diacylglycerol, the side chain; fMet only the N-terminus
    TEST_EQUAL((int)db.occupancy()[5][0].second, 3)
    TEST_EQUAL((int)db.occupancy()[4][0].second, 3)
    TEST_EQUAL((int)db.occupancy()[3][0].second, 2)
    TEST_EQUAL((int)db.occupancy()[1][0].second, 2)
    TEST_EQUAL(db.entries[5].proforma, "C[+238.229666|MOD:00069][+576.511761|MOD:00116]SSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNMATKYRK")
    TEST_EQUAL(db.proteins[0].forms.size(), 2)
    TEST_EQUAL(db.proteins[0].forms[1].first, 21)
    TEST_EQUAL(db.proteins[0].forms[1].second, 78)
    // the decoy: right after the target, its entries mirror the target's
    TEST_EQUAL(db.proteins[0].partner, 1)
    TEST_EQUAL(db.proteins[1].partner, 0)
    TEST_EQUAL(db.proteins[1].entries.size(), 6)
    for (Size k = 0; k < 6; k++)
    {
      const Size t = k, d = 6 + k;
      TEST_EQUAL(db.entries[d].protein, 1)
      TEST_EQUAL(db.entries[d].offset, db.entries[t].offset)
      TEST_EQUAL(sorted(fasta[d].sequence), sorted(fasta[t].sequence))  // composition per form
      TEST_EQUAL(fasta[d].sequence[0], fasta[t].sequence[0])            // first and last residue of a form
      TEST_EQUAL(fasta[d].sequence.back(), fasta[t].sequence.back())
      TEST_REAL_SIMILAR(entryMass(db, fasta, d), entryMass(db, fasta, t))
      TEST_EQUAL(db.entries[d].fixed.size(), db.entries[t].fixed.size())
      for (Size m = 0; m < db.entries[d].fixed.size(); m++) // moved with their residues
      {
        TEST_EQUAL(fasta[d].sequence[db.entries[d].fixed[m].position], fasta[t].sequence[db.entries[t].fixed[m].position])
      }
    }
    // residues 1, 20 (end of nothing), 21 and 78 are form boundaries: 1, 21 and 78 stay
    TEST_EQUAL(db.proteins[1].sequence[0], 'M')
    TEST_EQUAL(db.proteins[1].sequence[20], 'C')
    TEST_EQUAL(db.proteins[1].sequence[77], 'K')
    // deterministic
    ProteoformDatabase db2;
    std::vector<FASTAFile::FASTAEntry> fasta2;
    db2.build({lpp}, 1, default_mod_masses, fasta2);
    TEST_EQUAL(db2.proteins[1].sequence, db.proteins[1].sequence)
    TEST_EQUAL(db.proteins[1].sequence, "MAALITVSKGALLTGLKGVLCDMQVASADRQSTNSVAAVNLSQKDNVNSRAKQSNDKDNDLKQLAARALMAYITDRDK") // fixed seed
    // k = 0: the forms only; k = 2: the same (no form has two units)
    db2.build({lpp}, 0, default_mod_masses, fasta2);
    TEST_EQUAL(db2.entries.size(), 4)
    db2.build({lpp}, 2, default_mod_masses, fasta2);
    TEST_EQUAL(db2.entries.size(), 12)
  }

  // IF-3: forms 1-180, 2-180 (initiator Met), 7-180; N-methyl-Met only on 1-180
  {
    ProteoformDatabase db;
    std::vector<FASTAFile::FASTAEntry> fasta;
    db.build({if3}, 1, default_mod_masses, fasta);
    TEST_EQUAL(db.entries.size(), 8)
    TEST_EQUAL(db.entries[1].fixed.size(), 1)
    TEST_EQUAL(db.entries[2].offset, 1)
    TEST_EQUAL(db.entries[3].offset, 6)
    TEST_EQUAL(db.entries[3].fixed.size(), 0)
    TEST_EQUAL(db.proteins[0].description, "Translation initiation factor IF-3")
  }

  // DsbA: a disulfide bond only with its pair (Option B), all bonds of a form one unit
  {
    ProteoformDatabase db;
    std::vector<FASTAFile::FASTAEntry> fasta;
    db.build({{"sp:P0AEG4", dsba_keys + " \\DisulfideBond=(4:0,1)", dsba_seq}}, 1, default_mod_masses, fasta);
    TEST_EQUAL(db.entries.size(), 8) // [1,208] {} {bond}; [20,208] {} {bond}; the decoy's
    TEST_EQUAL(db.entries[1].fixed.size(), 2)
    TEST_EQUAL(db.entries[1].fixed[0].position, 48)
    TEST_EQUAL(db.entries[1].fixed[1].position, 51)
    TEST_EQUAL(db.entries[1].fixed[0].bond, 51) // each half cystine knows the other
    TEST_EQUAL(db.entries[1].fixed[1].bond, 48)
    TEST_EQUAL(db.entries[3].fixed[0].bond, 32)
    TEST_REAL_SIMILAR(db.entries[1].fixed[0].mass, -1.007825)
    TEST_EQUAL(db.entries[3].offset, 19)
    TEST_EQUAL(db.entries[3].fixed[0].position, 29)
    TEST_EQUAL(db.getLoadStatistics().disulfide_bonds, 1)
    // without the pairs: the half cystines are left out
    db.build({{"sp:P0AEG4", dsba_keys, dsba_seq}}, 1, default_mod_masses, fasta);
    TEST_EQUAL(db.entries.size(), 4)
    TEST_EQUAL(db.getLoadStatistics().unpaired_half_cystines, 2)
  }

  // alternatives on one residue: a nonempty subset on distinct sites is one unit; a C-terminal residue's modification is kind 3
  {
    ProteoformDatabase db;
    std::vector<FASTAFile::FASTAEntry> fasta;
    db.build({{"sp:P00002", "\\ID=ALT_ECOLI \\ModResPsi=(0:3|MOD:00064|N6-acetyl-L-lysine)(1:3|MOD:01819|N6-succinyl-L-lysine)(2:6|MOD:00064|N6-acetyl-L-lysine)",
               "MAKLLK"}},
             1, default_mod_masses, fasta);
    // K3: {ac} {su} (both on the epsilon-amino group: never together); K6 (the last residue): {ac}; and the empty set
    TEST_EQUAL(db.proteins[0].entries.size(), 4)
    TEST_EQUAL(db.entries[1].fixed.size(), 1)
    TEST_REAL_SIMILAR(db.entries[2].fixed[0].mass, 100.016044)
    TEST_EQUAL(db.entries[3].fixed[0].kind, 3)
    TEST_EQUAL(db.entries[1].fixed[0].kind, 1)
    // N6-acetyl-K6 at the C-terminus: the C-terminus and the side chain (TermSpec none) occupied
    TEST_EQUAL((int)db.occupancy()[3][0].second, 5)
    TEST_EQUAL((int)db.occupancy()[1][0].second, 1)
    TEST_EQUAL(db.entries[3].proforma, "MAKLLK[+42.010565|MOD:00064]")
    db.build({{"sp:P00002", "\\ID=ALT_ECOLI \\ModResPsi=(0:3|MOD:00064|N6-acetyl-L-lysine)(1:3|MOD:01819|N6-succinyl-L-lysine)(2:6|MOD:00064|N6-acetyl-L-lysine)",
               "MAKLLK"}},
             2, default_mod_masses, fasta);
    TEST_EQUAL(db.proteins[0].entries.size(), 6) // + {K3 choice} x {K6}
    // MoaD-like C-terminal Gly: 1-thioglycine (C-term) and glycyl adenylate (not in the bundled PSI-MOD: side chain and C-terminus at a
    // form end) are alternatives on the carboxyl group
    db.build({{"sp:P00003", "\\ID=MOAD_ECOLI \\ModResPsi=(0:6|MOD:00270|1-thioglycine)(1:6|MOD:01614|glycyl adenylate)", "MAKLGG"}}, 1,
             default_mod_masses, fasta);
    TEST_EQUAL(db.proteins[0].entries.size(), 3)
    TEST_EQUAL(db.entries[1].fixed[0].sites, 4)
    TEST_EQUAL(db.entries[2].fixed[0].sites, 5)
    // an annotation on a residue that is not its modification's origin is left out
    db.build({{"sp:P00004", "\\ID=WRONG_ECOLI \\ModResPsi=(0:2|MOD:00064|N6-acetyl-L-lysine)(1:3|MOD:00064|N6-acetyl-L-lysine)", "MAKLLK"}}, 1,
             default_mod_masses, fasta);
    TEST_EQUAL(db.proteins[0].entries.size(), 2)
    TEST_EQUAL(db.getLoadStatistics().wrong_residue.at("MOD:00064"), 1)
  }

  // decoys per protein: a protein with a decoy in the file gets none, the others one each, right after them
  {
    ProteoformDatabase db;
    std::vector<FASTAFile::FASTAEntry> fasta;
    db.build({lpp, {"sp|Q1|PLAIN", "plain", "MPEPTIDEKR"}, {"DECOY_sp|Q1|PLAIN", "plain", "MEDITPEKPR"}}, 1, default_mod_masses, fasta);
    TEST_EQUAL(db.proteins.size(), 4)
    TEST_EQUAL(db.proteins[1].identifier, "DECOY_sp|P69776|LPP_ECOLI")
    TEST_EQUAL(db.proteins[2].identifier, "sp|Q1|PLAIN")
    TEST_EQUAL(db.proteins[3].identifier, "DECOY_sp|Q1|PLAIN")
    TEST_EQUAL(db.proteins[3].sequence, "MEDITPEKPR")
    TEST_EQUAL(db.proteins[2].partner, 3)
    TEST_EQUAL(db.getLoadStatistics().generated_decoys, 1)
    TEST_EQUAL(db.getLoadStatistics().unmirrored_decoys, 0)
  }

  // decoys of the file are paired by the identifier as rebuilt, else by the accession; a target with one gets none generated
  {
    ProteoformDatabase db;
    std::vector<FASTAFile::FASTAEntry> fasta;
    // a FASTA decoy of the PEFF target sp:P69776 (\\ID=LPP_ECOLI): named after sp|P69776|LPP_ECOLI
    db.build({lpp, {"DECOY_sp|P69776|LPP_ECOLI", "", "MKATKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNMATKYRK"}}, 1, default_mod_masses,
             fasta);
    TEST_EQUAL(db.proteins.size(), 2)
    TEST_EQUAL(db.getLoadStatistics().generated_decoys, 0)
    TEST_EQUAL(db.proteins[0].partner, 1)
    TEST_EQUAL(db.proteins[1].partner, 0)
    TEST_EQUAL(db.getLoadStatistics().unmirrored_decoys, 1) // one entry against six
    // named after the raw identifier, without PEFF keys: by its accession
    db.build({lpp, {"DECOY_sp:P69776", "", "MKATKLVLGAVILGSTLLAGCSSNAKIDQLSSDVQTLNAKVDQLSNDVNAMRSDVQAAKDDAARANQRLDNMATKYRK"}}, 1, default_mod_masses, fasta);
    TEST_EQUAL(db.proteins.size(), 2)
    TEST_EQUAL(db.proteins[0].partner, 1)
    TEST_EQUAL(db.proteins[1].partner, 0)
    // a second decoy of a target (paired with it, the target with the first, as in a FASTA search), a decoy without a target, and an
    // entry without residues
    db.build({{"sp|Q1|PLAIN", "plain", "MPEPTIDEKR"}, {"DECOY_sp|Q1|PLAIN", "plain", "MEDITPEKPR"}, {"DECOYsp|Q1|PLAIN", "plain", "MEDITPEKPR"},
              {"DECOY_sp|Q2|OTHER", "other", "MEDITPEKPR"}, {"sp|Q3|EMPTY", "empty", ""}},
             1, default_mod_masses, fasta);
    TEST_EQUAL(db.proteins.size(), 4)
    TEST_EQUAL(db.getLoadStatistics().unpaired_decoys, 1)
    TEST_EQUAL(db.getLoadStatistics().empty_sequences, 1)
    TEST_EQUAL(db.getLoadStatistics().generated_decoys, 0)
    TEST_EQUAL(db.proteins[0].partner, 1)
    TEST_EQUAL(db.proteins[2].partner, 0)
    TEST_EQUAL(db.proteins[3].partner, -1)
    TEST_EQUAL(db.isIdentity(), true)
  }
}
END_SECTION

START_SECTION(static bool isPEFF(const String& file))
{
  TEST_EQUAL(ProteoformDatabase::isPEFF("db.peff"), true)
  TEST_EQUAL(ProteoformDatabase::isPEFF(OPENMS_GET_TEST_DATA_PATH("FastaIterator_test.fasta")), false)
  TEST_EQUAL(ProteoformDatabase::isPEFF(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta")), true) // its first line: # PEFF
}
END_SECTION

// the components that read the database entries

START_SECTION(FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs with fixed modifications)
{
  // N6-acetyl-K at position 1: the flanking residue mass of a tag after it includes the acetyl
  const String seq = "AKPEPTLDE";
  const double ak = residueSum("AK"), acetyl = 42.010565;
  const FLASHHelperClasses::Tag tag("PEP", ak + acetyl, -1, {100.0, 197.0, 326.0, 423.0}, {1, 1, 1, 1}, 1);
  std::vector<int> positions;
  std::vector<double> masses;
  FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(positions, masses, -1, seq, tag);
  TEST_EQUAL(positions.size(), 1)
  TEST_REAL_SIMILAR(masses[0], acetyl)
  const std::vector<std::pair<int, double>> fixed {{1, acetyl}};
  positions.clear();
  masses.clear();
  FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(positions, masses, -1, seq, tag, &fixed);
  TEST_EQUAL(positions.size(), 1)
  TOLERANCE_ABSOLUTE(1e-6)
  TEST_REAL_SIMILAR(masses[0], 0.0)
  TOLERANCE_ABSOLUTE(1e-5)
}
END_SECTION

START_SECTION(FLASHTaggerAlgorithm::scoreShifts with a density mass and the best count)
{
  // residue sums 0, 100, ..., 1000; masses that land on 300, 500 and 700 at shift 0, on nothing at shift 50
  std::vector<double> sums;
  for (int k = 0; k <= 10; k++) sums.push_back(100.0 * k);
  const std::vector<double> masses {0.0, 300.0, 500.0, 700.0};
  const std::vector<double> shifts {50.0, 0.0}, flanks {0.0, 0.0};
  int count = -1;
  const double score = FLASHTaggerAlgorithm::scoreShifts(masses, sums, shifts.data(), flanks.data(), 2, 10.0, 0.0, &count);
  TEST_EQUAL(count, 4) // the terminus (0) too
  TEST_REAL_SIMILAR(FLASHTaggerAlgorithm::scoreShifts(masses, sums, shifts.data(), flanks.data(), 2, 10.0), score)
  // the same sums spread over less mass (the residues without fixed masses): a match is likelier by chance
  TEST_EQUAL(FLASHTaggerAlgorithm::scoreShifts(masses, sums, shifts.data(), flanks.data(), 2, 10.0, 900.0, &count) < score, true)
  TEST_EQUAL(count, 4)
}
END_SECTION

START_SECTION(FLASHTaggerAlgorithm::scoreShifts with a floor)
{
  // random proteins, spectra and shifts: with a floor the score is the same, or -infinity (with a best count of 0) if it is below the floor
  std::mt19937 rng(5);
  std::uniform_real_distribution<double> residue(57.0, 186.0), unit(0.0, 1.0);
  Size below = 0, same = 0;
  for (int trial = 0; trial < 300; trial++)
  {
    std::vector<double> sums {0.0};
    const int length = 50 + (int)(rng() % 300);
    for (int k = 0; k < length; k++) sums.push_back(sums.back() + residue(rng));
    std::vector<std::uint32_t> index;
    FLASHTaggerAlgorithm::indexResidueSums(sums, index);
    // masses: some residue sums of the protein at one shift, and random ones
    const double shift = (double)(rng() % 3) * 1000.0 * unit(rng);
    std::vector<double> masses;
    const int planted = (int)(rng() % 12);
    for (int k = 0; k < planted; k++) masses.push_back(sums[1 + rng() % (sums.size() - 1)] - shift);
    for (int k = 0; k < 40; k++) masses.push_back(300.0 + 20000.0 * unit(rng));
    masses.erase(std::remove_if(masses.begin(), masses.end(), [](double m) { return m <= 0; }), masses.end());
    masses.push_back(0.0); // the terminus
    std::sort(masses.begin(), masses.end());
    std::vector<double> shifts {shift}, flanks {0.0};
    for (int k = 0; k < 6; k++)
    {
      shifts.push_back(2000.0 * unit(rng) - 1000.0);
      flanks.push_back(100.0 * unit(rng));
    }
    int count = -1;
    const double score = FLASHTaggerAlgorithm::scoreShifts(masses, sums, index, shifts.data(), flanks.data(), shifts.size(), 10.0, 0.0, &count);
    for (double score_floor : {0.0, score - 1.0, score, score + 1e-3, score + 0.5, score + 3.0, 1000.0})
    {
      int floor_count = -1, plain_count = -1;
      const double with_floor = FLASHTaggerAlgorithm::scoreShifts(masses, sums, index, shifts.data(), flanks.data(), shifts.size(), 10.0, 0.0,
                                                                  &floor_count, score_floor);
      const double plain = FLASHTaggerAlgorithm::scoreShifts(masses, sums, shifts.data(), flanks.data(), shifts.size(), 10.0, 0.0, &plain_count, score_floor);
      TEST_EQUAL(with_floor == plain || (std::isinf(with_floor) && std::isinf(plain)), true)
      if (with_floor == -std::numeric_limits<double>::infinity())
      {
        below++;
        TEST_EQUAL(score < score_floor, true)
        TEST_EQUAL(floor_count, 0)
      }
      else
      {
        same++;
        TEST_EQUAL(with_floor == score, true) // the same double
        TEST_EQUAL(floor_count, count)
      }
    }
  }
  TEST_EQUAL(below > 300, true)
  TEST_EQUAL(same > 300, true)
}
END_SECTION

START_SECTION(MassLadderAlgorithm precursor windows of several entries per protein)
{
  // the windows of the processed forms of IF-3 (1-180, 2-180, 7-180) are one hypothesis each, as in the FASTA database
  const double h2o = Residue::getInternalToFull().getMonoWeight();
  String reversed = if3_seq;
  std::reverse(reversed.begin(), reversed.end());
  const FASTAFile::FASTAEntry decoy("DECOY_sp:P0A707", "", reversed);
  const std::vector<FASTAFile::FASTAEntry> fasta {{"sp:P0A707", "", if3_seq}, decoy};
  const double mass1 = residueSum(if3_seq.substr(10, 20)) + h2o;                // residues 11-30
  const double mass2 = residueSum(if3_seq.substr(4, 20)) + 42.010565 + h2o;      // residues 5-24, K5 acetylated
  auto windows = [&](const ProteoformDatabase& db, const std::vector<FASTAFile::FASTAEntry>& entries, double precursor_mass) {
    std::vector<std::string> seqs;
    for (const auto& e : entries) seqs.push_back(e.sequence);
    MassLadderAlgorithm ladders;
    ladders.buildIndex(seqs, {"b", "y"}, false, true, &db);
    std::vector<MassLadderAlgorithm::PrecursorLadderResult> results;
    ladders.searchPrecursorLadders({500.0, 1000.0}, {precursor_mass}, 10.0, 10.0, 5, results);
    return results[0].windows;
  };
  const auto identity = ProteoformDatabase::fromFASTA(fasta);
  const Size identity1 = windows(identity, fasta, mass1), identity2 = windows(identity, fasta, mass2);
  TEST_EQUAL(identity1 >= 1, true)

  const String forms = "\\ID=IF3_ECOLI \\Processed=(1:1|180|PEFF:0001020|mature protein)(2:1|1|PEFF:0001035|initiator methionine)(4:7|180|PEFF:0001020|mature protein)";
  ProteoformDatabase db;
  std::vector<FASTAFile::FASTAEntry> entries;
  db.build({{"sp:P0A707", forms, if3_seq}, decoy}, 1, default_mod_masses, entries);
  TEST_EQUAL(entries.size(), 4) // three forms and the decoy of the file
  TEST_EQUAL(windows(db, entries, mass1), identity1)
  // with N6-acetyl-K5 (in 1-180 and 2-180): one more hypothesis at the acetylated mass, the window with the acetyl
  db.build({{"sp:P0A707", forms + " \\ModResPsi=(5:5|MOD:00064|N6-acetyl-L-lysine)", if3_seq}, decoy}, 1, default_mod_masses, entries);
  TEST_EQUAL(entries.size(), 6)
  TEST_EQUAL(windows(db, entries, mass1), identity1)
  TEST_EQUAL(windows(db, entries, mass2), identity2 + 1)
}
END_SECTION

START_SECTION(FLASHTnTFile::getProFormaString with fixed modifications)
{
  ProteinHit hit;
  hit.setSequence("MAKLLKQ");
  hit.setMetaValue("StartPosition", 1);
  hit.setMetaValue("EndPosition", 7);
  hit.setMetaValue("Modifications", std::vector<double> {42.010565, 14.999666});
  hit.setMetaValue("ModificationStarts", std::vector<int> {2, 6});
  hit.setMetaValue("ModificationEnds", std::vector<int> {2, 6});
  hit.setMetaValue("ModificationIDs", std::vector<String> {"N6-acetyl-L-lysine", "L-glutamic acid 5-methyl ester (Gln)"});
  hit.setMetaValue("ModificationACCs", std::vector<String> {"MOD:00064", "MOD:00657"});
  hit.setMetaValue("ModificationVariable", std::vector<int> {2, 2});
  // the mass searched, with the accession as an alternative (MOD:00657 gives 0.984 in the bundled PSI-MOD)
  TEST_EQUAL(FLASHTnTFile::getProFormaString(hit), "MAK[+42.010565|MOD:00064]LLKQ[+14.999666|MOD:00657]")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
