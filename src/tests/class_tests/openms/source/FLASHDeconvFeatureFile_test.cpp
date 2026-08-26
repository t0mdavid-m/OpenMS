// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Locks the TopFD/TopPIC *_ms2.feature layout that FLASHDeconv emits. TopPIC parses this file
// POSITIONALLY and without bounds checks (SpecFeature::SpecFeature), so a wrong column count is
// an access violation rather than an error message, and a shifted column is a silently wrong
// result. See docs/adr/0034 for why 1.8.x is the pinned target.
//
// The fixture is built to FAIL against the pre-fix writer, which is the only reason it is worth
// having. Two properties do that work:
//
//   1. mass_features carries a feature that is NOT an MS1 target (f3, ms_level 2). Synthetic
//      feature IDs used to be minted from the MS1-TARGET maximum, so they landed inside that
//      occupied index space and the ms_level == 2 branch -- which tells minted from real on
//      `index < mass_features.size() + 1` -- could not distinguish them. With only MS1 targets in
//      the vector the old and new writers agree byte for byte and the test proves nothing.
//
//   2. The writer is called for ms_level 1 BEFORE ms_level 2, as FLASHDeconv.cpp does. Minting
//      happens in the MS1 pass and mutates deconvolved_spectra; run the MS2 pass first and no
//      precursor has been minted an ID at all.
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FLASHDeconvFeatureFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/METADATA/Precursor.h>

#include <fstream>
#include <map>
#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
  /// Split on tabs WITHOUT dropping trailing empties: N tabs must always yield N+1 fields, or a
  /// truncated row would tokenize one column short and read as well-formed.
  std::vector<String> splitTabs(const String& line)
  {
    std::vector<String> out;
    std::string cur;
    for (char c : line)
    {
      if (c == '\t')
      {
        out.push_back(cur);
        cur.clear();
      }
      else
        cur += c;
    }
    out.push_back(cur);
    return out;
  }

  std::vector<String> readLines(const String& path)
  {
    std::vector<String> lines;
    std::ifstream in(path.c_str());
    std::string l;
    while (std::getline(in, l))
    {
      if (! l.empty() && l.back() == '\r') l.pop_back();
      if (! l.empty()) lines.push_back(String(l));
    }
    return lines;
  }

  FLASHHelperClasses::MassFeature makeFeature(uint index, uint ms_level, double mono_mass, int rep_charge, double rt_start)
  {
    std::vector<Peak2D> peaks;
    for (int i = 0; i < 3; ++i)
    {
      Peak2D p;
      p.setRT(rt_start + 10.0 * i);
      p.setMZ(mono_mass); // FLASHDeconv traces run in mass space, not m/z
      p.setIntensity(1000.0f * (i == 1 ? 2.0f : 1.0f)); // apex in the middle
      peaks.push_back(p);
    }
    MassTrace mt(peaks);
    mt.updateWeightedMeanMZ();

    FLASHHelperClasses::MassFeature f;
    f.index = index;
    f.ms_level = ms_level;
    f.is_decoy = false;
    f.mt = mt;
    f.iso_offset = 0;
    f.avg_mass = mono_mass + 0.6; // a plausible averagine delta; only its sign matters here
    f.rep_charge = rep_charge;
    // mirrors MassFeatureTrace: the representative charge's AVERAGE m/z, protons included
    f.rep_mz = f.avg_mass / rep_charge + Constants::PROTON_MASS_U;
    f.min_charge = rep_charge - 1;
    f.max_charge = rep_charge + 1;
    f.charge_count = 3;
    f.scan_number = 1;
    f.min_scan_number = 1;
    f.max_scan_number = 3;
    f.isotope_score = 0.9;
    f.qscore = 0.8;
    return f;
  }

  /// A precursor peak group: two peaks at the representative charge, the +1 isotope the more
  /// intense, so the observed rep m/z is distinguishable from the monoisotopic m/z.
  PeakGroup makePrecursorPG(double mono_mass, int z, int ms1_scan, uint feature_index)
  {
    PeakGroup pg(z, z, true);
    const double mono_mz = mono_mass / z + Constants::PROTON_MASS_U;
    for (int iso = 0; iso < 2; ++iso)
    {
      Peak1D p1d(mono_mz + iso * Constants::ISOTOPE_MASSDIFF_55K_U / z, iso == 0 ? 500.0f : 900.0f);
      FLASHHelperClasses::LogMzPeak lp(p1d, true);
      lp.abs_charge = z;
      lp.isotopeIndex = iso;
      pg.push_back(lp);
    }
    pg.setAbsChargeRange(z, z);
    pg.setRepAbsCharge(z);
    pg.setMonoisotopicMass(mono_mass);
    pg.setScanNumber(ms1_scan);
    pg.setFeatureIndex(feature_index);
    return pg;
  }
} // namespace

START_TEST(FLASHDeconvFeatureFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/// ---- fixture -------------------------------------------------------------------------------
// f3 is the trap: an MSn feature occupying index 3, which the old minting scheme would collide with.
std::vector<FLASHHelperClasses::MassFeature> mass_features;
mass_features.push_back(makeFeature(1, 1, 5000.0, 5, 100.0));  // traced MS1 target
mass_features.push_back(makeFeature(2, 1, 6000.0, 7, 200.0));  // traced MS1 target
mass_features.push_back(makeFeature(3, 2, 9999.0, 9, 3000.0)); // MSn feature -- never in _ms1.feature

const int kMs1Scan = 1;
const double kMs1Rt = 60.0; // seconds -> 1.000000 min

MSSpectrum ms1_spec;
ms1_spec.setRT(kMs1Rt);
ms1_spec.setMSLevel(1);
DeconvolvedSpectrum ms1(kMs1Scan);
ms1.setOriginalSpectrum(ms1_spec);

Precursor prec;
prec.setIntensity(4242.0f);

MSSpectrum ms2_spec;
ms2_spec.setMSLevel(2);

// MS2 whose precursor WAS traced -> branch A, payload from mass_features[0]
ms2_spec.setRT(66.0);
DeconvolvedSpectrum ms2_traced(2);
ms2_traced.setOriginalSpectrum(ms2_spec);
ms2_traced.setPrecursorPeakGroup(makePrecursorPG(5000.0, 5, kMs1Scan, 1));
ms2_traced.setPrecursor(prec);

// MS2 whose precursor was NEVER traced -> gets a minted ID -> branch B, payload from the pg itself
ms2_spec.setRT(72.0);
DeconvolvedSpectrum ms2_minted(3);
ms2_minted.setOriginalSpectrum(ms2_spec);
ms2_minted.setPrecursorPeakGroup(makePrecursorPG(2450.0, 3, kMs1Scan, 0));
ms2_minted.setPrecursor(prec);

std::vector<DeconvolvedSpectrum> spectra;
spectra.push_back(ms1);
spectra.push_back(ms2_traced);
spectra.push_back(ms2_minted);

std::map<int, double> scan_rt_map;
scan_rt_map[kMs1Scan] = kMs1Rt;

/// ---- drive the writer, MS1 first (that pass mints) ------------------------------------------
String ms1_out;
NEW_TMP_FILE(ms1_out)
{
  std::fstream fs;
  fs.open(ms1_out.c_str(), std::fstream::out);
  FLASHDeconvFeatureFile::writeTopFDFeatureHeader(fs, 1);
  FLASHDeconvFeatureFile::writeTopFDFeatures(spectra, mass_features, scan_rt_map, "in.mzML", fs, 1);
  fs.close();
}

String ms2_out;
NEW_TMP_FILE(ms2_out)
{
  std::fstream fs;
  fs.open(ms2_out.c_str(), std::fstream::out);
  FLASHDeconvFeatureFile::writeTopFDFeatureHeader(fs, 2);
  FLASHDeconvFeatureFile::writeTopFDFeatures(spectra, mass_features, scan_rt_map, "in.mzML", fs, 2);
  fs.close();
}

std::vector<String> ms1_lines = readLines(ms1_out);
std::vector<String> ms2_lines = readLines(ms2_out);

START_SECTION((static void writeTopFDFeatureHeader(std::fstream& fs, uint ms_level)))
{
  ABORT_IF(ms2_lines.empty())
  std::vector<String> h2 = splitTabs(ms2_lines[0]);
  // TopPIC 1.8.x reads strs[16]; 16 columns is the crash.
  TEST_EQUAL(h2.size(), 17)
  ABORT_IF(h2.size() != 17)
  TEST_STRING_EQUAL(h2[12], "Precursor_neutral_monoisotopic_mass")
  TEST_STRING_EQUAL(h2[13], "Precursor_monoisotopic_mz")
  TEST_STRING_EQUAL(h2[14], "Precursor_average_mz")
  TEST_STRING_EQUAL(h2[15], "Precursor_charge")
  TEST_STRING_EQUAL(h2[16], "Precursor_intensity")

  // the _ms1.feature layout is unchanged at 18
  ABORT_IF(ms1_lines.empty())
  TEST_EQUAL(splitTabs(ms1_lines[0]).size(), 18)
}
END_SECTION

START_SECTION((static void writeTopFDFeatures(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const std::vector<FLASHHelperClasses::MassFeature>& mass_features, const std::map<int, double>& scan_rt_map, const String& file_name, std::fstream& fs, uint ms_level)))
{
  // 2 traced MS1 targets + 1 minted synthetic row; f3 is MSn and is never written here
  TEST_EQUAL(ms1_lines.size(), 1 + 3)
  // one row per MS2 spectrum with a precursor
  TEST_EQUAL(ms2_lines.size(), 1 + 2)

  for (Size i = 1; i < ms2_lines.size(); ++i)
  {
    TEST_EQUAL(splitTabs(ms2_lines[i]).size(), 17)
  }
  for (Size i = 1; i < ms1_lines.size(); ++i)
  {
    TEST_EQUAL(splitTabs(ms1_lines[i]).size(), 18)
  }
}
END_SECTION

START_SECTION([EXTRA] a minted precursor reports its OWN payload rather than mass_features[id - 1])
{
  ABORT_IF(ms2_lines.size() != 3)
  std::vector<String> traced = splitTabs(ms2_lines[1]);
  std::vector<String> minted = splitTabs(ms2_lines[2]);
  ABORT_IF(traced.size() != 17 || minted.size() != 17)

  // Minting starts above EVERY feature (3), not above the MS1-target max (2).
  TEST_EQUAL(traced[6].toInt(), 1)
  TEST_EQUAL(minted[6].toInt(), 4)

  // The payload must come from the precursor peak group itself. The old writer took branch A and
  // reported f3 here: mass 9999.0 and apex time 3010/60 = 50.166667 min.
  TEST_REAL_SIMILAR(minted[12].toDouble(), 2450.0)
  TEST_REAL_SIMILAR(minted[11].toDouble(), kMs1Rt / 60.0)
  TEST_EQUAL(minted[15].toInt(), 3)

  // The traced row still resolves through mass_features[0].
  TEST_REAL_SIMILAR(traced[12].toDouble(), 5000.0)
  TEST_REAL_SIMILAR(traced[11].toDouble(), 110.0 / 60.0)
  TEST_EQUAL(traced[15].toInt(), 5)
}
END_SECTION

START_SECTION([EXTRA] precursor mass/mz identity and isotope ordering)
{
  ABORT_IF(ms2_lines.size() != 3)
  for (Size i = 1; i < ms2_lines.size(); ++i)
  {
    std::vector<String> f = splitTabs(ms2_lines[i]);
    ABORT_IF(f.size() != 17)
    const int z = f[15].toInt();
    TEST_NOT_EQUAL(z, 0)
    // TopPIC's own invariant: col 14 is col 13 expressed as an m/z at col 16's charge.
    TEST_REAL_SIMILAR(f[13].toDouble(), f[12].toDouble() / z + Constants::PROTON_MASS_U)
  }

  // Traced row: Precursor_average_mz is the feature's average mass as an m/z, protons included,
  // and therefore sits ABOVE the monoisotopic m/z. Without the proton it sat below it on every row.
  std::vector<String> traced = splitTabs(ms2_lines[1]);
  TEST_REAL_SIMILAR(traced[14].toDouble(), (5000.0 + 0.6) / 5 + Constants::PROTON_MASS_U)
  TEST_EQUAL(traced[14].toDouble() > traced[13].toDouble(), true)

  // _ms1.feature's Rep_average_mz is the same quantity and moves with it.
  ABORT_IF(ms1_lines.size() < 2)
  std::vector<String> f1_row = splitTabs(ms1_lines[1]);
  ABORT_IF(f1_row.size() != 18)
  TEST_REAL_SIMILAR(f1_row[15].toDouble(), (5000.0 + 0.6) / 5 + Constants::PROTON_MASS_U)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
