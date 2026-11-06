// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

// ADR-0061 d.4: the real-time search is FLASHTnT's own. The fixture is 200 consecutive MS2s of a lab E. coli LC-MS
// run (test-data-private/ecoli_parity/, in the PRIVATE parent repo) and OpenMS's public K-12 FASTA. Every MS2 goes
// through the engine's MS2 deconvolution (FLASHIda/Deconvolution::deconvolveMSn) and is handed to FLASHTnT the
// way its TOPP tool reads FLASHDeconv's output (DeconvolvedSpectrum::toSpectrum), with its raw peaks for the
// annotation. FLASHTnT's parameters are its defaults with exactly two diffs: one_proteoform_per_feature false and
// ex:skip_precursor_inference true (d.13).
//
// This first version holds only the capture: FLASHTnTAlgorithm::run on the whole slice BEFORE the prepared-database
// edit, written to ecoli_parity_capture.tsv (cwd OpenMS/build), which CI keeps as a local artifact. It is promoted
// to test-data-private/ecoli_parity/tnt_reference.tsv and pins run() across the edit.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ProteinHit.h>

#include "FLASHIda_TestHelpers.h"  // buildJsonWithLogDir

#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  const std::string kSlice = "../../test-data-private/ecoli_parity/ms2_slice.mzML";

  // The converted MS2s (FLASHTnT's input, in slice order) and every MS2's raw spectrum by scan (the raw loader's source).
  struct ParityInput
  {
    MSExperiment deconvolved;
    std::map<int, MSSpectrum> raw;
    int ms2_count = 0;
  };

  int scanOf(const String& native_id)
  {
    const auto p = native_id.rfind("scan=");
    return p == String::npos ? -1 : std::stoi(native_id.substr(p + 5));
  }

  // The engine's MS2 deconvolution of every MS2: the deconvolution section of buildJsonWithLogDir, which is
  // method_ms3_cytc_real.json's (charges 4-50, masses 500-50000, 10 ppm at every level; its tqscore_threshold gates
  // MS1 exclusion only). The precursor mass is the selected ion's (m/z - proton) x charge. A spectrum without a mass
  // is left out of both sides: toSpectrum reads its first mass.
  ParityInput buildInput()
  {
    MSExperiment slice;
    MzMLFile().load(kSlice, slice);
    Config config(buildJsonWithLogDir(""));
    Deconvolution deconv(config, config.toleranceList());
    const double ms2_tol = config.toleranceList()[1];
    ParityInput in;
    for (const auto& spec : slice)
    {
      if (spec.getMSLevel() != 2) continue;
      in.ms2_count++;
      const int scan = scanOf(spec.getNativeID());
      in.raw[scan] = spec;
      std::vector<double> mzs, ints;
      for (const auto& p : spec)
      {
        mzs.push_back(p.getMZ());
        ints.push_back(p.getIntensity());
      }
      double mass = 0;
      int z = 0;
      if (! spec.getPrecursors().empty() && spec.getPrecursors()[0].getCharge() > 0)
      {
        z = spec.getPrecursors()[0].getCharge();
        mass = (spec.getPrecursors()[0].getMZ() - Constants::PROTON_MASS_U) * z;
      }
      deconv.deconvolveMSn(mzs.data(), ints.data(), (int)mzs.size(), spec.getRT(), mass, z);
      DeconvolvedSpectrum ds = deconv.storedMS2();
      if (ds.empty()) continue;
      ds.sort();  // toSpectrum sorts the peaks but not its per-mass lists
      MSSpectrum out = ds.toSpectrum(0, ms2_tol, false);
      out.setNativeID(spec.getNativeID());
      out.setMSLevel(2);
      out.setRT(spec.getRT());
      in.deconvolved.addSpectrum(out);
    }
    return in;
  }

  FLASHTnTAlgorithm::RawSpectrumLoader loaderOf(const ParityInput& in)
  {
    return [&in](const std::set<int>& scans, std::map<int, MSSpectrum>& out) {
      for (int s : scans)
      {
        const auto it = in.raw.find(s);
        if (it != in.raw.end()) out[s] = it->second;
      }
    };
  }

  Param parityParams()
  {
    Param p = FLASHTnTAlgorithm().getDefaults();
    p.setValue("one_proteoform_per_feature", "false");
    p.setValue("ex:skip_precursor_inference", "true");
    return p;
  }

  std::vector<FASTAFile::FASTAEntry> reference()
  {
    std::vector<FASTAFile::FASTAEntry> fasta;
    FASTAFile().load(OPENMS_GET_TEST_DATA_PATH("uniprot-mg1655-filtered-reviewed_yes.fasta"), fasta);
    return fasta;
  }

  // One PrSM as a capture row: every value a parity comparison reads, doubles at full precision.
  const std::vector<std::string> kPrsmMeta = {"Scan", "StartPosition", "EndPosition", "Modifications", "ModificationStarts",
                                              "ModificationEnds", "ExtensionScore", "RawScore"};

  std::string cellOf(const ProteinHit& h, const std::string& key)
  {
    if (! h.metaValueExists(key)) return "NA";
    const DataValue& v = h.getMetaValue(key);
    std::ostringstream o;
    o << std::setprecision(17);
    switch (v.valueType())
    {
      case DataValue::DOUBLE_VALUE: o << (double)v; break;
      case DataValue::DOUBLE_LIST: {
        const auto l = v.toDoubleList();
        for (Size i = 0; i < l.size(); i++) o << (i ? "," : "") << l[i];
        break;
      }
      case DataValue::INT_LIST: {
        const auto l = v.toIntList();
        for (Size i = 0; i < l.size(); i++) o << (i ? "," : "") << l[i];
        break;
      }
      default: o << v.toString(); break;
    }
    return o.str();
  }

  std::string prsmHeader()
  {
    std::string s = "accession\tscore";
    for (const auto& k : kPrsmMeta) s += "\t" + k;
    return s;
  }

  std::string prsmRow(const ProteinHit& h)
  {
    std::ostringstream o;
    o << std::setprecision(17) << h.getAccession() << "\t" << h.getScore();
    for (const auto& k : kPrsmMeta) o << "\t" << cellOf(h, k);
    return o.str();
  }
} // namespace

START_TEST(FLASHIda_ProteomeSearch, "$Id$")

START_SECTION(capture_the_whole_run_before_the_prepared_database_edit)
{
  const ParityInput in = buildInput();
  const auto fasta = reference();
  FLASHTnTAlgorithm tnt;
  tnt.setParameters(parityParams());
  tnt.setRawSpectrumLoader(loaderOf(in));
  tnt.run(in.deconvolved, fasta);

  std::ofstream out("ecoli_parity_capture.tsv");
  out << prsmHeader() << "\n";
  std::set<int> identified;
  for (const auto& h : tnt.getProteoforms())
  {
    out << prsmRow(h) << "\n";
    identified.insert((int)h.getMetaValue("Scan"));
  }
  std::cout << "[PARITY] ms2=" << in.ms2_count << " deconvolved=" << in.deconvolved.size() << " prsms=" << tnt.getProteoforms().size()
            << " scans_with_prsm=" << identified.size() << " fasta=" << fasta.size() << std::endl;
  TEST_EQUAL(in.ms2_count, 200)
  TEST_EQUAL(identified.size() >= 20, true)
}
END_SECTION

END_TEST
