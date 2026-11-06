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
// annotation. FLASHTnT's parameters are ProteomeSearch::defaultParameters(): its defaults with exactly two diffs,
// one_proteoform_per_feature false and ex:skip_precursor_inference true (d.13).
//
// Three things are pinned to each other:
//   - FLASHTnTAlgorithm::run on the whole slice, nothing prepared, to tnt_reference.tsv: the same run captured
//     BEFORE the prepared-database edit (so the edit changed no PrSM of offline FLASHTnT);
//   - ProteomeSearch's per-MS2 search (a prepared database, run() on one spectrum) to that whole run;
//   - the per-MS2 search on 4 threads, each with its own Worker, to the serial one.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteomeSearch.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ProteinHit.h>

#include "FLASHIda_TestHelpers.h"  // buildJsonWithLogDir

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;

namespace
{
  const std::string kSlice = "../../test-data-private/ecoli_parity/ms2_slice.mzML";
  const std::string kReference = "../../test-data-private/ecoli_parity/tnt_reference.tsv";

  // The converted MS2s (FLASHTnT's input, in slice order), and per scan the engine's deconvolution and the raw spectrum.
  struct ParityInput
  {
    MSExperiment deconvolved;
    std::vector<int> scans;  // the scans of `deconvolved`, in order
    std::map<int, DeconvolvedSpectrum> ds;
    std::map<int, MSSpectrum> raw;
    double ms2_tol = 0;
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
    ParityInput in;
    in.ms2_tol = config.toleranceList()[1];
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
      in.ds.emplace(scan, ds);
      ds.sort();  // toSpectrum sorts the peaks but not its per-mass lists
      MSSpectrum out = ds.toSpectrum(0, in.ms2_tol, false);
      out.setNativeID(spec.getNativeID());
      out.setMSLevel(2);
      out.setRT(spec.getRT());
      in.deconvolved.addSpectrum(out);
      in.scans.push_back(scan);
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

  std::vector<FASTAFile::FASTAEntry> reference()
  {
    std::vector<FASTAFile::FASTAEntry> fasta;
    FASTAFile().load(OPENMS_GET_TEST_DATA_PATH("uniprot-mg1655-filtered-reviewed_yes.fasta"), fasta);
    return fasta;
  }

  // One PrSM as a row: every value a parity comparison reads, doubles at full precision.
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

  std::vector<std::string> split(const std::string& s, char d)
  {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream is(s);
    while (std::getline(is, cur, d)) out.push_back(cur);
    if (! s.empty() && s.back() == d) out.emplace_back();
    return out;
  }

  bool asNumber(const std::string& s, double& x)
  {
    if (s.empty()) return false;
    char* end = nullptr;
    x = std::strtod(s.c_str(), &end);
    return end == s.c_str() + s.size();
  }

  // equal as text, or as numbers (each comma-separated piece) to relative 1e-9 (ADR-0061 d.4)
  bool sameCell(const std::string& a, const std::string& b)
  {
    if (a == b) return true;
    const auto pa = split(a, ','), pb = split(b, ',');
    if (pa.size() != pb.size()) return false;
    for (Size i = 0; i < pa.size(); i++)
    {
      double x, y;
      if (! asNumber(pa[i], x) || ! asNumber(pb[i], y)) return false;
      if (! (std::fabs(x - y) <= 1e-9 * std::max(std::fabs(x), std::fabs(y)))) return false;
    }
    return true;
  }

  bool sameRow(const std::string& a, const std::string& b, const std::string& what)
  {
    const auto ca = split(a, '\t'), cb = split(b, '\t');
    bool same = ca.size() == cb.size();
    for (Size i = 0; same && i < ca.size(); i++) same = sameCell(ca[i], cb[i]);
    if (! same) std::cout << "    [PARITY-MISMATCH] " << what << "\n      " << a << "\n      " << b << std::endl;
    return same;
  }

  // FLASHTnT's rank of a spectrum's PrSMs (FLASHTnTAlgorithm.cpp, score then coverage); the earlier one on a tie
  const ProteinHit* bestOf(const std::vector<ProteinHit>& hits)
  {
    const ProteinHit* best = nullptr;
    for (const auto& h : hits)
    {
      if (best == nullptr || h.getScore() > best->getScore() || (h.getScore() == best->getScore() && h.getCoverage() > best->getCoverage()))
        best = &h;
    }
    return best;
  }

  std::map<int, std::vector<ProteinHit>> byScan(const std::vector<ProteinHit>& hits)
  {
    std::map<int, std::vector<ProteinHit>> out;
    for (const auto& h : hits) out[(int)h.getMetaValue("Scan")].push_back(h);
    return out;
  }

  void setTeamSize()
  {
#ifdef _OPENMP
    omp_set_num_threads(4);  // as SearchPool's workers do: a new thread does not inherit the team size
#endif
  }
} // namespace

START_TEST(FLASHIda_ProteomeSearch, "$Id$")

const ParityInput in = buildInput();
const auto fasta = reference();
const Param params = ProteomeSearch::defaultParameters();
std::map<int, std::vector<ProteinHit>> whole;   // the whole run's PrSMs by scan
std::map<int, std::vector<ProteinHit>> serial;  // the serial per-MS2 search's

START_SECTION(the_whole_run_equals_the_capture_from_before_the_prepared_database_edit)
{
  FLASHTnTAlgorithm tnt;
  tnt.setParameters(params);
  tnt.setRawSpectrumLoader(loaderOf(in));
  tnt.run(in.deconvolved, fasta);
  whole = byScan(tnt.getProteoforms());

  std::ifstream ref(kReference);
  TEST_EQUAL(ref.good(), true)
  std::vector<std::string> rows;
  for (std::string line; std::getline(ref, line);)
  {
    if (! line.empty() && line.back() == '\r') line.pop_back();  // the parent checkout may carry CRLF
    rows.push_back(line);
  }
  TEST_EQUAL(rows.empty(), false)
  TEST_EQUAL(rows.empty() ? std::string() : rows[0], prsmHeader())
  const auto& hits = tnt.getProteoforms();
  TEST_EQUAL(hits.size() + 1, rows.size())
  Size same = 0;
  for (Size i = 0; i < hits.size() && i + 1 < rows.size(); i++) same += sameRow(prsmRow(hits[i]), rows[i + 1], "whole run vs capture, row " + std::to_string(i + 1));
  TEST_EQUAL(same, hits.size())
  std::cout << "[PARITY] ms2=" << in.ms2_count << " deconvolved=" << in.deconvolved.size() << " prsms=" << hits.size()
            << " scans_with_prsm=" << whole.size() << " fasta=" << fasta.size() << std::endl;
  TEST_EQUAL(in.ms2_count, 200)
}
END_SECTION

START_SECTION(the_per_ms2_search_equals_the_whole_run)
{
  ProteomeSearch search(fasta, params);
  TEST_EQUAL(search.decoysDropped(), 0)
  ProteomeSearch::Worker worker(search);
  std::vector<double> ms;
  for (int scan : in.scans)
  {
    const auto t0 = std::chrono::steady_clock::now();
    auto hits = worker.search(in.ds.at(scan), in.ms2_tol, in.raw.at(scan), scan);
    ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
    if (! hits.empty()) serial[scan] = std::move(hits);
  }
  std::vector<int> whole_scans, serial_scans;
  for (const auto& [s, h] : whole) whole_scans.push_back(s);
  for (const auto& [s, h] : serial) serial_scans.push_back(s);
  TEST_EQUAL(serial_scans == whole_scans, true)  // the same spectra without a PrSM
  Size same = 0;
  for (const auto& [scan, hits] : whole)
  {
    const auto it = serial.find(scan);
    if (it == serial.end()) continue;
    same += sameRow(prsmRow(*bestOf(it->second)), prsmRow(*bestOf(hits)), "best PrSM, per MS2 vs whole run, scan " + std::to_string(scan));
  }
  TEST_EQUAL(same, whole.size())
  TEST_EQUAL(whole.size() >= 20, true)  // an empty pair cannot pass
  std::sort(ms.begin(), ms.end());
  std::cout << "[PARITY] per-MS2 search: " << ms.size() << " spectra, " << serial.size() << " with a PrSM, median "
            << ms[ms.size() / 2] << " ms, p95 " << ms[(ms.size() * 95) / 100] << " ms, max " << ms.back() << " ms" << std::endl;
}
END_SECTION

START_SECTION(the_per_ms2_search_on_4_threads_equals_the_serial_search)
{
  ProteomeSearch search(fasta, params);
  std::vector<std::map<int, std::vector<ProteinHit>>> found(4);
  std::vector<std::thread> threads;
  for (int t = 0; t < 4; t++)
  {
    threads.emplace_back([&, t]() {
      setTeamSize();
      ProteomeSearch::Worker worker(search);
      for (Size k = t; k < in.scans.size(); k += 4)
      {
        const int scan = in.scans[k];
        auto hits = worker.search(in.ds.at(scan), in.ms2_tol, in.raw.at(scan), scan);
        if (! hits.empty()) found[t][scan] = std::move(hits);
      }
    });
  }
  for (auto& th : threads) th.join();
  std::map<int, std::vector<ProteinHit>> parallel;
  for (auto& f : found) parallel.insert(f.begin(), f.end());
  TEST_EQUAL(parallel.size(), serial.size())
  Size same = 0, compared = 0;
  for (const auto& [scan, hits] : serial)
  {
    const auto it = parallel.find(scan);
    if (it == parallel.end() || it->second.size() != hits.size())
    {
      std::cout << "    [PARITY-MISMATCH] scan " << scan << ": " << hits.size() << " PrSMs serially, "
                << (it == parallel.end() ? 0 : it->second.size()) << " on 4 threads" << std::endl;
      continue;
    }
    for (Size i = 0; i < hits.size(); i++)
    {
      compared++;
      same += sameRow(prsmRow(it->second[i]), prsmRow(hits[i]), "4 threads vs serial, scan " + std::to_string(scan));
    }
  }
  Size serial_prsms = 0;
  for (const auto& [scan, hits] : serial) serial_prsms += hits.size();
  TEST_EQUAL(compared, serial_prsms)
  TEST_EQUAL(same, serial_prsms)
}
END_SECTION

START_SECTION(a_decoy_entry_is_dropped_at_load)
{
  auto with_decoy = fasta;
  with_decoy.emplace_back("DECOY_sp|P0A7V0|RS2_ECOLI", "a decoy", "MATVSMRDMLKAGVHFGHQTRYWNPKMKPFIFGARNKVHIINLEKTVPMFNEALAELNKIASRKGKILFVGTKRAASEAVKDAALSCDQFFVNHRWLGGMLTNWKTVRQSIKRLKDLETQSQDGTFDKLTKKEALMRTRELEKLENSLGGIKDMGGLPDALFVIDADHEHIAIKEANNLGIPVFAIVDTNSDPDGVDFVIPGNDDAIRAVTLYLGAVAATVREGRSQDLASQAEESFVEAE");
  ProteomeSearch search(with_decoy, params);
  TEST_EQUAL(search.decoysDropped(), 1)
  TEST_EQUAL(search.reference().size(), fasta.size())
  ProteomeSearch::Worker worker(search);
  Size same = 0, compared = 0;
  for (Size k = 0; k < in.scans.size() && k < 20; k++)  // results unchanged on the first 20 spectra
  {
    const int scan = in.scans[k];
    const auto hits = worker.search(in.ds.at(scan), in.ms2_tol, in.raw.at(scan), scan);
    const auto it = serial.find(scan);
    const Size expected = it == serial.end() ? 0 : it->second.size();
    TEST_EQUAL(hits.size(), expected)
    for (Size i = 0; i < hits.size() && i < expected; i++)
    {
      compared++;
      same += sameRow(prsmRow(hits[i]), prsmRow(it->second[i]), "with a decoy entry vs without, scan " + std::to_string(scan));
    }
  }
  TEST_EQUAL(same, compared)
}
END_SECTION

END_TEST
