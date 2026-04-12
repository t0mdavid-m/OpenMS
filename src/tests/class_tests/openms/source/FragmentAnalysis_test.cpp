// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>

#include <fstream>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  struct ScanData
  {
    std::vector<double> mzs;
    std::vector<double> ints;
    double rt;
    std::string scan_id;
  };

  std::vector<ScanData> loadTsvScans(const std::string& path)
  {
    std::vector<ScanData> scans;
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line))
    {
      if (line.substr(0, 4) == "Spec")
      {
        scans.emplace_back();
        auto tab = line.find('\t');
        scans.back().scan_id = line.substr(10, tab - 10);
        scans.back().rt = std::stod(line.substr(tab + 1));
      }
      else if (!scans.empty())
      {
        auto tab = line.find('\t');
        if (tab != std::string::npos)
        {
          scans.back().mzs.push_back(std::stod(line.substr(0, tab)));
          scans.back().ints.push_back(std::stod(line.substr(tab + 1)));
        }
      }
    }
    return scans;
  }

  // Horse heart cytochrome c sequence (104 AA)
  const char* cytochrome_c_seq = "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE";

  const std::string ms2_cytc_path = "../../FlashIDA/test-data/spectra/ms2_cytc_scan149.txt";

  // JSON config for fragment analysis tests
  const char* fragment_test_config = R"({
    "deconvolution": {
      "score_threshold": 0.0,
      "tqscore_threshold": 0.9,
      "min_charge": 4,
      "max_charge": 50,
      "min_mass": 500,
      "max_mass": 50000,
      "tol": [10, 10]
    },
    "precursor_selection": {
      "RT_window": 180,
      "target_mode": 0,
      "IDScore": false,
      "AllCharges": false,
      "HCDEnergy": 29,
      "strict_inclusion": false,
      "tie_threshold": 0.1
    },
    "tagging": {
      "min_tag_length": 3,
      "max_tag_length": 8,
      "max_ptm_count": 3,
      "max_flanking_mass_diff": 50000
    },
    "quantification": {
      "enabled": false,
      "reporter_mz_tol": 0.002,
      "fold_change_threshold": 1.4
    },
    "faims": {
      "cv_values": [-50],
      "max_cv_skip": 0,
      "cv_precursor_threshold": 15
    },
    "ms_settings": {
      "ms1": {
        "analyzer": "Orbitrap",
        "first_mass": 500,
        "last_mass": 2000,
        "resolution": 120000,
        "agc_target": 800000,
        "max_it": 246
      },
      "ms2": [
        {
          "analyzer": "Orbitrap",
          "activation": "HCD",
          "collision_energy": 29,
          "resolution": 120000
        }
      ]
    },
    "scheduling": {
      "cycle_time": { "enabled": false, "value_ms": 60000 },
      "scan_timeout": { "enabled": false, "value_ms": 30000 }
    },
    "files": {
      "target_logs": [],
      "fasta": "",
      "inclusion_list": "",
      "ptm_list": ""
    },
    "ms3": {
      "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE"
    },
    "conditional_ms2": false,
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_targets": 3 },
      "ms2": { "selection": "intensity", "max_targets": 3 },
      "ms3": { "selection": "intensity", "max_targets": 3 }
    }
  })";

} // anonymous namespace


START_TEST(FragmentAnalysis, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(getTopFragmentMatches_cytochrome_c)
{
  auto scans = loadTsvScans(ms2_cytc_path);
  ABORT_IF(scans.empty())

  Config cfg{std::string(fragment_test_config)};
  Deconvolution deconv(cfg);

  // Deconvolve scan 149 (MS2, no precursor info)
  int pg_count = deconv.deconvolveMSn(scans[0].mzs.data(), scans[0].ints.data(),
                                       (int)scans[0].mzs.size(), scans[0].rt, 0.0, 0);
  TEST_EQUAL(pg_count > 0, true)

  FragmentAnalysis fragments(cfg);

  double masses[20], qscores[20], wstarts[20], wends[20];
  int charges[20], frag_indices[20];
  char ion_types[20];

  int found = fragments.getTopFragmentMatches(cytochrome_c_seq, 20,
                                               masses, qscores, charges,
                                               wstarts, wends, ion_types,
                                               frag_indices, deconv.storedMS2());
  TEST_EQUAL(found > 0, true)

  for (int i = 0; i < found; ++i)
  {
    TEST_EQUAL(masses[i] > 0, true)
    TEST_EQUAL(charges[i] > 0, true)
    TEST_EQUAL(frag_indices[i] > 0, true)
    TEST_EQUAL(ion_types[i] == 'b' || ion_types[i] == 'y', true)
  }
}
END_SECTION

START_SECTION(getTerminalFragmentIons_cytochrome_c)
{
  auto scans = loadTsvScans(ms2_cytc_path);
  ABORT_IF(scans.empty())

  Config cfg{std::string(fragment_test_config)};
  Deconvolution deconv(cfg);

  deconv.deconvolveMSn(scans[0].mzs.data(), scans[0].ints.data(),
                        (int)scans[0].mzs.size(), scans[0].rt, 0.0, 0);

  FragmentAnalysis fragments(cfg);

  double masses[20], qscores[20], wstarts[20], wends[20];
  int charges[20], frag_indices[20];
  char ion_types[20];

  int found = fragments.getTerminalFragmentIons(cytochrome_c_seq, 10,
                                                 masses, qscores, charges,
                                                 wstarts, wends, ion_types,
                                                 frag_indices, deconv.storedMS2());
  TEST_EQUAL(found > 0, true)

  for (int i = 0; i < found; ++i)
  {
    TEST_EQUAL(masses[i] > 0, true)
    TEST_EQUAL(ion_types[i] == 'b' || ion_types[i] == 'y', true)
  }
}
END_SECTION

START_SECTION(getTopFragmentMatches_empty_spectrum_returns_zero)
{
  Config cfg{std::string(fragment_test_config)};
  FragmentAnalysis fragments(cfg);

  DeconvolvedSpectrum ds(1);

  double masses[20], qscores[20], wstarts[20], wends[20];
  int charges[20], frag_indices[20];
  char ion_types[20];

  int found = fragments.getTopFragmentMatches(cytochrome_c_seq, 20,
                                               masses, qscores, charges,
                                               wstarts, wends, ion_types,
                                               frag_indices, ds);
  TEST_EQUAL(found, 0)
}
END_SECTION

START_SECTION(getBestMS2Masses_returns_deconvolved_peaks)
{
  auto scans = loadTsvScans(ms2_cytc_path);
  ABORT_IF(scans.empty())

  Config cfg{std::string(fragment_test_config)};
  Deconvolution deconv(cfg);

  deconv.deconvolveMSn(scans[0].mzs.data(), scans[0].ints.data(),
                        (int)scans[0].mzs.size(), scans[0].rt, 0.0, 0);

  FragmentAnalysis fragments(cfg);

  double masses[10], qscores[10], wstarts[10], wends[10];
  int charges[10];

  int found = fragments.getBestMS2Masses(10, masses, qscores, charges,
                                          wstarts, wends, deconv.storedMS2());
  TEST_EQUAL(found > 0, true)

  for (int i = 0; i < found; ++i)
  {
    TEST_EQUAL(masses[i] >= 500.0, true)
    TEST_EQUAL(wstarts[i] < wends[i], true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
