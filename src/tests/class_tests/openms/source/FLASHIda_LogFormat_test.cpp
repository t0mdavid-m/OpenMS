// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// ADR-0079: the W gate reads its inputs "as logged" -- through LogFormat.h, the one formatter the logger writes with. These
// sections pin that the cells the logger writes ARE the LogFormat strings, so the logger and the gate cannot drift apart,
// and that composeMS3 + assignMS3Id is buildMS3 byte for byte (a command is composed before the gate decides and issued
// only when sent).

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/LogFormat.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include "FLASHIda_TestHelpers.h"

#include <cmath>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  // An MS2 context with awkward values: more digits than any log prints, so a second formatting would show.
  ScanCommand awkwardMs2()
  {
    ScanCommand c{};
    c.scan_id = 77;
    c.msn_level = 2;
    c.num_stages = 1;
    c.stages[0].precursor_mz = 1234.5678912345;
    c.stages[0].isolation_width = 1.23456789;
    c.stages[0].charge_state = 11;
    c.stages[0].collision_energy = 29.0;
    std::strncpy(c.stages[0].activation_type, "HCD", sizeof(c.stages[0].activation_type) - 1);
    c.mono_mass = 13569.123456789;
    c.qscore = 0.123456789;
    c.charge_cos = 0.987654321;
    c.charge_snr = 12.3456789;
    c.iso_cos = 0.99912345;
    c.snr = 3.14159265;
    c.charge_score = 0.55555555;
    c.ppm_error = -1.23456789;
    c.precursor_intensity = 1234567.891;
    c.peakgroup_intensity = 7654321.987;
    std::strncpy(c.scan_description, "!!#R13.6k@11", sizeof(c.scan_description) - 1);
    return c;
  }

  FragmentAnalysis::FragmentScores awkwardFragment()
  {
    FragmentAnalysis::FragmentScores s;
    s.qscore = 0.876543219;
    s.mono_mass = 2010.13098765;
    s.charge_cos = 0.950716789;
    s.charge_snr = 35.7401234;
    s.iso_cos = 0.984057891;
    s.snr = 99.9142345;
    s.charge_score = 1.0;
    s.ppm_error = 0.756906789;
    s.precursor_intensity = 7269.0234;
    s.peakgroup_intensity = 8123.4567;
    return s;
  }
}

START_TEST(FLASHIda_LogFormat, "$Id$")

START_SECTION(formatters_are_the_stream_defaults)
{
  TEST_STRING_EQUAL(formatSig6(1234567.891), "1.23457e+06")
  TEST_STRING_EQUAL(formatSig6(0.123456789), "0.123457")
  TEST_STRING_EQUAL(formatSig6(-1.0), "-1")
  TEST_STRING_EQUAL(formatFixed4(13569.123456789), "13569.1235")
  TEST_STRING_EQUAL(formatFixed4(-1.0), "-1.0000")
  TEST_EQUAL(parseLogged("1.23457e+06"), 1234570.0)
  TEST_EQUAL(std::isnan(parseLogged("")), true)
  TEST_EQUAL(std::isnan(parseLogged("x")), true)
}
END_SECTION

START_SECTION(compose_then_assign_is_build)
{
  const std::string dir = freshLogDir("lfmt_compose");
  Config cfg{buildJsonWithLogDir(dir, true)};
  ScanCommandQueue built_q(cfg), composed_q(cfg);
  const ScanCommand ms2 = awkwardMs2();
  const FragmentAnalysis::FragmentScores fs = awkwardFragment();
  const ScanConfig& ms3_cfg = cfg.level(3).scans.empty() ? cfg.level(2).scans[0] : cfg.level(3).scans[0];

  const ScanCommand built = built_q.buildMS3(ms2, ms3_cfg, 671.05, 3, 1.5, ms2.scan_id, 'y', 18, 1, fs);
  ScanCommandQueue::ComposedMS3 composed = composed_q.composeMS3(ms2, ms3_cfg, 671.05, 3, 1.5, ms2.scan_id, 'y', 18, 1, fs);
  TEST_EQUAL(composed.cmd.scan_id, 0)   // not issued yet
  const int id = composed_q.assignMS3Id(composed);
  TEST_EQUAL(id, built.scan_id)
  // The whole 2048-byte command, field for field (both start from a zero-initialised ScanCommand).
  TEST_EQUAL(std::memcmp(&composed.cmd, &built, sizeof(ScanCommand)), 0)
  TEST_STRING_EQUAL(std::string(composed.cmd.scan_description), std::string(built.scan_description))
}
END_SECTION

START_SECTION(formatters_equal_logged_cells)
{
  const std::string dir = freshLogDir("lfmt_cells");
  Config cfg{buildJsonWithLogDir(dir, true)};
  ScanCommandQueue q(cfg);
  const ScanCommand ms2 = awkwardMs2();
  const ScanConfig& ms3_cfg = cfg.level(3).scans.empty() ? cfg.level(2).scans[0] : cfg.level(3).scans[0];
  const ScanCommand ms3 = q.buildMS3(ms2, ms3_cfg, 671.05123, 3, 1.5, ms2.scan_id, 'y', 18, 1, awkwardFragment());
  {
    IdaLogger logger(cfg);
    logger.writeScanCommandRow(ms2, 1, "");
    logger.writeScanCommandRow(ms3, 1, "PEPTIDEK");

    IdaLogger::ScanRowDescriptor r;
    r.tracking_id = "!!#";
    r.ms_level = 2;
    r.precursor_share = 0.0139028765;
    r.tic_coverage = 0.321568765f;
    logger.writeScanResultRow(r);

    IdaLogger::IdRowDescriptor i;
    i.tracking_id = "!!#";
    i.ms_level = 2;
    i.scan_mode = 'R';
    i.match.proteoform_sequence = "PEPTIDEK";
    FragmentAnalysis::ProteoformMatch::FragmentMatch f;
    f.ion_type = "b";
    f.ion_index = 3;
    f.observed_mass = 300.5;
    i.match.fragments.push_back(f);
    i.match.identification_score = 179.644187654;
    logger.writeIdentificationRow(i);
  }

  auto c = TSVFile::parse(dir + "/scan_commands.tsv");
  TEST_EQUAL(c.rows.size(), 2)
  const std::vector<ScanCommand> cmds{ms2, ms3};
  for (size_t k = 0; k < cmds.size() && k < c.rows.size(); ++k)
  {
    const ScanCommandCells e = scanCommandCells(cmds[k]);
    const auto& row = c.rows[k];
    TEST_STRING_EQUAL(cell(c, row, "mono_mass"), e.mono_mass)
    TEST_STRING_EQUAL(cell(c, row, "charge"), e.charge)
    TEST_STRING_EQUAL(cell(c, row, "precursor_mz"), e.precursor_mz)
    TEST_STRING_EQUAL(cell(c, row, "isolation_width"), e.isolation_width)
    TEST_STRING_EQUAL(cell(c, row, "qscore"), e.qscore)
    TEST_STRING_EQUAL(cell(c, row, "charge_cos"), e.charge_cos)
    TEST_STRING_EQUAL(cell(c, row, "charge_snr"), e.charge_snr)
    TEST_STRING_EQUAL(cell(c, row, "iso_cos"), e.iso_cos)
    TEST_STRING_EQUAL(cell(c, row, "snr"), e.snr)
    TEST_STRING_EQUAL(cell(c, row, "charge_score"), e.charge_score)
    TEST_STRING_EQUAL(cell(c, row, "collision_energy"), e.collision_energy)
    TEST_STRING_EQUAL(cell(c, row, "ppm_error"), e.ppm_error)
    TEST_STRING_EQUAL(cell(c, row, "precursor_intensity"), e.precursor_intensity)
    TEST_STRING_EQUAL(cell(c, row, "peakgroup_intensity"), e.peakgroup_intensity)
    TEST_STRING_EQUAL(cell(c, row, "first_mass"), e.first_mass)
    TEST_STRING_EQUAL(cell(c, row, "last_mass"), e.last_mass)
  }
  // The MS3 row's stage-1 cells are the fragment's values at the log's precision, not a re-rounding of them.
  if (c.rows.size() == 2)
  {
    TEST_STRING_EQUAL(cell(c, c.rows[1], "mono_mass"), formatFixed4(ms2.mono_mass) + ";" + formatFixed4(2010.13098765))
    TEST_STRING_EQUAL(cell(c, c.rows[1], "precursor_intensity"), formatSig6(ms2.precursor_intensity) + ";" + formatSig6(7269.0234))
  }

  auto r = TSVFile::parse(dir + "/scan_results.tsv");
  TEST_EQUAL(r.rows.size(), 1)
  if (r.rows.size() == 1)
  {
    TEST_STRING_EQUAL(cell(r, r.rows[0], "precursor_share"), formatSig6(0.0139028765))
    TEST_STRING_EQUAL(cell(r, r.rows[0], "tic_coverage"), formatSig6(0.321568765f))
  }

  auto i = TSVFile::parse(dir + "/identification.tsv");
  TEST_EQUAL(i.rows.size(), 1)
  if (i.rows.size() == 1) TEST_STRING_EQUAL(cell(i, i.rows[0], "identification_score"), formatFixed4(179.644187654))
}
END_SECTION

END_TEST
