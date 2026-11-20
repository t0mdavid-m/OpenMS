// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MassOffset.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// ADR-0081: one mass offset per run and MS level. The engine's state, matcher evidence and fallback vote against the
// analysis' test vectors (task v23 / v24 of the MS3-gate analysis), copied to the PRIVATE parent repo
// (test-data-private/mass_offset/: they hold masses of recorded E. coli runs, which the public forks must not carry):
//   - v24 (tools/v24/off24_ref.py, reading A): every scan's estimate used, restart and estimate after, exact;
//   - v23 (tools/v23/off_ref_impl.py): every MS3's matches, bonds, errors fed, the vote counter and every vote, with the
//     offset the vector used for that scan injected. The v23 estimates follow its retired median rule and are not
//     asserted; its MS2 sequence is not matched here either (the engine's MS2 match is FLASHIda's own PrSM match, not
//     the closest-unused rule) -- the MS2 state itself is pinned by the v24 MS2 sequence.

using namespace OpenMS;
using json = nlohmann::json;

namespace
{
  const std::string kDir = "../../test-data-private/mass_offset/";

  json loadJson(const std::string& path)
  {
    std::ifstream in(path);
    if (!in) return json();
    std::stringstream ss;
    ss << in.rdbuf();
    return json::parse(ss.str());
  }

  bool near(double a, double b, double tol = 1e-9) { return std::abs(a - b) <= tol * std::max(1.0, std::abs(b)); }

  bool sameList(const std::vector<double>& a, const json& b)
  {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
      if (!near(a[i], b[i].get<double>())) return false;
    return true;
  }

  // The analysis' result names, as the engine writes them.
  std::string engineResult(const std::string& v23)
  {
    static const std::map<std::string, std::string> m = {{"too few errors", "too_few_errors"},
                                                         {"too few scans", "too_few_scans"},
                                                         {"runner-up too close", "runner_up"},
                                                         {"confirms the estimate", "confirms"},
                                                         {"jump", "jump"}};
    const auto it = m.find(v23);
    return it == m.end() ? v23 : it->second;
  }
} // namespace

START_TEST(FLASHIda_MassOffset, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION(reading_a_estimate_is_the_posterior_mean)
{
  // Prior N(-25, 0.6^2); one error of -24 at sigma 0.5787868575652204: (-25/0.36 + -24/s^2) / (1/0.36 + 1/s^2)
  MassOffsetState st(-25.0, kMassOffsetSigmaMs3Ppm);
  TEST_REAL_SIMILAR(st.estimate(), -25.0)
  st.feed({-24.0});
  const double w0 = 1.0 / 0.36, w = 1.0 / (kMassOffsetSigmaMs3Ppm * kMassOffsetSigmaMs3Ppm);
  TEST_REAL_SIMILAR(st.estimate(), (-25.0 * w0 - 24.0 * w) / (w0 + w))
  TEST_EQUAL(st.count(), 1)
  st.restart(-185.0);   // a jump: the prior at the vote's centre, every observation dropped
  TEST_REAL_SIMILAR(st.estimate(), -185.0)
  TEST_EQUAL(st.count(), 0)
}
END_SECTION

START_SECTION(the_vectors_are_present)
{
  TEST_EQUAL(loadJson(kDir + "v24_offset24_test_vectors.json").is_object(), true)
  TEST_EQUAL(loadJson(kDir + "v23_engine_offset_test_vectors.json").is_object(), true)
}
END_SECTION

START_SECTION(v24_reading_a_reproduces_every_scan)
{
  const json d = loadJson(kDir + "v24_offset24_test_vectors.json");
  ABORT_IF(!d.is_object())
  TEST_EQUAL(d["vectors"].size(), 5)
  // The spreads the engine fixes are the ones the vectors were produced with.
  TEST_REAL_SIMILAR(d["sigma_ppm"]["ms3"].get<double>(), kMassOffsetSigmaMs3Ppm)
  TEST_REAL_SIMILAR(d["sigma_ppm"]["ms2"].get<double>(), kMassOffsetSigmaMs2Ppm)
  int scans = 0, mismatches = 0;
  for (const json& v : d["vectors"])
  {
    const double sigma = v["level"].get<int>() == 3 ? kMassOffsetSigmaMs3Ppm : kMassOffsetSigmaMs2Ppm;
    MassOffsetState st(v["start_ppm"].get<double>(), sigma);
    for (const json& sc : v["scans"])
    {
      ++scans;
      bool ok = near(st.estimate(), sc["estimate_used"].get<double>());
      std::vector<double> fed;
      for (const json& e : sc["errors_fed"]) fed.push_back(e.get<double>());
      st.feed(fed);
      if (!sc["restart_to"].is_null()) st.restart(sc["restart_to"].get<double>());
      ok = ok && near(st.estimate(), sc["estimate_after"].get<double>());
      if (!ok)
      {
        ++mismatches;
        std::cout << "v24 " << v["name"].get<std::string>() << " k=" << sc["k"] << " estimate=" << st.estimate()
                  << " expected_after=" << sc["estimate_after"] << std::endl;
      }
    }
  }
  TEST_EQUAL(scans > 100, true)
  TEST_EQUAL(mismatches, 0)
}
END_SECTION

START_SECTION(v23_matching_and_vote_reproduce_every_ms3)
{
  const json d = loadJson(kDir + "v23_engine_offset_test_vectors.json");
  ABORT_IF(!d.is_object())
  int scans = 0, votes = 0, mismatches = 0;
  for (const json& v : d["vectors"])
  {
    if (v["level"].get<int>() != 3) continue;   // see the header: the MS2 match is not this rule
    MassOffsetVote vote;
    for (const json& sc : v["scans"])
    {
      ++scans;
      const json& in = sc["inputs"];
      const json& ex = sc["expected"];
      DeconvolvedSpectrum spec(0);
      for (const json& m : in["masses"])
      {
        PeakGroup pg(1, 1, true);
        pg.setMonoisotopicMass(m.get<double>());
        spec.push_back(pg);
      }
      // The ladder as the vector gives it: each sorted mass with the bond it cleaves, read as a prefix ion at that
      // position so the engine's bond rule returns the vector's bond.
      std::vector<MS3FragmentMatcher::TheoreticalMass> theo;
      int max_bond = 1;
      for (size_t i = 0; i < in["ladder_sorted"].size(); ++i)
      {
        MS3FragmentMatcher::TheoreticalMass t;
        t.mass = in["ladder_sorted"][i].get<double>();
        t.position = in["bond_of_ladder_mass"][i].get<int>();
        t.ion_type = "b";
        max_bond = std::max(max_bond, t.position);
        theo.push_back(t);
      }
      const double used = ex["offset_used"].get<double>();
      const MS3FragmentMatcher::OffsetEvidence ev =
        MS3FragmentMatcher::offsetEvidence(spec, theo, used, ex["tolerance_ppm"].get<double>(), max_bond + 1);
      const bool confident = ev.bonds >= kMassOffsetConfidentMs3Bonds;
      bool ok = ev.bonds == ex["n"].get<int>() && confident == ex["confident"].get<bool>();
      ok = ok && sameList(confident ? ev.raw_errors : std::vector<double> {}, ex["errors_fed"]);

      MassOffsetVoteResult vr;
      const bool voted = vote.step(ev.n_masses, confident, ev.wide_errors, used, vr);
      ok = ok && vote.counter() == ex["counter_after"].get<int>();
      ok = ok && voted == !ex["vote"].is_null();
      if (voted && !ex["vote"].is_null())
      {
        ++votes;
        const json& xv = ex["vote"];
        ok = ok && vr.result == engineResult(xv["result"].get<std::string>());
        if (xv.contains("centre"))
          ok = ok && vr.has_centre && near(vr.centre, xv["centre"].get<double>()) && vr.errors == xv["errors"].get<int>()
               && vr.scans == xv["scans"].get<int>() && vr.runner_up == xv["runner_up"].get<int>();
      }
      if (!ok)
      {
        ++mismatches;
        std::cout << "v23 " << v["name"].get<std::string>() << " k=" << sc["k"] << " bonds=" << ev.bonds
                  << " counter=" << vote.counter() << " voted=" << voted << " result=" << vr.result << std::endl;
      }
    }
  }
  TEST_EQUAL(scans > 900, true)
  TEST_EQUAL(votes >= 14, true)   // every vote of the eight MS3 sequences, the decoy false jump included
  TEST_EQUAL(mismatches, 0)
}
END_SECTION

START_SECTION(the_vote_is_armed_only_by_an_authored_ms3_start)
{
  // 20 MS3 with five masses and nothing confident: armed, the counter reaches the limit and a vote is taken; not armed,
  // MassOffsets never steps it.
  MassOffsets armed(0.0, 0.0, true), unarmed(0.0, 0.0, false);
  const std::vector<double> wide = {-185.1, -185.3, -184.9, -185.2, -185.0};
  for (int i = 0; i < 20; ++i)
  {
    armed.onMs3("t", armed.ms3.estimate(), 5, 0, 0, {}, wide);
    unarmed.onMs3("t", unarmed.ms3.estimate(), 5, 0, 0, {}, wide);
  }
  // armed: 100 errors around -185.1 in 20 scans, no runner-up, more than 3 ppm from 0 -- a jump to their median
  TEST_REAL_SIMILAR(armed.ms3.estimate(), -185.1)
  TEST_EQUAL(armed.ms3.count(), 0)   // the restart drops every observation
  TEST_REAL_SIMILAR(unarmed.ms3.estimate(), 0.0)
}
END_SECTION

END_TEST
