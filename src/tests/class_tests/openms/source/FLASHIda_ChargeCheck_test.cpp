// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ChargeCheck.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;

// ADR-0078: the charge check's measurement, on synthetic windows and on replayed decision-survey windows.
//
// The synthetic envelopes are an averagine-like Poisson (lambda = M / 1800) at the averagine spacing; the same envelopes, the same
// deterministic noise and the same co-isolation set were run beforehand through the prototype's Python (flashida-paper
// w08 / w14d), so every bound below is the prototype's own behaviour, not a guess. The replay fixtures are the w14 truth windows
// (FlashIDA/test-data/charge_check, written by flashida-paper w17_export_fixtures.py); their expected counts are spec 10's.

namespace
{
  using CC = ChargeCheck::Status;
  constexpr double ISO = 1.00235;
  constexpr double PROTON = 1.007276;
  constexpr int Z_MIN = 4;   // deconvolution.min_charge / max_charge of the CA and myoglobin methods
  constexpr int Z_MAX = 50;

  using Points = std::vector<std::pair<double, double>>;

  // The averagine-like envelope of mass M at charge z: isotope k at (M + k * ISO + z * PROTON) / z, Poisson(M / 1800) heights,
  // every isotope above 1e-3 of the apex.
  Points envelope(double M, int z, double scale = 1.0)
  {
    const double lam = M / 1800.0;
    const int kmax = static_cast<int>(lam + 6.0 * std::sqrt(lam) + 6.0);
    std::vector<double> p(kmax);
    double top = 0.0;
    for (int k = 0; k < kmax; ++k)
    {
      p[k] = std::exp(-lam + k * std::log(lam) - std::lgamma(k + 1.0));
      top = std::max(top, p[k]);
    }
    Points out;
    for (int k = 0; k < kmax; ++k)
      if (p[k] > 1e-3 * top) out.emplace_back((M + k * ISO + z * PROTON) / z, scale * 1e6 * p[k]);
    return out;
  }

  double apexMz(double M, int z)
  {
    return (M + std::floor(M / 1800.0) * ISO + z * PROTON) / z;
  }

  MSSpectrum spectrumOf(Points pts)
  {
    std::sort(pts.begin(), pts.end());
    MSSpectrum s;
    for (const auto& [mz, it] : pts)
    {
      Peak1D p;
      p.setMZ(mz);
      p.setIntensity(static_cast<float>(it));
      s.push_back(p);
    }
    return s;
  }

  // The window the replay and the prototype use: the commanded m/z +- 0.7 Th.
  ChargeCheck::Measurement check(const Points& pts, double mz, int z_cmd, double r200)
  {
    return ChargeCheck::measure(spectrumOf(pts), mz - 0.7, mz + 0.7, mz, z_cmd, r200, Z_MIN, Z_MAX, ChargeCheck::Settings{});
  }

  bool confident(CC s)
  {
    return s == CC::Confirmed || s == CC::Changed || s == CC::Inconsistent;
  }

  // The prototype's deterministic generator (x <- (1103515245 x + 12345) mod 2^31), so C++ and Python draw the same numbers.
  struct Lcg
  {
    std::uint64_t x;
    explicit Lcg(std::uint64_t seed) : x(seed) {}
    double next()
    {
      x = (1103515245ULL * x + 12345ULL) % 2147483648ULL;
      return static_cast<double>(x) / 2147483648.0;
    }
  };

  struct ReplayWindow
  {
    std::string id;
    double mz = 0, lo = 0, hi = 0, r200 = 0;
    int z_cmd = 0, z_true = 0;
    Points pts;
  };

  // A replay fixture: '#' provenance lines, then "W id mz lo hi z_cmd r200 z_true" rows, each followed by its "P mz intensity" rows.
  std::vector<ReplayWindow> readReplay(const std::string& path)
  {
    std::vector<ReplayWindow> out;
    std::ifstream in(path);
    std::string line;
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::string tag;
      ss >> tag;
      if (tag == "W")
      {
        ReplayWindow w;
        ss >> w.id >> w.mz >> w.lo >> w.hi >> w.z_cmd >> w.r200 >> w.z_true;
        out.push_back(w);
      }
      else if (tag == "P" && !out.empty())
      {
        double mz = 0, it = 0;
        ss >> mz >> it;
        out.back().pts.emplace_back(mz, it);
      }
    }
    return out;
  }

  struct ReplayCounts
  {
    int windows = 0, errors = 0, right = 0, recall = 0, false_alarms = 0, wrong = 0, changed = 0;
    std::vector<double> ms;
  };

  // One pass over a fixture, with the decision the engine makes: a Changed is a change of the commanded charge.
  ReplayCounts replay(const std::vector<ReplayWindow>& ws)
  {
    ReplayCounts c;
    for (const auto& w : ws)
    {
      const MSSpectrum s = spectrumOf(w.pts);
      const auto t0 = std::chrono::steady_clock::now();
      const auto m = ChargeCheck::measure(s, w.lo, w.hi, w.mz, w.z_cmd, w.r200, Z_MIN, Z_MAX, ChargeCheck::Settings{});
      c.ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
      ++c.windows;
      const bool changed = m.status == CC::Changed;
      c.changed += changed;
      if (w.z_true <= 0) continue;   // no reference charge
      const bool err = w.z_true != w.z_cmd;
      c.errors += err;
      c.right += !err;
      c.recall += changed && err && m.z_meas == w.z_true;
      c.false_alarms += changed && !err;
      c.wrong += changed && err && m.z_meas != w.z_true;
    }
    return c;
  }

  double percentile(std::vector<double> v, double q)
  {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t i = std::min(v.size() - 1, static_cast<size_t>(q * (v.size() - 1) + 0.5));
    return v[i];
  }

  const std::string FI_MYO = "../../FlashIDA/test-data/charge_check/myo_ox_windows.tsv";
  const std::string FI_CA = "../../FlashIDA/test-data/charge_check/ca2_windows.tsv";
} // namespace

START_TEST(FLASHIda_ChargeCheck, "$Id$")

/////////////////////////////////////////////////////////////

// The status vocabulary scan_commands.tsv logs (Q7), and the resolvability the gate reads.
// FAILS UNDER: a renamed status, or a resolvability not (1.00235 / z) / (m/z / (R200 sqrt(200 / (m/z)))).
START_SECTION(vocabulary_and_resolvability)
{
  TEST_STRING_EQUAL(ChargeCheck::toString(CC::NotChecked), "not_checked")
  TEST_STRING_EQUAL(ChargeCheck::toString(CC::Confirmed), "confirmed")
  TEST_STRING_EQUAL(ChargeCheck::toString(CC::Changed), "changed")
  TEST_STRING_EQUAL(ChargeCheck::toString(CC::UnverifiedResolution), "unverified_resolution")
  TEST_STRING_EQUAL(ChargeCheck::toString(CC::UnverifiedLowScore), "unverified_low_score")
  TEST_STRING_EQUAL(ChargeCheck::toString(CC::UnverifiedNoSpecies), "unverified_no_species")
  TEST_STRING_EQUAL(ChargeCheck::toString(CC::Inconsistent), "inconsistent")
  TEST_STRING_EQUAL(ChargeCheck::toString(CC::Empty), "empty")
  TEST_STRING_EQUAL(ChargeCheck::toString(CC::Swapped), "swapped")
  TEST_STRING_EQUAL(ChargeCheck::toString(CC::Remeasured), "remeasured")
  TEST_REAL_SIMILAR(ChargeCheck::resolvability(850.0, 20, 120000.0), (1.00235 / 20) / (850.0 / (120000.0 * std::sqrt(200.0 / 850.0))))
  TEST_REAL_SIMILAR(ChargeCheck::resolvability(850.0, 0, 120000.0), 0.0)
  TEST_REAL_SIMILAR(ChargeCheck::resolvability(850.0, 20, 0.0), 0.0)
}
END_SECTION

// A1. Clean single envelopes, 10-30 kDa, z 5-40, at 120k and 240k, commanded at z-2..z+2 (1800 windows): no confident outcome
// (confirmed, changed, inconsistent) ever names a wrong charge; below the resolvability gate the status is unverified_resolution;
// at the true charge the comb reads z wherever the resolvability is >= 3, and never 2z where it is >= 2.5. (High charges of
// light envelopes measure the right z below the 0.30 score, so they are unverified_low_score, not confirmed: spec 4's comb
// shares one lag span across charges, and those envelopes are narrower than it.)
// FAILS UNDER: a comb whose lag span differs per charge (picks 2z / 3z), a half-charge rule picking the lowest charge, a gate
// order that lets a low-resolution window through.
START_SECTION(clean_envelopes_never_confidently_wrong)
{
  int windows = 0, wrong_confident = 0, gate_missed = 0, wrong_at_res3 = 0, doubled = 0;
  for (double r200 : {120000.0, 240000.0})
  {
    for (double M : {10000.0, 15000.0, 20000.0, 25000.0, 30000.0})
    {
      for (int z = 5; z <= 40; ++z)
      {
        const double mz = apexMz(M, z);
        const Points env = envelope(M, z);
        for (int z_cmd = z - 2; z_cmd <= z + 2; ++z_cmd)
        {
          const auto m = check(env, mz, z_cmd, r200);
          ++windows;
          if (confident(m.status) && m.z_meas != z) ++wrong_confident;
          if (ChargeCheck::resolvability(mz, z_cmd, r200) < 2.5 && m.status != CC::UnverifiedResolution) ++gate_missed;
          if (z_cmd != z) continue;
          if (m.resolvability >= 3.0 && m.z_meas != z) ++wrong_at_res3;
          if (m.resolvability >= 2.5 && m.z_meas == 2 * z) ++doubled;
        }
      }
    }
  }
  std::cout << "A1: " << windows << " windows, wrong confident " << wrong_confident << ", gate missed " << gate_missed
            << ", wrong z at res>=3 " << wrong_at_res3 << ", 2z at res>=2.5 " << doubled << std::endl;
  TEST_EQUAL(windows, 1800)
  TEST_EQUAL(wrong_confident, 0)
  TEST_EQUAL(gate_missed, 0)
  TEST_EQUAL(wrong_at_res3, 0)
  TEST_EQUAL(doubled, 0)
}
END_SECTION

// A2. The test proteins' own ranges (spec 6): myoglobin 16951 Da at 120k, z 12-24; CA II 29000 Da at 240k, z 16-30. Wherever the
// commanded charge's resolvability clears 2.5, the true charge is confirmed and a charge one off is changed to the true one.
// FAILS UNDER: anything that loses recall where the spec's data put it: the comb, the refinement, the score.
START_SECTION(myoglobin_and_ca_ranges_confirm_and_correct)
{
  struct Range { double M; double r200; int z_lo; int z_hi; };
  int cases = 0, misses = 0;
  for (const Range& rg : {Range{16951.0, 120000.0, 12, 24}, Range{29000.0, 240000.0, 16, 30}})
  {
    for (int z = rg.z_lo; z <= rg.z_hi; ++z)
    {
      const double mz = apexMz(rg.M, z);
      const Points env = envelope(rg.M, z);
      for (int z_cmd : {z, z - 1, z + 1})
      {
        if (ChargeCheck::resolvability(mz, z_cmd, rg.r200) < 2.5) continue;
        const auto m = check(env, mz, z_cmd, rg.r200);
        ++cases;
        const bool ok = z_cmd == z ? m.status == CC::Confirmed : (m.status == CC::Changed && m.z_meas == z);
        if (!ok)
        {
          ++misses;
          std::cout << "A2 miss: M=" << rg.M << " z=" << z << " z_cmd=" << z_cmd << " -> " << ChargeCheck::toString(m.status)
                    << " z_meas=" << m.z_meas << " score=" << m.score << std::endl;
        }
      }
    }
  }
  TEST_EQUAL(cases > 60, true)
  TEST_EQUAL(misses, 0)
}
END_SECTION

// A3. Co-isolation: a main envelope (12/17/22 kDa, z 12/16/20/24, 120k) with one interferer at 30 % (9/15/25 kDa, z 8/15/22/30)
// whose apex sits -0.4/0/+0.4 Th off (432 windows). The prototype reads a wrong confident charge in 19 of them, all at the
// resolvability edge (22 kDa z 24, res 2.55, read as z 25): bound at 8 %, and most confident outcomes are right.
// FAILS UNDER: an interferer taking over the comb broadly (more than the resolvability-edge cases).
START_SECTION(co_isolated_envelopes_mostly_measure_the_strongest)
{
  int windows = 0, confident_right = 0, confident_wrong = 0;
  for (double M : {12000.0, 17000.0, 22000.0})
  {
    for (int z : {12, 16, 20, 24})
    {
      const double mz = apexMz(M, z);
      for (double M2 : {9000.0, 15000.0, 25000.0})
      {
        for (int z2 : {8, 15, 22, 30})
        {
          for (double off : {-0.4, 0.0, 0.4})
          {
            Points pts = envelope(M, z);
            const double shift = mz + off - apexMz(M2, z2);
            for (const auto& [m2, a2] : envelope(M2, z2, 0.3)) pts.emplace_back(m2 + shift, a2);
            const auto m = check(pts, mz, z, 120000.0);
            ++windows;
            if (confident(m.status)) (m.z_meas == z ? confident_right : confident_wrong)++;
          }
        }
      }
    }
  }
  std::cout << "A3: " << windows << " windows, confident right " << confident_right << ", confident wrong " << confident_wrong
            << std::endl;
  TEST_EQUAL(windows, 432)
  TEST_EQUAL(confident_wrong <= 34, true)   // 8 % of 432; the prototype: 19
  TEST_EQUAL(confident_right >= 240, true)  // the prototype: 269
}
END_SECTION

// A4. Windows with no envelope: nothing is empty; a lone centroid is never confirmed or changed at any commanded charge; random
// noise centroids (40 windows of 20/60/150 centroids, each checked at z_cmd 4..50 step 2: 960 checks) are changed rarely: the
// prototype changes 12, bound at 2 %. (The engine checks only admitted PeakGroups, so a pure-noise window is not its input.)
// FAILS UNDER: a check that acts on a window with no signal, or noise that reads as a comb broadly.
START_SECTION(windows_without_an_envelope)
{
  // nothing in the window
  const auto empty = check({}, 1000.0, 10, 120000.0);
  TEST_EQUAL(empty.status == CC::Empty, true)
  TEST_EQUAL(empty.z_meas, 0)

  // one centroid, every commanded charge
  int lone_acted = 0;
  for (int z_cmd = Z_MIN; z_cmd <= Z_MAX; ++z_cmd)
  {
    const auto m = check({{1000.0, 1e6}}, 1000.0, z_cmd, 120000.0);
    if (m.status == CC::Confirmed || m.status == CC::Changed) ++lone_acted;
  }
  TEST_EQUAL(lone_acted, 0)

  // random noise centroids; the draws in the prototype's order (m/z, then intensity, per centroid)
  int checks = 0, changed = 0;
  for (std::uint64_t seed = 1; seed <= 40; ++seed)
  {
    Lcg g(seed);
    const double mz = 700.0 + 800.0 * g.next();
    const int n = (seed % 3 == 0) ? 20 : (seed % 3 == 1) ? 60 : 150;
    Points pts;
    for (int i = 0; i < n; ++i)
    {
      const double at = mz - 0.75 + 1.5 * g.next();
      const double it = (0.5 + g.next()) * 1e5;
      pts.emplace_back(at, it);
    }
    for (int z_cmd = 4; z_cmd <= 50; z_cmd += 2)
    {
      ++checks;
      changed += check(pts, mz, z_cmd, 120000.0).status == CC::Changed;
    }
  }
  std::cout << "A4: noise checks " << checks << ", changed " << changed << std::endl;
  TEST_EQUAL(checks, 960)
  TEST_EQUAL(changed <= 19, true)   // 2 % of 960; the prototype: 12
}
END_SECTION

// A5. Intensity noise: every A1 envelope at resolvability >= 3, each centroid scaled by 1 + 0.1 u, u uniform in [-1, 1]
// (deterministic draws), still reads its own charge (the prototype: 234 of 234).
// FAILS UNDER: a comb that needs exact Poisson heights.
START_SECTION(intensity_noise_keeps_the_charge)
{
  Lcg g(99);
  int windows = 0, right = 0;
  for (double r200 : {120000.0, 240000.0})
  {
    for (double M : {10000.0, 15000.0, 20000.0, 25000.0, 30000.0})
    {
      for (int z = 5; z <= 40; ++z)
      {
        const double mz = apexMz(M, z);
        if (ChargeCheck::resolvability(mz, z, r200) < 3.0) continue;
        Points pts;
        for (const auto& [m, a] : envelope(M, z))
        {
          const double u = g.next();
          pts.emplace_back(m, a * (1.0 + 0.1 * (2.0 * u - 1.0)));
        }
        ++windows;
        right += check(pts, mz, z, r200).z_meas == z;
      }
    }
  }
  TEST_EQUAL(windows, 234)
  TEST_EQUAL(right, windows)
}
END_SECTION

// R1, R2. The replay (spec 10): the decision-survey centroids of the oxidized-myoglobin (156 windows, 81 charge errors, 120k) and
// CA II (56 windows, 27 errors, 240k) commanded precursors with a resolved truth. The prototype: myoglobin 75 / 81 corrected,
// CA 9 / 27, no false alarm, no correction to a wrong charge. Timing over both (R3): the spec's budget is 2 ms median and 10 ms
// at the 99th percentile on the acquisition PC; CI shares this host with two other runners and ctest -j 6, so twice that.
// FAILS UNDER: any drift from the prototype's measurement on real windows.
START_SECTION(replay_myoglobin_and_ca)
{
  const auto myo_w = readReplay(FI_MYO);
  const auto ca_w = readReplay(FI_CA);
  TEST_EQUAL(myo_w.size(), 156)
  TEST_EQUAL(ca_w.size(), 56)

  const ReplayCounts myo = replay(myo_w);
  const ReplayCounts ca = replay(ca_w);
  std::cout << "R1 myoglobin: recall " << myo.recall << "/" << myo.errors << ", false alarms " << myo.false_alarms << "/" << myo.right
            << ", wrong " << myo.wrong << std::endl;
  std::cout << "R2 CA II: recall " << ca.recall << "/" << ca.errors << ", false alarms " << ca.false_alarms << "/" << ca.right
            << ", wrong " << ca.wrong << std::endl;
  TEST_EQUAL(myo.errors, 81)
  TEST_EQUAL(myo.recall >= 73, true)   // >= 90 % of 81; the prototype: 75
  TEST_EQUAL(myo.false_alarms, 0)
  TEST_EQUAL(myo.wrong, 0)
  TEST_EQUAL(ca.errors, 27)
  TEST_EQUAL(ca.recall >= 9, true)     // >= 30 % of 27; the prototype: 9
  TEST_EQUAL(ca.false_alarms, 0)
  TEST_EQUAL(ca.wrong, 0)

  std::vector<double> ms = myo.ms;
  ms.insert(ms.end(), ca.ms.begin(), ca.ms.end());
  const double median = percentile(ms, 0.5);
  const double p99 = percentile(ms, 0.99);
  std::cout << "R3 timing over " << ms.size() << " windows: median " << median << " ms, p99 " << p99 << " ms" << std::endl;
  TEST_EQUAL(median <= 4.0, true)
  TEST_EQUAL(p99 <= 20.0, true)
}
END_SECTION

// The E. coli LC replay (spec 10: at most 1.0 % of 2283 commanded MS2 changed, at most 2 of those with a right reference charge),
// run by hand against a built test exe: CHARGE_CHECK_EXTRA_FIXTURE=<flashida-paper>/charge_check_ecoli.tsv. Too large for CI; the
// section passes without the variable.
START_SECTION(replay_ecoli_when_supplied)
{
  const char* path = std::getenv("CHARGE_CHECK_EXTRA_FIXTURE");
  if (path == nullptr || *path == '\0')
  {
    std::cout << "E. coli replay: CHARGE_CHECK_EXTRA_FIXTURE not set, skipped" << std::endl;
    TEST_EQUAL(true, true)
  }
  else
  {
    const auto ws = readReplay(path);
    const ReplayCounts c = replay(ws);
    std::cout << "E. coli replay: " << c.windows << " windows, changed " << c.changed << " ("
              << (c.windows > 0 ? 100.0 * c.changed / c.windows : 0.0) << " %); with a reference " << (c.right + c.errors)
              << ", false alarms " << c.false_alarms << ", to the reference " << c.recall << "; median "
              << percentile(c.ms, 0.5) << " ms, p99 " << percentile(c.ms, 0.99) << " ms" << std::endl;
    TEST_EQUAL(c.windows > 0, true)
    TEST_EQUAL(c.changed <= c.windows / 100, true)
    TEST_EQUAL(c.false_alarms <= 2, true)
  }
}
END_SECTION


/////////////////////////////////////////////////////////////
END_TEST
