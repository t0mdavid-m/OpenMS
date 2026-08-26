// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Shared helpers for the FLASHIda logging tests:
//   * TSV / IDA-log loaders and parsers
//   * full MS1->MS2->MS3 cycle driver and single-exploration-group driver
//   * JSON config builders (with runtime log paths injected)
//   * plausibility predicates used by FLASHIda_LoggingFields_test
//
// These live in an anonymous namespace so each ClassTest translation unit (every
// *_test.cpp is its own executable) gets its own internal-linkage copy -- no ODR
// concerns across the separate test binaries.

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>  // encode/decode for distinct MS1 ids
#include <OpenMS/SYSTEM/File.h>  // freshLogDir: makeDir / removeDirRecursively
#include "FLASHIda_TestAccess.h"  // FLASHIdaTestAccess::explorationActive (private-state access)

#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Test-only header: pull OpenMS into scope so includers need not order the include
// relative to their own `using namespace OpenMS;` (and so the unqualified FLASHIda /
// ScanCommand references below resolve regardless of include position).
using namespace OpenMS;

namespace
{
  // --- Config bounds for the cytC fixture (must match the JSON configs below) ---
  const double FI_MIN_MASS = 500.0;
  const double FI_MAX_MASS = 50000.0;
  const int    FI_MIN_CHARGE = 4;
  const int    FI_MAX_CHARGE = 50;
  const double FI_RT_WINDOW = 180.0;

  // ----------------------------------------------------------------------------
  // Plausibility predicates (drift-stable: ranges / sets / signs, not exact scalars)
  // ----------------------------------------------------------------------------

  inline bool finiteVal(double x) { return std::isfinite(x); }
  inline bool posFinite(double x) { return std::isfinite(x) && x > 0.0; }
  inline bool nonNegFinite(double x) { return std::isfinite(x) && x >= -1e-9; }
  inline bool inUnit(double x) { return std::isfinite(x) && x >= -1e-6 && x <= 1.0 + 1e-6; }
  inline bool inMass(double x) { return std::isfinite(x) && x >= FI_MIN_MASS - 1.0 && x <= FI_MAX_MASS + 1.0; }
  inline bool inChargeD(double z)
  {
    return std::isfinite(z) && z >= FI_MIN_CHARGE && z <= FI_MAX_CHARGE && std::floor(z) == z;
  }
  // Fragment / MSn-precursor charge: integer in [1, parent charge] (and <= max_charge). Unlike a
  // precursor charge it is NOT bounded below by min_charge — a fragment can be charge 1..3.
  inline bool inFragCharge(double zf, double zp)
  {
    return std::isfinite(zf) && std::floor(zf) == zf && zf >= 1.0 && zf <= zp && zf <= FI_MAX_CHARGE;
  }
  inline bool inDuration(double t) { return std::isfinite(t) && t >= 0.0 && t < 3600000.0; }
  inline bool inActivationSet(const std::string& s)
  {
    return s == "HCD" || s == "ETD" || s == "EThcD" || s == "CID";
  }
  // tracking_id: exactly 3 chars, each printable ASCII 0x21-0x7E (the base-94 alphabet).
  inline bool isTrackingId(const std::string& s)
  {
    if (s.size() != 3) return false;
    for (char c : s)
    {
      unsigned char u = static_cast<unsigned char>(c);
      if (u < 0x21 || u > 0x7E) return false;
    }
    return true;
  }

  // Parse a TSV cell to double; NaN on malformed input (NaN fails every finite predicate).
  inline double toD(const std::string& s)
  {
    try { return std::stod(s); }
    catch (...) { return std::numeric_limits<double>::quiet_NaN(); }
  }
  // Split a cell into ';'-delimited (or arbitrary-delimited) tokens.
  inline std::vector<std::string> splitTokens(const std::string& s, char delim)
  {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string t;
    while (std::getline(ss, t, delim)) out.push_back(t);
    return out;
  }

  // ----------------------------------------------------------------------------
  // Spectrum / TSV loaders (verbatim from FLASHIda_Logging_test / FLASHIda_exploration_test)
  // ----------------------------------------------------------------------------

  struct ScanData
  {
    std::vector<double> mzs;
    std::vector<double> ints;
    double rt = 0.0;
    std::string scan_id;
  };

  // The INSTRUMENT scan number a fed spectrum claims, for processScan's trailing argument.
  //
  // Fixtures carry it in their "Spec scan=N" header, which loadTsvScans already parses into
  // ScanData::scan_id -- a field that was populated and then never read by anything. It is the same
  // quantity the C# twin sources from Header["Scan"] (MockMsScan.cs:267,340), so both mirrors of the
  // one drive contract feed the engine the same number from the same place.
  //
  // 0 for a fixture that names none, which the engine reads as "not supplied" (ADR-0035).
  inline int instrumentScanNumberOf(const ScanData& s)
  {
    if (s.scan_id.empty()) return 0;
    try { return std::stoi(s.scan_id); } catch (...) { return 0; }
  }

  // Parse multi-scan TSV: "Spec scan=N\tRT" headers, "mz\tintensity" data lines
  inline std::vector<ScanData> loadTsvScans(const std::string& path)
  {
    std::vector<ScanData> scans;
    std::ifstream f(path);
    if (!f.good()) return scans;
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

  // Push all scans through processScan (ms_level=1), return total command count.
  // Each MS1 gets a DISTINCT, valid 3-char tracking-id prefix (engine decodes desc[0:3] as the MS1 id),
  // drawn from a high base (800000+) so it never collides with the engine's low command-id counter — so
  // every MS2's parent_tracking_id resolves to a distinct fed MS1 id (was "scan_"+id → all decode to "sca").
  // When fed_ms1_ids != nullptr, records each fed id (without the 'S' marker) for parent-resolution checks.
  inline int pushAllScans(FLASHIda* ida, const std::vector<ScanData>& scans,
                          std::vector<std::string>* fed_ms1_ids = nullptr)
  {
    int total = 0;
    int i = 0;
    for (const auto& scan : scans)
    {
      std::string id = ScanCommandQueue::encode(800000 + (i++));
      total += ida->processScan(scan.mzs.data(), scan.ints.data(),
                                (int)scan.mzs.size(), scan.rt, 1, (id + "S").c_str(),
                                0.0, instrumentScanNumberOf(scan));
      if (fed_ms1_ids) fed_ms1_ids->push_back(id);
    }
    return total;
  }

  // Parse a TSV file into rows of string vectors (first row = header)
  struct TSVFile
  {
    std::vector<std::string> headers;
    std::vector<std::vector<std::string>> rows;

    static TSVFile parse(const std::string& path)
    {
      TSVFile result;
      std::ifstream f(path);
      std::string line;
      bool first = true;
      while (std::getline(f, line))
      {
        // Manual split that emits N+1 fields for N tabs, preserving a trailing empty field — matching
        // C# String.Split('\t'). The prior std::getline(iss, col, '\t') loop DROPPED a trailing empty
        // field, so scan_results rows whose new last column (winner_tracking_id) is empty tokenized to
        // one fewer column than the header and failed row.size()==headers.size().
        std::vector<std::string> cols;
        size_t start = 0;
        while (true)
        {
          size_t t = line.find('\t', start);
          if (t == std::string::npos) { cols.push_back(line.substr(start)); break; }
          cols.push_back(line.substr(start, t - start));
          start = t + 1;
        }

        if (first) { result.headers = cols; first = false; }
        else        { result.rows.push_back(cols); }
      }
      return result;
    }

    int colIndex(const std::string& name) const
    {
      for (size_t i = 0; i < headers.size(); i++)
        if (headers[i] == name) return static_cast<int>(i);
      return -1;
    }
  };

  // Safe cell access: "" when the column is missing or the row is short.
  inline std::string cell(const TSVFile& tsv, const std::vector<std::string>& row, const std::string& name)
  {
    int c = tsv.colIndex(name);
    if (c < 0 || c >= (int)row.size()) return std::string();
    return row[c];
  }

  // ----------------------------------------------------------------------------
  // Full MS1->MS2->MS3 cycle driver (verbatim from FLASHIda_Logging_test)
  // ----------------------------------------------------------------------------

  struct CycleResult
  {
    std::vector<ScanCommand> ms2_cmds;
    std::vector<ScanCommand> ms3_cmds;
    std::vector<ScanCommand> ms1_cmds;
    int total_dequeued = 0;
  };

  // Decode the trailing precursor ion of an MS3 scan_description per the descriptor contract:
  //   {id}R{mass}k@{charge}{ion_type}{ion_index}  (the 'k' is part of the mass token, before '@').
  // Find the LAST '@', skip the run of fragment-charge digits after it, require an ion_type char in
  // {a,b,c,x,y,z}, then an all-digit ion_index >= 1. Returns false on the no-ion form {id}R{mass}k@{charge}
  // (nothing after the charge digits) -- tolerated by callers. On success, ion_key = ion_type + ion_index
  // (e.g. "b44"). Kept identical, byte-for-byte, to decodeTrailingIon in FLASHIda_LoggingFields_test.cpp.
  // [ION-DECODE C#<->C++ — see docs/kb/test-harness] byte-for-byte twin of C# DecodeIonFromScanDescription
  inline bool decodeTrailingIonKey(const std::string& d, std::string& ion_key)
  {
    auto at = d.rfind('@');
    if (at == std::string::npos) return false;
    size_t i = at + 1;
    while (i < d.size() && std::isdigit(static_cast<unsigned char>(d[i]))) ++i;  // skip fragment charge digits
    if (i >= d.size()) return false;                                             // no-ion form -> tolerated
    char t = d[i];
    if (!(t == 'a' || t == 'b' || t == 'c' || t == 'x' || t == 'y' || t == 'z')) return false;
    std::string idx = d.substr(i + 1);
    if (idx.empty() || idx.find_first_not_of("0123456789") != std::string::npos) return false;
    if (std::atoi(idx.c_str()) < 1) return false;
    ion_key = std::string(1, t) + idx;
    return true;
  }

  // [DRAIN-CONTRACT C#<->C++ — see docs/kb/test-harness] AcqResult + the one canonical interleaved driver.
  struct AcqResult
  {
    std::vector<ScanCommand> ms1_cmds, ms2_cmds, ms3_cmds;
    // Cross-level capability (C++-only, no C# twin): EVERY dequeued command — workload AND idle AGC/MS1
    // ticks — pushed in raw dequeue order, so a caller can assert cross-level interleave (e.g. a prio-0
    // CV-transition MS1 drained BEFORE the prio-2 MS2s) that the per-level buckets above lose. Additive;
    // the drain contract is unchanged.
    std::vector<ScanCommand> all_cmds;
    // Parallel to all_cmds (|all_active| == |all_cmds|): the engine's exploration_active_ flag captured AT
    // each command's dequeue (FLASHIda::explorationActive()). Lets a caller scope an assertion to exactly
    // the exploration-ACTIVE window (e.g. cycle-time MS1 suppression holds only while a group is active).
    std::vector<char> all_active;
    int total_dequeued = 0;
    // single_group_only bookkeeping (see runInterleaved): the processScan return of the first MS1 that
    // forms a group (== # commands it pushed), and the running sum of MS2-feed processScan returns. Both
    // stay 0 for the default full-drive callers; driveOneExplorationGroup reads them to build ExplResult.
    int first_group_commands = 0;
    int ms2_feed_returns = 0;
  };

  // runInterleaved: the C++ implementation of the ground-truth interleaved engine-id-echo drive contract
  // (docs/kb/test-harness/README.md); twin of C# ContinuityTestHarness.PushScanAndDrainFull. Pull one command
  // at a time and feed exactly one response per requested command, stamped with the engine's own
  // scan_description; the engine paces the surveys. Terminate on idle>=3 (queue drained) or max_iters.
  //   ms1_scans   : fed one per survey command, in order (nMs1 = size); further MS1 surveys are idle ticks.
  //   ms2_scans   : ms2_scans[0] fed for every MS2 command.
  //   ms3_ion_map : per-ion MS3 fixtures (decode trailing ion -> map; absent/empty => SKIP, never fabricate).
  //                 When null, the legacy MS2-as-MS3 shortcut is used (C++ plausibility only — not the contract).
  //   single_group_only : C++-side convenience (no C# twin) — once the FIRST fed MS1 forms a group (processScan
  //                 returns >0), stop feeding further MS1 surveys (they become idle ticks) but keep draining that
  //                 group's MS2 variants + MS3. Records r.first_group_commands; used by driveOneExplorationGroup.
  //                 Default false leaves the core pull->classify->feed->idle>=3 contract byte-identical.
  //   rt_offset : added to every fed MS1 survey rt (default 0.0). Lets a caller drive the SAME persistent
  //               engine a second time with the surveys shifted in retention time (e.g. CBE two-pass
  //               exclusion: 2nd pass inside vs outside rt_window). MS2/MS3 feeds reuse the source rt as-is.
  //   ms2_ce_map : optional CE-keyed MS2-spectrum source (Task E.2-4) — twin of the C# PushScanAndDrainFull
  //               ms2CeMap. When non-null, an MS2 command's fed spectrum is selected by its stage-0 collision
  //               energy (round(cmd.stages[0].collision_energy)) so each CE variant of the exploration sweep
  //               gets its own energy-resolved fixture. NO FALLBACK: an unmapped CE hard-fails (throws) — a
  //               silent fall-through to ms2_scans[0] would collapse the sweep. When null, ms2_scans[0] is fed
  //               for every MS2 command (existing behaviour). The MS2-bucket recording is unchanged.
  inline AcqResult runInterleaved(FLASHIda* ida,
                                  const std::vector<ScanData>& ms1_scans,
                                  const std::vector<ScanData>& ms2_scans,
                                  const std::map<std::string, std::vector<ScanData>>* ms3_ion_map = nullptr,
                                  int max_iters = 600,
                                  bool single_group_only = false,
                                  double rt_offset = 0.0,
                                  const std::map<int, ScanData>* ms2_ce_map = nullptr)
  {
    AcqResult r;
    const int n_ms1 = static_cast<int>(ms1_scans.size());
    int idle = 0, ms1_fed = 0;
    bool group_formed = false;
    ScanCommand cmd{};
    for (int it = 0; it < max_iters && idle < 3; ++it)
    {
      if (ida->getNextScanCommand(cmd) != 1) break;
      ++r.total_dequeued;
      // Cross-level capability (C++-only, no C# twin): record EVERY dequeued command — workload AND idle
      // AGC/MS1 ticks — in raw dequeue order, BEFORE the idle-tick continue and the per-level bucketing,
      // so callers can assert cross-level interleave the per-level buckets lose.
      r.all_cmds.push_back(cmd);
      r.all_active.push_back(FLASHIdaTestAccess::explorationActive(*ida) ? 1 : 0);   // engine's exploration-active flag AT dequeue
      // Idle tick: AGC, empty descriptor, an MS1 re-survey after all ms1_scans have been fed, or (in
      // single_group_only mode) any MS1 survey once the first group has already formed.
      if (cmd.is_agc || cmd.scan_description[0] == '\0' || (cmd.msn_level <= 1 && ms1_fed >= n_ms1)
          || (single_group_only && group_formed && cmd.msn_level <= 1))
      {
        ++idle;
        cmd = ScanCommand{};
        continue;
      }
      idle = 0;
      if (cmd.msn_level <= 1)
      {
        r.ms1_cmds.push_back(cmd);
        const ScanData& s = ms1_scans[ms1_fed++];
        // F7: echo the command's FAIMS CV back to the engine on the re-fed MS1 (the C++ channel for the
        // CV is the processScan faims_cv argument — the C# twin reads it from the IMsScan "FAIMS CV"
        // trailer in PushScanAndDrainFull). Without it the re-fed MS1 carried CV 0.0 and FAIMS cycling
        // never observed the commanded CV, so downstream CV binding was lost.
        int ret = ida->processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt + rt_offset, 1,
                                   cmd.scan_description, cmd.faims_cv, instrumentScanNumberOf(s));
        if (single_group_only && ret > 0 && !group_formed) { r.first_group_commands = ret; group_formed = true; }
      }
      else if (cmd.msn_level == 2)
      {
        r.ms2_cmds.push_back(cmd);
        // CE-keyed MS2 spectrum (Task E.2-4): select the energy-resolved fixture for THIS variant's stage-0
        // collision energy. NO FALLBACK — an unmapped CE hard-fails (mirrors C# PushScanAndDrainFull).
        const ScanData* ms2_src = ms2_scans.empty() ? nullptr : &ms2_scans[0];
        if (ms2_ce_map != nullptr)
        {
          int ce = (int)std::round(cmd.num_stages > 0 ? cmd.stages[0].collision_energy : 0.0);
          auto cit = ms2_ce_map->find(ce);
          if (cit == ms2_ce_map->end())
            throw std::runtime_error("runInterleaved: MS2 command collision energy " + std::to_string(ce) +
                                     " has no CE-map fixture (no fallback).");
          ms2_src = &cit->second;
        }
        if (ms2_src != nullptr)
          r.ms2_feed_returns += ida->processScan(ms2_src->mzs.data(), ms2_src->ints.data(),
                           (int)ms2_src->mzs.size(), ms2_src->rt, 2, cmd.scan_description,
                           0.0, instrumentScanNumberOf(*ms2_src));
      }
      else  // msn_level >= 3
      {
        r.ms3_cmds.push_back(cmd);
        const std::vector<ScanData>* spectra = nullptr;
        if (ms3_ion_map != nullptr)
        {
          std::string ion_key;
          if (decodeTrailingIonKey(cmd.scan_description, ion_key))
          {
            auto mit = ms3_ion_map->find(ion_key);
            if (mit != ms3_ion_map->end() && !mit->second.empty()) spectra = &mit->second;  // real per-ion MS3
            // ion absent / no-ion descriptor -> skip (never fabricate)
          }
        }
        else
        {
          spectra = ms2_scans.empty() ? nullptr : &ms2_scans;  // legacy MS2-as-MS3 shortcut (no manifest)
        }
        if (spectra != nullptr && !spectra->empty())
          ida->processScan((*spectra)[0].mzs.data(), (*spectra)[0].ints.data(),
                           (int)(*spectra)[0].mzs.size(), (*spectra)[0].rt, 3, cmd.scan_description,
                           0.0, instrumentScanNumberOf((*spectra)[0]));
      }
      cmd = ScanCommand{};
    }
    return r;
  }

  // runFullCycle: legacy 2-pass cycle, now a THIN WRAPPER over runInterleaved (same CycleResult fields).
  // When ms3_ion_map is non-null, each level-3 command's trailing ion is decoded and fed from the manifest
  // (absent/no-ion -> skipped); when null, the MS2 spectrum is fed back as the MS3 scan (shortcut).
  inline CycleResult runFullCycle(FLASHIda* ida,
                                  const std::vector<ScanData>& ms1_scans,
                                  const std::vector<ScanData>& ms2_scans,
                                  const std::map<std::string, std::vector<ScanData>>* ms3_ion_map = nullptr)
  {
    // Budget the drive so EVERY input MS1 is fed (idle>=3 ends it, not the iteration cap). Sized generously
    // above the worst-case per-cycle command count (survey + MS2s + MS3s + AGCs) summed over all scans.
    const int budget = 256 + 64 * static_cast<int>(ms1_scans.size() + ms2_scans.size());
    AcqResult a = runInterleaved(ida, ms1_scans, ms2_scans, ms3_ion_map, budget);
    CycleResult result;
    result.ms1_cmds = a.ms1_cmds;
    result.ms2_cmds = a.ms2_cmds;
    result.ms3_cmds = a.ms3_cmds;
    result.total_dequeued = a.total_dequeued;
    return result;
  }

  // Full-acquisition id-chaining driver — now a THIN WRAPPER over runInterleaved (the canonical contract).
  // Feeds the SINGLE ms1 spectrum n_ms1 times with an RT step (so each survey re-selects its precursor and
  // forms a distinct log group); MS3 uses the MS2-as-MS3 shortcut (no manifest).
  inline AcqResult runFullAcquisition(FLASHIda* ida, const ScanData& ms1, const ScanData& ms2,
                                      int max_iters = 300, int n_ms1 = 1, double ms1_rt_step = 1000.0)
  {
    std::vector<ScanData> ms1_scans;
    ms1_scans.reserve(n_ms1);
    for (int i = 0; i < n_ms1; ++i)
    {
      ScanData s = ms1;
      s.rt = ms1.rt + (double)i * ms1_rt_step;
      // A real instrument stamps a FRESH, increasing scan number on every scan it acquires, so
      // replaying ONE fixture n times has to step the id alongside the RT. Left unstepped, all n
      // surveys claim the same instrument scan and their ida.log entries collapse onto a single map
      // key -- which is what ida_log_multi_scan_distinct_keys exists to catch.
      if (!ms1.scan_id.empty())
        s.scan_id = std::to_string(instrumentScanNumberOf(ms1) + i);
      ms1_scans.push_back(s);
    }
    return runInterleaved(ida, ms1_scans, std::vector<ScanData>{ms2}, nullptr, max_iters);
  }

  // Hard parent-resolution check: every non-empty parent_tracking_id in the commands TSV must resolve to
  // a known id = {fed MS1 input ids} U {command tracking_ids}. Returns false + sets err on the first miss.
  inline bool validateParentTrackingIds(const TSVFile& cmds, const std::vector<std::string>& fed_ms1_ids,
                                        std::string& err)
  {
    std::set<std::string> known(fed_ms1_ids.begin(), fed_ms1_ids.end());
    int id_col = cmds.colIndex("tracking_id");
    int par_col = cmds.colIndex("parent_tracking_id");
    if (id_col < 0 || par_col < 0) { err = "missing tracking_id/parent_tracking_id column"; return false; }
    for (const auto& row : cmds.rows)
      if (id_col < (int)row.size() && !row[id_col].empty()) known.insert(row[id_col]);
    for (const auto& row : cmds.rows)
    {
      if (par_col >= (int)row.size()) continue;
      const std::string& parent = row[par_col];
      if (parent.empty()) continue;
      if (known.count(parent) == 0) { err = "unresolved parent_tracking_id: '" + parent + "'"; return false; }
    }
    return true;
  }

  // ----------------------------------------------------------------------------
  // Give a test section its own empty log folder.
  //
  // Isolation used to come from a unique FILENAME per stream in the one shared ctest CWD.
  // The engine now writes fixed basenames under runtime.log_dir, so isolation moves to the
  // FOLDER and every engine instance needs its own -- including the several sections that
  // construct two or three engines in a row.
  //
  // IdaLogger deliberately does NOT create directories (the host owns the filesystem), and a
  // missing folder leaves every stream silently closed. So this must be called before building
  // the config, and it is the ONLY place a test should name a log directory.
  // ----------------------------------------------------------------------------
  inline std::string freshLogDir(const std::string& tag)
  {
    const std::string dir = "testlogs/" + tag;
    File::removeDirRecursively(dir);  // streams open in append mode; start from nothing
    File::makeDir(dir);               // recursive, idempotent, no-throw
    return dir;
  }

  // ----------------------------------------------------------------------------
  // JSON config builder writing all five log streams into one folder.  enable_ms3 emits the
  // proven inclusion-pinned cytC MS3 recipe (target_mode=1 + inclusion_cytc.txt + the
  // M-starting proteoform + ms3 stage); callers feed ms1_cytc + ms2_cytc_fresh_scan57.
  //
  // RENAMED from buildJsonWithRuntime, and the rename is load-bearing rather than cosmetic.
  // Under the new parameter list a stale 3-or-more-argument call would still COMPILE -- the
  // second argument is now `bool enable_ms3`, and a `const char*` converts to `true` -- so it
  // would silently force MS3 on and blank the activation string. Renaming makes every stale
  // call a named compile error instead.
  // ----------------------------------------------------------------------------
  inline std::string buildJsonWithLogDir(const std::string& log_dir,
                                         bool enable_ms3 = false,
                                         const std::string& ms2_activation = "HCD",
                                         int ms3_max_targets = 0)
  {
    std::string targeting_val     = enable_ms3 ? "inclusion" : "none";
    std::string inclusion_list_val = enable_ms3 ? "../../FlashIDA/test-data/configs/inclusion_cytc.txt" : "";
    std::string ms3_settings = enable_ms3
      ? R"(,
        "ms3": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 35, "resolution": 120000 })"
      : "";

    // characterization.mode is now the whole MS3 gate -- it replaces what used to be two separate
    // selection strings plus the presence of an ms3 scan block.
    //
    // BUDGET, and why the number is what it is. `ms3_max_targets` used to be appended to
    // selection_strategy.ms3.max_targets, which the engine NEVER READ (the budget is taken from
    // level 2). So §T9's "exhaustive MS3" request has always been silently ignored and every caller
    // has actually run the level-2 default. For hand-written bridge JSON like this that default was
    // the C++ literal 10, NOT the C# 3 -- so 10 is the behaviour-preserving value and is now stated
    // explicitly, because characterization.max_targets defaults to 3.
    //
    // Honouring ms3_max_targets for real would CHANGE §T9's behaviour, so it is deliberately not
    // done here: this rewrite preserves behaviour, and turning that request on is a separate,
    // golden-moving decision.
    const int effective_budget = 10;
    std::string ms3_max_targets_json = ", \"max_targets\": " + std::to_string(effective_budget);
    (void)ms3_max_targets;  // inert since before this rewrite; see above

    std::string ms3_block = enable_ms3
      ? R"("characterization": { "mode": "ambiguity", "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE")"
        + ms3_max_targets_json + R"( },)"
      : R"("characterization": { "mode": "off" },)";
    // reaction_time is an ETD-family parameter; only emit it for ETD/EThcD so the HCD
    // default stays byte-identical to the original FLASHIda_Logging_test config.
    std::string ms2_rt = (ms2_activation == "ETD" || ms2_activation == "EThcD")
                         ? ", \"reaction_time\": 10.0" : "";

    std::ostringstream oss;
    oss << R"({
      "deconvolution": {
        "score_threshold": 0.0, "tqscore_threshold": 0.9,
        "min_charge": 4, "max_charge": 50,
        "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
      },
      "precursor_selection": {
        "rt_window": 180, "targeting": ")" << targeting_val << R"(",
        "consider_all_charges": false,
        "strict_inclusion": false, "tie_threshold": 0.1,
        "rank_by": "qscore", "max_precursors": 5
      },
      "flashtnt": { "min_length": 3, "max_length": 8 },
      "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
      "faims": { "cv_values": [], "max_cv_skip": 0 },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
        "ms2": { "analyzer": "Orbitrap", "activation": ")" << ms2_activation << R"(", "collision_energy": 29)" << ms2_rt << R"(, "resolution": 120000 })" << ms3_settings << R"(
      },
      "scheduling": {
        "cycle_time": { "enabled": false, "value_ms": 60000 },
        "scan_timeout": { "enabled": false, "value_ms": 30000 },
        "agc_interval_seconds": 999999
      },
      )" << ms3_block << R"(
      "files": { "target_logs": [], "fasta": "", "inclusion_list": ")" << inclusion_list_val << R"(", "ptm_list": "" },
      "runtime": { "log_dir": ")" << log_dir << R"(" }
    })";
    return oss.str();
  }

  // Insert a "runtime" block naming the log folder into an arbitrary config JSON that lacks
  // one, immediately before its closing brace.  Used to drive the exploration configs (which
  // have no runtime section) and capture their log files.
  //
  // Keeps its name: every existing call passes 3-4 arguments and this takes 2, so a stale call
  // is an arity error rather than a silent conversion.
  inline std::string injectRuntime(std::string cfg, const std::string& log_dir)
  {
    auto last = cfg.rfind('}');
    if (last == std::string::npos) return cfg;
    cfg.insert(last, ", \"runtime\": {\"log_dir\": \"" + log_dir + "\"}");
    return cfg;
  }

  // ----------------------------------------------------------------------------
  // MS2 exploration config (mass_count, CE 20-40 step 5) -- cytC, no runtime block.
  // Mirrors FLASHIda_exploration_test::exploration_config.
  // ----------------------------------------------------------------------------
  inline std::string explorationConfig()
  {
    return R"({
      "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10] },
      "precursor_selection": { "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 3,
        "exploration": { "metric": "mass_count", "ce_min": 20.0, "ce_max": 40.0, "ce_step": 5.0 } },
      "flashtnt": { "min_length": 3, "max_length": 8 },
      "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
      "faims": { "cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
        "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        "ms3": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 35, "resolution": 120000 }
      },
      "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 }, "agc_interval_seconds": 999999 },
      "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
      "characterization": { "mode": "off", "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE", "max_targets": 3 },
      "conditional_ms2": false
    })";
  }

  // MS2 exploration (mass_count) + MS3 exploration (fragment_count, CID 15-35) — cytC, no runtime block.
  // SHARED ground truth for FLASHIda_exploration_test::ms3_exploration_config and the FLASHIda_LoggingFields
  // MS3-exploration golden-column section. inclusionPinCytc() turns it into the matching cytC recipe.
  inline std::string ms3ExplorationConfig()
  {
    return R"({
      "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10] },
      "precursor_selection": { "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 3,
        "exploration": { "metric": "mass_count", "ce_min": 20.0, "ce_max": 40.0, "ce_step": 5.0 } },
      "flashtnt": { "min_length": 3, "max_length": 8 },
      "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
      "faims": { "cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
        "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        "ms3": { "analyzer": "Orbitrap", "activation": "CID", "collision_energy": 25, "resolution": 120000 }
      },
      "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 } },
      "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
      "characterization": { "mode": "ambiguity", "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE", "max_targets": 3,
        "exploration": { "metric": "fragment_count", "ce_min": 15.0, "ce_max": 35.0, "ce_step": 5.0 } },
      "conditional_ms2": false
    })";
  }

  // MS2 exploration (mass_count) + MS3 SELECTION only (intensity, no MS3 exploration) — cytC, no runtime block.
  // SHARED ground truth for FLASHIda_exploration_test::ms3_selection_only_config and the FLASHIda_LoggingFields
  // exploration-follow-up golden-column section (add a non-tolerance override to make the winner re-acquire).
  inline std::string ms3SelectionOnlyConfig()
  {
    return R"({
      "deconvolution": { "score_threshold": 0.0, "tqscore_threshold": 0.9, "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10] },
      "precursor_selection": { "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false, "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 3,
        "exploration": { "metric": "mass_count", "ce_min": 20.0, "ce_max": 40.0, "ce_step": 5.0 } },
      "flashtnt": { "min_length": 3, "max_length": 8 },
      "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
      "faims": { "cv_values": [], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
        "ms2": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 },
        "ms3": { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 35, "resolution": 120000 }
      },
      "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 } },
      "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
      "characterization": { "mode": "ambiguity", "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE", "max_targets": 3 },
      "conditional_ms2": false
    })";
  }

  // Derive an inclusion-pinned, MS3-capable variant of an exploration config (verbatim
  // from FLASHIda_exploration_test): pin the cytC precursor + the validatable M-starting
  // proteoform so real ms2_cytc_fresh_scan57 b/y fragments match -> MS3 fires.
  inline std::string inclusionPinCytc(std::string cfg)
  {
    auto rep = [&cfg](const std::string& from, const std::string& to) {
      auto p = cfg.find(from);
      if (p != std::string::npos) cfg.replace(p, from.size(), to);
    };
    rep("\"targeting\": \"none\"", "\"targeting\": \"inclusion\"");
    rep("\"inclusion_list\": \"\"",
        "\"inclusion_list\": \"../../FlashIDA/test-data/configs/inclusion_cytc.txt\"");
    rep("GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
        "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE");
    rep("\"mode\": \"off\"", "\"mode\": \"ambiguity\"");
    return cfg;
  }

  // Drive a single exploration group end-to-end through FLASHIda — a THIN WRAPPER over the one canonical
  // driver runInterleaved (single_group_only=true): feed MS1 with the engine's own survey id until one forms
  // a group, then drain + feed each MS2 variant back so the winner fires MS3, and project AcqResult->ExplResult.
  struct ExplResult
  {
    bool found_ms3 = false;
    bool found_production_ms2 = false;
    int total_returns = 0;
    int ms3_num_stages = 0;
    int group_commands = 0;
  };

  inline ExplResult driveOneExplorationGroup(FLASHIda* ida,
                                             const std::vector<ScanData>& ms1_scans,
                                             const ScanData& ms2)
  {
    // single_group_only: feed engine-emitted survey MS1 ids (the always-on gate rejects fabricated ones) until
    // the first selects a precursor and forms the exploration group; runInterleaved then drains that group's MS2
    // variants (feeding the ms2 spectrum back so the winner fires MS3) and the MS3, all via the one contract.
    const int budget = 256 + 64 * static_cast<int>(ms1_scans.size() + 1);
    AcqResult a = runInterleaved(ida, ms1_scans, std::vector<ScanData>{ms2}, nullptr, budget,
                                 /*single_group_only=*/true);

    ExplResult r;
    r.group_commands = a.first_group_commands;                       // # variants the forming MS1 pushed
    r.total_returns  = a.ms2_feed_returns;                           // sum of MS2-variant feed returns
    r.found_ms3      = !a.ms3_cmds.empty();                          // an MS3 command was emitted
    r.ms3_num_stages = a.ms3_cmds.empty() ? 0 : a.ms3_cmds[0].num_stages;
    for (const auto& c : a.ms2_cmds)                                 // a production-MS2 (winner) re-acquisition
    {
      std::string d(c.scan_description);
      if (d.size() >= 4 && d[3] != 'E') { r.found_production_ms2 = true; break; }
    }
    return r;
  }

  // Non-draining exploration-group former (C++-side convenience, no C# twin — like single_group_only): feed
  // engine-emitted survey MS1 ids (the always-on gate rejects fabricated ones) until ONE forms a group, then
  // RETURN that survey's processScan count (the # of variant commands it pushed) LEAVING the variants QUEUED
  // for the caller to inspect (driveOneExplorationGroup DRAINS the group; this does not). Bounded idle>=3/max_iters.
  inline int bootstrapExplorationGroup(FLASHIda* ida, const std::vector<ScanData>& scans, int max_iters = 4000)
  {
    const int n_ms1 = static_cast<int>(scans.size());
    int ms1_fed = 0;
    int idle = 0;
    ScanCommand cmd{};
    for (int it = 0; it < max_iters && idle < 3; ++it)
    {
      if (ida->getNextScanCommand(cmd) != 1) break;
      if (cmd.msn_level == 1 && !cmd.is_agc && ms1_fed < n_ms1)
      {
        const ScanData& s = scans[ms1_fed++];
        int n = ida->processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 1,
                                 cmd.scan_description, 0.0, instrumentScanNumberOf(s));
        if (n > 0) return n;  // group created; leave its variants queued for the caller
        idle = 0;
      }
      else { ++idle; }       // AGC / re-survey-after-exhausted -> idle tick
      cmd = ScanCommand{};
    }
    return 0;
  }

  // Project an AcqResult's MS2 commands to (charge, mz, isolation_width) rows (stage-0 of each MS2 command).
  // For asserting which precursors the engine actually acquired — reads the harness result, so the
  // drive still goes through runInterleaved.
  struct AcquisitionRow { int charge; double mz; double width; };
  inline std::vector<AcquisitionRow> ms2AcquisitionRows(const AcqResult& a)
  {
    std::vector<AcquisitionRow> rows;
    for (const auto& c : a.ms2_cmds)
      if (c.num_stages >= 1)
        rows.push_back({c.stages[0].charge_state, c.stages[0].precursor_mz, c.stages[0].isolation_width});
    return rows;
  }

  // ----------------------------------------------------------------------------
  // Common spectrum paths (relative to the CTest working dir, OpenMS/build/)
  // ----------------------------------------------------------------------------
  const std::string FI_MS1_STD   = "../../FlashIDA/test-data/spectra/ms1_standard.txt";
  const std::string FI_MS2_HCD   = "../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt";
  const std::string FI_MS1_CYTC  = "../../FlashIDA/test-data/spectra/ms1_cytc.txt";
  const std::string FI_MS2_CYTC  = "../../FlashIDA/test-data/spectra/ms2_cytc_fresh_scan57.txt";
  const std::string FI_MS2_CYTC149 = "../../FlashIDA/test-data/spectra/ms2_cytc_scan149.txt";
} // anonymous namespace
