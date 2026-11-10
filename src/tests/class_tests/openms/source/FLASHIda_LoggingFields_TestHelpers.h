// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// The file-local helpers of FLASHIda_LoggingFields_test and FLASHIda_LoggingFields2_test (one suite split over two
// executables for ctest -j wall time, 2026-11-10). Moved here unchanged from FLASHIda_LoggingFields_test.cpp; an
// anonymous namespace, because each test exe is its own translation unit.

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include "FLASHIda_TestHelpers.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace OpenMS;

namespace
{
  // tracking_id -> ms_level, from a parsed scan_commands.tsv (used to classify result rows,
  // which carry no ms_level column: a result row's id is in this map iff it is an MS2/MS3 scan;
  // absent => MS1 survey input).
  std::map<std::string, int> commandLevels(const TSVFile& cmds)
  {
    std::map<std::string, int> m;
    int id_col = cmds.colIndex("tracking_id");
    int lvl_col = cmds.colIndex("ms_level");
    if (id_col < 0 || lvl_col < 0) return m;
    for (const auto& r : cmds.rows)
      if (id_col < (int)r.size() && lvl_col < (int)r.size())
        m[r[id_col]] = std::atoi(r[lvl_col].c_str());
    return m;
  }

  bool ionTypeOk(const std::string& s)
  {
    if (s.size() != 1) return false;
    char c = s[0];
    return c == 'a' || c == 'b' || c == 'c' || c == 'x' || c == 'y' || c == 'z';
  }

  // Decode the trailing precursor ion of an MS3 scan_description per the descriptor contract:
  //   {id}R{mass}k@{charge}{ion_type}{ion_index}  (the 'k' is part of the mass token, before '@').
  // Find the LAST '@', skip the run of fragment-charge digits after it, require an ion_type char in
  // {a,b,c,x,y,z} (via ionTypeOk), then an all-digit ion_index >= 1. Returns false for the no-ion form
  // {id}R{mass}k@{charge} (nothing after the charge digits) -- that branch is tolerated by callers.
  // Kept byte-for-byte identical to decodeTrailingIonKey in FLASHIda_TestHelpers.h.
  bool decodeTrailingIon(const std::string& d, char& t, int& idx)
  {
    std::size_t at = d.rfind('@');
    if (at == std::string::npos) return false;
    std::size_t i = at + 1;
    while (i < d.size() && std::isdigit(static_cast<unsigned char>(d[i]))) ++i;  // skip fragment charge digits
    if (i >= d.size() || !ionTypeOk(std::string(1, d[i]))) return false;          // no-ion branch -> false
    t = d[i++];
    std::string n = d.substr(i);
    if (n.empty() || n.find_first_not_of("0123456789") != std::string::npos) return false;  // non-truncated index
    idx = std::atoi(n.c_str());
    return idx >= 1;
  }

  // Build ion -> [ScanData...] from real MS3 fixtures ms3_cytc_<ion>_scan<N>.txt in a spectra dir.
  // The ion key (e.g. "b44") is the substring between "ms3_cytc_" and "_scan" (per the fixture naming
  // contract). loadTsvScans each matching file and group by ion. Returns empty if the dir is absent or
  // holds no matching fixtures (caller then skips cleanly -- never fabricates data).
  std::map<std::string, std::vector<ScanData>> buildMs3IonManifest(const std::string& dir)
  {
    std::map<std::string, std::vector<ScanData>> manifest;
    std::error_code ec;
    if (!std::filesystem::is_directory(dir, ec)) return manifest;
    const std::string prefix = "ms3_cytc_";
    for (std::filesystem::directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec))
    {
      if (!it->is_regular_file(ec)) continue;
      std::string name = it->path().filename().string();
      if (name.size() < prefix.size() || name.compare(0, prefix.size(), prefix) != 0) continue;
      // The ion token sits between "ms3_cytc_" and the LAST "_scan"; require both present and a non-empty ion.
      std::size_t scan_pos = name.rfind("_scan");
      if (scan_pos == std::string::npos || scan_pos <= prefix.size()) continue;  // e.g. legacy ms3_cytc_scan341.txt
      std::string ion = name.substr(prefix.size(), scan_pos - prefix.size());
      if (ion.empty()) continue;
      auto scans = loadTsvScans(it->path().string());
      if (scans.empty()) continue;
      auto& bucket = manifest[ion];
      for (auto& s : scans) bucket.push_back(std::move(s));
    }
    return manifest;
  }

  // Strict ProForma validator (per the plan glossary): every RESIDUE char (alphabetic, outside [...]
  // mod groups and the () region delimiters) must be UPPERCASE A-Z (no lowercase residues); modifications
  // are bracketed [..]; an ambiguity region '(' RESIDUES ')' must be immediately followed by a '[' mod;
  // brackets balanced; >=1 residue. e.g. GDVEK, PEP[+10.0000]TIDE, (GDVEKGKK)[+42.0157]IFVQK valid;
  // gdvek[+42] (lowercase), GDVEK(IFVQ)KIFVQK (region not followed by a mod) invalid.
  bool isProForma(const std::string& s)
  {
    if (s.empty()) return false;
    bool in_mod = false;             // inside [...]
    bool expect_mod = false;         // just closed a ')', a '[' mod must follow
    bool any_residue = false;
    for (char ch : s)
    {
      if (in_mod) { if (ch == ']') in_mod = false; continue; }
      if (expect_mod) { if (ch != '[') return false; expect_mod = false; in_mod = true; continue; }
      if (ch == '[') { in_mod = true; continue; }
      if (ch == '(') continue;
      if (ch == ')') { expect_mod = true; continue; }
      if (std::isalpha(static_cast<unsigned char>(ch)))
      {
        if (std::isupper(static_cast<unsigned char>(ch)) == 0) return false;  // lowercase residue
        any_residue = true;
      }
    }
    return any_residue && !in_mod && !expect_mod;
  }

  // Parse every PTM's (start,end) 1-based RESIDUE range from a ProForma proteoform, left-to-right.
  //   Localized  RES[+m]         -> [p,p]   (p = residue index the '[' immediately follows)
  //   Ambiguous  (RES..RES)[+m]  -> [s,e]   (s = first residue in parens, e = last)
  // Residues are the uppercase A-Z outside [..] mod content and () region delimiters; bracket bodies are
  // skipped. Mods appear in residue order and narrowing never adds/drops/reorders a mod, so index k of an
  // identification proteoform and index k of the same tracking_id's scan_results proteoform are the SAME
  // modification -- comparable as (start,end) pairs. Used by §T9/§F3 to assert identification ⊆ scan_results.
  std::vector<std::pair<int,int>> parseModRanges(const std::string& s)
  {
    std::vector<std::pair<int,int>> ranges;
    int residues = 0;              // residues read so far (1-based index of the last residue)
    int region_start = -1;         // >0 while inside (...)
    int pending_end = -1;          // end residue captured at ')', awaiting the '[' mod
    bool in_mod = false;           // inside [...]
    bool just_closed_region = false;
    for (char ch : s)
    {
      if (in_mod) { if (ch == ']') in_mod = false; continue; }
      if (ch == '[')
      {
        in_mod = true;
        if (just_closed_region) { ranges.emplace_back(region_start, pending_end); just_closed_region = false; region_start = -1; }
        else                    { ranges.emplace_back(residues, residues); }   // localized at current residue
        continue;
      }
      if (ch == '(') { region_start = residues + 1; continue; }   // next residue is the region's first
      if (ch == ')') { pending_end = residues; just_closed_region = true; continue; }
      if (ch >= 'A' && ch <= 'Z') ++residues;
    }
    return ranges;
  }

  const std::string CYTC_MS1 = FI_MS1_CYTC;
  const std::string CYTC_MS2 = FI_MS2_CYTC;
}
