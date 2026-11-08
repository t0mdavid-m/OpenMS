// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/CandidateAdmission.h>   // sameSpeciesMass
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

/// Per-proteoform MS3 re-planning (ADR-0064): the plain data and the pure rules the queue's registry and the
/// tracker's plan gate share. Header-only and free of `Config`, in the shape of `PlanUnit.h`.
///
/// With identification.characterization.per_proteoform on, a proteoform has ONE live MS3 plan: a newer
/// identification of it replaces its still-unsent MS3s; a sent MS3 commits its fragment (cleavage + charge + the
/// MS2 intensity it was chosen at); a sent fragment is re-acquired only at >= 2x that intensity; at most
/// max_targets distinct cleavages per proteoform per run.
namespace OpenMS
{
  /// (ion type, ion index) -- the same type as ProteoformTracker's FragmentKey. An unassigned raw mass (exhaustive
  /// 'u', a discriminating mass) is ("u", its nominal mass).
  using Ms3Cleavage = std::pair<std::string, int>;

  /// What identifies a proteoform across selections: the protein, its region (0-based, exclusive end, -1 already
  /// normalised to 0 / P) and the COMMANDED precursor mass, matched by sameProteoform.
  struct Ms3ProteoformKey
  {
    std::string accession;
    int region_start = 0;
    int region_end = 0;
    double commanded_mass = 0.0;
  };

  /// Same protein, same region, and masses equal within @p ms1_tol_ppm of the larger -- or the same species seen
  /// twice by sameSpeciesMass (1..3 isotope steps, a x2/x3/:2/:3 harmonic), the rule the close-out uses.
  inline bool sameProteoform(const Ms3ProteoformKey& a, const Ms3ProteoformKey& b, double ms1_tol_ppm)
  {
    if (a.accession != b.accession || a.region_start != b.region_start || a.region_end != b.region_end) return false;
    if (a.commanded_mass <= 0.0 || b.commanded_mass <= 0.0) return false;
    const double tol = ms1_tol_ppm * 1e-6 * std::max(a.commanded_mass, b.commanded_mass);
    return std::abs(a.commanded_mass - b.commanded_mass) <= tol
           || sameSpeciesMass(a.commanded_mass, b.commanded_mass, ms1_tol_ppm);
  }

  /// One replaceable unit of a proteoform's live MS3 plan: one MS3 command, or the N CE variants of one MS3
  /// exploration group. It commits as a whole on its first command's dequeue.
  struct Ms3PlanItem
  {
    Ms3Cleavage cleavage;
    std::vector<std::pair<int, double>> charges;  ///< (charge, MS2 intensity) it acquires: the anchor + stage-1 notches
    std::vector<ScanCommand> commands;
    int sweep_group_id = 0;        ///< the Exploration group to drop when withdrawn (0 = a single MS3)
    int owner_precursor_id = 0;    ///< the model that planned it
    double fragment_mass = 0.0;    ///< a discriminating record's lookup key
    int nominal = 0;               ///< the exhaustive / discriminating dispatch stamp it set
    bool committed = false;        ///< queue-owned
    int commands_left = 0;         ///< queue-owned
  };

  /// What the analysis thread plans against: the proteoform's sent fragments (cleavage -> charge -> the highest
  /// intensity it was sent at). @c generation counts commits, so a swap can name the ones made after the snapshot.
  struct Ms3RegistrySnapshot
  {
    bool found = false;
    uint64_t generation = 0;
    std::map<Ms3Cleavage, std::map<int, double>> sent;
  };

  /// The one eligibility rule for a candidate (cleavage, charge, intensity), used at planning and again at the
  /// swap: a cleavage never sent is admitted; a charge already sent needs >= 2x its sent intensity; another charge
  /// of a sent cleavage needs >= 2x the cleavage's highest sent intensity. (The cleavage cap is counted separately.)
  inline bool ms3Admits(const Ms3RegistrySnapshot& snap, const Ms3Cleavage& c, int z, double intensity)
  {
    const auto it = snap.sent.find(c);
    if (it == snap.sent.end() || it->second.empty()) return true;
    const auto zt = it->second.find(z);
    double ref = 0.0;
    if (zt != it->second.end()) ref = zt->second;
    else
      for (const auto& kv : it->second) ref = std::max(ref, kv.second);
    return intensity >= 2.0 * ref;
  }

  /// How many cleavages never sent a proteoform may still get this run.
  inline int ms3NewCleavageRoom(const Ms3RegistrySnapshot& snap, int max_targets)
  {
    return std::max(0, max_targets - static_cast<int>(snap.sent.size()));
  }

  /// What one replaceMs3Plan did. The caller undoes the withdrawn and reconciled items after the queue lock.
  struct Ms3ReplanResult
  {
    int entry_id = 0;
    std::vector<int> inserted_scan_ids;
    std::vector<int> withdrawn_scan_ids;
    std::vector<int> reconciled_scan_ids;
    std::vector<Ms3PlanItem> withdrawn;     ///< unsent items of the previous plan
    std::vector<Ms3PlanItem> reconciled;    ///< new items refused at the swap (a commit in the gap, the cap)
    std::vector<Ms3Cleavage> committed_in_gap;
    std::set<std::pair<int, int>> committed_owner_nominals;   ///< (precursor id, nominal) of committed items
  };
} // namespace OpenMS
