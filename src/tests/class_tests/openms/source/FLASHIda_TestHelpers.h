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

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
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

  // Push all scans through processScan (ms_level=1), return total command count
  inline int pushAllScans(FLASHIda* ida, const std::vector<ScanData>& scans)
  {
    int total = 0;
    for (const auto& scan : scans)
    {
      total += ida->processScan(scan.mzs.data(), scan.ints.data(),
                                (int)scan.mzs.size(), scan.rt, 1,
                                ("scan_" + scan.scan_id).c_str());
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
        std::vector<std::string> cols;
        std::istringstream iss(line);
        std::string col;
        while (std::getline(iss, col, '\t')) cols.push_back(col);

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

  inline CycleResult runFullCycle(FLASHIda* ida,
                                  const std::vector<ScanData>& ms1_scans,
                                  const std::vector<ScanData>& ms2_scans)
  {
    CycleResult result;

    pushAllScans(ida, ms1_scans);

    ScanCommand cmd{};
    while (ida->getNextScanCommand(cmd) > 0)
    {
      result.total_dequeued++;
      if (cmd.msn_level == 2) result.ms2_cmds.push_back(cmd);
      else if (cmd.msn_level == 1) result.ms1_cmds.push_back(cmd);
      if (cmd.is_agc) break;  // idle-cycle AGC => real queue drained
    }

    for (const auto& ms2_cmd : result.ms2_cmds)
    {
      if (!ms2_scans.empty())
        ida->processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                         (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 2,
                         ms2_cmd.scan_description);
    }

    while (ida->getNextScanCommand(cmd) > 0)
    {
      result.total_dequeued++;
      if (cmd.msn_level == 3) result.ms3_cmds.push_back(cmd);
      if (cmd.is_agc) break;
    }

    for (const auto& ms3_cmd : result.ms3_cmds)
    {
      if (!ms2_scans.empty())
        ida->processScan(ms2_scans[0].mzs.data(), ms2_scans[0].ints.data(),
                         (int)ms2_scans[0].mzs.size(), ms2_scans[0].rt, 3,
                         ms3_cmd.scan_description);
    }

    return result;
  }

  // ----------------------------------------------------------------------------
  // JSON config builder with runtime log paths (extends FLASHIda_Logging_test's
  // builder with an identification_log_path slot).  enable_ms3 emits the proven
  // inclusion-pinned cytC MS3 recipe (target_mode=1 + inclusion_cytc.txt + the
  // M-starting proteoform + ms3 stage); callers feed ms1_cytc + ms2_cytc_fresh_scan57.
  // ----------------------------------------------------------------------------
  inline std::string buildJsonWithRuntime(const std::string& ida_log_path,
                                          const std::string& commands_path,
                                          const std::string& results_path,
                                          bool enable_ms3 = false,
                                          const std::string& identification_path = "",
                                          const std::string& ms2_activation = "HCD")
  {
    std::string target_mode_val   = enable_ms3 ? "1" : "0";
    std::string inclusion_list_val = enable_ms3 ? "../../FlashIDA/test-data/configs/inclusion_cytc.txt" : "";
    std::string ms2_selection      = enable_ms3 ? "\"intensity\"" : "\"none\"";
    std::string ms3_settings = enable_ms3
      ? R"(,
        "ms3": [ { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 35, "resolution": 120000 } ])"
      : "";
    std::string ms3_block = enable_ms3
      ? R"("ms3": { "protein_sequence": "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE" },)"
      : "";
    std::string ms3_selection = enable_ms3 ? "\"intensity\"" : "\"none\"";
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
        "RT_window": 180, "target_mode": )" << target_mode_val << R"(,
        "IDScore": false, "AllCharges": false,
        "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1
      },
      "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
      "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
      "faims": { "cv_values": [-50], "max_cv_skip": 0 },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
        "ms2": [
          { "analyzer": "Orbitrap", "activation": ")" << ms2_activation << R"(", "collision_energy": 29)" << ms2_rt << R"(, "resolution": 120000 }
        ])" << ms3_settings << R"(
      },
      "scheduling": {
        "cycle_time": { "enabled": false, "value_ms": 60000 },
        "scan_timeout": { "enabled": false, "value_ms": 30000 },
        "agc_interval_seconds": 999999
      },
      "exploration": { "enabled": false, "max_depth": 1, "max_variants": 5 },
      )" << ms3_block << R"(
      "files": { "target_logs": [], "fasta": "", "inclusion_list": ")" << inclusion_list_val << R"(", "ptm_list": "" },
      "selection_strategy": {
        "ms1": { "selection": "qscore", "max_targets": 5 },
        "ms2": { "selection": )" << ms2_selection << R"( },
        "ms3": { "selection": )" << ms3_selection << R"( }
      },
      "runtime": {
        "ida_log_path": ")" << ida_log_path << R"(",
        "scan_commands_path": ")" << commands_path << R"(",
        "scan_results_path": ")" << results_path << R"(",
        "identification_log_path": ")" << identification_path << R"("
      }
    })";
    return oss.str();
  }

  // Insert a "runtime" block (the four log paths) into an arbitrary config JSON that
  // lacks one, immediately before its closing brace.  Used to drive the exploration
  // configs (which have no runtime section) and capture their log files.
  inline std::string injectRuntime(std::string cfg,
                                   const std::string& commands_path,
                                   const std::string& results_path,
                                   const std::string& identification_path = "",
                                   const std::string& ida_log_path = "")
  {
    auto last = cfg.rfind('}');
    if (last == std::string::npos) return cfg;
    std::ostringstream rt;
    rt << ", \"runtime\": {"
       << "\"ida_log_path\": \"" << ida_log_path << "\", "
       << "\"scan_commands_path\": \"" << commands_path << "\", "
       << "\"scan_results_path\": \"" << results_path << "\", "
       << "\"identification_log_path\": \"" << identification_path << "\"}";
    cfg.insert(last, rt.str());
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
      "precursor_selection": { "RT_window": 180, "target_mode": 0, "IDScore": false, "AllCharges": false, "HCDEnergy": 29, "strict_inclusion": false, "tie_threshold": 0.1 },
      "tagging": { "min_tag_length": 3, "max_tag_length": 8, "max_ptm_count": 3, "max_flanking_mass_diff": 50000 },
      "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
      "faims": { "cv_values": [-50], "max_cv_skip": 0, "cv_precursor_threshold": 15 },
      "ms_settings": {
        "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
        "ms2": [ { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 } ],
        "ms3": [ { "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 35, "resolution": 120000 } ]
      },
      "scheduling": { "cycle_time": { "enabled": false, "value_ms": 60000 }, "scan_timeout": { "enabled": false, "value_ms": 30000 }, "agc_interval_seconds": 999999 },
      "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
      "ms3": { "protein_sequence": "GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE" },
      "conditional_ms2": false,
      "selection_strategy": {
        "ms1": { "selection": "qscore", "max_targets": 3 },
        "ms2": { "selection": "intensity", "max_targets": 3, "exploration": { "metric": "mass_count", "ce_min": 20.0, "ce_max": 40.0, "ce_step": 5.0 } },
        "ms3": { "selection": "none" }
      }
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
    rep("\"target_mode\": 0", "\"target_mode\": 1");
    rep("\"inclusion_list\": \"\"",
        "\"inclusion_list\": \"../../FlashIDA/test-data/configs/inclusion_cytc.txt\"");
    rep("GDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFSYTDANKNKGITWGEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE",
        "MGDVEKGKKIFVQKCAQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGFTYTDANKNKGITWKEETLMEYLENPKKYIPGTKMIFAGIKKKTEREDLIAYLKKATNE");
    rep("\"selection\": \"none\"", "\"selection\": \"intensity\"");
    return cfg;
  }

  // Drive a single exploration group end-to-end through FLASHIda (verbatim from
  // FLASHIda_exploration_test): push MS1 until the first scan creates a group, then
  // drain + feed each MS2 variant back so the winner fires MS3.
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
    ExplResult r;
    for (const auto& s : ms1_scans)
    {
      int n = ida->processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 1,
                               ("scan_" + s.scan_id).c_str());
      if (n > 0) { r.group_commands = n; break; }
    }
    if (r.group_commands == 0) return r;

    int idle = 0;
    for (int i = 0; i < 100 && !r.found_ms3; ++i)
    {
      ScanCommand next{};
      if (ida->getNextScanCommand(next) != 1) break;
      if (next.is_agc || next.msn_level == 1) { if (++idle > 3) break; continue; }
      idle = 0;
      if (next.msn_level == 3) { r.found_ms3 = true; r.ms3_num_stages = next.num_stages; break; }
      if (next.msn_level == 2)
      {
        std::string d(next.scan_description);
        if (d.size() >= 4 && d[3] != 'E') r.found_production_ms2 = true;
        r.total_returns += ida->processScan(ms2.mzs.data(), ms2.ints.data(),
                                            (int)ms2.mzs.size(), ms2.rt, 2, next.scan_description);
      }
    }
    return r;
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
