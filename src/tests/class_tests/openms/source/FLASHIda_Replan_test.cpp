// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// Survey re-planning: a survey's MS2 decisions stay provisional while they wait in FLASHIda's queue,
// every newer survey at the same CV re-plans them, and only the drain's dequeue commits.
//
// The survey hands the queue its whole plan for one CV in one call (replacePlan), one PLAN UNIT per
// species, and prints `[REPLAN] cv= survey= kept= dropped= added= adjusted=`. These sections drive the
// real engine and read that line and the commands the drain hands out.

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>

#include "FLASHIda_TestAccess.h"
#include "FLASHIda_TestHelpers.h"  // ground-truth harness: ScanData/loadTsvScans/runInterleaved, TSVFile

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace OpenMS;

namespace
{
  /// One `[REPLAN]` line, as the survey printed it.
  struct ReplanLine
  {
    std::string survey;
    int kept = -1, dropped = -1, added = -1, adjusted = -1;
  };

  /// Every `[REPLAN]` line of @p out, in print order.
  std::vector<ReplanLine> replanLines(const std::string& out)
  {
    std::vector<ReplanLine> lines;
    std::istringstream in(out);
    std::string line;
    while (std::getline(in, line))
    {
      if (line.rfind("[REPLAN] ", 0) != 0) continue;
      ReplanLine r;
      std::istringstream tok(line.substr(9));
      std::string kv;
      while (tok >> kv)
      {
        const auto eq = kv.find('=');
        if (eq == std::string::npos) continue;
        const std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
        if (k == "survey") r.survey = v;
        else if (k == "kept") r.kept = std::atoi(v.c_str());
        else if (k == "dropped") r.dropped = std::atoi(v.c_str());
        else if (k == "added") r.added = std::atoi(v.c_str());
        else if (k == "adjusted") r.adjusted = std::atoi(v.c_str());
      }
      lines.push_back(r);
    }
    return lines;
  }

  /// The C++-schema twin of test-data/configs/method_inclusion_charge_set.json (the C++ Config cannot read
  /// the C# schema): inclusion of cytC at the authored charge set 10;13;16, `single`, five precursors.
  std::string inclusionChargeSetConfig(const std::string& log_dir,
                                       const std::string& list = "inclusion_cytc_charges.txt",
                                       int max_precursors = 5, bool strict_inclusion = false)
  {
    return std::string(R"JSON({
  "deconvolution": {
    "score_threshold": 0.0, "tqscore_threshold": 0.0,
    "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
  },
  "flashtnt": { "min_length": 3, "max_length": 8 },
  "faims": { "cv_values": [], "max_cv_skip": 0 },
  "scheduling": {
    "cycle_time": { "enabled": false, "value_ms": 60000 },
    "scan_timeout": { "enabled": false, "value_ms": 30000 },
    "agc_interval_seconds": 9999999
  },
  "files": { "target_logs": [], "fasta": "", "inclusion_list": "../../FlashIDA/test-data/configs/)JSON") + list + R"JSON(", "ptm_list": "" },
  "precursor_selection": {
    "rt_window": 180, "targeting": "inclusion", "precursor_charges": "single",
    "consider_all_charges": false, "strict_inclusion": )JSON" + (strict_inclusion ? "true" : "false") + R"JSON(,
    "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": )JSON" + std::to_string(max_precursors) + R"JSON(
  },
  "characterization": { "mode": "off", "max_targets": 3 },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
    "ms2": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 500000,
             "max_it": 100, "microscans": 4, "activation": "ETD", "reaction_time": 7,
             "reagent_agc_target": 700000, "reagent_max_it": 200 }
  },
  "tagging": {},
  "quantification": { "enabled": false, "reporter_mz_tol": 0.002, "fold_change_threshold": 1.4 },
  "runtime": { "log_dir": ")JSON" + log_dir + R"JSON(" }
})JSON";
  }

  /// Run @p fn with std::cout captured; returns what it printed.
  template <typename Fn>
  std::string captureCout(Fn fn)
  {
    std::ostringstream captured;
    std::streambuf* const real_cout = std::cout.rdbuf(captured.rdbuf());
    try
    {
      fn();
    }
    catch (...)
    {
      std::cout.rdbuf(real_cout);
      throw;
    }
    std::cout.rdbuf(real_cout);
    return captured.str();
  }

  /// One survey's plan as the run shows it: the `[REPLAN]` line's counts, and the nominal masses of the
  /// MS2 commands the survey parented.
  struct SurveyPlan
  {
    int added = -1;
    int dropped = -1;
    size_t commands = 0;
    std::set<int> nominals;
  };

  /// Drive @p cfg over cytC and pair every `[REPLAN]` line with the survey's MS2 commands. The drive is
  /// depth 1, so every MS2 a survey planned is dequeued before the next survey: the survey's plan is
  /// exactly the MS2 commands whose parent it is.
  std::vector<SurveyPlan> plansOf(const std::string& cfg)
  {
    auto scans = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
    AcqResult a;
    const std::string out = captureCout([&] {
      FLASHIda ida(const_cast<char*>(cfg.c_str()));
      a = runInterleaved(&ida, scans, std::vector<ScanData>{}, nullptr, 4000);
    });
    std::map<std::string, std::vector<const ScanCommand*>> by_survey;
    for (const auto& c : a.ms2_cmds) by_survey[std::string(c.parent_scan_id)].push_back(&c);

    std::vector<SurveyPlan> plans;
    for (const auto& r : replanLines(out))
    {
      SurveyPlan p;
      p.added = r.added;
      p.dropped = r.dropped;
      const auto it = by_survey.find(r.survey);
      if (it != by_survey.end())
      {
        p.commands = it->second.size();
        for (const ScanCommand* c : it->second) p.nominals.insert(SpectralDeconvolution::getNominalMass(c->mono_mass));
      }
      plans.push_back(p);
    }
    return plans;
  }

  /// Split @p s on tabs, keeping empty fields.
  std::vector<std::string> splitTabs(const std::string& s)
  {
    std::vector<std::string> out;
    std::string cur;
    for (char ch : s)
    {
      if (ch == '\t')
      {
        out.push_back(cur);
        cur.clear();
      }
      else
      {
        cur += ch;
      }
    }
    out.push_back(cur);
    return out;
  }

  /// scan_commands.tsv rows with the one column that carries wall-clock time (enqueue_ts) removed.
  std::vector<std::vector<std::string>> stableCommandRows(const std::string& path, std::vector<std::string>& headers)
  {
    const TSVFile tsv = TSVFile::parse(path);
    headers.clear();
    std::vector<int> keep;
    for (size_t i = 0; i < tsv.headers.size(); ++i)
    {
      if (tsv.headers[i] == "enqueue_ts") continue;
      keep.push_back(static_cast<int>(i));
      headers.push_back(tsv.headers[i]);
    }
    std::vector<std::vector<std::string>> rows;
    for (const auto& row : tsv.rows)
    {
      std::vector<std::string> cells;
      for (int i : keep) cells.push_back(i < static_cast<int>(row.size()) ? row[i] : std::string());
      rows.push_back(cells);
    }
    return rows;
  }

  /// Cells agree: byte-identical, or -- token by token over the ';' and ',' sub-cell axes -- identical
  /// strings or numbers within the golden comparer's relative tolerance. gcc and MSVC differ in FP
  /// reduction order (CLAUDE.md), so a printed last digit is not a behaviour; a different command is.
  bool sameCell(const std::string& a, const std::string& b)
  {
    if (a == b) return true;
    auto tokens = [](const std::string& s) {
      std::vector<std::string> t;
      std::string cur;
      for (char ch : s)
      {
        if (ch == ';' || ch == ',')
        {
          t.push_back(cur);
          cur.clear();
        }
        else
        {
          cur += ch;
        }
      }
      t.push_back(cur);
      return t;
    };
    const auto ta = tokens(a), tb = tokens(b);
    if (ta.size() != tb.size()) return false;
    for (size_t i = 0; i < ta.size(); ++i)
    {
      if (ta[i] == tb[i]) continue;
      if (ta[i].empty() || tb[i].empty()) return false;
      char* ea = nullptr;
      char* eb = nullptr;
      const double x = std::strtod(ta[i].c_str(), &ea);
      const double y = std::strtod(tb[i].c_str(), &eb);
      if (*ea != '\0' || *eb != '\0') return false;
      if (std::abs(x - y) > std::max(1e-6, 1e-3 * std::max(std::abs(x), std::abs(y)))) return false;
    }
    return true;
  }

  // scan_commands.tsv of the depth-1 drive below, captured at the Push 1 tip (OpenMS 3e5e58437b + Task 7b/7,
  // before any survey called replacePlan) without enqueue_ts. Rows in dequeue order.
  const std::vector<std::string> kDepth1PushOneTip = {
    "!!\"\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!\"S\t0\t0\t500.0000\t2000."
      "0000",
    "!!#\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!#S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!$\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!$S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!%\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!%S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!&\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!&S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!'\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!'S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!(\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!(S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!)\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!)S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!*\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!*S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!+\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!+S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!,\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!,S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!-\trecording\t2\t!!,\t1\t2\t2063.6058\t4\t517.034\t1.05037\t0.476161\t0.860556\t1.43069\t0.860556\t1.18869\t1\tHCD\t29"
      "\t29\t0\t0\t0\t0.23633\t5459.27\t5459.27\t\t0\t\t!!-R2.063606k@4\t0\t0\t0.0000\t0.0000",
    "!!.\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!.S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!/\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!/S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!0\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!0S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!1\trecording\t2\t!!0\t2\t2\t2277.2541\t4\t570.572\t1.29744\t0.771507\t0.93761\t4.96672\t0.93761\t4.96672\t1\tHCD\t29\t"
      "29\t0\t0\t0\t2.99553\t5213.45\t5213.45\t\t0\t\t!!1R2.277254k@4\t0\t0\t0.0000\t0.0000",
    "!!2\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!2S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!3\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!3S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!4\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!4S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!5\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!5S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!6\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!6S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!7\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!7S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!8\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!8S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!9\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!9S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!:\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!:S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!;\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!;S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!<\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!<S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!=\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!=S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!>\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!>S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!?\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!?S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!@\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!@S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!A\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!AS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!B\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!BS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!C\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!CS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!D\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!DS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!E\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!ES\t0\t0\t500.0000\t2000.00"
      "00",
    "!!F\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!FS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!G\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!GS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!H\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!HS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!I\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!IS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!J\trecording\t2\t!!I\t3\t2\t5315.1289\t9\t591.911\t1.02839\t0.350165\t0.858807\t3.53716\t0.858807\t3.53716\t1\tHCD\t29"
      "\t29\t0\t0\t0\t3.67734\t5876.7\t5876.7\t\t0\t\t!!JR5.315129k@9\t0\t0\t0.0000\t0.0000",
    "!!K\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!KS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!L\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!LS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!M\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!MS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!N\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!NS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!O\trecording\t2\t!!N\t4\t2\t2256.1789\t4\t565.177\t1.05073\t0.419769\t0.853377\t1.30749\t0.853377\t0.528162\t1\tHCD\t2"
      "9\t29\t0\t0\t0\t0.161604\t3336.41\t3336.41\t\t0\t\t!!OR2.256179k@4\t0\t0\t0.0000\t0.0000",
    "!!P\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!PS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!Q\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!QS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!R\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!RS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!S\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!SS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!T\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!TS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!U\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!US\t0\t0\t500.0000\t2000.00"
      "00",
    "!!V\trecording\t2\t!!U\t5\t2\t4297.1774\t7\t615.249\t0.941357\t0.393986\t0.57079\t2.8727\t0.881239\t4.42826\t1\tHCD\t29"
      "\t29\t0\t0\t0\t3.00522\t3821.83\t5851.61\t\t0\t\t!!VR4.297177k@7\t0\t0\t0.0000\t0.0000",
    "!!W\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!WS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!X\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!XS\t0\t0\t500.0000\t2000.00"
      "00",
    "!!Y\trecording\t2\t!!X\t6\t2\t2135.7138\t4\t535.061\t1.0536\t0.3612\t0.872528\t12.5045\t0.872528\t0.352325\t1\tHCD\t29\t"
      "29\t0\t0\t0\t2.79475\t3519.72\t3519.72\t\t0\t\t!!YR2.135714k@4\t0\t0\t0.0000\t0.0000",
    "!!Z\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!ZS\t0\t0\t500.0000\t2000.00"
      "00",
    "!![\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!![S\t0\t0\t500.0000\t2000.00"
      "00",
    "!!\\\tsurvey\t1\t\t0\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tnone\t0\t0\t0\t0\t0\t0\t0\t0\t\t0\t\t!!\\S\t0\t0\t500.0000\t2000."
      "0000",
  };
  // ---- The bought R (Task 10): a quantification verdict buys an identification R that waits ----------------

  /// The C++-schema twin of test-data/configs/method_quant_welch.json (@p charges "separate", the quant_welch
  /// golden's method) or of method_quant_multiplexed.json ("multiplexed"): one species per survey, two to three
  /// charges, a Q per charge or one Q for the set, and the verdict buys an ETD identification R.
  std::string quantTwinConfig(const std::string& charges)
  {
    const bool separate = charges == "separate";
    return std::string(R"JSON({
  "deconvolution": {
    "score_threshold": 0.0, "tqscore_threshold": 0.9,
    "min_charge": 4, "max_charge": 50, "min_mass": 500, "max_mass": 50000, "tol": [10, 10, 10]
  },
  "flashtnt": { "min_length": 3, "max_length": 8 },
  "conditional_ms2": false,
  "faims": { "cv_values": [], "max_cv_skip": 0 },
  "scheduling": {
    "cycle_time": { "enabled": false, "value_ms": 60000 },
    "scan_timeout": { "enabled": false, "value_ms": 30000 },
    "agc_interval_seconds": 9999999
  },
  "files": { "target_logs": [], "fasta": "", "inclusion_list": "", "ptm_list": "" },
  "precursor_selection": {
    "rt_window": 180, "targeting": "none", "consider_all_charges": false, "strict_inclusion": false,
    "tie_threshold": 0.1, "rank_by": "qscore", "max_precursors": 1,
    "precursor_charges": ")JSON") + charges + R"JSON(", "min_charge_states": 2, "max_charge_states": 3
  },
  "characterization": { "mode": "off", "max_targets": 3 },
  "ms_settings": {
    "ms1": { "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000, "resolution": 120000, "agc_target": 800000, "max_it": 246 },
    "ms2": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 500000,
             "max_it": 100, "microscans": 4, "activation": "ETD", "reaction_time": 7,
             "reagent_agc_target": 700000, "reagent_max_it": 200 },
    "ms2_quant": { "analyzer": "Orbitrap", "first_mass": 100, "last_mass": 2000, "resolution": 120000, "agc_target": 500000,
                   "max_it": 100, "microscans": 4, "activation": "HCD", "collision_energy": 30 }
  },
  "tagging": {},
  "quantification": {
    "enabled": true, "labelling": "tmt6plex", "reporter_mz_tol": 0.002,)JSON"
           + (separate ? std::string(R"JSON(
    "test": "welch_and_fold_change", "p_value_threshold": 0.05, "fold_change_threshold": 1.5,
    "empty_condition": "incomplete", "aggregate": "median", "min_usable_scans": 1,)JSON")
                       : std::string(R"JSON(
    "fold_change_threshold": 1.4,)JSON"))
           + R"JSON(
    "conditions": [ { "name": "control", "channels": ["126", "127", "128"] },
                    { "name": "treated", "channels": ["129", "130", "131"] } ]
  }
})JSON";
  }

  /// A bought R, caught while it waits: the command as queued, the charges it was registered with, the
  /// survey that decided its Q, and every charge the Qs fed back to the engine isolated (sorted, unique).
  struct WaitingR
  {
    bool bought = false;
    ScanCommand r{};
    std::vector<int> covered;
    int q_survey = 0;
    std::vector<int> q_charges;
  };

  /// Feed @p a as the engine's first survey and return its 'Q' scans the two real iTMT spectra of one species by
  /// charge (odd -> z23, even -> z22, as the quant_welch golden does) until the verdict buys an R. Stops BEFORE
  /// the drain could take the R: it is read off the queue, never dequeued.
  WaitingR buyR(FLASHIda& ida, const ScanData& a, const ScanData& q_odd, const ScanData& q_even)
  {
    WaitingR out;
    ScanCommandQueue& q = FLASHIdaTestAccess::queue(ida);
    ScanCommand cmd{};
    ida.getNextScanCommand(cmd); // the idle survey
    if (cmd.msn_level != 1) return out;
    out.q_survey = cmd.scan_id;
    ida.processScan(a.mzs.data(), a.ints.data(), (int)a.mzs.size(), a.rt, 1, cmd.scan_description, cmd.faims_cv,
                    instrumentScanNumberOf(a));
    for (int it = 0; it < 20; ++it)
    {
      for (const ScanCommand& c : FLASHIdaTestAccess::queuedCommands(q))
      {
        const auto covered = FLASHIdaTestAccess::adjustableResolvedCharges(q, c.scan_id);
        if (!covered) continue;
        out.bought = true;
        out.r = c;
        out.covered = *covered;
        std::sort(out.q_charges.begin(), out.q_charges.end());
        out.q_charges.erase(std::unique(out.q_charges.begin(), out.q_charges.end()), out.q_charges.end());
        return out;
      }
      cmd = ScanCommand{};
      ida.getNextScanCommand(cmd);
      if (cmd.msn_level != 2 || cmd.scan_description[3] != 'Q') return out; // the plan ran out before a buy
      out.q_charges.push_back(std::abs(cmd.stages[0].charge_state));
      for (int j = 0; j < cmd.stage0_notch_count; ++j) out.q_charges.push_back(std::abs(cmd.notches[j].charge_state));
      const ScanData& s = std::abs(cmd.stages[0].charge_state) % 2 == 1 ? q_odd : q_even;
      ida.processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), s.rt, 2, cmd.scan_description, 0.0,
                      instrumentScanNumberOf(s));
    }
    return out;
  }

  /// A survey that reaches the engine while the R still waits -- a cycle-time or CV-transition MS1 the drain
  /// sent ahead of it. Minted and registered as the drain would; returns its tracking id.
  int feedSurvey(FLASHIda& ida, const ScanData& s, double rt)
  {
    ScanCommandQueue& q = FLASHIdaTestAccess::queue(ida);
    ScanCommand ms1 = q.makeMS1();
    ms1.scan_id = q.nextTrackingId();
    ScanCommandQueue::writeMs1Description(ms1, ScanRole::Survey);
    q.registerPending(ms1);
    ida.processScan(s.mzs.data(), s.ints.data(), (int)s.mzs.size(), rt, 1, ms1.scan_description, 0.0, 0);
    return ms1.scan_id;
  }

  /// @p s without the isotope envelope of the species of mono mass @p mass at each charge of @p charges: every peak
  /// between the monoisotopic m/z and 25 Da above it, padded by 0.3 Th. The other charges' windows lie tens of Th
  /// away, so the species stays in the spectrum at every charge not named.
  ScanData withoutCharges(const ScanData& s, double mass, const std::vector<int>& charges)
  {
    ScanData out = s;
    out.mzs.clear();
    out.ints.clear();
    for (size_t i = 0; i < s.mzs.size(); ++i)
    {
      bool inside = false;
      for (int z : charges)
      {
        const double lo = mass / z + 1.007276 - 0.3;
        const double hi = (mass + 25.0) / z + 1.007276 + 0.3;
        if (s.mzs[i] >= lo && s.mzs[i] <= hi) inside = true;
      }
      if (inside) continue;
      out.mzs.push_back(s.mzs[i]);
      out.ints.push_back(s.ints[i]);
    }
    return out;
  }

  /// The species' charges on deconvolution @p d as charge_selection "snr" (the default) ranks them, computed here
  /// from the PeakGroups rather than through the engine: within ±1 nominal mass, isolatable, charge SNR >= 1,
  /// restricted to @p allowed when non-empty; SNR descending, then intensity, then the lower charge.
  std::vector<int> rankedCharges(const DeconvolvedSpectrum& d, int nominal, const std::vector<int>& allowed)
  {
    struct C { int z; double snr; double intensity; };
    std::vector<C> cs;
    for (const PeakGroup& pg : d)
    {
      if (std::abs(SpectralDeconvolution::getNominalMass(pg.getMonoMass()) - nominal) > 1) continue;
      const auto [min_z, max_z] = pg.getAbsChargeRange();
      for (int z = min_z; z <= max_z; ++z)
      {
        const auto [lo, hi] = pg.getMzRange(z);
        if (hi <= lo || pg.getChargeSNR(z) < 1.0) continue;
        if (!allowed.empty() && std::find(allowed.begin(), allowed.end(), z) == allowed.end()) continue;
        cs.push_back({z, pg.getChargeSNR(z), pg.getChargeIntensity(z)});
      }
    }
    std::stable_sort(cs.begin(), cs.end(), [](const C& a, const C& b) {
      if (a.snr != b.snr) return a.snr > b.snr;
      if (a.intensity != b.intensity) return a.intensity > b.intensity;
      return a.z < b.z;
    });
    std::vector<int> out;
    for (const C& c : cs)
      if (std::find(out.begin(), out.end(), c.z) == out.end()) out.push_back(c.z);
    return out;
  }

  /// The anchor first, then the stage-0 notches: every charge @p c isolates.
  std::vector<int> isolatedCharges(const ScanCommand& c)
  {
    std::vector<int> zs = {std::abs(c.stages[0].charge_state)};
    for (int j = 0; j < c.stage0_notch_count; ++j) zs.push_back(std::abs(c.notches[j].charge_state));
    return zs;
  }

  /// The queued command with tracking id @p id, read without dequeuing; scan_id 0 when it is not queued.
  ScanCommand queued(FLASHIda& ida, int id)
  {
    for (const ScanCommand& c : FLASHIdaTestAccess::queuedCommands(FLASHIdaTestAccess::queue(ida)))
      if (c.scan_id == id) return c;
    return ScanCommand{};
  }
} // namespace

START_TEST(FLASHIda_Replan, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// Two PeakGroups of one species in one survey (ADR-0036 siblings) are ONE plan unit: the unit is keyed by
// nominal mass, so a survey adds exactly as many units as it selected distinct nominal masses.
START_SECTION(siblings_form_one_unit)
{
  // method_inclusion_charge_set: cytC at 10;13;16. Its winner carries all three named charges in every
  // survey of this fixture that acquires it, so no sibling is ever admitted -- the assertion holds, but it
  // cannot tell one unit per species from one unit per PeakGroup.
  const auto plans = plansOf(inclusionChargeSetConfig(freshLogDir("replan_siblings")));
  TEST_EQUAL(plans.empty(), false)
  for (const auto& p : plans)
  {
    TEST_EQUAL(p.added, static_cast<int>(p.nominals.size()))
    TEST_EQUAL(p.dropped, 0)
  }

  // So the same assertion again on the WIDE row 12..18 (FLASHIda_ChargeModes_test CM-09's config), whose
  // envelope arrives split: there a survey admits a sibling, and two commands of one nominal mass in one
  // survey must still be one unit.
  const auto split = plansOf(inclusionChargeSetConfig(freshLogDir("replan_siblings_split"), "inclusion_cytc_charges_wide.txt",
                                                      /*max_precursors=*/2, /*strict_inclusion=*/true));
  TEST_EQUAL(split.empty(), false)
  int surveys_with_siblings = 0;
  for (const auto& p : split)
  {
    TEST_EQUAL(p.added, static_cast<int>(p.nominals.size()))
    TEST_EQUAL(p.dropped, 0)
    if (p.commands > p.nominals.size()) ++surveys_with_siblings;
  }
  TEST_EQUAL(surveys_with_siblings > 0, true) // non-vacuity: a sibling was admitted
}
END_SECTION

// At depth 1 every MS2 a survey planned is dequeued before the next survey arrives, so no survey ever meets
// a waiting unit of its CV, and re-planning must not change a single command. The drive is
// FLASHIda_Logging_test's (method_default's C++ twin, ms1_standard + ms2_hcd_fragment), compared with a
// capture of the same drive taken at the Push 1 tip.
START_SECTION(depth1_drive_is_unchanged_without_a_stale_survey)
{
  const std::string log_dir = freshLogDir("replan_depth1");
  const std::string cfg = buildJsonWithLogDir(log_dir);
  auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_standard.txt");
  auto ms2 = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_hcd_fragment.txt");
  ABORT_IF(ms1.empty() || ms2.empty())
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    runFullCycle(&ida, ms1, ms2);
  });
  const auto lines = replanLines(out);
  TEST_EQUAL(lines.empty(), false)
  for (const auto& r : lines)
  {
    TEST_EQUAL(r.dropped, 0)
    TEST_EQUAL(r.kept, 0)
  }

  std::vector<std::string> headers;
  const auto rows = stableCommandRows(log_dir + "/scan_commands.tsv", headers);
  TEST_EQUAL(rows.empty(), false)
  TEST_EQUAL(rows.size(), kDepth1PushOneTip.size())
  ABORT_IF(rows.size() != kDepth1PushOneTip.size())
  for (size_t r = 0; r < rows.size(); ++r)
  {
    const std::vector<std::string> expected = splitTabs(kDepth1PushOneTip[r]);
    TEST_EQUAL(rows[r].size(), expected.size())
    if (rows[r].size() != expected.size()) continue;
    for (size_t c = 0; c < expected.size(); ++c)
    {
      if (sameCell(rows[r][c], expected[c])) continue;
      std::cerr << "row " << r << " column " << headers[c] << std::endl;
      TEST_EQUAL(rows[r][c], expected[c])
    }
  }
}
END_SECTION

// ---- The bought R is adjustable, never droppable --------------------------------------------------------------
//
// ms1_cytc.txt[8] fed as the first survey: a `separate` quantification group of cytC (three Qs, decided on the
// two real iTMT spectra) buys an identification R at its purest window, z13 -- while charge_selection ("snr")
// ranks z14 first on that survey. The R is caught WAITING, a newer survey reaches the engine ahead of it, and the
// drain then takes it. Every expected charge is ranked by rankedCharges from the PeakGroups, not by the engine.
// The R may move only among the charges its group's Q scans covered (owner 2026-10-23): the anchors and notches of
// the Qs that returned, registered with the R.

// A newer survey that ranks another Q's charge first moves the R there: the R row's charge is the newer survey's
// choice, that survey decided it, its parent is still the winning Q, and its queue time restarts.
START_SECTION(bought_r_adjusted_to_better_charge)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("separate");
  const ScanData& a = ms1[8];

  WaitingR w;
  int b_id = 0, expected = 0;
  ScanCommand after{}, sent{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, a, q_odd[0], q_even[0]);
    if (!w.bought) return;
    std::this_thread::sleep_for(std::chrono::milliseconds(5)); // a queue-time reset is then a strictly later stamp
    b_id = feedSurvey(ida, a, a.rt + 2.0);                      // the same spectrum, two seconds on
    const auto ranked = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), SpectralDeconvolution::getNominalMass(w.r.mono_mass),
                                      w.covered);
    if (!ranked.empty()) expected = ranked.front();
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
    ida.getNextScanCommand(sent);
  });
  ABORT_IF(!w.bought)
  TEST_EQUAL(std::abs(w.r.stages[0].charge_state), 13)
  TEST_EQUAL(w.covered == w.q_charges, true) // registered with exactly the charges the group's Qs isolated
  TEST_EQUAL(w.covered.size() > 1, true)
  TEST_EQUAL(expected, 14) // non-vacuity: the newer survey prefers another charge, and a Q covered it
  TEST_EQUAL(std::find(w.covered.begin(), w.covered.end(), expected) != w.covered.end(), true)

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(std::abs(after.stages[0].charge_state), expected)
  const std::string desc(after.scan_description);
  TEST_EQUAL(desc.substr(desc.rfind('@') + 1), std::to_string(expected))
  TEST_EQUAL(std::string(after.parent_scan_id), std::string(w.r.parent_scan_id)) // lineage: still the winning Q
  TEST_EQUAL(after.enqueue_timestamp_ms > w.r.enqueue_timestamp_ms, true)
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, b_id)
  TEST_EQUAL(decision->isolation_purity >= 0.0 && decision->isolation_purity <= 1.0, true)
  const auto lines = replanLines(out);
  ABORT_IF(lines.size() != 2)
  TEST_EQUAL(lines.back().survey, ScanCommandQueue::encode(b_id))
  TEST_EQUAL(lines.back().adjusted, 1)

  // ...and that is the command the drain hands out.
  TEST_EQUAL(sent.scan_id, w.r.scan_id)
  TEST_EQUAL(std::abs(sent.stages[0].charge_state), expected)
}
END_SECTION

// A newer survey that ranks the R's own charge first re-confirms it (owner 2026-10-22): the charge stays, but the
// newer survey is now the one that decided it, and the queue time restarts. The newer survey is A without whichever
// charges outrank the R's on it, so the R's charge leads among those the Qs covered.
START_SECTION(bought_r_reconfirmed_same_charge_resets_decision)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("separate");
  const ScanData& a = ms1[8];

  WaitingR w;
  int b_id = 0, expected = 0;
  std::vector<int> ranked_on_a;
  ScanCommand after{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, a, q_odd[0], q_even[0]);
    if (!w.bought) return;
    const int nominal = SpectralDeconvolution::getNominalMass(w.r.mono_mass);
    const int own = std::abs(w.r.stages[0].charge_state);
    ranked_on_a = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, w.covered);
    // Take out whatever charge outranks the R's, one at a time, until the R's leads: removing an envelope moves the
    // others' SNR, so the charges to remove are found on a scratch engine's deconvolution of each candidate.
    std::vector<int> removed;
    ScanData b = a;
    for (int i = 0; i < 6; ++i)
    {
      std::vector<int> ranked;
      {
        FLASHIda scratch(const_cast<char*>(cfg.c_str()));
        feedSurvey(scratch, b, a.rt + 2.0);
        ranked = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(scratch), nominal, w.covered);
      }
      if (ranked.empty() || ranked.front() == own) break;
      removed.push_back(ranked.front());
      b = withoutCharges(a, w.r.mono_mass, removed);
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    b_id = feedSurvey(ida, b, a.rt + 2.0);
    const auto ranked = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, w.covered);
    if (!ranked.empty()) expected = ranked.front();
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
  });
  ABORT_IF(!w.bought || ranked_on_a.empty())
  const int bought_z = std::abs(w.r.stages[0].charge_state);
  TEST_EQUAL(ranked_on_a.front() != bought_z, true) // non-vacuity: on the Q's survey another charge led
  TEST_EQUAL(expected, bought_z)                    // ...on the newer survey the R's own leads

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(std::abs(after.stages[0].charge_state), bought_z)
  TEST_EQUAL(std::string(after.scan_description), std::string(w.r.scan_description))
  TEST_EQUAL(after.enqueue_timestamp_ms > w.r.enqueue_timestamp_ms, true)
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, b_id)
  const auto lines = replanLines(out);
  ABORT_IF(lines.empty())
  TEST_EQUAL(lines.back().adjusted, 1)
}
END_SECTION

// A newer survey that resolves the species only at charges no Q covered cannot move the R: the newer survey is A
// without every covered charge, so it still resolves cytC, at charges the group never measured. The R goes out as
// decided, at its charge, still attributed to the Q's survey.
START_SECTION(bought_r_keeps_charge_when_no_covered_charge_resolves)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("separate");
  const ScanData& a = ms1[8];

  WaitingR w;
  std::vector<int> ranked_on_b, ranked_on_both;
  ScanCommand after{}, sent{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, a, q_odd[0], q_even[0]);
    if (!w.bought) return;
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    feedSurvey(ida, withoutCharges(a, w.r.mono_mass, w.covered), a.rt + 2.0);
    const int nominal = SpectralDeconvolution::getNominalMass(w.r.mono_mass);
    ranked_on_b = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, {});
    ranked_on_both = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, w.covered);
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
    ida.getNextScanCommand(sent);
  });
  ABORT_IF(!w.bought || ranked_on_b.empty())
  // Non-vacuity: the newer survey resolves cytC, at no charge a Q covered.
  TEST_EQUAL(std::find(w.covered.begin(), w.covered.end(), ranked_on_b.front()) == w.covered.end(), true)
  TEST_EQUAL(ranked_on_both.empty(), true)

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(std::string(after.scan_description), std::string(w.r.scan_description))
  TEST_EQUAL(after.stages[0].charge_state, w.r.stages[0].charge_state)
  TEST_REAL_SIMILAR(after.stages[0].precursor_mz, w.r.stages[0].precursor_mz)
  TEST_EQUAL(after.enqueue_timestamp_ms, w.r.enqueue_timestamp_ms) // nothing re-chose it
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, w.q_survey)
  const auto lines = replanLines(out);
  ABORT_IF(lines.empty())
  TEST_EQUAL(lines.back().adjusted, 0)
  TEST_EQUAL(sent.scan_id, w.r.scan_id)
  TEST_EQUAL(sent.stages[0].charge_state, w.r.stages[0].charge_state)
}
END_SECTION

// A multiplexed bought R re-chooses its whole SET (owner 2026-10-22) within the charges its Q covered -- the Q's
// anchor and notches (owner 2026-10-23): anchor and notches become those charges in the newer survey's order, and
// its scan_description names the new anchor. The newer survey is A without the last charge of the R's set: it
// still resolves cytC at charges the Q never isolated, and none of them may fill the freed slot.
START_SECTION(bought_r_multiplexed_set_rechosen)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("multiplexed");
  const ScanData& a = ms1[8];

  WaitingR w;
  int b_id = 0;
  std::vector<int> expected, ranked_on_b;
  ScanCommand after{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, a, q_odd[0], q_even[0]);
    if (!w.bought) return;
    b_id = feedSurvey(ida, withoutCharges(a, w.r.mono_mass, {isolatedCharges(w.r).back()}), a.rt + 2.0);
    const int nominal = SpectralDeconvolution::getNominalMass(w.r.mono_mass);
    ranked_on_b = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, {});
    expected = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, w.covered);
    if (expected.size() > 3) expected.resize(3);
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
  });
  ABORT_IF(!w.bought)
  TEST_EQUAL(isolatedCharges(w.r).size(), 3)
  std::vector<int> isolated_sorted = isolatedCharges(w.r);
  std::sort(isolated_sorted.begin(), isolated_sorted.end());
  TEST_EQUAL(w.covered == isolated_sorted, true) // the one Q's anchor and notches, which the R inherited
  TEST_EQUAL(expected.size(), 2)                 // the freed slot stays empty ...
  TEST_EQUAL(ranked_on_b.size() > expected.size(), true) // ... though the newer survey resolves more charges
  TEST_EQUAL(expected != isolatedCharges(w.r), true)     // non-vacuity: the newer survey's set differs

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(isolatedCharges(after) == expected, true)
  const std::string desc(after.scan_description);
  TEST_EQUAL(desc.substr(desc.rfind('@') + 1), std::to_string(expected.front()))
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, b_id)
  const auto lines = replanLines(out);
  ABORT_IF(lines.empty())
  TEST_EQUAL(lines.back().adjusted, 1)
}
END_SECTION

// The covered set is a FILTER, not a preference: a newer survey that ranks a charge no Q covered first and a
// covered one after it moves the R to the covered one. The Q's survey is A without A's own top charge, so that
// charge is not covered; the newer survey is A without the R's charge, so the top charge leads again and the R
// cannot stay where it is.
START_SECTION(bought_r_moves_only_among_covered_charges)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("separate");
  const ScanData& a = ms1[8];
  const double cytc = 12351.3; // the species ms1_cytc.txt[8]'s verdict buys an R for (the sections above)

  int top_on_a = 0;
  {
    FLASHIda scratch(const_cast<char*>(cfg.c_str()));
    captureCout([&] { feedSurvey(scratch, a, a.rt); });
    const auto ranked = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(scratch), SpectralDeconvolution::getNominalMass(cytc), {});
    if (!ranked.empty()) top_on_a = ranked.front();
  }
  ABORT_IF(top_on_a == 0)

  WaitingR w;
  int b_id = 0, expected = 0;
  std::vector<int> ranked_on_b;
  ScanCommand after{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, withoutCharges(a, cytc, {top_on_a}), q_odd[0], q_even[0]);
    if (!w.bought) return;
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    b_id = feedSurvey(ida, withoutCharges(a, w.r.mono_mass, {std::abs(w.r.stages[0].charge_state)}), a.rt + 2.0);
    const int nominal = SpectralDeconvolution::getNominalMass(w.r.mono_mass);
    ranked_on_b = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, {});
    const auto ranked = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), nominal, w.covered);
    if (!ranked.empty()) expected = ranked.front();
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
  });
  ABORT_IF(!w.bought || ranked_on_b.empty() || expected == 0)
  const int bought_z = std::abs(w.r.stages[0].charge_state);
  // Non-vacuity: the newer survey's best charge is one no Q covered, and a covered charge other than the R's own
  // follows it.
  TEST_EQUAL(std::find(w.covered.begin(), w.covered.end(), top_on_a) == w.covered.end(), true)
  TEST_EQUAL(std::find(w.covered.begin(), w.covered.end(), ranked_on_b.front()) == w.covered.end(), true)
  TEST_EQUAL(expected != bought_z, true)

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(std::abs(after.stages[0].charge_state), expected)
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, b_id)
  const auto lines = replanLines(out);
  ABORT_IF(lines.empty())
  TEST_EQUAL(lines.back().adjusted, 1)
}
END_SECTION

// A multiplexed R the newer survey could re-isolate at fewer than min_charge_states (2) of its covered charges is
// kept as decided (owner 2026-10-23, ADR-0040's floor): its set, its description and its decision survey stay. The
// newer survey is A without every covered charge but one.
START_SECTION(bought_r_multiplexed_below_floor_is_kept)
{
  const auto ms1 = loadTsvScans("../../FlashIDA/test-data/spectra/ms1_cytc.txt");
  const auto q_odd = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z23_18295.txt");
  const auto q_even = loadTsvScans("../../FlashIDA/test-data/spectra/ms2_quant_itmt_z22_18295.txt");
  ABORT_IF(ms1.size() < 9 || q_odd.empty() || q_even.empty())
  const std::string cfg = quantTwinConfig("multiplexed");
  const ScanData& a = ms1[8];

  WaitingR w;
  std::vector<int> ranked_on_b;
  ScanCommand after{};
  std::optional<DecisionInfo> decision;
  const std::string out = captureCout([&] {
    FLASHIda ida(const_cast<char*>(cfg.c_str()));
    w = buyR(ida, a, q_odd[0], q_even[0]);
    if (!w.bought || w.covered.empty()) return;
    std::this_thread::sleep_for(std::chrono::milliseconds(5));
    const std::vector<int> others(w.covered.begin() + 1, w.covered.end());
    feedSurvey(ida, withoutCharges(a, w.r.mono_mass, others), a.rt + 2.0);
    ranked_on_b = rankedCharges(FLASHIdaTestAccess::surveyDeconvolution(ida), SpectralDeconvolution::getNominalMass(w.r.mono_mass),
                                w.covered);
    after = queued(ida, w.r.scan_id);
    decision = FLASHIdaTestAccess::decisionInfo(FLASHIdaTestAccess::queue(ida), w.r.scan_id);
  });
  ABORT_IF(!w.bought)
  TEST_EQUAL(isolatedCharges(w.r).size(), 3)
  TEST_EQUAL(ranked_on_b.size(), 1) // non-vacuity: one covered charge resolves, so there was something to move to

  TEST_EQUAL(after.scan_id, w.r.scan_id)
  TEST_EQUAL(isolatedCharges(after) == isolatedCharges(w.r), true)
  TEST_EQUAL(std::string(after.scan_description), std::string(w.r.scan_description))
  TEST_EQUAL(after.enqueue_timestamp_ms, w.r.enqueue_timestamp_ms) // nothing re-chose it
  ABORT_IF(!decision)
  TEST_EQUAL(decision->decision_survey_id, w.q_survey)
  const auto lines = replanLines(out);
  ABORT_IF(lines.empty())
  TEST_EQUAL(lines.back().adjusted, 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
