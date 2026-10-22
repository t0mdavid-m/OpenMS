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
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
