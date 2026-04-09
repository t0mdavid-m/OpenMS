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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>

#include <iostream>
#include <stdexcept>
#include <nlohmann/json.hpp>

namespace OpenMS
{

  // Static default level config (selection=None, exploration=None)
  const MSLevelConfig Config::default_level_ = {
    {},                           // scans (empty)
    SelectionMetric::None,        // selection
    10,                           // max_targets
    ExplorationMetric::None,      // exploration
    20.0, 40.0, 5.0,              // ce_min, ce_max, ce_step
    "HCD",                        // exploration_activation
    {},                           // overrides
    10.0                          // tolerance_ppm
  };

  Config::Config(const std::string& json_str)
  {
    // Reject non-JSON input
    if (json_str.empty() || json_str[0] != '{')
    {
      throw std::invalid_argument("Config: input must be JSON (starts with '{')");
    }

    using json = nlohmann::json;
    json config = json::parse(json_str);

    // --- deconvolution section ---
    auto deconv = config.value("deconvolution", json::object());
    targeting_.qscore_threshold = deconv.value("score_threshold", 0.0);
    targeting_.tqscore_threshold = deconv.value("tqscore_threshold", 0.9);
    deconv_.min_charge = deconv.value("min_charge", 4);
    deconv_.max_charge = deconv.value("max_charge", 50);
    deconv_.min_mass = deconv.value("min_mass", 500.0);
    deconv_.max_mass = deconv.value("max_mass", 50000.0);

    // Tolerance values: tol[0] -> MS1, tol[1] -> MS2+
    std::vector<double> tol_values;
    if (deconv.contains("tol") && deconv["tol"].is_array())
    {
      for (const auto& v : deconv["tol"])
        tol_values.push_back(v.get<double>());
    }
    if (tol_values.empty())
      tol_values = {10.0, 10.0};
    if (tol_values.size() == 1)
      tol_values.push_back(tol_values[0]);

    // Store tolerance per level (will be set on MSLevelConfig after level parsing)
    double tol_ms1 = tol_values[0];
    double tol_ms2 = tol_values.size() >= 2 ? tol_values[1] : tol_values[0];

    // Use MS2 tolerance for tag matching
    targeting_.tag_matching_tolerance_ppm = tol_ms2;

    // --- precursor_selection section ---
    auto ps = config.value("precursor_selection", json::object());
    targeting_.rt_window = ps.value("RT_window", 180.0);
    targeting_.mode = ps.value("target_mode", 0);
    targeting_.use_idscore = ps.value("IDScore", false);
    targeting_.consider_all_charges = ps.value("AllCharges", false);
    targeting_.ms3_all_charges = ps.value("MS3AllCharges", false);
    targeting_.hcd_energy = ps.value("HCDEnergy", -1);
    targeting_.strict_inclusion = ps.value("strict_inclusion", false);
    targeting_.tie_threshold = ps.value("tie_threshold", 0.1);

    if (targeting_.mode == 1)
      std::cout << "Inclusion mode: " << (targeting_.strict_inclusion ? "strict" : "non-strict") << "\n";

    // --- tagging section ---
    auto tagging = config.value("tagging", json::object());
    targeting_.min_tag_length = tagging.value("min_tag_length", 3);
    targeting_.max_tag_length = tagging.value("max_tag_length", 8);
    targeting_.max_total_ptm_count = tagging.value("max_ptm_count", 3);
    targeting_.max_flanking_mass_diff = tagging.value("max_flanking_mass_diff", 50000.0);

    // --- files section (paths only; loading stays in FLASHIda) ---
    auto files = config.value("files", json::object());

    if (files.contains("target_logs") && files["target_logs"].is_array())
    {
      for (const auto& f : files["target_logs"])
        targeting_.target_log_files.push_back(f.get<std::string>());
    }
    if (files.contains("fasta") && !files["fasta"].get<std::string>().empty())
      targeting_.fasta_file = files["fasta"].get<std::string>();
    if (files.contains("inclusion_list") && !files["inclusion_list"].get<std::string>().empty())
      targeting_.inclusion_list_file = files["inclusion_list"].get<std::string>();
    if (files.contains("ptm_list") && !files["ptm_list"].get<std::string>().empty())
      targeting_.ptm_list_file = files["ptm_list"].get<std::string>();

    // --- ms3 section ---
    auto ms3 = config.value("ms3", json::object());
    targeting_.ms3_mode = ms3.value("mode", 0);
    targeting_.protein_sequence = ms3.value("protein_sequence", "");
    // ms3 max_per_ms2 -> levels_[3].max_targets (set after selection_strategy parsing)
    int ms3_max_per_ms2 = ms3.value("max_per_ms2", 4);

    // --- conditional_ms2 (derived from ms2 scan count, stored for backward compat) ---
    // conditional_ms2 flag from JSON is now redundant; actual conditional MS2
    // support is derived from levels_[2].scans.size() >= 2.
    // We still read it silently so unknown keys don't cause parse errors.
    (void)config.value("conditional_ms2", false);

    // --- quantification ---
    auto quant = config.value("quantification", json::object());
    quant_.enabled = quant.value("enabled", false);
    quant_.reporter_mz_tol = quant.value("reporter_mz_tol", 0.002);
    quant_.fold_change_threshold = quant.value("fold_change_threshold", 1.4);

    // --- faims ---
    auto faims_section = config.value("faims", json::object());
    if (faims_section.contains("cv_values") && faims_section["cv_values"].is_array())
    {
      for (const auto& v : faims_section["cv_values"])
        faims_.cv_values.push_back(v.get<double>());
    }
    faims_.max_cv_skip = faims_section.value("max_cv_skip", 0);
    faims_.precursor_threshold = faims_section.value("cv_precursor_threshold", 15);
    faims_.enabled = (faims_.cv_values.size() > 1);

    // --- ms_settings: populate levels_[1] (MS1) and levels_[2] (MS2) scan configs ---
    auto ms_settings = config.value("ms_settings", json::object());

    // MS1 scan config -> levels_[1].scans[0]
    auto ms1_json = ms_settings.value("ms1", json::object());
    ScanConfig ms1_scan;
    ms1_scan.analyzer = ms1_json.value("analyzer", "");
    ms1_scan.first_mass = ms1_json.value("first_mass", 0.0);
    ms1_scan.last_mass = ms1_json.value("last_mass", 0.0);
    ms1_scan.resolution = ms1_json.value("resolution", 0);
    ms1_scan.agc_target = ms1_json.value("agc_target", 0);
    ms1_scan.max_it = ms1_json.value("max_it", 0.0);

    // Ensure levels_[1] exists before populating scans
    if (levels_.find(1) == levels_.end())
      levels_[1] = MSLevelConfig{};
    levels_[1].scans.push_back(ms1_scan);

    // MS2 scan configs -> levels_[2].scans[0..N]
    if (ms_settings.contains("ms2") && ms_settings["ms2"].is_array())
    {
      if (levels_.find(2) == levels_.end())
        levels_[2] = MSLevelConfig{};
      for (const auto& m : ms_settings["ms2"])
      {
        ScanConfig ms2_scan;
        ms2_scan.analyzer = m.value("analyzer", "");
        ms2_scan.activation = m.value("activation", "");
        ms2_scan.collision_energy = m.value("collision_energy", 0);
        ms2_scan.resolution = m.value("resolution", 0);
        levels_[2].scans.push_back(ms2_scan);
      }
    }

    // --- scheduling ---
    auto sched = config.value("scheduling", json::object());
    auto ct = sched.value("cycle_time", json::object());
    scheduling_.cycle_time_enabled = ct.value("enabled", false);
    scheduling_.cycle_time_ms = ct.value("value_ms", 60000.0);
    auto to = sched.value("scan_timeout", json::object());
    scheduling_.timeout_enabled = to.value("enabled", false);
    scheduling_.timeout_ms = to.value("value_ms", 30000.0);
    double agc_interval_sec = sched.value("agc_interval_seconds", 30.0);
    scheduling_.agc_interval_ms = static_cast<uint64_t>(agc_interval_sec * 1000.0);

    // --- selection_strategy (required) ---
    if (!config.contains("selection_strategy"))
    {
      throw std::runtime_error("Config: missing required 'selection_strategy' in JSON config");
    }
    const auto& sel_strategy = config["selection_strategy"];
    for (auto it = sel_strategy.begin(); it != sel_strategy.end(); ++it)
    {
      std::string ms_key = it.key();
      if (ms_key.substr(0, 2) == "ms" && ms_key.size() > 2)
      {
        int level_num = std::stoi(ms_key.substr(2));
        const auto& level_obj = it.value();

        // Ensure the level exists in the map
        if (levels_.find(level_num) == levels_.end())
          levels_[level_num] = MSLevelConfig{};
        MSLevelConfig& cfg = levels_[level_num];

        // Selection metric
        std::string sel_str = level_obj.value("selection",
            level_num == 1 ? std::string("qscore") : std::string("intensity"));
        if (sel_str == "intensity") cfg.selection = SelectionMetric::Intensity;
        else if (sel_str == "qscore") cfg.selection = SelectionMetric::QScore;
        else if (sel_str == "none") cfg.selection = SelectionMetric::None;
        else cfg.selection = SelectionMetric::Intensity;

        // Max targets (with aliases)
        cfg.max_targets = level_obj.value("max_targets",
            level_obj.value("max_precursors",
            level_obj.value("max_fragments", 10)));

        // Exploration (optional, MS2+ only; guard against JSON null)
        if (level_obj.contains("exploration") && !level_obj["exploration"].is_null() && level_num > 1)
        {
          const auto& expl_obj = level_obj["exploration"];
          std::string met_str = expl_obj.value("metric", std::string("none"));
          if (met_str == "mass_count") cfg.exploration = ExplorationMetric::MassCount;
          else if (met_str == "remaining_precursor") cfg.exploration = ExplorationMetric::RemainingPrecursor;
          else if (met_str == "fragment_count") cfg.exploration = ExplorationMetric::FragmentCount;
          else cfg.exploration = ExplorationMetric::None;
          cfg.ce_min = expl_obj.value("ce_min", 20.0);
          cfg.ce_max = expl_obj.value("ce_max", 40.0);
          cfg.ce_step = expl_obj.value("ce_step", 5.0);
          cfg.exploration_activation = expl_obj.value("activation", std::string("HCD"));
          if (expl_obj.contains("overrides"))
          {
            const auto& ov_obj = expl_obj["overrides"];
            for (auto ov_it = ov_obj.begin(); ov_it != ov_obj.end(); ++ov_it)
              cfg.overrides[ov_it.key()] = ov_it.value().get<std::string>();
          }
        }
      }
    }

    // Apply ms3 max_per_ms2 to levels_[3].max_targets if level 3 exists
    if (levels_.find(3) != levels_.end())
    {
      // Only override if ms3.max_per_ms2 was set and level 3's max_targets
      // is still the default (10), otherwise selection_strategy takes precedence
      if (ms3_max_per_ms2 != 4 || levels_[3].max_targets == 10)
        levels_[3].max_targets = ms3_max_per_ms2;
    }

    // Set per-level tolerance values
    for (auto& [lvl, cfg] : levels_)
    {
      cfg.tolerance_ppm = (lvl <= 1) ? tol_ms1 : tol_ms2;
    }

    // Compute convenience boolean
    exploration_enabled_ = false;
    for (const auto& [lvl, cfg] : levels_)
    {
      if (cfg.exploration != ExplorationMetric::None)
      {
        exploration_enabled_ = true;
        break;
      }
    }

    // SNR threshold (hardcoded in original parseJSONConfig_)
    targeting_.snr_threshold = 1.0;
  }

  const MSLevelConfig& Config::level(int msn_level) const
  {
    auto it = levels_.find(msn_level);
    return (it != levels_.end()) ? it->second : default_level_;
  }

  bool Config::hasExploration(int msn_level) const
  {
    auto it = levels_.find(msn_level);
    if (it == levels_.end()) return false;
    return it->second.exploration != ExplorationMetric::None;
  }

  int Config::getInt(const std::string& key) const
  {
    if (key == "targeting_mode") return targeting_.mode;
    if (key == "hcd_energy") return targeting_.hcd_energy;
    return 0;
  }

  double Config::getDouble(const std::string& key) const
  {
    if (key == "rt_window") return targeting_.rt_window;
    return 0.0;
  }

} // namespace OpenMS
