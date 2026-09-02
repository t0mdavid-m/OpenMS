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

// For the isobaric labelling vocabulary and the channel-name table, so an authored condition
// channel resolves to an ordinal HERE, at load. Quantification.h includes Config.h, not the other
// way round, so this is a one-way dependency and not a cycle.
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Quantification.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <set>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

namespace
{
  // Strict schema: throw naming every key in `obj` that is not in `allowed`. Dynamic objects
  // (e.g. exploration.overrides) are simply not validated (never passed here).
  void rejectUnknownKeys(const nlohmann::json& obj, const std::set<std::string>& allowed,
                         const std::string& path)
  {
    if (!obj.is_object()) return;
    std::vector<std::string> bad;
    for (auto it = obj.begin(); it != obj.end(); ++it)
      if (allowed.find(it.key()) == allowed.end()) bad.push_back(it.key());
    if (bad.empty()) return;
    std::string joined;
    for (size_t i = 0; i < bad.size(); ++i) { if (i) joined += ", "; joined += bad[i]; }
    throw std::invalid_argument(
        "Config: unknown key(s) in '" + path + "': " + joined +
        ". Keys are case-sensitive snake_case; see FlashIDA/test-data/config_schema_reference.json.");
  }

  // Shared parse for the two charge-acquisition keys (precursor_selection.precursor_charges and
  // characterization.fragment_charges). One helper because the value sets are identical and must
  // stay identical -- two hand-written if/else chains is how "Coverage" once silently meant
  // ambiguity. Hard-rejected, not defaulted, per ADR-0007.
  OpenMS::ChargeAcquisitionMode parseChargeMode(const nlohmann::json& obj, const std::string& key,
                                                const std::string& path)
  {
    const std::string v = obj.value(key, std::string("single"));
    if (v == "single")      return OpenMS::ChargeAcquisitionMode::Single;
    if (v == "separate")    return OpenMS::ChargeAcquisitionMode::Separate;
    if (v == "multiplexed") return OpenMS::ChargeAcquisitionMode::Multiplexed;
    throw std::invalid_argument(
        "Config: " + path + "." + key + " must be one of \"single\", \"separate\", \"multiplexed\"; "
        "got \"" + v + "\" (values are case-sensitive).");
  }

  // The authored spelling of a characterization mode, for error messages that quote the user's
  // config back at them.
  //
  // A switch over every enumerator with no `default:`, deliberately: this replaced a two-way
  // ternary (`mode == Coverage ? "coverage" : "ambiguity"`) which, the moment a third on-value
  // existed, printed "ambiguity" for an exhaustive config -- a factually wrong mode name inside a
  // message telling the user what their config said. With the switch, the next value is a
  // compiler warning rather than a wrong string.
  const char* characterizationModeName(OpenMS::CharacterizationMode m)
  {
    switch (m)
    {
      case OpenMS::CharacterizationMode::Off: return "off";
      case OpenMS::CharacterizationMode::Ambiguity: return "ambiguity";
      case OpenMS::CharacterizationMode::Coverage: return "coverage";
      case OpenMS::CharacterizationMode::Exhaustive: return "exhaustive";
    }
    return "off";  // unreachable for a valid enumerator; keeps return-path checks quiet
  }

  // The authored spelling of a quantification objective (ADR-0039), for the [CONFIG-WARN] line that
  // quotes the user's config back at them. A switch over every enumerator with NO `default:`, for
  // the same reason as characterizationModeName above: the next value must be a compiler warning
  // rather than a wrong mode name inside a message telling the user what their config said.
  const char* quantIdentifyName(OpenMS::QuantIdentify i)
  {
    switch (i)
    {
      case OpenMS::QuantIdentify::Differential: return "differential";
      case OpenMS::QuantIdentify::Quantified: return "quantified";
      case OpenMS::QuantIdentify::All: return "all";
      case OpenMS::QuantIdentify::None: return "none";
    }
    return "differential";  // unreachable for a valid enumerator; keeps return-path checks quiet
  }

  // One lenient scan-config allowlist for every scan object (ms1/ms2/ms3/follow_up_scan): the union
  // of MS1- and MS2/MS3-level keys. Rejects non-schema scan keys such as the removed 'IsolationMode'.
  const std::set<std::string> kScanKeys = {
    "analyzer", "activation", "collision_energy", "resolution", "agc_target", "max_it",
    "first_mass", "last_mass", "microscans", "data_type", "rf_lens", "source_cid",
    "source_cid_scaling", "scan_rate", "reaction_time", "reagent_max_it", "reagent_agc_target"};

  // Read EVERY key kScanKeys admits, so validation and parsing share one source of truth. When the
  // two lists drifted apart, nine keys passed validation and were then silently discarded --
  // notably reaction_time on a follow_up_scan, which made an ETD follow-up unconfigurable.
  //
  // A ScanConfig fully determines its scan's instrument parameters (ADR-0009); an unset value means
  // "use the instrument method default", never "inherit from another scan".
  void parseScanConfig(const nlohmann::json& j, OpenMS::ScanConfig& sc,
                       const std::string& analyzer_default)
  {
    // is_string() guard, not .value(): a present-but-null string (the generated
    // config_schema_reference.json carries "data_type": null) makes .value() throw type_error.302.
    auto str = [&](const char* k, const std::string& d) {
      auto it = j.find(k);
      return (it != j.end() && it->is_string()) ? it->get<std::string>() : d;
    };

    sc.analyzer           = str("analyzer", analyzer_default);
    sc.activation         = str("activation", "");
    sc.data_type          = str("data_type", "");
    sc.scan_rate          = str("scan_rate", "");
    sc.collision_energy   = j.value("collision_energy", 0);
    sc.resolution         = j.value("resolution", 0);
    sc.agc_target         = j.value("agc_target", 0);
    sc.max_it             = j.value("max_it", 0.0);
    sc.first_mass         = j.value("first_mass", 0.0);
    sc.last_mass          = j.value("last_mass", 0.0);
    sc.microscans         = j.value("microscans", 0);
    sc.rf_lens            = j.value("rf_lens", 0.0);
    sc.source_cid         = j.value("source_cid", 0.0);
    sc.source_cid_scaling = j.value("source_cid_scaling", 0.0);
    sc.reaction_time      = j.value("reaction_time", 0.0);
    sc.reagent_max_it     = j.value("reagent_max_it", 0.0);
    sc.reagent_agc_target = j.value("reagent_agc_target", 0);
  }

  // --- scan-name resolution -------------------------------------------------------------
  // File-local, not members: none of these touches Config state, and their signatures need
  // nlohmann::json -- which Config.h deliberately does not know about, since the public
  // constructor takes raw JSON as a std::string precisely to keep the library an
  // implementation detail. Declaring them in the header broke that and would not compile.
  bool isValidScanName(const std::string& name)
  {
    if (name.empty() || name.size() > 32) return false;
    if (name[0] < 'a' || name[0] > 'z') return false;
    for (char ch : name)
      if (!((ch >= 'a' && ch <= 'z') || (ch >= '0' && ch <= '9') || ch == '_')) return false;
    return true;
  }

  std::string knownScanNames(const std::map<std::string, OpenMS::ScanConfig>& additional_ms2)
  {
    if (additional_ms2.empty())
      return "ms_settings.additional_ms2 defines no scan configs.";
    std::string s = "Defined in ms_settings.additional_ms2:";
    for (const auto& kv : additional_ms2) s += " " + kv.first;
    return s + ".";
  }

  void resolveFollowUp_(const nlohmann::json& config, const std::string& section,
                                const std::map<std::string, OpenMS::ScanConfig>& additional_ms2,
                                OpenMS::ScanConfig& out, bool& has_follow_up, std::string& out_name)
  {
    const auto sec = config.value(section, nlohmann::json::object());
    if (!sec.contains("follow_up_scan") || sec["follow_up_scan"].is_null())
      return;

    // A follow-up is just another MS2, so it names an additional_ms2 entry rather than repeating a
    // 17-key block. An object here is someone carrying the old shape forward.
    if (sec["follow_up_scan"].is_object())
      throw std::invalid_argument(
          "Config: " + section + ".follow_up_scan is no longer an inline scan object. It is now the "
          "NAME of an ms_settings.additional_ms2 entry, e.g. \"" + section + "_follow_up\".");
    if (!sec["follow_up_scan"].is_string())
      throw std::invalid_argument("Config: " + section + ".follow_up_scan must be a scan-config name.");

    const std::string name = sec["follow_up_scan"].get<std::string>();
    if (name.empty()) return;

    auto found = additional_ms2.find(name);
    if (found == additional_ms2.end())
      throw std::invalid_argument(
          "Config: " + section + ".follow_up_scan references unknown MS2 scan config '" + name
          + "'. " + knownScanNames(additional_ms2));

    out = found->second;
    has_follow_up = true;
    out_name = name;
  }
}

namespace OpenMS
{

  // Activation-coupled parameters: a parameter that is meaningful only for particular activation
  // types, and must therefore travel with the activation whenever it is set (CONTEXT.md).
  //
  // Declared in Config.h and defined at namespace scope rather than in this file's anonymous
  // namespace, because Exploration::buildVariants_ asks the same question -- which axis of a sweep
  // varies, and which one that activation's baseline zeroes. It used to re-inline these literals,
  // so the answer had two definitions that could disagree without anything noticing.
  bool needsCollisionEnergy(const std::string& act)
  {
    return act == "HCD" || act == "CID" || act == "EThcD";
  }
  bool needsReactionTime(const std::string& act)
  {
    return act == "ETD" || act == "EThcD";
  }

  // Static default level config (selection=None, exploration=None)
  const MSLevelConfig Config::default_level_ = {
    {},                           // scans (empty)
    SelectionMetric::None,        // selection
    10,                           // max_targets
    0,                            // min_charge (no filter)
    ExplorationMetric::None,      // exploration
    20.0, 40.0, 5.0,              // ce_min, ce_max, ce_step
    0.0, 0.0, 1.0,                // rt_min, rt_max, rt_step
    {},                           // activations
    {},                           // overrides
    10.0,                         // tolerance_ppm
    10.0,                         // exploration_tolerance_ppm
    0.1                           // remaining_precursor_target
  };

  const std::set<std::string>& Config::scanKeys()
  {
    return kScanKeys;
  }

  void ScanConfig::applyOverrides(const std::unordered_map<std::string, std::string>& overrides)
  {
    for (const auto& [key, val] : overrides)
    {
      if (key == "analyzer") analyzer = val;
      else if (key == "activation") activation = val;
      else if (key == "collision_energy") collision_energy = static_cast<int>(std::stod(val));
      else if (key == "resolution") resolution = static_cast<int>(std::stod(val));
      else if (key == "agc_target") agc_target = static_cast<int>(std::stod(val));
      else if (key == "first_mass") first_mass = std::stod(val);
      else if (key == "last_mass") last_mass = std::stod(val);
      else if (key == "max_it") max_it = std::stod(val);
      else if (key == "microscans") microscans = static_cast<int>(std::stod(val));
      else if (key == "rf_lens") rf_lens = std::stod(val);
      else if (key == "source_cid") source_cid = std::stod(val);
      else if (key == "source_cid_scaling") source_cid_scaling = std::stod(val);
      else if (key == "data_type") data_type = val;
      else if (key == "scan_rate") scan_rate = val;
      else if (key == "reaction_time") reaction_time = std::stod(val);
      else if (key == "reagent_max_it") reagent_max_it = std::stod(val);
      else if (key == "reagent_agc_target") reagent_agc_target = static_cast<int>(std::stod(val));
    }
  }

  Config::Config(const std::string& json_str)
  {
    // Reject non-JSON input
    if (json_str.empty() || json_str[0] != '{')
    {
      throw std::invalid_argument("Config: input must be JSON (starts with '{')");
    }

    using json = nlohmann::json;
    json config = json::parse(json_str);

    // --- strict schema: reject unknown keys (top-level + every section below) ---
    if (config.contains("ms3"))
      throw std::invalid_argument(
          "Config: 'ms3' is no longer a top-level section. "
          "MS3 is configured under 'characterization' (mode/protein_sequence/max_targets) with its "
          "scan parameters in 'ms_settings.ms3'.");
    // Dedicated migration error rather than the generic unknown-key message, which would say only
    // "unknown key 'selection_strategy'" and leave the reader to find all seven destinations.
    if (config.contains("selection_strategy"))
      throw std::invalid_argument(
          "Config: 'selection_strategy' has been removed. Its keys moved to the two decision "
          "sections:\n"
          "  ms1.selection    -> precursor_selection.rank_by\n"
          "  ms1.max_targets  -> precursor_selection.max_precursors   (it is the MS2 count)\n"
          "  ms1.min_charge   -> precursor_selection.min_precursor_charge\n"
          "  ms2.exploration  -> precursor_selection.exploration\n"
          "  ms2.max_targets  -> characterization.max_targets         (it is the MS3 budget)\n"
          "  ms2.min_charge   -> characterization.min_fragment_charge\n"
          "  ms3.exploration  -> characterization.exploration\n"
          "  ms2.selection and ms3.selection are replaced by characterization.mode "
          "(off|ambiguity|coverage|exhaustive); ms3.max_targets and ms3.min_charge were never read and are "
          "deleted.");
    rejectUnknownKeys(config,
        {"global", "deconvolution", "precursor_selection", "flashtnt", "tagging", "conditional_ms2",
         "quantification", "faims", "ms_settings", "scheduling",
         "characterization", "files", "runtime"},
        "(root)");
    rejectUnknownKeys(config.value("global", json::object()),
        {"method_name", "method_description", "duration"}, "global");

    // --- deconvolution section ---
    auto deconv = config.value("deconvolution", json::object());
    rejectUnknownKeys(deconv,
        {"score_threshold", "tqscore_threshold", "min_charge", "max_charge", "min_mass", "max_mass", "tol"},
        "deconvolution");
    targeting_.qscore_threshold = deconv.value("score_threshold", 0.0);
    targeting_.tqscore_threshold = deconv.value("tqscore_threshold", 0.9);
    deconv_.min_charge = deconv.value("min_charge", 4);
    deconv_.max_charge = deconv.value("max_charge", 50);
    deconv_.min_mass = deconv.value("min_mass", 500.0);
    deconv_.max_mass = deconv.value("max_mass", 50000.0);

    // Tolerance values: one entry per MS level, indexed by level-1
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

    // Use MS2 tolerance for tag matching (index 1, guaranteed to exist)
    targeting_.tag_matching_tolerance_ppm = tol_values.size() >= 2 ? tol_values[1] : tol_values[0];

    // --- precursor_selection section: WHICH intact species do we fragment? ---
    // Absorbs what used to be selection_strategy.ms1. Keys are snake_case throughout and named for
    // what they PRODUCE: max_precursors is the MS2 count. rank_by/max_precursors/min_precursor_charge
    // are read into levels_[1] by applyCharacterizationMode_() once every section has been parsed.
    auto ps = config.value("precursor_selection", json::object());
    // Dedicated migration error rather than the generic unknown-key message: the flag's only
    // user-visible effect was a fan-out that is now a named acquisition mode, and a reader hitting
    // "unknown key" would have no way to know that (ADR-0021).
    if (ps.contains("charge_based_exclusion"))
      throw std::invalid_argument(
          "Config: precursor_selection.charge_based_exclusion has been removed (ADR-0021). It was a "
          "developer flag that keyed exclusion per (mass, charge); as a side effect it was also the "
          "only thing that made precursor_charges: \"separate\" fan out, so acquisition geometry "
          "depended on an exclusion flag. Acquiring several charge states of one species is now "
          "requested directly:\n"
          "  precursor_selection.precursor_charges: \"separate\"     one MS2 per charge state\n"
          "  precursor_selection.precursor_charges: \"multiplexed\"  one MS2 co-isolating them\n"
          "Exclusion is mass-keyed. There is no replacement for re-selecting one mass at a "
          "different charge on a LATER survey.");

    rejectUnknownKeys(ps,
        {"rt_window", "targeting", "consider_all_charges",
         "precursor_charges", "strict_inclusion", "tie_threshold",
         "rank_by", "max_precursors", "min_precursor_charge", "additional_scans", "exploration",
         "tag_expansion"},
        "precursor_selection");
    targeting_.rt_window = ps.value("rt_window", 180.0);
    targeting_.consider_all_charges = ps.value("consider_all_charges", false);
    targeting_.precursor_charges = parseChargeMode(ps, "precursor_charges", "precursor_selection");
    targeting_.strict_inclusion = ps.value("strict_inclusion", false);
    targeting_.tie_threshold = ps.value("tie_threshold", 0.1);

    // targeting: int -> string enum. Values map to the CODE, not the doc comments -- 2 is in-depth
    // and 3 is exclusion (PrecursorSelection.cpp:136-139 logs exactly that), while MethodConfig.cs,
    // Config.h and PrecursorSelection.cpp:544-546 all had 2 and 3 the wrong way round.
    // Rejected rather than defaulted: an unknown value used to fall through silently.
    {
      const std::string t = ps.value("targeting", std::string("none"));
      if (t == "none") targeting_.mode = 0;
      else if (t == "inclusion") targeting_.mode = 1;
      else if (t == "in_depth") targeting_.mode = 2;
      else if (t == "exclusion_masses") targeting_.mode = 3;
      else
        throw std::invalid_argument(
            "Config: precursor_selection.targeting must be one of \"none\", \"inclusion\", "
            "\"in_depth\", \"exclusion_masses\"; got \"" + t + "\" (values are case-sensitive).");
    }

    if (targeting_.mode == 1)
      std::cout << "Inclusion mode: " << (targeting_.strict_inclusion ? "strict" : "non-strict") << "\n";

    // tag_expansion: these two used to sit in `flashtnt`, which was a misnomer -- neither is a
    // FLASHTagger/FLASHExtender Param. max_ptm_count is read only by
    // PrecursorSelection::generatePTMCombinations_, and max_flanking_mass_diff is a call argument
    // FLASHIda passes to the static FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs.
    // Both belong to FLASHIda's tag-based target expansion, so they are authored where that feature
    // lives. STORAGE IS DELIBERATELY UNCHANGED (still TargetingConfig): this is a parse-path move
    // only, so every read site keeps working untouched. The old flashtnt placement now fails through
    // the ordinary unknown-key path -- no migration message, by design.
    auto te = ps.value("tag_expansion", json::object());
    rejectUnknownKeys(te, {"max_ptm_count", "max_flanking_mass_diff"}, "precursor_selection.tag_expansion");
    targeting_.max_total_ptm_count = te.value("max_ptm_count", 3);
    targeting_.max_flanking_mass_diff = te.value("max_flanking_mass_diff", 50000.0);

    // --- flashtnt section (FLASHTagger/FLASHExtender tuning) ---
    auto flashtnt = config.value("flashtnt", json::object());
    rejectUnknownKeys(flashtnt,
        {"min_length", "max_length", "allow_gap", "max_aa_in_gap", "max_blind_mod_count", "max_mod_mass", "fixed_mod"},
        "flashtnt");
    targeting_.min_tag_length = flashtnt.value("min_length", 3);
    targeting_.max_tag_length = flashtnt.value("max_length", 8);
    targeting_.allow_gap = flashtnt.value("allow_gap", false);
    targeting_.max_aa_in_gap = flashtnt.value("max_aa_in_gap", 2);
    targeting_.max_blind_mod_count = flashtnt.value("max_blind_mod_count", 2);
    targeting_.max_mod_mass = flashtnt.value("max_mod_mass", 700.0);  // 700 preserves current MS2 behavior; NOT the extender's own 500 default
    targeting_.fixed_mod.clear();
    if (flashtnt.contains("fixed_mod") && flashtnt["fixed_mod"].is_array())
    {
      for (const auto& m : flashtnt["fixed_mod"])
        targeting_.fixed_mod.push_back(m.get<std::string>());
    }

    // --- tagging section (acquisition-workflow controls; algorithm knobs live in flashtnt) ---
    auto tagging = config.value("tagging", json::object());
    rejectUnknownKeys(tagging, {"follow_up_scan"}, "tagging");

    if (tagging.contains("follow_up_scan") && tagging["follow_up_scan"].is_object())
    {
      auto fus = tagging["follow_up_scan"];
      rejectUnknownKeys(fus, kScanKeys, "tagging.follow_up_scan");
      parseScanConfig(fus, targeting_.tagging_follow_up_scan, "Orbitrap");
    }

    // --- files section (paths only; loading stays in FLASHIda) ---
    auto files = config.value("files", json::object());
    rejectUnknownKeys(files, {"target_logs", "fasta", "inclusion_list", "ptm_list"}, "files");

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

    // --- characterization section: WHETHER and HOW we characterize ---
    // `mode` is the single MS3 switch. It absorbs the old `objective` key AND the two
    // selection gates (selection_strategy.ms2.selection / .ms3.selection), which were booleans in
    // disguise. Decisions only -- the MS3 scan's instrument parameters stay in ms_settings.ms3.
    auto charact = config.value("characterization", json::object());
    // Migration, checked BEFORE the allowlist so the message names the replacement instead of the
    // generic unknown-key error. ms3_all_charges was a bool; its two states are now values of a
    // three-valued mode whose third value (multiplexed) it could not express.
    if (charact.contains("ms3_all_charges"))
      throw std::invalid_argument(
          "Config: characterization.ms3_all_charges was replaced by characterization.fragment_charges "
          "(ADR-0016). Use \"single\" for the old false (the fragment's best-MS2 charge) or "
          "\"separate\" for the old true (one MS3 per observed charge); \"multiplexed\" co-isolates "
          "them into one MS3.");
    rejectUnknownKeys(charact,
        {"mode", "protein_sequence", "max_targets", "min_fragment_charge", "min_target_mass",
         "fragment_charges", "exploration"},
        "characterization");
    {
      // Hard-rejected, not defaulted. The old parse was
      //     if (obj_str == "coverage") Coverage; else Ambiguity;
      // so "Coverage", "covrage" or any typo silently selected ambiguity, and with mode now
      // carrying the on/off bit a typo'd "Off" would silently ENABLE MS3.
      const std::string m = charact.value("mode", std::string("off"));
      if (m == "off")
        characterization_.mode = CharacterizationMode::Off;
      else if (m == "ambiguity")
      {
        characterization_.mode = CharacterizationMode::Ambiguity;
        characterization_.objective = CharacterizationObjective::Ambiguity;
      }
      else if (m == "coverage")
      {
        characterization_.mode = CharacterizationMode::Coverage;
        characterization_.objective = CharacterizationObjective::Coverage;
      }
      else if (m == "exhaustive")
      {
        characterization_.mode = CharacterizationMode::Exhaustive;
        // The engine branches on `objective` and never reads `mode` -- it has no read site outside
        // this file. Assigning only the mode here would ship a mode byte-identical to "ambiguity":
        // accepted, green, and inert (ADR-0023 D-a).
        characterization_.objective = CharacterizationObjective::Exhaustive;
      }
      else
        throw std::invalid_argument(
            "Config: characterization.mode must be one of \"off\", \"ambiguity\", \"coverage\", "
            "\"exhaustive\"; "
            "got \"" + m + "\" (values are case-sensitive).");

      characterization_.protein_sequence = charact.value("protein_sequence", "");
      characterization_.max_targets = charact.value("max_targets", 3);
      characterization_.min_fragment_charge = charact.value("min_fragment_charge", 0);
      // Read beside its sibling floors. Inert unless mode == Exhaustive; no other objective
      // consults it, because no other objective targets a raw deconvolved mass.
      characterization_.min_target_mass = charact.value("min_target_mass", 0.0);
      characterization_.fragment_charges = parseChargeMode(charact, "fragment_charges", "characterization");
    }

    // --- conditional_ms2 (top-level only) ---
    targeting_.conditional_ms2_enabled = config.value("conditional_ms2", false);

    // --- quantification (ADR-0038) ---
    auto quant = config.value("quantification", json::object());

    // Two retired keys, both worth their own message: silently rejecting them as "unknown" would
    // send an author looking for a typo when the model underneath changed.
    if (quant.contains("follow_up_scan"))
      throw std::invalid_argument(
          "Config: 'quantification.follow_up_scan' is retired (ADR-0038). The scan it named was "
          "acquired but never measured, while the scan that WAS measured -- the base MS2 -- could "
          "not release the reporter ion. The roles are now explicit: 'ms_settings.ms2_quant' is the "
          "quantification scan (rostered per precursor, measured), and 'ms_settings.ms2' is the "
          "identification scan a differential verdict buys.");
    if (quant.contains("only_one_condition"))
      throw std::invalid_argument(
          "Config: 'quantification.only_one_condition' is retired (ADR-0038). It was never "
          "reachable -- no emit DTO carried it -- and its intent is now unconditional: a condition "
          "that is WHOLLY absent reports 'differential', because a species present in one condition "
          "and absent in the other is the strongest result the experiment can produce.");

    rejectUnknownKeys(quant,
        {"enabled", "labelling", "reporter_mz_tol", "fold_change_threshold",
         "conditions", "correction_matrix",
         "identify", "enriched_in"},  // ADR-0039
        "quantification");
    quant_.enabled = quant.value("enabled", false);
    quant_.labelling = quant.value("labelling", std::string("tmt6plex"));
    quant_.reporter_mz_tol = quant.value("reporter_mz_tol", 0.002);
    quant_.fold_change_threshold = quant.value("fold_change_threshold", 1.4);

    // ADR-0039. Which verdicts buy the identification scan. Hard-rejected rather than defaulted,
    // for exactly characterization.mode's reason: a typo'd "Differential" must not silently select
    // a different acquisition policy, and here the wrong branch changes what the run acquires.
    {
      // NOT quant.value("identify", ...). ToCppJson uses the stock JavaScriptSerializer, which
      // EMITS NULLS -- value() on a present-but-null string throws a type_error, which is exactly
      // how `conditions` broke all 41 committed configs once (see the explicit-nulls block below
      // and its regression test). Absent and null must both mean "unauthored".
      const std::string id = (quant.contains("identify") && !quant["identify"].is_null())
                                 ? quant["identify"].get<std::string>()
                                 : "differential";
      if (id == "differential")     quant_.identify = QuantIdentify::Differential;
      else if (id == "quantified")  quant_.identify = QuantIdentify::Quantified;
      else if (id == "all")         quant_.identify = QuantIdentify::All;
      else if (id == "none")        quant_.identify = QuantIdentify::None;
      else
        throw std::invalid_argument(
            "Config: quantification.identify must be one of \"differential\", \"quantified\", "
            "\"all\", \"none\"; got \"" + id + "\" (values are case-sensitive). It selects which "
            "quantification verdicts buy the identification scan ms_settings.ms2 (ADR-0039).");
    }

    const auto& valid_labellings = Quantification::labellingNames();
    if (std::find(valid_labellings.begin(), valid_labellings.end(), quant_.labelling)
        == valid_labellings.end())
    {
      std::string known;
      for (const auto& n : valid_labellings) { if (!known.empty()) known += ", "; known += n; }
      throw std::invalid_argument(
          "Config: unknown 'quantification.labelling' value '" + quant_.labelling
          + "'. Valid: " + known + ". ('none' is deliberately not accepted -- "
            "'quantification.enabled' is the switch.)");
    }

    // NOTE the is_null() arm on both of these, and it is not defensive padding. ToCppJson uses the
    // STOCK JavaScriptSerializer, which EMITS nulls -- so every config that authors neither key
    // sends `"conditions": null, "correction_matrix": null`, and a bare contains() check would
    // throw "must be an ARRAY" for all 41 of them. Same idiom resolveFollowUp_ already uses above.
    if (quant.contains("correction_matrix") && !quant["correction_matrix"].is_null())
    {
      if (!quant["correction_matrix"].is_array())
        throw std::invalid_argument(
            "Config: 'quantification.correction_matrix' must be an array of strings, one per "
            "channel (Thermo data-sheet spelling, e.g. \"0.0/0.0/5.09/0.0\"). Omit it to use the "
            "method's stock matrix; an all-zero matrix turns correction off.");
      for (const auto& row : quant["correction_matrix"])
        quant_.correction_matrix.push_back(row.get<std::string>());
    }

    // Conditions are an ARRAY, never an object: nlohmann's object_t is a std::map, so an object
    // would sort the two alphabetically and silently decide which is the numerator. The array order
    // IS the ratio direction. Channel NAMES resolve to ordinals here, at load, so an unknown
    // channel fails loudly instead of reading the wrong intensity at acquisition time.
    if (quant.contains("conditions") && !quant["conditions"].is_null())
    {
      if (!quant["conditions"].is_array())
        throw std::invalid_argument(
            "Config: 'quantification.conditions' must be an ARRAY of exactly two objects "
            "{name, channels}. The array order is the ratio direction: "
            "fold_change = mean(conditions[0]) / mean(conditions[1]).");

      const std::vector<std::string> channel_names = Quantification::channelNames(quant_.labelling);
      std::string known_channels;
      for (const auto& c : channel_names)
      { if (!known_channels.empty()) known_channels += ", "; known_channels += c; }

      for (const auto& cond : quant["conditions"])
      {
        rejectUnknownKeys(cond, {"name", "channels"}, "quantification.conditions[]");
        QuantConfig::Condition parsed;
        parsed.name = cond.value("name", std::string{});
        // ADR-0039. "either" is quantification.enriched_in's sentinel for "either direction", so a
        // condition may not claim the name -- otherwise `enriched_in: "either"` is ambiguous and
        // whichever reading loses is silently wrong about which way the experiment ran.
        if (parsed.name == "either")
          throw std::invalid_argument(
              "Config: a quantification condition may not be named \"either\". That is the "
              "quantification.enriched_in sentinel meaning 'either direction' (ADR-0039); a "
              "condition of that name would make enriched_in ambiguous. Rename the condition.");
        if (!cond.contains("channels") || !cond["channels"].is_array() || cond["channels"].empty())
          throw std::invalid_argument(
              "Config: quantification condition '" + parsed.name
              + "' must name a non-empty 'channels' array. Valid channels for "
              + quant_.labelling + ": " + known_channels + ".");
        for (const auto& ch : cond["channels"])
        {
          const std::string cname = ch.get<std::string>();
          auto it = std::find(channel_names.begin(), channel_names.end(), cname);
          if (it == channel_names.end())
            throw std::invalid_argument(
                "Config: quantification condition '" + parsed.name + "' names unknown channel '"
                + cname + "'. Valid channels for " + quant_.labelling + ": " + known_channels + ".");
          parsed.channels.push_back(static_cast<size_t>(it - channel_names.begin()));
        }
        quant_.conditions.push_back(std::move(parsed));
      }
    }

    // ADR-0039. Direction, named by CONDITION rather than as up/down, and resolved to an ordinal
    // here -- exactly the treatment channel names get above, and for the same reason: an unknown
    // name fails loudly at load instead of silently reading the wrong side of the ratio.
    //
    // OUTSIDE the conditions block deliberately, so `enriched_in` authored without conditions is
    // rejected rather than ignored. Null-tolerant for the JavaScriptSerializer reason above.
    {
      const std::string ei = (quant.contains("enriched_in") && !quant["enriched_in"].is_null())
                                 ? quant["enriched_in"].get<std::string>()
                                 : "either";
      if (ei != "either")
      {
        for (size_t i = 0; i < quant_.conditions.size(); ++i)
          if (quant_.conditions[i].name == ei) { quant_.enriched_in = static_cast<int>(i); break; }

        if (quant_.enriched_in < 0)
        {
          std::string known;
          for (const auto& c : quant_.conditions)
          { if (!known.empty()) known += ", "; known += "\"" + c.name + "\""; }
          if (known.empty()) known = "(none authored)";
          throw std::invalid_argument(
              "Config: quantification.enriched_in names unknown condition \"" + ei
              + "\". Valid: " + known + ", or \"either\" for either direction. It names the "
                "condition a species must be ENRICHED IN for a differential verdict to buy the "
                "identification scan (ADR-0039).");
        }
      }
    }

    // --- faims ---
    auto faims_section = config.value("faims", json::object());
    rejectUnknownKeys(faims_section, {"cv_values", "max_cv_skip", "cv_precursor_threshold"}, "faims");
    if (faims_section.contains("cv_values") && faims_section["cv_values"].is_array())
    {
      for (const auto& v : faims_section["cv_values"])
        faims_.cv_values.push_back(v.get<double>());
    }
    faims_.max_cv_skip = faims_section.value("max_cv_skip", 0);
    faims_.precursor_threshold = faims_section.value("cv_precursor_threshold", 15);
    // An empty cv_values is the ONLY way to say "no FAIMS" (ADR-0012). This used to be
    // `size() > 1`, which conflated "FAIMS is off" with "there is nothing to cycle between": a
    // single-CV method -- an ordinary fixed-CV FAIMS run -- silently acquired at whatever CV the
    // instrument method carried, and FLASHIda had no way to say so. Cycling is a separate question,
    // answered by FAIMS::isCycling().
    faims_.enabled = !faims_.cv_values.empty();

    // --- ms_settings: a library of scan configs. NOTHING here decides whether a scan happens. ---
    // ms1/ms2/ms3 are bare objects (they were arrays); extra MS2 blocks live in the name-keyed
    // additional_ms2 map and reach the dispatch roster only by being REFERENCED.
    auto ms_settings = config.value("ms_settings", json::object());
    // ms2_quant is a bare slot beside ms1/ms2/ms3, not a name reference (ADR-0038) -- the same
    // shape characterization already uses, where the decision section holds no scan config and
    // ms_settings.ms3 supplies the parameters with no key pointing at it.
    rejectUnknownKeys(ms_settings, {"ms1", "ms2", "ms2_quant", "ms3", "additional_ms2"}, "ms_settings");

    if (ms_settings.contains("ms2") && ms_settings["ms2"].is_array())
      throw std::invalid_argument(
          "Config: 'ms_settings.ms2' is no longer an array. It is now a single scan object; "
          "additional MS2 configs go in 'ms_settings.additional_ms2' as a name->object map and are "
          "referenced from precursor_selection.additional_scans or tagging.follow_up_scan.");
    if (ms_settings.contains("ms3") && ms_settings["ms3"].is_array())
      throw std::invalid_argument(
          "Config: 'ms_settings.ms3' is no longer an array. It is now a single scan object. "
          "There is no additional_ms3: every level-3 consumer reads scans[0], so a second MS3 "
          "config was never reachable.");

    // levels_ is materialised unconditionally for {1,2,3}. This used to hold only by accident --
    // selection_strategy was required and every config named all three levels. toleranceList() and
    // explorationToleranceList() walk levels_ POSITIONALLY to build the DoubleLists that construct
    // Deconvolution, so a level going missing would silently shift every tols_[ms_level-1] index.
    for (int lvl : {1, 2, 3})
      if (levels_.find(lvl) == levels_.end())
        levels_[lvl] = MSLevelConfig{};

    auto parseNamedScan = [&](const json& node, const std::string& path) {
      rejectUnknownKeys(node, kScanKeys, path);
      ScanConfig sc;
      parseScanConfig(node, sc, "");
      return sc;
    };

    auto ms1_json = ms_settings.value("ms1", json::object());
    rejectUnknownKeys(ms1_json, kScanKeys, "ms_settings.ms1");
    ScanConfig ms1_scan;
    parseScanConfig(ms1_json, ms1_scan, "");
    levels_[1].scans.push_back(ms1_scan);

    // Definitions first; the roster is assembled afterwards from the reference list.
    ScanConfig primary_ms2;
    bool has_primary_ms2 = false;
    if (ms_settings.contains("ms2") && ms_settings["ms2"].is_object())
    {
      primary_ms2 = parseNamedScan(ms_settings["ms2"], "ms_settings.ms2");
      has_primary_ms2 = true;
    }
    ScanConfig quant_scan;
    bool has_quant_scan = false;
    if (ms_settings.contains("ms2_quant") && ms_settings["ms2_quant"].is_object())
    {
      quant_scan = parseNamedScan(ms_settings["ms2_quant"], "ms_settings.ms2_quant");
      has_quant_scan = true;
    }
    if (ms_settings.contains("ms3") && ms_settings["ms3"].is_object())
      levels_[3].scans.push_back(parseNamedScan(ms_settings["ms3"], "ms_settings.ms3"));

    // additional_ms2: user-authored KEYS, so they cannot be allowlisted -- validated as identifiers
    // instead. The VALUES go through the normal 17-key scan allowlist.
    std::map<std::string, ScanConfig> additional_ms2;
    // `&& !is_null()` is load-bearing, not defensive noise: JavaScriptSerializer does not omit null
    // properties, so ToCppJson emits `"additional_ms2": null` for the ~30 configs that have no
    // extras. contains() is true for an explicit null while is_object() is false, so checking only
    // contains() rejected almost every real config. Same idiom the exploration parse already uses.
    if (ms_settings.contains("additional_ms2") && !ms_settings["additional_ms2"].is_null())
    {
      const auto& add = ms_settings["additional_ms2"];
      if (!add.is_object())
        throw std::invalid_argument("Config: ms_settings.additional_ms2 must be a name->object map.");
      static const std::set<std::string> kReservedScanNames = {"ms1", "ms2", "ms3", "none", "off", "all"};
      for (auto it = add.begin(); it != add.end(); ++it)
      {
        const std::string& name = it.key();
        if (!isValidScanName(name))
          throw std::invalid_argument(
              "Config: scan-config name '" + name + "' in ms_settings.additional_ms2 is invalid. "
              "Names are snake_case identifiers matching ^[a-z][a-z0-9_]{0,31}$.");
        if (kReservedScanNames.count(name) != 0)
          throw std::invalid_argument(
              "Config: '" + name + "' is a reserved word and cannot be used as a scan-config name. "
              "Reserved: ms1, ms2, ms3, none, off, all.");
        additional_ms2[name] = parseNamedScan(it.value(), "ms_settings.additional_ms2." + name);
      }
    }

    // The DISPATCH ROSTER: ms_settings.ms2 first, then precursor_selection.additional_scans in the
    // order the array gives. Nothing else fires unconditionally -- a block that only backs a
    // follow-up reference is simply absent here. Order comes from the ARRAY, never from map
    // iteration: nlohmann's object_t is a std::map, so iterating additional_ms2 would sort the
    // names alphabetically and silently reorder dispatch.
    // ADR-0038 inverts which of the two MS2 configs is unconditional. With quantification enabled
    // the QUANTIFICATION scan takes the roster's primary slot -- it is the screen, so it must be
    // acquired to decide anything -- and ms_settings.ms2 becomes the IDENTIFICATION scan a
    // differential verdict buys, held on quant_ rather than dispatched. `ms_settings.ms2` therefore
    // means "the identification MS2" in every mode; only WHEN it fires changes.
    // Copy both scan configs unconditionally; only the ROSTER below is conditional.
    quant_.has_quant_scan = has_quant_scan;
    // ADR-0039. Carried onto quant_ because validate() is const and separate from this constructor,
    // and identification_scan cannot answer "was one authored" -- it is default-constructed when
    // ms_settings.ms2 is absent, which is precisely the state validate() has to reject.
    quant_.has_identification_scan = has_primary_ms2;
    if (has_quant_scan)  quant_.quantification_scan = quant_scan;
    if (has_primary_ms2) quant_.identification_scan = primary_ms2;

    if (quant_.enabled && has_quant_scan)
    {
      levels_[2].scans.push_back(quant_scan);
    }
    else if (has_primary_ms2)
    {
      levels_[2].scans.push_back(primary_ms2);
    }
    {
      std::set<std::string> seen;
      const auto& add_scans = ps.value("additional_scans", json::array());
      if (!add_scans.is_array())
        throw std::invalid_argument("Config: precursor_selection.additional_scans must be an array of names.");
      for (const auto& n : add_scans)
      {
        const std::string name = n.get<std::string>();
        if (!seen.insert(name).second)
          throw std::invalid_argument(
              "Config: precursor_selection.additional_scans lists '" + name + "' more than once.");
        auto found = additional_ms2.find(name);
        if (found == additional_ms2.end())
          throw std::invalid_argument(
              "Config: precursor_selection.additional_scans references unknown MS2 scan config '"
              + name + "'. " + knownScanNames(additional_ms2));
        levels_[2].scans.push_back(found->second);
      }
      additional_scan_names_ = seen;
    }

    // Follow-ups are name references now, resolved here so no downstream consumer learns that
    // names exist -- FLASHIda.cpp keeps reading targeting_.tagging_follow_up_scan verbatim.
    // Tagging is the only remaining name-referenced follow-up. Quantification's went away with
    // ADR-0038: its bought scan is ms_settings.ms2, a bare slot, so there is no name to resolve.
    resolveFollowUp_(config, "tagging", additional_ms2, targeting_.tagging_follow_up_scan,
                     targeting_.has_tagging_follow_up, targeting_.tagging_follow_up_name);

    // additional_ms2 is ONE flat namespace serving two roles, and the roles are mutually exclusive:
    // additional_scans dispatches unconditionally once per precursor, a follow_up_scan fires
    // conditionally off a returning MS2. A name in both does BOTH -- the same scan config acquired
    // twice per precursor, at two different priorities (2 for the roster entry, 0 for the
    // follow-up), with no diagnostic. Nothing else in the schema catches it: both references
    // resolve, so the dangling-name check passes, and the block is referenced, so the unreferenced
    // warning below stays quiet.
    auto rejectDoubleDuty = [this](const std::string& section, const std::string& name) {
      if (name.empty() || additional_scan_names_.count(name) == 0) return;
      throw std::invalid_argument(
          "Config: ms_settings.additional_ms2." + name + " is listed in "
          "precursor_selection.additional_scans AND referenced by " + section
          + ".follow_up_scan. It would fire unconditionally once per precursor and again as a "
            "conditional follow-up. Define two entries, or drop it from additional_scans.");
    };
    rejectDoubleDuty("tagging", targeting_.tagging_follow_up_name);

    // A definition nobody references never fires. That is legal but almost always a mistake, and it
    // is the only check that catches a typo on the DEFINITION side (a typo on the reference side is
    // caught above). Warn rather than throw -- commenting a scan out of the roster while tuning is
    // a normal thing to do.
    // Quantification no longer appears here: its two scans are the bare ms_settings.ms2_quant and
    // ms_settings.ms2 slots, so neither is an additional_ms2 definition that could go unreferenced.
    for (const auto& kv : additional_ms2)
      if (additional_scan_names_.count(kv.first) == 0
          && targeting_.tagging_follow_up_name != kv.first)
        std::cout << "[CONFIG-WARN] ms_settings.additional_ms2." << kv.first
                  << " is defined but never referenced; it will never be acquired.\n";

    // --- scheduling ---
    auto sched = config.value("scheduling", json::object());
    // target_depth is ACCEPTED AND UNUSED, deliberately. It sizes the instrument's custom-scan
    // queue, which only the host can do -- Flash.cs:ProcessSpectrum is the sole submitter, and the
    // engine has no visibility into how many commands the instrument is holding (docs/adr/0033).
    // It is listed here because the schema is strict on both sides (ADR-0007) and
    // test-data/config_schema_reference.json is generated from the C# model, so a host-only key
    // that C++ rejected would fail ConfigSchemaParity_test rather than any real config error.
    // Do not "wire it up" engine-side; there is nothing here for it to control.
    rejectUnknownKeys(sched, {"cycle_time", "scan_timeout", "agc_interval_seconds", "target_depth"}, "scheduling");
    auto ct = sched.value("cycle_time", json::object());
    rejectUnknownKeys(ct, {"enabled", "value_ms"}, "scheduling.cycle_time");
    scheduling_.cycle_time_enabled = ct.value("enabled", false);
    scheduling_.cycle_time_ms = ct.value("value_ms", 60000.0);
    auto to = sched.value("scan_timeout", json::object());
    rejectUnknownKeys(to, {"enabled", "value_ms"}, "scheduling.scan_timeout");
    scheduling_.timeout_enabled = to.value("enabled", false);
    scheduling_.timeout_ms = to.value("value_ms", 30000.0);
    double agc_interval_sec = sched.value("agc_interval_seconds", 1.0);
    scheduling_.agc_interval_ms = static_cast<uint64_t>(agc_interval_sec * 1000.0);

    // --- exploration: one block per decision section, each sweeping the scans that section
    //     dispatches. precursor_selection.exploration -> level 2, characterization.exploration ->
    //     level 3. They must stay separate: a single config legitimately sweeps different ranges at
    //     the two levels (method_exploration_ms3_followup sweeps HCD 20-40 at MS2, CID 15-35 at MS3).
    //     There is no level-1 exploration: it was parsed nowhere and discarded on both sides.
    auto parseExploration = [&](const json& parent, const std::string& path, int level) {
      if (!parent.contains("exploration") || parent["exploration"].is_null())
        return;
      const auto& e = parent["exploration"];
      if (!e.is_object())
        throw std::invalid_argument("Config: " + path + ".exploration must be an object.");
      rejectUnknownKeys(e,
          {"metric", "ce_min", "ce_max", "ce_step", "overrides", "remaining_precursor_target",
           "reaction_time_min", "reaction_time_max", "reaction_time_step", "activations",
           "tolerance_ppm"},
          path + ".exploration");

      MSLevelConfig& cfg = levels_[level];

      // Hard-rejected. The old parse fell through to None, so a typo silently collapsed an
      // N-variant sweep to a single scan -- the exact opposite direction from the selection typo,
      // which silently ENABLED a level. Neither should guess.
      const std::string m = e.value("metric", std::string("none"));
      if (m == "none") cfg.exploration = ExplorationMetric::None;
      else if (m == "mass_count") cfg.exploration = ExplorationMetric::MassCount;
      else if (m == "remaining_precursor") cfg.exploration = ExplorationMetric::RemainingPrecursor;
      else if (m == "fragment_count") cfg.exploration = ExplorationMetric::FragmentCount;
      else
        throw std::invalid_argument(
            "Config: " + path + ".exploration.metric must be one of \"none\", \"mass_count\", "
            "\"remaining_precursor\", \"fragment_count\"; got \"" + m + "\".");

      cfg.ce_min = e.value("ce_min", 20.0);
      cfg.ce_max = e.value("ce_max", 40.0);
      cfg.ce_step = e.value("ce_step", 5.0);
      cfg.remaining_precursor_target = e.value("remaining_precursor_target", 0.1);
      cfg.rt_min = e.value("reaction_time_min", 0.0);
      cfg.rt_max = e.value("reaction_time_max", 0.0);
      cfg.rt_step = e.value("reaction_time_step", 1.0);
      // First-class now. It used to be smuggled through the overrides map and then ERASED from it
      // at this point -- before Exploration.cpp:749 tested that same map for emptiness to decide
      // whether to acquire the production scan, so a tolerance-only map silently suppressed a scan.
      // 0 means "use deconvolution.tol for this level"; resolved just below.
      cfg.exploration_tolerance_ppm = e.value("tolerance_ppm", 0.0);

      if (e.contains("activations") && e["activations"].is_array())
        for (const auto& a : e["activations"])
          cfg.activations.push_back(a.get<std::string>());

      if (e.contains("overrides") && e["overrides"].is_object())
      {
        const auto& ov = e["overrides"];
        for (auto ov_it = ov.begin(); ov_it != ov.end(); ++ov_it)
        {
          if (ov_it.key() == "tolerance_ppm")
            throw std::invalid_argument(
                "Config: " + path + ".exploration.overrides no longer carries 'tolerance_ppm'. "
                "It is a first-class key: move it to " + path + ".exploration.tolerance_ppm. "
                "(applyOverrides has no branch for it, so leaving it here would drop it silently.)");
          cfg.overrides[ov_it.key()] = ov_it.value().get<std::string>();
        }
      }
    };
    parseExploration(ps, "precursor_selection", 2);
    parseExploration(charact, "characterization", 3);

    // --- THE PROJECTION: derive every level's selection state from characterization.mode ---
    // Runs after all sections are parsed and BEFORE the tol loop below, which sizes itself from
    // levels_. This is what makes `mode` the single MS3 switch.
    applyCharacterizationMode_(ps.value("rank_by", std::string("qscore")),
                               ps.value("max_precursors", 10),
                               ps.value("min_precursor_charge", 0));

    // Validate tol array length covers all configured MS levels. levels_ is now always {1,2,3}, so
    // this is effectively a fixed ">= 3" -- kept level-derived so it stays correct if that changes.
    int max_level = 0;
    for (const auto& [lvl, unused_cfg] : levels_)
      max_level = std::max(max_level, lvl);
    if (static_cast<int>(tol_values.size()) < max_level)
      throw std::invalid_argument("deconvolution.tol must have at least "
        + std::to_string(max_level) + " entries (one per MS level), got "
        + std::to_string(tol_values.size()));

    // Per-level tolerances. An unset exploration tolerance falls back to the level's base tol --
    // the same value the old overrides path produced, so this is value-preserving.
    for (auto& [lvl, cfg] : levels_)
    {
      cfg.tolerance_ppm = tol_values[lvl - 1];
      if (cfg.exploration_tolerance_ppm <= 0.0)
        cfg.exploration_tolerance_ppm = tol_values[lvl - 1];
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

    // --- runtime section (log folder, optional) ---
    // The json::object() default is what makes an absent section and "runtime": {} both mean
    // "open nothing" -- do not replace it with a required section.
    auto rt_section = config.value("runtime", json::object());
    // Dedicated migration error: the generic unknown-key message would name the offending key
    // without saying that ALL FIVE collapsed into one, or that the value is now a folder rather
    // than a file path.
    for (const char* dead : {"ida_log_path", "scan_commands_path", "scan_results_path",
                             "identification_log_path", "pooled_identification_log_path"})
    {
      if (rt_section.contains(dead))
        throw std::invalid_argument(
            std::string("Config: runtime.") + dead + " has been removed. The five per-stream log "
            "paths are replaced by a single 'runtime.log_dir' naming the FOLDER that receives all "
            "of them, under their fixed basenames (ida.log, scan_commands.tsv, scan_results.tsv, "
            "identification.tsv, pooled_identification.tsv).");
    }
    rejectUnknownKeys(rt_section, {"log_dir"}, "runtime");
    runtime_.log_dir = rt_section.value("log_dir", std::string{});

    // SNR threshold (hardcoded in original parseJSONConfig_)
    targeting_.snr_threshold = 1.0;

    validate();
  }

  void Config::applyCharacterizationMode_(const std::string& rank_by, int max_precursors,
                                          int min_precursor_charge)
  {
    // ---- LEVEL 1 ----
    // Do not remove. MSLevelConfig::selection defaults to None, so if level 1 is left unassigned
    // FLASHIda.cpp:169-170 short-circuits before MS1 selection and the run acquires NOTHING at all,
    // silently. This is the only level whose selection VALUE is read (PrecursorSelection.cpp:246
    // picks sortByIntensity vs sortByQscore); at levels 2 and 3 only None-vs-non-None matters.
    if (rank_by == "qscore") levels_[1].selection = SelectionMetric::QScore;
    else if (rank_by == "intensity") levels_[1].selection = SelectionMetric::Intensity;
    else if (rank_by == "none") levels_[1].selection = SelectionMetric::None;
    else
      throw std::invalid_argument(
          "Config: precursor_selection.rank_by must be one of \"qscore\", \"intensity\", \"none\"; "
          "got \"" + rank_by + "\" (values are case-sensitive).");
    levels_[1].max_targets = max_precursors;   // the MS2 count
    levels_[1].min_charge = min_precursor_charge;

    // ---- LEVELS 2 AND 3 ----
    // MS3 requires BOTH gates non-None: FLASHIda.cpp:368 and Exploration.cpp:903 test level 2,
    // Exploration.cpp:905 tests level 3. Driving both from one enum is what makes the incoherent
    // states (MS3 on with MS2 off) unrepresentable rather than merely discouraged.
    //
    // Intensity is the value used when on: intensity and qscore share a case in
    // ProteoformTracker::selectNextLevelTargets, and every MS3-enabled config in the corpus used
    // intensity, so this reproduces today's matcher exactly.
    // Keep this an INEQUALITY rather than an enumeration of the on-values: it is what makes every
    // future mode project levels 2 and 3 correctly without touching this function. Level 1 is
    // assigned from rank_by above, unconditionally, for the reason stated there.
    const bool on = characterization_.mode != CharacterizationMode::Off;
    levels_[2].selection = on ? SelectionMetric::Intensity : SelectionMetric::None;
    levels_[3].selection = on ? SelectionMetric::Intensity : SelectionMetric::None;

    // The MS3 budget and fragment-charge floor are read off level 2
    // (ProteoformTracker.cpp:455, Exploration.cpp:908/:981) because level-2 selection is what
    // produces level-3 targets. They are AUTHORED in characterization, where they belong.
    levels_[2].max_targets = characterization_.max_targets;
    levels_[2].min_charge = characterization_.min_fragment_charge;
  }

  void Config::validate() const
  {
    // Re-keyed off the REFERENCE rather than the resolved block's activation string. The old
    // sentinel (`tagging_follow_up_scan.activation.empty()`) was a proxy for "nobody configured
    // one"; now that the block is defined elsewhere and merely named, presence of the name is the
    // direct answer.
    if (targeting_.conditional_ms2_enabled && !targeting_.has_tagging_follow_up)
      throw std::invalid_argument(
          "conditional_ms2 is true but tagging.follow_up_scan is not set. Name an "
          "ms_settings.additional_ms2 entry, or set conditional_ms2 to false.");

    // --- quantification (ADR-0038) --------------------------------------------------------------
    // Three structural rejections. These are NOT chemistry judgements -- there is deliberately no
    // guard that the quantification scan's activation can release a reporter ion, and none that
    // reporter_mz_tol is narrower than half the scheme's channel spacing. The user is trusted on
    // both, which is exactly why quant_verdict is a four-way enum: it is now the only route by
    // which a misconfigured screen is discovered.
    if (quant_.enabled)
    {
      // 1. No quantification scan => nothing is rostered, ms_settings.ms2 stays unconditional, and
      //    the run is plain DDA that silently quantifies nothing.
      if (!quant_.has_quant_scan)
        throw std::invalid_argument(
            "Config: quantification.enabled is true but ms_settings.ms2_quant is not set. "
            "ms2_quant is the quantification scan -- the one rostered per precursor and actually "
            "measured (ADR-0038); ms_settings.ms2 is the identification scan a differential verdict "
            "buys. Without ms2_quant nothing is measured and no identification scan is ever bought.");

      // 2. No ratio without exactly two conditions. There is deliberately nothing behind this: the
      //    positional 3-vs-3 split it replaces was correct only for six-plex and silently wrong for
      //    every other scheme this ADR enables.
      if (quant_.conditions.size() != 2)
        throw std::invalid_argument(
            "Config: quantification.conditions must name EXACTLY TWO conditions (got "
            + std::to_string(quant_.conditions.size())
            + "). fold_change = mean(conditions[0]) / mean(conditions[1]) is a two-group ratio; a "
              "time course or dose series needs a different statistic, not a ratio of the first "
              "two groups.");

      // 3. Exploration at level 2 dispatches CE-sweep variants INSTEAD of the roster, so the
      //    quantification scan is never acquired and every variant is labelled 'E' -- never
      //    measured, never buying anything. The feature dies completely, and the one guard that
      //    might have caught it (exactly-one-scan-config at a swept level) is SATISFIED here,
      //    because the inverted roster has exactly one entry. Incompatible by construction.
      if (hasExploration(2))
        throw std::invalid_argument(
            "Config: quantification.enabled and precursor_selection.exploration are incompatible. "
            "Exploration replaces the level-2 roster with CE-sweep variants, so the quantification "
            "scan would never be dispatched and nothing would ever be measured. Turn one off.");

      // 4. ADR-0039. ms_settings.ms2 is the scan the screen buys. Without it, identification_scan is
      //    DEFAULT-CONSTRUCTED and buildFollowUp builds the bought scan out of defaults. That was
      //    latent while only a Differential verdict reached it; `identify: "all"` fires it on EVERY
      //    precursor, so the gap gets a guard in the same change that widens it.
      //
      //    Required whenever enabled, INERT under identify: "none" rather than optional -- so all
      //    four identify values stay interchangeable and flipping the objective can never invalidate
      //    a config. That is ADR-0013's ms_settings.ms3 rule, and the same reason the enriched_in
      //    mismatch below warns instead of throwing.
      if (!quant_.has_identification_scan)
        throw std::invalid_argument(
            "Config: quantification.enabled is true but ms_settings.ms2 is not set. ms2 is the "
            "IDENTIFICATION scan the quantification screen buys (ADR-0038); ms2_quant is the scan "
            "that does the measuring. It is required whenever quantification is enabled -- inert "
            "under quantification.identify: \"none\", so the objective stays a one-word edit.");

      // ADR-0039. A WARNING, not a throw. enriched_in is live under `differential` and merely inert
      // under the other three -- that is ms_settings.ms3-under-mode-off, NOT only_one_condition,
      // which was unreachable in every possible config. Throwing would force a second edit every
      // time `identify` is flipped and would invalidate a template that sets enriched_in once, so
      // the inertness is announced instead -- the same [CONFIG-WARN] treatment an unreferenced
      // ms_settings.additional_ms2 block gets in the constructor above.
      if (quant_.enriched_in >= 0 && quant_.identify != QuantIdentify::Differential)
        std::cout << "[CONFIG-WARN] quantification.enriched_in is set but quantification.identify "
                     "is \""
                  << quantIdentifyName(quant_.identify)
                  << "\", so direction is not applied; it restricts only the \"differential\" "
                     "objective.\n";
    }

    // NO ACTIVATION-COUPLING THROW HERE, deliberately (ADR-0030). This used to reject an authored
    // scan config pairing ETD with reaction_time <= 0, or HCD/CID with collision_energy <= 0. Its
    // stated purpose was to stop a SILENT failure: ScanFactory dropped a zero reaction_time
    // entirely, so the instrument fell back to its own method default and nothing anywhere said so.
    //
    // ScanFactory now gates that key on the stage's ACTIVATION rather than on its value, so a zero
    // reaction time is commanded as a real 0 and the logged value equals the commanded one. The
    // silence the guard existed to prevent is gone, and with it the guard -- on both branches, so
    // ETD and HCD stay symmetric. A zero on either axis is now simply "do not fragment", which is
    // exactly what an exploration baseline asks for.

    // Restated on the ROSTER, not the definition map: a CE/RT sweep varies one base scan config, so
    // the level it sweeps must dispatch exactly one.
    for (const auto& [lvl, cfg] : levels_)
    {
      if (cfg.exploration != ExplorationMetric::None && cfg.scans.size() != 1)
        throw std::invalid_argument(
            std::string(lvl == 2 ? "precursor_selection" : "characterization")
            + ".exploration is enabled, so level " + std::to_string(lvl)
            + " must dispatch exactly one scan config; it dispatches "
            + std::to_string(cfg.scans.size())
            + ". Remove entries from precursor_selection.additional_scans.");
    }

    for (const auto& [lvl, cfg] : levels_)
    {
      if (cfg.exploration == ExplorationMetric::FragmentCount && characterization_.protein_sequence.empty())
        throw std::invalid_argument(
            "ExplorationMetric::FragmentCount at level " + std::to_string(lvl) +
            " requires a non-empty characterization.protein_sequence.");
    }

    // RemainingPrecursor scores a variant from the intensity inside its isolation window alone
    // (Exploration::precursorWindowIntensity_) and its winner is ALWAYS re-acquired by a production
    // scan built from the un-overridden config, so its pre-scans are throwaway measurements --
    // which is why ADR-0026 narrows their scan range to exactly that window. The narrowing leaves an
    // MS2 pre-scan with no fragments in it, so the winning variant's deconvolved spectrum is a
    // useless MS3 target list.
    //
    // This rejection does NOT protect the MS2->MS3 cascade -- it REPLACES its source. Exploration's
    // post-winner handling is ONE if/else chain: `!level_config.overrides.empty() ||
    // measuring_ms3_sweep` takes the production re-acquisition, and the cascade is that same chain's
    // `else if (group.msn_level < 3)` arm, so the two are mutually exclusive. Forcing every MS2
    // remaining_precursor config to carry non-empty overrides makes the first condition
    // unconditionally true at MS2, which makes the cascade arm structurally UNREACHABLE for them.
    // That is the intended outcome, not collateral: rather than cascading off a window-only
    // spectrum, MS3 is dispatched off the FULL-RANGE production re-acquisition -- it is not an
    // exploration variant, so it returns on the regular MS2 path and cascades from the stored MS2
    // spectrum like any other MS2. A better cascade source, one scan later.
    //
    // Without the rejection, empty overrides fail BOTH conditions -- measuring_ms3_sweep is level-3
    // only -- so the cascade arm does run, off the narrowed spectrum, and yields ZERO MS3 targets
    // with no throw and no warning: only `[MS3-PLAN] no_containing_fragment` and a user concluding
    // their protein did not fragment. Requiring non-empty overrides is the STATIC form of ADR-0020
    // gate #1: the production re-acquisition that makes narrowing safe is guaranteed by the schema
    // rather than by the author's habit of writing an analyzer override anyway.
    for (const auto& [lvl, cfg] : levels_)
    {
      // levels_ always holds {1,2,3} (Config.cpp:755), but an exploration block is parsed for levels 2
      // and 3 only (Config.cpp:686, :745-746), so level 1 is structurally None today. Skipped rather
      // than left to fall through, because `sect` below resolves everything that is not level 2 to
      // "characterization": a future level-1 sweep would be told to edit a key that does not govern
      // it. A guard, not dead code.
      if (lvl < 2) continue;

      // The level-3 arm applies only when MS3 actually runs. parseExploration returns early when a
      // section carries no "exploration" key, so levels_[3].exploration is written whenever that block
      // is PRESENT (Config.cpp:746) -- and nothing afterwards clears it, because
      // applyCharacterizationMode_ resets only `selection`. A leftover characterization.exploration
      // block therefore outlives characterization.mode == "off". Rejecting on it would make that
      // leftover a LOAD ERROR for a run that emits no MS3 and narrows no pre-scan, breaking ADR-0013's
      // promise that toggling MS3 off stays a one-word edit -- under "off" the MS3 keys are carried and
      // never read. Level 2 is deliberately not guarded: an MS2 sweep runs whatever
      // characterization.mode says.
      //
      // That promise is NOT upheld tree-wide, and the asymmetry is known rather than an oversight: the
      // FragmentCount / protein_sequence rejection a few loops above carries no such guard, so
      // `characterization: { mode: "off", exploration: { metric: "fragment_count" } }` with an empty
      // protein_sequence still throws at load. The guard is applied HERE because THIS rejection is new
      // with ADR-0026 and gets to be born correct; the FragmentCount check predates it, and widening an
      // existing rejection's guard is a behaviour change to configs that load today -- out of scope for
      // ADR-0026, not a claim that it is right.
      if (lvl == 3 && characterization_.mode == CharacterizationMode::Off) continue;

      const std::string sect = (lvl == 2 ? "precursor_selection" : "characterization");
      if (cfg.exploration == ExplorationMetric::RemainingPrecursor && cfg.overrides.empty())
        throw std::invalid_argument(
            sect + ".exploration.metric is \"remaining_precursor\" but " + sect
            + ".exploration.overrides is empty. A remaining_precursor sweep never keeps its "
              "pre-scans -- they are scanned over their isolation window only, and the winner is "
              "re-acquired -- so it must declare the settings they run at. Add an overrides block, "
              "e.g. \"overrides\": { \"analyzer\": \"IonTrap\" } (ADR-0026).");
    }

    // Level-matched multiplexing. [first_mass, last_mass] is ONE interval and a notch set is not,
    // so ADR-0026's binding cannot express a multiplexed readout: charge states 10-16 of a ~12 kDa
    // protein scatter their 2 Th windows across ~463 Th, and binding to the anchor alone would
    // isolate seven charge states while reading one, while spanning them all would cut the speed
    // win from ~900x to ~4x. Two shapes stay LEGAL and are deliberately not caught here:
    //   - `separate` at BOTH levels: it fans out to one anchor per scan (buildMS2's notch guard at
    //     ScanCommandQueue.cpp:314 tests `== Multiplexed` only), so every readout is a single
    //     interval and each gets its own correct range.
    //   - the CROSS-LEVEL case -- an MS3 remaining_precursor sweep under
    //     precursor_selection.precursor_charges == multiplexed -- because stage-0 notches change
    //     WHICH precursors are fragmented, not where the MS3 readout sits; the sub-fragment scan is
    //     still one contiguous stage-1 window. Hence two pair-specific checks rather than one
    //     "multiplexed anywhere" check.
    if (level(2).exploration == ExplorationMetric::RemainingPrecursor
        && targeting_.precursor_charges == ChargeAcquisitionMode::Multiplexed)
      throw std::invalid_argument(
          "precursor_selection.exploration.metric is \"remaining_precursor\" but "
          "precursor_selection.precursor_charges is \"multiplexed\". A multiplexed scan reads "
          "several non-contiguous isolation windows, which cannot be expressed as the one scan "
          "range such a sweep's pre-scans are bound to. Set precursor_charges to \"single\" or "
          "\"separate\", or pick a different exploration metric (ADR-0026).");

    // Guarded on `mode` for the reason spelled out at the level-3 arm above: levels_[3].exploration
    // outlives characterization.mode == "off" (applyCharacterizationMode_ resets only `selection`),
    // and under "off" no MS3 pre-scan is ever emitted or narrowed, so there is no bound readout for a
    // multiplexed fragment isolation to conflict with (ADR-0013).
    if (characterization_.mode != CharacterizationMode::Off
        && level(3).exploration == ExplorationMetric::RemainingPrecursor
        && characterization_.fragment_charges == ChargeAcquisitionMode::Multiplexed)
      throw std::invalid_argument(
          "characterization.exploration.metric is \"remaining_precursor\" but "
          "characterization.fragment_charges is \"multiplexed\". A multiplexed scan reads several "
          "non-contiguous isolation windows, which cannot be expressed as the one scan range such "
          "a sweep's pre-scans are bound to. Set fragment_charges to \"single\" or \"separate\", "
          "or pick a different exploration metric (ADR-0026).");

    // Re-keyed onto `mode`. It used to fire off the UPSTREAM gate (any level >= 2 selecting), which
    // is why 17 test configs that run no MS3 at all had to carry a placeholder "SEQUENCE": their
    // ms2.selection defaulted to intensity and nothing downstream ever read the sequence. Now the
    // requirement tracks the thing that actually consumes it, and those 17 placeholders are blanked --
    // an empty sequence is the honest encoding of "no protein was supplied", and a non-empty one is
    // load-bearing rather than decorative.
    if (characterization_.mode != CharacterizationMode::Off
        && characterization_.protein_sequence.empty())
      throw std::invalid_argument(
          "characterization.mode is not \"off\" but characterization.protein_sequence is empty. "
          "MS3 characterization matches fragments against that sequence.");


    // The converse of "mode: off does not forbid ms_settings.ms3" -- and the direction that
    // segfaults. Exploration::initiateNextLevel reads next_cfg.scans[0] unguarded, so MS3 being
    // reachable with no level-3 scan config is an OOB read, not a no-op.
    if (characterization_.mode != CharacterizationMode::Off && level(3).scans.empty())
      throw std::invalid_argument(
          "characterization.mode is \"" +
          std::string(characterizationModeName(characterization_.mode))
          + "\" but ms_settings.ms3 is not defined. An MS3 scan config is required to build the "
            "MS3 command into.");

    // MS2 needs somewhere to dispatch into whenever MS1 selects at all.
    if (level(1).selection != SelectionMetric::None && level(2).scans.empty())
      throw std::invalid_argument(
          "precursor_selection.rank_by is not \"none\" but ms_settings.ms2 is not defined.");

    // Validate that each exploration activation type has its required sweep range
    for (const auto& [lvl, cfg] : levels_)
    {
      if (cfg.exploration == ExplorationMetric::None) continue;

      // If activations is empty, default to base scan config activation
      std::vector<std::string> acts = cfg.activations;
      if (acts.empty() && !cfg.scans.empty())
        acts.push_back(cfg.scans[0].activation);

      // Step guards. These were validated on NEITHER side, and a non-positive step is not a
      // no-op: Exploration.cpp's `for (ce = ce_min; ce <= ce_max + 1e-9; ce += ce_step)` never
      // advances, so it spins forever INSIDE processScan -- on the C# TPL ActionBlock thread, with
      // the instrument still waiting for commands. A hang, not an error.
      const std::string sect = (lvl == 2 ? "precursor_selection" : "characterization");
      if (cfg.ce_step <= 0.0)
        throw std::invalid_argument(
            sect + ".exploration.ce_step must be > 0; got " + std::to_string(cfg.ce_step)
            + ". A non-positive step never terminates the sweep loop.");
      if (cfg.rt_max > cfg.rt_min && cfg.rt_step <= 0.0)
        throw std::invalid_argument(
            sect + ".exploration.reaction_time_step must be > 0 when a reaction-time range is set; got "
            + std::to_string(cfg.rt_step) + ". A non-positive step never terminates the sweep loop.");

      for (const auto& act : acts)
      {
        bool needs_ce = needsCollisionEnergy(act);
        bool needs_rt = needsReactionTime(act);

        if (needs_ce && cfg.ce_max <= cfg.ce_min)
          throw std::invalid_argument(
              "Exploration activation '" + act + "' at level " + std::to_string(lvl) +
              " requires ce_min < ce_max for CE sweep.");

        if (needs_rt && cfg.rt_max <= cfg.rt_min)
          throw std::invalid_argument(
              "Exploration activation '" + act + "' at level " + std::to_string(lvl) +
              " requires rt_min < rt_max for RT sweep.");

        // The sweep GRID is authored, not floored (ADR-0029). Exploration raises only its own
        // synthesized baseline to MIN_REACTION_TIME_MS, so a grid that starts below the floor would
        // put a scan the instrument REJECTS into every sweep -- and, because the baseline would no
        // longer coincide with the grid's first point, would also resurrect the duplicate reference
        // scan the suppression rule exists to remove. One wrong scan per group, visible only on the
        // hardware. Reject it at load instead, and name the value that works.
        //
        // Gated on needs_rt: a CE-only sweep leaves reaction_time_min at its 0 default and must not
        // be rejected for it, which is every committed config except the ETD one.
        if (needs_rt && cfg.rt_min < MIN_REACTION_TIME_MS)
          throw std::invalid_argument(
              sect + ".exploration.reaction_time_min must be >= 0.03 when an ETD-family activation "
              "is swept ('" + act + "' at level " + std::to_string(lvl) + "); got "
              + std::to_string(cfg.rt_min) + ". The instrument rejects a reaction time of 0 -- use "
              "0.03 for a near-zero sweep point.");
      }
    }
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

  DoubleList Config::toleranceList() const
  {
    DoubleList tol;
    for (const auto& [lvl, cfg] : levels_)
      tol.push_back(cfg.tolerance_ppm);
    return tol;
  }

  DoubleList Config::explorationToleranceList() const
  {
    DoubleList tol;
    for (const auto& [lvl, cfg] : levels_)
      tol.push_back(cfg.exploration_tolerance_ppm);
    return tol;
  }

} // namespace OpenMS
