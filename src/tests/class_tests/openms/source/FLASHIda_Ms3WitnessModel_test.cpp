// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// ADR-0076: the witness probability W (Ms3WitnessModel.h) against the ms3_esm package's acceptance test (ENGINE.md
// section 7): the 300 reference MS3 and the 200 whose inputs sit exactly on a split threshold, from the raw values to W,
// to 1e-6, and the decisions of the package's two threshold variants. Plus the shipped file's structure, the refusals
// of a model the engine cannot use, the drift guard between the two share copies, and the input transforms.
//
// Reads (cwd OpenMS/build): ../share/OpenMS/FLASHIDA/ms3_w_gate.json, ../../FlashIDA/share/OpenMS/FLASHIDA/ms3_w_gate.json,
// ../../FlashIDA/test-data/ms3_w_gate/reference_rows{,_ties}.tsv. Writes nothing.

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3WitnessModel.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3WitnessGate.h>   // the inputs read off the logs

#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

using namespace OpenMS;

namespace
{
  const std::string kModel = "../share/OpenMS/FLASHIDA/ms3_w_gate.json";
  const std::string kModelFlashIda = "../../FlashIDA/share/OpenMS/FLASHIDA/ms3_w_gate.json";
  const std::string kRows = "../../FlashIDA/test-data/ms3_w_gate/reference_rows.tsv";
  const std::string kTies = "../../FlashIDA/test-data/ms3_w_gate/reference_rows_ties.tsv";
  // MANIFEST.json gate.variants (neither is chosen; they are the package's acceptance points).
  const double kW50 = 0.5;
  const double kW10 = 0.835;

  std::string slurp(const std::string& path)
  {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream buf;
    buf << in.rdbuf();
    return buf.str();
  }

  // A TSV as header -> column, rows of cells; CR stripped (the fixtures are CRLF on disk).
  struct Table
  {
    std::map<std::string, size_t> col;
    std::vector<std::vector<std::string>> rows;
    const std::string& at(size_t r, const std::string& name) const { return rows[r][col.at(name)]; }
  };

  Table readTsv(const std::string& path)
  {
    Table t;
    std::ifstream in(path);
    std::string line;
    bool header = true;
    while (std::getline(in, line))
    {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty()) continue;
      std::vector<std::string> cells;
      std::stringstream ss(line);
      std::string cell;
      while (std::getline(ss, cell, '\t')) cells.push_back(cell);
      if (line.back() == '\t') cells.push_back("");
      if (header)
      {
        for (size_t i = 0; i < cells.size(); ++i) t.col[cells[i]] = i;
        header = false;
      }
      else
        t.rows.push_back(cells);
    }
    return t;
  }

  // One reference row's raw values, as w_gate.model_inputs reads them: frag_type N/C, every other value by name.
  WitnessRaw rawOf(const Table& t, size_t r)
  {
    WitnessRaw raw;
    const std::string& ft = t.at(r, "frag_type");
    raw[kWFragN] = (ft == "N") ? 1.0 : (ft == "C" ? 0.0 : std::numeric_limits<double>::quiet_NaN());
    for (int i = 1; i < kWitnessRawCount; ++i)
    {
      const std::string& s = t.at(r, witnessRawName(i));
      raw[i] = (s.empty() || s == "nan") ? std::numeric_limits<double>::quiet_NaN() : std::stod(s);
    }
    return raw;
  }

  bool decide(double w, const std::string& n_uncovered, double threshold)
  {
    return std::stod(n_uncovered) > 0 && witnessAdmits(w, threshold);
  }

  // The shipped model's JSON with one change, for the refusal cases.
  std::string mutated(void (*change)(nlohmann::json&))
  {
    nlohmann::json j = nlohmann::json::parse(slurp(kModel));
    change(j);
    return j.dump();
  }

  bool refuses(const std::string& bytes, const std::string& needle)
  {
    try
    {
      WitnessModel::parse(bytes, "probe.json");
    }
    catch (const std::invalid_argument& e)
    {
      const std::string what = e.what();
      return what.find("probe.json") != std::string::npos && what.find(needle) != std::string::npos;
    }
    return false;
  }
} // namespace

START_TEST(FLASHIda_Ms3WitnessModel, "$Id$")

START_SECTION(shipped_model_loads)
{
  const WitnessModel m = WitnessModel::load(kModel);
  TEST_EQUAL(m.treeCount(), 84)
  TEST_EQUAL(m.sha8(), "47f0babe")   // SHA-256 47f0babed8f4... of the package's w_ext_v3.json
  // Through OPENMS_DATA_PATH, as the engine finds it.
  const WitnessModel shipped = WitnessModel::loadShipped();
  TEST_EQUAL(shipped.sha8(), "47f0babe")
  TEST_EQUAL(shipped.origin().find("FLASHIDA") != std::string::npos, true)
}
END_SECTION

START_SECTION(reference_rows_reproduce_w)
{
  const WitnessModel m = WitnessModel::load(kModel);
  const Table t = readTsv(kRows);
  TEST_EQUAL(t.rows.size(), 300)
  double worst = 0.0;
  int decisions_differ = 0;
  for (size_t r = 0; r < t.rows.size(); ++r)
  {
    const double w = m.w(rawOf(t, r));
    worst = std::max(worst, std::fabs(w - std::stod(t.at(r, "w_expected"))));
    const std::string& unc = t.at(r, "n_uncovered_bonds");
    decisions_differ += decide(w, unc, kW50) != (t.at(r, "send_w50") == "1");
    decisions_differ += decide(w, unc, kW10) != (t.at(r, "send_w10") == "1");
  }
  STATUS("largest |W - w_expected| over 300 rows: " << worst)
  TEST_EQUAL(worst < 1e-6, true)
  TEST_EQUAL(decisions_differ, 0)
}
END_SECTION

START_SECTION(tie_rows_reproduce_w)
{
  const WitnessModel m = WitnessModel::load(kModel);
  const Table t = readTsv(kTies);
  TEST_EQUAL(t.rows.size(), 200)
  double worst = 0.0;
  int decisions_differ = 0;
  int wrong_precision_agrees = 0;
  for (size_t r = 0; r < t.rows.size(); ++r)
  {
    const double w = m.w(rawOf(t, r));
    worst = std::max(worst, std::fabs(w - std::stod(t.at(r, "w_expected"))));
    const std::string& unc = t.at(r, "n_uncovered_bonds");
    decisions_differ += decide(w, unc, kW50) != (t.at(r, "send_w50") == "1");
    decisions_differ += decide(w, unc, kW10) != (t.at(r, "send_w10") == "1");
    // Every tie row was chosen because the 64-bit decimal comparison gives it another W: ours must not be that one.
    wrong_precision_agrees += std::fabs(w - std::stod(t.at(r, "w_if_thresholds_are_compared_as_64_bit_decimals"))) <= 1e-6;
  }
  STATUS("largest |W - w_expected| over 200 tie rows: " << worst)
  TEST_EQUAL(worst < 1e-6, true)
  TEST_EQUAL(decisions_differ, 0)
  TEST_EQUAL(wrong_precision_agrees, 0)
}
END_SECTION

START_SECTION(a_wrong_model_is_refused)
{
  TEST_EXCEPTION(std::invalid_argument, WitnessModel::load("w_model_that_does_not_exist.json"))
  const bool not_json = refuses("{not json", "not JSON");
  const bool objective = refuses(mutated([](nlohmann::json& j) { j["learner"]["objective"]["name"] = "reg:squarederror"; }),
                                 "binary:logistic");
  const bool features = refuses(mutated([](nlohmann::json& j) { j["learner"]["learner_model_param"]["num_feature"] = "35"; }),
                                "num_feature");
  const bool tree_count = refuses(mutated([](nlohmann::json& j) {
                                    j["learner"]["gradient_booster"]["model"]["gbtree_model_param"]["num_trees"] = "83";
                                  }),
                                  "num_trees");
  const bool split_index = refuses(mutated([](nlohmann::json& j) {
                                     j["learner"]["gradient_booster"]["model"]["trees"][0]["split_indices"][0] = 36;
                                   }),
                                   "splits on input 36");
  const bool base_score = refuses(mutated([](nlohmann::json& j) { j["learner"]["learner_model_param"].erase("base_score"); }),
                                  "base_score");
  const bool categorical = refuses(mutated([](nlohmann::json& j) {
                                     j["learner"]["gradient_booster"]["model"]["trees"][3]["categories_nodes"] = nlohmann::json::array({0});
                                   }),
                                   "categorical");
  TEST_EQUAL(not_json, true)
  TEST_EQUAL(objective, true)
  TEST_EQUAL(features, true)
  TEST_EQUAL(tree_count, true)
  TEST_EQUAL(split_index, true)
  TEST_EQUAL(base_score, true)
  TEST_EQUAL(categorical, true)
}
END_SECTION

START_SECTION(share_copies_are_identical)
{
  // The engine's ctest reads OpenMS/share, every C# process FlashIDA/share: one model, two copies, byte for byte.
  const std::string a = slurp(kModel);
  const std::string b = slurp(kModelFlashIda);
  TEST_EQUAL(a.empty(), false)
  TEST_EQUAL(a == b, true)
}
END_SECTION

START_SECTION(transforms_follow_the_training)
{
  WitnessRaw raw;
  TEST_EQUAL(raw.firstMissing(), 0)
  for (int i = 0; i < kWitnessRawCount; ++i) raw[i] = 2.0;
  TEST_EQUAL(raw.firstMissing(), -1)
  raw[kWFragN] = 1.0;
  raw[kWZMs3] = 7.0;
  raw[kWL] = 44.0;
  raw[kWIntensity] = 5974.8;
  raw[kWIdScore] = 0.5;            // below 1: log10(max(v, 1)) = 0
  raw[kWPrecIntensity] = 250163.0;
  raw[kWMs1Plan] = 0.0;            // absent from the last MS1: 0, logged as 0
  raw[kWMs2MassCount] = 22.0;
  raw[kWRtFrac] = 0.0573674492;
  const std::array<float, kWitnessInputCount> x = witnessModelInputs(raw);
  TEST_EQUAL(x[0], 1.0f)
  TEST_EQUAL(x[3], toFloat32(7.0 / 44.0))                       // z_ms3_per_res
  TEST_EQUAL(x[6], toFloat32(std::log10(5974.8)))               // log10 of the raw intensity
  TEST_EQUAL(x[7], 44.0f)
  TEST_EQUAL(x[8], 0.0f)                                        // id_score 0.5 -> log10(1)
  TEST_EQUAL(x[17], toFloat32(std::log10(static_cast<double>(toFloat32(250163.0)))))   // 32 bits first, then log10
  TEST_EQUAL(x[21], 0.0f)                                       // x_ms1_plan 0 -> log10(1)
  TEST_EQUAL(x[26], toFloat32(std::log10(22.0)))                // x_ms2_mass_count
  TEST_EQUAL(x[35], toFloat32(0.0573674492))                    // x_rt_frac, last column
  // A missing L leaves z_ms3_per_res missing; a missing value stays NaN through every transform.
  raw[kWL] = std::numeric_limits<double>::quiet_NaN();
  raw[kWIdScore] = std::numeric_limits<double>::quiet_NaN();
  const std::array<float, kWitnessInputCount> y = witnessModelInputs(raw);
  TEST_EQUAL(std::isnan(y[3]), true)
  TEST_EQUAL(std::isnan(y[7]), true)
  TEST_EQUAL(std::isnan(y[8]), true)
  TEST_EQUAL(raw.firstMissing(), static_cast<int>(kWL))
  TEST_EQUAL(std::isinf(toFloat32(1e40)), true)
  TEST_EQUAL(std::string(witnessRawName(kWRtFrac)), "x_rt_frac")
  TEST_EQUAL(std::string(witnessRawName(kWFragN)), "frag_N")

  // The inputs read off the logs, by the training's definitions (inputs_from_log.py, ms3_esm 8.2).
  // has_shift: a '[' after the leading <...> global modification; the <...> itself is not a shift.
  TEST_EQUAL(proformaHasShift("<[Carbamidomethyl]@C>PEPCK"), false)
  TEST_EQUAL(proformaHasShift("<[Carbamidomethyl]@C>PEP[+15.995]CK"), true)
  TEST_EQUAL(proformaHasShift("[Acetyl]-GDVEK"), true)
  TEST_EQUAL(proformaHasShift("GDVEK"), false)
  // The letters: the <...> removed, every [...] removed (nested ones too), every other non-letter dropped.
  TEST_EQUAL(proformaLetters("<[Carbamidomethyl]@C>[Acetyl]-GD(VE)[+1.2]K"), "GDVEK")
  TEST_EQUAL(proformaLetters("PEP[Phospho[+79.97]]TK"), "PEPTK")
  TEST_EQUAL(proformaLetters(""), "")
  // Placement: an N-terminal fragment is the first L letters of the parent's stretch, a C-terminal one the last L.
  TEST_EQUAL(witnessPlaced("GDVEKGKK", "GDV", true), true)
  TEST_EQUAL(witnessPlaced("GDVEKGKK", "GKK", false), true)
  TEST_EQUAL(witnessPlaced("GDVEKGKK", "GKK", true), false)
  TEST_EQUAL(witnessPlaced("GDVEKGKK", "DVE", true), false)     // inside the stretch, at neither terminus
  TEST_EQUAL(witnessPlaced("GDV", "GDVEK", true), false)        // longer than the stretch
  TEST_EQUAL(witnessPlaced("GDVEKGKK", "", true), false)
  // Printed precision: a stage-scored column at 6 significant digits, a mass (and id_score) at 4 decimals.
  TEST_EQUAL(asLoggedSig6(250163.37), 250163.0)
  TEST_EQUAL(asLoggedSig6(0.123456789), 0.123457)
  TEST_EQUAL(asLoggedSig6(1234567.0), 1.23457e6)
  TEST_EQUAL(asLoggedFixed4(12357.123456), 12357.1235)
  TEST_EQUAL(std::isnan(asLoggedSig6(std::numeric_limits<double>::quiet_NaN())), true)
  // x_ms1_plan: the largest summed intensity within 10 ppm of the mass at 0 or +-1.00335 Da, else 0.
  WitnessLastMs1 ms1;
  ms1.masses = {{12357.1, 4.0e6}, {12358.10335, 9.0e6}, {12355.5, 7.0e6}, {12357.1 * (1 + 2e-5), 8.0e6}};
  TEST_EQUAL(witnessMs1Plan(ms1, 12357.1), 9.0e6)                // the +1.00335 isotope is the larger
  ms1.masses.erase(ms1.masses.begin() + 1);
  TEST_EQUAL(witnessMs1Plan(ms1, 12357.1), 4.0e6)                // 20 ppm off does not count, nor 1.6 Da off
  TEST_EQUAL(witnessMs1Plan(WitnessLastMs1{}, 12357.1), 0.0)
}
END_SECTION

END_TEST
