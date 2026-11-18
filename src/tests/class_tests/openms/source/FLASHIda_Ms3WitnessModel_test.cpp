// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------
//
// ADR-0079: the shipped W model (Ms3WitnessModel.h, package 17's w_v13.json): its structure and SHA-256, the refusals of a
// model the engine cannot use, the drift guard between the two share copies, and the input transforms and the helpers that
// read the inputs as the logs print them. The contract values of package 17 (12,440 recorded MS3) are tested by
// FLASHIda_Ms3WGateV17_test.
//
// Reads (cwd OpenMS/build): ../share/OpenMS/FLASHIDA/ms3_w_gate.json, ../../FlashIDA/share/OpenMS/FLASHIDA/ms3_w_gate.json.
// Writes nothing.

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3WitnessModel.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3WitnessGate.h>   // the inputs read off the logs

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <nlohmann/json.hpp>

using namespace OpenMS;

namespace
{
  const std::string kModel = "../share/OpenMS/FLASHIDA/ms3_w_gate.json";
  const std::string kModelFlashIda = "../../FlashIDA/share/OpenMS/FLASHIDA/ms3_w_gate.json";

  std::string slurp(const std::string& path)
  {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream buf;
    buf << in.rdbuf();
    return buf.str();
  }

  // The shipped model's JSON with one change, for the refusal cases.
  std::string mutated(void (*change)(nlohmann::json&))
  {
    nlohmann::json j = nlohmann::json::parse(slurp(kModel));
    change(j);
    return j.dump();
  }

  bool refuses(const std::string& bytes, const std::string& needle, bool checked = false)
  {
    try
    {
      if (checked) WitnessModel::parseChecked(bytes, "probe.json");
      else WitnessModel::parse(bytes, "probe.json");
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
  TEST_EQUAL(m.treeCount(), 658)
  TEST_EQUAL(m.numFeatures(), 38)
  TEST_EQUAL(m.sha8(), "5d3f44a9")   // SHA-256 5d3f44a9... of package 17's w_v13.json
  TEST_EQUAL(m.sha256(), std::string(kWitnessModelSha256))
  // Through OPENMS_DATA_PATH, as the engine finds it, with its SHA-256 checked.
  const WitnessModel shipped = WitnessModel::loadShipped();
  TEST_EQUAL(shipped.sha8(), "5d3f44a9")
  TEST_EQUAL(shipped.origin().find("FLASHIDA") != std::string::npos, true)
}
END_SECTION

START_SECTION(a_wrong_model_is_refused)
{
  TEST_EXCEPTION(std::invalid_argument, WitnessModel::load("w_model_that_does_not_exist.json"))
  const bool not_json = refuses("{not json", "not JSON");
  const bool objective = refuses(mutated([](nlohmann::json& j) { j["learner"]["objective"]["name"] = "reg:squarederror"; }),
                                 "binary:logistic");
  const bool features = refuses(mutated([](nlohmann::json& j) { j["learner"]["learner_model_param"]["num_feature"] = "37"; }),
                                "num_feature");
  const bool tree_count = refuses(mutated([](nlohmann::json& j) {
                                    j["learner"]["gradient_booster"]["model"]["gbtree_model_param"]["num_trees"] = "657";
                                  }),
                                  "num_trees");
  const bool split_index = refuses(mutated([](nlohmann::json& j) {
                                     j["learner"]["gradient_booster"]["model"]["trees"][0]["split_indices"][0] = 38;
                                   }),
                                   "splits on input 38");
  const bool base_score = refuses(mutated([](nlohmann::json& j) { j["learner"]["learner_model_param"].erase("base_score"); }),
                                  "base_score");
  const bool categorical = refuses(mutated([](nlohmann::json& j) {
                                     j["learner"]["gradient_booster"]["model"]["trees"][3]["categories_nodes"] = nlohmann::json::array({0});
                                   }),
                                   "categorical");
  // ADR-0079: a structurally sound model that is not the package's file is refused by the checked loader.
  const bool other_file = refuses(mutated([](nlohmann::json& j) {
                                    j["learner"]["gradient_booster"]["model"]["trees"][0]["split_conditions"][0] = 0.123;
                                  }),
                                  "SHA-256", true);
  TEST_EQUAL(not_json, true)
  TEST_EQUAL(objective, true)
  TEST_EQUAL(features, true)
  TEST_EQUAL(tree_count, true)
  TEST_EQUAL(split_index, true)
  TEST_EQUAL(base_score, true)
  TEST_EQUAL(categorical, true)
  TEST_EQUAL(other_file, true)
  TEST_EQUAL(WitnessModel::parseChecked(slurp(kModel), "probe.json").treeCount(), 658)   // the file itself passes
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
  const double nan = std::numeric_limits<double>::quiet_NaN();
  std::array<double, kW17RawCount> raw{};
  raw.fill(2.0);
  raw[kW17FragN] = 1.0;
  raw[kW17ZMs3] = 7.0;
  raw[kW17L] = 44.0;
  raw[kW17Intensity] = 5974.8;
  raw[kW17IdScore] = 0.5;            // below 1: log10(max(v, 1)) = 0
  raw[kW17PrecIntensity] = 250163.0;
  raw[kW17Ms1Plan] = 0.0;            // absent from the last MS1: 0, logged as 0
  raw[kW17Ms2MassCount] = 22.0;
  raw[kW17Ms2Ions] = 9.0;
  const std::array<double, 8> pk3 = {12.5, 3.25, 0.5, 1.0, 0.75, 4.125, -0.3333333333333333, 6.0};
  const std::array<float, kW17InputCount> x = w17ModelInputs(raw, pk3);
  TEST_EQUAL(x[0], 1.0f)
  TEST_EQUAL(x[3], toFloat32(7.0 / 44.0))                       // z_ms3_per_res, in 64 bits then 32
  TEST_EQUAL(x[6], toFloat32(std::log10(5974.8)))               // log10 of the raw intensity
  TEST_EQUAL(x[7], 44.0f)
  TEST_EQUAL(x[8], 0.0f)                                        // id_score 0.5 -> log10(1)
  TEST_EQUAL(x[15], toFloat32(std::log10(static_cast<double>(toFloat32(250163.0)))))   // 32 bits first, then log10
  TEST_EQUAL(x[19], 0.0f)                                       // x_ms1_plan 0 -> log10(1)
  TEST_EQUAL(x[24], toFloat32(std::log10(22.0)))                // x_ms2_mass_count
  TEST_EQUAL(x[29], 9.0f)                                       // x_ms2_ions, the last log input
  TEST_EQUAL(x[30], 12.5f)                                      // pk3_since_start
  TEST_EQUAL(x[36], toFloat32(-0.3333333333333333))             // pk3_trend, rounded to 32 bits at the trees
  TEST_EQUAL(x[37], 6.0f)                                       // pk3_n_seen, the last column
  TEST_EQUAL(std::string(w17InputName(30)), "pk3_since_start")
  TEST_EQUAL(std::string(w17InputName(37)), "pk3_n_seen")
  TEST_EQUAL(std::string(w17RawName(kW17Ms1Plan)), "x_ms1_plan")
  // A missing L leaves z_ms3_per_res missing; a missing value stays NaN through every transform.
  raw[kW17L] = nan;
  raw[kW17IdScore] = nan;
  const std::array<float, kW17InputCount> y = w17ModelInputs(raw, pk3);
  TEST_EQUAL(std::isnan(y[3]), true)
  TEST_EQUAL(std::isnan(y[7]), true)
  TEST_EQUAL(std::isnan(y[8]), true)
  TEST_EQUAL(std::isinf(toFloat32(1e40)), true)

  // The inputs read off the logs, by the package's definitions (inputs_from_log_v17.py).
  // has_shift: a '[' after the leading <...> global modification; the <...> itself is not a shift.
  TEST_EQUAL(proformaHasShift("<[Carbamidomethyl]@C>PEPCK"), false)
  TEST_EQUAL(proformaHasShift("<[Carbamidomethyl]@C>PEP[+15.995]CK"), true)
  TEST_EQUAL(proformaHasShift("[Acetyl]-GDVEK"), true)
  TEST_EQUAL(proformaHasShift("GDVEK"), false)
  // L: the capital letters after removing <...> and every [...].
  TEST_EQUAL(proformaResidueCount("<[Carbamidomethyl]@C>[Acetyl]-GD(VE)[+1.2]K"), 5)
  TEST_EQUAL(proformaResidueCount("PEPTK"), 5)
  TEST_EQUAL(proformaResidueCount(""), 0)
  // Stage 0 / stage 1 of a logged cell: the first / last ';'-group, its first ','-value (the anchor of notches).
  TEST_EQUAL(loggedStage("11,12;3", 0), 11.0)
  TEST_EQUAL(loggedStage("11,12;3", -1), 3.0)
  TEST_EQUAL(loggedStage("29;35", 0), 29.0)
  TEST_EQUAL(loggedStage("13569.1235;2010.1310", -1), 2010.131)
  TEST_EQUAL(std::isnan(loggedStage("", 0)), true)
  // Printed precision: a stage-scored column at 6 significant digits, a mass (and id_score) at 4 decimals.
  TEST_EQUAL(asLoggedSig6(250163.37), 250163.0)
  TEST_EQUAL(asLoggedSig6(0.123456789), 0.123457)
  TEST_EQUAL(asLoggedSig6(1234567.0), 1.23457e6)
  TEST_EQUAL(asLoggedFixed4(12357.123456), 12357.1235)
  TEST_EQUAL(std::isnan(asLoggedSig6(nan)), true)
  // x_ms1_plan: the largest summed intensity within 10 ppm of the mass (of the precursor mass) at 0 or +-1.00335 Da.
  WitnessLastMs1 ms1;
  ms1.masses = {{12357.1, 4.0e6}, {12358.10335, 9.0e6}, {12355.5, 7.0e6}, {12357.1 * (1 + 2e-5), 8.0e6}};
  TEST_EQUAL(witnessMs1Plan(ms1, 12357.1), 9.0e6)                // the +1.00335 isotope is the larger
  ms1.masses.erase(ms1.masses.begin() + 1);
  TEST_EQUAL(witnessMs1Plan(ms1, 12357.1), 4.0e6)                // 20 ppm off does not count, nor 1.6 Da off
  TEST_EQUAL(witnessMs1Plan(WitnessLastMs1{}, 12357.1), 0.0)
}
END_SECTION

END_TEST
