// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms3WitnessModel.h>

#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QByteArray>
#include <QtCore/QCryptographicHash>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <nlohmann/json.hpp>

namespace OpenMS
{
  namespace
  {
    const double kNaN = std::numeric_limits<double>::quiet_NaN();

    /// log10(max(value, 1)) of the value rounded to 32 bits, rounded to 32 bits (w_gate.py log10_floor1).
    float log10Floor1(double v)
    {
      const float f = toFloat32(v);
      if (std::isnan(f)) return std::numeric_limits<float>::quiet_NaN();
      return toFloat32(std::log10(std::max(static_cast<double>(f), 1.0)));
    }

    [[noreturn]] void refuse(const std::string& origin, const std::string& what)
    {
      throw std::invalid_argument("FLASHIda: the W model " + origin + " is not usable: " + what
                                  + " (ADR-0079; share/OpenMS/FLASHIDA/ms3_w_gate.json is a byte copy of the ms3_esm "
                                    "package 17's w_v13.json).");
    }

    /// The member @p key of the object @p j; refuses when @p j is not an object or has no such member.
    const nlohmann::json& child(const nlohmann::json& j, const char* key, const std::string& origin)
    {
      if (!j.is_object() || !j.contains(key)) refuse(origin, std::string("no ") + key);
      return j.at(key);
    }

    /// An integer the XGBoost JSON writes as a decimal string ("36").
    long stringInt(const nlohmann::json& j, const std::string& origin, const char* what)
    {
      if (!j.is_string()) refuse(origin, std::string(what) + " is not a string");
      const std::string s = j.get<std::string>();
      long v = 0;
      bool ok = false;
      try
      {
        size_t used = 0;
        v = std::stol(s, &used);
        ok = (used == s.size());
      }
      catch (const std::exception&)
      {
        ok = false;
      }
      if (!ok) refuse(origin, std::string(what) + " is not an integer: " + s);
      return v;
    }
  } // namespace

  float toFloat32(double v)
  {
    if (std::isnan(v)) return std::numeric_limits<float>::quiet_NaN();
    if (std::fabs(v) > static_cast<double>(std::numeric_limits<float>::max()))
    {
      // Round-to-nearest still maps values just past FLT_MAX (below the halfway point to 2^128) onto FLT_MAX; anything
      // beyond overflows, as struct.pack('<f') does in w_gate.py. No input of the fit comes near.
      const double halfway = static_cast<double>(std::numeric_limits<float>::max())
                             + std::ldexp(1.0, std::numeric_limits<float>::max_exponent - std::numeric_limits<float>::digits - 1);
      if (std::fabs(v) < halfway) return v < 0 ? -std::numeric_limits<float>::max() : std::numeric_limits<float>::max();
      return v < 0 ? -std::numeric_limits<float>::infinity() : std::numeric_limits<float>::infinity();
    }
    return static_cast<float>(v);
  }

  namespace
  {
    const char* const kW17InputNames[kW17InputCount] = {
        "frag_N", "z_ms2", "z_ms3", "z_ms3_per_res", "ce1", "cid3", "log10_intensity", "L", "id_score", "iso_cos", "qscore", "snr",
        "has_shift", "mass_ms2", "mass_ms3", "x_prec_intensity", "x_prec_qscore", "x_prec_iso_cos", "x_prec_snr", "x_ms1_plan",
        "x_charge_snr", "x_charge_cos", "x_charge_score", "x_ppm_error", "x_ms2_mass_count", "x_ms2_tic_coverage",
        "x_ms2_precursor_share", "x_ms2_fragment_count", "x_ms2_tag_count", "x_ms2_ions", "pk3_since_start", "pk3_since_apex",
        "pk3_since_last", "pk3_missed", "pk3_rel_max", "pk3_log_max", "pk3_trend", "pk3_n_seen"};
    const char* const kW17RawNames[kW17RawCount] = {
        "frag_N", "z_ms2", "z_ms3", "ce1", "cid3", "intensity", "L", "id_score", "iso_cos", "qscore", "snr", "has_shift", "mass_ms2",
        "mass_ms3", "x_prec_intensity", "x_prec_qscore", "x_prec_iso_cos", "x_prec_snr", "x_ms1_plan", "x_charge_snr", "x_charge_cos",
        "x_charge_score", "x_ppm_error", "x_ms2_mass_count", "x_ms2_tic_coverage", "x_ms2_precursor_share", "x_ms2_fragment_count",
        "x_ms2_tag_count", "x_ms2_ions"};
  } // namespace

  const char* w17InputName(int index) { return (index >= 0 && index < kW17InputCount) ? kW17InputNames[index] : "?"; }
  const char* w17RawName(int index) { return (index >= 0 && index < kW17RawCount) ? kW17RawNames[index] : "?"; }

  std::array<float, kW17InputCount> w17ModelInputs(const std::array<double, kW17RawCount>& raw, const std::array<double, 8>& pk3)
  {
    std::array<float, kW17InputCount> x{};
    const double z3 = raw[kW17ZMs3];
    const double L = raw[kW17L];
    const double inten = raw[kW17Intensity];
    x[0] = toFloat32(raw[kW17FragN]);
    x[1] = toFloat32(raw[kW17ZMs2]);
    x[2] = toFloat32(z3);
    x[3] = (!std::isnan(L) && L != 0.0) ? toFloat32(z3 / L) : std::numeric_limits<float>::quiet_NaN();
    x[4] = toFloat32(raw[kW17Ce1]);
    x[5] = toFloat32(raw[kW17Cid3]);
    // log10_intensity: of the raw value as given, not rounded to 32 bits first (w_gate_v17.model_inputs).
    x[6] = std::isnan(inten) ? std::numeric_limits<float>::quiet_NaN() : toFloat32(std::log10(std::max(inten, 1.0)));
    x[7] = toFloat32(L);
    // Raw slots kW17IdScore .. kW17Ms2Ions are model columns 8 .. 29, in the same order.
    for (int r = kW17IdScore; r < kW17RawCount; ++r)
    {
      const bool log10_input = (r == kW17IdScore || r == kW17PrecIntensity || r == kW17Ms1Plan || r == kW17Ms2MassCount);
      x[static_cast<size_t>(r + 1)] = log10_input ? log10Floor1(raw[static_cast<size_t>(r)]) : toFloat32(raw[static_cast<size_t>(r)]);
    }
    for (size_t k = 0; k < 8; ++k) x[30 + k] = toFloat32(pk3[k]);
    return x;
  }

  WitnessModel WitnessModel::loadShipped()
  {
    String path;
    try
    {
      path = File::find("FLASHIDA/ms3_w_gate.json");
    }
    catch (const std::exception&)
    {
      refuse("FLASHIDA/ms3_w_gate.json", "not found under OPENMS_DATA_PATH (" + File::getOpenMSDataPath() + ")");
    }
    std::ifstream in(path, std::ios::binary);
    if (!in) refuse(path, "the file cannot be opened");
    std::ostringstream buf;
    buf << in.rdbuf();
    return parseChecked(buf.str(), path);
  }

  WitnessModel WitnessModel::parseChecked(const std::string& bytes, const std::string& origin)
  {
    WitnessModel m = parse(bytes, origin, kWitnessInputCount);
    if (m.sha256() != kWitnessModelSha256)
      refuse(origin, "its SHA-256 is " + m.sha256() + ", not " + std::string(kWitnessModelSha256)
                         + " (package 17's w_v13.json)");
    return m;
  }

  WitnessModel WitnessModel::load(const std::string& path) { return load(path, kWitnessInputCount); }

  WitnessModel WitnessModel::load(const std::string& path, int num_features)
  {
    std::ifstream in(path, std::ios::binary);
    if (!in) refuse(path, "the file cannot be opened");
    std::ostringstream buf;
    buf << in.rdbuf();
    return parse(buf.str(), path, num_features);
  }

  WitnessModel WitnessModel::parse(const std::string& bytes, const std::string& origin, int num_features)
  {
    nlohmann::json j;
    try
    {
      j = nlohmann::json::parse(bytes);
    }
    catch (const std::exception& e)
    {
      refuse(origin, std::string("not JSON: ") + e.what());
    }
    const nlohmann::json& learner = child(j, "learner", origin);

    // The structure checks; the shipped loader also checks the file's SHA-256 (ADR-0079).
    const nlohmann::json& objective_name = child(child(learner, "objective", origin), "name", origin);
    if (!objective_name.is_string() || objective_name.get<std::string>() != "binary:logistic")
      refuse(origin, "the objective is not binary:logistic");
    const nlohmann::json& lmp = child(learner, "learner_model_param", origin);
    if (stringInt(child(lmp, "num_feature", origin), origin, "num_feature") != num_features)
      refuse(origin, "num_feature is not " + std::to_string(num_features));
    const nlohmann::json& base_json = child(lmp, "base_score", origin);
    if (!base_json.is_string()) refuse(origin, "base_score is not a string");
    std::string base_text = base_json.get<std::string>();   // "[4.1593963E-1]" in XGBoost 3
    base_text.erase(std::remove(base_text.begin(), base_text.end(), '['), base_text.end());
    base_text.erase(std::remove(base_text.begin(), base_text.end(), ']'), base_text.end());
    double base = 0.0;
    try
    {
      base = std::stod(base_text);
    }
    catch (const std::exception&)
    {
      refuse(origin, "base_score is not a number: " + base_text);
    }
    if (!(base > 0.0 && base < 1.0)) refuse(origin, "base_score is not in (0, 1)");

    const nlohmann::json& booster = child(learner, "gradient_booster", origin);
    const nlohmann::json& booster_name = child(booster, "name", origin);
    if (!booster_name.is_string() || booster_name.get<std::string>() != "gbtree") refuse(origin, "the booster is not gbtree");
    const nlohmann::json& model = child(booster, "model", origin);
    const nlohmann::json& trees = child(model, "trees", origin);
    if (!trees.is_array()) refuse(origin, "trees is not an array");
    const long num_trees = stringInt(child(child(model, "gbtree_model_param", origin), "num_trees", origin), origin, "num_trees");
    if (num_trees != static_cast<long>(trees.size()) || trees.empty())
      refuse(origin, "num_trees " + std::to_string(num_trees) + " does not match the " + std::to_string(trees.size())
                         + " trees in the file");

    WitnessModel m;
    m.origin_ = origin;
    m.num_features_ = num_features;
    m.intercept_ = std::log(base / (1.0 - base));
    m.trees_.reserve(trees.size());
    for (size_t t = 0; t < trees.size(); ++t)
    {
      const nlohmann::json& tj = trees[t];
      const std::string where = "tree " + std::to_string(t) + ": ";
      if (!tj.is_object()) refuse(origin, where + "not an object");
      if (tj.contains("categories_nodes") && !tj.at("categories_nodes").empty())
        refuse(origin, where + "categorical splits are not supported");
      Tree tr;
      try
      {
        tr.left = tj.at("left_children").get<std::vector<int32_t>>();
        tr.right = tj.at("right_children").get<std::vector<int32_t>>();
        tr.feature = tj.at("split_indices").get<std::vector<int32_t>>();
        tr.default_left = tj.at("default_left").get<std::vector<uint8_t>>();
        // Thresholds and leaf values as 32-bit floats: the file prints each in its shortest decimal form, which reads
        // back to the same float (ENGINE.md section 4).
        for (double c : tj.at("split_conditions").get<std::vector<double>>()) tr.value.push_back(toFloat32(c));
      }
      catch (const std::exception& e)
      {
        refuse(origin, where + "malformed: " + e.what());
      }
      const size_t n = tr.left.size();
      if (n == 0 || tr.right.size() != n || tr.feature.size() != n || tr.value.size() != n || tr.default_left.size() != n)
        refuse(origin, where + "the node arrays differ in length");
      for (size_t i = 0; i < n; ++i)
      {
        const bool leaf = tr.left[i] == -1;
        if (leaf != (tr.right[i] == -1)) refuse(origin, where + "node " + std::to_string(i) + " has one child");
        if (leaf) continue;
        if (tr.left[i] <= static_cast<int32_t>(i) || tr.right[i] <= static_cast<int32_t>(i)
            || tr.left[i] >= static_cast<int32_t>(n) || tr.right[i] >= static_cast<int32_t>(n))
          refuse(origin, where + "node " + std::to_string(i) + " has a child out of order or out of range");
        if (tr.feature[i] < 0 || tr.feature[i] >= num_features)
          refuse(origin, where + "node " + std::to_string(i) + " splits on input " + std::to_string(tr.feature[i])
                             + ", not one of the " + std::to_string(num_features));
      }
      m.trees_.push_back(std::move(tr));
    }

    QCryptographicHash hash(QCryptographicHash::Sha256);
    hash.addData(QByteArray::fromRawData(bytes.data(), static_cast<int>(bytes.size())));
    m.sha256_ = hash.result().toHex().toStdString();
    m.sha8_ = m.sha256_.substr(0, 8);
    return m;
  }

  double WitnessModel::wOf(const float* x) const
  {
    double sum = intercept_;
    for (const Tree& t : trees_)
    {
      size_t i = 0;
      while (t.left[i] != -1)
      {
        const float v = x[static_cast<size_t>(t.feature[i])];
        if (std::isnan(v))
          i = static_cast<size_t>(t.default_left[i] ? t.left[i] : t.right[i]);
        else
          i = static_cast<size_t>(v < t.value[i] ? t.left[i] : t.right[i]);
      }
      sum += static_cast<double>(t.value[i]);
    }
    return 1.0 / (1.0 + std::exp(-sum));
  }
} // namespace OpenMS
