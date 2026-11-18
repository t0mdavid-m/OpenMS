// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief The witness probability W of an MS3 (ADR-0079): the chance that the MS3 shows at least one bond of its target
    fragment, from the ms3_esm package 17's gradient-boosted trees (w_v13.json, XGBoost binary:logistic, 38 inputs, 658
    trees).

    The trees are read at runtime from share/OpenMS/FLASHIDA/ms3_w_gate.json, a byte copy of the package's w_v13.json,
    whose SHA-256 the shipped loader checks. The evaluation follows the package's SPEC section 9: every input and every
    split threshold and leaf value is a 32-bit float, a node goes left iff input < threshold (strictly) and a missing
    input goes where the node's default_left says, the leaves are summed in 64 bits, and
    W = 1 / (1 + exp(-(ln(b / (1 - b)) + sum))). The thresholds sit at values of the training data, so an input is often
    exactly equal to one: compared in any other precision, such an MS3 goes the wrong way.
  */

  /// The number of model inputs (columns of the trees) of package 17.
  constexpr int kWitnessInputCount = 38;
  /// The SHA-256 of the shipped model file (package 17's w_v13.json).
  constexpr const char* kWitnessModelSha256 = "5d3f44a99ed580c6557605ee6433ea44f4a8b6cfeac857072ca09f8224f755d3";

  /// @p v rounded to the nearest 32-bit float; NaN stays NaN, a value past the float range becomes +-infinity.
  OPENMS_DLLAPI float toFloat32(double v);

  /// ADR-0079 (package 17): the 29 raw values of one candidate, in the order of the package's w_gate_v17.RAW. Slot 0 is
  /// frag_N (1 for an a, b or c ion, 0 for x, y or z). NaN = missing.
  enum W17Raw : int
  {
    kW17FragN = 0,
    kW17ZMs2,
    kW17ZMs3,
    kW17Ce1,
    kW17Cid3,
    kW17Intensity,
    kW17L,
    kW17IdScore,
    kW17IsoCos,
    kW17Qscore,
    kW17Snr,
    kW17HasShift,
    kW17MassMs2,
    kW17MassMs3,
    kW17PrecIntensity,
    kW17PrecQscore,
    kW17PrecIsoCos,
    kW17PrecSnr,
    kW17Ms1Plan,
    kW17ChargeSnr,
    kW17ChargeCos,
    kW17ChargeScore,
    kW17PpmError,
    kW17Ms2MassCount,
    kW17Ms2TicCoverage,
    kW17Ms2PrecursorShare,
    kW17Ms2FragmentCount,
    kW17Ms2TagCount,
    kW17Ms2Ions,
    kW17RawCount   ///< 29
  };
  /// The number of model inputs of package 17: the 29 raw values with z_ms3_per_res derived, and the eight peak features.
  constexpr int kW17InputCount = 38;
  /// The package's name of model input @p index (w_gate_v17.INPUTS), "?" out of range.
  OPENMS_DLLAPI const char* w17InputName(int index);
  /// The package's name of raw value @p index (w_gate_v17.RAW; frag_N for slot 0), "?" out of range.
  OPENMS_DLLAPI const char* w17RawName(int index);

  /// The 38 model inputs (SPEC 7): frag_N; z_ms3_per_res = z_ms3 / L in 64 bits; log10_intensity = log10(max(v, 1)) in
  /// 64 bits; id_score, x_prec_intensity, x_ms1_plan, x_ms2_mass_count rounded to 32 bits, then log10(max(v, 1)) in 64
  /// bits; the eight peak features as computed; every value rounded to 32 bits at the end. A missing value stays NaN.
  OPENMS_DLLAPI std::array<float, kW17InputCount> w17ModelInputs(const std::array<double, kW17RawCount>& raw,
                                                               const std::array<double, 8>& peak_features);

  class OPENMS_DLLAPI WitnessModel
  {
  public:
    /// The shipped model: File::find("FLASHIDA/ms3_w_gate.json") through OPENMS_DATA_PATH, refused unless its SHA-256
    /// is kWitnessModelSha256. Throws std::invalid_argument naming the path and the failed check.
    static WitnessModel loadShipped();

    /// parse(), then the SHA-256 of @p bytes against kWitnessModelSha256 (what loadShipped reads with).
    static WitnessModel parseChecked(const std::string& bytes, const std::string& origin);

    /// Read and check the model file at @p path. Throws std::invalid_argument naming @p path and the failed check.
    static WitnessModel load(const std::string& path);

    /// Check and build a model from the bytes of an XGBoost JSON model file; @p origin names it in an error.
    /// @p num_features: the inputs the file must declare (ADR-0079: 38 for package 17).
    static WitnessModel parse(const std::string& bytes, const std::string& origin, int num_features = kWitnessInputCount);

    /// Read and check the model file at @p path with @p num_features inputs.
    static WitnessModel load(const std::string& path, int num_features);

    /// W for the numFeatures() model inputs at @p x.
    double wOf(const float* x) const;
    /// W for the 38 inputs of package 17 (the model must have 38 features).
    double w(const std::array<float, kW17InputCount>& x) const { return wOf(x.data()); }

    int numFeatures() const { return num_features_; }
    /// The SHA-256 of the bytes the model was built from, 64 hex digits.
    const std::string& sha256() const { return sha256_; }

    size_t treeCount() const { return trees_.size(); }
    /// The first eight hex digits of the SHA-256 of the bytes the model was built from (the [MS3-W] line's sha=).
    const std::string& sha8() const { return sha8_; }
    /// Where the model came from (its path), for the [MS3-W] line and errors.
    const std::string& origin() const { return origin_; }

  private:
    struct Tree
    {
      std::vector<int32_t> left;
      std::vector<int32_t> right;
      std::vector<int32_t> feature;
      std::vector<float> value;          ///< split threshold at an inner node, leaf value at a leaf
      std::vector<uint8_t> default_left;
    };
    std::vector<Tree> trees_;
    double intercept_ = 0.0;             ///< ln(b / (1 - b)), b = base_score
    int num_features_ = kWitnessInputCount;
    std::string sha256_;
    std::string sha8_;
    std::string origin_;
  };

  /// The gate's test: an MS3 is admitted when W >= threshold.
  inline bool witnessAdmits(double w, double threshold) { return w >= threshold; }
} // namespace OpenMS
