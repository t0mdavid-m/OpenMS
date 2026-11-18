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
    @brief The witness probability W of an MS3 (ADR-0076): the chance that the MS3 shows at least one bond of its
    target fragment, from the ms3_esm package's gradient-boosted trees (w_ext_v3, XGBoost binary:logistic, 36 inputs).

    The trees are read at runtime from share/OpenMS/FLASHIDA/ms3_w_gate.json, a byte copy of the package's
    w_ext_v3.json. The evaluation follows the package's ENGINE.md exactly: every input and every split threshold is a
    32-bit float, a node goes left iff input < threshold (strictly) and a missing input goes where the node's
    default_left says, and W = 1 / (1 + exp(-(ln(b / (1 - b)) + sum of the leaf values))). The thresholds sit at values
    of the training data, so an input is often exactly equal to one: compared in any other precision, such an MS3 goes
    the wrong way (1,470 of the 18,266 recorded MS3 would, reference_rows_ties.tsv).
  */

  /// The 35 raw values of one candidate MS3, in the order of the package's w_gate.RAW. Slot 0 is frag_N (1 for an a, b
  /// or c ion, 0 for x, y or z) rather than the package's frag_type letter. NaN = missing.
  enum WitnessRawIndex : int
  {
    kWFragN = 0,
    kWZMs2,
    kWZMs3,
    kWCe1,
    kWCid3,
    kWIntensity,
    kWL,
    kWIdScore,
    kWIsoCos,
    kWQscore,
    kWSnr,
    kWHasShift,
    kWNMs2Prior,
    kWNMs3Prior,
    kWMassMs2,
    kWMassMs3,
    kWPrecIntensity,
    kWPrecQscore,
    kWPrecIsoCos,
    kWPrecSnr,
    kWMs1Plan,
    kWChargeSnr,
    kWChargeCos,
    kWChargeScore,
    kWPpmError,
    kWMs2MassCount,
    kWMs2TicCoverage,
    kWMs2PrecursorShare,
    kWMs2FragmentCount,
    kWMs2TagCount,
    kWMs2Ions,
    kWQueueAhead,
    kWRank,
    kWPlanAge,
    kWRtFrac,
    kWitnessRawCount   ///< 35
  };

  /// The package's name of each raw value (inputs.tsv's raw_value_to_give; frag_N for slot 0), for the [MS3-W] line.
  OPENMS_DLLAPI const char* witnessRawName(int index);

  /// The number of model inputs (columns of the trees): the 35 raw values with z_ms3_per_res derived.
  constexpr int kWitnessInputCount = 36;

  struct OPENMS_DLLAPI WitnessRaw
  {
    std::array<double, kWitnessRawCount> v;   ///< NaN-filled by the constructor
    WitnessRaw();
    double& operator[](int i) { return v[static_cast<size_t>(i)]; }
    double operator[](int i) const { return v[static_cast<size_t>(i)]; }
    /// The first raw value that is missing (NaN), or -1 when all 35 are present.
    int firstMissing() const;
  };

  /// The 36 model inputs of @p raw, as ENGINE.md section 3 builds them: frag_N; z_ms3_per_res = z_ms3 / L;
  /// log10_intensity = log10(max(intensity, 1)); id_score, x_prec_intensity, x_ms1_plan and x_ms2_mass_count rounded to
  /// 32 bits first, then log10(max(value, 1)); every value rounded to 32 bits. A missing value stays NaN.
  OPENMS_DLLAPI std::array<float, kWitnessInputCount> witnessModelInputs(const WitnessRaw& raw);

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
    /// The shipped model: File::find("FLASHIDA/ms3_w_gate.json") through OPENMS_DATA_PATH. Throws
    /// std::invalid_argument naming the path and the failed check when it cannot be found, read or accepted.
    static WitnessModel loadShipped();

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

    /// The predicted W for the 36 model inputs @p x.
    double w(const std::array<float, kWitnessInputCount>& x) const;

    /// W of a candidate's raw values: w(witnessModelInputs(raw)).
    double w(const WitnessRaw& raw) const { return w(witnessModelInputs(raw)); }

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
