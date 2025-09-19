// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qscore.h>
#include <iomanip>

namespace OpenMS
{
//-28.0219, -0.1854, -0.1084, 0.0312, 0, 25.2087
  std::vector<double> Qscore::weight_centroid_ { -21.0476, 1.5045, -0.1303, 0.183, 0.1834, 17.804};
  // Att0                21.0476
  // Att1                -1.5045
  // Att2                 0.1303
  // Att3                 -0.183
  // Att4                -0.1834
  // Intercept           -17.804
  std::unordered_map<int, std::vector<double>> weight_idscore_ {  
    // {8,  { 26.1052, -19.5365, 9.4123, -17.0423 }},
    // {9,  { 11.2583,  0.5577, -3.1241, -9.1871 }},
    // {10, { 11.9844, -1.5740, -0.5244, -9.6836 }},
    // {11, { 10.4434, -1.6515, -0.5515, -8.2053 }},
    // {12, { 6.9904,  -0.2193, -1.0180, -5.5812 }},
    // {13, { 7.3458,  -0.2778, -0.3595, -6.1223 }},
    // {14, { 4.9498,   0.0001, -0.1290, -4.2176 }},
    // {15, { 5.6324,   0.7917,  0.0000, -5.3158 }},
    // {16, { 6.0133,   0.7975,  0.4735, -5.8266 }},
    // {17, { 5.4421,   1.4079,  0.6455, -5.7315 }},
    // {18, { 5.6485,   1.9504,  0.3745, -6.0669 }},
    // {19, { 4.9349,   0.8475,  1.5304, -5.2819 }},
    // {20, { 4.2456,   2.3796,  0.8399, -5.2344 }},
    {21, { 4.4034,   2.2993,  1.1397, -5.3708  }},
    {22, { 5.6630,   1.0324,  2.6398, -6.3672  }},
    {23, { 5.2754,   2.6164,  1.4306, -6.3910  }},
    {24, { 7.4469,   1.3985,  2.0493, -7.7982  }},
    {25, { 8.3606,   0.5524,  3.0693, -8.5154  }},
    {26, { 7.0476,   1.6421,  2.4276, -7.6368  }},
    {27, { 8.4782,   0.5909,  3.4890, -8.6889  }},
    {28, { 8.4931,   2.0384,  1.4288, -8.5946  }},
    {29, { 8.7550,   1.2820,  2.5195, -8.8471  }},
    {30, { 9.6638,   0.5185,  3.1293, -9.4531  }},
    {31, { 9.9730,   1.1585,  3.5695, -10.0997 }},
    {32, { 9.5586,   1.2271,  3.1828, -9.6668  }},
    {33, { 9.7409,   0.1352,  5.0665, -9.9739  }},
    {34, { 10.4065,  0.9047,  3.9907, -10.5302 }},
    {35, { 10.6173,  0.7364,  4.6384, -10.8352 }},
    {36, { 11.0791,  0.2346,  5.0102, -11.1553 }},
    {37, { 11.0333,  0.9526,  3.7254, -11.1153 }},
    {38, { 11.4484,  1.1964,  2.8130, -11.3256 }},
    {39, { 11.6916,  1.6224,  4.0290, -12.2005 }},
    {40, { 10.3172,  0.0000,  5.9899, -10.8953 }},
    {41, { 10.3514,  0.4176,  5.5409, -11.1420 }},
    {42, { 11.5446,  1.2718,  4.6245, -12.2802 }},
    {43, { 9.6587,   0.4470,  6.4660, -10.9487 }},
    {44, { 10.3758,  0.0000,  6.3461, -11.4003 }},
    {45, { 10.3325,  0.0000,  8.6158, -12.1381 }},
    {46, { 10.4593,  0.5347,  6.2897, -11.7876 }},
    {47, { 10.6613,  0.1740,  8.0770, -12.4947 }},
    {48, { 9.1491,   2.2513,  5.2683, -11.2414 }},
    {49, { 10.0220,  0.5149,  7.0649, -11.8878 }},
    {50, { 9.7405,   0.9225,  6.8569, -11.7534 }},
    {51, { 9.6670,   0.3868,  6.2585, -11.2815 }},
    {52, { 8.0511,  -0.1232,  6.4464, -9.6824  }},
    {53, { 8.0043,  -1.0276,  6.8444, -9.3480  }},
    {54, { 9.2859,   0.3974,  7.4335, -11.4685 }},
    {55, { 8.7304,   0.7053,  6.8483, -10.9784 }}
  };
  // Qscore, Monoisotopic Mass, Charge, Intercept
  
  std::unordered_map<int, std::unordered_map<int, float>> Qscore::getIDscores(const PeakGroup* pg) {
    std::unordered_map<int, std::unordered_map<int, float>> idscores;

    auto qscores = pg->getAllQscores();
    float mass = pg->getMonoMass();

    for (const auto& [charge, qscore] : qscores) {
      for (const auto& [hcd, weights] : weight_idscore_) {
        auto affine = weights[0] * qscore + weights[1] * (mass / 20000) + weights[2] * (charge / 30) + weights[3];
        auto idscore = 1.0 / (1.0 + exp(-affine));
        idscores[charge][hcd] = idscore;
      }        
    }
    return idscores;
  }


  /// calculate Qscore using PeakGroup attributes
  std::unordered_map<int, float> Qscore::getQscores(const PeakGroup* pg)
  {
    std::unordered_map<int, float> qscores;
    if (pg->empty())
    { // all zero
      return qscores;
    }

    auto [min_charge, max_charge] = pg->getAbsChargeRange();

    for (int c = min_charge; c <= max_charge; c++) {
          if (pg->getChargeIntensity(c) <= 0) {
            // Skip computation if charge is not present
            continue;
          }

          double score = weight_centroid_.back() + .5;
          auto fv = toFeatureVector_(pg, c);

          for (Size i = 0; i < weight_centroid_.size() - 1; i++)
          {
            score += fv[i] * weight_centroid_[i];
          }
          double qscore = 1.0 / (1.0 + exp(score));
          qscores[c] = qscore;
    }

    return qscores;

  }

  /// convert PeakGroup into feature (attribute) vector
  std::vector<double> Qscore::toFeatureVector_(const PeakGroup* pg, int charge)
  {
    std::vector<double> fvector(5, .0); // length of weights vector - 1, excluding the intercept weight.
    if (pg->empty())
      return fvector;
    int index = 0;
    fvector[index++] = pg->getIsotopeCosine(); // (log2(a + d));

    fvector[index++] = pg->getIsotopeCosine() - pg->getChargeIsotopeCosine(charge); // (log2(d + a / (d + a)));

    fvector[index++] = log2(1 + pg->getChargeSNR(charge)); //(log2(d + a / (d + a)));

    fvector[index++] = log2(1 + pg->getChargeSNR(charge)) - log2(1 + pg->getSNR()); //(log2(a + d));

    fvector[index++] = pg->getAvgPPMError();

    return fvector;
  }

  /// to write down training csv file header.
  void Qscore::writeAttCsvForQscoreTrainingHeader(std::fstream& f)
  {
    PeakGroup pg;
    Size att_count = toFeatureVector_(&pg, 1).size();
    for (Size i = 0; i < att_count; i++)
      f << "Att" << i << ",";
    f << "Class\n";
  }

  /// to write down training csv file rows
  void Qscore::writeAttCsvForQscoreTraining(const DeconvolvedSpectrum& deconvolved_spectrum, std::fstream& f)
  {
    DeconvolvedSpectrum dspec;
    dspec.reserve(deconvolved_spectrum.size());
    for (auto& pg : deconvolved_spectrum)
    {
      dspec.push_back(pg);
    }

    if (dspec.empty())
      return;

    for (auto& pg : dspec)
    {
      bool target = pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::target;
      auto fv = toFeatureVector_(&pg, pg.getRepAbsCharge());

      for (auto& item : fv)
      {
        f << item << ",";
      }
      f << (target ? "T" : "F") << "\n";
    }
  }
} // namespace OpenMS