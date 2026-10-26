// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHQuantAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <QFileInfo>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------
/**
  @page TOPP_FLASHQuantConsensus TOPP_FLASHQuantConsensus

  @brief Quantifies the MS1 features of several runs together: FLASHDeconv features matched across the runs (retention times aligned),
  and every matched feature measured the same way in every run from its raw signal, with match between runs.

  Input: the MS1 mzML files and the FLASHDeconv feature tsv files (-out) of the runs. The retention times of each run are aligned to the
  first run by the features most runs share; features of the same mass (10 ppm, or one isotope apart) and aligned apex (1 min) form a
  consensus feature. A consensus feature detected in at least -min_run_fraction of the runs is measured in every run: the raw mass
  traces at its charges (the median range of its detections) and core isotopes within 10 ppm whose apex lies in its window (the aligned
  apex mapped into the run, +- half its median width), their bulk areas summed. In a run where it was not detected (match between runs),
  two decoys are measured the same way (mass + 5.5 Da, between the isotopes of the true mass; window 10 min later), and the value is kept
  when it exceeds -mbr_min_ratio times the larger decoy (else NaN).
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFLASHQuantConsensus : public TOPPBase
{
public:
  TOPPFLASHQuantConsensus() : TOPPBase("FLASHQuantConsensus", "Quantification of FLASHDeconv features across runs, with match between runs.", false)
  {
  }

protected:
  struct Feature
  {
    double mass, apex, start, end, aligned; // minutes
    int min_charge, max_charge;
  };

  struct Consensus
  {
    double mass, aligned;
    std::map<Size, Size> runs; // run -> feature
  };

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", {}, "MS1 mzML files of the runs (as given to FLASHDeconv), centroided.", true);
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFileList_("in_feature", "<files>", {}, "FLASHDeconv feature tsv files (-out) of the runs, in the order of -in.", true);
    setValidFormats_("in_feature", ListUtils::create<String>("tsv"));
    registerOutputFile_("out", "<file>", "", "Consensus feature tsv: per feature its mass, aligned apex and charges, and per run its abundance, "
                        "status (detected, mbr, or nan) and, for a match between runs, the ratio of the value to the larger decoy.", true);
    setValidFormats_("out", ListUtils::create<String>("tsv"));
    registerDoubleOption_("min_run_fraction", "<fraction>", 0.5, "Consensus features detected in at least this fraction of the runs are quantified in every run.", false);
    setMinFloat_("min_run_fraction", 0.0);
    setMaxFloat_("min_run_fraction", 1.0);
    registerDoubleOption_("mbr_min_ratio", "<ratio>", 3.0, "Match between runs: a value is kept when it exceeds this many times the larger of its two "
                          "decoys; 0 keeps every value, a negative value switches match between runs off (NaN where not detected).", false);
  }

  /// the FLASHDeconv MS1 features of a feature tsv
  std::vector<Feature> readFeatures_(const String& file) const
  {
    std::ifstream in(file);
    std::string line;
    std::getline(in, line);
    std::vector<String> header;
    String(line).split('\t', header);
    std::map<String, Size> col;
    for (Size i = 0; i < header.size(); ++i) col[header[i]] = i;
    for (const String& c : {"MSLevel", "MonoisotopicMass", "StartRetentionTime", "EndRetentionTime", "ApexRetentionTime", "MinCharge", "MaxCharge"})
    {
      if (col.find(c) == col.end()) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file, String("no column ") + c);
    }
    std::vector<Feature> features;
    while (std::getline(in, line))
    {
      std::vector<String> v;
      String(line).split('\t', v);
      if (v.size() < header.size() || v[col["MSLevel"]] != "1") continue;
      features.push_back({v[col["MonoisotopicMass"]].toDouble(), v[col["ApexRetentionTime"]].toDouble(), v[col["StartRetentionTime"]].toDouble(),
                          v[col["EndRetentionTime"]].toDouble(), v[col["ApexRetentionTime"]].toDouble(), v[col["MinCharge"]].toInt(),
                          v[col["MaxCharge"]].toInt()});
    }
    return features;
  }

  static bool sameMass_(double m1, double m2)
  {
    const double iso = Constants::ISOTOPE_MASSDIFF_55K_U, tol = 10e-6 * m1;
    return std::abs(m1 - m2) <= tol || std::abs(m1 - m2 - iso) <= tol || std::abs(m1 - m2 + iso) <= tol;
  }

  /// features of the same mass and retention time (by @p rt) across runs, greedily in mass order
  static std::vector<Consensus> match_(const std::vector<std::vector<Feature>>& runs, bool aligned, double rt_tol)
  {
    std::vector<std::tuple<double, Size, Size>> all;
    for (Size r = 0; r < runs.size(); ++r)
    {
      for (Size i = 0; i < runs[r].size(); ++i) all.emplace_back(runs[r][i].mass, r, i);
    }
    std::sort(all.begin(), all.end());
    std::vector<Consensus> cons;
    for (const auto& [m, r, i] : all)
    {
      const double rt = aligned ? runs[r][i].aligned : runs[r][i].apex;
      bool placed = false;
      for (Size k = cons.size() > 300 ? cons.size() - 300 : 0; k < cons.size(); ++k)
      {
        auto& c = cons[k];
        if (c.runs.count(r) == 0 && sameMass_(c.mass, m) && std::abs(c.aligned - rt) <= rt_tol)
        {
          c.runs[r] = i;
          placed = true;
          break;
        }
      }
      if (!placed) cons.push_back({m, rt, {{r, i}}});
    }
    return cons;
  }

  static double median_(std::vector<double> v)
  {
    if (v.empty()) return 0;
    std::nth_element(v.begin(), v.begin() + v.size() / 2, v.end());
    return v[v.size() / 2];
  }

  /// the retention time offset of a run to the first one at retention time @p rt: the median offset of the anchor features within 2 min
  static double offset_(const std::vector<std::pair<double, double>>& anchors, double rt)
  {
    const auto lo = std::lower_bound(anchors.begin(), anchors.end(), std::make_pair(rt - 2.0, -1e300));
    const auto hi = std::upper_bound(anchors.begin(), anchors.end(), std::make_pair(rt + 2.0, 1e300));
    std::vector<double> v;
    for (auto it = lo; it != hi; ++it) v.push_back(it->second);
    if (v.empty())
    {
      for (const auto& a : anchors) v.push_back(a.second);
    }
    return median_(v);
  }

  ExitCodes main_(int, const char**) override
  {
    const StringList in = getStringList_("in"), in_feature = getStringList_("in_feature");
    const String out = getStringOption_("out");
    const double min_fraction = getDoubleOption_("min_run_fraction"), mbr_min_ratio = getDoubleOption_("mbr_min_ratio");
    if (in.size() != in_feature.size() || in.empty())
    {
      OPENMS_LOG_ERROR << "-in and -in_feature need the same number of files" << endl;
      return ILLEGAL_PARAMETERS;
    }
    const Size run_count = in.size();
    std::vector<std::vector<Feature>> runs(run_count);
    for (Size r = 0; r < run_count; ++r) runs[r] = readFeatures_(in_feature[r]);

    // retention time alignment to the first run: anchors are the features matched (3 min, unaligned) in >= 80% of the runs and the first
    std::vector<std::vector<std::pair<double, double>>> anchors(run_count); // per run: (its apex, offset to the first run), by apex
    for (const auto& c : match_(runs, false, 3.0))
    {
      if (c.runs.size() < 0.8 * run_count || c.runs.count(0) == 0) continue;
      const double ref = runs[0][c.runs.at(0)].apex;
      for (const auto& [r, i] : c.runs) anchors[r].emplace_back(runs[r][i].apex, runs[r][i].apex - ref);
    }
    for (Size r = 0; r < run_count; ++r)
    {
      std::sort(anchors[r].begin(), anchors[r].end());
      for (auto& f : runs[r]) f.aligned = anchors[r].empty() ? f.apex : f.apex - offset_(anchors[r], f.apex);
    }
    std::vector<Consensus> cons = match_(runs, true, 1.0);
    std::vector<Consensus> kept;
    for (auto& c : cons)
    {
      if (c.runs.size() >= min_fraction * run_count && !c.runs.empty()) kept.push_back(std::move(c));
    }
    OPENMS_LOG_INFO << "Consensus features: " << cons.size() << ", detected in at least " << min_fraction * 100 << "% of the " << run_count
                    << " runs: " << kept.size() << endl;

    // per consensus feature: median mass, aligned apex, half width and charge range of its detections
    struct Summary { double mass, aligned, half; int zmin, zmax; };
    std::vector<Summary> summary;
    double max_mass = 0;
    for (const auto& c : kept)
    {
      std::vector<double> m, a, h, z1, z2;
      for (const auto& [r, i] : c.runs)
      {
        const auto& f = runs[r][i];
        m.push_back(f.mass);
        a.push_back(f.aligned);
        h.push_back((f.end - f.start) / 2);
        z1.push_back(f.min_charge);
        z2.push_back(f.max_charge);
      }
      summary.push_back({median_(m), median_(a), std::max(0.5, median_(h)), (int)median_(z1), (int)median_(z2)});
      max_mass = std::max(max_mass, summary.back().mass);
    }
    SpectralDeconvolution sd;
    Param sd_param = sd.getDefaults();
    sd_param.setValue("max_mass", std::max(1000.0, max_mass * 1.1 + 100));
    sd.setParameters(sd_param);
    sd.calculateAveragine(false);
    const auto& avg = sd.getAveragine();

    // every run: the features and, where not detected, their two decoys
    std::vector<std::vector<FLASHQuantAlgorithm::MassAbundance>> value(kept.size(), std::vector<FLASHQuantAlgorithm::MassAbundance>(run_count));
    std::vector<std::vector<double>> ratio(kept.size(), std::vector<double>(run_count, -1));
    for (Size r = 0; r < run_count; ++r)
    {
      MzMLFile mzml;
      mzml.setLogType(log_type_);
      mzml.getOptions().setMSLevels({1});
      PeakMap map;
      mzml.load(in[r], map);
      map.sortSpectra(true);
      std::vector<FLASHQuantAlgorithm::QuantTarget> targets;
      std::vector<std::pair<Size, int>> owner; // (consensus, 0 target / 1 mass decoy / 2 time decoy)
      for (Size k = 0; k < kept.size(); ++k)
      {
        const auto& s = summary[k];
        const double center = (s.aligned + (anchors[r].empty() ? 0.0 : offset_(anchors[r], s.aligned))) * 60.0, half = s.half * 60.0;
        targets.push_back({s.mass, center - half, center + half, s.zmin, s.zmax});
        owner.emplace_back(k, 0);
        if (kept[k].runs.count(r) > 0 || mbr_min_ratio < 0) continue;
        targets.push_back({s.mass + 5.5, center - half, center + half, s.zmin, s.zmax});
        owner.emplace_back(k, 1);
        targets.push_back({s.mass, center + 600 - half, center + 600 + half, s.zmin, s.zmax});
        owner.emplace_back(k, 2);
      }
      FLASHQuantAlgorithm fq;
      std::vector<FLASHQuantAlgorithm::MassAbundance> abundances;
      fq.quantifyTargets(map, targets, avg, abundances);
      std::vector<double> decoy(kept.size(), 0);
      for (Size t = 0; t < targets.size(); ++t)
      {
        const auto [k, kind] = owner[t];
        if (kind == 0) value[k][r] = abundances[t];
        else decoy[k] = std::max(decoy[k], abundances[t].first);
      }
      for (Size k = 0; k < kept.size(); ++k)
      {
        if (kept[k].runs.count(r) > 0) continue;
        ratio[k][r] = decoy[k] > 0 ? value[k][r].first / decoy[k] : (value[k][r].first > 0 ? std::numeric_limits<double>::infinity() : 0.0);
      }
      OPENMS_LOG_INFO << "Run " << (r + 1) << " of " << run_count << ": " << targets.size() << " targets measured" << endl;
    }

    // output
    std::ofstream os(out);
    os << "ConsensusIndex\tMonoisotopicMass\tApexRetentionTime\tMinCharge\tMaxCharge\tDetectedRuns";
    std::vector<String> names;
    for (const auto& f : in) names.push_back(QFileInfo(f.toQString()).completeBaseName());
    for (const auto& n : names) os << "\tAbundance_" << n << "\tStatus_" << n << "\tMatchRatio_" << n;
    os << "\n";
    Size mbr_kept = 0, mbr_total = 0;
    for (Size k = 0; k < kept.size(); ++k)
    {
      const auto& s = summary[k];
      os << k << "\t" << std::to_string(s.mass) << "\t" << s.aligned << "\t" << s.zmin << "\t" << s.zmax << "\t" << kept[k].runs.size();
      for (Size r = 0; r < run_count; ++r)
      {
        const double v = value[k][r].first;
        if (kept[k].runs.count(r) > 0)
        {
          os << "\t" << std::to_string(v) << "\tdetected\t";
          continue;
        }
        mbr_total++;
        const bool ok = mbr_min_ratio >= 0 && v > 0 && ratio[k][r] > mbr_min_ratio;
        mbr_kept += ok;
        os << "\t" << (ok ? std::to_string(v) : String("nan")) << "\t" << (ok ? "mbr" : "nan") << "\t" << (mbr_min_ratio >= 0 ? std::to_string(ratio[k][r]) : String(""));
      }
      os << "\n";
    }
    OPENMS_LOG_INFO << "Match between runs: " << mbr_kept << " of " << mbr_total << " values kept (ratio to the larger decoy > " << mbr_min_ratio << ")" << endl;
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPFLASHQuantConsensus tool;
  return tool.main(argc, argv);
}

/// @endcond
