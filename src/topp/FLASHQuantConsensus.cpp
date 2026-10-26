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
#include <array>
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
  its score is the cosine of its (charge, isotope) slot areas to the feature's template (its slot pattern where detected), if its matched
  traces pass three gates set from the detected features: |m/z error| at most twice their median, apex offset at most their 90%
  quantile, and co-elution (each trace's elution profile against the sum of the others) >= 0.5; else 0. The null: a decoy mass below the
  envelope by its width + 5 Da (a shift within the envelope, as in bottom-up, falls on the species' own isotopes), measured, gated and
  scored the same way. The FDR of a score threshold is the number of null scores above it over the number of target
  scores above it; values whose q-value exceeds -mbr_fdr are NaN.
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
                        "status (detected, mbr, or nan) and, for a match between runs, the cosine of its slot pattern to the template and its q-value.", true);
    setValidFormats_("out", ListUtils::create<String>("tsv"));
    registerDoubleOption_("min_run_fraction", "<fraction>", 0.5, "Consensus features detected in at least this fraction of the runs are quantified in every run.", false);
    setMinFloat_("min_run_fraction", 0.0);
    setMaxFloat_("min_run_fraction", 1.0);
    registerOutputFile_("out_mbr_detail", "<file>", "", "Match between runs, per transfer: the consensus feature, run, value, the cosine to the template, "
                        "its decoy value and cosine, and the median abundance of the feature where detected.", false, true);
    setValidFormats_("out_mbr_detail", ListUtils::create<String>("tsv"));
    registerDoubleOption_("mbr_fdr", "<fdr>", 0.01, "Match between runs: values are kept up to this FDR (target-decoy, see above); 1 keeps every "
                          "value, a negative value switches match between runs off (NaN where not detected).", false);
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
    const double min_fraction = getDoubleOption_("min_run_fraction"), mbr_fdr = getDoubleOption_("mbr_fdr");
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
    // the decoy mass of a feature: below its isotope envelope (core isotopes, >= 10% of the most abundant) by the envelope's width + 5 Da,
    // so that no isotope of the species at any charge falls on it
    std::vector<double> decoy_shift(kept.size());
    for (Size k = 0; k < kept.size(); ++k)
    {
      const auto iso = avg.get(summary[k].mass);
      double iso_max = 0;
      for (Size i = 0; i < iso.size(); ++i) iso_max = std::max(iso_max, (double)iso[i].getIntensity());
      int first = -1, last = -1;
      for (Size i = 0; i < iso.size(); ++i)
      {
        if (iso[i].getIntensity() < 0.1 * iso_max) continue;
        if (first < 0) first = (int)i;
        last = (int)i;
      }
      decoy_shift[k] = (last - first + 1) * Constants::ISOTOPE_MASSDIFF_55K_U + 5.0;
    }

    // every run: the features and, where not detected, their two decoys
    std::vector<std::vector<FLASHQuantAlgorithm::MassAbundance>> value(kept.size(), std::vector<FLASHQuantAlgorithm::MassAbundance>(run_count));
    std::vector<std::vector<double>> ratio(kept.size(), std::vector<double>(run_count, -1));
    std::vector<std::vector<std::array<double, 2>>> null_cos(kept.size(), std::vector<std::array<double, 2>>(run_count, std::array<double, 2> {0, 0}));
    // per feature and run: the (charge, isotope) slot areas of the target and, where not detected, of its two decoys
    using Slots = std::vector<std::tuple<int, int, double>>;
    std::vector<std::vector<std::array<Slots, 3>>> slot_patterns(kept.size(), std::vector<std::array<Slots, 3>>(run_count)); // (not "slots": a Qt macro)
    std::vector<std::vector<std::array<double, 3>>> measured(kept.size(), std::vector<std::array<double, 3>>(run_count, std::array<double, 3> {0, 0, 0}));
    // per feature, run and kind (target, decoy): apex offset and apex SD (in half windows), |m/z error| (ppm) and co-elution of the matched traces
    std::vector<std::vector<std::array<std::array<double, 4>, 3>>> evidence(kept.size(), std::vector<std::array<std::array<double, 4>, 3>>(run_count));
    for (Size r = 0; r < run_count; ++r)
    {
      MzMLFile mzml;
      mzml.setLogType(log_type_);
      mzml.getOptions().setMSLevels({1});
      PeakMap map;
      mzml.load(in[r], map);
      map.sortSpectra(true);
      std::vector<FLASHQuantAlgorithm::QuantTarget> targets;
      std::vector<std::pair<Size, int>> owner; // (consensus, 0 target / 1 mass decoy / 2 time decoy / 3, 4 decoys of the mass decoy)
      for (Size k = 0; k < kept.size(); ++k)
      {
        const auto& s = summary[k];
        const double center = (s.aligned + (anchors[r].empty() ? 0.0 : offset_(anchors[r], s.aligned))) * 60.0, half = s.half * 60.0;
        targets.push_back({s.mass, center - half, center + half, s.zmin, s.zmax});
        owner.emplace_back(k, 0);
        if (kept[k].runs.count(r) == 0 && mbr_fdr >= 0)
        {
          targets.push_back({s.mass - decoy_shift[k], center - half, center + half, s.zmin, s.zmax});
          owner.emplace_back(k, 1);
        }
      }
      FLASHQuantAlgorithm fq;
      std::vector<FLASHQuantAlgorithm::MassAbundance> abundances;
      std::vector<std::vector<std::tuple<int, int, double>>> slot_areas;
      std::vector<FLASHQuantAlgorithm::TargetEvidence> target_evidence;
      fq.quantifyTargets(map, targets, avg, abundances, &slot_areas, &target_evidence);
      for (Size t = 0; t < targets.size(); ++t)
      {
        const auto [k, kind] = owner[t];
        if (kind == 0) value[k][r] = abundances[t];
        slot_patterns[k][r][kind] = std::move(slot_areas[t]);
        measured[k][r][kind] = abundances[t].first;
        const double center = (targets[t].rt_start + targets[t].rt_end) / 2, half = (targets[t].rt_end - targets[t].rt_start) / 2;
        evidence[k][r][kind] = {std::abs(target_evidence[t].apex - center) / half, target_evidence[t].apex_sd / half, target_evidence[t].ppm,
                                target_evidence[t].coelution};
      }
      OPENMS_LOG_INFO << "Run " << (r + 1) << " of " << run_count << ": " << targets.size() << " targets measured" << endl;
    }

    // The template of a feature: its (charge, isotope) slot pattern where detected (the mean of the unit-length slot vectors). A transfer
    // scores the cosine of its slot pattern to the template; the null: the pattern of its decoy (below the envelope) against the same
    // template. Random matches do not follow the charge and isotope pattern of the species.
    std::vector<double> target_scores, null_scores;
    // Gates on the evidence of the matched traces, from the features where detected (measured the same way): |m/z error| at most twice
    // the median, apex offset at most the 90% quantile, and co-elution (the isotopes and charges aligned in retention time) >= 0.5. A
    // decoy passes the same gates; the q-values are over the cosines of those passing (PXD042169: at 1% FDR, ~1.2-1.5k transfers in
    // SpikeIn and ~0.9k in ProteomeMix instead of none on the cosine alone).
    double ppm_gate = 0, offset_gate = 0;
    {
      std::vector<double> ppms, offsets;
      for (Size k = 0; k < kept.size(); ++k)
      {
        for (Size r = 0; r < run_count; ++r)
        {
          if (kept[k].runs.count(r) == 0 || measured[k][r][0] <= 0) continue;
          offsets.push_back(evidence[k][r][0][0]);
          ppms.push_back(evidence[k][r][0][2]);
        }
      }
      ppm_gate = 2 * median_(ppms);
      if (!offsets.empty())
      {
        std::sort(offsets.begin(), offsets.end());
        offset_gate = offsets[(Size)(0.9 * (offsets.size() - 1))];
      }
      OPENMS_LOG_INFO << "Match between runs: gates from the detected features: |m/z error| <= " << ppm_gate << " ppm, apex offset <= " << offset_gate
                      << " half windows, co-elution >= 0.5" << endl;
    }
    const auto gated = [&](Size k, Size r, int kind) {
      const auto& e = evidence[k][r][kind];
      return measured[k][r][kind] > 0 && e[2] <= ppm_gate && e[0] <= offset_gate && e[3] >= 0.5;
    };
    for (Size k = 0; k < kept.size() && mbr_fdr >= 0; ++k)
    {
      std::map<std::pair<int, int>, double> tmpl;
      for (Size r = 0; r < run_count; ++r)
      {
        if (kept[k].runs.count(r) == 0) continue;
        double norm = 0;
        for (const auto& [z, i, a] : slot_patterns[k][r][0]) norm += a * a;
        if (norm <= 0) continue;
        for (const auto& [z, i, a] : slot_patterns[k][r][0]) tmpl[{z, i}] += a / std::sqrt(norm);
      }
      double tmpl_norm = 0;
      for (const auto& [key, a] : tmpl) tmpl_norm += a * a;
      const auto cosine_to = [&](const Slots& v, const std::map<std::pair<int, int>, double>& with) {
        double dot = 0, norm = 0;
        for (const auto& [z, i, a] : v)
        {
          norm += a * a;
          const auto it = with.find({z, i});
          if (it != with.end()) dot += a * it->second;
        }
        return (norm > 0 && tmpl_norm > 0) ? dot / std::sqrt(norm * tmpl_norm) : 0.0;
      };
      const auto cosine = [&](const Slots& v) { return cosine_to(v, tmpl); };
      for (Size r = 0; r < run_count; ++r)
      {
        if (kept[k].runs.count(r) > 0) continue;
        ratio[k][r] = cosine(slot_patterns[k][r][0]);
        null_cos[k][r] = {cosine(slot_patterns[k][r][1]), 0.0};
        // the score: the cosine of a transfer that passes the gates, else 0 (never kept)
        if (!gated(k, r, 0)) ratio[k][r] = 0;
        if (!gated(k, r, 1)) null_cos[k][r][0] = 0;
        target_scores.push_back(ratio[k][r]);
        null_scores.push_back(null_cos[k][r][0]);
      }
    }
    // q-values: the FDR of a threshold T is (null scores >= T) / (target scores >= T) (one null per transfer); the q-value of a
    // score is the lowest FDR of the thresholds at or below it
    std::sort(target_scores.begin(), target_scores.end());
    std::sort(null_scores.begin(), null_scores.end());
    std::map<double, double> qvalue; // distinct target score -> q-value
    {
      double lowest = std::numeric_limits<double>::max();
      for (Size i = 0; i < target_scores.size(); ++i)
      {
        if (i > 0 && target_scores[i] == target_scores[i - 1]) continue;
        const double t = target_scores[i];
        const double targets_above = (double)(target_scores.size() - i);
        const double nulls_above = (double)(null_scores.end() - std::lower_bound(null_scores.begin(), null_scores.end(), t));
        lowest = std::min(lowest, std::min(1.0, nulls_above / targets_above));
        qvalue[t] = lowest;
      }
    }
    const auto q_of = [&qvalue](double t) { const auto it = qvalue.find(t); return it == qvalue.end() ? 1.0 : it->second; };
    double threshold = std::numeric_limits<double>::infinity();
    for (const auto& [t, q] : qvalue)
    {
      if (t > 0 && q <= mbr_fdr) { threshold = t; break; }
    }
    if (mbr_fdr >= 0)
    {
      const Size null_above = null_scores.end() - std::lower_bound(null_scores.begin(), null_scores.end(), threshold);
      const Size target_above = target_scores.end() - std::lower_bound(target_scores.begin(), target_scores.end(), threshold);
      OPENMS_LOG_INFO << "Match between runs: " << target_scores.size() << " transfers; score threshold at FDR " << mbr_fdr << ": " << threshold
                      << " (targets at or above it " << target_above << ", nulls " << null_above << ")" << endl;
    }

    const String detail_out = getStringOption_("out_mbr_detail");
    if (!detail_out.empty())
    {
      std::ofstream ds(detail_out);
      ds << "ConsensusIndex\tRun\tValue\tDecoyValue\tCosine\tDecoyCosine\tDetectedMedian\tApexOffset\tApexSD\tPPM\tCoelution\tDecoyApexOffset\tDecoyApexSD\tDecoyPPM\tDecoyCoelution\n";
      for (Size k = 0; k < kept.size(); ++k)
      {
        std::vector<double> det;
        for (Size r = 0; r < run_count; ++r) { if (kept[k].runs.count(r) > 0) det.push_back(value[k][r].first); }
        for (Size r = 0; r < run_count; ++r)
        {
          if (kept[k].runs.count(r) > 0) continue;
          ds << k << "\t" << r << "\t" << std::to_string(measured[k][r][0]) << "\t" << std::to_string(measured[k][r][1])
             << "\t" << std::to_string(ratio[k][r]) << "\t" << std::to_string(null_cos[k][r][0])
             << "\t" << std::to_string(median_(det));
          for (int kind = 0; kind < 2; ++kind) for (double e : evidence[k][r][kind]) ds << "\t" << std::to_string(e);
          ds << "\n";
        }
      }
    }

    // output
    std::ofstream os(out);
    os << "ConsensusIndex\tMonoisotopicMass\tApexRetentionTime\tMinCharge\tMaxCharge\tDetectedRuns";
    std::vector<String> names;
    for (const auto& f : in) names.push_back(QFileInfo(f.toQString()).completeBaseName());
    for (const auto& n : names) os << "\tAbundance_" << n << "\tStatus_" << n << "\tMatchCosine_" << n << "\tMatchQvalue_" << n;
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
          os << "\t" << std::to_string(v) << "\tdetected\t\t";
          continue;
        }
        mbr_total++;
        const double q = (mbr_fdr >= 0 && v > 0) ? q_of(ratio[k][r]) : 1.0;
        const bool ok = mbr_fdr >= 0 && v > 0 && q <= mbr_fdr;
        mbr_kept += ok;
        os << "\t" << (ok ? std::to_string(v) : String("nan")) << "\t" << (ok ? "mbr" : "nan") << "\t" << (mbr_fdr >= 0 ? std::to_string(ratio[k][r]) : String(""))
           << "\t" << (mbr_fdr >= 0 ? std::to_string(q) : String(""));
      }
      os << "\n";
    }
    OPENMS_LOG_INFO << "Match between runs: " << mbr_kept << " of " << mbr_total << " values kept (q-value <= " << mbr_fdr << ")" << endl;
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPFLASHQuantConsensus tool;
  return tool.main(argc, argv);
}

/// @endcond
