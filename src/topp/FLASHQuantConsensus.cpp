// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHQuantAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
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
#include <tuple>

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
  FLASHDeconv looks for its mass in every MS1 scan of the window: the peak group of the mass as FLASHDeconv forms and scores one (raw
  peaks at its charges +- 2, isotope offset search, charge fit, isotope cosine, charge SNR, Qscore), kept if its mass stays within the
  tolerance of the target (or one isotope off), and the scans combined as FLASHDeconv combines the scans of a feature (Qscore2D: one
  minus the lowest product of (1 - Qscore) over consecutive scans, so the whole envelope must recur scan after scan). The score of the
  transfer is that Qscore2D if FLASHDeconv's mass of it (intensity-weighted over the scans) lies within three times the median mass
  error of the detected features of the feature's own mass (measured the same way where it was detected), and it has a raw abundance;
  else 0. The null: a decoy mass below the envelope by its
  width + 5 Da (a shift within the envelope, as in bottom-up, falls on the species' own isotopes), measured and scored the same way. The
  FDR of a score threshold is the number of null scores above it over the number of target scores above it; values whose q-value
  exceeds -mbr_fdr are NaN.
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
                        "status (detected, mbr, or nan) and, for a match between runs, its score (FLASHDeconv's Qscore2D of the mass, see above) and q-value.", true);
    setValidFormats_("out", ListUtils::create<String>("tsv"));
    registerDoubleOption_("min_run_fraction", "<fraction>", 0.25, "Consensus features detected in at least this fraction of the runs are quantified in every run "
                          "(where not detected, by match between runs at -mbr_fdr). Lower values quantify more features, but their many weak transfers "
                          "raise the score threshold of all transfers (PXD042169: 0.25 doubled the quantified features over 0.5 at the same accuracy; "
                          "features in >= 2 runs lost transfers of the strong features).", false);
    setMinFloat_("min_run_fraction", 0.0);
    setMaxFloat_("min_run_fraction", 1.0);
    registerOutputFile_("out_mbr_detail", "<file>", "", "Match between runs, per consensus feature and run (detected or transfer): the value, FLASHDeconv's "
                        "evidence for the mass (Qscore2D, best scan Qscore, scans found, mass error), the same for its decoy where not detected, "
                        "the feature's own mass error where detected (the reference of the mass errors) and its median abundance where detected.", false, true);
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

  /// FLASHDeconv's evidence for a mass in a retention time window
  struct FDEvidence
  {
    double qscore2d = 0, max_qscore = 0, ppm = 0; // ppm: intensity-weighted signed error of the scans' masses to the mass (nearest isotope)
    int scans = 0;                                 // scans where the mass was found
  };

  /// FLASHDeconv's evidence for @p target: in every MS1 scan of its window, the peak group of its mass as FLASHDeconv forms and scores
  /// one (raw peaks recruited at its charges +- 2, isotope offset search, charge fit, isotope cosine >= @p min_cos, charge SNR >= @p min_snr,
  /// Qscore), found if its mass stays within @p tol of the target's (or one isotope off); the scans combined as FLASHDeconv combines the
  /// scans of a feature (Qscore2D: one minus the lowest product of (1 - Qscore) over consecutive scans)
  static FDEvidence fdEvidence_(const PeakMap& map, const std::vector<double>& rts, const FLASHQuantAlgorithm::QuantTarget& target,
                                const FLASHHelperClasses::PrecalculatedAveragine& avg, double tol, double min_cos, double min_snr)
  {
    FDEvidence e;
    const double iso = Constants::ISOTOPE_MASSDIFF_55K_U;
    const int z1 = std::max(1, target.min_charge - 2), z2 = target.max_charge + 2;
    double chain = 1, lowest = 1, weights = 0, ppm_sum = 0;
    bool in_chain = false;
    std::vector<FLASHHelperClasses::LogMzPeak> noisy;
    for (Size s = std::lower_bound(rts.begin(), rts.end(), target.rt_start) - rts.begin(); s < map.size() && rts[s] <= target.rt_end; ++s)
    {
      const MSSpectrum& spec = map[s];
      PeakGroup pg(z1, z2, true);
      pg.setIsotopeDaDistance(iso);
      double next = target.mass, prev_cos = 0;
      int offset = 0;
      bool determined = false;
      for (int k = 0; k < 10; ++k) // as SpectralDeconvolution: recruit and score until the isotope offset settles
      {
        pg.recruitAllPeaksInSpectrum(spec, tol, avg, next, noisy);
        if (pg.empty()) break;
        auto [a, b] = pg.getAbsChargeRange();
        offset = pg.updateQscore(noisy, avg, min_cos, tol, (a + b) < 20, {}, false);
        if (prev_cos > pg.getIsotopeCosine() || offset == 0 || k == 9)
        {
          if (pg.getChargeSNR(pg.getRepAbsCharge()) < min_snr) break;
          if (offset != 0) pg.recruitAllPeaksInSpectrum(spec, tol, avg, pg.getMonoMass() + offset * iso, noisy);
          std::tie(a, b) = pg.getAbsChargeRange();
          pg.updateQscore(noisy, avg, min_cos, tol, (a + b) < 20, {}, true);
          determined = true;
          break;
        }
        prev_cos = pg.getIsotopeCosine();
        next = pg.getMonoMass() + offset * iso;
      }
      double error = 1e300;
      for (int i = -1; i <= 1 && determined; ++i)
      {
        const double d = (pg.getMonoMass() - target.mass - i * iso) / target.mass * 1e6;
        if (std::abs(d) < std::abs(error)) error = d;
      }
      const double q = (determined && std::abs(error) <= tol * 1e6) ? pg.getQscore() : 0.0;
      if (q <= 0)
      {
        in_chain = false;
        continue;
      }
      chain = in_chain ? chain * (1 - q) : 1 - q;
      in_chain = true;
      lowest = std::min(lowest, chain);
      e.scans++;
      e.max_qscore = std::max(e.max_qscore, q);
      weights += pg.getIntensity();
      ppm_sum += pg.getIntensity() * error;
    }
    e.qscore2d = 1 - lowest;
    e.ppm = weights > 0 ? ppm_sum / weights : 0;
    return e;
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
    // FLASHDeconv's MS1 tolerance, minimum isotope cosine and minimum charge SNR (its defaults)
    const double fd_tol = sd_param.getValue("tol").toDoubleVector()[0] * 1e-6, fd_min_cos = sd_param.getValue("min_cos").toDoubleVector()[0],
                 fd_min_snr = sd_param.getValue("min_snr").toDoubleVector()[0];
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

    // every run: the features and, where not detected, their decoys; their raw abundance and FLASHDeconv's evidence for their masses
    std::vector<std::vector<FLASHQuantAlgorithm::MassAbundance>> value(kept.size(), std::vector<FLASHQuantAlgorithm::MassAbundance>(run_count));
    std::vector<std::vector<double>> decoy_value(kept.size(), std::vector<double>(run_count, 0));
    std::vector<std::vector<std::array<FDEvidence, 2>>> fd_evidence(kept.size(), std::vector<std::array<FDEvidence, 2>>(run_count)); // target, decoy
    for (Size r = 0; r < run_count; ++r)
    {
      MzMLFile mzml;
      mzml.setLogType(log_type_);
      mzml.getOptions().setMSLevels({1});
      PeakMap map;
      mzml.load(in[r], map);
      map.sortSpectra(true);
      std::vector<FLASHQuantAlgorithm::QuantTarget> targets;
      std::vector<std::pair<Size, bool>> owner; // (consensus, decoy)
      for (Size k = 0; k < kept.size(); ++k)
      {
        const auto& s = summary[k];
        const double center = (s.aligned + (anchors[r].empty() ? 0.0 : offset_(anchors[r], s.aligned))) * 60.0, half = s.half * 60.0;
        targets.push_back({s.mass, center - half, center + half, s.zmin, s.zmax});
        owner.emplace_back(k, false);
        if (kept[k].runs.count(r) == 0 && mbr_fdr >= 0)
        {
          targets.push_back({s.mass - decoy_shift[k], center - half, center + half, s.zmin, s.zmax});
          owner.emplace_back(k, true);
        }
      }
      FLASHQuantAlgorithm fq;
      std::vector<FLASHQuantAlgorithm::MassAbundance> abundances;
      fq.quantifyTargets(map, targets, avg, abundances);
      std::vector<FDEvidence> evidence(targets.size());
      if (mbr_fdr >= 0)
      {
        std::vector<double> rts(map.size());
        for (Size s = 0; s < map.size(); ++s) rts[s] = map[s].getRT();
#pragma omp parallel for schedule(dynamic, 16)
        for (SignedSize t = 0; t < (SignedSize)targets.size(); ++t) evidence[t] = fdEvidence_(map, rts, targets[t], avg, fd_tol, fd_min_cos, fd_min_snr);
      }
      for (Size t = 0; t < targets.size(); ++t)
      {
        const auto [k, decoy] = owner[t];
        if (decoy) decoy_value[k][r] = abundances[t].first;
        else value[k][r] = abundances[t];
        fd_evidence[k][r][decoy ? 1 : 0] = evidence[t];
      }
      OPENMS_LOG_INFO << "Run " << (r + 1) << " of " << run_count << ": " << targets.size() << " targets measured" << endl;
    }

    // The score of a transfer: FLASHDeconv's Qscore2D of its mass in its window, if FLASHDeconv's mass of it is within three times the
    // median mass error of the features where detected and it has a raw abundance; else 0 (never kept). Its decoy is scored the same way.
    // FLASHDeconv's mass, from all the peaks of the envelope in a scan, is far more precise than any single trace. The error is taken
    // against the feature's own mass measured the same way where detected (the median signed error of those runs; for a detected run, of
    // the other ones): the mass of the FLASHDeconv feature table is another estimate, off by a fixed amount for a quarter of the features.
    // PXD042169: median |error| 0.12 ppm where detected, 2-3 ppm for the decoys; twice, three and four times the median kept ~70, ~80
    // and ~82% of the detected features at 1% FDR.
    std::vector<double> reference(kept.size(), 0); // per feature: the median signed error (ppm) of its mass where detected
    double mass_gate = 0;
    {
      std::vector<double> errors;
      for (Size k = 0; k < kept.size(); ++k)
      {
        std::vector<std::pair<Size, double>> own; // (run, signed error) where detected
        for (const auto& [r, i] : kept[k].runs)
        {
          if (fd_evidence[k][r][0].scans > 0) own.emplace_back(r, fd_evidence[k][r][0].ppm);
        }
        std::vector<double> all;
        for (const auto& [r, e] : own) all.push_back(e);
        reference[k] = median_(all);
        for (const auto& [r, e] : own)
        {
          std::vector<double> others;
          for (const auto& [r2, e2] : own) { if (r2 != r) others.push_back(e2); }
          if (!others.empty()) errors.push_back(std::abs(e - median_(others)));
        }
      }
      mass_gate = 3 * median_(errors);
      if (mbr_fdr >= 0) OPENMS_LOG_INFO << "Match between runs: |mass error| gate from the detected features: " << mass_gate << " ppm" << endl;
    }
    const auto score_of = [&](const FDEvidence& e, double v, double ref) {
      return (v > 0 && e.scans > 0 && std::abs(e.ppm - ref) <= mass_gate) ? e.qscore2d : 0.0;
    };
    std::vector<std::vector<double>> score(kept.size(), std::vector<double>(run_count, 0));
    std::vector<double> target_scores, null_scores;
    for (Size k = 0; k < kept.size() && mbr_fdr >= 0; ++k)
    {
      for (Size r = 0; r < run_count; ++r)
      {
        if (kept[k].runs.count(r) > 0) continue;
        score[k][r] = score_of(fd_evidence[k][r][0], value[k][r].first, reference[k]);
        target_scores.push_back(score[k][r]);
        null_scores.push_back(score_of(fd_evidence[k][r][1], decoy_value[k][r], reference[k]));
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
      ds << "ConsensusIndex\tRun\tStatus\tValue\tQscore2D\tMaxQscore\tScans\tMassErrorPPM\tDecoyValue\tDecoyQscore2D\tDecoyMaxQscore\tDecoyScans\tDecoyMassErrorPPM"
            "\tReferenceErrorPPM\tDetectedMedian\n";
      for (Size k = 0; k < kept.size(); ++k)
      {
        std::vector<double> det;
        for (const auto& [r, i] : kept[k].runs) det.push_back(value[k][r].first);
        const double det_median = median_(det);
        for (Size r = 0; r < run_count; ++r)
        {
          const bool detected = kept[k].runs.count(r) > 0;
          ds << k << "\t" << r << "\t" << (detected ? "detected" : "transfer") << "\t" << std::to_string(value[k][r].first);
          for (int kind = 0; kind < 2; ++kind)
          {
            const auto& e = fd_evidence[k][r][kind];
            if (kind == 1) ds << "\t" << std::to_string(decoy_value[k][r]);
            ds << "\t" << std::to_string(e.qscore2d) << "\t" << std::to_string(e.max_qscore) << "\t" << e.scans << "\t" << std::to_string(e.ppm);
          }
          ds << "\t" << std::to_string(reference[k]) << "\t" << std::to_string(det_median) << "\n";
        }
      }
    }

    // output
    std::ofstream os(out);
    os << "ConsensusIndex\tMonoisotopicMass\tApexRetentionTime\tMinCharge\tMaxCharge\tDetectedRuns";
    std::vector<String> names;
    for (const auto& f : in) names.push_back(QFileInfo(f.toQString()).completeBaseName());
    for (const auto& n : names) os << "\tAbundance_" << n << "\tStatus_" << n << "\tMatchScore_" << n << "\tMatchQvalue_" << n;
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
        const double q = (mbr_fdr >= 0 && v > 0) ? q_of(score[k][r]) : 1.0;
        const bool ok = mbr_fdr >= 0 && v > 0 && q <= mbr_fdr;
        mbr_kept += ok;
        os << "\t" << (ok ? std::to_string(v) : String("nan")) << "\t" << (ok ? "mbr" : "nan") << "\t" << (mbr_fdr >= 0 ? std::to_string(score[k][r]) : String(""))
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
