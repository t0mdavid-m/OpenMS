// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim, Jaekwan Kim $
// $Authors: Kyowon Jeong, Jihyung Kim, Jaekwan Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandJoin.h>
#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/Qvalue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/TopDownIsobaricQuantification.h>
#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>

#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{
inline const Size max_peak_count_for_centroid_ = 3e4;
inline const Size max_peak_count_for_profile_ = 1e5;

FLASHDeconvAlgorithm::FLASHDeconvAlgorithm(): DefaultParamHandler("FLASHDeconvAlgorithm"), ProgressLogger()
{
  // ADR-0046: a commanded MS2 is LOCATED through the acquisition's scan_commands.tsv, not searched for.
  defaults_.setValue("scan_commands", "", "Provide the scan_commands.tsv written by FLASHIda (in the run folder <raw name>_<yyyy-MM-dd-HH-mm-ss>) to couple this mzML with its FLASHIda acquisition. Needs an mzML that keeps each spectrum's 'scan description' (msconvert does). Replaces the former ida_log parameter.");
  defaults_.setValue(
    "report_FDR", "false",
    "Generate q-values (approximate FDR) for deconvolved masses. Decoy masses are also reported which were used for FDR calculation. (Beta version)");
  defaults_.setValidStrings("report_FDR", {"true", "false"});

  defaults_.setValue("allowed_isotope_error", 0,
                     "Tolerance for isotope index errors when calculating FDR. For instance, setting a value of 2 permits the inclusion of up to 2 isotope errors as valid matches. Beta version.");
  defaults_.addTag("allowed_isotope_error", "advanced");

  defaults_.setValue("use_RNA_averagine", "false", "Use the RNA (nucleotide) averagine model for deconvolution.");
  defaults_.setValidStrings("use_RNA_averagine", {"true", "false"});
  defaults_.addTag("use_RNA_averagine", "advanced");

//  defaults_.setValue(
//    "precursor_MS1_window", 1,
//    "Number of MS1 spectra around each MS2 spectrum to search for precursor peaks when determining the MS2 precursors. For MS2 spectrum, the mass of precursor ion should be determined for better deconvolution and reliable identification. "
//    "If the mass of precursor ion is not found in the immediately preceding MS1 spectrum, previous or next MS1 spectra may be used instead. "
//    "This parameter determines up to how many MS1 spectra around each MS2 spectrum will be searched.");
//  defaults_.setMinInt("precursor_MS1_window", 1);
//  defaults_.addTag("precursor_MS1_window", "advanced");

  defaults_.setValue(
    "isolation_window", 5.0,
    "Specify the isolation window width for precursor determination. Used when this information is absent in the mzML file.");
  defaults_.addTag("isolation_window", "advanced");

  defaults_.setValue("precursor_rescue_window", 2,
                     "Number of MS1 spectra before and after the precursor (preceding MS1) spectrum of an MS2 spectrum that are searched for its "
                     "precursor mass when the precursor spectrum has no deconvolved mass in the isolation window: in an LC-MS run the precursor "
                     "elutes over several MS1 scans, so a mass missing from one scan is often in the neighbours (the best charge SNR in the "
                     "window wins). 0 turns the rescue off, e.g. for native or targeted data.");
  defaults_.setMinInt("precursor_rescue_window", 0);
  defaults_.addTag("precursor_rescue_window", "advanced");

  defaults_.setValue("merging_method", 0,
                     "Method for merging spectra before deconvolution. 0: No merging  1: Gaussian averaging per MS level, effective for Q-TOF datasets. For MSn (n > 1), only the spectra from the same precursor mass "
                     "(subject to tolerance set by SD:tol) are averaged. 2: Block merging, combining all spectra into one per MS level (e.g., for NativeMS datasets).");
  defaults_.setMinInt("merging_method", 0);
  defaults_.setMaxInt("merging_method", 2);

  defaults_.setValue("merging_min_ms_level", 1, "Min MS level for merging");
  defaults_.setValue("merging_max_ms_level", 2, "Max MS level for merging");
  defaults_.setMinInt("merging_min_ms_level", 1);
  defaults_.setMinInt("merging_max_ms_level", 1);

  auto sd_defaults = SpectralDeconvolution().getDefaults();
  sd_defaults.remove("allowed_isotope_error");
  defaults_.insert("SD:", sd_defaults);

  Param mf_defaults = MassFeatureTrace().getDefaults();
  mf_defaults.setValue("min_cos", -1.0,
                       "Cosine similarity threshold between avg. and observed isotope pattern. When negative, MS1 cosine threshold for spectral "
                       "deconvolution (set by -SD:min_cos will be used ");
  mf_defaults.setValue(
    "mass_error_ppm", -1.0,
    "Specifies the mass error tolerance for feature tracing in ppm. When negative, the MS1 tolerance for deconvolution is used (e.g., 16 ppm is used when -SD:tol 16).");
  mf_defaults.addTag("min_cos", "advanced");
  mf_defaults.addTag("mass_error_ppm", "advanced");

  mf_defaults.remove("noise_threshold_int");
  mf_defaults.remove("reestimate_mt_sd");
  mf_defaults.remove("trace_termination_criterion");
  mf_defaults.remove("trace_termination_outliers");
  mf_defaults.remove("chrom_peak_snr");

  defaults_.insert("ft:", mf_defaults);
  defaults_.insert("iq:", TopDownIsobaricQuantification().getDefaults());

  defaultsToParam_();
}

void FLASHDeconvAlgorithm::updateMembers_()
{
  tols_ = param_.getValue("SD:tol");
  min_cos_ = param_.getValue("SD:min_cos");
  //precursor_MS1_window_ = param_.getValue("precursor_MS1_window");
  use_RNA_averagine_ = param_.getValue("use_RNA_averagine") != "false";
  report_decoy_ = param_.getValue("report_FDR") != "false";

  merge_spec_ = param_.getValue("merging_method");
  isolation_window_size_ = param_.getValue("isolation_window");
  scan_commands_file_ = param_.getValue("scan_commands").toString();
  precursor_rescue_window_ = param_.getValue("precursor_rescue_window");
  current_min_ms_level_ = max_ms_level_;
  current_max_ms_level_ = 0;
}

void FLASHDeconvAlgorithm::updateMSLevels_(MSExperiment& map)
{
  // read input dataset once to count spectra
  for (auto& it : map)
  {
    // if forced_ms_level > 0, force MS level of all spectra to 1.
    if (forced_ms_level_ > 0) { it.setMSLevel(forced_ms_level_); }

    if (it.empty()) { continue; }
    if (it.getMSLevel() > max_ms_level_) { continue; }

    uint ms_level = it.getMSLevel();
    current_max_ms_level_ = current_max_ms_level_ < ms_level ? ms_level : current_max_ms_level_;
    current_min_ms_level_ = current_min_ms_level_ > ms_level ? ms_level : current_min_ms_level_;
  }
  // Max MS Level is adjusted according to the input dataset
  current_max_ms_level_ = current_max_ms_level_ > max_ms_level_ ? max_ms_level_ : current_max_ms_level_;
}

void FLASHDeconvAlgorithm::filterLowPeaks_(MSExperiment& map)
{
  OPENMS_LOG_INFO << "Filtering low peaks in spectra ... ";
  ThresholdMower threshold_mower_filter;                         // threshold
  Param t_filter_param = threshold_mower_filter.getParameters(); //"threshold", .00001
  t_filter_param.setValue("threshold", 1e-6);
  threshold_mower_filter.setParameters(t_filter_param);
  threshold_mower_filter.filterPeakMap(map);

#pragma omp parallel for default(none), shared(map)
  for (int i = 0; i < (int) map.size(); i++)
  {
    auto& it = map[i];
    if (it.empty()) continue;
    Size count = it.getType(false) == SpectrumSettings::CENTROID ? max_peak_count_for_centroid_ : max_peak_count_for_profile_;
    it.sortByIntensity(true);
    double threshold = it.size() < count ? 0 : it[count].getIntensity();
    threshold = std::max(threshold, (double)it.begin()->getIntensity() / 1000);
    // pop back the low intensity peaks using threshold

    while (! it.empty() && it.back().getIntensity() <= threshold)
    {
      it.pop_back();
    }
    it.sortByPosition();
  }
  OPENMS_LOG_INFO << "Done" << std::endl;
}

void FLASHDeconvAlgorithm::mergeSpectra_(MSExperiment& map, uint ms_level)
{
  SpectraMerger merger;
  merger.setLogType(CMD);
  Param sm_param = merger.getDefaults();
  sm_param.setValue("mz_binning_width", tols_[ms_level - 1] / 2.5);
  sm_param.setValue("mz_binning_width_unit", "ppm");
  uint min_ms_level = param_.getValue("merging_min_ms_level");
  uint max_ms_level = param_.getValue("merging_max_ms_level");
  if (merge_spec_ == 1 && ms_level >= min_ms_level && ms_level <= max_ms_level)
  {
    if (ms_level == 1)
    {
      OPENMS_LOG_INFO << "Gaussian averaging MS1 spectra ... " << std::endl;
      merger.setParameters(sm_param);
      map.sortSpectra();
      merger.average(map, "gaussian", (int)ms_level);
    }
    else
    {
      // For ms n, first find precursors for all ms n. then make a tmp map having the precursor masses as precursor
      std::map<String, std::vector<Precursor>> original_precursor_map;
//#pragma omp parallel for default(none), shared(map, ms_level, original_precursor_map)
      for (int i = 0; i < (int) map.size(); i++)
      {
        const auto& spec = map[i];
        if (spec.getMSLevel() != ms_level) continue;

        const auto& native_id = spec.getNativeID();
        original_precursor_map[native_id] = spec.getPrecursors();

        if (! spec.getPrecursors().empty() && native_id_precursor_peak_group_map_.find(native_id) != native_id_precursor_peak_group_map_.end())
        {
          auto precursor_pg = native_id_precursor_peak_group_map_[native_id];
          auto precursor = spec.getPrecursors()[0];
          precursor.setCharge(1);
          precursor.setMZ(precursor_pg.getMonoMass());
          precursor.setIntensity(precursor_pg.getIntensity());
          map[i].setPrecursors(std::vector<Precursor> {precursor});
        }
      }
      // merge MS n using precursor method
      OPENMS_LOG_INFO << "Merging MS" << ms_level << " spectra from the same deconvolved precursor masses... " << std::endl;
      sm_param.setValue("precursor_method:mz_tolerance", 0.2);
      sm_param.setValue("precursor_method:rt_tolerance", 30.0); // TODO make this as a user input?
      merger.setParameters(sm_param);
      map.sortSpectra();
      merger.mergeSpectraPrecursors(map);

      for (auto& mspec : map)
      {
        auto native_id_str = mspec.getNativeID();
        std::vector<String> native_ids;
        native_id_str.split(",", native_ids);

        for (auto& native_id : native_ids)
        {
          if (original_precursor_map.find(native_id) == original_precursor_map.end()) continue;
          mspec.setPrecursors(original_precursor_map[native_id]);
        }
      }
    }
  }
  else if (merge_spec_ == 2)
  {
    OPENMS_LOG_INFO << "Merging spectra into a single spectrum for MS" << ms_level << std::endl;
    sm_param.setValue("block_method:rt_block_size", map.size() + 1);
    map.sortSpectra();
    sm_param.setValue("block_method:ms_levels", IntList {(int)ms_level});
    merger.setParameters(sm_param);
    merger.mergeSpectraBlockWise(map);
  }
  filterLowPeaks_(map);
}

int FLASHDeconvAlgorithm::getScanNumber(const MSExperiment& map, Size index)
{
  auto native_id_str = map[index].getNativeID();
  std::vector<String> native_ids;
  native_id_str.split(",", native_ids);
  String type_accession = "MS:1000768";

  if (!map.getSourceFiles().empty())
  {
    type_accession = map.getSourceFiles()[0].getNativeIDTypeAccession();
    if (type_accession.empty()) type_accession = "MS:1000768";
  }
  int scan_number = SpectrumLookup::extractScanNumber(native_ids.back(), type_accession);
  if (scan_number < 0) { scan_number = (int)index + 1; }

  return scan_number;
}

std::vector<double> FLASHDeconvAlgorithm::getTolerances() const
{
  return tols_;
}

void FLASHDeconvAlgorithm::runSpectralDeconvolution_(MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra)
{
  startProgress(0, (SignedSize)map.size(), "running FLASHDeconv");
  std::map<double, int> rt_scan_map;
  for (Size index = 0; index < map.size(); index++)
  {
    int scan_number = getScanNumber(map, index);
    rt_scan_map[map[index].getRT()] = scan_number;
  }

  for (uint ms_level = 1; ms_level <= current_max_ms_level_; ms_level++)
  {
    if (ms_level > 1)
    {
      // here, register precursor peak groups to the ms2 spectra.
      findPrecursorPeakGroupsForMSnSpectra_(map, deconvolved_spectra, ms_level);
    }
    if (merge_spec_ > 0)
    {
      mergeSpectra_(map, ms_level);
      // repeat precursor peak group registration for merged/averaged spectra
      if (ms_level > 1 && merge_spec_ > 0)
      {
        // here, register precursor peak groups to the ms2 spectra.
        findPrecursorPeakGroupsForMSnSpectra_(map, deconvolved_spectra, ms_level);
      }
    }

    // run Spectral deconvolution
    for (Size index = 0; index < map.size(); index++)
    {
      int scan_number = merge_spec_ == 0? getScanNumber(map, index) :
                                         (rt_scan_map.find(map[index].getRT()) == rt_scan_map.end()? getScanNumber(map, index) :
                                                                                                    rt_scan_map[map[index].getRT()]);//getScanNumber(map, index);
      const auto& spec = map[index];

      if (ms_level != spec.getMSLevel()) { continue; }

      nextProgress();
      if (spec.empty()) { continue; }
      String native_id = spec.getNativeID();

      PeakGroup precursor_pg;
      if (native_id_precursor_peak_group_map_.find(native_id) != native_id_precursor_peak_group_map_.end())
        precursor_pg = native_id_precursor_peak_group_map_[native_id];

      // now do it
      sd_.performSpectrumDeconvolution(spec, scan_number, precursor_pg);

      auto& deconvolved_spectrum = sd_.getDeconvolvedSpectrum();

      // add precursor information
      if (native_id_precursor_peak_group_map_.find(native_id) != native_id_precursor_peak_group_map_.end())
      {
        deconvolved_spectrum.setDeconvolvedPrecursor(native_id_precursor_peak_map_[native_id]);
        auto decoy_it = native_id_decoy_precursor_mass_map_.find(native_id);
        if (decoy_it != native_id_decoy_precursor_mass_map_.end())
        {
          deconvolved_spectrum.setDecoyPrecursorMass(decoy_it->second);
          // a noise decoy counts with the noise decoy weight, known after the q-values below
          deconvolved_spectrum.setDecoyPrecursorWeight(native_id_decoy_precursor_noise_map_[native_id] ? -1.0 : 1.0);
        }
      }

      if (native_id_precursor_mass_intensity_map_.find(native_id) != native_id_precursor_mass_intensity_map_.end())
      {
        deconvolved_spectrum.setPrecursorMassIntensityMap(native_id_precursor_mass_intensity_map_[native_id]);
      }
      // add precursor scan number information when a precursor peak group is not found
      if (ms_level > 1 && precursor_pg.empty())
      {
        for (int p_index = (int)index - 1; p_index >= 0; p_index--)
        {
          if(map[p_index].getMSLevel() == ms_level - 1)
          {
            deconvolved_spectrum.setPrecursorScanNumber(getScanNumber(map, p_index));
            break;
          }
        }
      }

      if (report_decoy_ && !deconvolved_spectrum.empty())
      {
#pragma omp parallel sections default(none) shared(spec, scan_number, precursor_pg, deconvolved_spectrum)
        {
#pragma omp section
          sd_noise_decoy_.performSpectrumDeconvolution(spec, scan_number, precursor_pg);
#pragma omp section
          sd_signal_decoy_.performSpectrumDeconvolution(spec, scan_number, precursor_pg);
        }
        deconvolved_spectrum.sortByQscore();

        deconvolved_spectrum.reserve(deconvolved_spectrum.size() + sd_signal_decoy_.getDeconvolvedSpectrum().size()
                                     + sd_noise_decoy_.getDeconvolvedSpectrum().size());

        for (const auto& pg : sd_signal_decoy_.getDeconvolvedSpectrum())
            deconvolved_spectrum.push_back(pg);

        for (const auto& pg : sd_noise_decoy_.getDeconvolvedSpectrum())
            deconvolved_spectrum.push_back(pg);

        deconvolved_spectrum.sort();
      }
      deconvolved_spectra.push_back(deconvolved_spectrum);
    }
    std::sort(deconvolved_spectra.begin(), deconvolved_spectra.end());
  }
  endProgress();
}


const FLASHHelperClasses::PrecalculatedAveragine& FLASHDeconvAlgorithm::getAveragine()
{
  return sd_.getAveragine();
}

const FLASHHelperClasses::PrecalculatedAveragine& FLASHDeconvAlgorithm::getDecoyAveragine()
{
  if (!report_decoy_) return getAveragine();
  return sd_noise_decoy_.getAveragine();
}

std::vector<int> FLASHDeconvAlgorithm::getHistogram_(const std::vector<double>& data, double min_range, double max_range, double bin_size)
{
  int num_bins = static_cast<int>((max_range - min_range) / bin_size) + 1;
  std::vector<int> bins(num_bins, 0);

  // Populate the bins
  for (double value : data) {
    if (value >= min_range && value <= max_range) {
      int bin_index = static_cast<int>((value - min_range) / bin_size);
      bins[bin_index]++;
    }
  }
  return bins;
}

void FLASHDeconvAlgorithm::determineTolerance_(const MSExperiment& map, const Param& sd_param, const FLASHHelperClasses::PrecalculatedAveragine& avg, const uint ms_level)
{
  OPENMS_LOG_INFO << "Determining tolerance for MS" << ms_level << " ... ";
  auto sd = SpectralDeconvolution();
  auto sd_param_t = sd_param;
  sd.setAveragine(avg);
  tols_[ms_level - 1] = 200; // maximum tolerance
  sd_param_t.setValue("min_charge", 1); // better to include charge 1 to determine ppm error.
  sd_param_t.setValue("min_mass", 50.0); // better to include small masses to determine ppm error.
  sd_param_t.setValue("tol", tols_);
  sd.setParameters(sd_param_t);
  sd.setToleranceEstimation();

  int sample_rate = 100; //
  int count = 0;
  std::vector<double> sampled_tols;
  for (const auto& spec : map)
  {
    if (ms_level != (uint)spec.getMSLevel()) { continue; }
    if (spec.empty()) { continue; }
    if (count++ % sample_rate != 0) continue;
    PeakGroup precursor_pg;
    sd.performSpectrumDeconvolution(spec, 0, precursor_pg);

    const auto& deconvolved_spectrum = sd.getDeconvolvedSpectrum(); // estimate both.

    for (const auto& pg : deconvolved_spectrum)
    {
      if (pg.getQscore() < .9) // TODO automatically find the good mass threshold
        continue;
      for (auto error : pg.getMassErrors())
      {
        sampled_tols.push_back(error);
      }
    }
  }

  if (sampled_tols.size() < 6)
  {
    OPENMS_LOG_INFO << "failed. Cannot be determined - no MS" << ms_level << " spectrum. Set to 10ppm (default tolerance)." << std::endl;
    tols_[ms_level - 1] = 10;
    return;
  }

  try
  {
    const int l = -100, r = 100, bin_size = 1;
    const auto& bins = getHistogram_(sampled_tols, l, r, bin_size);
    // Calculate mean using std::accumulate

    std::vector<DPosition<2>> points;
    for (int b = l; b <= r; b += bin_size)
    {
      points.emplace_back(b, bins[b - l]);
    }
    Math::GaussFitter fitter;
    const auto fit = fitter.fit(points);

    tols_[ms_level - 1] = round(fit.sigma * 2.17 * 2); // 97% area under curve in gaussian
    OPENMS_LOG_INFO << "done. Determined tolerance: " << tols_[ms_level - 1] << " ppm. You may test around this tolerance for better results."
                    << std::endl;
  }
  catch (const Exception::UnableToFit& e)
  {
    // Handle the exception and issue a warning
    OPENMS_LOG_INFO << "failed. Cannot be determined - Gaussian fitting failure. Set to 10ppm (default tolerance)." << std::endl;
    tols_[ms_level - 1] = 10;
  }
  catch (const std::exception& e)
  {
    // Handle other exceptions if necessary
    OPENMS_LOG_INFO << "failed. Cannot be determined - Gaussian fitting failure. Set to 10ppm (default tolerance)." << std::endl;
    tols_[ms_level - 1] = 10;
  }
}

void FLASHDeconvAlgorithm::run(MSExperiment& map,
                               std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                               std::vector<FLASHHelperClasses::MassFeature>& deconvolved_features)
{
  // ADR-0046: join this data file to its acquisition BEFORE any deconvolution, so a scan_commands file that
  // does not belong to it aborts in seconds. Metadata only -- nothing below has touched the spectra yet.
  located_by_scan_.clear();
  if (!scan_commands_file_.empty())
  {
    std::vector<ScanCommandJoin::Scan> scans;
    scans.reserve(map.size());
    for (Size index = 0; index < map.size(); index++)
    {
      ScanCommandJoin::Scan s;
      s.scan_number = getScanNumber(map, index);
      s.ms_level = (int)map[index].getMSLevel();
      if (map[index].metaValueExists("scan description")) { s.scan_description = map[index].getMetaValue("scan description").toString(); }
      for (const auto& precursor : map[index].getPrecursors())
      {
        s.isolation_targets.push_back(precursor.getMZ());
        // MzMLHandler overwrites getMZ() with the selected-ion m/z when the two differ, and keeps the target here
        if (precursor.metaValueExists("isolation window target m/z")) { s.isolation_targets.push_back((double)precursor.getMetaValue("isolation window target m/z")); }
      }
      scans.push_back(s);
    }
    located_by_scan_ = ScanCommandJoin::join(scans, ScanCommandJoin::parse(scan_commands_file_));
    OPENMS_LOG_INFO << "FLASHIda scan commands used: " << scan_commands_file_ << " (" << located_by_scan_.size() << " commanded MS2 scans located)" << std::endl;
  }

  updateMSLevels_(map);
  filterLowPeaks_(map);

  sd_ = SpectralDeconvolution();
  Param sd_param = param_.copy("SD:", true);
  sd_param.setValue("allowed_isotope_error", param_.getValue("allowed_isotope_error"));
  OPENMS_LOG_INFO<< "Calculating Averagines ... " << std::flush;
  sd_.setParameters(sd_param);
  sd_.calculateAveragine(use_RNA_averagine_);
  OPENMS_LOG_INFO<< "Done" << std::endl;
  const auto& avg = sd_.getAveragine();

  // determine tolerance in case tolerance input is negative
  for (uint ms_level = 1; ms_level <= current_max_ms_level_; ms_level++)
  {
    if (tols_[ms_level - 1] > 0) continue;
    determineTolerance_(map, sd_param, avg, ms_level);
  }

  sd_param.setValue("tol", tols_);
  sd_.setParameters(sd_param);

  if (report_decoy_)
  {
    sd_noise_decoy_.setParameters(sd_param);
    sd_noise_decoy_.setTargetDecoyType(PeakGroup::TargetDecoyType::noise_decoy, sd_.getDeconvolvedSpectrum()); // noise
    sd_noise_decoy_.calculateAveragine(use_RNA_averagine_); // for noise, averagine needs to be calculated differently.

    sd_signal_decoy_.setParameters(sd_param);
    sd_signal_decoy_.setTargetDecoyType(PeakGroup::TargetDecoyType::signal_decoy, sd_.getDeconvolvedSpectrum()); // isotope
    sd_signal_decoy_.setAveragine(avg);
  }

  setLogType(CMD);
  deconvolved_spectra.reserve(map.size() * 4);

  // run spectral deconvolution here and get deconvolved spectra
  runSpectralDeconvolution_(map, deconvolved_spectra);

  // feature tracing here and update FeatureQScores
  runFeatureFinding_(deconvolved_spectra, deconvolved_features);

  noise_decoy_weight_ = Qvalue::updatePeakGroupQvalues(deconvolved_spectra);
  for (auto& dspec : deconvolved_spectra)
  {
    if (dspec.getDecoyPrecursorMass() > 0 && dspec.getDecoyPrecursorWeight() < 0) dspec.setDecoyPrecursorWeight(noise_decoy_weight_);
  }

  TopDownIsobaricQuantification quantifier;
  Param quant_param = param_.copy("iq:", true);
  quantifier.setParameters(quant_param);
  // Isobaric quant run
  quantifier.quantify(map, deconvolved_spectra, deconvolved_features);
}

// ADR-0046 decision 3: FLASHDeconv's deconvolution of the commanded survey holds NOTHING in the isolation window,
// so the precursor is rebuilt from the command itself -- one peak at the anchor charge, scores as the engine
// logged them. The construction the ida.log path used, over doubles instead of floats, and with no charge range
// to claim. A window that holds something the command did not name takes the charge-SNR pick instead: a measured
// species carries an isotope envelope, and this one does not -- which also bounds the FRAGMENT deconvolution
// (SpectralDeconvolution.cpp:197-218), so a rebuilt precursor costs the MS2 as well as the mass.
PeakGroup FLASHDeconvAlgorithm::peakGroupFromCommand_(const ScanCommandJoin::Located& located) const
{
  const auto& row = located.row;
  FLASHHelperClasses::LogMzPeak precursor_log_mz_peak;
  precursor_log_mz_peak.abs_charge = std::abs(row.charge);
  precursor_log_mz_peak.is_positive = row.charge > 0;
  precursor_log_mz_peak.isotopeIndex = 0;
  precursor_log_mz_peak.mass = row.mono_mass;
  precursor_log_mz_peak.intensity = (float)row.peakgroup_intensity;

  PeakGroup precursor_pg(precursor_log_mz_peak.abs_charge, precursor_log_mz_peak.abs_charge, true);
  precursor_pg.push_back(precursor_log_mz_peak);
  precursor_pg.setChargeIsotopeCosine(precursor_log_mz_peak.abs_charge, (float)row.charge_cos);
  precursor_pg.setChargeSNR(precursor_log_mz_peak.abs_charge, (float)row.charge_snr);
  precursor_pg.setIsotopeCosine((float)row.iso_cos);
  precursor_pg.setSNR((float)row.snr);
  precursor_pg.setChargeScore((float)row.charge_score);
  precursor_pg.setAvgPPMError((float)row.ppm_error);
  precursor_pg.setQscore(row.qscore);
  precursor_pg.setRepAbsCharge(precursor_log_mz_peak.abs_charge);
  precursor_pg.updateMonoMassAndIsotopeIntensities(tols_[0]);
  precursor_pg.setScanNumber(located.parent_scan_number);
  return precursor_pg;
}

void FLASHDeconvAlgorithm::findPrecursorPeakGroupsForMSnSpectra_(const MSExperiment& map,
                                                                 const std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                                                 uint ms_level)
{
  for (Size index = 0; index < map.size(); index++)
  {
    auto spec = map[index]; // MS2 index
    if (spec.getMSLevel() != ms_level) { continue; }

    int scan_number = getScanNumber(map, index);
    String native_id = spec.getNativeID();

    // ADR-0046: a commanded MS2 is located by its command. Null for an uncommanded scan, for MSn > 2, and when
    // the survey the command names is not in this data file (an RT crop) -- those keep the search below unchanged.
    const ScanCommandJoin::Located* located = nullptr;
    if (ms_level == 2)
    {
      auto lit = located_by_scan_.find(scan_number);
      if (lit != located_by_scan_.end() && lit->second.parent_scan_number > 0) { located = &lit->second; }
    }

    auto cv_of = [](const MSSpectrum& sp) {
      auto fs = sp.getMetaValue("filter string").toString();
      Size p = fs.find("cv=");
      if (p == String::npos) return 1e5;
      Size e = fs.find(" ", p);
      if (e == String::npos) e = fs.length() - 1;
      return std::stod(fs.substr(p + 3, e - p));
    };
    const double cv = cv_of(spec);

    // The survey (MS n-1) spectra searched for the precursor. A located scan reads the survey its COMMAND names, however many newer
    // surveys were acquired in between (ADR-0046; 88 % of commanded MS2 have at least one), and only that one. Otherwise the preceding
    // survey of the same CV, then the neighbours by distance (one before, one after, two before, ...) up to precursor_rescue_window_
    // away: an LC-MS precursor elutes over several survey scans, so a mass missing from the preceding scan is often in a neighbour
    // (TopFD takes its precursors from LC-MS features).
    std::vector<Size> survey_indices;
    if (located != nullptr)
    {
      for (int i = (int)index - 1; i >= 0; i--)
      {
        if (map[i].getMSLevel() == ms_level - 1 && getScanNumber(map, (Size)i) == located->parent_scan_number)
        {
          survey_indices.push_back((Size)i);
          break;
        }
      }
    }
    else
    {
      std::vector<Size> before, after;
      for (int i = (int)index - 1; i >= 0 && (int)before.size() <= precursor_rescue_window_; i--)
      {
        if (map[i].getMSLevel() == ms_level - 1 && cv_of(map[i]) == cv) before.push_back((Size)i);
      }
      for (Size i = index + 1; i < map.size() && (int)after.size() < precursor_rescue_window_; i++)
      {
        if (map[i].getMSLevel() == ms_level - 1 && cv_of(map[i]) == cv) after.push_back(i);
      }
      if (! before.empty())
      {
        survey_indices.push_back(before[0]);
        for (int k = 1; k <= precursor_rescue_window_; k++)
        {
          if ((Size)k < before.size()) survey_indices.push_back(before[k]);
          if ((Size)(k - 1) < after.size()) survey_indices.push_back(after[k - 1]);
        }
      }
    }
    if (survey_indices.empty()) continue;

    // the deconvolved spectrum of a survey spectrum, by its RT (deconvolved_spectra are sorted by scan number)
    auto deconvolved_of = [&](Size survey_index) -> const DeconvolvedSpectrum* {
      const double rt = map[survey_index].getRT();
      auto it = std::lower_bound(deconvolved_spectra.begin(), deconvolved_spectra.end(), DeconvolvedSpectrum(getScanNumber(map, survey_index) + 1));
      while (it != deconvolved_spectra.begin())
      {
        --it;
        if (std::abs(it->getOriginalSpectrum().getRT() - rt) < 1e-5) return &*it;
        if (it->getOriginalSpectrum().getRT() < rt - 1e-5 && it->getOriginalSpectrum().getMSLevel() == ms_level - 1) break;
      }
      return nullptr;
    };

    // the isolation window
    double start_mz = 0;
    double end_mz = 0;
    for (auto& precursor : spec.getPrecursors())
    {
      double loffset = precursor.getIsolationWindowLowerOffset();
      double uoffest = precursor.getIsolationWindowUpperOffset();
      loffset = loffset <= 0 ? isolation_window_size_ / 2.0 : loffset;
      uoffest = uoffest <= 0 ? isolation_window_size_ / 2.0 : uoffest;

      // mzML isolation window offsets are relative to the isolation window target m/z. MzMLHandler stores the selected ion m/z as the
      // precursor m/z and, when the two differ, keeps the target as a meta value.
      double target_mz = precursor.metaValueExists("isolation window target m/z") ? (double)precursor.getMetaValue("isolation window target m/z")
                                                                                  : precursor.getMZ();
      start_mz = loffset > 100.0 ? loffset : -loffset + target_mz;
      end_mz = uoffest > 100.0 ? uoffest : uoffest + target_mz;
    }

    // The best precursor of one survey spectrum: among the target masses whose peaks lie in the window, the one of the highest charge
    // SNR (at the charge the window implies); for a located scan the COMMANDED mass wins when the survey holds it (ADR-0046: same
    // isotope only, so an isotope-off call competes as an ordinary candidate). With commit, it is registered for this MSn spectrum
    // together with the masses in the window and their intensities, and the best decoy mass of the window (with decoy reporting) as
    // the decoy precursor. Returns the SNR of the best target (0 or more), or -1 if the window holds none.
    auto scan_survey = [&](Size survey_index, bool commit) -> double {
      const DeconvolvedSpectrum* ds = deconvolved_of(survey_index);
      if (ds == nullptr) return -1.0;
      const MSSpectrum& precursor_raw_spec = map[survey_index];
      double sum_intensity = .0;
      if (! precursor_raw_spec.empty())
      {
        const Size k = precursor_raw_spec.findNearest(start_mz);
        for (Size l = k + 1; l < precursor_raw_spec.size(); l++)
        {
          if (precursor_raw_spec[l].getMZ() < start_mz) continue;
          if (precursor_raw_spec[l].getMZ() > end_mz) break;
          sum_intensity += precursor_raw_spec[l].getIntensity();
        }
        for (Size l = k;; l--)
        {
          if (precursor_raw_spec[l].getMZ() < start_mz) break;
          if (precursor_raw_spec[l].getMZ() <= end_mz) { sum_intensity += precursor_raw_spec[l].getIntensity(); }
          if (l == 0) break;
        }
      }

      double max_snr = -1.0, max_decoy_snr = -1.0;
      bool found_commanded = false; // ADR-0046: the commanded mass outranks every other species in this window
      double best_residual = 0;
      std::map<double, double> mass_to_intensity;
      mass_to_intensity[.0] = sum_intensity; // total intensity
      const PeakGroup* best_pg = nullptr;
      const FLASHHelperClasses::LogMzPeak* best_peak = nullptr;
      double decoy_mass = 0;
      bool decoy_noise = false;

      for (auto& pg : *ds)
      {
        if (pg[0].mz > end_mz || pg[pg.size() - 1].mz < start_mz) { continue; }
        // With decoy reporting the survey spectrum also holds decoy masses (noise and signal decoys). The best of them in the window,
        // chosen as the target precursor is, becomes the decoy precursor of the MSn spectrum (a decoy spectrum with that mass instead
        // of the precursor mass, e.g. for FLASHTnT -decoy_precursor); only target masses are precursors.
        if (pg.getTargetDecoyType() != PeakGroup::target)
        {
          int c = int(round(pg.getMonoMass() / start_mz));
          const FLASHHelperClasses::LogMzPeak* top = nullptr;
          for (auto& tmp_peak : pg)
          {
            if (tmp_peak.abs_charge != c || tmp_peak.mz < start_mz || tmp_peak.mz > end_mz) { continue; }
            if (top == nullptr || tmp_peak.intensity > top->intensity) top = &tmp_peak;
          }
          if (top == nullptr) { continue; }
          const double snr = pg.getChargeSNR(top->abs_charge);
          if (snr < max_decoy_snr) { continue; }
          max_decoy_snr = snr;
          decoy_mass = pg.getMonoMass();
          decoy_noise = pg.getTargetDecoyType() == PeakGroup::TargetDecoyType::noise_decoy;
          continue;
        }

        double max_intensity = .0;
        const FLASHHelperClasses::LogMzPeak* tmp_precursor = nullptr;

        int c = int(round(pg.getMonoMass() / start_mz));
        double intensity = 0;
        for (auto& tmp_peak : pg)
        {
          if (tmp_peak.abs_charge != c) { continue; }
          if (tmp_peak.mz < start_mz || tmp_peak.mz > end_mz) { continue; }
          intensity += tmp_peak.intensity;

          if (tmp_peak.intensity < max_intensity) { continue; }
          max_intensity = tmp_peak.intensity;
          tmp_precursor = &tmp_peak;
        }

        if (tmp_precursor == nullptr) { continue; }
        mass_to_intensity[pg.getMonoMass()] = intensity;

        // ADR-0046: for a located scan the COMMANDED mass wins when this survey holds it -- SAME isotope only, so an
        // isotope-off call no longer stands in for it and competes as an ordinary candidate. When the window holds
        // no commanded mass the loudest species is the better answer, so this falls through to the charge-SNR rule
        // below and peakGroupFromCommand_ is left for a window that holds nothing at all.
        bool commanded = false;
        if (located != nullptr)
        {
          double residual = 0;
          commanded = ScanCommandJoin::massRank(pg.getMonoMass(), located->row.mono_mass, residual) == 0;
          if (commanded)
          {
            if (found_commanded && residual >= best_residual) { continue; }
            found_commanded = true;
            best_residual = residual;
          }
        }
        if (! commanded)
        {
          if (found_commanded) { continue; } // never displace the commanded mass with a louder neighbour
          auto snr = pg.getChargeSNR(tmp_precursor->abs_charge); // the highest snr one should determine the mass
          if (snr < max_snr) { continue; }
          max_snr = snr;
        }
        best_pg = &pg;
        best_peak = tmp_precursor;
      }

      if (! commit) return best_pg == nullptr ? -1.0 : std::max(max_snr, 0.0);
      if (max_decoy_snr >= 0)
      {
        native_id_decoy_precursor_mass_map_[native_id] = decoy_mass;
        native_id_decoy_precursor_noise_map_[native_id] = decoy_noise;
      }
      if (best_pg != nullptr)
      {
        native_id_precursor_peak_group_map_[native_id] = *best_pg;
        Precursor precursor;
        precursor.setCharge(best_pg->getRepAbsCharge());
        precursor.setMZ(best_peak->mz);
        precursor.setIntensity(best_peak->intensity);
        native_id_precursor_peak_map_[native_id] = precursor;
      }
      auto& mi = native_id_precursor_mass_intensity_map_[native_id];
      mi.clear();
      for (const auto& [m, j] : mass_to_intensity) mi.emplace_back(m, j);
      return best_pg == nullptr ? -1.0 : std::max(max_snr, 0.0);
    };

    // the first survey spectrum as before; the neighbours only when it holds no precursor (uncommanded scans)
    double snr = scan_survey(survey_indices[0], true);
    if (snr < 0 && survey_indices.size() > 1)
    {
      double best_snr = -1.0;
      Size best_index = 0;
      for (Size k = 1; k < survey_indices.size(); k++)
      {
        const double s = scan_survey(survey_indices[k], false);
        if (s > best_snr) // ties keep the nearer scan
        {
          best_snr = s;
          best_index = survey_indices[k];
        }
      }
      if (best_snr >= 0)
      {
        native_id_decoy_precursor_mass_map_.erase(native_id);
        native_id_decoy_precursor_noise_map_.erase(native_id);
        scan_survey(best_index, true);
        rescued_precursors_++;
      }
    }
    if (located != nullptr && native_id_precursor_peak_group_map_.find(native_id) == native_id_precursor_peak_group_map_.end())
    {
      native_id_precursor_peak_group_map_[native_id] = peakGroupFromCommand_(*located);
    }
  }
  if (rescued_precursors_ > 0)
  {
    OPENMS_LOG_INFO << "Precursor masses taken from a neighbouring MS" << (ms_level - 1) << " spectrum (none in the preceding one): " << rescued_precursors_
                    << std::endl;
  }
}

void FLASHDeconvAlgorithm::updatePrecursorQScores_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, int ms_level)
{
  // update precursor feature QScores and qvalues
  std::map<int, DeconvolvedSpectrum> scan_fullscan;

  for (auto& dspec : deconvolved_spectra)
  {
    if ((int)dspec.getOriginalSpectrum().getMSLevel() != ms_level - 1) continue;
    int scan = dspec.getScanNumber();
    scan_fullscan[scan] = dspec;
  }

  for (auto& dspec : deconvolved_spectra)
  {
    if ((int)dspec.getOriginalSpectrum().getMSLevel() != ms_level) continue;
    if (dspec.getPrecursorPeakGroup().empty()) continue;

    auto precursor_pg = dspec.getPrecursorPeakGroup();

    int pscan = precursor_pg.getScanNumber();
    if (scan_fullscan.find(pscan) == scan_fullscan.end()) continue;

    auto fullscan = scan_fullscan[pscan];

    auto iter = std::lower_bound(fullscan.begin(), fullscan.end(), precursor_pg);

    if (iter == fullscan.end()) continue;

    if (precursor_pg.getMonoMass() == iter->getMonoMass())
    {
      precursor_pg.setFeatureIndex(iter->getFeatureIndex());
      precursor_pg.setQscore(iter->getQscore());
      if (iter->getFeatureIndex() > 0) precursor_pg.setQscore2D(iter->getQscore2D());
    }
    else { precursor_pg.setFeatureIndex(0); }
    dspec.setPrecursorPeakGroup(precursor_pg);
  }
}

void FLASHDeconvAlgorithm::runFeatureFinding_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                              std::vector<FLASHHelperClasses::MassFeature>& deconvolved_features)
{
  if (merge_spec_ == 2) return;

  auto mass_tracer = MassFeatureTrace();
  auto decoy_mass_tracer = MassFeatureTrace();

  Param mf_param = param_.copy("ft:", true);

  if (((double)mf_param.getValue("mass_error_ppm")) < 0) { mf_param.setValue("mass_error_ppm", tols_[0]); }
  if (((double)mf_param.getValue("min_cos")) < 0) { mf_param.setValue("min_cos", min_cos_[0]); }

  mf_param.setValue("noise_threshold_int", .0);
  mf_param.setValue("reestimate_mt_sd", "false");
  mf_param.setValue("trace_termination_criterion", "outlier");
  mf_param.setValue("trace_termination_outliers", 20);
  mf_param.setValue("chrom_peak_snr", .0);

  mass_tracer.setParameters(mf_param); // maybe go to set param
  // Find features for MS1 or the minimum MS level in the dataset.
  deconvolved_features = mass_tracer.findFeaturesAndUpdateQscore2D(sd_.getAveragine(), deconvolved_spectra, (int)current_min_ms_level_, false);

  if (report_decoy_)
  {
    const auto& decoy_deconvolved_features = mass_tracer.findFeaturesAndUpdateQscore2D(sd_.getAveragine(), deconvolved_spectra, (int)current_min_ms_level_, true);
    deconvolved_features.insert(deconvolved_features.end(), decoy_deconvolved_features.begin(), decoy_deconvolved_features.end());
  }

  mf_param.setValue("min_trace_length", 1e-5); // allow all traces for MSn
  mass_tracer.setParameters(mf_param);
  // Find features for MSn
  for (int ms_level = (int)current_min_ms_level_ + 1; ms_level <= (int)current_max_ms_level_; ms_level++)
  {
    updatePrecursorQScores_(deconvolved_spectra, ms_level);
    std::map<uint, std::vector<Size>> feature_index_set;
    for (Size i = 0; i < deconvolved_spectra.size(); i++)
    {
      const auto& dspec = deconvolved_spectra[i];
      if ((int)dspec.getOriginalSpectrum().getMSLevel() != ms_level) continue;
      if (dspec.getPrecursorPeakGroup().empty()) continue;
      uint findex = dspec.getPrecursorPeakGroup().getFeatureIndex();
      if (findex == 0) continue;
      feature_index_set[findex].push_back(i);
    }

    for (const auto& element : feature_index_set)
    {
      std::vector<DeconvolvedSpectrum> tmp_dspec;
      tmp_dspec.reserve(element.second.size());
      for (Size i : element.second)
        tmp_dspec.push_back(deconvolved_spectra[i]);

      const auto& df = mass_tracer.findFeaturesAndUpdateQscore2D(sd_.getAveragine(), tmp_dspec, ms_level, false);
      deconvolved_features.insert(deconvolved_features.end(), df.begin(), df.end());

      if (report_decoy_)
      {
        const auto& df_decoy = mass_tracer.findFeaturesAndUpdateQscore2D(sd_.getAveragine(), tmp_dspec, ms_level, true);
        deconvolved_features.insert(deconvolved_features.end(), df_decoy.begin(), df_decoy.end());
      }

      Size j = 0;
      for (Size i : element.second)
        deconvolved_spectra[i] = tmp_dspec[j++];
    }
  }
}
} // namespace OpenMS