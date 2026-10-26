// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHQuantAlgorithm.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/ML/NNLS/NonNegativeLeastSquaresSolver.h>
#include <OpenMS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#ifdef _OPENMP
  #include <omp.h>
#endif
#include <chrono>
#include <functional>
#include <memory>
#include <numeric>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <tuple>

namespace OpenMS
{
  FLASHQuantAlgorithm::FLASHQuantAlgorithm() :
      ProgressLogger(),
      DefaultParamHandler("FLASHQuantAlgorithm")
  {
    defaults_.setValue("min_charge", 5, "Minimum charge state to consider");
    defaults_.setValue("max_charge", 50, "Maximum charge state to consider");
    defaults_.setValue("min_mass", 10000, "Minimum mass (Da)");
    defaults_.setValue("max_mass", 70000, "Maximum mass (Da)");
    defaults_.setValue("mz_tol", 10, "Ppm tolerance for m/z values in deconvolution");
    defaults_.setValue("mass_tol", 3, "Mass tolerance in Dalton for integrating similar feature groups into a single one");
    defaults_.setValue("min_isotope_cosine", 0.85, "Cosine threshold between averagine and observed isotope pattern. Note that 0.7 is used for deconvolution", {"advanced"});

    defaults_.setValue("use_smoothed_intensities",
                       "true",
                       "Use smoothed intensities (Savitzky-Golay) instead of raw intensities.",
                       {"advanced"});
    defaults_.setValidStrings("use_smoothed_intensities", {"false", "true"});
    defaults_.setValue("out_shared_details", "false", "Outputs a tsv file including detailed information about the resolved signals (filename = <out_file_name>_shared.tsv", {"advanced"});
    defaults_.setValidStrings("out_shared_details", {"false", "true"});
    defaults_.setValue("resolving_shared_signal", "true", "Resolve shared signals between feature groups (i.e., co-elution)", {"advanced"});
    defaults_.setValidStrings("resolving_shared_signal", {"true", "false"});
    defaultsToParam_();

    this->setLogType(CMD);
  }

  void FLASHQuantAlgorithm::updateMembers_()
  {
    charge_lower_bound_ = (Size) param_.getValue("min_charge");
    charge_upper_bound_ = (Size) param_.getValue("max_charge");

    min_mass_ = (double) param_.getValue("min_mass");
    max_mass_ = (double) param_.getValue("max_mass");

    mz_tolerance_ = (double) param_.getValue("mz_tol");
    mass_tolerance_da_ = (double) param_.getValue("mass_tol");
    min_isotope_cosine_ = (double) param_.getValue("min_isotope_cosine");

    use_smoothed_intensities_ = param_.getValue("use_smoothed_intensities").toBool();
    shared_output_requested_ = param_.getValue("out_shared_details").toBool();
    resolving_shared_signal_ = param_.getValue("resolving_shared_signal").toBool();
  }

  Param FLASHQuantAlgorithm::getFLASHDeconvParams_() const
  {
    Param fd_defaults = SpectralDeconvolution().getDefaults();
    // overwrite algorithm default so we export everything (important for copying back MSstats results)
    fd_defaults.setValue("min_charge", (int) charge_lower_bound_);
    fd_defaults.setValue("max_charge", (int) charge_upper_bound_);
    fd_defaults.setValue("min_mass", min_mass_);
    fd_defaults.setValue("max_mass", max_mass_);
    fd_defaults.setValue("min_cos", DoubleList{.8, .8});
    fd_defaults.setValue("tol", DoubleList{mz_tolerance_, 10.0});
    return fd_defaults;
  }

  void FLASHQuantAlgorithm::run(std::vector<MassTrace> &input_mtraces, std::vector<FeatureGroup>& out_fgs)
  {
    // *********************************************************** //
    // Step 1 deconvolute mass traces
    // *********************************************************** //
    // initialize input & output
    std::vector<FeatureSeed> input_seeds;
    input_seeds.reserve(input_mtraces.size());
    std::vector<MassTrace> updated_masstraces;
    updated_masstraces.reserve(input_mtraces.size());
    Size index = 0;
    double trace_min_length = 100.0;
    for (auto iter = input_mtraces.begin(); iter != input_mtraces.end(); ++iter) {
      if (iter->getSize() < min_nr_peaks_in_mtraces_)
      {
        continue;
      }
      if (iter->getFWHM() == 0)
      {
        continue;
      }
      if (iter->getTraceLength() < trace_min_length)
      {
        trace_min_length = iter->getTraceLength();
      }
      updated_masstraces.push_back(*iter);

      FeatureSeed tmp_seed(*iter, use_smoothed_intensities_);
      tmp_seed.setTraceIndex(index);
      input_seeds.push_back(tmp_seed);
      ++index;
    }
    input_mtraces.swap(updated_masstraces);
    updated_masstraces.clear();
    OPENMS_LOG_INFO << "# final input mass traces : " << input_mtraces.size()  << " (min length " << trace_min_length << " s)" << endl;

    std::vector<FeatureGroup> features;
    features.reserve(input_mtraces.size());

    // run deconvolution
    buildMassTraceGroups_(input_seeds, features);
    features.shrink_to_fit();
    input_seeds.clear();

    // filter out peaks lower than noise level
    filterOutNoisePeaks(features, input_mtraces);

    // *********************************************************** //
    // Step 2 mass artifact removal & post processing...
    // *********************************************************** //
    refineFeatureGroups_(features);
    OPENMS_LOG_INFO << "#Detected feature groups : " << features.size() << endl;

    // *********************************************************** //
    // Step 3 clustering features
    // *********************************************************** //
    if (resolving_shared_signal_)
    {
      if (shared_output_requested_)
      {
        setOptionalDetailedOutput_();
      }

      clusterFeatureGroups_(features, input_mtraces);
      if (shared_output_requested_)
      {
        shared_out_stream_.close();
      }
    }
    OPENMS_LOG_INFO << "#Final feature groups: " << features.size() << endl;

    // output
    setFeatureGroupMembersForResultWriting_(features);
    out_fgs.swap(features);
    features.clear();
  }

   void FLASHQuantAlgorithm::makeMSSpectrum_(std::vector<FeatureSeed *> &local_traces, MSSpectrum &spec, const double &rt) const
  {
    for (auto &tmp_trace : local_traces)
    {
      spec.emplace_back(tmp_trace->getCentroidMz(), tmp_trace->getIntensity());
    }
    spec.setMSLevel(1);
    spec.setName("");
    spec.setRT(rt);
    spec.sortByPosition();
  }

  void FLASHQuantAlgorithm::setFeatureGroupMembersForResultWriting_(std::vector<FeatureGroup> &f_groups) const
  {
    // this cannot be done in FeatureGroup (in FLASHQuantHelper) due to some methods to be used only in here

    for (auto &fgroup : f_groups)
    {
      // initialize
      auto per_charge_int = std::vector<float>(1 + fgroup.getMaxCharge(), .0);
      auto per_charge_cos = std::vector<float>(1 + fgroup.getMaxCharge(), .0);

      std::vector<std::vector<float>> per_cs_isos(1 + fgroup.getMaxCharge(), std::vector<float>(fgroup.getIsotopeIntensities().size(), .0));

      // iterate all FeatureSeeds to collect per feature values
      for(auto& seed : fgroup)
      {
        if (seed.getIsotopeIndex() < 0)
        {
          continue;
        }
        per_charge_int[seed.getCharge()] += seed.getIntensity();
        per_cs_isos[seed.getCharge()][seed.getIsotopeIndex()] += seed.getIntensity();
      }

      for(auto& cs : fgroup.getChargeSet())
      {
        auto &this_cs_isos = per_cs_isos[cs];
        int min_isotope_index = std::distance( this_cs_isos.begin(), std::find_if( this_cs_isos.begin(), this_cs_isos.end(), [](auto &x) { return x != 0; }));
        int max_isotope_index = std::distance( this_cs_isos.begin(), (std::find_if( this_cs_isos.rbegin(), this_cs_isos.rend(), [](auto &x) { return x != 0; }) + 1).base());

        auto iso_dist = iso_model_.get(fgroup.getMonoisotopicMass());
        float cos_score = SpectralDeconvolution::getCosine(per_cs_isos[cs], min_isotope_index, max_isotope_index + 1,
                                                           iso_dist, 0, 2); // min_iso_size_: based on FLASHDeconv
        per_charge_cos[cs] = cos_score;
      }
      // calculate average mass
      auto avg_mass = iso_model_.getAverageMassDelta(fgroup.getMonoisotopicMass()) + fgroup.getMonoisotopicMass();

      // setting values
      fgroup.setPerChargeIntensities(per_charge_int);
      fgroup.setPerChargeCosineScore(per_charge_cos);
      fgroup.setAverageMass(avg_mass);
    }
  }

  void FLASHQuantAlgorithm::getFeatureFromSpectrum_(std::vector<FeatureSeed *> &local_traces,
                                                          std::vector<FeatureGroup> &local_fgroup,
                                                          const double &rt)
  {
     // convert local_traces_ to MSSpectrum
     MSSpectrum spec;
     makeMSSpectrum_(local_traces, spec, rt);

     // run deconvolution
     PeakGroup empty_pg; // empty one, since only MS1s are considered.
     deconv_.performSpectrumDeconvolution(spec, 0, empty_pg);
     DeconvolvedSpectrum& deconv_spec = deconv_.getDeconvolvedSpectrum();

     if (deconv_spec.empty()) // if no result was found
     {
       return;
     }

     // convert deconvolved result into FeatureGroup
     for (auto &deconv : deconv_spec)
     {
       // filter out if deconv results are not sufficient
       if (deconv.size() < min_nr_mtraces_)
       {
         continue;
       }

       FeatureGroup fg(deconv);
       fg.setMaxIsotopeIndex(iso_model_.get(deconv.getMonoMass()).size());

       // Add individual FeatureSeeds to FeatureGroup
       for (auto &peak: deconv)
       {
         // if isotope index of this peak is out of threshold, don't include this
         if (peak.isotopeIndex >= (int) fg.getMaxIsotopeIndex())
         {
           continue;
         }

         // find seed index
         auto it = std::find_if(local_traces.begin(), local_traces.end(),
                                [peak] (FeatureSeed* const& f) -> bool {
                                  return (f->getCentroidMz() == peak.mz); // intensity changes in performSpectrumDeconvolution, thus cannot be used for filtering
                                });
         FeatureSeed tmp_seed(**it);
         tmp_seed.setCharge(peak.abs_charge);
         tmp_seed.setIsotopeIndex(peak.isotopeIndex);
         tmp_seed.setMass(peak.mass);

         fg.push_back(tmp_seed);
       }
       local_fgroup.push_back(fg);
     }
  }

  bool FLASHQuantAlgorithm::doFWHMbordersOverlap_(const std::pair<double, double> &border1,
                                                        const std::pair<double, double> &border2) const
  {
    if ((border1.first > border2.second) || (border2.first > border1.second))
      return false;

    const double overlap_length = std::min(border1.second, border2.second) - std::max(border1.first, border2.first);
    if ((overlap_length / (border1.second - border1.first) < 0.5) &&
        (overlap_length / (border2.second - border2.first) < 0.5))
      return false;

    return true;
  }

  bool FLASHQuantAlgorithm::doMassTraceIndicesOverlap(const FeatureGroup &fg1, const FeatureGroup &fg2,
                                                            const double overlap_percentage_threshold, const bool charge_specific) const
  {
    // get overlapping charge states
    int min_overlapping_charge = std::max(fg1.getMinCharge(), fg2.getMinCharge());
    int max_overlapping_charge = std::min(fg1.getMaxCharge(), fg2.getMaxCharge());

    if (charge_specific && (min_overlapping_charge > max_overlapping_charge)) // no overlapping charge
    {
      return true;
    }

    // collect possible overlapping mass_traces based on charges
    std::vector<Size> mt_indices_1;
    std::vector<Size> mt_indices_2;
    mt_indices_1.reserve(fg1.size());
    mt_indices_2.reserve(fg2.size());
    for (Size fg1_idx = 0; fg1_idx < fg1.size(); ++fg1_idx)
    {
      if (!charge_specific || (fg1[fg1_idx].getCharge() >= min_overlapping_charge && fg1[fg1_idx].getCharge() <= max_overlapping_charge))
      {
        mt_indices_1.push_back(fg1[fg1_idx].getTraceIndex());
      }
    }
    if (mt_indices_1.empty())
    {
      return false;
    }
    for (Size fg2_idx = 0; fg2_idx < fg2.size(); ++fg2_idx)
    {
      if (!charge_specific || (fg2[fg2_idx].getCharge() >= min_overlapping_charge && fg2[fg2_idx].getCharge() <= max_overlapping_charge))
      {
        mt_indices_2.push_back(fg2[fg2_idx].getTraceIndex());
      }
    }
    if (mt_indices_2.empty())
    {
      return false;
    }

    std::sort(mt_indices_1.begin(), mt_indices_1.end());
    std::sort(mt_indices_2.begin(), mt_indices_2.end());

    Size min_vec_size = std::min(mt_indices_1.size(), mt_indices_2.size());
    std::vector<Size> inters_vec;
    inters_vec.reserve(min_vec_size);
    std::set_intersection(mt_indices_1.begin(),
                          mt_indices_1.end(),
                          mt_indices_2.begin(),
                          mt_indices_2.end(),
                          std::back_inserter(inters_vec));

    double overlap_percentage = static_cast<double>(inters_vec.size()) / static_cast<double>(min_vec_size);
    // TODO : change this to overlapping only major cs?
    if (overlap_percentage < overlap_percentage_threshold)
    {
      return false;
    }
    return true;
  }

  bool FLASHQuantAlgorithm::rescoreFeatureGroup_(FeatureGroup &fg) const
  {
    if (fd_seeded_)
    {
      fg.updateMembersForScoring();
      if (!fg.empty()) fg.updateMembers(use_smoothed_intensities_);
      return true;
    }
    if ((!scoreAndFilterFeatureGroup_(fg)))
    {
      return false;
    }
    // update private members in FeatureGroup based on the changed FeatureSeeds
    fg.updateMembers(use_smoothed_intensities_);
    return true;
  }

  bool FLASHQuantAlgorithm::scoreAndFilterFeatureGroup_(FeatureGroup &fg, double min_iso_score) const
  {
    /// return false when scoring is not done (filtered out)
    // update monoisotopic mass, isotope_intensities_ and charge vector to use for filtering
    fg.updateMembersForScoring();

    /// if this FeatureGroup is within the target, pass any filter
    bool isNotTarget = true;
//    if (with_target_masses_ && isThisMassOneOfTargets_(fg.getMonoisotopicMass(), fg.getRtOfMostAbundantMT()))
//    {
//      isNotTarget = false;
//    }

    /// filter if the number of charges are not enough
    if (isNotTarget && (fg.getChargeSet().size() < min_nr_charges_))
    {
      return false;
    }

    /// filter if the mass is out of range
    if (isNotTarget && (fg.getMonoisotopicMass() < min_mass_ || fg.getMonoisotopicMass() > max_mass_))
    {
      return false;
    }

    /// isotope cosine calculation
    /// based on: SpectralDeconvolution::scoreAndFilterPeakGroups_()
    if (min_iso_score == -1)
    {
      min_iso_score = min_isotope_cosine_;
    }
    int offset = 0;
    float isotope_score = SpectralDeconvolution::getIsotopeCosineAndIsoOffset(fg.getMonoisotopicMass(),
                                                                              fg.getIsotopeIntensities(),
                                                                              offset, iso_model_, 1, 0, std::vector<double> {});
    fg.setIsotopeCosine(isotope_score);
    if (isNotTarget && isotope_score < min_iso_score)
    {
      return false;
    }
    // update values based on the calculated offset
    if (offset > 0)
    {
      fg.setMonoisotopicMass(fg.getMonoisotopicMass() + iso_da_distance_ * offset);
      fg.updateIsotopeIndices(offset); // TODO: change?
    }

    return true;
  }

  void FLASHQuantAlgorithm::refineFeatureGroups_(std::vector<FeatureGroup> &in_features)
  {
    // change min, max charges based on built FeatureGroups (for later use in scoring)
    int min_abs_charge = INT_MAX;
    int max_abs_charge = INT_MIN;

    // minimum isotope cosine score for the refining FeatureGroups
    double min_iso_score = 0.5;
    if (!resolving_shared_signal_) // if not resolving shared signal, use the final min_isotope_cosine_ (because this is the last step)
    {
      min_iso_score = min_isotope_cosine_;
    }

    // output features
    std::vector<FeatureGroup> out_feature;
    out_feature.reserve(in_features.size());

    // sort features by masses
    std::sort(in_features.begin(), in_features.end());

    // setting pointer vector (to reduce time complexity for push_back & erase)
    std::vector<FeatureGroup*> in_feature_ptrs;
    in_feature_ptrs.reserve(in_features.size());

    // set variables according to the detected features
    for (auto &f: in_features)
    {
      in_feature_ptrs.push_back(&f);

      min_abs_charge =
          min_abs_charge < f.getMinCharge() ? min_abs_charge : f.getMinCharge();
      max_abs_charge =
          max_abs_charge > f.getMaxCharge() ? max_abs_charge : f.getMaxCharge();
    }
    charge_lower_bound_ = min_abs_charge;
    charge_upper_bound_ = max_abs_charge;

    Size initial_size = in_features.size();
    this->startProgress(0, initial_size, "refining feature groups");
    // insert FeatureGroup with the highest score to out_features, merge if other FeatureGroup exist within mass_tol
    while (!in_feature_ptrs.empty())
    {
      this->setProgress(initial_size - in_feature_ptrs.size());

      // get a feature with the highest Intensity
      auto candidate_fg_pointer = std::max_element(in_feature_ptrs.begin(), in_feature_ptrs.end(), FLASHQuantHelper::CmpFeatureGroupByIntensity());
      auto candidate_fg = *candidate_fg_pointer;

      // get all features within mass_tol from candidate FeatureGroup
      std::vector<FeatureGroup*>::iterator low_it, up_it;

      // open up the search space (10 Da) to check the mass trace overlap. mass_tolerance_da_ will be checked later
      FeatureGroup lower_fg(candidate_fg->getMonoisotopicMass() - 10);
      FeatureGroup upper_fg(candidate_fg->getMonoisotopicMass() + 10);

      low_it = std::lower_bound(in_feature_ptrs.begin(), in_feature_ptrs.end(), &lower_fg, FLASHQuantHelper::CmpFeatureGroupPointersByMass());
      up_it = std::upper_bound(in_feature_ptrs.begin(), in_feature_ptrs.end(), &upper_fg, FLASHQuantHelper::CmpFeatureGroupPointersByMass());

      // no matching in features (found only candidate itself)
      if (up_it - low_it == 1)
      {
        // save it to out_features
        if (scoreAndFilterFeatureGroup_(*candidate_fg, min_iso_score)) // rescoring is needed to set scores in FeatureGroup
        {
          candidate_fg->updateMembers(use_smoothed_intensities_);
          out_feature.push_back(*candidate_fg);
        }

        // remove candidate from features
        in_feature_ptrs.erase(candidate_fg_pointer);
        continue;
      }

      // check if found features are overlapping with the candidate feature
      std::vector<int> v_indices_to_remove;
      v_indices_to_remove.reserve(up_it - low_it);
      std::set<Size> mt_indices_to_add;
      std::vector<FeatureSeed *> mts_to_add; // One unique mt can be included in different FGs -> different FeatureSeed
      mts_to_add.reserve((up_it - low_it) * candidate_fg->size());

      for (; low_it != up_it; ++low_it)
      {
        // if low_it is candidate feature, ignore
        if (candidate_fg_pointer == low_it)
        {
          v_indices_to_remove.push_back(low_it - in_feature_ptrs.begin());
          continue;
        }

        // if the mass difference is larger than mass_tolerance_da_, check the mass trace overlap
        if ((std::abs(candidate_fg->getMonoisotopicMass() - (*low_it)->getMonoisotopicMass()) > mass_tolerance_da_)
            && !doMassTraceIndicesOverlap(**low_it, *candidate_fg))
        {
          continue;
        }

        // check if fwhm overlaps
        if (!doFWHMbordersOverlap_((*low_it)->getFwhmRange(), candidate_fg->getFwhmRange()))
        {
          continue;
        }

        // merge found feature to candidate feature
        auto &trace_indices = candidate_fg->getTraceIndices();
        for (auto &new_mt: **low_it)
        {
          // if this mass trace is not used in candidate_fg
          if (std::find(trace_indices.begin(), trace_indices.end(), new_mt.getTraceIndex()) == trace_indices.end())
          {
            mt_indices_to_add.insert(new_mt.getTraceIndex());
            mts_to_add.push_back(&new_mt);
          }
        }
        // add index of found feature to "to_be_removed_vector"
        v_indices_to_remove.push_back(low_it - in_feature_ptrs.begin());
      }

      // sort mts_to_add by abundance
      std::sort(mts_to_add.begin(), mts_to_add.end(), FLASHQuantHelper::CmpFeatureSeedByIntensity());

      // add extra masstraces to candidate_feature
      FeatureGroup final_candidate_fg(*candidate_fg); // copy of candidate_feature
      double mono_mass = final_candidate_fg.getMonoisotopicMass();
      int max_isotope = (int)iso_model_.getLastIndex(mono_mass);
      int min_isotope = (int)(iso_model_.getApexIndex(mono_mass) - iso_model_.getLeftCountFromApex(mono_mass));
      FeatureSeed* apex_of_this_fg = final_candidate_fg.getApexLMT();
      double rt_of_apex = (apex_of_this_fg->getFwhmEnd() + apex_of_this_fg->getFwhmStart())/2;
      for (auto &new_mt: mts_to_add)
      {
        // to skip duplicated masstraces that are included
        if (mt_indices_to_add.find(new_mt->getTraceIndex()) == mt_indices_to_add.end())
        {
          continue;
        }
        mt_indices_to_add.erase(new_mt->getTraceIndex());

        // if this mt is introducing new charge
        if (!final_candidate_fg.doesThisChargeExist(new_mt->getCharge()))
        {
          final_candidate_fg.push_back(*new_mt);
          continue;
        }

        /// re-calculate isotope index (from PeakGroup::recruitAllPeaksInSpectrum)
        double cmz = (mono_mass) / new_mt->getCharge() + Constants::PROTON_MASS_U; // mono mz
        const double iso_delta = iso_da_distance_ / (double)new_mt->getCharge();
        int iso_index = (int)round((new_mt->getCentroidMz() - cmz) / iso_delta);
        if (iso_index > max_isotope)
        {
          break;
        }
        if (iso_index < min_isotope)
        {
          continue;
        }

        new_mt->setIsotopeIndex(iso_index);

        // checking if this mt is colliding with the existing mt (iso index)
        bool introduce_new_mt = true;
        for (Size i = 0; i < final_candidate_fg.size(); ++i)
        {
          const FeatureSeed& seed = final_candidate_fg[i];
          if (seed.getCharge() != new_mt->getCharge())
          {
            continue;
          }
          if (seed.getIsotopeIndex() == iso_index) // if this isotope exist in feature group,
          {
            // the one with closer apex retention time to the max peak of feature group wins.
            double diff_to_new_mt = abs((new_mt->getFwhmEnd()+new_mt->getFwhmStart())/2 - rt_of_apex);
            double diff_to_org_mt = abs((seed.getFwhmEnd()+seed.getFwhmStart())/2 - rt_of_apex);
            if (diff_to_new_mt >= diff_to_org_mt)
            {
              introduce_new_mt = false;
            }
            else
            {
              final_candidate_fg.erase(final_candidate_fg.begin() + i);
            }
            break;
          }
        }
        if (introduce_new_mt)
        {
          final_candidate_fg.push_back(*new_mt);
        }
      }

      // don't merge when it failed to exceed filtering threshold
      if (!scoreAndFilterFeatureGroup_(final_candidate_fg, min_iso_score))
      {
        if (scoreAndFilterFeatureGroup_(*candidate_fg, min_iso_score))
        {
          candidate_fg->updateMembers(use_smoothed_intensities_);
          out_feature.push_back(*candidate_fg);
        }
      }
      else // if to be merged, save the updated one to out_feature
      {
        final_candidate_fg.updateMembers(use_smoothed_intensities_);
        out_feature.push_back(final_candidate_fg);
      }

      // remove candidate from features
      std::sort(v_indices_to_remove.begin(), v_indices_to_remove.end(), std::greater<>());

      for (auto &idx : v_indices_to_remove)
      {
        in_feature_ptrs.erase(in_feature_ptrs.begin() + idx);
      }
    }
    this->endProgress();

    in_features.swap(out_feature);
  }

  void FLASHQuantAlgorithm::buildMassTraceGroups_(std::vector<FeatureSeed> &mtraces, std::vector<FeatureGroup> &features)
  {
    /// FLASHDeconvAlgorithm setting
    Param fd_defaults = getFLASHDeconvParams_();
    deconv_.setParameters(fd_defaults);
    deconv_.calculateAveragine(false);
    iso_model_ = deconv_.getAveragine();
//    std::vector<double> target_masses_; // monoisotope
//    fd_.setTargetMasses(target_masses_, ms_level);

    /// group mass traces to spectrum
    std::vector<std::pair<double, FeatureSeed*>> mt_rt_starts;
    std::vector<std::pair<double, FeatureSeed*>> mt_rt_ends;
    mt_rt_starts.reserve(mtraces.size());
    mt_rt_ends.reserve(mtraces.size());
    int counter = 0;

    // sort input mass traces in RT
    std::sort(mtraces.begin(), mtraces.end(), FLASHQuantHelper::CmpFeatureSeedByRT());

    // collect rt information from mtraces to generate spectrum
    double rt_binning_size = numeric_limits<double>::max();
    for (auto &trace: mtraces)
    {
      mt_rt_starts.push_back(std::make_pair(trace.getFwhmStart(), &trace));
      mt_rt_ends.push_back(std::make_pair(trace.getFwhmEnd(), &trace));
      if (trace.getMassTrace().getFWHM() < rt_binning_size)
      {
        rt_binning_size = trace.getMassTrace().getFWHM();
      }
    }
    // minimum bin size is 1
    if (rt_binning_size < 1)
    {
      rt_binning_size = 1;
    }

    // sorting mass traces in rt
    std::sort(mt_rt_starts.begin(), mt_rt_starts.end());
    std::sort(mt_rt_ends.begin(), mt_rt_ends.end());

    std::vector<std::pair<double, FeatureSeed *>>::const_iterator rt_s_iter = mt_rt_starts.begin();
    std::vector<std::pair<double, FeatureSeed *>>::const_iterator rt_e_iter = mt_rt_ends.begin();
    auto end_of_iter = mt_rt_starts.end();
//    double end_of_current_rt_window = mt_rt_starts[0].first;
    double end_of_current_rt_window = rts_from_org_scans_[0];
    std::vector<double>::const_iterator scan_rt_iter = rts_from_org_scans_.begin();
    double last_rt = mt_rt_ends[mt_rt_ends.size() - 1].first;

    // mass traces to be added in a spectrum
    std::vector<FeatureSeed *> local_traces;
    local_traces.reserve(mtraces.size());

//    int possible_spec_size = int((mt_rt_starts[mt_rt_starts.size() - 1].first - end_of_current_rt_window) / rt_binning_size);
    int possible_spec_size = rts_from_org_scans_.size();
    this->startProgress(0, possible_spec_size, "assembling mass traces to features");

    while (rt_s_iter != end_of_iter && end_of_current_rt_window < last_rt)
    {
      this->setProgress(counter);

      // initial rt binning is 1 sec (for generating spectrum)
      end_of_current_rt_window = *(++scan_rt_iter);
//      end_of_current_rt_window += rt_binning_size;

      // add mass traces within rt range
      bool is_new_mt_added = false;
      for (; rt_s_iter != end_of_iter && rt_s_iter->first <= end_of_current_rt_window; ++rt_s_iter)
      {
        local_traces.push_back(rt_s_iter->second);
        is_new_mt_added = true;
      }

      // if nothing is added, increase current_rt
      if (!is_new_mt_added)
      {
        continue;
      }

      // remove mass traces out of rt range
      for (; rt_e_iter != mt_rt_ends.end() && rt_e_iter->first < end_of_current_rt_window; ++rt_e_iter)
      {
        local_traces.erase(std::remove_if(local_traces.begin(), local_traces.end(),
                                          [&rt_e_iter](auto const &p) { return rt_e_iter->second == p; }));
      }
      if (local_traces.empty())
      {
        continue;
      }

      // sort local traces in mz
      sort(local_traces.begin(), local_traces.end(), FLASHQuantHelper::CmpFeatureSeedByMZ());

      std::vector<FeatureGroup> local_fgroup;
      getFeatureFromSpectrum_(local_traces, local_fgroup, end_of_current_rt_window);
      ++counter; // to track the number of generated spectra

      // no feature has been detected
      if (local_fgroup.size() == 0)
      {
        continue;
      }

      for (auto &tmp_fg: local_fgroup)
      {
        sort(tmp_fg.begin(), tmp_fg.end());

        tmp_fg.updateMembersForScoring();
        tmp_fg.updateMembers(use_smoothed_intensities_);
      }

      features.insert(features.end(), local_fgroup.begin(), local_fgroup.end());
    }

    this->endProgress();
    OPENMS_LOG_INFO << "# generated spec from mass traces : " << counter << endl;
    OPENMS_LOG_INFO << "# generated feature groups from mass traces : " << features.size() << endl;
  }

  /// cluster FeatureGroups with shared mass traces, and resolve the shared ones. If not, report as output
  void FLASHQuantAlgorithm::clusterFeatureGroups_(std::vector<FeatureGroup> &fgroups,
                                                  std::vector<MassTrace> &input_mtraces)
  {
    // *********************************************************** //
    // Step 1 preparation for hypergraph : collect feature idx with shared mass traces
    // *********************************************************** //
    std::vector<std::vector<Size>> shared_m_traces(input_mtraces.size(), std::vector<Size>());
    for (Size fg_index = 0; fg_index < fgroups.size(); ++fg_index)
    {
      for (auto &mt_i: fgroups[fg_index].getTraceIndices())
      {
        shared_m_traces[mt_i].push_back(fg_index);
      }
    }

    // *********************************************************** //
    // Step 2 constructing hypergraph from featuregroups
    //        node = mass traces
    //        hyperedge = feature groups
    // *********************************************************** //
    Size num_nodes = shared_m_traces.size();
    std::vector<bool> bfs_visited;
    bfs_visited.resize(num_nodes, false);
    std::queue<Size> bfs_queue;
    Size search_pos = 0; // keeping track of mass trace index to look for seed

    std::vector<FeatureGroup> out_features;
    out_features.reserve(fgroups.size());
    std::vector<std::set<Size>> clusters; // FLASHDeconv masses: the clusters, resolved after the search

    // BFS
    while (true)
    {
      // finding a seed 'shared_mass_trace' to start with (for constructing a cluster)
      bool finished = true;
      for (Size i = search_pos; i < num_nodes; ++i)
      {
        if (!bfs_visited[i])
        {
          // check if this mass_trace is used to any FeatureGroup
          if (shared_m_traces[i].size() == 0)
          {
            bfs_visited[i] = true;
            continue;
          }

          bfs_queue.push(i);
          bfs_visited[i] = true;
          finished = false;
          search_pos = i + 1;
          break;
        }
      }
      if (finished) // if no possible seed is left
        break;

      std::set<Size> fg_indices_in_current_cluster;

      while (!bfs_queue.empty())
      {
        Size i = bfs_queue.front(); // index of seed
        bfs_queue.pop();

        // get FeatureGroup indices sharing this seed
        for (vector<Size>::const_iterator it = shared_m_traces[i].begin();
             it != shared_m_traces[i].end();
             ++it)
        {
          // if this FeatureGroup was visited before
          if (fg_indices_in_current_cluster.find(*it) != fg_indices_in_current_cluster.end())
          {
            continue;
          }

          fg_indices_in_current_cluster.insert(*it);

          FeatureGroup &current_fg = fgroups[*it];
          for (const auto &mt_index: current_fg.getTraceIndices())
          {
            if (!bfs_visited[mt_index])
            {
              bfs_queue.push(mt_index);
              bfs_visited[mt_index] = true;
            }
          }
        }
      }

      if (fd_seeded_) // resolved below, cluster by cluster in parallel
      {
        clusters.push_back(std::move(fg_indices_in_current_cluster));
        continue;
      }

      // this feature is not sharing any mass traces with others
      if (fg_indices_in_current_cluster.size() == 1)
      {
        // re-scoring (score threshold changed from this method)
        auto &out_feat = fgroups[*(fg_indices_in_current_cluster.begin())];
        if (rescoreFeatureGroup_(out_feat))
        {
          out_features.push_back(out_feat);
        }
        continue;
      }

      // resolve the conflict among feature groups
      resolveConflictInCluster_(fgroups, input_mtraces, shared_m_traces, fg_indices_in_current_cluster, out_features);
    }

    if (fd_seeded_)
    {
      // A cluster's groups and traces belong to no other cluster, so the clusters are resolved in parallel (the largest first), and their
      // groups are collected in cluster order: the same groups as one after the other.
      std::vector<std::vector<FeatureGroup>> cluster_out(clusters.size());
      std::vector<Size> order(clusters.size());
      std::iota(order.begin(), order.end(), 0);
      std::stable_sort(order.begin(), order.end(), [&clusters](Size a, Size b) { return clusters[a].size() > clusters[b].size(); });
      std::vector<std::exception_ptr> errors(clusters.size());
      const auto resolve_cluster = [&](const Size c) {
        try
        {
          auto& cluster = clusters[c];
          if (cluster.size() == 1)
          {
            auto& out_feat = fgroups[*(cluster.begin())];
            if (rescoreFeatureGroup_(out_feat)) cluster_out[c].push_back(out_feat);
            return;
          }
          resolveConflictInCluster_(fgroups, input_mtraces, shared_m_traces, cluster, cluster_out[c]);
        }
        catch (...)
        {
          errors[c] = std::current_exception();
        }
      };
      // the large clusters one after the other, each with its regions' fits and traces in parallel; then the others in parallel
      const Size large = 500;
      Size first_small = 0;
      while (first_small < order.size() && clusters[order[first_small]].size() > large)
      {
        resolve_cluster(order[first_small++]);
      }
#pragma omp parallel for schedule(dynamic, 1)
      for (SignedSize o = (SignedSize)first_small; o < (SignedSize)order.size(); ++o)
      {
        resolve_cluster(order[o]);
      }
      for (const auto& error : errors)
      {
        if (error) std::rethrow_exception(error);
      }
      for (auto& groups_of_cluster : cluster_out)
      {
        out_features.insert(out_features.end(), groups_of_cluster.begin(), groups_of_cluster.end());
      }
      // the groups set aside (in the order of their FLASHDeconv groups, as the clusters ran in parallel)
      std::stable_sort(fd_set_aside_.begin(), fd_set_aside_.end(), [](const FeatureGroup& a, const FeatureGroup& b) { return a.getFDGroup() < b.getFDGroup(); });
      out_features.insert(out_features.end(), fd_set_aside_.begin(), fd_set_aside_.end());
      fd_set_aside_.clear();
    }
    out_features.shrink_to_fit();
    std::swap(out_features, fgroups);
  }

  void FLASHQuantAlgorithm::filterOutNoisePeaks(std::vector<FeatureGroup> &fgroups, std::vector<MassTrace> &mtraces) const
  {
    /// 1. collect excluded traces from fgroups (consider them as noises)
    std::vector<bool> traces_used(mtraces.size(), false);
    Size max_num_of_peaks_in_seed = 0;
    for (const auto &fgroup : fgroups)
    {
      for (const auto &seed : fgroup)
      {
        if (!traces_used[seed.getTraceIndex()])
        { // when this trace has not been visited
          Size num_of_peaks = seed.getMassTrace().getSize();
          if (num_of_peaks > max_num_of_peaks_in_seed)
          {
            max_num_of_peaks_in_seed = num_of_peaks;
          }
        }
        traces_used[seed.getTraceIndex()] = true;
      }
    }

    /// 2. get median peak intensity of noises
    std::vector<double> intensities;
    int num_of_used_traces = std::count(traces_used.begin(), traces_used.end(), true);
    // the peaks of the unused traces (was num_of_unused * max peaks of a seed: gigabytes on a lysate run)
    Size unused_peaks = 0;
    for (Size i = 0; i < mtraces.size(); ++i)
    {
      if (!traces_used[i]) unused_peaks += mtraces[i].getSize();
    }
    intensities.reserve(unused_peaks);
    (void) num_of_used_traces;
    Size trace_size = mtraces.size();
    for (Size i = 0; i < trace_size; ++i)
    {
      if (traces_used[i]) continue;

      if (use_smoothed_intensities_)
      {
        for (const auto &p : mtraces[i].getSmoothedIntensities())
        {
          intensities.push_back(p);
        }
      }else
      {
        for (const auto &p : mtraces[i])
        {
          intensities.push_back(p.getIntensity());
        }
      }
    }
    intensities.shrink_to_fit();
    std::sort(intensities.begin(), intensities.end());

    // median intensity of unused mass traces for deconvolution
    double noise_intensity = intensities.empty() ? .0 : intensities[std::floor(intensities.size()/2)];
    intensities.clear();

    /// 3. update features & input_mtraces by reducing peak intensities based on noise intensity
    Size num_of_traces = traces_used.size();
    for (Size trace_idx = 0; trace_idx < num_of_traces; ++trace_idx)
    {
      if (!traces_used[trace_idx])
      {
        continue;
      }
      auto &trace = mtraces[trace_idx];
      for (auto &p : trace)
      {
        p.setIntensity(std::max(.0, p.getIntensity() - noise_intensity));
      }
      if (use_smoothed_intensities_)
      {
        std::vector<double> updated_smoothed_intys;
        updated_smoothed_intys.reserve(trace.getSize());
        for (auto &p : trace.getSmoothedIntensities())
        {
          updated_smoothed_intys.push_back(std::max(.0, p - noise_intensity));
        }
        trace.setSmoothedIntensities(updated_smoothed_intys);
      }
      trace.estimateFWHM(use_smoothed_intensities_);
    }

    std::vector<FeatureGroup> updated_fgs;
    updated_fgs.reserve(fgroups.size());
    for (auto &this_fg : fgroups)
    {
      auto seed_iter = this_fg.begin();
      while (seed_iter != this_fg.end())
      {
        auto &updated_trace = mtraces[seed_iter->getTraceIndex()];
        Size number_of_nonzeros = 0;
        for (auto & p : updated_trace)
        {
          number_of_nonzeros += p.getIntensity() > 0? 1 : 0;
        }

        if (number_of_nonzeros < min_nr_peaks_in_mtraces_ ||
            updated_trace.getFWHM() == 0 ||
            updated_trace.getIntensity(use_smoothed_intensities_) == 0)
        {
          seed_iter = this_fg.erase(seed_iter);
        }
        else
        {
          seed_iter->setMassTrace(updated_trace);
          seed_iter->setIntensity(updated_trace.getIntensity(use_smoothed_intensities_));
          ++seed_iter;
        }
      }
      if (this_fg.size() < min_nr_mtraces_ && !fd_seeded_)
      {
        continue;
      }
      this_fg.updateMembersForScoring();
      if (!this_fg.empty()) this_fg.updateMembers(use_smoothed_intensities_);
      updated_fgs.push_back(this_fg);
    }
    fgroups.swap(updated_fgs);
    updated_fgs.clear();
  }

  void FLASHQuantAlgorithm::setOptionalDetailedOutput_()
  {
    String out_path = output_file_path_.substr(0, output_file_path_.find_last_of(".")) + "_shared.tsv";
    shared_out_stream_.open(out_path, std::fstream::out);
    shared_out_stream_ << "FeatureGroupID\tTraceType\tMass\tCharge\tIsotopeIndex\tQuantValue\tCentroidMz\tTraceLabel\tRTs\tMZs\tIntensities\n"; // header
    // TraceType = 0 (not shared), 1 (true, before resolution), 2(true, after resolution), 3(theoretical shape)
  }

  void FLASHQuantAlgorithm::writeMassTracesOfFeatureGroup_(const FeatureGroup &fgroup,
                                                           const Size &fgroup_idx,
                                                           const std::vector<std::vector<Size>> &shared_m_traces_indices,
                                                           const bool &is_before_resolution)
  {
    /// shared_type = 0 (before resolution - shared), 1(unique), 2 (after resolution - unique & shared)
    auto mt_idxs = fgroup.getTraceIndices();
    Size shared_tag_for_output = is_before_resolution? 0 : 2;

    for (const auto &mt: fgroup)
    {
      // check if current mt is shared with other features
      if (is_before_resolution) // before resolution
      {
        shared_tag_for_output = (shared_m_traces_indices[mt.getTraceIndex()].size() == 1)? 0 : 1; // shared mass trace = 0; unique mass trace = 1;
      }

      stringstream rts;
      stringstream mzs;
      stringstream intys;

      for (const auto &peak: mt.getMassTrace())
      {
        mzs << std::to_string(peak.getMZ()) << ",";
        rts << std::to_string(peak.getRT()) << ",";
        intys << std::to_string(peak.getIntensity()) << ",";
      }
      std::string peaks = rts.str();
      peaks.pop_back();
      peaks = peaks + "\t" + mzs.str();
      peaks.pop_back();
      peaks = peaks + "\t" + intys.str();
      peaks.pop_back();

      shared_out_stream_ << fgroup_idx << "\t"
                         << shared_tag_for_output << "\t"
                         << std::to_string(fgroup.getMonoisotopicMass()) << "\t"
                         << mt.getCharge() << "\t"
                         << mt.getIsotopeIndex() << "\t"
                         << std::to_string(mt.getIntensity()) << "\t"
                         << std::to_string(mt.getCentroidMz()) << "\t"
                         << mt.getMassTrace().getLabel() << "\t"
                         << peaks + "\n";
    }
  }

  void FLASHQuantAlgorithm::writeTheoreticalShapeForConflictResolution_(const Size &fgroup_idx,
                                                                              const FeatureSeed &shared_mt,
                                                                              const std::vector<double> &theo_intensities,
                                                                              const double &calculated_ratio)
  {
    stringstream rts;
    stringstream mzs;
    stringstream intys;

    for (const auto &peak: shared_mt.getMassTrace())
    {
      mzs << std::to_string(peak.getMZ()) << ",";
      rts << std::to_string(peak.getRT()) << ",";
    }
    for (const auto &inty : theo_intensities)
    {
      intys << std::to_string(inty) << ",";
    }

    std::string peaks = rts.str();
    peaks.pop_back();
    peaks = peaks + "\t" + mzs.str();
    peaks.pop_back();
    peaks = peaks + "\t" + intys.str();
    peaks.pop_back();

    shared_out_stream_ << fgroup_idx << "\t"
                       << 3 << "\t" // shared tag = 3
                       << std::to_string(shared_mt.getMass()) << "\t"
                       << shared_mt.getCharge() << "\t"
                       << shared_mt.getIsotopeIndex() << "\t"
                       << std::to_string(calculated_ratio) << "\t" // ratio
                       << std::to_string(shared_mt.getCentroidMz()) << "\t"
                       << shared_mt.getMassTrace().getLabel() << "\t"
                       << peaks + "\n";
  }

  void FLASHQuantAlgorithm::filterOutIneligibleFeatureGroupsInCluster(std::vector<FeatureGroup>& feature_groups,
                                                                      std::vector<std::vector<Size>>& shared_m_traces_indices,
                                                                      std::set<Size>& candidate_fg_indices) const
  {
    /// Ineligible feature groups = harmonics OR the ones without unique feature groups
    std::set<Size> fg_idx_to_remove;
    // the groups of the cluster by mass: a group is compared only with the later ones whose mass is near one of its harmonic masses (the
    // others pass the harmonic test without effect), in the same order as before
    std::vector<std::pair<double, Size>> by_mass;
    by_mass.reserve(candidate_fg_indices.size());
    for (const Size idx : candidate_fg_indices)
    {
      by_mass.emplace_back(feature_groups[idx].getMonoisotopicMass(), idx);
    }
    std::sort(by_mass.begin(), by_mass.end());
    for (auto fg_idx = candidate_fg_indices.begin(); fg_idx != candidate_fg_indices.end(); ++fg_idx)
    {
      if (fg_idx_to_remove.count(*fg_idx) > 0) // if this is already removed
      {
        continue;
      }

      auto& fgroup = feature_groups[*fg_idx];

      // check if any masstrace is unique
      bool is_this_all_sharing = true;
      for (auto& seed_iter : fgroup)
      {
        auto trace_idx = seed_iter.getTraceIndex();
        if (shared_m_traces_indices[trace_idx].size() == 1) // this masstrace is unique
        {
          is_this_all_sharing = false;
          break;
        }
      }
      if (is_this_all_sharing)
      {
        fg_idx_to_remove.insert(*fg_idx);
        continue;
      }

      // check if this is harmonics to the other
      std::vector<double> harmonic_masses; // harmonic mass candidates
      const double fg_mass = fgroup.getMonoisotopicMass();
      for (int harmonic_ratio = 2; harmonic_ratio < 8; ++harmonic_ratio)
      {
        for (int iso = -2; iso < 3; ++iso)
        {
          // high harmonics
          double h_mass = (fg_mass + iso*iso_da_distance_) / harmonic_ratio;
          if (h_mass > min_mass_)
          {
            harmonic_masses.push_back(h_mass);
          }
          // low harmonic
          h_mass = (fg_mass + iso*iso_da_distance_) * harmonic_ratio;
          if (h_mass < max_mass_)
          {
            harmonic_masses.push_back(h_mass);
          }
        }
      }
      std::sort(harmonic_masses.begin(), harmonic_masses.end());

      // the later groups whose mass m can be within m * 2e-5 of a harmonic mass h, i.e. h / (1 + 2e-5) <= m <= h / (1 - 2e-5) (with margin)
      std::vector<Size> to_compare;
      for (const double h : harmonic_masses)
      {
        auto it = std::lower_bound(by_mass.begin(), by_mass.end(), std::make_pair(h / (1 + 2.1e-5), (Size)0));
        for (; it != by_mass.end() && it->first <= h / (1 - 2.1e-5); ++it)
        {
          if (it->second > *fg_idx) to_compare.push_back(it->second);
        }
      }
      std::sort(to_compare.begin(), to_compare.end());
      to_compare.erase(std::unique(to_compare.begin(), to_compare.end()), to_compare.end());

      for (const Size compare_idx : to_compare)
      {
        const auto fg_idx_to_compare = &compare_idx; // as the iterator over the cluster before
        if (fg_idx_to_remove.count(*fg_idx_to_compare) > 0) // if this is already removed
        {
          continue;
        }

        auto& fgroup_to_compare = feature_groups[*fg_idx_to_compare];
        const double mass_tolerance = fgroup_to_compare.getMonoisotopicMass() * 2 *  1e-5; // 20 ppm

        // is harmonic or not
        auto low_it = std::lower_bound(harmonic_masses.begin(), harmonic_masses.end(), fgroup_to_compare.getMonoisotopicMass() - mass_tolerance); // inclusive
        auto up_it = std::upper_bound(harmonic_masses.begin(), harmonic_masses.end(), fgroup_to_compare.getMonoisotopicMass() + mass_tolerance); // exclusive

        if (low_it == up_it) // if nothing is found, continue
        {
          continue;
        }

        Size harmonic_fg_idx = fgroup.getIntensity() < fgroup_to_compare.getIntensity() ? *fg_idx : *fg_idx_to_compare;
        if (doMassTraceIndicesOverlap(fgroup, fgroup_to_compare, 0.5, false))
        {
          fg_idx_to_remove.insert(harmonic_fg_idx);
        }
        else
        {
          // how many mass traces are shared?
          Size num_of_shared_mt = 0;
          auto &fg_to_remove = feature_groups[harmonic_fg_idx];
          for (auto &mt : fg_to_remove)
          {
            auto trace_idx = mt.getTraceIndex();
            if (shared_m_traces_indices[trace_idx].size() > 1) // this masstrace is shared
            {
              num_of_shared_mt++;
            }
          }

          // if the number of unique mt are less than the range of charges, remove this feature group
          if ((fg_to_remove.size()-num_of_shared_mt) < fg_to_remove.getChargeSet().size())
          {
            fg_idx_to_remove.insert(harmonic_fg_idx);
          }
        }
      }
    }

    // remove feature groups for further work & update the mass trace link
    for (auto &idx : fg_idx_to_remove)
    {
      if (fd_seeded_)
      {
        // the group keeps the traces it does not share; the shared ones go to the other groups
        FeatureGroup kept = feature_groups[idx];
        kept.clear();
        for (const auto& seed : feature_groups[idx])
        {
          if (shared_m_traces_indices[seed.getTraceIndex()].size() <= 1) kept.push_back(seed);
        }
        kept.updateMembersForScoring();
        if (!kept.empty()) kept.updateMembers(use_smoothed_intensities_);
#pragma omp critical (FLASHQuant_set_aside)
        const_cast<FLASHQuantAlgorithm*>(this)->fd_set_aside_.push_back(kept);
      }
      // remove feature groups for further work
      candidate_fg_indices.erase(idx);

      // update the mass link
      for (auto& seed_iter : feature_groups[idx])
      {
        auto trace_idx = seed_iter.getTraceIndex();
        if (shared_m_traces_indices[trace_idx].size() > 1) // update the link
        {
          auto &link = shared_m_traces_indices[trace_idx];
          link.erase(std::find(link.begin(), link.end(), idx));
        }
      }
      feature_groups[idx].clear(); // remove
    }
  }

  void FLASHQuantAlgorithm::resolveConflictInCluster_(std::vector<FeatureGroup> &feature_groups,
                                                      std::vector<MassTrace> &input_masstraces,
                                                      std::vector<std::vector<Size>> &shared_m_traces_indices,
                                                      std::set<Size> &fg_indices_in_this_cluster,
                                                      std::vector<FeatureGroup> &out_featuregroups)
  {
    /// conflict resolution is done in feature level (not feature group level), starting from the most abundant shared signal
    /// input masstraces needs to be changed, if any resolution is done per FeatureGroup -> add new input_masstrace with modified intensities, and change the pointer to it.

    /// remove FeatureGroups without unique masstraces (+ harmonics)
    filterOutIneligibleFeatureGroupsInCluster(feature_groups, shared_m_traces_indices, fg_indices_in_this_cluster);
    if (fg_indices_in_this_cluster.size() == 1)
    {
      auto &fg = feature_groups[*(fg_indices_in_this_cluster.begin())];
      if (rescoreFeatureGroup_(fg))
      {
        out_featuregroups.push_back(fg);
      }
      return;
    }

    /// 1. find shared MTs and corresponding "Features" (not FeatureGroup)
    std::map<Size, std::vector<std::pair<Size, int>>> shared_mt_N_features; // key: mt_idx, value: pair<feature_group_idx, charge>
    std::vector<std::pair<Size, double>> shared_mt_idx_inty_pairs; // to sort shared MTs by their abundance; first: mt_idx, second: mt intensity
    for (auto it = fg_indices_in_this_cluster.begin(); it != fg_indices_in_this_cluster.end(); )
    {
      auto &fgroup = feature_groups[*it];
      bool fg_with_only_unique_traces = true;
      for (auto &seed_iter : fgroup)
      {
        auto trace_idx = seed_iter.getTraceIndex();
        // collect only shared one. Skip unique traces
        if (shared_m_traces_indices[trace_idx].size() < 2)
        {
          continue;
        }

        if (shared_mt_N_features.find(trace_idx) == shared_mt_N_features.end()) // if this seed is not added to the map
        {
          shared_mt_N_features.emplace(trace_idx, std::vector<std::pair<Size, int>>());
          shared_mt_idx_inty_pairs.push_back(std::make_pair(trace_idx, seed_iter.getIntensity()));
        }
        shared_mt_N_features[trace_idx].push_back(std::make_pair(*it, seed_iter.getCharge()));
        fg_with_only_unique_traces = false;
      }

      // double check if any feature group is without shared traces after filtering
      if (fg_with_only_unique_traces)
      {
        if (rescoreFeatureGroup_(fgroup))
        {
          out_featuregroups.push_back(fgroup);
        }
        it = fg_indices_in_this_cluster.erase(it);
        continue;
      }

      if(shared_output_requested_)
      {
        writeMassTracesOfFeatureGroup_(fgroup, *it, shared_m_traces_indices, true);
      }
      ++it;
    }
    if (shared_mt_idx_inty_pairs.empty())
    {
      return;
    }

    /// 2. sort shared MT by intensity (descending)
    std::sort(shared_mt_idx_inty_pairs.begin(), shared_mt_idx_inty_pairs.end(), [](const auto &left, const auto &right) {
      return left.second > right.second;
    });

    std::map<std::vector<Size>, std::shared_ptr<EGHTraceFitter>> fit_cache; // the elution models fitted in this cluster
    /// 3. start from the highest intensity, iterate shared MT and get conflicting region (collecting conflicting Features)
    std::multimap<Size, Feature> solved_features; // key: feature group index, value: solved Feature
    // The shared traces by intensity, each starting a region unless an earlier region took it (it was erased from the list, from its
    // front and per trace of a region: quadratic in the cluster size), and the traces of a region visited in the order they are found
    // (it was a scan of the hash map for the first unvisited one per step: quadratic in the region size; Eclip_1480 has a region of
    // ~50k features). The regions and their features are the same; only the order the features are collected in may differ.
    std::unordered_set<Size> traces_in_regions;
    for (Size next_pair = 0; next_pair < shared_mt_idx_inty_pairs.size(); ++next_pair)
    {
      const Size first_mt = shared_mt_idx_inty_pairs[next_pair].first;
      if (traces_in_regions.count(first_mt) > 0) continue;
      // double check if this trace is shared
      if (shared_m_traces_indices[first_mt].size() < 2)
      {
        continue;
      }

      /// 3.1 collecting all shared mts and features in this conflicting region
      std::unordered_map<Size, bool> shared_mts_in_this_region; // to track all detected shared mass traces. key = shared_mt_idx, value = visited?
      std::set<String> visited_features_in_this_region; // to track all detected features. key = [FeatGroupID]&[CS]
      std::vector<Feature> conflicting_features;
      std::vector<Feature> features_not_for_resolution;
      std::vector<MassTrace> conflicting_mts; // collect shared traces for resolution (not FeatureSeed, but original MassTrace. FeatureSeed has conflicting information)
      std::vector<Size> conflicting_mt_indices; // trace_index of collect shared traces for resolution

      shared_mts_in_this_region.emplace(first_mt, false); // add the first shared mt
      std::vector<Size> region_queue {first_mt};
      for (Size q = 0; q < region_queue.size(); ++q)
      {
        Size this_shared_mt_index = region_queue[q];
        shared_mts_in_this_region[this_shared_mt_index] = true; // mark as visited

        // per conflicting features, make Feature
        for (auto &feat : shared_mt_N_features[this_shared_mt_index])
        {
          FeatureGroup &org_featgroup = feature_groups[feat.first];
          int cs_of_feat = feat.second;

          // check if FeatureElement is not created yet
          String feat_label = std::to_string(feat.first) + '&' + std::to_string(cs_of_feat);
          if (visited_features_in_this_region.find(feat_label) != visited_features_in_this_region.end())
          {
            continue;
          }
          visited_features_in_this_region.insert(feat_label); // mark as this Feature is visited

          // initialize Feature
          Feature new_feature;
          new_feature.prepareVectors(org_featgroup.size());
          new_feature.charge = cs_of_feat;
          new_feature.feature_group_index = feat.first;

          // add seeds in element
          for (auto &feat_seed: org_featgroup)
          {
            if (feat_seed.getCharge() != cs_of_feat)
            {
              continue;
            }

            if (shared_m_traces_indices[feat_seed.getTraceIndex()].size() > 1) // this seed is shared
            {
              new_feature.shared_traces.push_back(feat_seed);
              new_feature.shared_trace_indices.push_back(feat_seed.getTraceIndex());
            }
            else // this seed is not shared
            {
              new_feature.unique_traces.push_back(feat_seed);
              new_feature.unique_trace_indices.push_back(feat_seed.getTraceIndex());
            }
          }

          // check if this Feature is eligible for resolution
          if (!isEligibleFeatureForConflictResolution_(new_feature, shared_m_traces_indices, org_featgroup))
          {
            // this feature will not be used
            features_not_for_resolution.push_back(new_feature);
            continue;
          }
          new_feature.shrinkVectors();

          // if any seed is shared, mark for visiting in the next round
          for (auto &mt_idx : new_feature.shared_trace_indices)
          {
            if (shared_mts_in_this_region.find(mt_idx) == shared_mts_in_this_region.end()) // if this shared trace is not visited yet
            {
              shared_mts_in_this_region.emplace(mt_idx, false);
              region_queue.push_back(mt_idx);
            }
          }
          conflicting_features.push_back(new_feature);
        }
      }

      /// 3.2 the conflicting MTs visited in this round start no other region
      traces_in_regions.insert(region_queue.begin(), region_queue.end());

      /// 3.3 remove invalid feature from FeatureGroup and then collect conflicting mass traces
      for (auto &feat : features_not_for_resolution)
      {
        auto &feat_group = feature_groups[feat.feature_group_index];
        auto &feat_cs = feat.charge;

        auto pair_to_remove = std::make_pair(feat.feature_group_index, feat_cs);
        for (auto &trace_idx : feat.shared_trace_indices)
        {
          // update the lookup table to collect shared masstrace to be updated later
          auto &feat_vec = shared_mt_N_features[trace_idx];
          feat_vec.erase(std::find(feat_vec.begin(), feat_vec.end(), pair_to_remove));

          // if this seed is from "ineligible feature" and shared, remove it
          feat_group.erase(std::remove_if(feat_group.begin(), feat_group.end(), [trace_idx](auto const &seed) {return seed.getTraceIndex() == trace_idx;}));
          auto &link = shared_m_traces_indices[trace_idx];
          link.erase(std::find(link.begin(), link.end(), feat.feature_group_index)); // update "look up table" for shared traces
        }

        // only update members for scoring. (FeatureGroup can be empty at this point)
        feat_group.updateMembersForScoring();
      }

      // collect shared mass trace to be resolved after removing ineligible features
      for (const Size shared_mt_idx : region_queue)
      {
        auto &conflict_feat_tags = shared_mt_N_features[shared_mt_idx];

        // if this is not used by any features
        if (conflict_feat_tags.empty())
        {
          // mark the shared mt as not shared
          shared_m_traces_indices[shared_mt_idx] = std::vector<Size>();
          continue;
        }
        // check if this mt is not shared one anymore
        else if (conflict_feat_tags.size() == 1)
        {
          Size fg_id = conflict_feat_tags[0].first; // feature group id
          int f_charge = conflict_feat_tags[0].second; // feature charge

          // mark this MT is not shared one
          shared_m_traces_indices[shared_mt_idx] = std::vector<Size>{fg_id};

          // if the corresponding feature is already added to conflict_features, update it
          auto found_feat = find_if(conflicting_features.begin(), conflicting_features.end(),
                                    [fg_id, f_charge](auto &&x) { return (x.feature_group_index == fg_id) && (x.charge == f_charge); });
          if (found_feat == conflicting_features.end())
          {
            continue;
          }
          // find the shared trace
          auto it = std::find(found_feat->shared_trace_indices.begin(), found_feat->shared_trace_indices.end(), shared_mt_idx);
          Size trace_loc = it - found_feat->shared_trace_indices.begin();

          // move the shared trace to unique
          found_feat->unique_trace_indices.push_back(shared_mt_idx);
          FeatureSeed the_trace = found_feat->shared_traces[trace_loc];
          found_feat->unique_traces.push_back(the_trace);

          // delete the trace from the shared trace vectors
          found_feat->shared_trace_indices.erase(found_feat->shared_trace_indices.begin()+trace_loc);
          found_feat->shared_traces.erase(found_feat->shared_traces.begin()+trace_loc);

          // check if any unique_trace is left
          if (found_feat->shared_traces.empty())
          {
            // remove from conflicting features
            conflicting_features.erase(found_feat);
          }
          continue;
        }
        // this trace is shared one -> add this shared mt as a conflicting one
        auto i = input_masstraces.begin() + shared_mt_idx;
        conflicting_mts.push_back(*i);
        conflicting_mt_indices.push_back(shared_mt_idx);
      }
      // if no resolving conflict is needed
      if (conflicting_features.size() < 2)
      {
        continue;
      }

//      /// 3.4 set isotope probabilities for features (to be used as theoretical intensity when modeling a theoretical shape)
//      double minimum_probability = 1.0;
//      for (auto &feat : conflicting_features)
//      {
//        auto &feat_group = feature_groups[feat.feature_group_index];
//        auto theo_isodist = iso_model_.get(feat_group.getMonoisotopicMass());
//        feat.isotope_probabilities.reserve(feat.unique_traces.size()); // isotope prob for each mass_trace, based on their iso position
//        for (auto &lmt: feat.unique_traces)
//        {
//          // if this trace is taken from different Feature
//          if (lmt.getCharge() != feat.charge)
//          {
//            feat.isotope_probabilities.push_back(-1); // to be set later (as the minimum value)
//            continue;
//          }
//
//          // if isotope index of lmt exceed tmp_iso length, give 0 (to not use for the modeling)
//          Size tmp_iso_idx = (Size) lmt.getIsotopeIndex();
//          if (tmp_iso_idx < theo_isodist.size() && tmp_iso_idx >= 0) // To avoid bad access error. If index is larger than the isodist length, no intensity to use.
//          {
//            // give weight per masstrace based on their isotope position, not the real intensity.
//            double tmp_prob = theo_isodist[tmp_iso_idx].getIntensity();
//            feat.isotope_probabilities.push_back(tmp_prob);
//            if (tmp_prob < minimum_probability)
//            {
//              minimum_probability = tmp_prob;
//            }
//          }
//          else
//          {
//            feat.isotope_probabilities.push_back(0.0);
//          }
//        }
//      }
//      // switch -1 probability to minimum prob
//      minimum_probability = minimum_probability < 0 ? 0 : minimum_probability;
//      for (auto &feat : conflicting_features)
//      {
//        for (auto &p: feat.isotope_probabilities)
//        {
//          if (p == -1)
//          {
//            p = minimum_probability;
//          }
//        }
//      }

      /// 3.5 resolve the conflict in this region
      resolveConflictRegion_(conflicting_features, conflicting_mts, conflicting_mt_indices, &fit_cache);

      /// 3.6 Update FeatureGroup (only FeatureSeed vectors, other members will be updated later)
      for (auto &feat : conflicting_features)
      {
        // get original FeatureGroup, and update it
        FeatureGroup &feature_group = feature_groups[feat.feature_group_index];
        std::vector<Size> solved_traces;
        solved_traces.reserve(feat.shared_trace_indices.size());
        for (Size t_index = 0; t_index < feat.shared_trace_indices.size(); ++t_index)
        {
          auto &updated_seed = feat.shared_traces[t_index];
          Size &index_of_updated_seed = feat.shared_trace_indices[t_index];

          // find corresponding FeatureSeed in FeatureGroup
          auto found_seed = find_if(feature_group.begin(), feature_group.end(),
                                    [index_of_updated_seed](auto &&x) { return x.getTraceIndex() == index_of_updated_seed; });
          // this seed is not changed
          if (updated_seed.getIntensity() == found_seed->getIntensity())
          {
            continue;
          }

          // if updated seed is empty, add remove found seed
          if (updated_seed.getIntensity() == 0)
          {
            feature_group.erase(found_seed);
            continue;
          }
          solved_traces.push_back(index_of_updated_seed);
        }
        if (!solved_traces.empty())
        {
          solved_traces.shrink_to_fit();
          feat.solved_trace_indices.swap(solved_traces);
          solved_features.insert({feat.feature_group_index, feat}); // FeatureGroups will be updated later
        }
      }
    }

    /// 4.1 double checking feature groups are still conflicting
    std::map<Size, std::vector<Feature*>> traces_feature_links; // key: mass trace index, value: pointer to solved_features' second
    std::vector<Size> updated_fg_indices_in_this_cluster;
    updated_fg_indices_in_this_cluster.reserve(updated_fg_indices_in_this_cluster.size());
    for (auto &idx: fg_indices_in_this_cluster)
    {
      auto &fgroup = feature_groups[idx];

      // remove 0 intensity mass traces
      auto fg_iter = fgroup.begin();
      while (fg_iter != fgroup.end())
      {
        if (fg_iter->getIntensity() == 0)
        {
          fg_iter = fgroup.erase(fg_iter);
        }
        else
        {
          ++fg_iter;
        }
      }
      // check if enough number of mass traces are left
      if (fgroup.size() < min_nr_mtraces_ && !fd_seeded_)
      {
        continue;
      }

      updated_fg_indices_in_this_cluster.push_back(idx);
      // save solved mass trace indices
      const auto solved_of_group = solved_features.equal_range(idx); // (was a scan of all solved features per group: quadratic)
      for (auto solved_it = solved_of_group.first; solved_it != solved_of_group.second; ++solved_it)
      {
        auto &solved = *solved_it;

        for(auto &solved_trace_idx : solved.second.solved_trace_indices)
        {
          auto map_itr = traces_feature_links.find(solved_trace_idx);
          if (map_itr==traces_feature_links.end())
          {
            traces_feature_links.insert({solved_trace_idx, std::vector<Feature*>()});
            map_itr = traces_feature_links.find(solved_trace_idx);
          }
          map_itr->second.push_back(&solved.second);
        }
      }
    }

    /// 4.2 update FeatureGroup members and quantities
    for (auto &idx: updated_fg_indices_in_this_cluster)
    {
      auto &fgroup = feature_groups[idx];
      std::vector<FeatureSeed> theoretical_shapes;
      theoretical_shapes.reserve(fgroup.size());

      // update trace
      for (auto &seed : fgroup)
      {
        // if this seed was updated
        Size trace_idx_of_this_seed = seed.getTraceIndex();
        auto map_itr = traces_feature_links.find(trace_idx_of_this_seed);
        if (map_itr==traces_feature_links.end() || map_itr->second.size() == 1)
        { // this seed is not shared with other FeatureGroup
          continue;
        }

        // solved features from this fgroup
        for (auto &feat_ptr : map_itr->second)
        {
          if (feat_ptr->feature_group_index != idx)
          { // this feature does not belong to this fgroup
            continue;
          }

          // update this seed with the solved trace
          // find corresponding FeatureSeed in FeatureGroup
          auto &updated_traces = feat_ptr->shared_traces;
          auto found_seed = find_if(updated_traces.begin(), updated_traces.end(),
                                    [trace_idx_of_this_seed](auto &&x) { return x.getTraceIndex() == trace_idx_of_this_seed; });
          seed = *found_seed;
          if (!fd_seeded_) // (with FLASHDeconv masses the resolved seed keeps its trace index: the clusters run in parallel, and nothing reads the new traces)
          {
            seed.setTraceIndex(input_masstraces.size()); // add new trace index
            input_masstraces.push_back(seed.getMassTrace()); // add new FeatureSeed
          }

          // add theoretical shapes
          for (auto &theos : feat_ptr->theoretical_shapes)
          {
            if (theos.getTraceIndex() != trace_idx_of_this_seed) // theoretical shape has original seed index
              continue;
            theoretical_shapes.push_back(theos);
            break; // already found
          }
        }
      }

      if (rescoreFeatureGroup_(fgroup))
      {
        if (shared_output_requested_)
        {
          writeMassTracesOfFeatureGroup_(fgroup, idx, shared_m_traces_indices, false);
          // add theoretical shape information
          fgroup.updateTheoreticalShapes(theoretical_shapes);
          for (auto &theo : theoretical_shapes)
          {
            std::vector<double> shape_inty;
            shape_inty.reserve(theo.getMassTrace().getSize());
            for (auto& p : theo.getMassTrace())
            {
              shape_inty.push_back(p.getIntensity());
            }
            writeTheoreticalShapeForConflictResolution_(idx, theo, shape_inty, theo.getIntensity());
          }
        }
        out_featuregroups.push_back(fgroup);
      }
    }
  }

  bool FLASHQuantAlgorithm::isEligibleFeatureForConflictResolution_(Feature &feat, std::vector<std::vector<Size>> &shared_m_traces_indices, FeatureGroup &feat_group) const
  {
    // check if this feature is composed of only shared mass traces
    if (feat.unique_traces.size() > 0)
    {
      return true;
    }

    // get the most abundant mass trace from FeatureGroup (except recruited ones)
    FeatureSeed *most_abundant_mt_ptr = nullptr;
    double max_intensity = .0; // maximum mt intensity in this feature group
    const int& skip_this_charge = feat.charge;
    for (auto lmt = feat_group.begin(); lmt != feat_group.end(); ++lmt)
    {
      if (skip_this_charge > 0 && lmt->getCharge() == skip_this_charge)
      {
        continue;
      }

      // check if this mt is shared with other FeatureGroup
      if (lmt->getTraceIndex() < shared_m_traces_indices.size() && shared_m_traces_indices[lmt->getTraceIndex()].size() > 1)
      {
        continue;
      }

      if (lmt->getIntensity() > max_intensity)
      {
        max_intensity = lmt->getIntensity();
        most_abundant_mt_ptr = (FeatureSeed *) &(*lmt);
      }
    }

    if (most_abundant_mt_ptr == nullptr)
    {
      return false;
    }

    // get the most abundant mass trace from FeatureGroup (except recruited ones)
    feat.unique_traces.push_back(*most_abundant_mt_ptr);
    feat.unique_trace_indices.push_back(most_abundant_mt_ptr->getTraceIndex());
    return true;
  }

  void FLASHQuantAlgorithm::resolveConflictRegion_(std::vector<Feature> &conflicting_features,
                                                         const std::vector<MassTrace> &conflicting_mts,
                                                         const std::vector<Size> &conflicting_mt_indices,
                                                         std::map<std::vector<Size>, std::shared_ptr<EGHTraceFitter>>* fit_cache)
  {
    // if only one feature has been passed, skip.
    if (conflicting_features.size() < 2)
    {
      return;
    }
    // a region can hold tens of thousands of features and shared traces (Eclip_1480): the fits and the per-trace resolutions are
    // independent and run in parallel unless the region is already resolved inside a parallel loop over clusters
#ifdef _OPENMP
    const bool nested = omp_in_parallel() != 0;
#else
    const bool nested = true;
#endif

    /// the elution model of each feature from its unique traces: the cached ones, then the others fitted, then cached
    const Size feature_count = conflicting_features.size();
    std::vector<std::shared_ptr<EGHTraceFitter>> models(feature_count);
    std::vector<std::vector<Size>> fit_keys(feature_count);
    std::vector<Size> to_fit;
    for (Size i_of_f = 0; i_of_f < feature_count; ++i_of_f)
    {
      const auto &tmp_feat = conflicting_features[i_of_f];
      if (fit_cache != nullptr)
      {
        auto &fit_key = fit_keys[i_of_f];
        fit_key.reserve(tmp_feat.unique_trace_indices.size() + 2);
        fit_key.push_back(tmp_feat.feature_group_index);
        fit_key.push_back((Size)tmp_feat.charge);
        fit_key.insert(fit_key.end(), tmp_feat.unique_trace_indices.begin(), tmp_feat.unique_trace_indices.end());
        const auto cached = fit_cache->find(fit_key);
        if (cached != fit_cache->end())
        {
          models[i_of_f] = cached->second;
          continue;
        }
      }
      to_fit.push_back(i_of_f);
    }
    std::vector<std::exception_ptr> errors(to_fit.size());
#pragma omp parallel for schedule(dynamic, 4) if (!nested && to_fit.size() >= 16)
    for (SignedSize i = 0; i < (SignedSize)to_fit.size(); ++i)
    {
      try
      {
        auto model = std::make_shared<EGHTraceFitter>();
        fitTraceModelFromUniqueTraces_(conflicting_features[to_fit[i]], model.get());
        models[to_fit[i]] = model;
      }
      catch (...)
      {
        errors[i] = std::current_exception();
      }
    }
    for (const auto &error : errors)
    {
      if (error) std::rethrow_exception(error);
    }
    if (fit_cache != nullptr)
    {
      for (const Size i_of_f : to_fit) fit_cache->emplace(fit_keys[i_of_f], models[i_of_f]);
    }

    /// per conflicting trace (row): the features sharing it, in feature order, each with its component, i.e. its model over the trace
    /// normalized to sum 1 (-1: no valid model, or one that is 0 over the whole trace)
    std::unordered_map<Size, Size> row_of_trace;
    row_of_trace.reserve(conflicting_mt_indices.size());
    for (Size row = 0; row < conflicting_mt_indices.size(); ++row)
    {
      row_of_trace.emplace(conflicting_mt_indices[row], row);
    }
    std::vector<std::vector<std::pair<Size, long>>> row_entries(conflicting_mts.size());
    std::vector<std::vector<double>> components;
    for (Size i_of_f = 0; i_of_f < feature_count; ++i_of_f)
    {
      const auto &tmp_feat = conflicting_features[i_of_f];
      EGHTraceFitter* fitted_model = models[i_of_f].get();
      const double area = fitted_model->getArea();
      const bool valid = (area == area); // x != x: test for NaN
      for (const Size trace_index : tmp_feat.shared_trace_indices)
      {
        const auto found = row_of_trace.find(trace_index);
        if (found == row_of_trace.end()) continue;
        const Size row = found->second;
        long component_index = -1;
        if (valid)
        {
          std::vector<double> fit_intensities;
          double summed_intensities = .0;
          for (auto &peak: conflicting_mts[row])
          {
            double fitted_value = fitted_model->getValue(peak.getRT());
            summed_intensities += fitted_value;
            fit_intensities.push_back(fitted_value);
          }
          // if fit_model's RT range doesn't overlap with the conflict_mt -> fit_intensities vec is all composed of 0.
          if (summed_intensities != 0)
          {
            std::vector<double> component;
            component.reserve(fit_intensities.size());
            for (auto &inty: fit_intensities)
            {
              component.push_back(inty / summed_intensities);
            }
            component_index = (long)components.size();
            components.push_back(std::move(component));
          }
        }
        row_entries[row].emplace_back(i_of_f, component_index);
      }
    }

    /// reconstruction of observed XICs (per conflicting mt). A row changes only the copies of its own trace in the features.
    std::vector<std::exception_ptr> row_errors(conflicting_mts.size());
#pragma omp parallel for schedule(dynamic, 16) if (!nested && !shared_output_requested_ && conflicting_mts.size() >= 64)
    for (SignedSize r = 0; r < (SignedSize)conflicting_mts.size(); ++r)
    {
      try
      {
        const Size row = (Size)r;
        const Size org_index_of_this_trace = conflicting_mt_indices[row];
        const auto &entries = row_entries[row];
        const Size column_size = std::count_if(entries.begin(), entries.end(), [](const auto &e) { return e.second >= 0; });

        // if no resolution is needed, skip the next part
        if (column_size < 2) // only one or no Feature has valid model
        {
          for (const auto &[i_of_f, component_index] : entries)
          {
            Feature &feat = conflicting_features[i_of_f];
            auto lmt_ptr = std::find_if(feat.shared_traces.begin(), feat.shared_traces.end(),
                                        [org_index_of_this_trace](auto &&x) { return x.getTraceIndex() == org_index_of_this_trace; });
            if (lmt_ptr == feat.shared_traces.end()) continue;
            if (component_index < 0) // if this feature is not the one with the valid model
            {
              lmt_ptr->setIntensity(0); // to be removed later
              if (shared_output_requested_)
              {
                std::vector<double> zero_vec = std::vector<double>(lmt_ptr->getMassTrace().getSize(), 0);
                feat.setTheoreticalShapes(*lmt_ptr, zero_vec, 0);
              }
            }
            else if (shared_output_requested_) // valid Feature
            {
              feat.setTheoreticalShapes(*lmt_ptr, components[component_index], 1);
            }
          }
          continue;
        }
        updateFeatureWithFitModel(conflicting_features, entries, conflicting_mts[row], org_index_of_this_trace, components);
      }
      catch (...)
      {
        row_errors[r] = std::current_exception();
      }
    }
    for (const auto &error : row_errors)
    {
      if (error) std::rethrow_exception(error);
    }
  }

  void FLASHQuantAlgorithm::updateFeatureWithFitModel(std::vector<Feature>& conflicting_features, const std::vector<std::pair<Size, long>>& entries,
                                                      const MassTrace& obs_masstrace, const Size& org_index_of_obs_mt,
                                                      const std::vector<std::vector<double>>& components)
  {
    // the features of this mt with a valid component, in feature order
    Size column_size = std::count_if(entries.begin(), entries.end(), [](const auto &e) { return e.second >= 0; });

    // prepare observed XIC vector
    Size mt_size = obs_masstrace.getSize();
    double obs_total_intensity = .0;
    for (auto &peak : obs_masstrace)
    {
      obs_total_intensity += peak.getIntensity();
    }
    Matrix<double> obs(mt_size, 1); // TODO: smooth intensity?
    for (Size i = 0; i < mt_size; ++i)
    {
      obs.setValue(i, 0, obs_masstrace[i].getIntensity()/obs_total_intensity); // save normalized value
    }

    // prepare theoretical matrix (include only related features)
    Matrix<double> theo_matrix(mt_size, column_size);
    Size col = 0;
    for (const auto &[i_of_f, component_index] : entries)
    {
      // skipping no-feature column
      if (component_index < 0)
        continue;

      auto &temp_comp = components[component_index];
      for (Size tmp_r = 0; tmp_r < temp_comp.size(); ++tmp_r)
      {
        theo_matrix.setValue(tmp_r, col, temp_comp[tmp_r]);
      }
      col++;
    }

    Matrix<double> out_quant(column_size, 1);
    NonNegativeLeastSquaresSolver::solve(theo_matrix, obs, out_quant);

    // if any out_quant is zero, give the other group all.
    vector<double> calculated_ratio;
    Size zero_ratio_counter = 0;
    for (Size i = 0; i < out_quant.rows(); ++i)
    {
      calculated_ratio.push_back(out_quant.getValue(i, 0));
      if (out_quant.getValue(i, 0) == 0)
      {
        ++zero_ratio_counter;
      }
    }
    if (calculated_ratio.size() - zero_ratio_counter == 1) // only one out_quant was non-zero
    {
      for (auto &ratio:calculated_ratio)
      {
        if (ratio != 0)
        {
          ratio = 1;
        }
      }
    }

    /// update Features based on the calculated ratio
    col = 0;
    for (const auto &[i_of_f, component_index] : entries)
    {
      // skipping not-this-feature column
      if (component_index < 0)
        continue;

      auto &feat = conflicting_features[i_of_f];
      auto lmt_ptr = std::find_if(feat.shared_traces.begin(), feat.shared_traces.end(),
                                  [org_index_of_obs_mt](auto &&x) { return x.getTraceIndex() == org_index_of_obs_mt; });

      // Kyowon's advice! ratio should be applied to real intensity, not theoretical one
      lmt_ptr->setIntensity(calculated_ratio[col] * obs_total_intensity);

      if (this->shared_output_requested_)
      {
        auto theo_vec = theo_matrix.getEigenMatrix().col(col);
        std::vector<double> this_col(theo_vec.begin(), theo_vec.end());
        feat.setTheoreticalShapes(*lmt_ptr, this_col, calculated_ratio[col]);
      }

      /// Update MassTrace itself
      MassTrace updated = this->updateMassTrace_(lmt_ptr->getMassTrace(), calculated_ratio[col]);
      lmt_ptr->setMassTrace(updated);

      col++;
    }
  }

  // modified ElutionModelFitter::fitElutionModels
  void FLASHQuantAlgorithm::fitTraceModelFromUniqueTraces_(FLASHQuantHelper::Feature const& tmp_feat, EGHTraceFitter* fitter) const
  {
    Size mass_traces_size = tmp_feat.unique_traces.size();

    // preparation for ElutionModelFit
    FeatureFinderAlgorithmPickedHelperStructs::MassTraces traces_for_fitting;
    traces_for_fitting.reserve(mass_traces_size);
    traces_for_fitting.max_trace = 0;
    vector<Peak1D> peaks;
    // reserve space once, to avoid copying and invalidating pointers:
    peaks.reserve(tmp_feat.getPeakSizes());

    // get theoretical information from unique mass traces
    for (Size idx = 0; idx < mass_traces_size; ++idx)
    {
      auto &mass_trace_ptr = tmp_feat.unique_traces[idx].getMassTrace();

      // prepare for fitting
      FeatureFinderAlgorithmPickedHelperStructs::MassTrace tmp_mtrace;
      Size m_trace_size = mass_trace_ptr.getSize();
      tmp_mtrace.peaks.reserve(m_trace_size);
      for (auto &p_2d: mass_trace_ptr)
      {
        auto tmp_inty = p_2d.getIntensity();
        if (tmp_inty > 0.0) // only use non-zero intensities for fitting
        {
          Peak1D peak;
          peak.setMZ(p_2d.getMZ());
          peak.setIntensity(tmp_inty);
          peaks.push_back(peak);
          tmp_mtrace.peaks.emplace_back(p_2d.getRT(), &peaks.back());
        }
      }
      tmp_mtrace.updateMaximum();
//      tmp_mtrace.theoretical_int = tmp_feat.isotope_probabilities[idx];
      traces_for_fitting.push_back(tmp_mtrace);
    }

    // ElutionModelFit
//    // TODO : is this necessary? giving isotope probability as a weight
//    Param params = fitter->getDefaults();
//    params.setValue("weighted", "true");
//    fitter->setParameters(params);
    runElutionModelFit_(traces_for_fitting, fitter);
  }

  // update MassTrace
  MassTrace FLASHQuantAlgorithm::updateMassTrace_(const MassTrace& ref_trace, const double &ratio) const
  {
    /// reference = ElutionPeakDetection::detectElutionPeaks_ (step3)
    // if ratio = 0, return empty MassTrace
    if (ratio == 0)
    {
      return MassTrace();
    }

    std::vector<PeakType> tmp_mt;
    std::vector<double> smoothed_tmp;
    double total_inty(0.0);
    for (Size p_index = 0; p_index < ref_trace.getSize(); ++p_index)
    {
      auto &peak = ref_trace[p_index];
      double updated_intensity = peak.getIntensity() * ratio;
      PeakType new_peak(peak);
      new_peak.setIntensity(updated_intensity);
      tmp_mt.push_back(new_peak);
      smoothed_tmp.push_back(ref_trace.getSmoothedIntensities()[p_index] * ratio);
      total_inty += updated_intensity;
    }

    if (total_inty < std::numeric_limits<double>::epsilon())
    {
      return MassTrace();
    }

    // create new mass trace
    MassTrace new_trace(tmp_mt);
    new_trace.setSmoothedIntensities(smoothed_tmp);
    // set label as if it's a sub-trace
    new_trace.setLabel(ref_trace.getLabel() + "_" + String(ratio));
    if (use_smoothed_intensities_)
      new_trace.updateSmoothedMaxRT();
    new_trace.updateWeightedMeanMZ();
    new_trace.updateWeightedMZsd();
    new_trace.setQuantMethod(ref_trace.getQuantMethod());
    new_trace.estimateFWHM(use_smoothed_intensities_);

    return new_trace;
  }

  // from ElutionModelFitter.
  void FLASHQuantAlgorithm::runElutionModelFit_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces &m_traces,
                                                      EGHTraceFitter *fitter) const
  {
    // find the trace with maximal intensity:
    Size max_trace = 0;
    double max_intensity = 0;
    for (Size i = 0; i < m_traces.size(); ++i)
    {
      if (m_traces[i].max_peak->getIntensity() > max_intensity)
      {
        max_trace = i;
        max_intensity = m_traces[i].max_peak->getIntensity();
      }
    }
    m_traces.max_trace = max_trace;
    m_traces.baseline = 0.0;

    // Fitting
    try
    {
      fitter->fit(m_traces);
    }
    catch (Exception::UnableToFit &except)
    {
      OPENMS_LOG_ERROR << "Error fitting model to feature '"
                       << except.getName()
                       << " - " << except.getMessage() << endl;
    }

    // record model parameters:
    //    double center = fitter->getCenter(), height = fitter->getHeight();
    //    double sigma = fitter->getSigma();
    //    double tau = fitter->getTau();
    //    double width = sigma * 0.6266571 + abs(tau);
    //    double asymmetry = abs(tau) / sigma;
    //
    //    double lower_rt_bound = fitter->getLowerRTBound();
    //    double upper_rt_bound = fitter->getUpperRTBound();
  }

  void FLASHQuantAlgorithm::getFLASHDeconvConsensusResult()
  {
    with_target_masses_ = true;
    std::fstream fin;
    fin.open("/Users/jeek/Documents/A4B/FDQ/Kiel-Human/FDQ_target_from_rep9.csv", std::ios::in);
//    fin.open("/Users/jeek/Documents/A4B/FDQ/Kiel-Human/ProSightPD_concensus.csv", std::ios::in);

    std::string line;
    std::getline(fin, line); // skip the first line
    while (std::getline(fin, line))
    {
      Size comma_loc = line.find(',');
      target_masses_.push_back(std::make_pair(stod(line.substr(0, comma_loc)), stod(line.substr(comma_loc + 1))));
    }

    std::sort(target_masses_.begin(), target_masses_.end());
    fin.close();
  }

  bool FLASHQuantAlgorithm::isThisMassOneOfTargets_(const double &candi_mass, const double &candi_rt) const
  {
    auto low_it = std::lower_bound(target_masses_.begin(), target_masses_.end(), std::make_pair(candi_mass - 1.5, .0));
    auto up_it = std::upper_bound(target_masses_.begin(), target_masses_.end(), std::make_pair(candi_mass + 1.5, .0));

    bool is_it_in_the_list = false;
    // if only one mass is found
    if (low_it == up_it)
    {
      // check if any mass is within range
      if (abs(low_it->first-candi_mass) <= 1.5 && abs(low_it->second-candi_rt) < 180)
      {
        return true;
      }
    }
    for (auto tmp_it = low_it; tmp_it != up_it; ++tmp_it)
    {
      // check if any mass is within range
      if (abs(tmp_it->first-candi_mass) <= 1.5 && abs(tmp_it->second-candi_rt) < 180)
      {
        is_it_in_the_list = true;
      }
    }
    return is_it_in_the_list;
  }

  void FLASHQuantAlgorithm::detectElutionPeaks_(MassTrace& mt, std::vector<MassTrace>& single_mtraces, const ElutionPeakDetection& epd) const
  {
    if (mt.getSize() < min_peaks_of_trace_()) return;
    const int polynomial_order = 3;
    const double scan_time = mt.getAverageMS1CycleTime();
    const auto smooth = [&mt, polynomial_order](Size win_size) {
      MSSpectrum spectrum;
      for (Size i = 0; i != mt.getSize(); ++i)
      {
        spectrum.push_back(Peak1D(mt[i].getRT(), mt[i].getIntensity()));
      }
      // one filter per frame length and thread: constructing a filter and setting its parameters computes its coefficients (two SVDs,
      // which were most of the elution peak detection)
      static thread_local std::map<int, SavitzkyGolayFilter> filters;
      const int frame_length = std::max(polynomial_order + 1, (int)win_size); // at least polynomial_order + 1, or SG fails
      auto [sg, inserted] = filters.try_emplace(frame_length);
      if (inserted)
      {
        Param param;
        param.setValue("polynomial_order", polynomial_order);
        param.setValue("frame_length", frame_length);
        sg->second.setParameters(param);
      }
      sg->second.filter(spectrum);
      std::vector<double> smoothed;
      smoothed.reserve(spectrum.size());
      for (const auto& p : spectrum)
      {
        smoothed.push_back(p.getIntensity());
      }
      mt.setSmoothedIntensities(smoothed);
    };

    // Step 1: smooth with a window of the trace's own FWHM, then take the window again from the smoothed FWHM
    Size win_size = std::ceil(mt.estimateFWHM(false) / scan_time); // mt is not smoothed yet
    smooth(win_size);
    win_size = std::ceil(mt.estimateFWHM(true) / scan_time);

    // Step 2: local maxima and minima; a minimum closer than min_nr_peaks_in_mtraces_ peaks to the previous split (or to the end) is no split
    std::vector<Size> maxes, mins;
    epd.findLocalExtrema(mt, std::max(1, (int)win_size / 2), maxes, mins);
    if (!mins.empty())
    {
      Size previous_loc = 0;
      auto iter = mins.begin();
      while (iter != mins.end())
      {
        if ((*iter) - previous_loc < min_nr_peaks_in_mtraces_)
        {
          iter = mins.erase(iter);
        }
        else
        {
          previous_loc = *iter;
          ++iter;
        }
      }
      if (!mins.empty() && (mt.getSize() - 1) - mins.back() < min_nr_peaks_in_mtraces_)
      {
        mins.pop_back();
      }
    }

    // Step 3: the trace itself, or its sub-traces between the minima (no peak width or S/N filter: FLASHQuant's "auto" width filtering)
    if (maxes.empty()) return;
    if (maxes.size() == 1)
    {
      mt.updateSmoothedMaxRT();
      mt.estimateFWHM(true);
      single_mtraces.push_back(mt);
      return;
    }
    MassTrace::const_iterator cp_it = mt.begin();
    Size last_idx = 0;
    mins.push_back(mt.getSize() - 1);
    for (Size min_idx = 0; min_idx < mins.size(); ++min_idx)
    {
      std::vector<Peak2D> tmp_mt;
      std::vector<double> smoothed_tmp;
      while (last_idx <= mins[min_idx])
      {
        tmp_mt.push_back(*cp_it);
        smoothed_tmp.push_back(mt.getSmoothedIntensities()[last_idx]);
        ++cp_it;
        ++last_idx;
      }
      if (tmp_mt.empty()) continue;
      MassTrace new_mt(tmp_mt);
      new_mt.setSmoothedIntensities(smoothed_tmp);
      new_mt.setCentroidIM(mt.getCentroidIM());
      new_mt.fwhm_mz_avg = mt.fwhm_mz_avg;
      new_mt.fwhm_im_avg = mt.fwhm_im_avg;
      new_mt.setLabel(mt.getLabel() + "." + String(min_idx + 1));
      new_mt.updateSmoothedMaxRT();
      new_mt.updateWeightedMeanMZ();
      new_mt.updateWeightedMZsd();
      new_mt.setQuantMethod(mt.getQuantMethod());
      new_mt.estimateFWHM(true);
      single_mtraces.push_back(new_mt);
    }
  }

  double FLASHQuantAlgorithm::groupQuantity_(const FeatureGroup& fg) const
  {
    // the sum over all traces (FLASHQuant). The three largest traces did better on PXD042169 SpikeIn but compressed the ratios of
    // ProteomeMix (human at 1/5: +0.60 log2), the isotope-scaled estimates raised the CVs.
    double quantity = .0;
    for (const auto& seed : fg)
    {
      if (seed.getIsotopeIndex() < 0) continue;
      quantity += seed.computeBulkPeakArea(use_smoothed_intensities_);
    }
    return quantity;
  }

  void FLASHQuantAlgorithm::runFromDeconvolvedMasses(const PeakMap& ms1_map, const std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                                     const FLASHHelperClasses::PrecalculatedAveragine& avg,
                                                     std::vector<std::vector<MassAbundance>>& mass_abundances,
                                                     std::map<uint, MassAbundance>& feature_abundances,
                                                     const std::vector<FLASHHelperClasses::MassFeature>* mass_features)
  {
    fd_seeded_ = true;
    fd_set_aside_.clear();
    iso_model_ = avg;
    // wall time of the steps, logged at the end
    auto clock_start = std::chrono::steady_clock::now();
    std::vector<std::pair<String, double>> step_seconds;
    const auto lap = [&clock_start, &step_seconds](const String& name) {
      const auto now = std::chrono::steady_clock::now();
      step_seconds.emplace_back(name, std::chrono::duration<double>(now - clock_start).count());
      clock_start = now;
    };
    mass_abundances.assign(deconvolved_spectra.size(), std::vector<MassAbundance>());
    feature_abundances.clear();

    // *********************************************************** //
    // Step 1 raw mass traces, as the TOPP tool FLASHQuant detects them
    // *********************************************************** //
    std::vector<MassTrace> input_mtraces;
    std::vector<FeatureSeed> input_seeds;
    detectRawTraces_(ms1_map, input_mtraces, input_seeds);
    lap("mass traces and elution peaks");

    // the trace of each raw peak, by the retention time of its spectrum
    std::unordered_map<double, std::vector<std::pair<double, Size>>> peaks_at_rt;
    for (Size t = 0; t < input_mtraces.size(); ++t)
    {
      for (const auto& p : input_mtraces[t])
      {
        peaks_at_rt[p.getRT()].emplace_back(p.getMZ(), t);
      }
    }
    for (auto& rt_peaks : peaks_at_rt)
    {
      std::sort(rt_peaks.second.begin(), rt_peaks.second.end());
    }
    const auto trace_of = [&peaks_at_rt](const double rt, const double mz) -> long {
      const auto it = peaks_at_rt.find(rt);
      if (it == peaks_at_rt.end()) return -1;
      const auto& rt_peaks = it->second;
      const auto lb = std::lower_bound(rt_peaks.begin(), rt_peaks.end(), mz - 1e-6,
                                       [](const std::pair<double, Size>& a, const double v) { return a.first < v; });
      if (lb == rt_peaks.end() || std::abs(lb->first - mz) > 1e-6) return -1;
      return (long)lb->second;
    };

    // *********************************************************** //
    // Step 2 the groups. A FLASHDeconv feature is one group; a mass without a feature joins the group of another mass (feature or not) when
    // a raw trace carries peaks of both and their masses agree (within 10 ppm, or one isotope apart): the masses FLASHDeconv reports scan
    // by scan for one species without tracing them into a feature are one species, not several that share its signal. Two features are
    // never joined.
    // *********************************************************** //
    struct MassEntry
    {
      Size spectrum, mass;
      uint feature;
      double mono_mass;
      float intensity;
      std::vector<std::pair<Size, int>> traces; // (trace, charge) of its raw peaks
    };
    std::vector<MassEntry> entries;
    for (Size si = 0; si < deconvolved_spectra.size(); ++si)
    {
      const auto& dspec = deconvolved_spectra[si];
      if (dspec.isDecoy() || dspec.getOriginalSpectrum().getMSLevel() != 1) continue;
      const double rt = dspec.getOriginalSpectrum().getRT();
      mass_abundances[si].assign(dspec.size(), MassAbundance(.0, 0));
      for (Size k = 0; k < dspec.size(); ++k)
      {
        const auto& pg = dspec[k];
        if (pg.getTargetDecoyType() != PeakGroup::TargetDecoyType::target) continue;
        MassEntry entry {si, k, pg.getFeatureIndex(), pg.getMonoMass(), pg.getIntensity(), {}};
        for (const auto& p : pg)
        {
          const long t = trace_of(rt, p.mz);
          if (t >= 0) entry.traces.emplace_back((Size)t, p.abs_charge);
        }
        entries.push_back(std::move(entry));
      }
    }

    // union-find over the masses; a root carries the feature of its set (0: none)
    std::vector<Size> parent(entries.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::vector<uint> set_feature(entries.size());
    for (Size e = 0; e < entries.size(); ++e) set_feature[e] = entries[e].feature;
    std::function<Size(Size)> root = [&](Size a) { return parent[a] == a ? a : parent[a] = root(parent[a]); };
    const auto join = [&](Size a, Size b) {
      a = root(a);
      b = root(b);
      if (a == b) return;
      if (set_feature[a] > 0 && set_feature[b] > 0 && set_feature[a] != set_feature[b]) return; // two features
      parent[b] = a;
      if (set_feature[a] == 0) set_feature[a] = set_feature[b];
    };
    {
      std::map<uint, Size> first_of_feature;
      for (Size e = 0; e < entries.size(); ++e)
      {
        if (entries[e].feature == 0) continue;
        const auto it = first_of_feature.emplace(entries[e].feature, e).first;
        join(it->second, e);
      }
      // the masses on each trace, by mass; neighbours within one isotope (+ 10 ppm) are compared
      std::vector<std::vector<Size>> masses_of_trace(input_mtraces.size());
      for (Size e = 0; e < entries.size(); ++e)
      {
        for (const auto& tc : entries[e].traces)
        {
          if (masses_of_trace[tc.first].empty() || masses_of_trace[tc.first].back() != e) masses_of_trace[tc.first].push_back(e);
        }
      }
      for (auto& on_trace : masses_of_trace)
      {
        if (on_trace.size() < 2) continue;
        std::sort(on_trace.begin(), on_trace.end(), [&entries](Size a, Size b) { return entries[a].mono_mass < entries[b].mono_mass; });
        for (Size i = 0; i < on_trace.size(); ++i)
        {
          const double mi = entries[on_trace[i]].mono_mass;
          const double tol = mi * 10e-6;
          for (Size j = i + 1; j < on_trace.size(); ++j)
          {
            const double diff = entries[on_trace[j]].mono_mass - mi;
            if (diff > iso_da_distance_ + tol) break;
            if (diff <= tol || std::abs(diff - iso_da_distance_) <= tol) join(on_trace[i], on_trace[j]);
          }
        }
      }
    }

    struct MassGroup
    {
      double mono_mass = .0;          // of its most intense mass
      float best_intensity = -1;
      double summed_intensity = .0;   // of its masses
      std::map<Size, std::map<int, int>> charges_of_trace; // trace -> charge -> number of its peaks with that charge
      std::vector<std::pair<Size, Size>> masses;            // (deconvolved spectrum, mass)
    };
    std::vector<MassGroup> groups;
    std::map<uint, Size> group_of_feature;
    {
      std::vector<long> group_of_root(entries.size(), -1);
      for (Size e = 0; e < entries.size(); ++e)
      {
        const Size r = root(e);
        if (group_of_root[r] < 0)
        {
          group_of_root[r] = (long)groups.size();
          groups.emplace_back();
        }
        const Size g = (Size)group_of_root[r];
        auto& group = groups[g];
        const auto& entry = entries[e];
        if (entry.feature > 0) group_of_feature[entry.feature] = g;
        group.masses.emplace_back(entry.spectrum, entry.mass);
        group.summed_intensity += entry.intensity;
        if (entry.intensity > group.best_intensity)
        {
          group.best_intensity = entry.intensity;
          group.mono_mass = entry.mono_mass;
        }
        for (const auto& tc : entry.traces)
        {
          group.charges_of_trace[tc.first][tc.second]++;
        }
      }
    }
    std::vector<MassEntry>().swap(entries);
    lap("mass groups");

    // FeatureGroups: per trace the charge most of its peaks have in the group, and its isotope index from the group's monoisotopic mass
    std::vector<FeatureGroup> features;
    features.reserve(groups.size());
    const double charge_mass = FLASHHelperClasses::getChargeMass(true);
    min_mass_ = std::numeric_limits<double>::max();
    max_mass_ = .0;
    for (Size g = 0; g < groups.size(); ++g)
    {
      const auto& group = groups[g];
      min_mass_ = std::min(min_mass_, group.mono_mass);
      max_mass_ = std::max(max_mass_, group.mono_mass);
      if (group.charges_of_trace.empty()) continue;
      FeatureGroup fg(group.mono_mass);
      fg.setFDGroup((int)g);
      const Size max_isotope_index = iso_model_.get(group.mono_mass).size();
      fg.setMaxIsotopeIndex(max_isotope_index);
      for (const auto& [t, charge_counts] : group.charges_of_trace)
      {
        int charge = 0, count = 0;
        for (const auto& [c, n] : charge_counts)
        {
          if (n > count)
          {
            count = n;
            charge = c;
          }
        }
        FeatureSeed seed = input_seeds[t];
        const double mass = (seed.getCentroidMz() - charge_mass) * charge;
        const int iso = (int)std::round((mass - group.mono_mass) / iso_da_distance_);
        if (iso < 0 || iso >= (int)max_isotope_index) continue;
        seed.setCharge(charge);
        seed.setIsotopeIndex(iso);
        seed.setMass(mass);
        fg.push_back(seed);
      }
      if (fg.empty()) continue;
      fg.sort();
      fg.updateMembersForScoring();
      fg.updateMembers(use_smoothed_intensities_);
      features.push_back(fg);
    }
    // A group of at least two traces also takes the traces no mass claims at the m/z of its charges (+-2) and core isotopes (>= 10% of
    // the most abundant) within 10 ppm and apex within half its median FWHM of its apex, if no other group wants them: FLASHDeconv
    // reports fewer peaks per mass at a lower amount (PXD042169: quantified features 936 -> 943 and 873 -> 881, ProtMix |log2 error|
    // 0.770 -> 0.737). Traces another mass claims are not taken: they made shared traces of most of them.
    {
      std::vector<std::pair<double, Size>> by_mz(input_mtraces.size());
      for (Size t = 0; t < input_mtraces.size(); ++t) by_mz[t] = {input_seeds[t].getCentroidMz(), t};
      std::sort(by_mz.begin(), by_mz.end());
      std::vector<char> claimed(input_mtraces.size(), 0);
      for (const auto& fg : features) for (const auto& seed : fg) claimed[seed.getTraceIndex()] = 1;
      std::vector<std::vector<std::pair<Size, FeatureSeed>>> wanted(features.size()); // per group: (trace, seed)
      std::vector<int> claimants(input_mtraces.size(), 0);
#pragma omp parallel for schedule(dynamic, 64)
      for (SignedSize f = 0; f < (SignedSize)features.size(); ++f)
      {
        const auto& fg = features[f];
        if (fg.size() < 2) continue;
        const double mono = groups[fg.getFDGroup()].mono_mass;
        int zmin = std::numeric_limits<int>::max(), zmax = 0;
        std::vector<double> fwhms;
        const FeatureSeed* apex_seed = nullptr;
        for (const auto& seed : fg)
        {
          zmin = std::min(zmin, seed.getCharge());
          zmax = std::max(zmax, seed.getCharge());
          fwhms.push_back(seed.getFwhmEnd() - seed.getFwhmStart());
          if (apex_seed == nullptr || seed.getIntensity() > apex_seed->getIntensity()) apex_seed = &seed;
        }
        std::nth_element(fwhms.begin(), fwhms.begin() + fwhms.size() / 2, fwhms.end());
        const double apex = apex_seed->getMassTrace().getCentroidRT(), half = std::max(1.0, fwhms[fwhms.size() / 2] / 2);
        const auto iso = iso_model_.get(mono);
        double iso_max = 0;
        for (Size k = 0; k < iso.size(); ++k) iso_max = std::max(iso_max, (double)iso[k].getIntensity());
        for (int z = std::max(1, zmin - 2); z <= zmax + 2; ++z)
        {
          for (Size i = 0; i < iso.size(); ++i)
          {
            if (iso[i].getIntensity() < 0.1 * iso_max) continue;
            const double mz = (mono + i * iso_da_distance_) / z + charge_mass, tol = mz * 10e-6;
            for (auto it = std::lower_bound(by_mz.begin(), by_mz.end(), std::make_pair(mz - tol, (Size)0)); it != by_mz.end() && it->first <= mz + tol; ++it)
            {
              const Size t = it->second;
              if (claimed[t]) continue;
              if (std::abs(input_mtraces[t].getCentroidRT() - apex) > half) continue;
              FeatureSeed seed = input_seeds[t];
              seed.setCharge(z);
              seed.setIsotopeIndex((int)i);
              seed.setMass((seed.getCentroidMz() - charge_mass) * z);
              wanted[f].emplace_back(t, seed);
#pragma omp atomic
              claimants[t]++;
            }
          }
        }
      }
      Size recruited = 0;
      for (Size f = 0; f < features.size(); ++f)
      {
        bool changed = false;
        for (auto& [t, seed] : wanted[f])
        {
          if (claimants[t] != 1) continue;
          features[f].push_back(seed);
          changed = true;
          recruited++;
        }
        if (!changed) continue;
        features[f].sort();
        features[f].updateMembersForScoring();
        features[f].updateMembers(use_smoothed_intensities_);
      }
      OPENMS_LOG_INFO << "FLASHQuant from FLASHDeconv masses: " << recruited << " unclaimed traces recruited" << std::endl;
    }
    OPENMS_LOG_INFO << "FLASHQuant from FLASHDeconv masses: " << input_mtraces.size() << " raw mass traces, " << groups.size() << " mass groups ("
                    << group_of_feature.size() << " features), " << features.size() << " with traces" << std::endl;

    // *********************************************************** //
    // Step 3 FLASHQuant: noise level, then the resolution of the traces the groups share
    // *********************************************************** //
    if (!features.empty())
    {
      // FLASHQuant's global noise level (the median peak intensity of the traces no group uses, subtracted from every peak) is not applied:
      // on a lysate run it is 5e5 and removed the whole elution profile of a third of the masses with traces
      if (resolving_shared_signal_)
      {
        clusterFeatureGroups_(features, input_mtraces);
      }
      lap("shared traces");
    }

    // *********************************************************** //
    // Step 4 abundances: FLASHQuant's quantity of each group; a group without a raw trace: NaN
    // *********************************************************** //
    std::vector<MassAbundance> group_abundance(groups.size(), MassAbundance(.0, 0));
    for (const auto& fg : features)
    {
      if (fg.getFDGroup() < 0) continue;
      Size trace_count = 0;
      for (const auto& seed : fg)
      {
        if (seed.getIsotopeIndex() >= 0) trace_count++;
      }
      auto& abundance = group_abundance[fg.getFDGroup()];
      abundance.first += groupQuantity_(fg);
      abundance.second += trace_count;
    }
    Size without_trace = 0, emptied = 0;
    for (Size g = 0; g < groups.size(); ++g)
    {
      auto& abundance = group_abundance[g];
      if (groups[g].charges_of_trace.empty()) // no raw trace of at least 4 peaks (a peak seen in one scan): no abundance
      {
        abundance = MassAbundance(std::numeric_limits<double>::quiet_NaN(), 0);
        without_trace++;
      }
      else if (abundance.second == 0)
      {
        emptied++;
      }
      for (const auto& [si, k] : groups[g].masses)
      {
        mass_abundances[si][k] = abundance;
      }
    }
    for (const auto& [feature_index, g] : group_of_feature)
    {
      feature_abundances[feature_index] = group_abundance[g];
    }
    // A feature without a raw abundance: its FeatureQuantity on the scale of the raw abundances (PXD042169, every feature quantified as
    // with FeatureQuantity: SpikeIn ProtMix |log2 error| 0.935 -> 0.764, ProteomeMix human 0.358 -> 0.306, against FeatureQuantity).
    Size filled = 0;
    if (mass_features != nullptr)
    {
      std::map<uint, double> feature_quantity;
      std::vector<double> ratios;
      for (const auto& mf : *mass_features)
      {
        if (mf.ms_level != 1 || mf.is_decoy) continue;
        const double quantity = mf.mt.computePeakArea();
        feature_quantity[mf.index] = quantity;
        const auto it = feature_abundances.find(mf.index);
        if (it != feature_abundances.end() && it->second.first > 0 && quantity > 0) ratios.push_back(it->second.first / quantity);
      }
      if (!ratios.empty())
      {
        std::nth_element(ratios.begin(), ratios.begin() + ratios.size() / 2, ratios.end());
        const double scale = ratios[ratios.size() / 2];
        for (const auto& [feature_index, g] : group_of_feature)
        {
          auto& abundance = feature_abundances[feature_index];
          if (abundance.first > 0) continue; // (NaN or 0)
          const auto q = feature_quantity.find(feature_index);
          if (q == feature_quantity.end() || !(q->second > 0)) continue;
          abundance = MassAbundance(q->second * scale, 0);
          for (const auto& [si, k] : groups[g].masses) mass_abundances[si][k] = abundance;
          filled++;
        }
      }
    }
    lap("abundances");
    OPENMS_LOG_INFO << "FLASHQuant from FLASHDeconv masses: " << (groups.size() - without_trace) << " groups with raw traces (" << emptied
                    << " left without a trace by the resolution of shared traces: 0), " << without_trace << " without (NaN); " << filled
                    << " features filled from their FeatureQuantity; wall time";
    for (const auto& [name, seconds] : step_seconds)
    {
      OPENMS_LOG_INFO << " " << name << " " << String::number(seconds, 1) << " s";
    }
    OPENMS_LOG_INFO << std::endl;
    fd_seeded_ = false;
  }

  void FLASHQuantAlgorithm::detectRawTraces_(const PeakMap& ms1_map, std::vector<MassTrace>& input_mtraces, std::vector<FeatureSeed>& input_seeds) const
  {
    Param p_mtd = MassTraceDetection().getDefaults();
    p_mtd.setValue("chrom_peak_snr", 3.0);
    p_mtd.setValue("mass_error_ppm", 10.0);
    p_mtd.setValue("trace_termination_outliers", 2);
    // The traces are detected in m/z slices in parallel (MassTraceDetection is sequential; 17 s on a lysate run): each slice with 100 ppm
    // of its neighbours on both sides, keeping the traces whose centroid lies in the slice itself. A trace stays within its mass error
    // (10 ppm) of its centroid, so the slices give the traces of the whole map, but for traces seeded differently near a slice border.
    // (MassTraceDetection on the whole map: 986,396 traces; in 64 slices: about as many.)
    std::vector<double> sample_mzs;
    for (Size i = 0; i < ms1_map.size(); i += std::max((Size)1, ms1_map.size() / 200))
    {
      for (const auto& peak : ms1_map[i]) sample_mzs.push_back(peak.getMZ());
    }
    std::sort(sample_mzs.begin(), sample_mzs.end());
    // a fixed number of slices, so that the result does not depend on the number of threads (the resolution of shared traces carries
    // small differences of the traces near a slice border on: Eclip_1480, 7 vs 48 slices: 80 of 9.3k features differ by > 10%)
    const int slice_count = sample_mzs.size() < 1000 ? 1 : 64;
    std::vector<double> borders {0.0};
    for (int k = 1; k < slice_count; ++k) borders.push_back(sample_mzs[sample_mzs.size() * k / slice_count]);
    borders.push_back(std::numeric_limits<double>::max());
    std::vector<std::vector<MassTrace>> slice_traces(slice_count);
#pragma omp parallel for schedule(dynamic, 1)
    for (int k = 0; k < slice_count; ++k)
    {
      const double lo = borders[k], hi = borders[k + 1], lo_margin = lo * (1 - 100e-6), hi_margin = hi * (1 + 100e-6);
      PeakMap slice;
      slice.reserve(ms1_map.size());
      for (const auto& spec : ms1_map)
      {
        MSSpectrum part;
        static_cast<SpectrumSettings&>(part) = spec;
        part.setRT(spec.getRT());
        part.setMSLevel(spec.getMSLevel());
        part.setNativeID(spec.getNativeID());
        for (auto it = spec.MZBegin(lo_margin); it != spec.end() && it->getMZ() < hi_margin; ++it) part.push_back(*it);
        slice.addSpectrum(std::move(part));
      }
      slice.updateRanges();
      MassTraceDetection slice_mtd;
      slice_mtd.setLogType(ProgressLogger::NONE);
      slice_mtd.setParameters(p_mtd);
      std::vector<MassTrace> found;
      slice_mtd.run(slice, found);
      for (auto& mt : found)
      {
        const double centroid = mt.getCentroidMZ();
        if (centroid >= lo && centroid < hi) slice_traces[k].push_back(std::move(mt));
      }
    }
    std::vector<MassTrace> traces;
    Size found_count = 0;
    for (const auto& found : slice_traces) found_count += found.size();
    traces.reserve(found_count);
    for (auto& found : slice_traces)
    {
      for (auto& mt : found) traces.push_back(std::move(mt));
      std::vector<MassTrace>().swap(found);
    }
    // a canonical order (by centroid m/z, then first retention time): the resolution of shared traces depends on the order of the traces
    // for some groups (Eclip_1480: 71 of 9.3k features changed by > 10% with the order of the slices)
    std::sort(traces.begin(), traces.end(), [](const MassTrace& a, const MassTrace& b) {
      if (a.getCentroidMZ() != b.getCentroidMZ()) return a.getCentroidMZ() < b.getCentroidMZ();
      return a.begin()->getRT() < b.begin()->getRT();
    });

    ElutionPeakDetection epd;
    std::vector<std::vector<MassTrace>> split_traces(traces.size());
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize i = 0; i < (SignedSize)traces.size(); ++i)
    {
      detectElutionPeaks_(traces[i], split_traces[i], epd);
    }
    std::vector<MassTrace>().swap(traces);
    Size split_count = 0;
    for (const auto& sub_traces : split_traces) split_count += sub_traces.size();
    input_mtraces.reserve(input_mtraces.size() + split_count);
    input_seeds.reserve(input_seeds.size() + split_count);
    for (auto& sub_traces : split_traces)
    {
      for (auto& mt : sub_traces)
      {
        if (mt.getSize() < min_peaks_of_trace_() || mt.getFWHM() == 0) continue;
        FeatureSeed seed(mt, use_smoothed_intensities_);
        seed.setTraceIndex(input_mtraces.size());
        input_seeds.push_back(std::move(seed));
        input_mtraces.push_back(std::move(mt));
      }
      std::vector<MassTrace>().swap(sub_traces);
    }
  }

  void FLASHQuantAlgorithm::quantifyTargets(const PeakMap& ms1_map, const std::vector<QuantTarget>& targets,
                                            const FLASHHelperClasses::PrecalculatedAveragine& avg, std::vector<MassAbundance>& abundances)
  {
    // per target, the raw traces at the m/z of its charges and core isotopes (>= 10% of the most abundant) within 10 ppm whose apex
    // lies in its retention time window, their bulk areas summed (no resolution of shared traces)
    const bool was_seeded = fd_seeded_;
    fd_seeded_ = true; // the trace detection of FLASHDeconv masses
    std::vector<MassTrace> input_mtraces;
    std::vector<FeatureSeed> input_seeds;
    detectRawTraces_(ms1_map, input_mtraces, input_seeds);
    fd_seeded_ = was_seeded;
    std::vector<std::pair<double, Size>> by_mz(input_mtraces.size());
    for (Size t = 0; t < input_mtraces.size(); ++t) by_mz[t] = {input_seeds[t].getCentroidMz(), t};
    std::sort(by_mz.begin(), by_mz.end());
    const double charge_mass = FLASHHelperClasses::getChargeMass(true);
    abundances.assign(targets.size(), MassAbundance(.0, 0));
#pragma omp parallel for schedule(dynamic, 64)
    for (SignedSize k = 0; k < (SignedSize)targets.size(); ++k)
    {
      const auto& target = targets[k];
      const auto iso = avg.get(target.mass);
      double iso_max = 0;
      for (Size i = 0; i < iso.size(); ++i) iso_max = std::max(iso_max, (double)iso[i].getIntensity());
      std::set<Size> used;
      double area = 0;
      for (int z = target.min_charge; z <= target.max_charge; ++z)
      {
        for (Size i = 0; i < iso.size(); ++i)
        {
          if (iso[i].getIntensity() < 0.1 * iso_max) continue;
          const double mz = (target.mass + i * iso_da_distance_) / z + charge_mass, tol = mz * 10e-6;
          for (auto it = std::lower_bound(by_mz.begin(), by_mz.end(), std::make_pair(mz - tol, (Size)0)); it != by_mz.end() && it->first <= mz + tol; ++it)
          {
            const Size t = it->second;
            const double apex = input_mtraces[t].getCentroidRT();
            if (apex < target.rt_start || apex > target.rt_end || !used.insert(t).second) continue;
            area += input_seeds[t].computeBulkPeakArea(use_smoothed_intensities_);
          }
        }
      }
      abundances[k] = MassAbundance(area, used.size());
    }
  }
}
