// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHQuantHelper.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FEATUREFINDER/EGHTraceFitter.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/FEATUREFINDER/TraceFitter.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <map>
#include <memory>
#include <algorithm>
#include <fstream>

using namespace std;
namespace OpenMS
{
  class OPENMS_DLLAPI FLASHQuantAlgorithm :
      public ProgressLogger,
      public DefaultParamHandler
{
  public:
    typedef FLASHQuantHelper::FeatureSeed FeatureSeed;
    typedef FLASHQuantHelper::FeatureGroup FeatureGroup;
    typedef FLASHQuantHelper::Feature Feature;

    /// Default constructor
    FLASHQuantAlgorithm();

    /// Default destructor
    ~FLASHQuantAlgorithm() = default;

    /// copy constructor
    FLASHQuantAlgorithm(const FLASHQuantAlgorithm& ) = delete;

    /// move constructor
    FLASHQuantAlgorithm(FLASHQuantAlgorithm&& other) = default;

    /// assignment operator
    FLASHQuantAlgorithm& operator=(const FLASHQuantAlgorithm& fd);

    /// main method of FeatureFindingMetabo
    void run(std::vector<MassTrace> &input_mtraces, std::vector<FeatureGroup> &out_fgroups);

    /// abundance of a FLASHDeconv mass and the number of raw mass traces it is summed over (0: no trace of at least 4 peaks, see below)
    typedef std::pair<double, Size> MassAbundance;

    /**
      @brief FLASHDeconv -quant: quantify every MS1 mass of FLASHDeconv from its raw signal. The masses are the starting point (they replace
      FLASHQuant's own grouping of mass traces by deconvolution): the masses of one FLASHDeconv feature form one feature group, any other
      mass one of its own, and each raw peak a mass is made of brings in the raw mass trace it lies on (mass traces as FLASHQuant detects
      them, with its charge and the isotope index relative to the group's monoisotopic mass). The traces two groups share are then resolved
      as in FLASHQuant (elution model of each group from its own traces, isotope pattern, non-negative least squares), and the abundance of
      a group is FLASHQuant's quantity: the area of each of its traces down to 10% of its apex, summed. A group without any trace of at least
      4 peaks gets the summed intensities of its masses times the MS1 cycle time (and a trace count of 0).
      @param ms1_map the MS1 spectra as loaded (before FLASHDeconv's peak filtering), each sorted by m/z
      @param deconvolved_spectra FLASHDeconv's deconvolved spectra; the MS1 ones are quantified
      @param avg FLASHDeconv's averagine
      @param mass_abundances out: per deconvolved spectrum (by position), per mass (by position): the abundance of its group
      @param feature_abundances out: per FLASHDeconv feature index: the abundance of its group
      @param mass_features FLASHDeconv's features: a feature without a raw abundance (no trace of at least 4 peaks, or none left by the
      resolution of shared traces) gets its FeatureQuantity times the median Abundance / FeatureQuantity of the run's features with both
      (trace count 0), and its masses with it
    */
    void runFromDeconvolvedMasses(const PeakMap& ms1_map, const std::vector<DeconvolvedSpectrum>& deconvolved_spectra,
                                  const FLASHHelperClasses::PrecalculatedAveragine& avg,
                                  std::vector<std::vector<MassAbundance>>& mass_abundances, std::map<uint, MassAbundance>& feature_abundances,
                                  const std::vector<FLASHHelperClasses::MassFeature>* mass_features = nullptr);

    // test purpose
    String output_file_path_;
    std::vector<double> rts_from_org_scans_;

  protected:
    void updateMembers_() override;

  private:
    /// FLASHQuant's elution peak detection (TOPP FLASHQuant): smoothing window from the trace's own FWHM, Savitzky-Golay of order 3, and
    /// sub-traces of at least min_nr_peaks_in_mtraces_ peaks; @p epd only lends its local extrema search
    void detectElutionPeaks_(MassTrace& mt, std::vector<MassTrace>& single_mtraces, const ElutionPeakDetection& epd) const;

    /// FLASHQuant's quantity of a feature group: the bulk peak area of each trace with an isotope index, summed
    double groupQuantity_(const FeatureGroup& fg) const;

    Param getFLASHDeconvParams_() const;

    // equivalent to FLASHDeconvAlgorithm::generatePeakGroupsFromSpectrum_
    void getFeatureFromSpectrum_(std::vector<FeatureSeed*> &local_traces, std::vector<FeatureGroup> &local_fgroup, const double &rt);

    void buildMassTraceGroups_(std::vector<FeatureSeed> &in_seeds, std::vector<FeatureGroup> &features);

    void filterOutNoisePeaks(std::vector<FeatureGroup> &fgroups, std::vector<MassTrace> &input_mtraces) const;

    bool scoreAndFilterFeatureGroup_(FeatureGroup& fg, double min_iso_score = -1) const;

    void refineFeatureGroups_(std::vector<FeatureGroup>& features);

    bool rescoreFeatureGroup_(FeatureGroup& fg) const;

    bool doFWHMbordersOverlap_(const std::pair<double, double>& border1, const std::pair<double, double>& border2) const;

    bool doMassTraceIndicesOverlap(const FeatureGroup& fg1, const FeatureGroup& fg2, const double overlap_percentage_threshold = 0.5, const bool charge_specific = true) const;

    void clusterFeatureGroups_(std::vector<FeatureGroup>& fgroups,
                               std::vector<MassTrace>& input_mtraces);

    void resolveConflictInCluster_(std::vector<FeatureGroup>& feature_groups,
                                   std::vector<MassTrace> &input_masstraces,
                                   std::vector<std::vector<Size> >& shared_m_traces_indices,
                                   std::set<Size>& hypo_indices,
                                   std::vector<FeatureGroup>& out_features);

    void filterOutIneligibleFeatureGroupsInCluster(std::vector<FeatureGroup>& feature_groups,
                                                   std::vector<std::vector<Size> >& shared_m_traces_indices,
                                                   std::set<Size>& hypo_indices) const;

    void setOptionalDetailedOutput_();

    void writeMassTracesOfFeatureGroup_(const FeatureGroup &fgroup,
                                       const Size &fgroup_idx,
                                       const std::vector<std::vector<Size>> &shared_m_traces_indices,
                                       const bool &is_before_resolution);

    void writeTheoreticalShapeForConflictResolution_(const Size &fgroup_idx,
                                                     const FeatureSeed &shared_mt,
                                                     const std::vector<double> &theo_intensities,
                                                     const double &calculated_ratio);

    /// @p fit_cache (one per cluster): the elution model of a feature by (group, charge, its unique traces in order), fitted once; the same
    /// input gives the same fit, so a feature met again in a later region of the cluster is not fitted again
    void resolveConflictRegion_(std::vector<Feature> &conflicting_features,
                                const std::vector<MassTrace> &conflicting_mts,
                                const std::vector<Size> &conflicting_mt_indices,
                                std::map<std::vector<Size>, std::shared_ptr<EGHTraceFitter>>* fit_cache = nullptr);

    MassTrace updateMassTrace_(const MassTrace& ref_trace, const double &ratio) const;

    void fitTraceModelFromUniqueTraces_(Feature const& tmp_feat, EGHTraceFitter* fitted_model) const;

    void runElutionModelFit_(FeatureFinderAlgorithmPickedHelperStructs::MassTraces &m_traces, EGHTraceFitter* fitter) const;

    /// @p entries: the features of the trace with their component (-1: none), in feature order
    void updateFeatureWithFitModel(std::vector<Feature>& conflicting_features, const std::vector<std::pair<Size, long>>& entries,
                                   const MassTrace& obs_masstrace, const Size& org_index_of_obs_mt,
                                   const std::vector<std::vector<double>>& components);

    void getFLASHDeconvConsensusResult();

    bool isThisMassOneOfTargets_(const double &candi_mass, const double &candi_rt) const;

    void makeMSSpectrum_(std::vector<FeatureSeed *> &local_traces, MSSpectrum &spec, const double &rt) const;

    void setFeatureGroupMembersForResultWriting_(std::vector<FeatureGroup> &f_groups) const;

    bool isEligibleFeatureForConflictResolution_(Feature &new_feature, std::vector<std::vector<Size>> &shared_m_traces_indices, FeatureGroup &feat_group) const;

    /// default parameters
    Size charge_lower_bound_;
    Size charge_upper_bound_;
    double min_mass_;
    double max_mass_;
    double mz_tolerance_; // ppm
    bool resolving_shared_signal_;

    double mass_tolerance_da_; // Da, for feature mass collection
    double iso_da_distance_ = Constants::ISOTOPE_MASSDIFF_55K_U;

    // advanced parameter?
    Size min_nr_mtraces_ = 3; //  minimum number of mass traces to support feature group
    Size min_nr_charges_ = 2; //  minimum number of consecutive bridges to support feature group
    Size min_nr_peaks_in_mtraces_ = 4; // at least 4 is needed for EGHTraceFitter
    /// the shortest trace taken: FLASHDeconv masses 3 (for its area; the elution models are fitted on all unique traces of a feature
    /// together). PXD042169: features filled from FeatureQuantity -11%, SpikeIn ProtMix |log2 error| 0.764 -> 0.739, ProteomeMix human
    /// 0.306 -> 0.297.
    Size min_peaks_of_trace_() const { return fd_seeded_ ? 3 : min_nr_peaks_in_mtraces_; }
    bool use_smoothed_intensities_;

    /// variables for internal use (not for user input)
    FLASHHelperClasses::PrecalculatedAveragine iso_model_;
    Size max_nr_traces_; // calculated from iso_model_ (setAveragineModel())

    /// cosine threshold between observed and theoretical isotope patterns for MS1
    double min_isotope_cosine_;

    /// SpectralDeconvolution class for deconvolution
    SpectralDeconvolution deconv_;

    // loop up table
    std::vector<std::pair<double, double>> target_masses_; // mass and rt
    bool with_target_masses_ = false;

    /// for detailed outputs (mostly test purpose)
    bool shared_output_requested_;
    std::fstream shared_out_stream_ = std::fstream();

    /// runFromDeconvolvedMasses: every group is kept and keeps its mass (no score filter, no isotope offset correction); a group the
    /// conflict resolution sets aside (all traces shared, or a harmonic) keeps its own traces in fd_set_aside_ and is reported with them
    bool fd_seeded_ = false;
    std::vector<FeatureGroup> fd_set_aside_;
};
}