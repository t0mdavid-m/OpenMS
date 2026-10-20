// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandJoin.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <boost/dynamic_bitset.hpp>
#include <iostream>

namespace OpenMS
{
  /**
  @brief FLASHDeconv algorithm: ultrafast mass deconvolution algorithm for top down mass spectrometry dataset
  From MSSpectrum, this class outputs DeconvolvedSpectrum.
  Deconvolution takes three steps:
   i) decharging and select candidate masses - speed up via binning
   ii) collecting isotopes from the candidate masses and deisotoping - peak groups are defined here
   iii) scoring and filter out low scoring masses (i.e., peak groups)
  @ingroup Topdown
*/

  class OPENMS_DLLAPI FLASHDeconvAlgorithm : public DefaultParamHandler, public ProgressLogger
  {
  public:
    /// default constructor
    FLASHDeconvAlgorithm();

    /// copy constructor
    FLASHDeconvAlgorithm(const FLASHDeconvAlgorithm&) = default;

    /// move constructor
    FLASHDeconvAlgorithm(FLASHDeconvAlgorithm&& other)  noexcept = default;

    /// assignment operator
    FLASHDeconvAlgorithm& operator=(const FLASHDeconvAlgorithm& fd) = default;

    /// move assignment operator
    FLASHDeconvAlgorithm& operator=(FLASHDeconvAlgorithm&& fd) = default;

    /// destructor
    ~FLASHDeconvAlgorithm() override = default;

    std::vector<double> getTolerances() const;

    /**
     * @brief Run FLASHDeconv algorithm for @p map and store @p deconvolved_spectra and @p deconvolved_feature
     * @param map the dataset
     * @param deconvolved_spectra the deconvolved spectra will be stored in here
     * @param deconvolved_feature the deconvolved features wll be strored in here
     */
    void run(MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<FLASHHelperClasses::MassFeature>& deconvolved_feature);

    /// get calculated averagine. Call after calculateAveragine is called.
    const FLASHHelperClasses::PrecalculatedAveragine& getAveragine();

    /// get calculated decoy averagine. Call after calculateAveragine is called.
    const FLASHHelperClasses::PrecalculatedAveragine& getDecoyAveragine();

    /// get noise decoy weight
    double getNoiseDecoyWeight() const
    {
      return noise_decoy_weight_;
    }
    /// get scan number of the spectrum in @p index -th in @p map
    static int getScanNumber(const MSExperiment& map, Size index);

  protected:
    void updateMembers_() override;

  private:
    /// SpectralDeconvolution  instances for spectral deconvolution for target and decoys
    SpectralDeconvolution sd_, sd_noise_decoy_, sd_signal_decoy_;

    /// to merge spectra.
    int merge_spec_ = 0;

    /// forced MS level
    int forced_ms_level_ = 0;

    /// maximum MS level, which is 4.
    UInt max_ms_level_ = 4;

    /// current maximum MS level - i.e., for MS2, this is the precursor charge
    UInt current_max_ms_level_ = 0;

    /// current minimum MS level.
    UInt current_min_ms_level_ = 0;

    /// MS1 spectra before and after the precursor spectrum searched for a precursor the precursor spectrum lacks (0: off)
    int precursor_rescue_window_ = 2;
    /// precursors taken from a neighbouring MS1 spectrum in this run
    Size rescued_precursors_ = 0;
    /// choose the precursor of an MSn spectrum by the complementary fragment pairs among the target masses in its isolation window and,
    /// with report_FDR, its decoy precursor among the decoy masses of the window against the chosen precursor (precursor_by_fragments)
    bool precursor_by_fragments_ = true;
    /// MSn spectra with several target masses in the isolation window and those whose precursor the fragments changed, in this run
    int spectra_with_several_window_masses_ = 0, precursors_by_fragments_ = 0;
    /// MSn spectra with a precursor and decoy masses in the isolation window and those in which a decoy mass beat the chosen precursor
    /// (and became the decoy precursor), in this run
    int spectra_with_window_decoys_ = 0, decoy_precursors_by_fragments_ = 0, decoy_precursors_ = 0;

    /// the FLASHIda acquisition's scan_commands.tsv (ADR-0046); empty = no coupling
    String scan_commands_file_;

    /// mass tolerances, and minimum cosine scores per MS level
    DoubleList tols_, min_cos_;

    /// to use RNA averagine model
    bool use_RNA_averagine_ = false;

    /// should decoy deconvolution be done?
    bool report_decoy_ = false;

    /// default precursor isolation window size.
    double isolation_window_size_{3.0};

    /// noise decoy weight determined with qvalue calcualtion.
    double noise_decoy_weight_ = 1;
    /// ADR-0046: the commanded MS2 scans by scan number, each naming the survey its command was decided from
    std::map<int, ScanCommandJoin::Located> located_by_scan_;
    /// a map from native ID to precursor peak
    std::map<String, Precursor> native_id_precursor_peak_map_;
    /// a map from native ID to precursor peak group
    std::map<String, PeakGroup> native_id_precursor_peak_group_map_;
    /// per MSn native id, every target mass in the isolation window: its peak group and its precursor ion (the most intense peak of the
    /// mass in the window, at the charge the window implies). The candidates among which precursor_by_fragments chooses; kept only with
    /// precursor_by_fragments and several masses in the window
    std::map<String, std::vector<std::pair<PeakGroup, Precursor>>> native_id_window_candidates_;
    /// per MSn native id, every MS1 decoy mass in the isolation window as FLASHDeconv generated it (with a peak at the charge the window
    /// implies), and whether it is a noise decoy: the challengers of the chosen precursor for the decoy precursor; kept only with
    /// report_FDR and precursor_by_fragments
    std::map<String, std::vector<std::pair<double, bool>>> native_id_window_decoys_;
    /// per MSn spectrum: the precursor score (2D Qscore) of each window decoy mass, as scored when the window was scanned
    std::map<String, std::map<double, double>> native_id_window_decoy_qscores_;
    /// precursor score of the decoy precursor of each MSn spectrum (refined after the feature tracing, updatePrecursorQScores_)
    std::map<String, double> native_id_decoy_precursor_qscore_map_;
    /// per MSn native id, the decoy precursor mass (with report_FDR). With precursor_by_fragments, the decoy mass of the isolation window
    /// that beat the chosen precursor by the fragments (no entry if none did); otherwise the decoy mass of the window of the highest charge
    /// SNR
    std::map<String, double> native_id_decoy_precursor_mass_map_;
    /// per MSn native id, whether that decoy precursor is a noise decoy (weighted by the noise decoy weight in an FDR)
    std::map<String, bool> native_id_decoy_precursor_noise_map_;
    /// a map from native ID to precursor mass - intensity pair (within isolation window)
    std::map<String, std::vector<std::tuple<double, double>>> native_id_precursor_mass_intensity_map_;

    /// read dataset to update ms level information
    void updateMSLevels_(MSExperiment& map);

    /// merge spectra
    void mergeSpectra_(MSExperiment& map, uint ms_level);

    /// run spectral deconvolution
    void runSpectralDeconvolution_(MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra);

    /// run feature finding to get deconvolved features
    void runFeatureFinding_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, std::vector<FLASHHelperClasses::MassFeature>& deconvolved_features);

    /// with found deconvolved features, update QScores for masses that are contained in features.
    static void updatePrecursorQScores_(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, int ms_level);

    /// ADR-0046: the one-peak precursor peak group of a commanded scan whose survey window holds no
    /// deconvolved mass at all
    PeakGroup peakGroupFromCommand_(const ScanCommandJoin::Located& located) const;

    /// register the precursor peak group (or mass) if possible for MSn (n>1) spectrum.
    void findPrecursorPeakGroupsForMSnSpectra_(const MSExperiment& map, const std::vector<DeconvolvedSpectrum>& deconvolved_spectra, uint ms_level);

    /**
     * Choose a precursor mass of an MSn spectrum among candidate masses by the fragments: a candidate mass M scores its complementary
     * fragment pairs (b + y = M, c + z = M + H; two masses less than 10 Da apart are no pair) minus the average count at shifted masses
     * (chance). The default stays unless another candidate has at least two pairs beyond chance and two more than the default has (the
     * highest such score wins, the first one on a tie). Used up to twice per MSn spectrum on the same fragment masses: among the target
     * masses of the isolation window with the registered precursor as the default (with several of them), and among the chosen precursor
     * (the default) and the decoy masses of the window for the decoy precursor (with decoy masses in the window)
     * @param fragment_masses the masses of the MSn spectrum deconvolved within the bounds of the target masses of the window (with a single
     * one, its own deconvolution), sorted
     * @param cand_masses the candidate masses
     * @param default_mass the mass of the default candidate (the first candidate of that mass)
     * @param tol_ppm the MSn mass tolerance in ppm
     * @return the index of the chosen candidate, or -1 if the default stays (or is no candidate)
     */
    static int choosePrecursorByFragments_(const std::vector<double>& fragment_masses, const std::vector<double>& cand_masses, double default_mass,
                                           double tol_ppm);

    /// determine tolerance
    void determineTolerance_(const MSExperiment& map, const Param& sd_param, const FLASHHelperClasses::PrecalculatedAveragine& avg, uint ms_level);

    /// get histogram
    static std::vector<int> getHistogram_(const std::vector<double>& data, double min_range, double max_range, double bin_size);

    /// filter low intensity peaks
    static void filterLowPeaks_(MSExperiment& map);
  };
} // namespace OpenMS