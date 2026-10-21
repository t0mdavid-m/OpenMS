// Copyright (c) 2002-2025, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <boost/dynamic_bitset.hpp>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>

namespace OpenMS
{
/**
@brief
@ingroup Topdown
*/

class OPENMS_DLLAPI FLASHTnTAlgorithm : public DefaultParamHandler, public ProgressLogger
{
public:
  /// constructor
  FLASHTnTAlgorithm();

  /// destructor
  ~FLASHTnTAlgorithm() override = default;

  /// copy constructor
  FLASHTnTAlgorithm(const FLASHTnTAlgorithm&);

  /// move constructor
  FLASHTnTAlgorithm(FLASHTnTAlgorithm&& other) = default;

  /// assignment operator
  FLASHTnTAlgorithm& operator=(const FLASHTnTAlgorithm& other);

  /// Find sequence tags from @p mzs and @p intensities then store them in @p tags.
  /**
    @brief
    Decoy or MS level 1 spectra are removed by this process.
    Overlapping PeakGroups in merged spectra are also removed.

    @param map spectra deconvolved by FLASHDeconv.
    @param fasta_entry fasta entry to searched against

  */
  void run(const MSExperiment& map, const std::vector<FASTAFile::FASTAEntry>& fasta_entry);
  void getProteoformHitsMatchedBy(const FLASHHelperClasses::Tag& tag, std::vector<ProteinHit>& hits) const;
  void getTags(std::vector<FLASHHelperClasses::Tag>& tags) const;
  void getProteoforms(std::vector<ProteinHit>& hits) const
  {
    for (const auto& hit : proteoform_hits_)
    {
      hits.push_back(hit);
    }
  }
  double getMaxTotalModificationMass() const
  {
    return max_total_mod_mass_;
  }

  /// a function that loads the raw MSn spectra of the given scan numbers (FLASHAnnotatorAlgorithm::loadSpectra on the raw mzML)
  using RawSpectrumLoader = std::function<void(const std::set<int>&, std::map<int, MSSpectrum>&)>;

  /// With a loader, every PrSM (targets and decoys alike) is annotated on its raw spectrum before the FDR (FLASHAnnotatorAlgorithm, the
  /// raw: parameters) and ranked by its extension score plus raw_score_weight times its raw score.
  void setRawSpectrumLoader(RawSpectrumLoader loader)
  {
    raw_loader_ = std::move(loader);
  }

  /// The protein-level hits, as they were ranked before proteoform characterization. Filled in
  /// both modes, so a run can be used to shortlist proteins even when the extension is performed.
  void getProteinHits(std::vector<ProteinHit>& hits) const
  {
    for (const auto& hit : protein_hits_)
    {
      hits.push_back(hit);
    }
  }

protected:
  void updateMembers_() override;
  /// implemented for DefaultParamHandler
  void setDefaultParams_();

private:

  Param tagger_param_, extender_param_, gapped_param_, annotator_param_;
  RawSpectrumLoader raw_loader_;
  /// weight of the raw-spectrum evidence score added to the extension score for the FDR
  double raw_score_weight_ = 8.0;
  /// the raw spectra of loadRawSpectra_, by position (raw_slot_: scan to position): their peaks centroided into raw_peaks_ and removed
  /// from raw_spectra_, which keep the precursor information
  std::vector<MSSpectrum> raw_spectra_;
  std::vector<FLASHAnnotatorAlgorithm::PeakArrays> raw_peaks_;
  std::map<int, int> raw_slot_;
  /// how many of the loaded raw spectra were profile spectra (centroided here)
  int raw_profile_count_ = 0;
  /// load the raw spectra of @p scans (raw_loader_) for annotateOnRawSpectra_
  void loadRawSpectra_(const std::set<int>& scans);
  /// free the raw spectra of loadRawSpectra_
  void releaseRawSpectra_();
  /// annotate every hit of @p hits on its raw spectrum (loaded by loadRawSpectra_) and set the combined score
  void annotateOnRawSpectra_(std::vector<ProteinHit>& hits) const;
  /// one PrSM per spectrum: the best by score (a determined mass and termini break a tie, then the coverage, then the earlier PrSM);
  /// @p label and @p score_name go into the log line
  static void keepBestPerSpectrum_(std::vector<ProteinHit>& hits, const String& label, const String& score_name);
  /// residual search rounds of chimeric spectra after the first search (0: off)
  int chimeric_rounds_ = 0;
  /// candidates from gapped tags: off, fallback (spectra without sequence-tag candidates), both, only
  String gapped_mode_ = "fallback";
  /// candidates from terminal ladders (residue sums from an intact protein terminus) in every spectrum
  bool terminal_ladders_ = true;
  /// candidates from precursor-anchored ladders (windows whose unmodified mass is the precursor mass) in every spectrum with a precursor mass
  bool precursor_ladders_ = true;
  /// learn variable modification penalties from target and decoy usage
  bool learn_var_mod_penalty_ = true;
  /// search spectra with MS1-based PrSMs again with decoy precursor masses, as decoys of the precursor in the FDR
  bool decoy_precursor_ = false;
  /// co-isolation: an MS1 mass reported by PrSMs of different proteins stays with the protein that explains it best
  bool one_proteoform_per_feature_ = true;
  /// native and targeted runs: candidates of every source chosen by the decoy margin of their spectrum, every candidate's PrSM reported
  bool native_ = false;
  /// with native_: the margin (-log10 units) over the best decoy candidate of the same source in the spectrum
  double native_margin_ = 2.0;
  double decoy_factor_ = 0, pro_fdr_ = 1;
  /// how many nodes (massses) will be used
  /// at most this many masses (best Qscore first) per spectrum are used. 300 cut key fragments of native spectra (440-570 masses from
  /// centroid input); 400 recovered them at about 1.8x the native search time, and lysate MS2 spectra rarely have this many masses.
  Size max_node_cntr_ = 400;
  double max_total_mod_mass_ = 0;

  bool keep_decoy_ = false;
  bool keep_underdetermined_ = true;
  bool multiple_hits_per_spec_ = false;
  bool protein_level_only_ = false;
  int max_protein_hits_ = 50;
  /// mass tolerance (ppm) of the protein-level search; <= 0: the MS2 deconvolution tolerance of the input
  double protein_tol_ = 0;
  int min_protein_hits_ = 5;
  double protein_fdr_ = 0.1;
  std::vector<ProteinHit> proteoform_hits_;
  std::vector<ProteinHit> protein_hits_;
  std::vector<FLASHHelperClasses::Tag> tags_;
  std::map<int, std::vector<int>> matching_hits_indices_;
  static bool areConsistent_(const ProteinHit& a, const ProteinHit& b, double tol) ;
  /// number of blind modifications of a proteoform hit (variable modifications excluded)
  static Size blindModCount_(const ProteinHit& hit);
  /// Groups the PrSMs of each protein into proteoforms (areConsistent_ with the representative) and marks the representatives
  /// (meta values Index, Representative, PrSMIndices). With @p by_qvalue, only PrSMs in an FDR class (@p hit_class >= 0) take part, and a
  /// determined PrSM with PrSM q-value (@p qvalues) <= @p rep_q, then the lowest q-value, comes first as representative; otherwise every
  /// PrSM takes part. Then a PrSM with a mass comes first, then the best scoring. Returns the proteoforms as positions in proteoform_hits_,
  /// the representative first, in the order of their representatives.
  std::vector<std::vector<Size>> markRepresentativeProteoformHits_(double tol, const std::vector<int>& hit_class, const std::vector<double>& qvalues,
                                                                   double rep_q, bool by_qvalue);
  /// Protein-level q-values (meta value ProteinQvalue on every hit). @p hit_class is the FDR class of each hit (-1: none; such a hit takes
  /// no part), @p qvalues its PrSM q-value, @p decoys and @p decoy_weights its decoy flag and weight, @p clusters the proteoforms
  /// (positions in proteoform_hits_, the representative first), @p rep_q the PrSM q-value of a passing PrSM.
  void setProteinQvalues_(const std::vector<int>& hit_class, const std::vector<double>& qvalues, const std::vector<bool>& decoys,
                          const std::vector<double>& decoy_weights, const std::vector<std::vector<Size>>& clusters, double rep_q);
  /// fasta indices of the proteins (targets and decoys) whose run-level q-value passes protein_fdr_, from each
  /// protein's best tag-based score in protein_hits_
  std::set<Size> getProteinsPassingFDR_(const std::vector<FASTAFile::FASTAEntry>& fasta_entry) const;
  static void vectorizeProteinSequence_(const std::vector<std::string>& cleaned_protein_seqs,
                            std::vector<std::unordered_set<int>>& vec_pro,
                            std::vector<std::unordered_set<int>>& rev_vec_pro);
};
} // namespace OpenMS