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
#include <OpenMS/ANALYSIS/TOPDOWN/ProteoformDatabase.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <boost/dynamic_bitset.hpp>
#include <cstdint>
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
  /// option protein_evidence proteoforms: whether a protein of @p proteoforms proteoforms takes part in the protein ranking (two or more,
  /// or @p strong: one of them is a determined tag-search proteoform without blind modifications of protein_evidence_min_mass)
  static bool proteinEvidenceSupported(Size proteoforms, bool strong);
  /// the grouping of the proteoform classes for the FDR (option merge_small_classes): in @p class_order, a group takes classes until it
  /// holds @p min_size proteoforms (targets and decoys) and a last smaller group joins the one before; with @p by_kind (kind) within
  /// the two kinds @p kind_of gives (0, 1) separately, a kind that is one small group joining the other kind's first group
  static std::vector<std::vector<int>> groupClasses(const std::vector<int>& class_order, const std::vector<Size>& class_size, Size min_size,
                                                    bool by_kind, const std::function<int(int)>& kind_of);
  /// the blind modifications of a PrSM: its BlindModCount, else the ModificationVariable entries of kind 0, else all its Modifications
  static Size blindModCount(const ProteinHit& hit);
  /// the fixed pseudo-random order of the database entries (splitmix64 of the index) that breaks the ties of the candidate cuts
  static std::uint64_t entryOrder(std::uint64_t index);
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
    @param fasta_entry fasta entry to searched against: the database entries (with @p db, its entries' letters in the same order)
    @param db the database of the entries (ProteoformDatabase::loadPEFF: several entries per protein, fixed modifications); nullptr: the
    entries are the proteins (ProteoformDatabase::fromFASTA)

  */
  void run(const MSExperiment& map, const std::vector<FASTAFile::FASTAEntry>& fasta_entry, const ProteoformDatabase* db = nullptr);
  void getProteoformHitsMatchedBy(const FLASHHelperClasses::Tag& tag, std::vector<ProteinHit>& hits) const;
  void getTags(std::vector<FLASHHelperClasses::Tag>& tags) const;
  void getProteoforms(std::vector<ProteinHit>& hits) const
  {
    for (const auto& hit : proteoform_hits_)
    {
      hits.push_back(hit);
    }
  }
  /// the PrSMs, as getProteoforms copies them
  const std::vector<ProteinHit>& getProteoforms() const
  {
    return proteoform_hits_;
  }
  /// the tags, as getTags copies them
  const std::vector<FLASHHelperClasses::Tag>& getTags() const
  {
    return tags_;
  }
  /// positions in getProteoforms() of the PrSMs matched by the tag of index @p tag_index; nullptr if none
  const std::vector<int>* getProteoformIndicesMatchedBy(int tag_index) const
  {
    auto it = matching_hits_indices_.find(tag_index);
    return it == matching_hits_indices_.end() ? nullptr : &it->second;
  }
  /// the protein-level hits, as getProteinHits copies them
  const std::vector<ProteinHit>& getProteinHits() const
  {
    return protein_hits_;
  }
  /// free the results (PrSMs, protein hits, tags) with all threads; they are empty afterwards
  void releaseResults();
  double getMaxTotalModificationMass() const
  {
    return max_total_mod_mass_;
  }

  /// a function that loads the raw MSn spectra of the given scan numbers (FLASHAnnotatorAlgorithm::loadSpectra on the raw mzML)
  using RawSpectrumLoader = std::function<void(const std::set<int>&, std::map<int, MSSpectrum>&)>;

  /// With a loader, every PrSM (targets and decoys alike) is annotated on its raw spectrum before the FDR (FLASHAnnotatorAlgorithm with its
  /// defaults and raw:tol) and ranked by its extension score plus 8 times its raw score.
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
  /// the raw-spectrum PrSM of @p hit (its proteoform, bounded by the entry where a terminus is open) and its given precursor mass (-1:
  /// none); false if its spectrum is not loaded
  bool toRawPrSM_(const ProteinHit& hit, FLASHAnnotatorAlgorithm::PrSM& prsm, double& precursor_mass) const;
  /// the raw-evidence score of a hit for the path alternatives of the extension (8 x raw score, as annotateOnRawSpectra_ adds it); the
  /// raw spectra must stay loaded while it is used
  std::function<double(const ProteinHit&)> alternativeRawScorer_() const;
  /// path alternatives per candidate protein and extension mode (native: 2, else 0)
  int path_alternatives_ = 0;
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
  /// search cascade: sequence tags and ladders in every spectrum, then gapped tags only in the spectra the first stage identifies nothing in
  bool search_cascade_ = true;
  /// search spectra with MS1-based PrSMs again with decoy precursor masses, as decoys of the precursor in the FDR
  bool decoy_precursor_ = false;
  /// co-isolation: an MS1 mass reported by PrSMs of different proteins stays with the protein that explains it best
  bool one_proteoform_per_feature_ = true;
  /// native and targeted runs: candidates of every source chosen by the decoy margin of their spectrum, every candidate's PrSM reported
  bool native_ = false;
  double decoy_factor_ = 0, pro_fdr_ = 1;
  /// how many nodes (massses) will be used
  /// at most this many masses (best Qscore first) per spectrum are used: 400, lysate MS2 spectra rarely have this many masses; no cap with
  /// native_ (native spectra have 440-570 masses from centroid input, and the cap cut key ladder masses)
  Size max_node_cntr_ = 400;
  double max_total_mod_mass_ = 0;

  bool keep_decoy_ = false;
  bool keep_underdetermined_ = true;
  /// terminal-ladder PrSMs as an FDR and proteoform class of their own
  bool terminal_ladder_class_ = true;
  /// a protein takes part in the protein ranking with two proteoforms or one strong one (proteoforms), or with any (none)
  bool protein_evidence_ = true;
  /// proteoform classes with few proteoforms counted together with the next ones within their kind (kind) or in class order (true)
  bool merge_by_kind_ = true;
  double protein_evidence_min_mass_ = 3000;
  /// several PrSMs per spectrum (native only)
  bool multiple_hits_per_spec_ = false;
  bool protein_level_only_ = false;
  int max_protein_hits_ = 50;
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
  /// (positions in proteoform_hits_, the representative first), @p rep_q the PrSM q-value of a passing PrSM, @p db the database (the
  /// processed forms of a protein, for a whole-protein PrSM).
  void setProteinQvalues_(const std::vector<int>& hit_class, const std::vector<double>& qvalues, const std::vector<bool>& decoys,
                          const std::vector<double>& decoy_weights, const std::vector<std::vector<Size>>& clusters, double rep_q,
                          const ProteoformDatabase& db);
  /// protein indices (into @p db's proteins; the fasta indices of a FASTA search) of the proteins (targets and decoys) whose run-level
  /// q-value passes protein_fdr_, from each protein's best tag-based score over its entries in protein_hits_, and of their partners
  std::set<Size> getProteinsPassingFDR_(const ProteoformDatabase& db) const;
  /// @p fixed: per sequence, (position, mass) of its fixed modifications in the sequence (nullptr: none)
  static void vectorizeProteinSequence_(const std::vector<std::string>& cleaned_protein_seqs,
                            std::vector<std::unordered_set<int>>& vec_pro,
                            std::vector<std::unordered_set<int>>& rev_vec_pro,
                            const std::vector<std::vector<std::pair<int, double>>>* fixed = nullptr);
};
} // namespace OpenMS