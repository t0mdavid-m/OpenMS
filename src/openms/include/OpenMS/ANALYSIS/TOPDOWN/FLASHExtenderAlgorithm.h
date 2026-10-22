// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

namespace OpenMS
{
class ProteoformDatabase;

/**
@brief Extend between tags found by FLASHTaggerAlgorithm. In practice, the proteoform characterization is
 performed here. Mass shifts are explained by variable modifications (a given list, residue-specific) or by blind
 modifications (any mass up to max_mod_mass); the two are counted separately and can occur in the same proteoform.
@ingroup Topdown
*/

class OPENMS_DLLAPI FLASHExtenderAlgorithm : public DefaultParamHandler, public ProgressLogger
{
public:
  /// constructor
  FLASHExtenderAlgorithm();

  /// destructor
  ~FLASHExtenderAlgorithm() override = default;

  /// copy constructor
  FLASHExtenderAlgorithm(const FLASHExtenderAlgorithm&) = default;

  /// move constructor
  FLASHExtenderAlgorithm(FLASHExtenderAlgorithm&& other) = default;

  /// assignment operator
  FLASHExtenderAlgorithm& operator=(const FLASHExtenderAlgorithm& other);

  /**
   * The main run function to perform extension algorithm. Take the candidate protein hits and perform extension for each protein with the tags.
   * @param hits the candidate protein hits
   * @param dspec the deconvolved spectrum
   * @param spec_vec
   * @param vec_pro
   * @param rev_vec_pro
   * @param tags the sequence tags from FLASHTaggerAlgorithm
   * @param ppm mass ppm tolerance
   * @param multiple_hits_per_spec should multiple proteins be considered per spectrum or only the best protein should be considered?
   */
  void run(std::vector<ProteinHit>& hits,
           const DeconvolvedSpectrum& dspec, const std::vector<int>& spec_vec, const std::vector<std::unordered_set<int>>& vec_pro,
           const std::vector<std::unordered_set<int>>& rev_vec_pro, const std::vector<FLASHHelperClasses::Tag>& tags, double ppm, bool multiple_hits_per_spec);

  /// names of the variable modifications, in the order of the var_mod_penalty parameter
  std::vector<String> getVariableModificationNames() const
  {
    std::vector<String> names;
    for (const auto& vm : var_mod_list_) names.push_back(vm.name);
    return names;
  }
  /// residues each variable modification can sit on (empty: any), in the same order
  std::vector<String> getVariableModificationResidues() const
  {
    std::vector<String> residues;
    for (const auto& vm : var_mod_list_) residues.push_back(vm.residues);
    return residues;
  }
  /// position of each variable modification: 0 anywhere, 1 proteoform N-terminal residue, 2 proteoform C-terminal residue
  std::vector<int> getVariableModificationPositions() const
  {
    std::vector<int> positions;
    for (const auto& vm : var_mod_list_) positions.push_back(vm.position);
    return positions;
  }

  /// fill the characterized proteoforms in @p hits
  void fillProteoforms(std::vector<ProteinHit>& hits) const
  {
    for (const auto& hit : proteoform_hits_)
    {
      hits.push_back(hit);
    }
  }

  /// set modification map to be considered for blind modification search. Note that they are not for variable modification search.
  /// The mass shifts found by blind modifications matching to the masses of these modifications are annotated.
  void setCandidateBlindModificationMap(const std::map<double, std::vector<ResidueModification>>& mod_map);

  /// The fixed modifications of the database entries (ProteoformDatabase; the hits' FastaIndex is the entry): their masses are part of
  /// the residue masses, and no variable modification sits where a fixed one does (a side-chain one on a residue with a fixed side-chain
  /// modification, an N- or C-terminal one on a terminal residue with a fixed one of that terminus). nullptr (the default) or a database
  /// without fixed modifications: none. The database must outlive the extension.
  void setFixedModifications(const ProteoformDatabase* db);

  /// With @p report, every proteoform hit carries the deconvolved masses its best path matched (meta value MatchedMassIndices: indices
  /// into the deconvolved spectrum given to run(), ascending): the masses of every node of the path whose protein position counts as
  /// matched (MatchedAA), every mass merged into such a node included (e.g. the b and c ion of one position, or a suffix mass moved to the
  /// prefix side by the precursor mass); the precursor node is no mass. Off by default; it changes no path, score or other meta value.
  void setReportMatchedMasses(bool report)
  {
    report_matched_masses_ = report;
  }

  /// Path alternatives. The extension keeps one path per candidate: in each mode, the best score, then the smaller modification mass, the
  /// longer protein and fewer modified residues. Among equal or close scores this choice ignores the raw spectrum (a native OMPF spectrum:
  /// two N-terminal paths of equal score, 1-362 and 205-362; the smaller shift chose 205-362, whose raw evidence was the weaker). With
  /// @p scorer (a hit's raw-evidence score on the scale of the extension score) and @p max_alternatives > 0, a candidate with distinct
  /// proteoforms (protein length, modification mass) among the paths of a mode within @p margin of that mode's best score is also extended
  /// with each of the next @p max_alternatives of them in that mode (the best one in the other modes), and the hit with the best extension
  /// score + scorer is kept (the usual one among equal ones). Targets and decoys alike. Without (the default): one path, as before.
  void setPathAlternatives(std::function<double(const ProteinHit&)> scorer, int max_alternatives, int margin)
  {
    alternative_scorer_ = std::move(scorer);
    max_path_alternatives_ = alternative_scorer_ ? std::max(0, max_alternatives) : 0;
    path_alternative_margin_ = std::max(0, margin);
  }

protected:
  void updateMembers_() override;
  /// implemented for DefaultParamHandler
  void setDefaultParams_();

private:
  /// The best score per (zero-score vertex, stretch key) while one DAG is built: lookups and stores only, on one flat open-addressing
  /// table (a map of hash maps allocated a node per store, and its stores took more than half of the extension time)
  class NodeScoreMap
  {
  public:
    /// the stored score of (@p vertex, @p key), or nullptr
    const int* find(Size vertex, Size key) const
    {
      if (table_.empty()) return nullptr;
      for (Size i = slot_(vertex, key); ; i = (i + 1) & mask_)
      {
        const auto& e = table_[i];
        if (e.vertex == empty_) return nullptr;
        if (e.vertex == vertex && e.key == key) return &e.score;
      }
    }

    /// store (or overwrite) the score of (@p vertex, @p key)
    void set(Size vertex, Size key, int score)
    {
      if ((count_ + 1) * 2 > table_.size()) grow_();
      for (Size i = slot_(vertex, key); ; i = (i + 1) & mask_)
      {
        auto& e = table_[i];
        if (e.vertex == empty_)
        {
          e = {vertex, key, score};
          ++count_;
          return;
        }
        if (e.vertex == vertex && e.key == key)
        {
          e.score = score;
          return;
        }
      }
    }

  private:
    struct Entry
    {
      Size vertex;
      Size key;
      int score;
    };
    static constexpr Size empty_ = std::numeric_limits<Size>::max(); // no vertex of a DAG (vertex < DAG size)
    std::vector<Entry> table_;
    Size mask_ = 0;
    Size count_ = 0;

    Size slot_(Size vertex, Size key) const
    {
      Size h = vertex * 0x9E3779B97F4A7C15ULL ^ (key + 0x632BE59BD9B4E019ULL) * 0xC2B2AE3D27D4EB4FULL;
      return (h ^ (h >> 29)) & mask_;
    }

    void grow_()
    {
      std::vector<Entry> old;
      old.swap(table_);
      table_.assign(old.empty() ? 1024 : old.size() * 2, Entry {empty_, 0, 0});
      mask_ = table_.size() - 1;
      count_ = 0;
      for (const auto& e : old)
      {
        if (e.vertex != empty_) set(e.vertex, e.key, e.score);
      }
    }
  };

  /// Protein hit information that should be retained throughout the run.
  struct OPENMS_DLLAPI HitInformation
  {
  public:
    /// vertices already reached when building the DAG. A set, not a bitset over the whole vertex space: that space (spectrum nodes x
    /// protein positions x modification counts x path scores) takes tens to hundreds of MB per hit while only ~10^2-10^4 vertices are
    /// reached.
    FLASHHelperClasses::VisitedVertexSet visited_;
    /// The DAG consisting of paths representing proteoforms
    FLASHHelperClasses::DAG dag_;
    /// We perform three times of extensions - for suffix, for prefix, and for both termini
    int mode_;
    /// For each mode, information on each mass should be kept. These retain the information of each mass
    std::map<int, MSSpectrum> node_spec_map_, tol_spec_map_;
    /// For each mode, protein suffix or prefix masses are stored here.
    std::map<int, std::vector<double>> pro_mass_map_;
    /// protein start and end position for protein truncation, calculated after extension (after proteoform characterization).
    int protein_start_position_ = -1, protein_end_position_ = -1;
    /// calculated precursor mass from fragment mass pairing (representing complementary ion pairs)
    double calculated_precursor_mass_ = -1;
    /// the precursor mass is taken as given (no isotope offsets tried): set when an inferred mass gave no path and the MS1 mass is tried
    bool precursor_as_given_ = false;
    /// the suffix and prefix node masses (residue sums) whose sum defined calculated_precursor_mass_ when it was inferred from fragments; -1 if not
    double inferred_suffix_node_ = -1, inferred_prefix_node_ = -1;
    /// in mode 2, the start position of the terminal (prefix) result bounds the start of the precursor-mass path only when unshifted prefix
    /// masses read from that start support it (see min_n_term_ladder)
    bool start_prior_ = true;
    /// n and c term best mass shift between spec and protein. For speed up
    int n_best_shift_ = 0, c_best_shift_ = 0;
    /// cached pointers to the current-mode (mode_) entries of the maps above; set once per mode to avoid re-resolving the std::map lookup on every recursion frame.
    const std::vector<double>* cur_pro_masses_ = nullptr;
    const MSSpectrum* cur_node_spec_ = nullptr;
    const MSSpectrum* cur_tol_spec_ = nullptr;
    /// the protein sequence in the orientation of the current mode's protein masses (reversed for mode 0)
    String cur_oriented_seq_;
    /// per mode, the fragment series of each node: 0 a/b/x/y, 1 c/z, -1 node 0, the precursor or a merged node of both
    std::map<int, std::vector<int>> node_series_map_;
    /// with report_matched_masses_, per mode: the deconvolved-spectrum index of each mass of the mode before merging (node_source_map_, in
    /// mass order) and, per node of node_spec_map_, the range [first, last] of those masses merged into it (first > last: node 0 and the
    /// precursor node, which are no masses)
    std::map<int, std::vector<int>> node_source_map_;
    std::map<int, std::vector<std::pair<int, int>>> node_range_map_;
    const std::vector<int>* cur_node_series_ = nullptr;
    /// per node k of the current mode, the smallest mass and the largest tolerance of the nodes k, k + 1, ... (index node count: none)
    std::vector<double> cur_tail_min_mass_, cur_tail_max_margin_;
    /// the protein masses of the current mode by 1 Da (bit b: a mass in [b, b + 1)), for the node filter of diagonal steps; empty (no
    /// filter) unless the masses ascend within [0, 2^24). The filter's rounding slack: 1e-9 of the largest magnitudes of the mode.
    std::vector<uint64_t> cur_mass_bits_;
    double cur_mass_slack_ = 0;
    /// per series, the cleavage weight of a matched mass at each protein index of the current mode (index k: the bond after k residues)
    std::vector<int> cur_bond_weight_[2];
    /// the fixed modifications of the hit's entry (setFixedModifications; nullptr: none), and their kinds per residue (bits 1 side chain,
    /// 2 N-terminal, 4 C-terminal) in the orientation of cur_oriented_seq_ (empty: none)
    const std::vector<std::pair<int, double>>* fixed_ = nullptr;
    const std::vector<std::pair<int, std::uint8_t>>* occupancy_ = nullptr;
    std::vector<std::uint8_t> cur_occupied_;
  };

  /// a variable modification: a residue-specific mass shift searched with its own count, independently of blind modifications
  struct VariableModification
  {
    String name;
    double mass = 0;
    String residues;       ///< one-letter codes of the residues it can sit on; empty or "*" for any residue
    int position = 0;      ///< 0 anywhere, 1 proteoform N-terminal residue, 2 proteoform C-terminal residue
    int unimod = -1;
    int sites = 1;         ///< residues it needs (2 for a residue written twice, e.g. "CC" for a disulfide)
  };

  /// one or more variable modifications explaining one mass shift between two matched masses
  struct VariableModificationCombination
  {
    double mass = 0;
    std::vector<int> mods; ///< indices into var_mod_list_, ascending
  };

  /// the candidate blind modifications (setCandidateBlindModificationMap) in flat arrays: the masses of the map (ascending) and, per mass k,
  /// the UniMod record ids and ids of its modifications in the map's order at [blind_mod_offsets_[k], blind_mod_offsets_[k + 1]). The
  /// map's nodes were looked up in the innermost loop of the extension, and copied for every extender.
  std::vector<double> blind_mod_masses_;
  std::vector<Size> blind_mod_offsets_;
  std::vector<int> blind_mod_unimods_;
  std::vector<String> blind_mod_ids_;
  /// per bucket of 1 / blind_mod_bucket_scale_ Da from blind_mod_masses_.front(), the first mass index in it: where blindModLowerBound_ starts
  std::vector<UInt> blind_mod_buckets_;
  double blind_mod_bucket_scale_ = 0;
  /// the index of the first mass of blind_mod_masses_ not below @p mass, as std::lower_bound (and std::map::lower_bound) give it
  Size blindModLowerBound_(double mass) const;
  /// blind_mod_masses_ followed by four +inf (blindModNear_ compares four masses from a bucket's first one without bounds checks)
  std::vector<double> blind_mod_padded_;
  /// whether the first mass of blind_mod_masses_ not below @p delta - @p tol exists and lies within @p tol of @p delta (the blind-shift
  /// test of findSubPathsBetweenTagEndPoints, called for most positions there): the result of blindModLowerBound_ and that test. Kept
  /// inline in that large function (the compiler called it otherwise).
#if defined(__GNUC__)
  __attribute__((always_inline))
#endif
  bool blindModNear_(double delta, double tol) const
  {
    const double mass = delta - tol;
    const auto& masses = blind_mod_masses_;
    Size i;
    if (blind_mod_buckets_.empty()) i = blindModLowerBound_(mass);
    else if (! (mass > masses.front())) i = 0; // no mass below it (or NaN): std::lower_bound gives the first
    else
    {
      // The masses before the first one of the bucket of `mass` are below it and those of later buckets are not (see
      // blindModLowerBound_), so the masses below it from that first one on are the first `below` of the next four, counted without
      // branches; only a crowded bucket can hold more.
      const double cell = (mass - masses.front()) * blind_mod_bucket_scale_;
      i = cell < (double)(blind_mod_buckets_.size() - 1) ? blind_mod_buckets_[(Size)cell] : blind_mod_buckets_.back();
      const double* m = blind_mod_padded_.data() + i;
      const Size below = (Size)(m[0] < mass) + (Size)(m[1] < mass) + (Size)(m[2] < mass) + (Size)(m[3] < mass);
      i += below;
      if (below == 4)
      {
        while (blind_mod_padded_[i] < mass) i++; // stops at the first +inf (index blind_mod_masses_.size()) at the latest
      }
    }
    return i < masses.size() && ! (std::abs(delta - blind_mod_padded_[i]) > tol);
  }
  std::vector<VariableModification> var_mod_list_;
  /// score penalty of each modification of var_mod_list_
  std::vector<int> var_mod_penalties_;
  /// all combinations of up to max_var_mod_count variable modifications, sorted by mass
  std::vector<VariableModificationCombination> var_mod_combinations_;

  /// read var_mod (modification names) and var_mod_file (TopPIC format: name,mass,residues,position,unimod)
  void loadVariableModifications_();

  /// can the modifications of @p combination sit on distinct residues among @p seq[from, to) (current mode orientation)?
  /// @p n_site / @p c_site: the index in @p seq of the proteoform's N- / C-terminal residue if this stretch can reach it, else -1;
  /// @p occupied: per residue of @p seq, the kinds of its fixed modifications (HitInformation::cur_occupied_), nullptr: none
  bool combinationFits_(const VariableModificationCombination& combination, const String& seq, int from, int to, int n_site, int c_site,
                        const std::uint8_t* occupied = nullptr) const;

  /// names of the variable modification combinations of mass @p mass within @p tol that fit @p seq[from, to); @p site is the only
  /// possible residue (index in @p seq) if a single one-site modification fits at a single residue, else -1
  String variableModificationNames_(double mass, double tol, const String& seq, int from, int to, int n_site, int c_site, String& accessions,
                                    int& site, const std::uint8_t* occupied = nullptr) const;

  /// whether a variable modification at @p position (0 anywhere, 1 N-terminal, 2 C-terminal) may sit on a residue with fixed kinds @p bits
  static bool freeFor_(std::uint8_t bits, int position)
  {
    return (bits & (position == 0 ? 1 : (position == 1 ? 2 : 4))) == 0;
  }

  /// the residues that can carry the shift between two matched masses of a path, (@p pre_node, @p pre_pro) and (@p node, @p pro), whose
  /// summed shift at (@p node, @p pro) is @p mass_shift: [@p pre_pro, returned end) in this mode's orientation, where the end is the first
  /// protein index in (@p pre_pro, @p pro) that a mass in between matches with the whole shift, else @p pro. Used for the DAG edge and for
  /// the report of the path, so that a variable modification is tested and named on the same residues.
  int shiftRangeEnd_(const MSSpectrum& node_spec, const MSSpectrum& tol_spec, const std::vector<double>& pro_masses, int pre_node, int pre_pro,
                     int node, int pro, double mass_shift) const;

  /// largest |mass| that up to @p count variable modifications can explain (0 without them)
  double varModRoom_(int count) const
  {
    if (count <= 0 || var_mod_room_.empty()) return 0;
    return var_mod_room_[std::min(count, (int)var_mod_room_.size() - 1)];
  }
  std::vector<double> var_mod_room_;

  /** get protein prefix or suffix masses (depending on the mode)
   * @param hit protein hit
   * @param pro_masses protein prefix or suffix masses to be filled
   * @param mode mode
   */
  static void getProMasses_(const ProteinHit& hit, std::vector<double>& pro_masses, int mode,
                            const std::vector<std::pair<int, double>>* fixed = nullptr);

  /***
   * calculate precursor mass with fragment mass pairs
   * @param hit protein hit
   * @param best_path_map best path(s) per mode (key: mode value: paths)
   * @param hi protein hit information
   */
  /// true if more complementary fragment pairs (prefix + suffix mass = @p precursor_mass, for the ion types in use) are found in
  /// @p sorted_masses than at any of 200 masses moved by up to +-300 Da (chance in this spectrum)
  bool hasComplementaryPairSupport_(const std::vector<double>& sorted_masses, double precursor_mass, double suffix_node = -1, double prefix_node = -1) const;

  void calculatePrecursorMass_(const ProteinHit& hit,
                               const std::map<int, std::vector<Size>>& best_path_map,
                               HitInformation& hi);

  /**
   * update HitInformation
   * @param dspec deconvolved spectrum
   * @param hi the HitInformation to be updated
   * @param max_mass maximum mass to consider
   */
  void updateHitInformation_(const DeconvolvedSpectrum& dspec, HitInformation& hi,
                   double max_mass);

  static int getBestMassShift_(const std::vector<int>& spec_vec, const std::unordered_set<int>& pro_vec);

  /// of the placements @p positions of @p tag on @p seq (by its letters), those whose residues, with the fixed modifications @p fixed
  /// (position in @p seq, mass), have the mass span of the tag's masses; a tag with lowercase letters is only kept off fixed residues
  void keepFixedConsistentPlacements_(const String& seq, const FLASHHelperClasses::Tag& tag, const std::vector<std::pair<int, double>>& fixed,
                                      std::vector<int>& positions) const;

  /**
   * The main function to perform extension algorithm for each hit.
   * @param hit the target protein hit
   * @param hi the HitInformation of the hit
   * @param matched_tags matched tags
   * @param tag_positions per matched tag, whether its protein positions are found yet and the positions (the same in every mode)
   * @param all_paths_per_mode all found paths per mode
   * @param max_mod_cntr_for_last_mode the maximum number of modifications for the last mode (mode == 2)
   * @param alternative per blind modification count, the rank of the distinct proteoform taken (0: the usual path; a count with fewer
   * distinct ones takes its last); see setPathAlternatives
   * @param alternatives set to the most distinct proteoforms of any blind modification count among the ranked paths
   */
  void run_(const ProteinHit& hit, HitInformation& hi,
            const std::vector<FLASHHelperClasses::Tag>& matched_tags, std::vector<std::pair<bool, std::vector<int>>>& tag_positions,
            std::map<int, std::vector<Size>>& all_paths_per_mode, int max_mod_cntr_for_last_mode, int alternative, int& alternatives); // per hit

  /**
   * get the vertex number from the values describing the position on the multidimensional graph for extension algorithm
   * @param node_index the index representing a deconvolved mass
   * @param pro_index the index respresnting a fragment mass
   * @param score the score of the node
   * @param num_blind_mod the number of blind modifications
   * @param num_var_mod the number of variable modifications
   * @param pro_mass_size the size of protein fragment masses
   * @return
   */
  Size getVertex_(int node_index, int pro_index, int score, int num_blind_mod, int num_var_mod, Size pro_mass_size) const;

  /// the inverse function of getVertex_(...). Get the values describing the position from the vertex number
  int getNodeIndex_(Size vertex, Size pro_mass_size) const;
  int getProIndex_(Size vertex, Size pro_mass_size) const;
  int getBlindModNumber_(Size vertex) const;
  int getVarModNumber_(Size vertex) const;
  int getScore_(Size vertex) const;

  /**
   * construct DAG on which the extension algorithm run.
   * @param sinks The sinks of the DAG
   * @param hi HitInformation
   * @param tag_edges edges corresponding to the matching tags
   * @param max_mod_cntr_for_last_mode maximum number of modifications for the last mode (mode == 2)
   * @param use_tags to use tags or not.
   */
  void constructDAG_(std::set<Size>& sinks,
                     HitInformation& hi,
                     const std::vector<std::vector<int>>& tag_edges,
                     int max_mod_cntr_for_last_mode,
                     bool use_tags);
  /**
   * a recursive function to find subpaths along vertices representing tag end points.
   * @param visited_tag_edges already visited tag edges
   * @param hi HitInformation
   * @param sinks the sinks of the DAG
   * @param vertex the current vertex
   * @param truncation_mass the proteoform truncation mass
   * @param cumulative_shift the summed mass shift of the modifications so far
   * @param node_max_score_map the  maximum score achieving vertex of the node
   * @param tag_edges edges corresponding to the matching tags
   * @param max_mod_cntr_for_last_mode maximum number of modifications for the last mode (mode == 2)
   * @param use_tags to use tags or not.
   */
  void findPathsAlongTagEndPoints(std::set<Size>& visited_tag_edges,
                           HitInformation& hi,
                           std::map<Size, std::set<std::pair<double, double>>>& sinks,
                           Size vertex,
                           double truncation_mass,
                           double cumulative_shift,
                           NodeScoreMap& node_max_score_map,
                           const std::vector<std::vector<int>>& tag_edges,
                           int max_mod_cntr_for_last_mode,
                           bool use_tags);

  /**
   * a recursive function to find subpaths between tag end points. It runs within findPathsAlongTagEndPoints(..).
   * @param sinks the sinks of the DAG
   * @param hi HitInformation
   * @param start_vertex the starting vertex (representing a tag end point)
   * @param end_node_index the node of the end vertex (representing a tag end point)
   * @param end_pro_index the protein position of the end vertex
   * @param diagonal_counter how many diagonal subpaths (representing non modified a.a. masses) are present so for?
   * @param truncation_mass the proteoform truncation mass
   * @param cumulative_mod_mass the summed mass shift of the modifications so far
   * @param node_max_score_map the  maximum score achieving vertex of the node
   * @param max_blind_mod_cntr_for_last_mode maximum number of modifications for the last mode (mode == 2)
   */
  void findSubPathsBetweenTagEndPoints(std::map<Size, std::set<std::pair<double, double>>>& sinks,
                          HitInformation& hi,
                          Size start_vertex,
                          int end_node_index,
                          int end_pro_index,
                          int diagonal_counter,
                          double truncation_mass,
                          double cumulative_mod_mass,
                          NodeScoreMap& node_max_score_map,
                          int max_blind_mod_cntr_for_last_mode);

  /// get proteoform specific information from the input @p path and relevant data.
  int getProteinLength_(const std::vector<Size>& path, const std::vector<double>& pro_masses) const;
  double getSpecMassSpan_(const std::vector<Size>& path, const MSSpectrum& node_spec, int pro_mass_size) const;
  double getProteinMassSpan_(const std::vector<Size>& path, const std::vector<double>& pro_masses) const;
  int getModifiedAACount_(const std::vector<Size>& path) const;

  std::vector<std::string> ion_types_str_;
  std::vector<double> prefix_shifts_;
  std::vector<int> prefix_series_, suffix_series_; ///< 0 a/b/x/y, 1 c/z, index-parallel to the shifts

  /// log2(cleavage rate / mean) of a bond by the residue before and after it, per series (0 b/y, 1 c/z), indexed by residue - 'A'
  double cleavage_before_[2][26] = {}, cleavage_after_[2][26] = {};
  double cleavage_weight_ = 0;
  /// fill the built-in cleavage tables, or read cleavage_table (series<TAB>before|after<TAB>residue<TAB>log2 ratio)
  void loadCleavageTables_();
  std::vector<double> suffix_shifts_;
  std::vector<ProteinHit> proteoform_hits_;
  double tol_;

  bool skip_precursor_inference_ = true;
  /// set MatchedMassIndices on every proteoform hit (setReportMatchedMasses)
  bool report_matched_masses_ = false;
  /// setPathAlternatives
  std::function<double(const ProteinHit&)> alternative_scorer_;
  int max_path_alternatives_ = 0;
  int path_alternative_margin_ = 0;
  /// per database entry, its fixed modifications (position, mass) and their kinds (setFixedModifications); nullptr: none
  const std::vector<std::vector<std::pair<int, double>>>* fixed_masses_ = nullptr;
  const std::vector<std::vector<std::pair<int, std::uint8_t>>>* occupancy_ = nullptr;
  int max_blind_mod_cntr_ = 1;
  int max_var_mod_cntr_ = 0;
  int allowed_isotope_error_ = 1;
  const int max_path_score_ = 1200;
  const int min_path_score_ = -20;
  const int max_extension_stretch_ = 50;
  double max_mod_mass_ = 500.0;
  /// largest absolute mass of a single variable modification (0 without them)
  double max_var_mod_mass_ = 0;
  double given_precursor_mass_ = -1;

  /// precomputed vertex-packing spans (set in updateMembers_): var, score, blind, and their products. Invariant during run().
  Size var_span_ = 1, score_span_ = 1221, blind_span_ = 2, var_score_span_ = 1221, var_score_blind_span_ = 2442;
};
} // namespace OpenMS