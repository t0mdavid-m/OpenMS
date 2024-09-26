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
#include <iomanip>
#include <iostream>

namespace OpenMS
{
/**
@brief
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

  void run(std::vector<ProteinHit>& hits, const std::vector<FLASHHelperClasses::Tag>& tags,
           const DeconvolvedSpectrum& dspec,
           const MSSpectrum& spec,
           double flanking_mass_tol, double ppm, bool multiple_hits_per_spec);

  void getProteoforms(std::vector<ProteinHit>& hits) const
  {
    for (const auto& hit : proteoform_hits_)
    {
      hits.push_back(hit);
    }
  }

  bool hasProteoforms() const
  {
    return ! proteoform_hits_.empty();
  }

  void setModificationMap(const std::map<double, std::vector<ResidueModification>>& mod_map)
  {
    mod_map_ = mod_map;
  }

  const static int multi_ion_score = 1;

protected:
  void updateMembers_() override;
  /// implemented for DefaultParamHandler
  void setDefaultParams_();

private:
  std::map<double, std::vector<ResidueModification>> mod_map_; // modification mass to modification index. To use find nearest function


  static void getProMasses(const ProteinHit& hit, std::vector<double>& pro_masses, int mode);
  static double calculatePrecursorMass_(const ProteinHit& hit,
                                          int protein_start_position,
                                          int protein_end_position,
                                          const std::vector<int>& mod_starts,
                                          const std::vector<int>& mod_ends,
                                          const std::vector<double>& mod_masses,
                                          double& total_mod_mass);
  void defineNodes(const DeconvolvedSpectrum& dspec,//const MSSpectrum& spec,
                   MSSpectrum& node_spec, MSSpectrum& tol_spec, double max_mass, double precursor_mass, int mode);
  void run_(const ProteinHit& hit,
            const std::vector<FLASHHelperClasses::Tag>& matched_tags,
            const MSSpectrum& node_spec,
            const MSSpectrum& tol_spec,
            const std::vector<double>& pro_masses,
            boost::dynamic_bitset<>& visited,
            const double precursor_mass,
            const double total_mod_mass,
            std::map<int, std::vector<Size>>& all_paths_per_mode, int max_mod_cntr_for_last_mode,
            int mode); // per hit
  Size getVertex_(int node_index, int pro_index, int score, int num_mod, Size pro_length) const;
  int getNodeIndex_(Size vertex, Size pro_length) const;
  int getProIndex_(Size vertex, Size pro_length) const;
  int getModNumber_(Size vertex) const;
  int getScore_(Size vertex) const;
  void constructDAG_(FLASHHelperClasses::DAG& dag,
                     std::set<Size>& sinks,
                     boost::dynamic_bitset<>& visited,
                     const MSSpectrum& node_spec,
                     const MSSpectrum& tol_spec,
                     const std::vector<double>& pro_masses,
                     const std::vector<std::vector<int>>& tag_edges,
                     const double total_mod_mass,
                     int max_mod_cntr_for_last_mode,
                     int mode,
                     bool use_tags);
  void connectBetweenTags_(FLASHHelperClasses::DAG& dag,
                           boost::dynamic_bitset<>& visited,
                           std::set<Size>& visited_tag_edges,
                             std::map<Size, double>& sinks,
                           Size vertex,
                           double cumulative_shift,
                           const MSSpectrum& node_spec,
                           const MSSpectrum& tol_spec,
                           const std::vector<double>& pro_masses,
                           const std::vector<std::vector<int>>& tag_edges,
                           const double total_mod_mass,
                           int max_mod_cntr_for_last_mode,
                           int mode,
                           bool use_tags);
  void extendBetweenTags_(FLASHHelperClasses::DAG& dag,
                          boost::dynamic_bitset<>& visited,
                          std::map<Size, double>& sinks,
                          Size start_vertex,
                          int end_node_index,
                          int end_pro_index,
                          int diagonal_counter,
                          double cumulative_shift,
                          const MSSpectrum& node_spec,
                          const MSSpectrum& tol_spec,
                          const std::vector<double>& pro_masses,
                          const double total_mod_mass,
                          int max_mod_cntr_for_last_mode,
                          int mode);
  Size getLength_(const std::vector<Size>& path, const std::vector<double>& pro_masses) const;
  Size getProteinSpan_(const std::vector<Size>& path, const std::vector<double>& pro_masses) const;

  std::vector<std::string> ion_types_str_;
  std::vector<double> prefix_shifts_;
  std::vector<double> suffix_shifts_;
  std::vector<ProteinHit> proteoform_hits_;
  std::vector<FLASHHelperClasses::Tag> tags_;
  double tol_, flanking_mass_tol_;
  int max_mod_cntr_ = 0;
  const int max_path_score_ = 400;
  const int min_path_score_ = -10;
  double max_mod_mass_ = 500.0;
  double precursor_mass_ = -1;
  //bool tmp = true;
};
} // namespace OpenMS