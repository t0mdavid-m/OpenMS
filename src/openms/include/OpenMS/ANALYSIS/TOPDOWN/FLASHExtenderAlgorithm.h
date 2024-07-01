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

    void run(const FLASHTaggerAlgorithm& tagger, double flanking_mass_tol, double ppm);

    void getProteoforms(std::vector<ProteinHit>& hits) const
    {
      for (const auto& hit : proteoform_hits_)
      {
        hits.push_back(hit);
      }
    }

  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();
  private:
    static void get_pro_masses_(const ProteinHit& hit, std::vector<double>& pro_masses, int mode);
    static double calculate_precursor_mass_(const ProteinHit& hit, int protein_start_position, int protein_end_position, const std::vector<int>& mod_starts, const std::vector<int>& mod_ends, const std::vector<double>& mod_masses) ;
    void define_nodes_(const FLASHTaggerAlgorithm& tagger, MSSpectrum& node_spec, MSSpectrum& tol_spec, double max_mass, double precursor_mass, int mode);
    void run_(const ProteinHit& hit,
              const MSSpectrum& node_spec, const MSSpectrum& tol_spec, const std::vector<double>& pro_masses,
              boost::dynamic_bitset<>& visited, double precursor_mass,
              std::map<int, std::vector<Size>>& all_paths_per_mode, int mode); // per hit
    Size getVertex_(int node_index, int pro_index, int score, int num_mod, Size pro_length) const;
    int getNodeIndex_(Size vertex, Size pro_length) const;
    int getProIndex_(Size vertex, Size pro_length) const;
    int getModNumber_(Size vertex) const;
    int getScore_(Size vertex) const;
    void constructDAG_(FLASHHelperClasses::DAG& dag, std::set<Size>& sinks, boost::dynamic_bitset<>& visited, const MSSpectrum& node_spec, const MSSpectrum& tol_spec, const std::vector<double>& pro_masses, const std::vector<std::vector<int>>& tag_edges, int mode);
    void connectBetweenTags_(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited, std::set<Size>& visited_tag_edges, std::set<Size>& sinks, Size vertex, const MSSpectrum& node_spec, const MSSpectrum& tol_spec, const std::vector<double>& pro_masses, const std::vector<std::vector<int>>& tag_edges, int mode);
    void extendBetweenTags_(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited, std::set<Size>& sinks,
                           Size vertex, int node_index, int pro_index, int diagonal_counter, const MSSpectrum& node_spec, const MSSpectrum& tol_spec, const std::vector<double>& pro_masses, int mode);

    std::vector<String> ion_types_str_;
    std::vector<double> prefix_shifts_;
    std::vector<double> suffix_shifts_;
    std::vector<ProteinHit> proteoform_hits_;
    std::vector<FLASHHelperClasses::Tag> tags_;

    double tol_, flanking_mass_tol_;
    int max_mod_cntr_ = 0;
    const int max_path_score_ = 400;
    const int min_path_score_ = -20;
    // double fdr_ = 1.0;
    // bool keep_decoy_ = false;
    double max_mod_mass_ = 500.0;
    double precursor_mass_ = -1;
  };
} // namespace OpenMS