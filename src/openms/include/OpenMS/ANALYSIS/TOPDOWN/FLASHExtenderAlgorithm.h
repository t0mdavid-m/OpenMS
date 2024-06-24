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

    void run(const FLASHTaggerAlgorithm& tagger);

  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();
  private:
    void get_pro_masses_(const ProteinHit& hit, std::vector<double>& pro_masses, int mode);
    void calcualte_precursor_mass_(const std::vector<MSSpectrum>& node_spectrum_list,
                                  const std::vector<std::vector<double>>& pro_masses_list, const std::map<int, std::vector<Size>>& best_paths);
    void define_nodes_(const FLASHTaggerAlgorithm& tagger, MSSpectrum& node_spec, MSSpectrum& tol_spec, double max_mass, int mode);
    void run_(const FLASHTaggerAlgorithm& tagger, const ProteinHit& hit,
              const MSSpectrum& node_spec, const MSSpectrum& tol_spec, const std::vector<double>& pro_masses,
              std::vector<std::vector<Size>>& all_paths, int mode); // per hit
    Size getVertex_(int node_index, int pro_index, int score, int num_mod, Size pro_length) const;
    int getNodeIndex_(Size vertex, Size pro_length) const;
    int getProIndex_(Size vertex, Size pro_length) const;
    int getModNumber_(Size vertex) const;
    int getScore_(Size vertex) const;
    void constructDAG_(FLASHHelperClasses::DAG& dag, std::set<Size>& sinks, const MSSpectrum& node_spec, const MSSpectrum& tol_spec, const std::vector<double>& pro_masses, const std::vector<std::vector<int>>& tag_edges, int mode);
    void connectBetweenTags_(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited, std::set<Size>& visited_tag_edges, std::set<Size>& sinks, Size vertex, const MSSpectrum& node_spec, const MSSpectrum& tol_spec, const std::vector<double>& pro_masses, const std::vector<std::vector<int>>& tag_edges, int mode);
    void extendBetweenTags_(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited, std::set<Size>& sinks,
                           Size vertex, int node_index, int pro_index, int diagonal_counter, const MSSpectrum& node_spec, const MSSpectrum& tol_spec, const std::vector<double>& pro_masses, int mode);

    std::vector<String> ion_types_str_;
    std::vector<double> prefix_shifts_;
    std::vector<double> suffix_shifts_;

    double tol_;
    double precursor_mass_ = .0;
    int max_mod_cntr_ = 0;
    int max_path_score_ = 0;
    int min_path_score_ = 0;
    double fdr_ = 1.0;
    double max_mod_mass_ = 500.0;
  };
} // namespace OpenMS