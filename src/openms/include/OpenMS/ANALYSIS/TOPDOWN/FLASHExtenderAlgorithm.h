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

    void run(const DeconvolvedSpectrum& dspec, const FLASHTaggerAlgorithm& tagger);

  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();
  private:

    int getVertex_(int node_index, int pro_index, int score, int num_mod) const;
    int getNodeIndex_(int vertex) const;
    int getProIndex_(int vertex) const;
    int getModNumber_(int vertex) const;
    int getScore_(int vertex) const;
    void constructDAG_(FLASHHelperClasses::DAG& dag, std::vector<int>& sinks, std::vector<int>& tag_node_starts, std::vector<int>& tag_pro_starts, std::vector<int>& tag_node_ends, std::vector<int>& tag_pro_ends);
    void connectBetweenTags(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited, std::vector<int>& sinks, int vertex, std::vector<int>& tag_node_starts, std::vector<int>& tag_pro_starts, std::vector<int>& tag_node_ends, std::vector<int>& tag_pro_ends);
    void extendBetweenTags(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited, std::vector<int>& sinks,
                          int vertex, int node_index, int pro_index, bool go_diagonal = false);

    std::map<String, Residue::ResidueType> res_type_names_;
    std::vector<int> node_scores_;
    std::vector<double> node_masses_;
    std::vector<double> pro_masses_;
    std::vector<String> ion_types_str_;

    double tol_;
    double precursor_mass_ = .0;
    int pro_length_ = 0;
    int max_mod_cntr_ = 0;
    int max_path_score_ = 0;
    int min_path_score_ = 0;
    double fdr_ = 1.0;
    double max_mod_mass_ = 500.0;
  };
} // namespace OpenMS