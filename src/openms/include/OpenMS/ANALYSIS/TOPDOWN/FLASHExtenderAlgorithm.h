// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong$
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
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

  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();
  private:

    int getVertex_(int peak_index, int pro_index, int score, int num_mod, int layer) const;
    int getPeakIndex_(int vertex) const;
    int getProIndex_(int vertex) const;
    int getModNumber_(int vertex) const;
    int getScore_(int vertex) const;
    void constructSubDAG_(FLASHHelperClasses::DAG& dag, boost::dynamic_bitset<>& visited,
                          int vertex1, int peak_index2, int pro_index2, int layer, bool allow_truncation = false);

    std::vector<int> peak_scores_;
    std::vector<int> peak_masses_;
    std::vector<int> pro_masses_;

    double tol_;
    int pro_length_ = 0;
    int max_mod_cntr_ = 0;
    int max_iso_in_tag_ = 0;
    int max_path_score_ = 0;
    int min_path_score_ = 0;
    int max_layer_ = 0;
    double fdr_ = 1.0;
    double max_edge_mass_ = 0;
    double max_mod_mass_ = 500.0;
  };
} // namespace OpenMS