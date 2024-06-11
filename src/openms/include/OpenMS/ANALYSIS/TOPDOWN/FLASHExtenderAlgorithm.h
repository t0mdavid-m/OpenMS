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
    FLASHExtenderAlgorithm(const FLASHExtenderAlgorithm&);

    /// move constructor
    FLASHExtenderAlgorithm(FLASHExtenderAlgorithm&& other) = default;

    /// assignment operator
    FLASHExtenderAlgorithm& operator=(const FLASHExtenderAlgorithm& other);

  protected:
    void updateMembers_() override;
    /// implemented for DefaultParamHandler
    void setDefaultParams_();
  private:

  };
} // namespace OpenMS