// Copyright (c) 2002-2024, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong$
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <utility>

namespace OpenMS
{
inline const Size max_node_cntr = 500;

FLASHExtenderAlgorithm::FLASHExtenderAlgorithm(): DefaultParamHandler("FLASHExtenderAlgorithm"), ProgressLogger()
{
  setDefaultParams_();
}

FLASHExtenderAlgorithm::FLASHExtenderAlgorithm(const FLASHExtenderAlgorithm& other): DefaultParamHandler(other), ProgressLogger(other)
{
}

FLASHExtenderAlgorithm& FLASHExtenderAlgorithm::operator=(const FLASHExtenderAlgorithm& rhs)
{
  if (this == &rhs) return *this;

  DefaultParamHandler::operator=(rhs);
  return *this;
}

void FLASHExtenderAlgorithm::setDefaultParams_()
{
  defaults_.setValue("max_tag_count", 500,
                     "Maximum number of the tags per length (lengths set by -min_length and -max_length options). The tags with different amino acid "
                     "combinations are all treated separately. E.g., "
                     "TII, TIL, TLI, TLL are distinct tags even though they have the same mass differences. "
                     "but are counted as four different tags. ");
  defaults_.setMinInt("max_tag_count", 0);

  defaults_.setValue(
    "min_length", 4,
    "Minimum length of a tag. Each mass gap contributes to a single length (even if a mass gap is represented by multiple amino acids). ");
  defaults_.setMaxInt("min_length", 30);
  defaults_.setMinInt("min_length", 3);

  defaults_.setValue(
    "max_length", 10,
    "Maximum length of a tag. Each mass gap contributes to a single length (even if a mass gap is represented by multiple amino acids). ");
  defaults_.setMaxInt("max_length", 30);
  defaults_.setMinInt("max_length", 3);

  defaults_.setValue("flanking_mass_tol", 700.0, "Flanking mass tolerance in Da.");
  defaults_.setValue("max_iso_error_count", 0, "Maximum isotope error count per tag.");
  defaults_.setMaxInt("max_iso_error_count", 2);
  defaults_.setMinInt("max_iso_error_count", 0);
  defaults_.addTag("max_iso_error_count", "advanced");
  defaults_.setValue("min_matched_aa", 5, "Minimum number of amino acids in matched proteins, covered by tags.");

  defaults_.setValue("fdr", 1.0, "Protein FDR threshold.");
  defaults_.setMaxFloat("fdr", 1.0);
  defaults_.setMinFloat("fdr", 0.01);

  defaults_.setValue("keep_decoy", "false", "Keep decoy proteins.");
  defaults_.addTag("keep_decoy", "advanced");
  defaults_.setValidStrings("keep_decoy", {"true", "false"});

  defaultsToParam_();
}

void FLASHExtenderAlgorithm::updateMembers_()
{
}

/*
int FLASHExtenderAlgorithm::getVertex_(int index, int path_score, int level, int iso_level) const
{
  return ((index * (max_tag_length_ + 1) + level) * (max_iso_in_tag_ + 1) + iso_level) * (max_path_score_ - min_path_score_ + 1)
         + (path_score - min_path_score_);
}

int FLASHExtenderAlgorithm::getIndex_(int vertex) const
{
  return ((vertex / (max_path_score_ - min_path_score_ + 1)) / (max_iso_in_tag_ + 1)) / (max_tag_length_ + 1);
}
*/


} // namespace OpenMS