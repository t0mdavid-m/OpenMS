// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteomeSearch.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <map>
#include <set>
#include <utility>

namespace OpenMS
{
  Param ProteomeSearch::defaultParameters()
  {
    Param p = FLASHTnTAlgorithm().getDefaults();
    p.setValue("one_proteoform_per_feature", "false");
    p.setValue("ex:skip_precursor_inference", "true");
    return p;
  }

  ProteomeSearch::ProteomeSearch(std::vector<FASTAFile::FASTAEntry> reference, const Param& tnt_params) : params_(tnt_params)
  {
    reference_.reserve(reference.size());
    for (auto& entry : reference)
    {
      if (String(entry.identifier).hasPrefix("DECOY"))  // FLASHTnT's decoy rule (FLASHTnTAlgorithm::run, ProteoformDatabase)
      {
        decoys_dropped_++;
        continue;
      }
      reference_.push_back(std::move(entry));
    }
    FLASHTnTAlgorithm builder;
    builder.setParameters(params_);
    prepared_ = builder.prepareDatabase(reference_);
  }

  ProteomeSearch::Worker::Worker(const ProteomeSearch& search) : search_(search)
  {
    // constructed and given the parameters, never copied: FLASHTnTAlgorithm's copy constructor copies only its Param
    tnt_.setParameters(search_.params_);
    tnt_.setPreparedDatabase(search_.prepared_);
  }

  std::vector<ProteinHit> ProteomeSearch::Worker::search(const DeconvolvedSpectrum& ms2, double tol_ppm, const MSSpectrum& raw, int scan)
  {
    if (ms2.empty()) return {};  // toSpectrum reads the first mass
    DeconvolvedSpectrum sorted = ms2;  // toSpectrum is non-const, and sorts its peaks but not its per-mass lists
    sorted.sort();
    MSSpectrum spec = sorted.toSpectrum(0, tol_ppm, false);
    spec.setNativeID("controllerType=0 controllerNumber=1 scan=" + String(scan));  // FLASHTnT reads the scan from the native id
    spec.setMSLevel(2);
    spec.setRT(raw.getRT());
    MSExperiment one;
    one.addSpectrum(spec);
    tnt_.setRawSpectrumLoader([&raw, scan](const std::set<int>& scans, std::map<int, MSSpectrum>& out) {
      if (scans.count(scan) > 0) out[scan] = raw;
    });
    tnt_.run(one, search_.reference_);
    std::vector<ProteinHit> hits = tnt_.getProteoforms();
    tnt_.releaseResults();  // run() keeps its results (and tags) until released; the next search must start from nothing
    tnt_.setRawSpectrumLoader(nullptr);
    return hits;
  }
} // namespace OpenMS
