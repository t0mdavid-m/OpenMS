// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteomeSearch.h>

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <map>
#include <set>
#include <stdexcept>
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

  const std::vector<std::pair<std::string, std::string>>& ProteomeSearch::knobNames()
  {
    // FLASHTnTAlgorithm's constructor inserts the Tagger's parameters under "tag:" and the Extender's under
    // "ex:"; max_blind_mod_count, max_mod_mass, var_mod and max_var_mod_count are its own top-level copies,
    // which it hands to the Extender.
    static const std::vector<std::pair<std::string, std::string>> names = {
      {"min_length", "tag:min_length"},
      {"max_length", "tag:max_length"},
      {"allow_gap", "tag:allow_gap"},
      {"max_aa_in_gap", "tag:max_aa_in_gap"},
      {"max_blind_mod_count", "max_blind_mod_count"},
      {"max_mod_mass", "max_mod_mass"},
      {"var_mod", "var_mod"},
      {"max_var_mod_count", "max_var_mod_count"},
      {"max_protein_hits", "max_protein_hits"},
      {"gapped_tags", "gapped_tags"},
      {"terminal_ladders", "terminal_ladders"},
      {"precursor_ladders", "precursor_ladders"},
      {"skip_precursor_inference", "ex:skip_precursor_inference"},
    };
    return names;
  }

  Param ProteomeSearch::parametersFor(const FlashTnTSettings& s)
  {
    const auto name = [](const std::string& key) -> const std::string& {
      for (const auto& kn : knobNames())
        if (kn.first == key) return kn.second;
      throw std::logic_error("ProteomeSearch::parametersFor: no FLASHTnT name for identification.flashtnt." + key);
    };
    const auto flag = [](bool b) { return std::string(b ? "true" : "false"); };

    Param p = defaultParameters();
    p.setValue(name("min_length"), s.min_length);
    p.setValue(name("max_length"), s.max_length);
    p.setValue(name("allow_gap"), flag(s.allow_gap));
    p.setValue(name("max_aa_in_gap"), s.max_aa_in_gap);
    p.setValue(name("max_blind_mod_count"), s.max_blind_mod_count);
    p.setValue(name("max_mod_mass"), s.max_mod_mass);
    p.setValue(name("var_mod"), s.var_mod);
    p.setValue(name("max_var_mod_count"), s.max_var_mod_count);
    p.setValue(name("max_protein_hits"), s.max_protein_hits);
    p.setValue(name("gapped_tags"), s.gapped_tags);
    p.setValue(name("terminal_ladders"), flag(s.terminal_ladders));
    p.setValue(name("precursor_ladders"), flag(s.precursor_ladders));
    p.setValue(name("skip_precursor_inference"), flag(s.skip_precursor_inference));
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
    // no progress lines per MS2: a prepared run() keeps this log type (FLASHTnTAlgorithm::run sets CMD only unprepared)
    tnt_.setLogType(ProgressLogger::NONE);
  }

  std::vector<ProteinHit> ProteomeSearch::Worker::search(const DeconvolvedSpectrum& ms2, double tol_ppm, const MSSpectrum& raw, int scan,
                                                         int* tag_count)
  {
    if (tag_count != nullptr) *tag_count = 0;
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
    if (tag_count != nullptr) *tag_count = static_cast<int>(tnt_.getTags().size());  // every tag this spectrum produced
    tnt_.releaseResults();  // run() keeps its results (and tags) until released; the next search must start from nothing
    tnt_.setRawSpectrumLoader(nullptr);
    return hits;
  }
} // namespace OpenMS
