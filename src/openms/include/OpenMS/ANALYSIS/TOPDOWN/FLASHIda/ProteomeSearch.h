// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom Mueller $
// $Authors: Tom Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTnTAlgorithm.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ProteinHit.h>

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
  class DeconvolvedSpectrum;
  struct FlashTnTSettings;

  /**
   * @brief FLASHTnT's own search of one MS2 against the Protein reference (ADR-0061 d.3).
   *
   * The reference is loaded once: entries with FLASHTnT's decoy prefix (DECOY) are dropped and counted, so
   * real time never searches a decoy, and the database side of FLASHTnTAlgorithm::run is prepared once
   * (FLASHTnTAlgorithm::prepareDatabase) and shared, read-only, by every Worker.
   *
   * A Worker is one search thread's: it owns its FLASHTnTAlgorithm (so its Tagger, Extender and Annotator) and
   * searches one deconvolved MS2 with its raw peaks by calling run() on a one-spectrum input. With a target-only
   * reference every run-level step of run() (protein FDR, cascade, q-values) has nothing to decide, so a
   * spectrum's PrSMs are those a whole run of MS2s gives it; FLASHIda_ProteomeSearch_test pins that.
   *
   * The calling thread sets OpenMP's team size (a new thread does not inherit it; SearchPool's on_thread_start
   * does it for the engine).
   */
  class OPENMS_DLLAPI ProteomeSearch
  {
  public:
    /// FLASHTnTAlgorithm's defaults with the real-time diffs: one_proteoform_per_feature false (it needs MS1
    /// input), ex:skip_precursor_inference true (ADR-0061 d.13)
    static Param defaultParameters();

    /// identification.flashtnt's knobs, each as (config key, FLASHTnTAlgorithm parameter name), in the order
    /// Config declares them. The one place the mapping lives: parametersFor writes through it and the tests
    /// read it. fixed_mod maps to FLASHTnTAlgorithm's own fixed_mod (ADR-0067).
    static const std::vector<std::pair<std::string, std::string>>& knobNames();

    /// defaultParameters() with every knob of @p settings written to its FLASHTnT name (booleans as the
    /// strings "true"/"false", as FLASHTnT declares them)
    static Param parametersFor(const FlashTnTSettings& settings);

    /// @p reference: the Protein reference's entries (decoys are dropped); @p tnt_params: FLASHTnTAlgorithm's
    /// parameters (ProteomeSearch::parametersFor; fixed modifications, the iodoTMT label of quantification.scan ms3/both
    /// included, are its fixed_mod, ADR-0067).
    ProteomeSearch(std::vector<FASTAFile::FASTAEntry> reference, const Param& tnt_params);

    /// the entries dropped at load for FLASHTnT's decoy prefix
    Size decoysDropped() const
    {
      return decoys_dropped_;
    }

    /// the entries searched (decoys dropped)
    const std::vector<FASTAFile::FASTAEntry>& reference() const
    {
      return reference_;
    }

    class OPENMS_DLLAPI Worker
    {
    public:
      explicit Worker(const ProteomeSearch& search);

      /// every PrSM FLASHTnT reports for @p ms2 (the engine's deconvolution of an MS2 at @p tol_ppm, whose raw
      /// spectrum is @p raw and scan number @p scan), in FLASHTnT's order; empty when @p ms2 has no mass. When
      /// @p tag_count is given it receives the number of tags the search produced (FLASHTnTAlgorithm::getTags, read
      /// before the results are released; 0 for a spectrum with no mass).
      std::vector<ProteinHit> search(const DeconvolvedSpectrum& ms2, double tol_ppm, const MSSpectrum& raw, int scan,
                                     int* tag_count = nullptr);

    private:
      const ProteomeSearch& search_;
      FLASHTnTAlgorithm tnt_;
    };

  private:
    std::vector<FASTAFile::FASTAEntry> reference_;
    Param params_;
    Size decoys_dropped_ = 0;
    std::shared_ptr<const FLASHTnTAlgorithm::PreparedDatabase> prepared_;
  };
} // namespace OpenMS
