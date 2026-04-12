// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Deconvolution.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS
{

  /**
   * @brief Manages exploration CE sweep groups, variant scoring, winner selection,
   * and next-level triggering for FLASHIda.
   *
   * Key design principle: Exploration returns commands, never pushes them.
   * The caller (orchestrator) enqueues.
   */
  class OPENMS_DLLAPI Exploration
  {
  public:
    /// Single variant in an exploration CE sweep
    struct ExplorationVariant
    {
      int variant_index = -1;
      double collision_energy = 0.0;
      std::string activation_type;
      std::string tracking_id;
      double score = -1.0;
      float tic_coverage = 0.0f;
      int fragment_count = 0;
      bool received = false;
      DeconvolvedSpectrum result{0};
    };

    /// Group of CE variants for one precursor at one MSn level
    struct ExplorationGroup
    {
      int group_id = 0;
      int msn_level = 2;
      ExplorationMetric exploration_metric = ExplorationMetric::MassCount;
      std::string parent_tracking_id;
      double precursor_mz = 0.0;
      double precursor_mass = 0.0;
      int precursor_charge = 0;
      PeakGroup precursor_pg;
      double isolation_width = 0.0;
      double faims_cv = 0.0;
      uint64_t start_ms = 0;
      std::vector<ExplorationVariant> variants;
      int winner_index = -1;
      bool complete = false;
    };

    /// Lookup reference from tracking ID to exploration group + variant
    struct VariantRef
    {
      int group_id;
      int variant_index;
    };

    /// Construct with a reference to the shared Config and Deconvolution engine
    explicit Exploration(const Config& config, Deconvolution& deconv);

    /// Create exploration group with CE variants. Returns commands for the caller to push.
    std::vector<ScanCommand> initiate(int msn_level, const PeakGroup& pg, int charge,
                                      double faims_cv, ScanCommandQueue& queue);

    /// Process returning exploration variant: deconvolve with correct precursor context,
    /// score, select winner, trigger next level. Returns commands for the caller to push.
    std::vector<ScanCommand> feedResult(int tracking_id,
                                        const double* mzs, const double* ints, int length,
                                        double rt, ScanCommandQueue& queue);

    /// Test-only: feed a pre-deconvolved result (bypasses deconvolution)
    std::vector<ScanCommand> feedResultForTest(int tracking_id,
                                               const DeconvolvedSpectrum& ms2_deconv,
                                               double rt, ScanCommandQueue& queue);

    /// Check whether a tracking_id belongs to an exploration variant
    bool isExplorationVariant(int tracking_id) const;

    /// Generic MSn+1 follow-up after MSn processing. Returns commands for the caller to push.
    std::vector<ScanCommand> initiateNextLevel(int msn_level, const DeconvolvedSpectrum& result,
                                               double faims_cv, ScanCommandQueue& queue);

    /// Number of currently active exploration groups
    int activeGroupCount() const;

    /// Get exploration group by ID (caller must ensure group exists)
    ExplorationGroup getGroup(int group_id) const;

  private:
    const Config& config_;
    Deconvolution& deconv_;

    /// Active exploration groups (group_id -> ExplorationGroup)
    std::unordered_map<int, ExplorationGroup> active_groups_;

    /// Maps tracking_id (int) -> {group_id, variant_index} for variant result routing
    std::unordered_map<int, VariantRef> variant_tracking_map_;

    /// Next group ID (monotonically increasing)
    int next_group_id_ = 1;

    /// Shared implementation: process a deconvolved spectrum for a tracked variant
    std::vector<ScanCommand> feedResultImpl_(int tracking_id, const DeconvolvedSpectrum& ms2_deconv,
                                             double rt, ScanCommandQueue& queue);

    /// Generate CE variant values from min/max/step
    std::vector<double> buildCEVariants_(double ce_min, double ce_max, double ce_step) const;

    /// Dispatch to metric-specific scorer
    double computeExplorationScore_(ExplorationMetric metric, const DeconvolvedSpectrum& spec) const;

    /// Score: number of deconvolved masses
    double computeMassCount_(const DeconvolvedSpectrum& spec) const;

    /// Score: fragmentation efficiency (higher = less remaining precursor)
    double computeRemainingPrecursorScore_(const DeconvolvedSpectrum& spec) const;

    /// Score: number of fragment ions
    double computeFragmentCount_(const DeconvolvedSpectrum& spec) const;

    /// Compute TIC coverage for metadata
    float computeTICCoverage_(const DeconvolvedSpectrum& spec) const;
  };

} // namespace OpenMS
