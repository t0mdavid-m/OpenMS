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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ProteoformTracker.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include <memory>
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
      double reaction_time = 0.0;
      std::string activation_type;
      std::string tracking_id;
      bool is_baseline = false;           ///< CE=0 reference scan for RemainingPrecursor
      double score = -1.0;
      float tic_coverage = 0.0f;
      int fragment_count = 0;
      bool received = false;
      DeconvolvedSpectrum result{0};
      ScanCommand cmd;
      FragmentAnalysis::ProteoformMatch identification_result;  ///< Per-fragment match details (populated at batch eval)
    };

    /// Cached MS2 context for identification.tsv output from MS3 results
    struct MS2Context
    {
      std::string proteoform_sequence;
      int start_pos = 0;
      int end_pos = 0;
      std::vector<FragmentAnalysis::PTMSite> ptm_sites;
      double ms1_precursor_mass = 0.0;
      double ms1_precursor_mz = 0.0;
      int ms1_precursor_charge = 0;
      char fragment_ion_type = '\0';
      int fragment_ion_index = 0;
      double fragment_mass = 0.0;
      double fragment_mz = 0.0;
      int fragment_charge = 0;
      int tag_count = 0;        ///< Sequence-tag count from the parent MS2 (for MS3 results-row tag_count)
      // I2: isolation-window reporting for identification.tsv. MS2 triplet = the parent MS2 precursor;
      // MS3 triplet = the MS3 fragment precursor (0 on MS2 rows). width/charge_intensity are the engine's
      // commanded values copied from the ScanCommand; window_snr is computed over the ACTUAL commanded
      // window (signal/noise incl. co-isolation, FragmentAnalysis::windowSnr), carried via the queue map.
      double ms2_isolation_width = 0.0;
      double ms2_window_snr = 0.0;
      double ms2_charge_intensity = 0.0;
      double ms3_isolation_width = 0.0;
      double ms3_window_snr = 0.0;
      double ms3_charge_intensity = 0.0;
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
      double baseline_intensity = 0.0;    ///< isolation-window intensity from CE=0 scan
      bool has_baseline = false;           ///< whether baseline result has arrived
      bool baseline_failed = false;        ///< RemainingPrecursor: CE=0 baseline had no in-window signal -> abort (no winner)
      uint64_t start_ms = 0;
      std::vector<ExplorationVariant> variants;
      int winner_index = -1;
      bool complete = false;
      ScanCommand originating_cmd{};  ///< MS2 context for buildMS3 (stage 0)
      char fragment_ion_type = '\0';   ///< Fragment ion type (e.g. 'b', 'y') for MS3 scan description
      int fragment_ion_index = 0;      ///< Fragment ion index for MS3 scan description
      MS3FragmentMatcher::ProteoformContext proteoform_ctx; ///< Cached MS2 proteoform for MS3 scoring
    };

    /// Lookup reference from tracking ID to exploration group + variant
    struct VariantRef
    {
      int group_id;
      int variant_index;
    };

    /// Result of initiateNextLevel: commands plus fragment matching metadata
    struct NextLevelResult
    {
      std::vector<ScanCommand> commands;
      std::vector<MS2Context> ms3_contexts;  ///< Parallel to commands: one MS2Context per MS3 command
      std::string matched_protein;
      std::string proteoform_sequence;
      float tic_coverage = 0.0f;
      int fragment_count = 0;
      FragmentAnalysis::ProteoformMatch proteoform_match;  ///< Full fragment detail for identification
    };

    /// Result of feedResult: commands plus exploration per-variant metadata

    // @ Claude we want to clean up this struct. It feels extremely messy.
    struct FeedResultInfo
    {
      struct IdentificationRowInfo
      {
        std::string tracking_id;
        FragmentAnalysis::ProteoformMatch identification_result;
        MS2Context ms2_context;
      };

      std::vector<ScanCommand> commands;
      int group_id = -1;
      int variant_index = -1;
      int total_variants = 0;
      double collision_energy = 0.0;
      std::string activation_type;
      double reaction_time = 0.0;
      double stage0_collision_energy = 0.0;  ///< MS2 isolation stage CE; 0.0 for non-MS3 groups
      std::string stage0_activation_type;    ///< MS2 isolation stage activation; "" for non-MS3 groups
      double stage0_reaction_time = 0.0;     ///< MS2 isolation stage reaction_time; 0.0 for non-MS3 groups
      double score = -1.0;
      float tic_coverage = 0.0f;
      int fragment_count = 0;
      int exploration_metric = 0;
      std::string matched_protein;
      std::string proteoform_sequence;
      double remaining_ratio = -1.0;  ///< Raw remaining_intensity / baseline_intensity (-1.0 = N/A)
      /// F5: encoded tracking id of the winning variant; set ONLY on the group-completing feedResult,
      /// empty ("") on every other (per-variant / non-exploration / no-winner) row. Logged as the trailing
      /// scan_results column so each variant row keeps its OWN metrics while the winner stays identifiable.
      std::string winner_tracking_id;
      char parent_scan_id[4]{};  ///< Parent's encoded tracking ID (from group's originating_cmd)
      FragmentAnalysis::ProteoformMatch identification_result;  ///< Per-fragment match details for identification.tsv
      MS2Context ms2_context;  ///< Cached MS2 context for this variant's group
      std::vector<IdentificationRowInfo> additional_identification_rows;  ///< Extra calibrated rows for other variants in same completed exploration group
      std::vector<std::string> child_ids;  ///< Encoded tracking ids of the pushed children, in info.commands order (pre-encoded by feedResultImpl_; replaces processScan's inline encode loop)
      /// (scan_id -> MS2Context) for production MS3 commands in `commands` that return on the REGULAR MS3
      /// path and need ms2_context_cache_ seeded by the caller. Mirrors the regular MS2->MS3 caching.
      std::vector<std::pair<int, MS2Context>> ms3_context_cache;
    };

    /// Construct with a reference to the shared Config and FragmentAnalysis
    explicit Exploration(const Config& config, FragmentAnalysis& fragments);

    /// Create exploration group with CE variants. Returns commands for the caller to push.
    /// The FAIMS CV is sourced from ms_ctx->faims_cv (Item 1: CV travels via the context, not a param).
    /// When source_spectrum != nullptr, each command's window_snr is computed over it (the survey MS1 for
    /// MS2 variants, the MS2 result for MS3 variants) and stamped on the command (Item 2: SNR travels with
    /// the command). nullptr leaves window_snr at its -1.0 default (used by direct-initiate unit tests).
    /// @param ms_ctx  Optional originating ScanCommand (carries faims_cv; needed for MS3 buildMS3 stage 0)
    std::vector<ScanCommand> initiate(int msn_level, const PeakGroup& pg, int charge,
                                      ScanCommandQueue& queue, const MSSpectrum* source_spectrum = nullptr,
                                      const ScanCommand* ms_ctx = nullptr,
                                      char ion_type = '\0', int frag_index = 0,
                                      const MS3FragmentMatcher::ProteoformContext& proto_ctx = {},
                                      // F2: stage-1 scores of the selected fragment, forwarded to buildMS3
                                      // so MS3 exploration variants carry real *_s1 scalars (default = empty).
                                      const FragmentAnalysis::FragmentScores& frag_scores = {},
                                      // 9b: per-ion best-MS2 params for MS3 stage[0] (ADR-0003). Defaulted so the
                                      // MS2-exploration caller and direct-initiate tests compile unchanged. Only
                                      // consulted on the MS3 (msn_level>=3) buildMS3 branch.
                                      const Ms2Params* stage0_params = nullptr);

    /// Process returning exploration variant: deconvolve with correct precursor context,
    /// score, select winner, trigger next level. Returns FeedResultInfo with commands and metadata.
    /// @param precursor_id The returning variant's per-MS1-selection identity (engine-side map lookup);
    ///        used as the ProteoformTracker model key so each model holds one MS1 selection (one charge).
    FeedResultInfo feedResult(int tracking_id,
                              const double* mzs, const double* ints, int length,
                              double rt, ScanCommandQueue& queue, ProteoformTracker* tracker = nullptr,
                              int precursor_id = 0);

    /// Test-only access (the feedResultImpl_ deconvolution-bypass + getGroup) lives in
    /// FLASHIda_TestAccess.h via this friend, so test scaffolding stays out of the production API.
    friend struct ExplorationTestAccess;

    /// Check whether a tracking_id belongs to an exploration variant
    bool isExplorationVariant(int tracking_id) const;

    /// Generic MSn+1 follow-up after MSn processing. Returns commands plus fragment matching metadata.
    /// @param ms_ctx  Optional originating MS2 ScanCommand (needed for buildMS3 when next_level >= 3)
    /// @param tracker Optional ProteoformTracker: when non-null and next_level >= 3 the model is the
    ///        dispatch authority — it selects the MS3 targets (planNextScans) and Exploration builds
    ///        them (single MS3, or the CE sweep when MS3 exploration is configured). nullptr keeps the
    ///        legacy getTopFragmentMatches direct path (used by tests).
    /// @param precursor_id The originating MS2's per-MS1-selection identity; used as the
    ///        ProteoformTracker model key for the regular MS2->MS3 one-shot model + planNextScans.
    NextLevelResult initiateNextLevel(int msn_level, const DeconvolvedSpectrum& result,
                                      double faims_cv, ScanCommandQueue& queue,
                                      const ScanCommand* ms_ctx = nullptr,
                                      ProteoformTracker* tracker = nullptr,
                                      int precursor_id = 0);

    /// Number of currently active exploration groups
    int activeGroupCount() const;

    /// Mass count of the stored exploration-deconv MS2 spectrum (0 if none). Encapsulates the
    /// exploration_deconv_ access that processScan used to perform inline.
    int explorationDeconvMassCount() const
    {
      return (exploration_deconv_ != nullptr && exploration_deconv_->hasStoredMS2())
                 ? static_cast<int>(exploration_deconv_->storedMS2().size()) : 0;
    }

    /// Pointer to the stored exploration-deconv MS2 spectrum, or nullptr if none. Targets an
    /// engine-owned member (valid until the next exploration deconvolution).
    const DeconvolvedSpectrum* explorationDeconvSpectrum() const
    {
      return (exploration_deconv_ != nullptr && exploration_deconv_->hasStoredMS2())
                 ? &exploration_deconv_->storedMS2() : nullptr;
    }

  private:
    /// Internal deconvolution engine for exploration-variant MS2 spectra. Was public; now reached
    /// only via the accessors above plus Exploration's own internals.
    std::unique_ptr<Deconvolution> exploration_deconv_;

    /// Get exploration group by ID (caller must ensure group exists). Test-only; reached via
    /// ExplorationTestAccess (was public).
    ExplorationGroup getGroup(int group_id) const;

    const Config& config_;
    FragmentAnalysis& fragments_;

    /// Active exploration groups (group_id -> ExplorationGroup)
    std::unordered_map<int, ExplorationGroup> active_groups_;

    /// Maps tracking_id (int) -> {group_id, variant_index} for variant result routing
    std::unordered_map<int, VariantRef> variant_tracking_map_;

    /// Next group ID (monotonically increasing)
    int next_group_id_ = 1;

    /// Shared implementation: process a deconvolved spectrum for a tracked variant.
    /// @p tracker may be nullptr (test-only bypass via ExplorationTestAccess); production always passes a
    /// real tracker.
    FeedResultInfo feedResultImpl_(int tracking_id, const DeconvolvedSpectrum& msn_deconv,
                                   const double* mzs, const double* ints, int length,
                                   ScanCommandQueue& queue, ProteoformTracker* tracker = nullptr,
                                   int precursor_id = 0);

    /// Parameters for one variant in a multi-activation sweep
    struct VariantParams
    {
      std::string activation;
      double collision_energy = 0.0;
      double reaction_time = 0.0;
    };

    /// Generate variants across activation types with CE and/or RT sweep
    std::vector<VariantParams> buildVariants_(const MSLevelConfig& cfg, const ScanConfig& base_config) const;

    /// Dispatch to metric-specific scorer
    double computeExplorationScore_(ExplorationMetric metric, const DeconvolvedSpectrum& spec,
                                    const ExplorationGroup& group,
                                    const double* mzs, const double* ints, int length,
                                    double* out_remaining_ratio,
                                    FragmentAnalysis::ProteoformMatch* out_frag,
                                    const std::string& activation_type) const;

    /// Score: number of deconvolved masses
    double computeMassCount_(const DeconvolvedSpectrum& spec) const;

    /// Score: fragmentation efficiency (higher = less remaining precursor)
    double computeRemainingPrecursorScore_(const ExplorationGroup& group,
                                           const double* mzs, const double* ints, int length,
                                           double* out_ratio = nullptr) const;

    /// Score: number of matched fragment ions + protein info
    FragmentAnalysis::ProteoformMatch computeFragmentMatch_(const DeconvolvedSpectrum& spec, int msn_level,
                                                           const std::string& activation_type) const;

    /// Compute TIC coverage for metadata
    float computeTICCoverage_(const DeconvolvedSpectrum& spec) const;
  };

} // namespace OpenMS
