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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Ms2Params.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/config.h>

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace OpenMS
{

  /// A single fragment observation, either from MS2 or MS3
  struct FragmentObservation
  {
    uint8_t ms_level = 2;       ///< 2 | 3
    double observed_mass = 0;   ///< MS2 frame (MS3 already converted via equiv/adjusted)
    double intensity = 0;
    int source_scan_id = 0;
    Ms2Params params;           ///< For an MS3 obs: the parent MS2's params
    double frag_mz = 0;         ///< Isolation m/z for targeting this fragment in MS3
    int frag_charge = 0;        ///< Charge state for targeting this fragment in MS3
  };

  /// Key identifying a fragment ion: (ion_type, ion_index)
  using FragmentKey = std::pair<std::string, int>;

  /// Hash for FragmentKey so it can be used as an unordered_map key
  struct FragmentKeyHash
  {
    std::size_t operator()(const FragmentKey& k) const noexcept
    {
      std::size_t h1 = std::hash<std::string>{}(k.first);
      std::size_t h2 = std::hash<int>{}(k.second);
      return h1 ^ (h2 * 2654435761ULL + 0x9e3779b9ULL + (h1 << 6) + (h1 >> 2));
    }
  };

  /// Accumulated observations and coverage for a single fragment ion across all scans
  struct MappedFragment
  {
    std::string ion_type;
    int ion_index = 0;
    bool is_prefix = false;
    int cover_start = 0, cover_end = 0;    ///< 1-based inclusive, proteoform frame
    double theoretical_mass = 0;
    std::optional<FragmentObservation> best_ms2, best_ms3;
    int n_ms2 = 0, n_ms3 = 0;
  };

  /// A detected or candidate modification with localization support bounds
  struct ModificationState
  {
    double mass_shift = 0;
    int candidate_start = 0, candidate_end = 0;   ///< start==end => localized
    double support_lower = 0, support_upper = 0;
  };

  /// Objective driving next-scan planning for a proteoform
  enum class CharacterizationObjective
  {
    Ambiguity,  ///< Resolve PTM site ambiguity
    Coverage    ///< Extend sequence coverage
  };

  /// Richer per-peak record so mapScanOntoModel_ can recover mz and charge for MS3 targeting
  struct PeakRecord
  {
    double mono_mass = 0;
    double mz = 0;       ///< Isolation centre m/z: (mz1+mz2)/2 from PeakGroup::getMzRange(charge)
    int charge = 0;      ///< getMaxIntensityAbsCharge() — same charge used for mz derivation
    double intensity = 0;
  };

  /// A deconvolved scan awaiting integration into a ProteoformModel
  struct PendingScan
  {
    int scan_id = 0;
    uint8_t ms_level = 2;
    Ms2Params params;
    std::vector<PeakRecord> peaks;   ///< One entry per deconvolved PeakGroup
    FragmentAnalysis::ProteoformMatch match;
    double id_score = -1;
  };

  /// Per-precursor accumulated state: fragments, modifications, pending scans
  struct ProteoformModel
  {
    int nominal_mass = 0;
    bool finalized = false;
    std::string proteoform_sequence;
    int region_start = -1, region_end = -1;
    double identification_score = -1;
    int winner_scan_id = 0;
    int update_index = 0;
    std::vector<PendingScan> pending;
    std::vector<ModificationState> modifications;
    std::unordered_map<FragmentKey, MappedFragment, FragmentKeyHash> fragments;
    ScanCommand ms2_ctx;        ///< MS2 command context for buildMS3 (captured once from the first feedScan)
    bool has_ms2_ctx = false;   ///< True once ms2_ctx has been set

    /// Fraction of proteoform residues covered by at least one fragment observation [0,1]
    double coveragePct() const;

    /// Union of all MS2-frame observed masses across MS2 and MS3 observations
    std::vector<double> combinedMs2FrameMasses() const;
  };

  /**
   * @brief Per-precursor proteoform tracking model for the FLASHIda characterization feature.
   *
   * Accumulates fragment observations from MS2 and MS3 scans belonging to the same
   * nominal-mass precursor, maps them onto a ProteoformModel, narrows modification
   * localization, and plans the next acquisition scans to maximize coverage or
   * resolve ambiguity.
   *
   * Key design principle: ProteoformTracker returns commands, never pushes them.
   * The caller (orchestrator) enqueues.
   */
  class OPENMS_DLLAPI ProteoformTracker
  {
  public:
    /// Construct with references to the shared Config and IdaLogger
    ProteoformTracker(const Config& cfg, IdaLogger& logger);

    /**
     * @brief Integrate a new deconvolved scan into the model for @p nominal_mass.
     *
     * Creates the model entry on first call. Appends a PendingScan to pending,
     * then triggers mapScanOntoModel_ + narrowModifications_ (skeleton: no-ops).
     */
    void feedScan(int nominal_mass, uint8_t ms_level, const Ms2Params& params, int scan_id,
                  const DeconvolvedSpectrum& deconv,
                  const FragmentAnalysis::ProteoformMatch& match, double id_score,
                  const ScanCommand& ms2_ctx);

    /// Mark the model for @p nominal_mass as finalized (no further scans expected).
    void finalize(int nominal_mass);

    /**
     * @brief Plan the next acquisition scans for @p nominal_mass.
     *
     * Inspects the current ProteoformModel and returns a list of ScanCommands
     * for the caller to push onto @p queue. Returns empty if no scans are needed
     * or the model is finalized.
     */
    std::vector<ScanCommand> planNextScans(int nominal_mass, ScanCommandQueue& queue);

    /// Return a pointer to the model for @p nominal_mass, or nullptr if absent.
    const ProteoformModel* model(int nominal_mass) const;

  private:
    /// Map all observations in @p scan onto @p mdl's fragments map (skeleton: no-op).
    void mapScanOntoModel_(ProteoformModel& mdl, const PendingScan& scan);

    /// Update ModificationState entries in @p mdl based on current fragment coverage (skeleton: no-op).
    void narrowModifications_(ProteoformModel& mdl);

    /// Emit a log/TSV row for the current state of @p mdl (skeleton: no-op).
    void emitRow_(const ProteoformModel& mdl);

    const Config& config_;
    IdaLogger& logger_;

    /// Active models keyed by nominal precursor mass
    std::unordered_map<int, ProteoformModel> models_;
  };

} // namespace OpenMS
