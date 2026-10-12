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

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/RefineSettings.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

/**
  @brief Refine -- the annotation method beside FLASHTnT (ADR-0049): a mass shift is detected before the path.

  Pure functions, no engine state. DETECT scans the mass-shift axis dm on a fixed grid and, at every
  grid point, counts the ladder ions of the supplied sequence that a deconvolved fragment mass explains
  once shifted by dm (the curve). Observed counts against all masses, Expected against the masses the
  unshifted ladder already explains, and the Signal S = Observed - Expected. The null is K shuffled
  sequences (the residues permuted: composition, length and fragment count kept, only the ladder goes)
  counted against the unexplained masses. FDR(h) = mean decoy peaks with >= h ions / Signal peaks with
  >= h ions; the cut is the smallest h >= min_ions with FDR(h) <= fdr, and the LEVELS are the regions of
  S at or above the cut. Three exclusions mark a level rather than drop it (ADR-0049 as amended 2026-10-12): a neutral
  loss (water/ammonia below 0 or below a stronger level), then a neutron satellite (1-3 x 1.003355 Da beside 0 or a
  stronger level), then the cap (of the rest, only the max_levels strongest are kept). ASSEMBLE builds its paths from
  the kept levels and charges every jump jump_penalty ions.

  The reference implementation is the owner's prototype (trunc_core.py / detect_fdr.py of the Refine
  publication notes); every rule here is that prototype's, one to one.
*/
namespace OpenMS
{
  // RefineSettings lives in RefineSettings.h (included above) so Config.h can hold it by value
  // without the Config.h -> Refine.h -> FragmentAnalysis.h -> Config.h include cycle.

  class OPENMS_DLLAPI Refine
  {
  public:
    /// Theoretical prefix/suffix ladders of one sequence for an ion alphabet (I->L, monoisotopic).
    /// prefix[n-1] = ion of cleavage site n (b_n or c_n), suffix[n-1] = y_{L-n} / z_{L-n}.
    /// residues = the per-residue masses in sequence order; sequence_mass = their sum + water.
    struct Ladder
    {
      std::vector<double> prefix, suffix;
      std::vector<double> residues;
      double sequence_mass = 0;
    };
    /// The prefix series is the first prefix type (a/b/c) of @p ion_types, the suffix series the first
    /// suffix type (x/y/z); b and y when a side names none. Masses come from
    /// FragmentAnalysis::computePTMAdjustedFragmentMasses, so the residue arithmetic is the engine's.
    static Ladder ladder(const String& sequence, const std::vector<std::string>& ion_types);

    /// One detected level: a region of the Signal at or above the cut.
    struct Level
    {
      /// Why a level is not one Assemble builds from; the first rule that matches, in this order.
      enum class Exclusion
      {
        None,             ///< kept
        NeutralLoss,      ///< water/ammonia below 0 or below a stronger level
        NeutronSatellite, ///< 1-3 neutrons (1.003355 Da) above or below 0 or a stronger level
        Capped            ///< beyond the max_levels strongest of the rest
      };
      double shift = 0;         ///< dm at the middle of the region's highest grid points
      int signal = 0;           ///< the region's maximum of S
      int n_prefix = 0;         ///< prefix-ladder ions that match at shift
      int n_suffix = 0;         ///< suffix-ladder ions that match at shift
      double lo = 0, hi = 0;    ///< dm of the region's first and last grid point
      Exclusion excluded = Exclusion::None;
    };
    struct Detect
    {
      int cut = 0;                       ///< smallest h >= min_ions with FDR(h) <= fdr; 0 = none passes
      std::vector<Level> levels;         ///< regions of S >= cut, ascending shift, every exclusion marked (none dropped)
      std::vector<double> fdr_by_height; ///< index h (1..60): mean decoy peaks / signal peaks; -1 where S has none
      int observed_at_zero = 0;          ///< masses the unshifted ladder explains (|M n T|)
      int n_kept = 0;                    ///< levels with excluded == None: what Lambda_b / Lambda_y are built from
    };
    /// Detect on one scan's deconvolved masses (ascending; @p intensities parallel -- the cap's tie-break reads them).
    static Detect detect(const std::vector<double>& masses_sorted, const std::vector<double>& intensities, const Ladder& ladder,
                         const RefineSettings& s);

    // --- building blocks, exposed for the unit sections ---
    /// n(f, h) for every h in [1, hmax]: maximal runs of consecutive grid points with f >= h.
    /// The result has hmax + 1 entries; index 0 is unused (0).
    static std::vector<int> regionCounts(const std::vector<int>& curve, int hmax);
    /// The curve: number of distinct ladder ions matching at least one mass, at every grid point lo + i*step.
    static std::vector<int> curve(const std::vector<double>& theo, const std::vector<double>& masses_sorted,
                                  double tol_ppm, double lo, double step, int n);
    /// Mask of masses the unshifted ladder explains within tolerance.
    static std::vector<bool> explained(const std::vector<double>& theo, const std::vector<double>& masses_sorted, double tol_ppm);
    /// Marks NeutralLoss (never drops) a level whose region holds r - v for a reference r in {0} u {shift of a level with
    /// a strictly higher signal} and a loss v in {water, ammonia}, with 0.005 Da slack on either side.
    static void flagNeutralLosses(std::vector<Level>& levels);
    /// Marks NeutronSatellite (never drops) a level not already excluded whose region holds r + s n 1.003355 for a reference
    /// r in {0} u {shift of any level with a strictly higher signal}, n in {1, 2, 3}, s in {+1, -1}, with 0.005 Da slack.
    /// Run after flagNeutralLosses: a level that is both stays a NeutralLoss.
    static void flagNeutronSatellites(std::vector<Level>& levels);
    /// The cap: if more than @p max_levels levels are not excluded, marks all but the max_levels strongest Capped -- by
    /// signal (desc), then @p tie_intensity (desc; parallel to @p levels: the summed intensity of the masses the level's ions
    /// match at its shift), then the smaller |shift|. Returns the number of levels left with excluded == None.
    static int capLevels(std::vector<Level>& levels, const std::vector<double>& tie_intensity, int max_levels);

    // --- ASSEMBLE: one column DP per lambda_k (column_dp_v85.py), the targets (targets_v88.py), the MS3 recount (ms3_recount_v79.py) ---

    /// The best path of one lambda_k: levels 0 = v_0, ..., v_k = lambda_k; a cleavage site n on v_j explains its prefix ion at
    /// t_n + v_j and its suffix ion at m_seq - t_n + lambda_k - v_j.
    struct Path
    {
      double total = 0;                        ///< lambda_k
      int score = 0;                           ///< S = ions - jump_penalty x jumps
      int ions = 0;                            ///< N, the explained ions (one mass may explain two)
      int jumps = 0;                           ///< mass shifts
      double intensity_sum = 0;                ///< sum of the intensities of the distinct explained masses (a tie-break)
      std::vector<double> levels;              ///< 0, ..., lambda_k in sequence order
      std::vector<std::pair<int, int>> windows; ///< per jump: first and last residue (1-based) it may sit on keeping N
      struct Ion
      {
        char type;          ///< the ladder's prefix or suffix ion type (b/c/a, y/z/x)
        int index;          ///< ion number: n for a prefix ion of site n, L - n for a suffix ion
        int mass_idx;       ///< index of the explaining observed mass
        double level;       ///< the shift on this ion's own series: v_j (prefix), lambda_k - v_j (suffix)
        double theoretical; ///< the unshifted ladder ion plus that shift
      };
      std::vector<Ion> ions_explained; ///< at the placement with every jump at its window's first residue, in site order
    };
    struct Assemble
    {
      std::vector<Path> candidates; ///< one per lambda_k; score desc, jumps asc, intensity_sum desc, then construction order
    };
    /// Assemble on one scan's deconvolved masses (ascending; @p intensities parallel). Lambda_b / Lambda_y = the levels of
    /// @p d that are not excluded with at least one prefix / suffix ion. A jump costs s.jump_penalty ions.
    static Assemble assemble(const std::vector<double>& masses_sorted, const std::vector<double>& intensities, const Ladder& ladder,
                             const Detect& d, const RefineSettings& s);

    /// An observed mass two paths explain differently: as different ions, or one path not at all.
    struct Discriminating
    {
      int mass_idx = -1;
      double mass = 0;
      double intensity = 0;
      std::optional<Path::Ion> ion_a, ion_b;
    };
    /// Masses the two paths explain differently (one or none, or different ions); intensity desc. One ion per mass and
    /// path: a mass both paths give the same ion is shared; otherwise a path's second ion of a mass goes for the one in the
    /// series of the other path's ion, else the prefix ion (targets_v88.py).
    static std::vector<Discriminating> discriminating(const Path& a, const Path& b, const std::vector<double>& masses,
                                                      const std::vector<double>& intensities);

    /// PTM sites a path implies (full-protein 1-based ranges = windows), for MS3FragmentMatcher.
    static std::vector<FragmentAnalysis::PTMSite> ptmSites(const Path& p);

    struct PieceVerdict
    {
      bool fits_a = false, fits_b = false;
      int shared_seen = 0;      ///< observed masses matching a shared piece at tight_ppm after the calibration
      double correction = 1.0;  ///< 1 / (1 + median ppm of the shared pieces at loose_ppm)
    };
    /// One MS3 scan of one discriminating mass: calibrate on shared pieces, then test discriminating,
    /// shift-carrying pieces per path (decision 8).
    static PieceVerdict pieceTest(const DeconvolvedSpectrum& ms3, const String& sequence, const Path& a, const Path& b,
                                  const Discriminating& target, double loose_ppm, double tight_ppm);

    struct Recount
    {
      int score_a = 0, score_b = 0;
      int confirmed_a = 0, confirmed_b = 0;
      bool a_wins = true;
    };
    /// S after MS3 = shared + confirmed targets of the path - jump_penalty x jumps (ms3_recount_v79.py, with the amended
    /// penalty); on equal S fewer jumps.
    /// shared = the DISTINCT observed masses both paths explain alike (Version 88: one ion per mass and path), counted as
    /// a path's explained masses minus the targets it explains -- masses, not ions (a path may explain one mass as two
    /// ions). Equal for both paths by construction; throws Exception::InvalidParameter when it is not.
    static Recount recount(const Path& a, const Path& b, const std::vector<Discriminating>& targets,
                           const std::vector<std::pair<bool, bool>>& confirmed /* per target: fits a, fits b */, int jump_penalty);
  };

  /// What Refine decided on one MS2 scan (FragmentAnalysis::ProteoformMatch::refine, ADR-0049): Detect, every candidate
  /// Assemble built (decision 5), and the precursor mass difference of the best path.
  struct RefineCandidates
  {
    Refine::Detect detect;
    Refine::Assemble assemble;
    double precursor_mass_difference = 0; ///< precursor mono mass - (sequence mass + best total); 0 without a precursor mass
    int cut = 0;                          ///< detect.cut
    int n_levels = 0;                     ///< detect.n_kept (logged as n_levels)
  };
} // namespace OpenMS
