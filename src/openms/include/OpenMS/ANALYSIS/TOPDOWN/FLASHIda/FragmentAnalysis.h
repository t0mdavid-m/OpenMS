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
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{

  /**
   * @brief Fragment analysis component for FLASHIda.
   *
   * Owns all MS2 fragment matching methods: direct mass queries, tag-based
   * fragment matching (FLASHTagger + FLASHExtender), terminal fragment ion
   * selection, and PTM ambiguity-enclosing ion selection.
   *
   * All query methods take the stored MS2 DeconvolvedSpectrum as a parameter
   * so that FragmentAnalysis is stateless with respect to deconvolution results.
   */
  class OPENMS_DLLAPI FragmentAnalysis
  {
  public:
    /// Structure representing a PTM site detected by FLASHExtender
    struct PTMSite
    {
      int position;        ///< Position in protein sequence (1-based, midpoint)
      int start_position;  ///< Start of the region where PTM could be localized (1-based)
      int end_position;    ///< End of the region where PTM could be localized (1-based)
      double mass_shift;   ///< Observed mass shift (modification mass)
    };

    /// Cached result from the last runTagBasedFragmentMatching_ call
    struct ProteoformInfo
    {
      int region_start = -1;  ///< 0-based start position in protein sequence (-1 = full sequence)
      int region_end = -1;    ///< 0-based exclusive end position (-1 = full sequence)
      std::vector<PTMSite> ptm_sites; ///< 1-based positions relative to proteoform
    };

    /// Get the proteoform region and PTM sites from the last tag-based matching call
    const ProteoformInfo& getLastProteoformInfo() const { return last_proteoform_info_; }

    /// Map fragmentation method name to ion type strings for FLASHTagger/FLASHExtender.
    /// Case-insensitive. Returns {"b","y"} for HCD/CID, {"c","z"} for ETD,
    /// {"b","c","y","z"} for EThcD/EtCID, {"a","b","c","x","y","z"} for UVPD.
    /// Defaults to {"b","y"} for unknown methods.
    static std::vector<std::string> getIonTypesForFragmentationMethod(const String& method);

    /// Constructor: takes a reference to the shared Config
    explicit FragmentAnalysis(const Config& config);

    // ---------------------------------------------------------------
    // Fragment query methods -- all take stored MS2 as parameter
    // ---------------------------------------------------------------

    /**
     * @brief Get the top N MS2 masses with isolation window info
     *
     * Requires prior MS2 deconvolution. Returns masses sorted by intensity.
     *
     * @param n maximum number of masses to return
     * @param masses output: monoisotopic masses
     * @param qscores output: quality scores / intensities
     * @param charges output: representative charges
     * @param window_starts output: isolation window start m/z
     * @param window_ends output: isolation window end m/z
     * @param stored_ms2 deconvolved MS2 spectrum
     * @return actual number of masses returned
     */
    int getBestMS2Masses(int n,
                         double* masses,
                         double* qscores,
                         int* charges,
                         double* window_starts,
                         double* window_ends,
                         DeconvolvedSpectrum& stored_ms2);

    /**
     * @brief Get top fragment ion matches against a protein sequence, sorted by qscore
     *
     * Performs tag-based fragment matching using stored MS2 deconvolution results.
     * Uses MS2 tolerance from config (config_.level(2).tolerance_ppm).
     *
     * @param protein_sequence the protein sequence to match against
     * @param n maximum number of matches to return
     * @param masses output: observed monoisotopic masses
     * @param qscores output: qscores of matched peaks
     * @param charges output: charges of matched peaks
     * @param window_starts output: isolation window start m/z values
     * @param window_ends output: isolation window end m/z values
     * @param ion_types output: ion type characters ('a','b','c','x','y','z')
     * @param fragment_indices output: 1-based fragment index
     * @param stored_ms2 deconvolved MS2 spectrum
     * @param fragmentation_method fragmentation method name (default "HCD")
     * @return number of matches found (may be less than n)
     */
    int getTopFragmentMatches(const String& protein_sequence,
                              int n,
                              double* masses,
                              double* qscores,
                              int* charges,
                              double* window_starts,
                              double* window_ends,
                              char* ion_types,
                              int* fragment_indices,
                              DeconvolvedSpectrum& stored_ms2,
                              const String& fragmentation_method = "HCD",
                              double tolerance_ppm = 0.0);

    /**
     * @brief Get terminal (innermost) fragment ions sorted by sequence position
     *
     * Returns fragment ions that extend furthest toward the opposite terminus:
     * - b-ions: sorted by highest fragment_index (rightmost, closest to C-terminus)
     * - y-ions: sorted by highest fragment_index (leftmost, closest to N-terminus)
     *
     * Output is interleaved: [top_b, top_y, 2nd_b, 2nd_y, ...]
     *
     * @param protein_sequence the protein sequence to match against
     * @param n maximum number of results to return
     * @param masses output: observed monoisotopic masses
     * @param qscores output: quality scores
     * @param charges output: charge states
     * @param window_starts output: isolation window start m/z
     * @param window_ends output: isolation window end m/z
     * @param ion_types output: ion type characters
     * @param fragment_indices output: 1-based fragment index
     * @param stored_ms2 deconvolved MS2 spectrum
     * @param fragmentation_method fragmentation method name (default "HCD")
     * @return number of results returned (may be < n)
     */
    int getTerminalFragmentIons(const String& protein_sequence,
                                int n,
                                double* masses,
                                double* qscores,
                                int* charges,
                                double* window_starts,
                                double* window_ends,
                                char* ion_types,
                                int* fragment_indices,
                                DeconvolvedSpectrum& stored_ms2,
                                const String& fragmentation_method = "HCD",
                                double tolerance_ppm = 0.0);

    /**
     * @brief Get unique fragment ions that enclose PTM ambiguity regions
     *
     * Identifies PTM sites and returns the best fragment ions that bracket
     * ambiguity regions. Returns deduplicated ions sorted by qscore.
     * Useful for MS3 targeting to resolve PTM localization.
     *
     * @param protein_sequence the protein sequence to analyze
     * @param n maximum number of ions to return
     * @param masses output: monoisotopic masses
     * @param qscores output: quality scores
     * @param charges output: representative charges
     * @param window_starts output: isolation window start m/z
     * @param window_ends output: isolation window end m/z
     * @param ion_types output: ion type characters
     * @param fragment_indices output: 1-based fragment index
     * @param stored_ms2 deconvolved MS2 spectrum
     * @param fragmentation_method fragmentation method name (default "HCD")
     * @return number of enclosing ions found
     */
    int getAmbiguityEnclosingIons(const String& protein_sequence,
                                  int n,
                                  double* masses,
                                  double* qscores,
                                  int* charges,
                                  double* window_starts,
                                  double* window_ends,
                                  char* ion_types,
                                  int* fragment_indices,
                                  DeconvolvedSpectrum& stored_ms2,
                                  const String& fragmentation_method = "HCD",
                                  double tolerance_ppm = 0.0);

    // ---------------------------------------------------------------
    // Python-friendly overloads
    // ---------------------------------------------------------------

    int getBestMS2MassesPy(int n,
                           std::vector<double>& masses,
                           std::vector<double>& qscores,
                           std::vector<int>& charges,
                           std::vector<double>& window_starts,
                           std::vector<double>& window_ends,
                           DeconvolvedSpectrum& stored_ms2);

    int getTopFragmentMatchesPy(const String& protein_sequence,
                                int n,
                                std::vector<double>& masses,
                                std::vector<double>& qscores,
                                std::vector<int>& charges,
                                std::vector<double>& window_starts,
                                std::vector<double>& window_ends,
                                std::vector<int>& is_b_ions,
                                std::vector<int>& fragment_indices,
                                DeconvolvedSpectrum& stored_ms2);

    int getTerminalFragmentIonsPy(const String& protein_sequence,
                                  int n,
                                  std::vector<double>& masses,
                                  std::vector<double>& qscores,
                                  std::vector<int>& charges,
                                  std::vector<double>& window_starts,
                                  std::vector<double>& window_ends,
                                  std::vector<int>& is_b_ions,
                                  std::vector<int>& fragment_indices,
                                  DeconvolvedSpectrum& stored_ms2);

    int getAmbiguityEnclosingIonsPy(const String& protein_sequence,
                                    int n,
                                    std::vector<double>& masses,
                                    std::vector<double>& qscores,
                                    std::vector<int>& charges,
                                    std::vector<double>& window_starts,
                                    std::vector<double>& window_ends,
                                    std::vector<int>& is_b_ions,
                                    std::vector<int>& fragment_indices,
                                    DeconvolvedSpectrum& stored_ms2);

  private:
    const Config& config_;
    ProteoformInfo last_proteoform_info_;

    /// Internal struct for fragment match results from tag-based matching
    struct TagBasedFragmentMatch
    {
      int peak_index;           ///< Index in stored MS2 spectrum
      double observed_mass;     ///< Monoisotopic mass
      double theoretical_mass;  ///< Theoretical mass
      double qscore;            ///< Quality score from PeakGroup
      int charge;               ///< Charge state
      int fragment_index;       ///< 1-based position in protein sequence
      char ion_type;            ///< Ion type: 'a', 'b', 'c', 'x', 'y', or 'z'
      double ppm_error;         ///< ppm error from match
    };

    /**
     * @brief Run FLASHTagger+FLASHExtender workflow to get matched fragments
     *
     * Uses the stored MS2 deconvolution to:
     * 1. Generate sequence tags via FLASHTagger
     * 2. Validate tags match the protein sequence
     * 3. Run FLASHExtender for path-based matching
     * 4. Extract matched fragments with qscores, charges, and positions
     *
     * @param protein_sequence the protein sequence to match against
     * @param matches output: vector of matched fragments sorted by qscore descending
     * @param stored_ms2 deconvolved MS2 spectrum
     * @param ptm_sites optional output: PTM sites detected by FLASHExtender
     * @param fragmentation_method fragmentation method name (default "HCD")
     * @return number of matches found
     */
    int runTagBasedFragmentMatching_(const String& protein_sequence,
                                    std::vector<TagBasedFragmentMatch>& matches,
                                    DeconvolvedSpectrum& stored_ms2,
                                    std::vector<PTMSite>* ptm_sites = nullptr,
                                    const String& fragmentation_method = "HCD",
                                    double tolerance_ppm = 0.0);
  };

} // namespace OpenMS
