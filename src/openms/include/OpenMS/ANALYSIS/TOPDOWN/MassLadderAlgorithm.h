// Copyright (c) 2002-2026, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace OpenMS
{
class ProteoformDatabase;

/**
@brief Candidate proteins from mass ladders: deconvolved fragment masses read as residue sums of a database protein, without any
sequence tag. Mass ladders are a candidate source of their own, next to sequence tags (FLASHTaggerAlgorithm) and gapped tags
(FLASHGappedTaggerAlgorithm).

- Terminal ladders: masses that are residue sums from an intact protein terminus (the N-terminus also without the initiator Met and
  with acetylation). A protein is scored by the chance probability of the number of masses on its ladder.
- Precursor-anchored ladders: the windows (stretches) of the proteins whose unmodified mass equals the precursor mass, scored by the
  chance probability of the masses on their prefix ladder (read from the window start) and their suffix ladder (read from the window
  end).
- Sliding ladders: masses read from any start of a protein (prefix ion types) or to any end (suffix ion types), without a precursor
  mass; a protein is scored by the chance probability of its best ladder, corrected for the starts or ends tested.

Only the standard residues (ACDEFGHIKLMNPQRSTVWYUO) make up a ladder: a terminal ladder ends at the first other residue, and no window
crosses one. buildIndex builds the database side once, from the protein sequences and the run's ion types; the searches of a spectrum
only read it, so spectra can be searched in parallel. A search returns plain candidates (protein, rank, score, and the matched masses of
the best ion type of each side as a ladder spanning residues of the protein), from which the caller makes its tags and protein hits.
@ingroup Topdown
*/
class OPENMS_DLLAPI MassLadderAlgorithm
{
public:
  /// precursor ladders: a window needs this many matched ladder positions to be a candidate (the prototype's calibrated setting)
  static constexpr int precursor_ladder_min_k = 3;
  /// precursor ladders: a window has at least this many residues. Shorter windows were chance hits: on Eclip_1480 against E. coli + B.
  /// theta (keep-decoy), the PrSMs of windows under 10 residues were 3 decoys, 1 B. theta and 2 E. coli PrSMs; from 10 residues on, 1
  /// decoy, 2 B. theta and 167 E. coli PrSMs (none of the B. theta or decoy ones at 1%). A short window has few ladder positions, its
  /// masses are where fragment masses are densest, and many windows fit a small precursor mass.
  static constexpr int precursor_ladder_min_length = 10;

  /// the matched masses of one ion type of a candidate, read from one end (a protein terminus or a window end)
  struct Ladder
  {
    bool prefix = true;             ///< a prefix ion type (read from the N-terminal end); otherwise a suffix ion type (from the C-terminal end)
    int start = 0;                  ///< first residue it spans (0-based): the residues between its lowest and its highest mass
    int end = 0;                    ///< one past the last residue it spans
    std::vector<Size> mass_indices; ///< the matched masses (indices into the masses searched), from the lowest; at least two
  };

  /// a candidate protein of a spectrum
  struct Candidate
  {
    Size protein = 0;            ///< index of the protein in the sequences given to buildIndex
    int rank = 0;                ///< rank among the best candidates of the spectrum (1-based); a ranked protein without a ladder is left out, its rank unused
    double score = 0;            ///< terminal ladders: -log10 P of the better side, corrected for its number of shifts; precursor ladders:
                                 ///< -log10 P of the window minus log10 of the number of windows that fit the precursor mass; sliding
                                 ///< ladders: -log10 P of the best ladder minus log10 of the anchors tested of its class (termini, others)
    std::vector<Ladder> ladders; ///< the prefix side first, then the suffix side (each side at most once)
    int window_start = -1;       ///< precursor ladders: the window is residues window_start..window_end-1 (0-based); -1 for terminal ladders
    int window_end = -1;
    int isotope = 0;             ///< precursor ladders: the isotope error (-1, 0, 1) at which the window mass equals the precursor mass
  };

  /// the precursor-ladder candidates of one precursor mass
  struct PrecursorLadderResult
  {
    Size windows = 0;                  ///< the windows whose mass equals the precursor mass (hypotheses), whatever their match count
    std::vector<Candidate> candidates; ///< the best windows with at least precursor_ladder_min_k matched positions, best first
  };

  /// precursor ladders: a run of standard residues a..b-1 of a protein, in which the windows lie. A window i..j-1 of it contains every
  /// fixed modification of its entry (ProteoformDatabase): i <= i_max and j >= j_min (no bound without fixed modifications).
  struct Segment
  {
    std::uint32_t protein, a, b;
    std::uint32_t i_max = UINT32_MAX, j_min = 0;
  };

  /// precursor ladders of a database with several entries per protein (ProteoformDatabase): the windows of entries of one protein that
  /// are the same residues with the same fixed modifications are one hypothesis, the window of the entry with the shortest form that
  /// contains it (then the lowest index). Empty: every window is its own.
  struct WindowOwnership
  {
    std::vector<int> offset, length;          ///< per entry: residues of its protein before it, and its residues
    std::vector<std::uint32_t> class_first;   ///< per entry: its class, members[class_first .. class_first + class_size - 1]
    std::vector<std::uint32_t> class_size;
    std::vector<std::uint32_t> members;       ///< per class: its entries by length, then index

    bool empty() const
    {
      return offset.empty();
    }

    /// whether the window i..j-1 of @p entry is that entry's hypothesis
    bool owns(std::uint32_t entry, std::size_t i, std::size_t j) const
    {
      const long s = (long)offset[entry] + (long)i, t = (long)offset[entry] + (long)j;
      for (std::uint32_t k = class_first[entry]; k < class_first[entry] + class_size[entry]; k++)
      {
        const std::uint32_t m = members[k];
        if (offset[m] <= s && t <= (long)offset[m] + length[m]) return m == entry;
      }
      return true;
    }
  };

  /// precursor ladders: a window found by preparePrecursorWindows, residues a + start .. a + start + length - 1 of a segment
  struct PreparedWindow
  {
    std::uint32_t segment;
    std::uint16_t start, length;
  };

  MassLadderAlgorithm() = default;
  /// not copyable: the window residue sums point into the object's own storage
  MassLadderAlgorithm(const MassLadderAlgorithm&) = delete;
  MassLadderAlgorithm& operator=(const MassLadderAlgorithm&) = delete;
  MassLadderAlgorithm(MassLadderAlgorithm&&) = default;
  MassLadderAlgorithm& operator=(MassLadderAlgorithm&&) = default;
  ~MassLadderAlgorithm() = default;

  /**
    @brief build the database side of the searches
    @param protein_seqs the protein sequences (the candidates refer to them by index)
    @param ion_types the run's ion types (a, b, c prefix; x, y, z, zp1, zp2 suffix)
    @param terminal_ladders build the index of searchTerminalLadders
    @param precursor_ladders build the segments of searchPrecursorLadders
    @param db the database the sequences are the entries of (their fixed modifications are added to the residue sums, and the searches
    rank one entry per protein); nullptr or an identity database: every sequence is a protein of its own
  */
  void buildIndex(const std::vector<std::string>& protein_seqs, const std::vector<std::string>& ion_types, bool terminal_ladders,
                  bool precursor_ladders, const ProteoformDatabase* db = nullptr);

  /**
    @brief the terminal-ladder candidates of a spectrum. Thread safe once the index is built.
    @param masses the deconvolved masses of the spectrum, ascending
    @param tol_ppm mass tolerance of a deconvolved mass in ppm
    @param top_n at most this many candidates (the best proteins by score; ties by protein index)
    @param excluded proteins (indices) that are not proposed
    @param candidates output, best first
  */
  void searchTerminalLadders(const std::vector<double>& masses, double tol_ppm, Size top_n, const std::set<Size>& excluded,
                             std::vector<Candidate>& candidates) const;

  /**
    @brief the precursor-ladder candidates of a spectrum for each of @p precursor_masses (e.g. its precursor mass and its decoy precursor
    mass; the masses are prepared once). Thread safe once the index is built.
    @param masses the deconvolved masses of the spectrum, ascending
    @param precursor_masses the precursor masses; one result per mass, in this order
    @param fragment_tol_ppm mass tolerance of a deconvolved mass in ppm
    @param precursor_tol_ppm mass tolerance of a window mass against a precursor mass in ppm
    @param top_n at most this many candidates per precursor mass
    @param results output
  */
  void searchPrecursorLadders(const std::vector<double>& masses, const std::vector<double>& precursor_masses, double fragment_tol_ppm,
                              double precursor_tol_ppm, Size top_n, std::vector<PrecursorLadderResult>& results) const;

  /**
    @brief precursor ladders: the windows that fit each of @p precursor_masses (e.g. those of all spectra of a run), found in one pass over
    the database. searchPrecursorLadders then reads the windows of a prepared mass instead of sweeping the database for it, with the same
    results; other masses (or another precursor tolerance) are swept as before, and so are all with fewer than 256 distinct masses. Not
    thread safe: call it before the searches. Needs the segments (buildIndex with @p precursor_ladders).
    @param precursor_masses the precursor masses (any order; repeats and masses <= 0 are ignored)
    @param precursor_tol_ppm mass tolerance of a window mass against a precursor mass in ppm, as given to searchPrecursorLadders
  */
  void preparePrecursorWindows(const std::vector<double>& precursor_masses, double precursor_tol_ppm);

  /// frees the windows of preparePrecursorWindows (the searches sweep the database again)
  void clearPrecursorWindows();

  /**
    @brief the sliding-ladder candidates of a spectrum: masses read as residue sums from any start of a protein (prefix ion types) or to
    any end (suffix ion types), so that proteoforms truncated at either end are proposed without a sequence tag or a precursor mass.
    Needs the segments (buildIndex with @p precursor_ladders). Thread safe once the index is built.
    @param masses the deconvolved masses of the spectrum, ascending
    @param tol_ppm mass tolerance of a deconvolved mass in ppm
    @param top_n at most this many candidates (the best proteins by score; ties by protein index)
    @param excluded proteins (indices) that are not proposed
    @param candidates output, best first
  */
  void searchSlidingLadders(const std::vector<double>& masses, double tol_ppm, Size top_n, const std::set<Size>& excluded,
                            std::vector<Candidate>& candidates) const;

  /// searchSlidingLadders for many spectra (@p masses and @p excluded per spectrum) at once: the database is read once per batch of 64
  /// spectra instead of once per spectrum; candidates[s] is what the search of spectrum s alone gives
  void searchSlidingLadders(const std::vector<const std::vector<double>*>& masses, double tol_ppm, Size top_n,
                            const std::vector<const std::set<Size>*>& excluded, std::vector<std::vector<Candidate>>& candidates) const;

  /// the number of runs of standard residues (segments) the precursor-ladder windows lie in
  Size getSegmentCount() const
  {
    return segments_.size();
  }

private:
  /// a residue sum from a terminus and its protein, in the mass-sorted terminal ladder index
  struct LadderEntry
  {
    double mass;
    std::uint32_t protein;
  };

  /// per ion type: whether it is a prefix ion, and the shift that makes its fragment mass a residue sum (mass + shift)
  std::vector<std::pair<bool, double>> ion_residue_shifts_;
  /// the protein lengths
  std::vector<Size> protein_lengths_;
  /// with a database of several entries per protein: per entry (sequence) its protein, its number of fixed modifications and its offset in
  /// the protein, by which the searches keep one entry per protein (bestEntryPerProtein_); empty: none
  std::vector<Size> entry_protein_, entry_fixed_count_;
  std::vector<int> entry_offset_;
  /// per entry, the fixed masses in its terminal ladder sums (ladder_pre_, ladder_suf_; the chance of a match is taken from the residues
  /// without them), and the running sum of its fixed masses by residue for the window sums (sliding ladders); empty: no fixed masses
  std::vector<double> ladder_pre_fixed_, ladder_suf_fixed_;
  std::vector<std::vector<double>> window_fixed_sums_;
  /// precursor ladders: the windows of several entries that are one hypothesis
  WindowOwnership window_owners_;
  /// the evidence of a ranked entry for bestEntryPerProtein_: the scores of both sides summed, and the matches of each side
  struct EntryEvidence
  {
    double score_sum = 0;
    int n = 0, c = 0;
  };
  /// Keep one candidate entry per protein of @p ranked (score, entry), as entry_protein_ says: the entry with the best scores of both
  /// sides summed (@p evidence, parallel to @p ranked; then fewer fixed modifications, then the lower index), or, of its processed form,
  /// the entry with as many matches on each side and the fewest fixed modifications (then the lowest index); the kept one carries the
  /// protein's best score (FLASHTnTAlgorithm's candidate collapse).
  void bestEntryPerProtein_(std::vector<std::pair<double, Size>>& ranked, const std::vector<EntryEvidence>& evidence) const;
  /// After the cut to the top proteins (@p ranked: one entry per protein, the first @p top kept): every entry of those proteins in
  /// @p all (the entries with evidence, before bestEntryPerProtein_) follows its protein's entry, so that each is extended; @p top grows
  /// to the number kept.
  void addEntriesOfTopProteins_(std::vector<std::pair<double, Size>>& ranked, Size& top, const std::vector<std::pair<double, Size>>& all) const;
  /// per protein, the residue sums from its N- and from its C-terminus (up to the first non-standard residue): the terminal ladders
  /// (ladder_suf_ only with terminal ladders); ladder_pre_ is also the window sums of a protein of standard residues only
  std::vector<std::vector<double>> ladder_pre_, ladder_suf_;
  /// terminal ladders: all of them sorted by mass (a sum of 0 left out)
  std::vector<LadderEntry> ladder_pre_index_, ladder_suf_index_;
  /// terminal ladders: per Da of mass from 0, the first entry of the index at or above it (empty: not built, e.g. beyond 10 MDa)
  std::vector<std::uint32_t> ladder_pre_first_, ladder_suf_first_;
  /// terminal ladders: observed mass + shift = residue sum from the terminus
  std::vector<double> n_ladder_shifts_, c_ladder_shifts_;
  /// precursor ladders: the runs of standard residues of the database
  std::vector<Segment> segments_;
  /// precursor ladders: per protein, its residue sums from the N-terminus, in which a non-standard residue adds 0 (no window crosses it)
  std::vector<const double*> window_sums_;
  /// precursor ladders: the storage of window_sums_ (the sums of a protein of standard residues only are its ladder_pre_ sums)
  std::vector<std::vector<double>> window_own_sums_;
  /// precursor ladders, preparePrecursorWindows: the precursor tolerance, the prepared masses (ascending, distinct) and their windows, by
  /// chunk of segments and in a chunk by mass (in the order of the sweep): those of mass n in chunk c start at
  /// prepared_offsets_[n * prepared_chunks_ + c] (a row more than masses for the ends)
  double prepared_tol_ppm_ = -1;
  std::vector<double> prepared_masses_;
  Size prepared_chunks_ = 0;
  std::unique_ptr<PreparedWindow[]> prepared_windows_; // not value-initialized: written once, by the threads that fill it
  Size prepared_window_count_ = 0;
  std::vector<std::uint32_t> prepared_offsets_;
};
} // namespace OpenMS
