// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <limits>
#include <vector>

namespace OpenMS
{
/**
@brief Annotates proteoform-spectrum matches (PrSMs) on the raw (centroided, not deconvolved) MS2 spectrum.

The theoretical isotope envelope of every terminal fragment of the proteoform is placed at every plausible charge on
the centroided spectrum, the way ExD Viewer annotates, and kept if it matches:
 - isotope peaks: averagine (FLASHHelperClasses::PrecalculatedAveragine, 25 Da bins) peaks of at least
   min_rel_isotope_intensity of the most abundant one, at mono + k * C13C12_MASSDIFF_U;
 - charges 1 .. min(precursor charge, floor(fragment mass / min_mass_per_charge));
 - each predicted peak is matched to the nearest observed peak within tol ppm (after calibration);
 - an envelope is accepted when at least min_peaks peaks, the most abundant one included, are matched and the cosine and
   the Pearson correlation between predicted and observed intensities (unmatched = 0) pass min_cos and min_pearson.

Blind modifications: a fragment carries a modification mass only if it contains the whole modification range;
fragments that split a range are not annotated.

Calibration (optional): a first pass at calibration_tol over the b/y fragments takes, for every envelope with at least
min_peaks matched peaks, the ppm error of the monoisotopic peak (if kept and matched) or else of the most abundant one;
the median error corrects the observed m/z. Fewer than 5 such envelopes: no correction.

Score: a raw match count grows with the number of fragments tried and with the density of the spectrum. The chance
match rate p0 is therefore measured per PrSM, in its own spectrum, with the same rule and charges on decoy fragment sets
(every target envelope, i.e. its isotope pattern and charges, moved by each of decoy_shifts):
p0 = (decoy matched hypotheses + 1) / (K * n + 2), with K the number of shifts and n the fragment hypotheses. A hypothesis
is an ion type and a length with at least min_peaks predicted isotope peaks and at least one charge whose most abundant
peak falls in the spectrum m/z range; ion types whose masses differ by one to three isotope spacings (z, z+1 and z+2) put
their envelopes on the same peaks and form one hypothesis. A hypothesis counts once, whatever the number of its matched
charges and ion types, and decoys are counted only for the n target hypotheses. With k matched hypotheses,
 - score = -log10 P(X >= k), X ~ Binomial(n, p0): the length enters only through n, so the score does not grow with the
   proteoform length. It is a p-value only as far as the decoy fragments match as often as random fragments do. The
   default shifts are whole multiples of 1.000495 Da, the spacing of the peptide mass bands, so that the decoy masses lie
   in the bands like the masses of random fragments; shifts between the bands underestimate p0;
 - llr = k ln(p1 / p0) + (n - k) ln((1 - p1) / (1 - p0)).
Both are evidence scores of one PrSM in its spectrum, not identification p-values: a candidate selected by a search
(tags, blind modifications) matches more than a random proteoform even when it is wrong. For identification, use
target-decoy competition, e.g. the scores of decoy-protein PrSMs.

@ingroup Topdown
*/
class OPENMS_DLLAPI FLASHAnnotatorAlgorithm : public DefaultParamHandler
{
public:
  /// a centroided spectrum as flat arrays, sorted by m/z, without zero-intensity peaks
  struct OPENMS_DLLAPI PeakArrays
  {
    std::vector<double> mzs;
    std::vector<double> intensities;
  };

  /// a blind modification: its mass and its 0-based inclusive residue range in the proteoform
  struct OPENMS_DLLAPI Modification
  {
    double mass = 0;
    int start = 0;
    int end = 0;
  };

  /// a PrSM to annotate
  struct OPENMS_DLLAPI PrSM
  {
    /// proteoform sequence (one-letter codes, without modifications)
    String sequence;
    std::vector<Modification> modifications;
    /// fragment charges are capped by this charge (<= 0: only by min_mass_per_charge)
    int precursor_charge = 0;
    /// if false (undetermined start or end position), prefix or suffix fragments are not annotated
    bool known_n_term = true;
    bool known_c_term = true;
    /// index of the spectrum in the spectra given to annotate(); negative if the spectrum is missing
    int spectrum_index = -1;
  };

  /// a terminal fragment hypothesis
  struct OPENMS_DLLAPI Fragment
  {
    /// index into getIonTypes()
    int ion_type = 0;
    /// number of residues
    int length = 0;
    /// the bond (1 .. L-1) the fragment supports
    int bond = 0;
    /// neutral monoisotopic mass
    double mono_mass = 0;
    /// highest charge tried
    int max_charge = 0;
  };

  /// an accepted isotope envelope
  struct OPENMS_DLLAPI Envelope
  {
    int ion_type = 0;
    int length = 0;
    int bond = 0;
    int charge = 0;
    double mono_mass = 0;
    /// theoretical m/z of the most abundant isotope peak
    double most_abundant_mz = 0;
    int matched_peaks = 0;
    int predicted_peaks = 0;
    double cosine = 0;
    double pearson = 0;
    /// summed intensity of the matched peaks
    double intensity = 0;
    /// ppm error of the most abundant peak after calibration
    double ppm_error = 0;
  };

  /**
    annotation of a PrSM. Values are -1 when not computed (no spectrum, unknown residue; for the scores also no decoy shifts
    or no fragment hypothesis), except llr and calibration_ppm, which can be negative and are then NaN.
  */
  struct OPENMS_DLLAPI Result
  {
    bool annotated = false;
    std::vector<Envelope> envelopes;
    int covered_bonds = -1;
    /// covered bonds / (L-1) * 100
    double coverage = -1;
    /// k: fragment hypotheses (ion type group, length) with an accepted envelope at any charge
    int matched_fragments = -1;
    /// n: fragment hypotheses (ion type group, length) that can be matched (enough predicted peaks, a charge in the spectrum m/z range)
    int fragment_hypotheses = -1;
    /// matched decoy hypotheses, summed over the decoy shifts
    int decoy_matched_fragments = -1;
    /// p0
    double chance_rate = -1;
    /// -log10 P(X >= k), X ~ Binomial(n, p0)
    double score = -1;
    double llr = std::numeric_limits<double>::quiet_NaN();
    double calibration_ppm = std::numeric_limits<double>::quiet_NaN();
    /// number of envelopes the calibration was computed from
    int calibration_envelopes = 0;
  };

  /// constructor
  FLASHAnnotatorAlgorithm();

  /// destructor
  ~FLASHAnnotatorAlgorithm() override = default;

  /// sorted flat arrays of the peaks of @p spec with positive intensity
  static void toPeakArrays(const MSSpectrum& spec, PeakArrays& peaks);

  /**
    @brief Annotate @p prsms on @p spectra (centroided) in parallel over the PrSMs.

    The averagine table is computed at the first call, up to the heaviest fragment (at most max_mass), and extended only
    if a later call needs heavier fragments.
    @param spectra spectra referred to by PrSM::spectrum_index
    @param prsms the PrSMs
    @param results one result per PrSM, in the same order
  */
  void annotate(const std::vector<PeakArrays>& spectra, const std::vector<PrSM>& prsms, std::vector<Result>& results);

  /// the terminal fragments of @p prsm that are annotated; false if the sequence contains an unknown residue
  bool getFragments(const PrSM& prsm, std::vector<Fragment>& fragments) const;

  /// theoretical m/z and relative intensities (most abundant = 1) of the kept isotope peaks of a fragment
  void getIsotopeEnvelope(double mono_mass, int charge, std::vector<double>& mzs, std::vector<double>& intensities);

  /// the ion types, in the order Fragment::ion_type and Envelope::ion_type refer to
  const std::vector<String>& getIonTypes() const;

  /// neutral monoisotopic mass of the proteoform (residues + water + modifications); -1 for an unknown residue
  static double getProteoformMass(const String& sequence, const std::vector<Modification>& modifications);

  /// -log10 P(X >= k) for X ~ Binomial(n, p), computed in log space and capped at max_score_
  static double binomialTailScore(int n, int k, double p);

protected:
  void updateMembers_() override;
  /// implemented for DefaultParamHandler
  void setDefaultParams_();

private:
  /// the flattened averagine table: the kept isotope peaks of every 25 Da mass bin
  struct IsotopeTable
  {
    /// bins cover masses min_mass_ + b * bin_width_
    double max_mass = 0;
    /// begin of the kept peaks of a bin in iso_index / iso_rel; size = number of bins + 1
    std::vector<int> offset;
    /// isotope index k of a kept peak
    std::vector<int> iso_index;
    /// intensity relative to the most abundant peak
    std::vector<double> iso_rel;
    /// per bin: position (in the kept peaks of the bin) of the most abundant peak and of the monoisotopic peak (-1 if not kept)
    std::vector<int> apex_pos;
    std::vector<int> mono_pos;
    /// per bin: sum and sum of squares of the relative intensities (cosine and Pearson)
    std::vector<double> rel_sum;
    std::vector<double> rel_sum_sq;
  };

  /// match statistics of one envelope
  struct EnvelopeMatch
  {
    int matched = 0;
    int apex_peak = -1;
    double cosine = 0;
    double pearson = 0;
    double intensity = 0;
  };

  /// make sure the isotope table covers @p mass (computed in parallel, once per needed mass range)
  void ensureIsotopeTable_(double mass);
  /// averagine mass bin of a monoisotopic mass
  int massToBin_(double mass) const;
  /**
    Test the envelope of mass @p mono_mass at charge @p z with the isotope peaks of bin @p bin, the predicted m/z scaled by
    @p cal_factor (1 + calibration ppm * 1e-6). True if accepted; @p match is then filled.
  */
  bool matchEnvelope_(const PeakArrays& spec, double mono_mass, int bin, int z, double cal_factor, EnvelopeMatch& match) const;
  /// the median ppm error of the calibration envelopes, or 0 if there are fewer than min_calibration_envelopes_
  double calibrate_(const PeakArrays& spec, const std::vector<Fragment>& fragments, int& count) const;
  /// annotate one PrSM (without the scores, which are computed afterwards)
  void annotateOne_(const PeakArrays& spec, const PrSM& prsm, Result& result) const;

  /// parameters
  double tol_ = 10.0;
  double cal_tol_ = 20.0;
  bool calibrate_on_ = true;
  int min_peaks_ = 3;
  double min_rel_iso_ = 0.1;
  double min_cos_ = 0.8;
  double min_pearson_ = 0.7;
  double min_mass_per_charge_ = 300.0;
  double max_mass_ = 200000.0;
  double p1_ = 0.3;
  std::vector<double> decoy_shifts_;
  std::vector<String> ion_types_;
  /// mass added to the residue sum of a fragment, per ion type; and whether the ion type is a prefix (N-terminal) one
  std::vector<double> ion_shifts_;
  std::vector<bool> ion_is_prefix_;
  /// used for calibration (b and y if present, else all)
  std::vector<bool> ion_for_calibration_;
  /// per ion type: its group of ion types whose envelopes fall on the same peaks (z, z+1, z+2), and the number of groups
  std::vector<int> ion_group_;
  int ion_group_count_ = 0;

  IsotopeTable table_;
  double table_min_rel_iso_ = -1;

  /// averagine bins as in SpectralDeconvolution::calculateAveragine
  static constexpr double min_mass_ = 50.0;
  static constexpr double bin_width_ = 25.0;
  /// fewer calibration envelopes than this: no correction
  static constexpr int min_calibration_envelopes_ = 5;
  /// cap of the score
  static constexpr double max_score_ = 1000.0;
  /// ion types of the same terminus whose shifts differ by 1 .. max_isotope_offset_ isotope spacings within this mass (Da)
  /// match the same envelopes (z, z+1 and z+2 differ by 4.5 mDa per spacing) and are one fragment hypothesis
  static constexpr double isotope_equivalence_tol_ = 0.01;
  static constexpr int max_isotope_offset_ = 3;
};
} // namespace OpenMS
