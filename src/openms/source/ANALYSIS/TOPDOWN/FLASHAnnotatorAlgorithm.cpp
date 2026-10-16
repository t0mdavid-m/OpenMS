// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHAnnotatorAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numeric>

namespace OpenMS
{
namespace
{
  /// internal (residue) monoisotopic masses by one-letter code; 0 for codes without a residue of positive mass
  const std::vector<double>& residueMasses()
  {
    static const std::vector<double> masses = [] {
      std::vector<double> m(256, 0.0);
      const ResidueDB* db = ResidueDB::getInstance();
      for (int c = 'A'; c <= 'Z'; c++)
      {
        const Residue* r = db->getResidue((unsigned char)c);
        if (r == nullptr) continue;
        double w = r->getMonoWeight(Residue::Internal);
        if (w <= 0) continue;
        m[c] = w;
        m[c - 'A' + 'a'] = w;
      }
      return m;
    }();
    return masses;
  }

  /**
    Index of the observed peak nearest to @p target within [target - w, target + w], or -1 (the first of equally near
    peaks wins). The first index with m/z >= target - w is searched in [lo, hi), which must contain it (or hi if it is
    hi), and returned in @p lb so that the next isotope peak can be searched in a narrower range.
  */
  inline int nearestPeak(const std::vector<double>& mzs, int lo, int hi, double target, double w, int& lb)
  {
    const double* base = mzs.data();
    const int size = (int)mzs.size();
    lb = (int)(std::lower_bound(base + lo, base + hi, target - w) - base);
    int best = -1;
    double best_d = 0;
    for (int i = lb; i < size && base[i] <= target + w; i++)
    {
      double d = std::abs(base[i] - target);
      if (best < 0 || d < best_d)
      {
        best = i;
        best_d = d;
      }
    }
    return best;
  }
} // namespace

FLASHAnnotatorAlgorithm::FLASHAnnotatorAlgorithm(): DefaultParamHandler("FLASHAnnotatorAlgorithm")
{
  setDefaultParams_();
}

void FLASHAnnotatorAlgorithm::setDefaultParams_()
{
  defaults_.setValue("tol", 10.0, "Tolerance (ppm) for matching isotope peaks to observed peaks, after calibration.");
  defaults_.setMinFloat("tol", 0.0);

  defaults_.setValue("calibrate", "true",
                     "Correct the observed m/z by the median ppm error of the fragment envelopes matched at calibration_tol "
                     "(computed per PrSM, from its b/y fragments).");
  defaults_.setValidStrings("calibrate", {"true", "false"});

  defaults_.setValue("calibration_tol", 20.0, "Tolerance (ppm) of the calibration pass.");
  defaults_.setMinFloat("calibration_tol", 0.0);

  defaults_.setValue("ion_type", std::vector<std::string> {"b", "y"}, "Ion types to annotate.");
  defaults_.setValidStrings("ion_type", {"b", "c", "a", "y", "z", "x", "zp1", "zp2"});

  defaults_.setValue("min_peaks", 3, "Minimum number of matched isotope peaks of an envelope, the most abundant one included.");
  defaults_.setMinInt("min_peaks", 1);

  defaults_.setValue("min_rel_isotope_intensity", 0.1,
                     "Isotope peaks of at least this intensity relative to the most abundant one are predicted.");
  defaults_.setMinFloat("min_rel_isotope_intensity", 0.0);
  defaults_.setMaxFloat("min_rel_isotope_intensity", 1.0);

  defaults_.setValue("min_cos", 0.8, "Minimum cosine between predicted and observed isotope intensities (unmatched peaks count as 0).");
  defaults_.setMinFloat("min_cos", -1.0);
  defaults_.setMaxFloat("min_cos", 1.0);

  defaults_.setValue("min_pearson", 0.7,
                     "Minimum Pearson correlation between predicted and observed isotope intensities (unmatched peaks count as "
                     "0). It rejects flat runs of peaks that pass the cosine.");
  defaults_.setMinFloat("min_pearson", -1.0);
  defaults_.setMaxFloat("min_pearson", 1.0);

  defaults_.setValue("min_mass_per_charge", 300.0, "Fragment charges are capped at floor(fragment mass / this value).");
  defaults_.setMinFloat("min_mass_per_charge", 1.0);

  defaults_.setValue("max_mass", 200000.0, "Fragments heavier than this are not annotated.");
  defaults_.setMinFloat("max_mass", 100.0);

  // whole multiples (9, 13, 21, 33, 39, 53) of 1.000495 Da, the spacing of the peptide mass bands
  defaults_.setValue("decoy_shifts",
                     std::vector<double> {-53.026235, -39.019305, -33.016335, -21.010395, -13.006435, -9.004455, 9.004455, 13.006435,
                                          21.010395, 33.016335, 39.019305, 53.026235},
                     "Mass shifts (Da) of the decoy fragment sets that measure the chance match rate of each PrSM in its own "
                     "spectrum. Whole multiples of 1.000495 Da keep the decoy masses in the peptide mass bands, where random "
                     "fragments lie; shifts between the bands match less often and underestimate the chance rate. Empty: no "
                     "chance rate, no RawScore and no RawLLR.");

  defaults_.setValue("p1", 0.3, "Probability that a fragment of a correct proteoform is observed, for the log-likelihood ratio.");
  defaults_.setMinFloat("p1", 1e-6);
  defaults_.setMaxFloat("p1", 1.0 - 1e-6);

  defaultsToParam_();
}

void FLASHAnnotatorAlgorithm::updateMembers_()
{
  tol_ = param_.getValue("tol");
  calibrate_on_ = param_.getValue("calibrate").toString() == "true";
  cal_tol_ = param_.getValue("calibration_tol");
  min_peaks_ = param_.getValue("min_peaks");
  min_rel_iso_ = param_.getValue("min_rel_isotope_intensity");
  min_cos_ = param_.getValue("min_cos");
  min_pearson_ = param_.getValue("min_pearson");
  min_mass_per_charge_ = param_.getValue("min_mass_per_charge");
  max_mass_ = param_.getValue("max_mass");
  decoy_shifts_ = param_.getValue("decoy_shifts").toDoubleVector();
  p1_ = param_.getValue("p1");

  // same fragment masses as FLASHTaggerAlgorithm and FLASHExtenderAlgorithm: internal residue sum + ion type shift
  ion_types_.clear();
  ion_shifts_.clear();
  ion_is_prefix_.clear();
  for (const auto& ion_str : param_.getValue("ion_type").toStringVector())
  {
    if (std::find(ion_types_.begin(), ion_types_.end(), ion_str) != ion_types_.end()) continue;
    double shift;
    bool prefix;
    if (ion_str == "a") { shift = Residue::getInternalToAIon().getMonoWeight(); prefix = true; }
    else if (ion_str == "b") { shift = Residue::getInternalToBIon().getMonoWeight(); prefix = true; }
    else if (ion_str == "c") { shift = Residue::getInternalToCIon().getMonoWeight(); prefix = true; }
    else if (ion_str == "x") { shift = Residue::getInternalToXIon().getMonoWeight(); prefix = false; }
    else if (ion_str == "y") { shift = Residue::getInternalToYIon().getMonoWeight(); prefix = false; }
    else if (ion_str == "z") { shift = Residue::getInternalToZIon().getMonoWeight(); prefix = false; }
    else if (ion_str == "zp1") { shift = Residue::getInternalToZp1Ion().getMonoWeight(); prefix = false; }
    else if (ion_str == "zp2") { shift = Residue::getInternalToZp2Ion().getMonoWeight(); prefix = false; }
    else { continue; }
    ion_types_.push_back(ion_str);
    ion_shifts_.push_back(shift);
    ion_is_prefix_.push_back(prefix);
  }

  // calibration uses the b and y fragments (as the prototype did); all ion types if neither is given
  bool has_by = std::find(ion_types_.begin(), ion_types_.end(), "b") != ion_types_.end()
                || std::find(ion_types_.begin(), ion_types_.end(), "y") != ion_types_.end();
  ion_for_calibration_.clear();
  for (const auto& ion : ion_types_)
  {
    ion_for_calibration_.push_back(! has_by || ion == "b" || ion == "y");
  }

  // ion types of the same terminus whose shifts differ by a few isotope spacings (z, z+1, z+2) match the same envelope one
  // or two isotope peaks apart: they form one fragment hypothesis for the score
  const int ion_count = (int)ion_types_.size();
  ion_group_.assign(ion_count, -1);
  ion_group_count_ = 0;
  for (int t = 0; t < ion_count; t++)
  {
    if (ion_group_[t] >= 0) continue;
    ion_group_[t] = ion_group_count_;
    // grow the group of t: an ion type joins it if it is equivalent to any member
    for (bool grown = true; grown;)
    {
      grown = false;
      for (int u = t + 1; u < ion_count; u++)
      {
        if (ion_group_[u] >= 0 || ion_is_prefix_[u] != ion_is_prefix_[t]) continue;
        for (int v = t; v < ion_count; v++)
        {
          if (ion_group_[v] != ion_group_count_) continue;
          const double diff = ion_shifts_[u] - ion_shifts_[v];
          const int offset = (int)std::round(diff / Constants::C13C12_MASSDIFF_U);
          if (offset == 0 || std::abs(offset) > max_isotope_offset_) continue;
          if (std::abs(diff - offset * Constants::C13C12_MASSDIFF_U) > isotope_equivalence_tol_) continue;
          ion_group_[u] = ion_group_count_;
          grown = true;
          break;
        }
      }
    }
    ion_group_count_++;
  }
}

const std::vector<String>& FLASHAnnotatorAlgorithm::getIonTypes() const
{
  return ion_types_;
}

void FLASHAnnotatorAlgorithm::toPeakArrays(const MSSpectrum& spec, PeakArrays& peaks)
{
  peaks.mzs.clear();
  peaks.intensities.clear();
  peaks.mzs.reserve(spec.size());
  peaks.intensities.reserve(spec.size());

  if (spec.isSorted())
  {
    for (const auto& p : spec)
    {
      if (p.getIntensity() <= 0) continue;
      peaks.mzs.push_back(p.getMZ());
      peaks.intensities.push_back(p.getIntensity());
    }
    return;
  }

  std::vector<std::pair<double, double>> sorted;
  sorted.reserve(spec.size());
  for (const auto& p : spec)
  {
    if (p.getIntensity() <= 0) continue;
    sorted.emplace_back(p.getMZ(), p.getIntensity());
  }
  std::sort(sorted.begin(), sorted.end());
  for (const auto& [mz, intensity] : sorted)
  {
    peaks.mzs.push_back(mz);
    peaks.intensities.push_back(intensity);
  }
}

double FLASHAnnotatorAlgorithm::getProteoformMass(const String& sequence, const std::vector<Modification>& modifications)
{
  static const double water = Residue::getInternalToFull().getMonoWeight();
  const auto& aa = residueMasses();
  double mass = water;
  for (const char c : sequence)
  {
    double m = aa[(unsigned char)c];
    if (m <= 0) return -1;
    mass += m;
  }
  for (const auto& mod : modifications)
  {
    mass += mod.mass;
  }
  return mass;
}

void FLASHAnnotatorAlgorithm::ensureIsotopeTable_(double mass)
{
  mass = std::min(mass, max_mass_);
  if (table_min_rel_iso_ == min_rel_iso_ && table_.max_mass >= mass && ! table_.apex_pos.empty()) return;

  // bins at min_mass_ + b * bin_width_ (50, 75, ... Da) as PrecalculatedAveragine(50, max mass, 25) in
  // SpectralDeconvolution::calculateAveragine, up to the bin nearest to mass
  const int bin_count = std::max(2, (int)std::floor((mass - min_mass_) / bin_width_) + 2);
  const double top_mass = min_mass_ + (bin_count - 1) * bin_width_;

  CoarseIsotopePatternGenerator generator(300);
  auto max_iso = generator.estimateFromPeptideWeight(top_mass);
  max_iso.trimRight(0.0001 * max_iso.getMostAbundant().getIntensity());
  const int max_isotope = std::max(200, (int)max_iso.size());

  std::vector<std::vector<int>> iso_index(bin_count);
  std::vector<std::vector<double>> iso_rel(bin_count);
  std::vector<int> apex_pos(bin_count, 0), mono_pos(bin_count, -1);

  // the distribution of a bin does not depend on the other bins, so chunks of bins are computed in parallel, each by its
  // own PrecalculatedAveragine. Each bin is then flattened into its kept peaks.
  const int chunk_size = 100;
  const int chunk_count = (bin_count + chunk_size - 1) / chunk_size;
#pragma omp parallel for schedule(dynamic)
  for (int c = 0; c < chunk_count; c++)
  {
    const int first = c * chunk_size;
    const int last = std::min(bin_count, first + chunk_size) - 1;
    CoarseIsotopePatternGenerator chunk_generator(max_isotope);
    FLASHHelperClasses::PrecalculatedAveragine avg(min_mass_ + first * bin_width_, min_mass_ + last * bin_width_, bin_width_, chunk_generator, false);
    for (int b = first; b <= last; b++)
    {
      const IsotopeDistribution iso = avg.get(min_mass_ + b * bin_width_);
      Size apex = 0;
      for (Size k = 1; k < iso.size(); k++)
      {
        if (iso[k].getIntensity() > iso[apex].getIntensity()) apex = k;
      }
      const double apex_intensity = iso.size() == 0 ? 0.0 : (double)iso[apex].getIntensity();
      if (apex_intensity <= 0) continue;
      for (Size k = 0; k < iso.size(); k++)
      {
        const double rel = iso[k].getIntensity() / apex_intensity;
        if (rel <= 0 || rel < min_rel_iso_) continue;
        if (k == apex) apex_pos[b] = (int)iso_rel[b].size();
        if (k == 0) mono_pos[b] = 0;
        iso_index[b].push_back((int)k);
        iso_rel[b].push_back(rel);
      }
    }
  }

  table_ = IsotopeTable();
  table_.max_mass = top_mass;
  table_.offset.resize(bin_count + 1, 0);
  table_.apex_pos = apex_pos;
  table_.mono_pos = mono_pos;
  table_.rel_sum.resize(bin_count, 0);
  table_.rel_sum_sq.resize(bin_count, 0);
  for (int b = 0; b < bin_count; b++)
  {
    table_.offset[b] = (int)table_.iso_index.size();
    for (Size j = 0; j < iso_rel[b].size(); j++)
    {
      table_.iso_index.push_back(iso_index[b][j]);
      table_.iso_rel.push_back(iso_rel[b][j]);
      table_.rel_sum[b] += iso_rel[b][j];
      table_.rel_sum_sq[b] += iso_rel[b][j] * iso_rel[b][j];
    }
  }
  table_.offset[bin_count] = (int)table_.iso_index.size();
  table_min_rel_iso_ = min_rel_iso_;
}

int FLASHAnnotatorAlgorithm::massToBin_(double mass) const
{
  // the same rounding as PrecalculatedAveragine::massToIndex_
  int bin = (int)std::round(std::max(.0, mass - min_mass_) / bin_width_);
  return std::min(bin, (int)table_.apex_pos.size() - 1);
}

void FLASHAnnotatorAlgorithm::getIsotopeEnvelope(double mono_mass, int charge, std::vector<double>& mzs, std::vector<double>& intensities)
{
  mzs.clear();
  intensities.clear();
  if (charge <= 0 || mono_mass <= 0) return;
  ensureIsotopeTable_(mono_mass + bin_width_);
  const int bin = massToBin_(mono_mass);
  for (int j = table_.offset[bin]; j < table_.offset[bin + 1]; j++)
  {
    mzs.push_back((mono_mass + table_.iso_index[j] * Constants::C13C12_MASSDIFF_U + charge * Constants::PROTON_MASS_U) / charge);
    intensities.push_back(table_.iso_rel[j]);
  }
}

bool FLASHAnnotatorAlgorithm::getFragments(const PrSM& prsm, std::vector<Fragment>& fragments) const
{
  fragments.clear();
  const auto& aa = residueMasses();
  const int length = (int)prsm.sequence.size();
  std::vector<double> prefix_mass(length + 1, .0);
  for (int i = 0; i < length; i++)
  {
    double m = aa[(unsigned char)prsm.sequence[i]];
    if (m <= 0) return false;
    prefix_mass[i + 1] = prefix_mass[i] + m;
  }

  // modification ranges, clamped into the proteoform
  std::vector<std::pair<int, int>> ranges;
  ranges.reserve(prsm.modifications.size());
  for (const auto& mod : prsm.modifications)
  {
    int s = std::max(0, std::min(length - 1, std::min(mod.start, mod.end)));
    int e = std::max(0, std::min(length - 1, std::max(mod.start, mod.end)));
    ranges.emplace_back(s, e);
  }

  fragments.reserve(2 * ion_types_.size() * (Size)std::max(0, length));
  for (int t = 0; t < (int)ion_types_.size(); t++)
  {
    const bool prefix = ion_is_prefix_[t];
    if ((prefix && ! prsm.known_n_term) || (! prefix && ! prsm.known_c_term)) continue;
    for (int len = 1; len < length; len++)
    {
      // residues [first, last] of the fragment
      const int first = prefix ? 0 : length - len;
      const int last = prefix ? len - 1 : length - 1;
      double mod_mass = 0;
      bool split = false;
      for (Size m = 0; m < ranges.size(); m++)
      {
        if (ranges[m].first >= first && ranges[m].second <= last) { mod_mass += prsm.modifications[m].mass; }
        else if (ranges[m].second >= first && ranges[m].first <= last)
        {
          split = true; // the fragment contains only part of the modification range: its mass is unknown
          break;
        }
      }
      if (split) continue;

      Fragment f;
      f.ion_type = t;
      f.length = len;
      f.bond = prefix ? len : length - len;
      f.mono_mass = prefix_mass[last + 1] - prefix_mass[first] + ion_shifts_[t] + mod_mass;
      if (f.mono_mass <= 0 || f.mono_mass > max_mass_) continue;
      f.max_charge = (int)std::floor(f.mono_mass / min_mass_per_charge_);
      if (prsm.precursor_charge > 0) f.max_charge = std::min(f.max_charge, prsm.precursor_charge);
      if (f.max_charge < 1) continue;
      fragments.push_back(f);
    }
  }
  return true;
}

bool FLASHAnnotatorAlgorithm::matchEnvelope_(const PeakArrays& spec, double mono_mass, int bin, int z, double cal_factor, EnvelopeMatch& match) const
{
  const int begin = table_.offset[bin];
  const int end = table_.offset[bin + 1];
  const int predicted = end - begin;
  if (predicted < min_peaks_) return false;

  const int size = (int)spec.mzs.size();
  const int apex = begin + table_.apex_pos[bin];
  // predicted m/z on the observed scale (observed / cal_factor is the calibrated m/z), and the tolerance with it
  const double scale = cal_factor / z;
  const double base = mono_mass + z * Constants::PROTON_MASS_U;
  const double tol = tol_ * 1e-6;

  // the most abundant peak first: it is required
  double target = (base + table_.iso_index[apex] * Constants::C13C12_MASSDIFF_U) * scale;
  int apex_lb;
  const int apex_peak = nearestPeak(spec.mzs, 0, size, target, target * tol, apex_lb);
  if (apex_peak < 0) return false;

  // intensities relative to the observed most abundant peak (cosine and Pearson are scale free)
  const double apex_intensity = spec.intensities[apex_peak];
  double sxy = table_.iso_rel[apex], sy = 1.0, syy = 1.0, intensity = apex_intensity;
  int matched = 1;
  int unsearched = predicted - 1;

  // isotope peaks above the most abundant one, in increasing m/z
  int lb = apex_lb;
  for (int j = apex + 1; j < end; j++)
  {
    target = (base + table_.iso_index[j] * Constants::C13C12_MASSDIFF_U) * scale;
    const int peak = nearestPeak(spec.mzs, lb, size, target, target * tol, lb);
    unsearched--;
    if (peak < 0)
    {
      if (matched + unsearched < min_peaks_) return false;
      continue;
    }
    const double y = spec.intensities[peak] / apex_intensity;
    matched++;
    sxy += table_.iso_rel[j] * y;
    sy += y;
    syy += y * y;
    intensity += spec.intensities[peak];
  }
  // and below it, in decreasing m/z
  lb = apex_lb;
  for (int j = apex - 1; j >= begin; j--)
  {
    target = (base + table_.iso_index[j] * Constants::C13C12_MASSDIFF_U) * scale;
    const int peak = nearestPeak(spec.mzs, 0, lb, target, target * tol, lb);
    unsearched--;
    if (peak < 0)
    {
      if (matched + unsearched < min_peaks_) return false;
      continue;
    }
    const double y = spec.intensities[peak] / apex_intensity;
    matched++;
    sxy += table_.iso_rel[j] * y;
    sy += y;
    syy += y * y;
    intensity += spec.intensities[peak];
  }
  if (matched < min_peaks_) return false;

  const double sx = table_.rel_sum[bin];
  const double sxx = table_.rel_sum_sq[bin];
  const double cosine = sxy / std::sqrt(sxx * syy);
  if (cosine < min_cos_) return false;

  const double n = predicted;
  const double var = (n * sxx - sx * sx) * (n * syy - sy * sy);
  const double pearson = var > 0 ? (n * sxy - sx * sy) / std::sqrt(var) : 0.0;
  if (pearson < min_pearson_) return false;

  match.matched = matched;
  match.apex_peak = apex_peak;
  match.cosine = cosine;
  match.pearson = pearson;
  match.intensity = intensity;
  return true;
}

double FLASHAnnotatorAlgorithm::calibrate_(const PeakArrays& spec, const std::vector<Fragment>& fragments, int& count) const
{
  count = 0;
  if (spec.mzs.empty()) return 0;
  const int size = (int)spec.mzs.size();
  const double min_mz = spec.mzs.front();
  const double max_mz = spec.mzs.back();
  const double tol = cal_tol_ * 1e-6;
  std::vector<double> errors;

  for (const auto& f : fragments)
  {
    if (! ion_for_calibration_[f.ion_type]) continue;
    const int bin = massToBin_(f.mono_mass);
    const int begin = table_.offset[bin];
    const int end = table_.offset[bin + 1];
    const int predicted = end - begin;
    if (predicted < min_peaks_) continue;
    const int apex = begin + table_.apex_pos[bin];
    const int mono = table_.mono_pos[bin] < 0 ? -1 : begin + table_.mono_pos[bin];

    for (int z = 1; z <= f.max_charge; z++)
    {
      const double base = f.mono_mass + z * Constants::PROTON_MASS_U;
      const double low = (base + table_.iso_index[begin] * Constants::C13C12_MASSDIFF_U) / z;
      const double high = (base + table_.iso_index[end - 1] * Constants::C13C12_MASSDIFF_U) / z;
      if (low * (1 - tol) > max_mz) continue;
      if (high * (1 + tol) < min_mz) break; // higher charges are further below

      // the error is measured on the monoisotopic peak if it is kept and matched, else on the most abundant one
      double target = (base + table_.iso_index[apex] * Constants::C13C12_MASSDIFF_U) / z;
      int apex_lb, lb;
      const int apex_peak = nearestPeak(spec.mzs, 0, size, target, target * tol, apex_lb);
      int mono_peak = -1;
      if (mono >= 0 && mono != apex)
      {
        const double mono_target = base / z;
        mono_peak = nearestPeak(spec.mzs, 0, apex_lb, mono_target, mono_target * tol, lb);
      }
      else if (mono == apex) { mono_peak = apex_peak; }
      if (apex_peak < 0 && mono_peak < 0) continue;

      int matched = (apex_peak >= 0 ? 1 : 0) + (mono >= 0 && mono != apex && mono_peak >= 0 ? 1 : 0);
      int unsearched = predicted - 1 - (mono >= 0 && mono != apex ? 1 : 0);
      bool fail = matched + unsearched < min_peaks_;
      lb = apex_lb;
      for (int j = apex + 1; j < end && ! fail; j++)
      {
        target = (base + table_.iso_index[j] * Constants::C13C12_MASSDIFF_U) / z;
        unsearched--;
        if (nearestPeak(spec.mzs, lb, size, target, target * tol, lb) >= 0) matched++;
        else fail = matched + unsearched < min_peaks_;
      }
      lb = apex_lb;
      for (int j = apex - 1; j >= begin && ! fail; j--)
      {
        if (j == mono) continue;
        target = (base + table_.iso_index[j] * Constants::C13C12_MASSDIFF_U) / z;
        unsearched--;
        if (nearestPeak(spec.mzs, 0, lb, target, target * tol, lb) >= 0) matched++;
        else fail = matched + unsearched < min_peaks_;
      }
      if (fail || matched < min_peaks_) continue;

      const int ref = mono_peak >= 0 ? mono : apex;
      const int ref_peak = mono_peak >= 0 ? mono_peak : apex_peak;
      const double expected = (base + table_.iso_index[ref] * Constants::C13C12_MASSDIFF_U) / z;
      errors.push_back((spec.mzs[ref_peak] - expected) / expected * 1e6);
    }
  }

  count = (int)errors.size();
  if (count < min_calibration_envelopes_) return 0;
  const Size mid = errors.size() / 2;
  std::nth_element(errors.begin(), errors.begin() + mid, errors.end());
  double median = errors[mid];
  if (errors.size() % 2 == 0)
  {
    median = (median + *std::max_element(errors.begin(), errors.begin() + mid)) / 2.0;
  }
  return median;
}

void FLASHAnnotatorAlgorithm::annotateOne_(const PeakArrays& spec, const PrSM& prsm, Result& result) const
{
  std::vector<Fragment> fragments;
  if (! getFragments(prsm, fragments)) return;

  const int length = (int)prsm.sequence.size();
  result.annotated = true;
  result.envelopes.clear();

  int cal_count = 0;
  const double cal_ppm = calibrate_on_ ? calibrate_(spec, fragments, cal_count) : 0.0;
  const double cal_factor = 1.0 + cal_ppm * 1e-6;
  result.calibration_ppm = cal_ppm;
  result.calibration_envelopes = cal_count;

  const bool empty = spec.mzs.empty();
  const double min_mz = empty ? 0 : spec.mzs.front();
  const double max_mz = empty ? 0 : spec.mzs.back();
  const double tol = tol_ * 1e-6;

  // fragment hypotheses (ion type group, length): ion types whose envelopes fall on the same peaks count once
  const int hypothesis_count = ion_group_count_ * std::max(1, length);
  std::vector<char> testable(hypothesis_count, 0), matched(hypothesis_count, 0);

  // target fragments: every accepted envelope is recorded
  std::vector<char> covered(std::max(1, length), 0);
  EnvelopeMatch match;
  for (const auto& f : fragments)
  {
    if (empty) break; // nothing can match; no hypothesis is testable
    const int bin = massToBin_(f.mono_mass);
    const int predicted = table_.offset[bin + 1] - table_.offset[bin];
    if (predicted < min_peaks_) continue; // can never be accepted, not a hypothesis
    const double apex_shift = table_.iso_index[table_.offset[bin] + table_.apex_pos[bin]] * Constants::C13C12_MASSDIFF_U;
    const int hypothesis = ion_group_[f.ion_type] * length + f.length;
    for (int z = 1; z <= f.max_charge; z++)
    {
      const double apex_mz = (f.mono_mass + apex_shift + z * Constants::PROTON_MASS_U) / z;
      const double observed = apex_mz * cal_factor;
      if (observed * (1 - tol) > max_mz) continue;
      if (observed * (1 + tol) < min_mz) break; // higher charges are further below
      testable[hypothesis] = 1;
      if (! matchEnvelope_(spec, f.mono_mass, bin, z, cal_factor, match)) continue;
      matched[hypothesis] = 1;
      covered[f.bond] = 1;

      Envelope e;
      e.ion_type = f.ion_type;
      e.length = f.length;
      e.bond = f.bond;
      e.charge = z;
      e.mono_mass = f.mono_mass;
      e.most_abundant_mz = apex_mz;
      e.matched_peaks = match.matched;
      e.predicted_peaks = predicted;
      e.cosine = match.cosine;
      e.pearson = match.pearson;
      e.intensity = match.intensity;
      e.ppm_error = (spec.mzs[match.apex_peak] / cal_factor - apex_mz) / apex_mz * 1e6;
      result.envelopes.push_back(e);
    }
  }
  const int hypotheses = (int)std::count(testable.begin(), testable.end(), 1);
  const int matched_hypotheses = (int)std::count(matched.begin(), matched.end(), 1);

  // decoy fragments: each target fragment (its isotope pattern and charges) moved by a mass shift; counts only, and only
  // for the testable target hypotheses, so that at most K * n decoy hypotheses match
  int decoy_matched = 0;
  std::vector<char> decoy_matched_hypothesis(hypothesis_count, 0);
  for (Size s = 0; s < decoy_shifts_.size() && hypotheses > 0; s++)
  {
    const double shift = decoy_shifts_[s];
    std::fill(decoy_matched_hypothesis.begin(), decoy_matched_hypothesis.end(), 0);
    for (const auto& f : fragments)
    {
      const int hypothesis = ion_group_[f.ion_type] * length + f.length;
      if (! testable[hypothesis] || decoy_matched_hypothesis[hypothesis]) continue;
      const double mono_mass = f.mono_mass + shift;
      if (mono_mass <= 0) continue;
      const int bin = massToBin_(f.mono_mass);
      if (table_.offset[bin + 1] - table_.offset[bin] < min_peaks_) continue;
      const double apex_shift = table_.iso_index[table_.offset[bin] + table_.apex_pos[bin]] * Constants::C13C12_MASSDIFF_U;
      for (int z = 1; z <= f.max_charge; z++)
      {
        const double observed = (mono_mass + apex_shift + z * Constants::PROTON_MASS_U) / z * cal_factor;
        if (observed * (1 - tol) > max_mz) continue;
        if (observed * (1 + tol) < min_mz) break;
        if (matchEnvelope_(spec, mono_mass, bin, z, cal_factor, match))
        {
          decoy_matched_hypothesis[hypothesis] = 1; // a hypothesis counts once, whatever the number of charges and ion types
          decoy_matched++;
          break;
        }
      }
    }
  }

  int covered_bonds = 0;
  for (int b = 1; b < length; b++)
  {
    covered_bonds += covered[b];
  }
  result.covered_bonds = covered_bonds;
  result.coverage = length > 1 ? 100.0 * covered_bonds / (length - 1) : 0.0;
  result.matched_fragments = matched_hypotheses;
  result.fragment_hypotheses = hypotheses;
  result.decoy_matched_fragments = decoy_matched;
}

double FLASHAnnotatorAlgorithm::binomialTailScore(int n, int k, double p)
{
  if (k <= 0 || n <= 0) return 0;
  if (k > n) k = n;
  if (p <= 0) return max_score_;
  if (p >= 1) return 0;

  const double log_p = std::log(p);
  const double log_q = std::log1p(-p);
  const double log_n_fact = std::lgamma(n + 1.0);
  auto log_pmf = [&](int j) { return log_n_fact - std::lgamma(j + 1.0) - std::lgamma(n - j + 1.0) + j * log_p + (n - j) * log_q; };

  // the terms of the upper tail grow up to the mode and decrease after it; sum them relative to the largest one
  const int mode = std::min(n, (int)std::floor((n + 1) * p));
  const int top = std::max(k, mode);
  const double log_top = log_pmf(top);
  double sum = 0;
  for (int j = k; j <= n; j++)
  {
    const double d = log_pmf(j) - log_top;
    if (j > top && d < -50) break; // the rest is negligible
    sum += std::exp(d);
  }
  const double score = -(log_top + std::log(sum)) / std::log(10.0);
  return std::max(0.0, std::min(max_score_, score));
}

void FLASHAnnotatorAlgorithm::annotate(const std::vector<PeakArrays>& spectra, const std::vector<PrSM>& prsms, std::vector<Result>& results)
{
  results.assign(prsms.size(), Result());

  // the isotope table must cover the heaviest fragment. It is computed here, once, before the parallel region.
  const auto& aa = residueMasses();
  double max_fragment_mass = 0;
  for (const auto& prsm : prsms)
  {
    double mass = 0;
    for (const char c : prsm.sequence)
    {
      mass += aa[(unsigned char)c];
    }
    for (const auto& mod : prsm.modifications)
    {
      mass += std::max(.0, mod.mass);
    }
    max_fragment_mass = std::max(max_fragment_mass, mass + 100.0); // + the largest ion type shift
  }
  ensureIsotopeTable_(max_fragment_mass);

  // longest proteoforms first, so that the dynamic schedule does not end with a long one
  std::vector<int> order(prsms.size());
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(), [&prsms](int a, int b) { return prsms[a].sequence.size() > prsms[b].sequence.size(); });

  const int spectrum_count = (int)spectra.size();
#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < (int)order.size(); i++)
  {
    const int index = order[i];
    const PrSM& prsm = prsms[index];
    if (prsm.spectrum_index < 0 || prsm.spectrum_index >= spectrum_count) continue;
    annotateOne_(spectra[prsm.spectrum_index], prsm, results[index]);
  }

  // scores (serial: std::lgamma need not be thread safe). Without a testable hypothesis there is nothing to score.
  const double decoy_sets = (double)decoy_shifts_.size();
  for (auto& r : results)
  {
    if (! r.annotated || decoy_shifts_.empty() || r.fragment_hypotheses <= 0) continue;
    const int n = r.fragment_hypotheses;
    const int k = r.matched_fragments;
    const double p0 = (r.decoy_matched_fragments + 1.0) / (decoy_sets * n + 2.0);
    r.chance_rate = p0;
    r.score = binomialTailScore(n, k, p0);
    r.llr = k * std::log(p1_ / p0) + (n - k) * std::log((1.0 - p1_) / (1.0 - p0));
  }
}
} // namespace OpenMS
