// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong, Jihyung Kim $
// $Authors: Kyowon Jeong, Jihyung Kim $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <cmath>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroupScoring.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <algorithm>
#include <functional>
#include <array>
#include <bit>
#include <cstdint>

namespace OpenMS
{
  PeakGroup::PeakGroup(const int min_abs_charge, const int max_abs_charge, const bool is_positive):
      min_abs_charge_(min_abs_charge),
      max_abs_charge_(max_abs_charge),
      is_positive_(is_positive)
  {
  }

  bool PeakGroup::operator<(const PeakGroup& a) const
  {
    if (this->monoisotopic_mass_ == a.monoisotopic_mass_) { return this->intensity_ < a.intensity_; }
    return this->monoisotopic_mass_ < a.monoisotopic_mass_;
  }

  bool PeakGroup::operator>(const PeakGroup& a) const
  {
    if (this->monoisotopic_mass_ == a.monoisotopic_mass_) { return this->intensity_ > a.intensity_; }
    return this->monoisotopic_mass_ > a.monoisotopic_mass_;
  }

  bool PeakGroup::operator==(const PeakGroup& a) const
  {
    return this->monoisotopic_mass_ == a.monoisotopic_mass_ && this->intensity_ == a.intensity_;
  }

  std::vector<float> PeakGroup::getMassErrors(bool ppm) const
  {
    std::vector<float> errors;
    int i = 0;
    float i_cntr = 0;
    float i_error = .0f;
    for (const auto& p : *this)
    {
      if (i != p.isotopeIndex)
      {
        if (i_cntr > 0)
        {
          errors.push_back(i_error / i_cntr);
        }
        i = p.isotopeIndex;
        i_cntr = 0;
        i_error = 0;

      }
      i_cntr ++;
      i_error += ppm ? getPPMError_(p) : getDaError_(p);
    }
    if (i_cntr > 0)
    {
      errors.push_back(i_error / i_cntr);
    }
    return errors;
  }


  void PeakGroup::updateAvgMassError_()
  {
    avg_ppm_error_ = averageAbsMassError_(true);
    avg_da_error_ = averageAbsMassError_(false);
  }

  float PeakGroup::averageAbsMassError_(bool ppm) const
  {
    // same arithmetic as averaging the absolute values of getMassErrors(ppm), without building the vector
    float sum = .0f;
    Size count = 0;
    int i = 0;
    float i_cntr = 0;
    float i_error = .0f;
    for (const auto& p : *this)
    {
      if (i != p.isotopeIndex)
      {
        if (i_cntr > 0)
        {
          sum += std::abs(i_error / i_cntr);
          ++count;
        }
        i = p.isotopeIndex;
        i_cntr = 0;
        i_error = 0;
      }
      i_cntr++;
      i_error += ppm ? getPPMError_(p) : getDaError_(p);
    }
    if (i_cntr > 0)
    {
      sum += std::abs(i_error / i_cntr);
      ++count;
    }
    return sum / (float)count;
  }

  float PeakGroup::getPPMError_(const LogMzPeak& p) const
  {
    auto mass = (float)(monoisotopic_mass_ + p.isotopeIndex * iso_da_distance_);
    return (float)((mass / (float)p.abs_charge + FLASHHelperClasses::getChargeMass(p.is_positive) - p.mz) / p.mz * 1e6);
  }

  float PeakGroup::getDaError_(const LogMzPeak& p) const
  {
    auto mass = (float)(monoisotopic_mass_ + p.isotopeIndex * iso_da_distance_);
    return (float)((mass - p.getUnchargedMass()));
  }

  int PeakGroup::getMinNegativeIsotopeIndex() const
  {
    return min_negative_isotope_index_;
  }

  void PeakGroup::updatePerChargeCos_(const FLASHHelperClasses::PrecalculatedAveragine& avg, double tol)
  {
    const auto& iso_dist = avg.getReference_(monoisotopic_mass_);
    static thread_local std::vector<float> current_per_isotope_intensities; // per-worker scratch, refilled by getPerIsotopeIntensities_
    current_per_isotope_intensities.assign(getIsotopeIntensities().size(), .0f);

    if (min_abs_charge_ == max_abs_charge_) setChargeIsotopeCosine(min_abs_charge_, getIsotopeCosine());
    else
    {
      for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_; abs_charge++)
      {
        int min_isotope_index, max_isotope_index;
        getPerIsotopeIntensities_(current_per_isotope_intensities, min_isotope_index, max_isotope_index, abs_charge, 0, tol);
        float cos_score
          = SpectralDeconvolution::getCosine(current_per_isotope_intensities, min_isotope_index, max_isotope_index, iso_dist, 0,
                                             SpectralDeconvolution::min_iso_size);
        setChargeIsotopeCosine(abs_charge, cos_score); //
      }
    }
  }

  int PeakGroup::updateQscore(const std::vector<LogMzPeak>& noisy_peaks,
                              const FLASHHelperClasses::PrecalculatedAveragine& avg,
                              const double min_cos,
                              const double tol,
                              const bool is_low_charge,
                              const std::vector<double>& excluded_masses,
                              const bool is_last)
  {
    qscore_ = 0;

    if (empty()) { return 0; }
    // In the last iteration, the noise power over all charges is needed only for the SNR, after the charge fit and cosine checks below that
    // most groups (e.g., harmonic masses) fail. It is calculated there, from the signal peaks and charge range of now (updateChargeRange_ changes them).
    const bool defer_total_noise = is_last && min_abs_charge_ > 0;
    const int charge_count = max_abs_charge_ - min_abs_charge_ + 1;
    static thread_local std::vector<LogMzPeak> signal_peaks; // per-worker scratch
    if (defer_total_noise) { signal_peaks.assign(logMzpeaks_.begin(), logMzpeaks_.end()); }
    updatePerChargeInformation_(noisy_peaks, tol, is_last);
    if (is_last && ! defer_total_noise)
    {
      per_charge_noise_pwr_[0] = getNoisePeakPower_(noisy_peaks.data(), noisy_peaks.data() + noisy_peaks.size(), logMzpeaks_.data(),
                                                    logMzpeaks_.data() + logMzpeaks_.size(), charge_count, 0, tol);
    }
    updateChargeRange_();
    updateChargeFitScoreAndChargeIntensities_(is_low_charge);
    if (charge_score_ < .7f) //
    {
      return 0;
    }
    updateMonoMassAndIsotopeIntensities(tol); //
    if (per_isotope_int_.empty() || max_abs_charge_ < min_abs_charge_) { return 0; }
    int h_offset;
    int window_width = is_last ? 0 : -1;

    isotope_cosine_score_ = SpectralDeconvolution::getIsotopeCosineAndIsoOffset(
      monoisotopic_mass_, per_isotope_int_, h_offset, avg,
      -min_negative_isotope_index_, // change if to select cosine calculation and if to get second best hits
      window_width, excluded_masses);
    if (!is_last && h_offset != 0) return h_offset;

    if (isotope_cosine_score_ < min_cos) { return 0; }
    if (defer_total_noise)
    {
      per_charge_noise_pwr_[0] = getNoisePeakPower_(noisy_peaks.data(), noisy_peaks.data() + noisy_peaks.size(), signal_peaks.data(),
                                                    signal_peaks.data() + signal_peaks.size(), charge_count, 0, tol);
    }
    updatePerChargeCos_(avg, tol);
    updateAvgMassError_();
    updateSNR_((float)avg.getSNRMultiplicationFactor(monoisotopic_mass_));
    //if (target_decoy_type_ == PeakGroup::signal_decoy) isotope_cosine_score_ = std::min(1.0f, isotope_cosine_score_ + .02f); // the smaller

    for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_; abs_charge++)
    {
      if (getChargeSNR(abs_charge) > getChargeSNR(max_snr_abs_charge_)) { max_snr_abs_charge_ = abs_charge; }
    }

    qscores_ = PeakGroupScoring::getQscores(this);
    qscore_ = PeakGroupScoring::getQscore(this);
    return h_offset;
  }

  float PeakGroup::getNoisePeakPower_(const LogMzPeak* noisy_begin,
                                      const LogMzPeak* noisy_end,
                                      const LogMzPeak* signal_begin,
                                      const LogMzPeak* signal_end,
                                      const int charge_count,
                                      const int z,
                                      const double tol) const
  {
    if (noisy_begin == noisy_end) return 0;
    const Size max_noisy_peak_number = z == 0? std::min(200, 40 * charge_count) :40; // too many noise peaks will slow down the process
    constexpr Size bin_number_margin = 8;
    constexpr Size max_bin_number = bin_number_margin + 12; // 12 bin + 8 extra bin
    float threshold = -1;
    static thread_local std::vector<std::pair<Peak1D, bool>> all_peaks; // peak + is signal? (reused across calls)
    all_peaks.clear();
    all_peaks.reserve(max_noisy_peak_number + (signal_end - signal_begin));

    // the ranges hold the peaks of charge z only (all charges for z == 0)
    const Size noise_peak_count = noisy_end - noisy_begin;
    // get intensity threshold: the (max_noisy_peak_number + 1)-th largest intensity, the element a full descending sort puts there
    if (noise_peak_count > max_noisy_peak_number)
    {
      static thread_local std::vector<float> intensities;
      intensities.clear();
      intensities.reserve(noise_peak_count);
      for (const LogMzPeak* p = noisy_begin; p != noisy_end; ++p) intensities.push_back(p->intensity);
      std::nth_element(intensities.begin(), intensities.begin() + max_noisy_peak_number, intensities.end(), std::greater<float>());
      threshold = intensities[max_noisy_peak_number];
    }

    // filter peaks and check which mzs are signal and which are noise.
    for (const LogMzPeak* p = noisy_begin; p != noisy_end; ++p)
    {
      if (p->intensity < threshold) continue;
      all_peaks.emplace_back(Peak1D(p->getUnchargedMass(), p->intensity), false);
    }

    if (all_peaks.empty()) return 0;

    for (const LogMzPeak* p = signal_begin; p != signal_end; ++p)
    {
      if (p->intensity < threshold) continue;
      all_peaks.emplace_back(Peak1D(p->getUnchargedMass(), p->intensity), true);
    }

    std::sort(all_peaks.begin(), all_peaks.end(),
              [](std::pair<Peak1D, bool>& left, std::pair<Peak1D, bool>& right) { return left.first.getMZ() < right.first.getMZ(); }); //

    float charge_noise_pwr = 0;
    const Size peak_count = all_peaks.size();

    // per-worker tables reused across calls: masses and intensities of the sorted peaks, and per peak (row) the edge of each bin (0 = none),
    // the bins in which it is the target of an edge, and the bins in which it heads a chain (has an edge but is no target)
    static thread_local std::vector<double> peak_masses;
    static thread_local std::vector<float> peak_intensities;
    static thread_local std::vector<unsigned> edge_table;
    static thread_local std::vector<std::uint32_t> target_bins;
    static thread_local std::vector<std::uint32_t> head_bins;
    peak_masses.resize(peak_count);
    peak_intensities.resize(peak_count);
    edge_table.resize(peak_count * max_bin_number);
    target_bins.assign(peak_count, 0);
    head_bins.resize(peak_count);
    for (Size i = 0; i < peak_count; i++)
    {
      peak_masses[i] = all_peaks[i].first.getMZ();
      peak_intensities[i] = all_peaks[i].first.getIntensity();
    }

    // (sum, bin) pairs; equal sums keep the last bin like a std::map would, visited by descending sum
    std::array<std::pair<float, Size>, max_bin_number> max_intensity_sum_to_bin;
    Size max_intensity_sum_count = 0;
    static constexpr std::array<double, 3> div_factors {1.0, 2.0, 3.0}; // allow two skips for each bin
    std::uint32_t used_bins = 0; // bins with at least one edge
    // first collect all possible edges. An edge means mass difference between two peaks.
    // largest normalized distance that can still produce an edge; beyond this even the largest divisor lands outside the bin range
    const double max_edge_normalized_dist = ((double)max_bin_number - 0.5) * div_factors.back() / (double)(max_bin_number - bin_number_margin);
    std::array<double, max_bin_number> per_bin_error {}; // of the edge of peak i in a bin, valid for the bins in edge_bins
    std::array<unsigned, max_bin_number> per_bin_edge {};
    for (Size i = 0; i < peak_count; i++)
    {
      const bool p1_signal = all_peaks[i].second;
      const double p1_mass = peak_masses[i];
      std::uint32_t edge_bins = 0;
      for (Size j = i + 1; j < peak_count; j++)
      {
        const double normalized_dist = (peak_masses[j] - p1_mass) / iso_da_distance_;
        if (normalized_dist >= max_edge_normalized_dist) break; // peaks sorted by mass: no further j can form an edge

        if (p1_signal && all_peaks[j].second
            && normalized_dist >= .75) // if both are signals, and they are different from each other by more than .75 isotope distance, do not connect.
                                       // Otherwise, connect as they may a part of consecutive other noisy peaks.
        {
          continue;
        }

        for (double d : div_factors)
        {
          double distance = normalized_dist / d * (max_bin_number - bin_number_margin);
          // == (Size)round(distance) for the non-negative distances of the mass-sorted peaks, without the libm call
          const Size bin_floor = (Size)distance;
          Size bin = bin_floor + (distance - (double)bin_floor >= 0.5 ? 1 : 0);
          if (bin == 0) { continue; }
          if (bin >= max_bin_number) { continue; } // larger divisors give smaller bins that may still be in range, so skip (not break) this d

          used_bins |= std::uint32_t(1) << bin;
          double current_error = d * max_bin_number + std::abs((double)bin - distance); // larger when d gets larger. For the same d, comparable.
          if ((edge_bins >> bin & 1) && per_bin_error[bin] < current_error) { continue; }
          per_bin_edge[bin] = (unsigned)j;
          per_bin_error[bin] = current_error;
          edge_bins |= std::uint32_t(1) << bin;
        }
      }
      unsigned* edges = edge_table.data() + i * max_bin_number;
      std::fill(edges, edges + max_bin_number, 0u);
      for (std::uint32_t b = edge_bins; b != 0; b &= b - 1)
      {
        const int bin = std::countr_zero(b);
        edges[bin] = per_bin_edge[bin];
        target_bins[per_bin_edge[bin]] |= std::uint32_t(1) << bin;
      }
      // edges point to later peaks only, so all edges into peak i are known here
      head_bins[i] = edge_bins & ~target_bins[i];
    }

    // then from each bin find the highest intensity path consisting of the same mass differences.
    // each chain is walked once from its head, in peak order per bin; the chain sum is maximized at the head (intensities are non-negative)
    std::array<float, max_bin_number> max_sum_intensity {};
    std::array<int, max_bin_number> per_bin_start_index; // -2 means bin is empty. -1 means bin is used. zero or positive = edge index
    for (Size k = 0; k < max_bin_number; k++)
    {
      per_bin_start_index[k] = (used_bins >> k & 1) ? -1 : -2;
    }
    for (Size i = 0; i < peak_count; i++)
    {
      for (std::uint32_t b = head_bins[i]; b != 0; b &= b - 1)
      {
        const int bin = std::countr_zero(b);
        float sum_intensity = peak_intensities[i];

        for (Size j = edge_table[i * max_bin_number + bin]; j != 0; j = edge_table[j * max_bin_number + bin])
        {
          sum_intensity += peak_intensities[j];
        }

        if (max_sum_intensity[bin] < sum_intensity) // at least two edges should be there.
        {
          max_sum_intensity[bin] = sum_intensity;
          per_bin_start_index[bin] = (int)i; //
        }
      }
    }
    for (Size k = 0; k < max_bin_number; k++)
    {
      if (per_bin_start_index[k] == -2) { continue; }
      // how to deal with profile peaks?
      bool replaced = false;
      for (Size m = 0; m < max_intensity_sum_count; ++m)
      {
        if (max_intensity_sum_to_bin[m].first == max_sum_intensity[k])
        {
          max_intensity_sum_to_bin[m].second = k;
          replaced = true;
          break;
        }
      }
      if (! replaced) { max_intensity_sum_to_bin[max_intensity_sum_count++] = {max_sum_intensity[k], k}; }
    }

    static thread_local std::vector<char> unused; // reused across calls
    unused.assign(peak_count, 1);
    std::sort(max_intensity_sum_to_bin.begin(), max_intensity_sum_to_bin.begin() + max_intensity_sum_count,
              [](const auto& l, const auto& r) { return l.first > r.first; });

    // Now from the highest intensity path to the lowest, sum up intensities excluding already used peaks or signal peaks.
    for (auto it = max_intensity_sum_to_bin.begin(); it != max_intensity_sum_to_bin.begin() + max_intensity_sum_count; ++it)
    {
      //if (signal_pwr / 100 > it->first * it->first) continue; // if the noise summed intensity is too small, skip
      Size bin = it->second;
      int index = per_bin_start_index[bin];
      if (index < 0) { continue; }

      const auto edge = [&](Size p) -> Size { return edge_table[p * max_bin_number + bin]; };
      const double ori_mass = peak_masses[index];
      const int ori_index = index;
      float sum_intensity = .0;

      while (index < (int)peak_count)
      {
        const double p_mass = peak_masses[index];
        if (p_mass - ori_mass > tol / 2.0 * p_mass) break;
        float intensity = // p_signal ? 0 :
          peak_intensities[index];

        Size j = edge(index);
        if (j == 0) { break; }

        if (unused[index])
        {
          sum_intensity += intensity;
          unused[index] = false;
        }
        else { break; }

        for (; j != 0; j = edge(j))
        {
          if (unused[j])
          {
            sum_intensity += peak_intensities[j];
            unused[j] = false;
          }
          else
            break;
        }
        index++;
      }

      index = ori_index - 1;
      while (index >= 0)
      {
        const double p_mass = peak_masses[index];
        if (ori_mass - p_mass > tol / 2.0 * ori_mass) break;
        float intensity = // p_signal ? 0 :
          peak_intensities[index];

        Size j = edge(index);
        if (j == 0) { break; }

        if (unused[index])
        {
          sum_intensity += intensity;
          unused[index] = false;
        }
        else { break; }

        for (; j != 0; j = edge(j))
        {
          if (unused[j])
          {
            sum_intensity += peak_intensities[j];
            unused[j] = false;
          }
          else
            break;
        }
        index--;
      }

      charge_noise_pwr += sum_intensity * sum_intensity;
    }

    // if still peaks are remaining, add their individual power.
    double prev_mass = .0;
    float tmp_int_sum = 0;
    for (Size index = 0; index < unused.size(); ++index)
    {
      if (! unused[index]) { continue; }
      if (! all_peaks[index].second)
      {
        const double p_mass = peak_masses[index];
        if (p_mass - prev_mass > p_mass * tol)
        {
          charge_noise_pwr += tmp_int_sum * tmp_int_sum;
          prev_mass = p_mass;
          tmp_int_sum = 0;
        }
        tmp_int_sum += peak_intensities[index];
      }
    }
    charge_noise_pwr += tmp_int_sum * tmp_int_sum;

    return charge_noise_pwr;
  }


  void PeakGroup::updatePerChargeInformation_(const std::vector<LogMzPeak>& noisy_peaks, const double tol, const bool is_last)
  {
    // assign reuses the capacity of the previous call (same sizes and values as new vectors)
    per_charge_sum_signal_squared_.assign(1 + max_abs_charge_, .0f);
    per_charge_int_.assign(1 + max_abs_charge_, .0f);
    int max_iso = 0;

    // calculate per charge intensity, and per charge sum of signal intensity squared
    for (const auto& p : logMzpeaks_)
    {
      per_charge_int_[p.abs_charge] += p.intensity;
      max_iso = std::max(max_iso, p.isotopeIndex);
    }
    // per-worker scratch replaces a fresh (charge x isotope) matrix per call; row-major, same summation order
    static thread_local std::vector<float> per_charge_isotope_int;
    static thread_local std::vector<std::uint64_t> touched_cells; // one bit per matrix cell that a peak was added to
    const Size iso_cols = 1 + (Size)max_iso;
    const Size cell_count = (1 + (Size)max_abs_charge_) * iso_cols;
    per_charge_isotope_int.assign(cell_count, .0f);
    touched_cells.assign((cell_count + 63) / 64, 0);
    for (const auto& p : logMzpeaks_)
    {
      const Size cell = p.abs_charge * iso_cols + p.isotopeIndex;
      per_charge_isotope_int[cell] += p.intensity;
      touched_cells[cell / 64] |= std::uint64_t(1) << (cell % 64);
    }

    // The other cells are +0, and adding their squares leaves the sums unchanged: only the touched cells are summed, in the same order
    for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
    {
      const Size row_begin = z * iso_cols;
      const Size row_end = row_begin + iso_cols;
      for (Size w = row_begin / 64; w * 64 < row_end; w++)
      {
        std::uint64_t bits = touched_cells[w];
        if (w * 64 < row_begin) { bits &= ~std::uint64_t(0) << (row_begin - w * 64); }
        if (row_end - w * 64 < 64) { bits &= (std::uint64_t(1) << (row_end - w * 64)) - 1; }
        for (; bits != 0; bits &= bits - 1)
        {
          float v = per_charge_isotope_int[w * 64 + std::countr_zero(bits)];
          per_charge_sum_signal_squared_[0] += v * v;
          per_charge_sum_signal_squared_[z] += v * v;
        }
      }
    }

    // for each charge calculate signal and noise power
    per_charge_noise_pwr_.assign(1 + max_abs_charge_, .0f);

    if (is_last) // the noise power over all charges (index 0) is set by updateQscore
    {
      const int charge_count = max_abs_charge_ - min_abs_charge_ + 1;
      // the noisy and the signal peaks grouped by charge once (stable, so each group keeps the original order) instead of rescanning both
      // lists for every charge
      static thread_local std::vector<LogMzPeak> noisy_by_charge, signal_by_charge; // per-worker scratch
      static thread_local std::vector<Size> noisy_offset, signal_offset, fill;
      const auto group_by_charge = [&](const std::vector<LogMzPeak>& peaks, std::vector<LogMzPeak>& grouped, std::vector<Size>& offset) {
        offset.assign(max_abs_charge_ + 2, 0);
        for (const auto& p : peaks)
        {
          if (p.abs_charge >= min_abs_charge_ && p.abs_charge <= max_abs_charge_) offset[p.abs_charge + 1]++;
        }
        for (int z = 1; z <= max_abs_charge_ + 1; z++) offset[z] += offset[z - 1];
        grouped.resize(offset[max_abs_charge_ + 1]);
        fill.assign(offset.begin(), offset.end());
        for (const auto& p : peaks)
        {
          if (p.abs_charge >= min_abs_charge_ && p.abs_charge <= max_abs_charge_) grouped[fill[p.abs_charge]++] = p;
        }
      };
      group_by_charge(noisy_peaks, noisy_by_charge, noisy_offset);
      group_by_charge(logMzpeaks_, signal_by_charge, signal_offset);
      for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
      {
        per_charge_noise_pwr_[z] = getNoisePeakPower_(noisy_by_charge.data() + noisy_offset[z], noisy_by_charge.data() + noisy_offset[z + 1],
                                                      signal_by_charge.data() + signal_offset[z], signal_by_charge.data() + signal_offset[z + 1],
                                                      charge_count, z, tol);
      }
    }
    else
    {
      for (const auto& p : noisy_peaks)
      {
        float pwr = p.intensity * p.intensity;
        if (p.abs_charge < min_abs_charge_ || p.abs_charge > max_abs_charge_) continue;
        per_charge_noise_pwr_[0] += pwr;
        per_charge_noise_pwr_[p.abs_charge] += pwr;
      }
    }
  }

  void PeakGroup::updateChargeRange_()
  {
    int max_sig_charge = 0;
    double max_sig = 0;

    // Noise floor: a charge without a detected noise peak does not have zero noise, its noise is below the detection limit of the spectrum.
    // Each isotope position examined for noise at a charge (recruitAllPeaksInSpectrum: isotope 0 up to the highest signal isotope) contributes
    // at most the power of a peak at the detection limit. Without the floor, a charge with no unmatched peak in its isotope window scored
    // signal^2 in raw intensity units (1e10 - 1e14) and became the anchor with one to three coincidental peaks, so the range collapsed around
    // it. The floor replaces the unit-less "1 +" and makes the per-charge SNR dimensionless (independent of the intensity scale).
    int max_isotope_index = 0;
    for (const auto& p : logMzpeaks_)
    {
      max_isotope_index = std::max(max_isotope_index, p.isotopeIndex);
    }
    const double noise_floor_pwr = (double)(max_isotope_index + 1) * detection_limit_ * detection_limit_;
    auto get_charge_snr = [&](const int z) -> double {
      if (noise_floor_pwr <= 0) { return per_charge_int_[z] * per_charge_int_[z] / (1 + per_charge_noise_pwr_[z]); } // no detection limit set
      return (double)per_charge_int_[z] * per_charge_int_[z] / (per_charge_noise_pwr_[z] + noise_floor_pwr);
    };

    // first, find the maximum snr charge.
    for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
    {
      double tmp_snr = get_charge_snr(z);
      if (max_sig < tmp_snr)
      {
        max_sig = tmp_snr;
        max_sig_charge = z;
      }
    }

    // determine the final charge ranges based on per charge power.
    // If more than two consecutive charges do not contain any signal peak, the charge range stops at that charge.

    int new_max_abs_charge;
    int new_min_abs_charge;

    new_max_abs_charge = new_min_abs_charge = max_sig_charge;
    double threshold = std::min(max_sig / 10, 1.0);
    for (int z = max_sig_charge; z <= max_abs_charge_; z++)
    {
      if (get_charge_snr(z) < threshold) { break; }
      new_max_abs_charge = z;
    }

    for (int z = max_sig_charge; z >= min_abs_charge_; z--)
    {
      if (get_charge_snr(z) < threshold) { break; }
      new_min_abs_charge = z;
    }

    // if the updated charge range is different from the original one, signal and noisy peaks are again updated
    if (max_abs_charge_ != new_max_abs_charge || min_abs_charge_ != new_min_abs_charge)
    {
      // now only take the signal and noise peaks within the new charge range (in place, order kept).
      logMzpeaks_.erase(std::remove_if(logMzpeaks_.begin(), logMzpeaks_.end(),
                                       [&](const LogMzPeak& p) { return p.abs_charge < new_min_abs_charge || p.abs_charge > new_max_abs_charge; }),
                        logMzpeaks_.end());
      max_abs_charge_ = new_max_abs_charge;
      min_abs_charge_ = new_min_abs_charge;
    }
    if (min_abs_charge_ > max_abs_charge_) { clear_(); }
    else { sort(); }
  }

  std::vector<FLASHHelperClasses::LogMzPeak> PeakGroup::recruitAllPeaksInSpectrum(const MSSpectrum& spec,
                                                                                  const double tol,
                                                                                  const FLASHHelperClasses::PrecalculatedAveragine& avg,
                                                                                  double mono_mass,
                                                                                  bool renew_signal_peaks)
  {
    std::vector<LogMzPeak> noisy_peaks;
    recruitAllPeaksInSpectrum(spec, tol, avg, mono_mass, noisy_peaks, renew_signal_peaks);
    return noisy_peaks;
  }

  void PeakGroup::recruitAllPeaksInSpectrum(const MSSpectrum& spec,
                                            const double tol,
                                            const FLASHHelperClasses::PrecalculatedAveragine& avg,
                                            double mono_mass,
                                            std::vector<LogMzPeak>& noisy_peaks,
                                            bool renew_signal_peaks,
                                            const std::vector<double>* log_mzs)
  {
    const double mul_tol = .8; // not all peaks within tolerance are considered as signal peaks.
    if (log_mzs != nullptr && log_mzs->size() != spec.size()) { log_mzs = nullptr; }
    // the log m/z of a peak: precomputed for the spectrum, or computed here
    const auto make_peak = [&](Size index) {
      return log_mzs != nullptr ? LogMzPeak(spec[index], is_positive_, (*log_mzs)[index]) : LogMzPeak(spec[index], is_positive_);
    };
    noisy_peaks.clear();
    // keep this group's own capacity instead of freeing it on every renewal
    if (renew_signal_peaks) logMzpeaks_.clear();

    if (mono_mass < 0) { return; }

    if (max_abs_charge_ - min_abs_charge_ < max_abs_charge_ / 20) // if charge range is too small ...
    {
      return;
    }

    monoisotopic_mass_ = mono_mass;

    // int iso_margin = 3; // how many isotopes do we want to scan before the monoisotopic mass?
    int max_isotope = (int)avg.getLastIndex(mono_mass);
    int min_isotope = (int)(avg.getApexIndex(mono_mass) - avg.getLeftCountFromApex(mono_mass) + min_negative_isotope_index_);
    int max_signal_isotope = 0;
    min_isotope = std::max(min_negative_isotope_index_, min_isotope);

    negative_iso_peaks_.clear();

    noisy_peaks.reserve(max_isotope * (max_abs_charge_ - min_abs_charge_ + 1) * 2);

    // scan from the largest to the smallest charges and recruit the raw peaks for this monoisotopic_mass
    for (int c = max_abs_charge_; c >= min_abs_charge_; c--)
    {
      if (c <= 0) { break; }
      double cmz = (mono_mass) / c + FLASHHelperClasses::getChargeMass(is_positive_);
      double left_mz = (mono_mass - (1 - min_negative_isotope_index_) * iso_da_distance_) / c + FLASHHelperClasses::getChargeMass(is_positive_);
      Size index = spec.findNearest(left_mz * (1 - tol * mul_tol));
      double iso_delta = iso_da_distance_ / c;

      for (; index < spec.size(); index++)
      {
        float pint = spec[index].getIntensity();
        if (pint <= 0) { continue; }
        double pmz = spec[index].getMZ();
        int iso_index = (int)round((pmz - cmz) / iso_delta);
        if (iso_index > max_isotope) { break; }

        if (iso_index < min_isotope) { continue; }

        if (std::abs(pmz - cmz - iso_index * iso_delta) <= pmz * tol * mul_tol)
        {
          auto p = make_peak(index);
          p.isotopeIndex = iso_index;
          p.abs_charge = c;
          if (iso_index < 0) { negative_iso_peaks_.push_back(p); }
          else
          {
            max_signal_isotope = std::max(max_signal_isotope, iso_index);
            if (renew_signal_peaks) logMzpeaks_.push_back(p);
          }
        }
        else if (iso_index >= 0)
        {
          auto p = make_peak(index);
          int noise_iso_index = (int)floor((pmz - cmz) / iso_delta);

          p.isotopeIndex = noise_iso_index;
          p.abs_charge = c;
          noisy_peaks.push_back(p);
        }
      }

      if (index >= spec.size()) { break; }
    }

    noisy_peaks.erase(std::remove_if(noisy_peaks.begin(), noisy_peaks.end(), [max_signal_isotope](const LogMzPeak& p) { return p.isotopeIndex > max_signal_isotope; }),
                      noisy_peaks.end());
  }

  void PeakGroup::updateChargeFitScoreAndChargeIntensities_(bool is_low_charge)
  {
    if (max_abs_charge_ == min_abs_charge_)
    {
      charge_score_ = 1;
      return;
    }
    float max_per_charge_intensity = .0;
    float summed_intensity = .0;
    int max_index = -1;
    int first_index = -1;
    int last_index = -1;

    for (int c = min_abs_charge_; c <= max_abs_charge_; c++)
    {
      summed_intensity += per_charge_int_[c];
      if (per_charge_int_[c] > 0)
      {
        if (first_index < 0) { first_index = c; }
        last_index = c;
      }

      if (max_per_charge_intensity > per_charge_int_[c]) { continue; }
      max_per_charge_intensity = per_charge_int_[c];
      max_index = c;
    }
    if (max_index < 0)
    {
      charge_score_ = 0;
      return;
    }

    first_index = first_index < 0 ? 0 : first_index;
    float p = .0f;
    const float factor = .3;
    for (int c = max_index; c < last_index; c++)
    {
      float diff = per_charge_int_[c + 1] - per_charge_int_[c];
      if (diff > 0) p += diff;
      else if (! is_low_charge && diff < -(per_charge_int_[c]) * factor)
        p -= diff + (per_charge_int_[c]) * factor;
    }
    for (int c = max_index; c > first_index; c--)
    {
      float diff = per_charge_int_[c - 1] - per_charge_int_[c];
      if (diff > 0) p += diff;
      else if (! is_low_charge && diff < -(per_charge_int_[c]) * factor)
        p -= diff + (per_charge_int_[c]) * factor;
    }
    charge_score_ = std::max(.0f, 1.0f - p / summed_intensity);
  }

  void PeakGroup::getPerIsotopeIntensities_(std::vector<float>& intensities,
                                            int& min_isotope_index,
                                            int& max_isotope_index,
                                            int abs_charge,
                                            int min_negative_isotope_index,
                                            double tol)
  {
    int da_tol = (int)(tol * logMzpeaks_[0].getUnchargedMass() / 2.0);

    min_isotope_index = 1e5;
    max_isotope_index = -1; // this is inclusive!!

    for (const auto& p : logMzpeaks_)
    {
      max_isotope_index = max_isotope_index < p.isotopeIndex ? p.isotopeIndex : max_isotope_index;
      min_isotope_index = min_isotope_index < p.isotopeIndex ? min_isotope_index : p.isotopeIndex;
    }
    intensities.assign(max_isotope_index + 1 + da_tol - min_negative_isotope_index_, .0f);

    double denom = 2.0 * std::pow(da_tol / 2.355, 2.0) / iso_da_distance_;

    for (const auto& peak : logMzpeaks_)
    {
      if (abs_charge > 0 && peak.abs_charge != abs_charge) { continue; }

      if (peak.isotopeIndex - min_negative_isotope_index >= (int)intensities.size()) { continue; }

      if (peak.isotopeIndex < 0) { continue; }

      for (int margin = -da_tol; margin <= da_tol; margin++)
      {
        int index = peak.isotopeIndex + margin - min_negative_isotope_index;
        if (index < 0 || index >= (int)intensities.size()) continue;
        intensities[index] += denom > 0 ? float(peak.intensity * exp(-margin * margin / denom)) : peak.intensity;
      }
    }

    if (min_negative_isotope_index != 0)
    {
      for (const auto& peak : negative_iso_peaks_)
      {
        if (abs_charge > 0 && peak.abs_charge != abs_charge) { continue; }

        if (peak.isotopeIndex - min_negative_isotope_index >= (int)intensities.size()) { continue; }

        if (peak.isotopeIndex - min_negative_isotope_index < 0) { continue; }

        for (int margin = -da_tol; margin <= da_tol; margin++)
        {
          int index = peak.isotopeIndex + margin - min_negative_isotope_index;
          if (index < 0 || index >= (int)intensities.size()) continue;
          intensities[index] += denom > 0 ? float(peak.intensity * exp(-margin * margin / denom)) : peak.intensity;
        }
      }
    }
  }

  void PeakGroup::updateMonoMassAndIsotopeIntensities(double tol)
  {
    if (logMzpeaks_.empty()) { return; }
    int min_isotope_index, max_isotope_index;
    sortPeaks_();

    intensity_ = .0;
    double nominator = .0;

    getPerIsotopeIntensities_(per_isotope_int_, min_isotope_index, max_isotope_index, 0, min_negative_isotope_index_, tol);

    for (const auto& p : logMzpeaks_)
    {
      float pi = p.intensity;
      if (p.isotopeIndex - min_negative_isotope_index_ < 0) { continue; }
      nominator += pi * (p.getUnchargedMass() - p.isotopeIndex * iso_da_distance_);
      intensity_ += pi;
    }

    monoisotopic_mass_ = nominator / intensity_;
  }

  void PeakGroup::setSNR(const float snr)
  {
    snr_ = snr;
  }

  void PeakGroup::setChargeSNR(const int abs_charge, const float c_snr)
  {
    if (max_abs_charge_ < abs_charge) { return; }
    if (per_charge_snr_.empty()) { per_charge_snr_ = std::vector<float>(1 + max_abs_charge_, .0); }
    per_charge_snr_[abs_charge] = c_snr;
  }

  void PeakGroup::setTargeted()
  {
    is_targeted_ = true;
  }

  void PeakGroup::setChargeIsotopeCosine(const int abs_charge, const float cos)
  {
    if (max_abs_charge_ < abs_charge) { return; }
    if (per_charge_cos_.empty()) { per_charge_cos_ = std::vector<float>(1 + max_abs_charge_, .0); }
    per_charge_cos_[abs_charge] = cos;
  }

  void PeakGroup::setAbsChargeRange(const int min_abs_charge, const int max_abs_charge)
  {
    min_abs_charge_ = min_abs_charge;
    max_abs_charge_ = max_abs_charge;
  }

  void PeakGroup::setScanNumber(const int sn)
  {
    scan_number_ = sn;
  }

  void PeakGroup::setIsotopeCosine(const float cos)
  {
    isotope_cosine_score_ = cos;
  }

  void PeakGroup::setMonoisotopicMass(double mono_mass)
  {
    monoisotopic_mass_ = mono_mass;
  }

  void PeakGroup::setRepAbsCharge(const int max_snr_abs_charge)
  {
    max_snr_abs_charge_ = max_snr_abs_charge;
  }

  void PeakGroup::setChargeScore(const float score)
  {
    charge_score_ = score;
  }

  void PeakGroup::setAvgPPMError(const float error)
  {
    avg_ppm_error_ = error;
  }

  void PeakGroup::setQscore(const double qscore)
  {
    qscore_ = std::min(1.0, qscore);
  }

  std::tuple<double, double> PeakGroup::getRepMzRange() const
  {
    return getMzRange(getRepAbsCharge());
  }

  std::tuple<double, double> PeakGroup::getMzRange(int abs_charge) const
  {
    double mz_start = -1;
    double mz_end = -10;
    if (! (abs_charge > max_abs_charge_ || abs_charge < min_abs_charge_))
    {
      for (const auto& tmp_p : logMzpeaks_)
      {
        if (tmp_p.abs_charge != abs_charge) { continue; }
        if (mz_start < 0) { mz_start = tmp_p.mz; }
        else { mz_start = mz_start < tmp_p.mz ? mz_start : tmp_p.mz; }
        mz_end = mz_end > tmp_p.mz ? mz_end : tmp_p.mz;
      }
    }
    return std::tuple<double, double> {mz_start, mz_end};
  }

  std::tuple<int, int> PeakGroup::getAbsChargeRange() const
  {
    return std::tuple<int, int> {min_abs_charge_, max_abs_charge_};
  }

  const std::vector<float>& PeakGroup::getIsotopeIntensities() const
  {
    return per_isotope_int_;
  }

  int PeakGroup::getScanNumber() const
  {
    return scan_number_;
  }

  double PeakGroup::getMonoMass() const
  {
    return monoisotopic_mass_;
  }

  float PeakGroup::getIntensity() const
  {
    return intensity_;
  }

  float PeakGroup::getMaxChargeIntensity() const
  {
    float max_int = 0.0f;
    for (int c = min_abs_charge_; c <= max_abs_charge_; c++)
    {
      if (c >= 0 && c < (int)per_charge_int_.size() && per_charge_int_[c] > max_int)
      {
        max_int = per_charge_int_[c];
      }
    }
    return max_int;
  }

  float PeakGroup::getIsotopeCosine() const
  {
    return isotope_cosine_score_;
  }

  int PeakGroup::getRepAbsCharge() const
  {
    return max_snr_abs_charge_;
  }

  int PeakGroup::getMaxIntensityAbsCharge() const
  {
    int max_int_charge = min_abs_charge_;
    float max_intensity = 0.0f;
    for (int c = min_abs_charge_; c <= max_abs_charge_; ++c)
    {
      float intensity = getChargeIntensity(c);
      if (intensity > max_intensity)
      {
        max_intensity = intensity;
        max_int_charge = c;
      }
    }
    return max_int_charge;
  }

  double PeakGroup::getQscore() const
  {
    return qscore_;
  }

  double PeakGroup::getQscore2D() const
  {
    return std::max(qscore_, qscore2D_);
  }

  std::unordered_map<int, float> PeakGroup::getAllQscores() const
  {
    return qscores_;
  }

  int PeakGroup::getBestQScoreCharge() const
  {
    if (qscores_.empty())
    {
      return -1;
    }

    float max_score = 0.0f;
    std::vector<int> best_charges;

    // Find maximum score and collect all charges with that score
    for (const auto& [charge, score] : qscores_)
    {
      if (score > max_score)
      {
        max_score = score;
        best_charges.clear();
        best_charges.push_back(charge);
      }
      else if (score == max_score)
      {
        best_charges.push_back(charge);
      }
    }

    if (best_charges.empty())
    {
      return -1;
    }

    // Handle ties by selecting charge closest to representative charge
    if (best_charges.size() == 1)
    {
      return best_charges[0];
    }

    int rep_charge = getRepAbsCharge();
    if (rep_charge <= 0)
    {
      return best_charges[0]; // Fallback if no representative charge available
    }

    int best_charge = best_charges[0];
    int min_distance = std::abs(best_charges[0] - rep_charge);

    for (int charge : best_charges)
    {
      int distance = std::abs(charge - rep_charge);
      if (distance < min_distance)
      {
        min_distance = distance;
        best_charge = charge;
      }
    }

    return best_charge;
  }

  float PeakGroup::getBestQScore() const
  {
    int best_charge = getBestQScoreCharge();
    if (best_charge == -1)
    {
      return 0.0f;
    }

    auto it = qscores_.find(best_charge);
    return (it != qscores_.end()) ? it->second : 0.0f;
  }

  void PeakGroup::setFeatureIndex(uint findex)
  {
    findex_ = findex;
  }

  bool PeakGroup::isTargeted() const
  {
    return is_targeted_;
  }

  void PeakGroup::updateSNR_(float mul_factor)
  {
    per_charge_snr_.assign(1 + max_abs_charge_, .0f);
    float total_nom = 0;
    float total_denom = 0;
    for (size_t c = min_abs_charge_; (int)c < 1 + max_abs_charge_; ++c)
    {
      if (per_charge_cos_.size() > c)
      {
        float per_charge_cos_squared = per_charge_cos_[c] * per_charge_cos_[c];
        float sig_pwr = per_charge_sum_signal_squared_[c] * per_charge_cos_squared;
        float nom = 1e-6f * mul_factor + mul_factor * sig_pwr;
        float denom = 1e-6f * mul_factor + per_charge_noise_pwr_[c] + (1 - per_charge_cos_squared) * per_charge_sum_signal_squared_[c];

        per_charge_snr_[c] = denom <= 0 ? .0f : (nom / denom);

        total_denom += (1 - per_charge_cos_squared) * per_charge_sum_signal_squared_[c];
        total_nom += sqrt(sig_pwr);
      }
    }

    snr_ = mul_factor * (1e-6f + total_nom * total_nom) / (1e-6f * mul_factor + total_denom + per_charge_noise_pwr_[0]);

    per_charge_sum_signal_squared_.clear();
    per_charge_noise_pwr_.clear();
  }

  float PeakGroup::getQvalue() const
  {
    return qvalue_;
  }


  float PeakGroup::getSNR() const
  {
    return snr_;
  }

  float PeakGroup::getChargeScore() const
  {
    return charge_score_;
  }

  float PeakGroup::getAvgPPMError() const
  {
    return avg_ppm_error_;
  }

  float PeakGroup::getAvgDaError() const
  {
    return avg_da_error_;
  }

  float PeakGroup::getChargeSNR(const int abs_charge) const
  {
    if (abs_charge < 0 || (int)per_charge_snr_.size() <= abs_charge) { return 0; }
    return per_charge_snr_[abs_charge];
  }

  float PeakGroup::getChargeIsotopeCosine(const int abs_charge) const
  {
    if (abs_charge < 0 || (int)per_charge_cos_.size() <= abs_charge) { return 0; }
    return per_charge_cos_[abs_charge];
  }

  float PeakGroup::getChargeIntensity(const int abs_charge) const
  {
    if (abs_charge < 0 || per_charge_int_.empty() || (int)per_charge_int_.size() <= abs_charge) { return 0; }
    return per_charge_int_[abs_charge];
  }

  bool PeakGroup::isPositive() const
  {
    return is_positive_;
  }

  PeakGroup::TargetDecoyType PeakGroup::getTargetDecoyType() const
  {
    return target_decoy_type_;
  }

  void PeakGroup::setTargetDecoyType(PeakGroup::TargetDecoyType index)
  {
    target_decoy_type_ = index;
  }

  void PeakGroup::setIsotopeDaDistance(const double d)
  {
    iso_da_distance_ = d;
  }

  double PeakGroup::getIsotopeDaDistance() const
  {
    return iso_da_distance_;
  }

  void PeakGroup::setDetectionLimit(const float intensity)
  {
    detection_limit_ = intensity;
  }

  void PeakGroup::setWindowBounds(const double max_mass, const int max_charge)
  {
    window_max_mass_ = max_mass;
    window_max_charge_ = max_charge;
  }

  double PeakGroup::getWindowMaxMass() const
  {
    return window_max_mass_;
  }

  int PeakGroup::getWindowMaxCharge() const
  {
    return window_max_charge_;
  }

  void PeakGroup::setUseWindowBounds(const bool use)
  {
    use_window_bounds_ = use;
  }

  bool PeakGroup::useWindowBounds() const
  {
    return use_window_bounds_;
  }

  void PeakGroup::setIndex(const uint i)
  {
    index_ = i;
  }

  uint PeakGroup::getIndex() const
  {
    return index_;
  }

  uint PeakGroup::getFeatureIndex() const
  {
    return findex_;
  }

  void PeakGroup::setQscore2D(double fqscore)
  {
    qscore2D_ = fqscore;
  }

  std::vector<FLASHHelperClasses::LogMzPeak>::const_iterator PeakGroup::begin() const noexcept
  {
    return logMzpeaks_.begin();
  }

  std::vector<FLASHHelperClasses::LogMzPeak>::const_iterator PeakGroup::end() const noexcept
  {
    return logMzpeaks_.end();
  }

  std::vector<FLASHHelperClasses::LogMzPeak>::iterator PeakGroup::begin() noexcept
  {
    return logMzpeaks_.begin();
  }

  std::vector<FLASHHelperClasses::LogMzPeak>::iterator PeakGroup::end() noexcept
  {
    return logMzpeaks_.end();
  }

  const FLASHHelperClasses::LogMzPeak& PeakGroup::operator[](const Size i) const
  {
    return logMzpeaks_[i];
  }

  void PeakGroup::push_back(const FLASHHelperClasses::LogMzPeak& pg)
  {
    logMzpeaks_.push_back(pg);
  }

  FLASHHelperClasses::LogMzPeak& PeakGroup::back()
  {
    return logMzpeaks_.back();
  }

  Size PeakGroup::size() const noexcept
  {
    return logMzpeaks_.size();
  }

  void PeakGroup::clear_()
  {
    std::vector<LogMzPeak>().swap(logMzpeaks_);
  }

  void PeakGroup::reserve(Size n)
  {
    logMzpeaks_.reserve(n);
  }

  bool PeakGroup::empty() const
  {
    return logMzpeaks_.empty();
  }

  void PeakGroup::swap(std::vector<FLASHHelperClasses::LogMzPeak>& x)
  {
    logMzpeaks_.swap(x);
  }

  void PeakGroup::sort()
  {
    sortPeaks_();
  }

  void PeakGroup::sortPeaks_()
  {
    // A strictly increasing range has exactly one sorted order, so skipping std::sort gives the
    // identical sequence. Ranges with ties (equal log m/z and intensity) are still sorted as before.
    const auto first_not_strict = std::adjacent_find(logMzpeaks_.begin(), logMzpeaks_.end(),
                                                     [](const LogMzPeak& a, const LogMzPeak& b) { return ! (a < b); });
    if (first_not_strict != logMzpeaks_.end()) { std::sort(logMzpeaks_.begin(), logMzpeaks_.end()); }
  }

  void PeakGroup::setQvalue(double q)
  {
    qvalue_ = (float)q;
  }

} // namespace OpenMS
