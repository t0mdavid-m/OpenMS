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
#include <array>
#include <bit>
#include <cstdint>

namespace OpenMS
{
  namespace
  {
    // per-worker scratch of PeakGroup::getNoisePeakPower_ (one thread_local, so one TLS lookup per call)
    struct NoisePowerScratch
    {
      std::vector<std::pair<Peak1D, bool>> all_peaks;
      std::vector<float> intensities;
      std::vector<double> peak_masses;
      std::vector<float> peak_intensities;
      std::vector<unsigned> edge_table;
      std::vector<double> chain_bounds;
      std::vector<std::uint32_t> edge_bins;
      std::vector<std::uint32_t> target_bins;
      std::vector<std::uint32_t> head_bins;
      std::vector<char> unused;      // per peak: not yet in a chain sum
    };
  } // namespace

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

  // The per-charge noise powers of the last scoring iteration, computed when first read (chargeNoisePower_) from the same noisy and signal
  // peaks as before. A lower bound of each: the square of the most intense noise peak of the charge. getNoisePeakPower_ always takes that
  // peak (it is above the intensity threshold), adds its intensity to one of the non-negative sums it squares, and adds the non-negative
  // squares: with finite non-negative intensities every float step is monotone, so the noise power is at least the bound.
  struct PeakGroup::ChargeNoise_
  {
    const std::vector<FLASHHelperClasses::LogMzPeak>* noisy_peaks = nullptr;
    const std::vector<FLASHHelperClasses::LogMzPeak>* signal_peaks = nullptr;
    int charge_count = 0;
    double tol = 0;
    int min_charge = 1, max_charge = 0;
    std::vector<float> lower_bounds; // per charge (index 0 and charges without noise peaks: 0)
    std::vector<char> computed;      // per charge

    /// set up for the charges [@p min_z, @p max_z]; false (nothing left to compute later) unless all intensities are finite, non-negative and
    /// sum to less than 1e30 (no NaN, no overflow of the float sums)
    bool prepare(const std::vector<FLASHHelperClasses::LogMzPeak>& noisy, const std::vector<FLASHHelperClasses::LogMzPeak>& signal, const int count,
                 const double t, const int min_z, const int max_z)
    {
      if (min_z < 1 || max_z < min_z) { return false; }
      double total = 0;
      for (const auto& p : signal)
      {
        if (! (p.intensity >= 0) || ! std::isfinite(p.intensity)) { return false; }
        total += p.intensity;
      }
      lower_bounds.assign(max_z + 1, .0f);
      for (const auto& p : noisy)
      {
        if (! (p.intensity >= 0) || ! std::isfinite(p.intensity)) { return false; }
        total += p.intensity;
        if (p.abs_charge >= min_z && p.abs_charge <= max_z) { lower_bounds[p.abs_charge] = std::max(lower_bounds[p.abs_charge], p.intensity); }
      }
      if (! (total < 1e30)) { return false; }
      for (int z = min_z; z <= max_z; z++)
      {
        lower_bounds[z] = lower_bounds[z] * lower_bounds[z];
      }
      computed.assign(max_z + 1, 0);
      noisy_peaks = &noisy;
      signal_peaks = &signal;
      charge_count = count;
      tol = t;
      min_charge = min_z;
      max_charge = max_z;
      return true;
    }
  };

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
    // The per-charge noise powers of the last iteration are also computed when first read, from those signal peaks: updateChargeRange_ reads
    // only the charges that can have the largest SNR or decide the new charge range, updateSNR_ those in the new range
    static thread_local ChargeNoise_ charge_noise; // per-worker scratch
    ChargeNoise_* const lazy_noise
      = defer_total_noise && charge_noise.prepare(noisy_peaks, signal_peaks, charge_count, tol, min_abs_charge_, max_abs_charge_) ? &charge_noise : nullptr;
    updatePerChargeInformation_(noisy_peaks, tol, is_last, lazy_noise);
    if (is_last && ! defer_total_noise) { per_charge_noise_pwr_[0] = getNoisePeakPower_(noisy_peaks, logMzpeaks_, charge_count, 0, tol); }
    updateChargeRange_(lazy_noise);
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
    if (defer_total_noise) { per_charge_noise_pwr_[0] = getNoisePeakPower_(noisy_peaks, signal_peaks, charge_count, 0, tol); }
    updatePerChargeCos_(avg, tol);
    updateAvgMassError_();
    updateSNR_((float)avg.getSNRMultiplicationFactor(monoisotopic_mass_), lazy_noise);
    //if (target_decoy_type_ == PeakGroup::signal_decoy) isotope_cosine_score_ = std::min(1.0f, isotope_cosine_score_ + .02f); // the smaller

    for (int abs_charge = min_abs_charge_; abs_charge <= max_abs_charge_; abs_charge++)
    {
      if (getChargeSNR(abs_charge) > getChargeSNR(max_snr_abs_charge_)) { max_snr_abs_charge_ = abs_charge; }
    }

    qscores_ = PeakGroupScoring::getQscores(this);
    qscore_ = PeakGroupScoring::getQscore(this);
    return h_offset;
  }

  float PeakGroup::getNoisePeakPower_(const std::vector<FLASHHelperClasses::LogMzPeak>& noisy_peaks,
                                      const std::vector<FLASHHelperClasses::LogMzPeak>& signal_peaks,
                                      const int charge_count,
                                      const int z,
                                      const double tol) const
  {
    if (noisy_peaks.empty()) return 0;
    const Size max_noisy_peak_number = z == 0? std::min(200, 40 * charge_count) :40; // too many noise peaks will slow down the process
    constexpr Size bin_number_margin = 8;
    constexpr Size max_bin_number = bin_number_margin + 12; // 12 bin + 8 extra bin
    float threshold = -1;
    static thread_local NoisePowerScratch scratch; // one per worker, reused across calls
    auto& all_peaks = scratch.all_peaks;           // peak + is signal?
    all_peaks.clear();
    all_peaks.reserve(max_noisy_peak_number + signal_peaks.size());

    auto noise_peak_count = std::count_if(noisy_peaks.begin(), noisy_peaks.end(), [&](const auto& noisy_peak) {
      return z == 0 || noisy_peak.abs_charge == z;
    });

    if (noise_peak_count == 0) return 0;
    // get intensity threshold
    if (noise_peak_count > (int)max_noisy_peak_number)
    {
      auto& intensities = scratch.intensities;
      intensities.clear();
      intensities.reserve(noise_peak_count);
      for (const auto& noisy_peak : noisy_peaks)
      {
        if (z > 0 && noisy_peak.abs_charge != z) continue;
        //if (noisy_peak.abs_charge < min_abs_charge_ || noisy_peak.abs_charge > max_abs_charge_) continue;
        intensities.push_back(noisy_peak.intensity);
      }

      std::sort(intensities.rbegin(), intensities.rend());
      threshold = intensities[max_noisy_peak_number];
    }

    // filter peaks and check which mzs are signal and which are noise.
    for (const auto& noisy_peak : noisy_peaks)
    {
      if ((z > 0 && noisy_peak.abs_charge != z) || noisy_peak.intensity < threshold) continue;
      //if (noisy_peak.abs_charge < min_abs_charge_ || noisy_peak.abs_charge > max_abs_charge_) continue;
      all_peaks.emplace_back(Peak1D(noisy_peak.getUnchargedMass(), noisy_peak.intensity), false);
    }

    if (all_peaks.empty()) return 0;

    for (const auto& peak : signal_peaks)
    {
      if ((z > 0 && peak.abs_charge != z) || peak.intensity < threshold) continue;
      //if (peak.abs_charge < min_abs_charge_ || peak.abs_charge > max_abs_charge_) continue;
      all_peaks.emplace_back(Peak1D(peak.getUnchargedMass(), peak.intensity), true);
    }

    std::sort(all_peaks.begin(), all_peaks.end(),
              [](std::pair<Peak1D, bool>& left, std::pair<Peak1D, bool>& right) { return left.first.getMZ() < right.first.getMZ(); }); //

    float charge_noise_pwr = 0;
    const Size peak_count = all_peaks.size();

    // per-worker tables (grown, never shrunk): masses and intensities of the sorted peaks, and per peak (row) the edge of each bin (0 = none),
    // the bins in which it has an edge, in which it is the target of an edge, and in which it heads a chain (has an edge but is no target)
    if (scratch.peak_masses.size() < peak_count)
    {
      scratch.peak_masses.resize(peak_count);
      scratch.peak_intensities.resize(peak_count);
      scratch.edge_table.resize(peak_count * max_bin_number);
      scratch.chain_bounds.resize(peak_count * max_bin_number);
      scratch.edge_bins.resize(peak_count);
      scratch.target_bins.resize(peak_count);
      scratch.head_bins.resize(peak_count);
      scratch.unused.resize(peak_count);
    }
    double* const peak_masses = scratch.peak_masses.data();
    float* const peak_intensities = scratch.peak_intensities.data();
    unsigned* const edge_table = scratch.edge_table.data();
    std::uint32_t* const edge_bins_of = scratch.edge_bins.data();
    std::uint32_t* const target_bins = scratch.target_bins.data();
    std::uint32_t* const head_bins = scratch.head_bins.data();
    std::fill_n(target_bins, peak_count, 0u);
    // bounded_chains: finite non-negative intensities whose float sums cannot overflow (for the chain sum bounds below)
    bool bounded_chains = peak_count < (Size(1) << 20);
    double intensity_total = 0;
    for (Size i = 0; i < peak_count; i++)
    {
      const float intensity = all_peaks[i].first.getIntensity();
      peak_masses[i] = all_peaks[i].first.getMZ();
      peak_intensities[i] = intensity;
      bounded_chains = bounded_chains && intensity >= 0 && std::isfinite(intensity);
      intensity_total += intensity;
    }
    bounded_chains = bounded_chains && intensity_total < 1e30;

    // (sum, bin) pairs; equal sums keep the last bin like a std::map would, visited by descending sum
    std::array<std::pair<float, Size>, max_bin_number> max_intensity_sum_to_bin;
    Size max_intensity_sum_count = 0;
    static constexpr std::array<double, 3> div_factors {1.0, 2.0, 3.0}; // allow two skips for each bin
    std::uint32_t used_bins = 0; // bins with at least one edge
    // first collect all possible edges. An edge means mass difference between two peaks.
    // largest normalized distance that can still produce an edge; beyond this even the largest divisor lands outside the bin range
    const double max_edge_normalized_dist = ((double)max_bin_number - 0.5) * div_factors.back() / (double)(max_bin_number - bin_number_margin);
    Size edge_count = 0, merged_edge_count = 0; // edges, and those into a peak that is already the target of an edge of the bin
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
          Size bin = (Size)round(distance);
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
      unsigned* edges = edge_table + i * max_bin_number;
      std::fill(edges, edges + max_bin_number, 0u);
      for (std::uint32_t b = edge_bins; b != 0; b &= b - 1)
      {
        const int bin = std::countr_zero(b);
        edges[bin] = per_bin_edge[bin];
        merged_edge_count += target_bins[per_bin_edge[bin]] >> bin & 1;
        target_bins[per_bin_edge[bin]] |= std::uint32_t(1) << bin;
        edge_count++;
      }
      edge_bins_of[i] = edge_bins;
      // edges point to later peaks only, so all edges into peak i are known here
      head_bins[i] = edge_bins & ~target_bins[i];
    }

    // then from each bin find the highest intensity path consisting of the same mass differences: per bin the largest chain sum (walked from
    // its head, float) and the first head in peak order that has it, as when every chain is walked in peak order with a strict "<" update
    std::array<float, max_bin_number> max_sum_intensity {};
    std::array<int, max_bin_number> per_bin_start_index; // -2 means bin is empty. -1 means bin is used. zero or positive = edge index
    for (Size k = 0; k < max_bin_number; k++)
    {
      per_bin_start_index[k] = (used_bins >> k & 1) ? -1 : -2;
    }
    const auto chain_sum = [&](Size head, int bin) {
      float sum_intensity = peak_intensities[head];
      for (Size j = edge_table[head * max_bin_number + bin]; j != 0; j = edge_table[j * max_bin_number + bin])
      {
        sum_intensity += peak_intensities[j];
      }
      return sum_intensity;
    };
    // (the bounds pay off when chains merge, so that heads share their walks; otherwise, or for a few peaks, the plain walks are cheaper)
    if (bounded_chains && peak_count >= 32 && merged_edge_count * 4 >= edge_count)
    {
      // Upper bounds of the chain sums, so that most chains need no walk: the chain sum from each peak in double, in reverse peak order (edges
      // point to later peaks), times a factor covering the rounding of a float sum of up to peak_count non-negative terms (each step at most
      // 1 + 2^-24) and of the double sums (1 - 2^-53 per step): bound >= float chain sum.
      double* const chain_bounds = scratch.chain_bounds.data();
      for (Size i = peak_count; i-- > 0;)
      {
        const unsigned* edges = edge_table + i * max_bin_number;
        for (std::uint32_t b = edge_bins_of[i]; b != 0; b &= b - 1)
        {
          const int bin = std::countr_zero(b);
          const Size j = edges[bin];
          const double rest = (edge_bins_of[j] >> bin & 1) ? chain_bounds[j * max_bin_number + bin] : (double)peak_intensities[j];
          chain_bounds[i * max_bin_number + bin] = (double)peak_intensities[i] + rest;
        }
      }
      const double bound_factor = 1.0 + (double)(peak_count + 1) * 0x1p-22;
      // per bin, the head with the largest bound is walked first; a head whose bound is below the best sum so far cannot reach it
      std::array<double, max_bin_number> best_bound {};
      std::array<int, max_bin_number> best_head;
      best_head.fill(-1);
      for (Size i = 0; i < peak_count; i++)
      {
        for (std::uint32_t b = head_bins[i]; b != 0; b &= b - 1)
        {
          const int bin = std::countr_zero(b);
          const double bound = chain_bounds[i * max_bin_number + bin];
          if (best_head[bin] < 0 || best_bound[bin] < bound)
          {
            best_bound[bin] = bound;
            best_head[bin] = (int)i;
          }
        }
      }
      std::array<float, max_bin_number> best_sum {};
      for (std::uint32_t b = used_bins; b != 0; b &= b - 1)
      {
        const int bin = std::countr_zero(b);
        if (best_head[bin] >= 0) { best_sum[bin] = chain_sum(best_head[bin], bin); }
      }
      for (Size i = 0; i < peak_count; i++)
      {
        for (std::uint32_t b = head_bins[i]; b != 0; b &= b - 1)
        {
          const int bin = std::countr_zero(b);
          if ((int)i == best_head[bin] || chain_bounds[i * max_bin_number + bin] * bound_factor < best_sum[bin]) { continue; }
          const float sum_intensity = chain_sum(i, bin);
          if (best_sum[bin] < sum_intensity || (best_sum[bin] == sum_intensity && (int)i < best_head[bin]))
          {
            best_sum[bin] = sum_intensity;
            best_head[bin] = (int)i;
          }
        }
      }
      for (std::uint32_t b = used_bins; b != 0; b &= b - 1)
      {
        const int bin = std::countr_zero(b);
        if (best_head[bin] >= 0 && max_sum_intensity[bin] < best_sum[bin]) // at least two edges should be there.
        {
          max_sum_intensity[bin] = best_sum[bin];
          per_bin_start_index[bin] = best_head[bin];
        }
      }
    }
    else
    {
      // each chain is walked once from its head, in peak order per bin
      for (Size i = 0; i < peak_count; i++)
      {
        for (std::uint32_t b = head_bins[i]; b != 0; b &= b - 1)
        {
          const int bin = std::countr_zero(b);
          const float sum_intensity = chain_sum(i, bin);
          if (max_sum_intensity[bin] < sum_intensity) // at least two edges should be there.
          {
            max_sum_intensity[bin] = sum_intensity;
            per_bin_start_index[bin] = (int)i; //
          }
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

    char* const unused = scratch.unused.data();
    std::fill_n(unused, peak_count, 1);
    // by descending sum. The sums are distinct (equal ones were merged above) and never NaN (a maximum starting at 0 with "<" updates), so
    // every sort gives this order: insertion sort of the at most max_bin_number pairs
    for (Size m = 1; m < max_intensity_sum_count; m++)
    {
      const auto pair = max_intensity_sum_to_bin[m];
      Size k = m;
      for (; k > 0 && max_intensity_sum_to_bin[k - 1].first < pair.first; k--)
      {
        max_intensity_sum_to_bin[k] = max_intensity_sum_to_bin[k - 1];
      }
      max_intensity_sum_to_bin[k] = pair;
    }

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
    for (Size index = 0; index < peak_count; ++index)
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


  void PeakGroup::updatePerChargeInformation_(const std::vector<LogMzPeak>& noisy_peaks, const double tol, const bool is_last, ChargeNoise_* lazy_noise)
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
    if (per_charge_isotope_int.size() < cell_count) { per_charge_isotope_int.resize(cell_count); }
    touched_cells.assign((cell_count + 63) / 64, 0);
    for (const auto& p : logMzpeaks_)
    {
      const Size cell = p.abs_charge * iso_cols + p.isotopeIndex;
      const std::uint64_t cell_bit = std::uint64_t(1) << (cell % 64);
      // only touched cells are read below: each is set to +0 when first touched instead of zeroing the whole matrix
      if (! (touched_cells[cell / 64] & cell_bit)) { per_charge_isotope_int[cell] = .0f; }
      per_charge_isotope_int[cell] += p.intensity;
      touched_cells[cell / 64] |= cell_bit;
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

    if (is_last && lazy_noise == nullptr) // the noise power over all charges (index 0) is set by updateQscore; lazy ones by chargeNoisePower_
    {
      const int charge_count = max_abs_charge_ - min_abs_charge_ + 1;
      for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
      {
        per_charge_noise_pwr_[z] = getNoisePeakPower_(noisy_peaks, logMzpeaks_, charge_count, z, tol);
      }
    }
    else if (! is_last)
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

  float PeakGroup::chargeNoisePower_(const int z, ChargeNoise_* lazy_noise)
  {
    if (lazy_noise != nullptr && z >= lazy_noise->min_charge && z <= lazy_noise->max_charge && ! lazy_noise->computed[z])
    {
      per_charge_noise_pwr_[z] = getNoisePeakPower_(*lazy_noise->noisy_peaks, *lazy_noise->signal_peaks, lazy_noise->charge_count, z, lazy_noise->tol);
      lazy_noise->computed[z] = 1;
    }
    return per_charge_noise_pwr_[z];
  }

  void PeakGroup::updateChargeRange_(ChargeNoise_* lazy_noise)
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
      const float noise_pwr = chargeNoisePower_(z, lazy_noise);
      if (noise_floor_pwr <= 0) { return per_charge_int_[z] * per_charge_int_[z] / (1 + noise_pwr); } // no detection limit set
      return (double)per_charge_int_[z] * per_charge_int_[z] / (noise_pwr + noise_floor_pwr);
    };
    // With lazy noise powers and a detection limit: an upper bound of the SNR of a charge, the same SNR with the lower bound of its noise power
    // (the SNR falls as the noise power grows; intensities are finite and non-negative). A charge whose bound is below what a decision needs
    // is decided without its noise power, as it would be with it.
    const bool bounded = lazy_noise != nullptr && noise_floor_pwr > 0;
    auto charge_snr_bound = [&](const int z) -> double {
      return (double)per_charge_int_[z] * per_charge_int_[z] / ((double)lazy_noise->lower_bounds[z] + noise_floor_pwr);
    };

    // first, find the maximum snr charge.
    if (bounded)
    {
      // the first charge with the largest SNR as in the loop below: the charge with the largest bound first, then the charges whose bound
      // reaches the best SNR so far (any other charge has a smaller SNR), equal SNRs keeping the lower charge
      int first_charge = -1;
      double first_bound = 0;
      for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
      {
        const double bound = charge_snr_bound(z);
        if (first_charge < 0 || first_bound < bound)
        {
          first_bound = bound;
          first_charge = z;
        }
      }
      if (first_charge >= 0)
      {
        double best_snr = get_charge_snr(first_charge);
        int best_charge = first_charge;
        for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
        {
          if (z == first_charge || charge_snr_bound(z) < best_snr) { continue; }
          const double tmp_snr = get_charge_snr(z);
          if (best_snr < tmp_snr || (best_snr == tmp_snr && z < best_charge))
          {
            best_snr = tmp_snr;
            best_charge = z;
          }
        }
        if (max_sig < best_snr)
        {
          max_sig = best_snr;
          max_sig_charge = best_charge;
        }
      }
    }
    else
    {
      for (int z = min_abs_charge_; z <= max_abs_charge_; z++)
      {
        double tmp_snr = get_charge_snr(z);
        if (max_sig < tmp_snr)
        {
          max_sig = tmp_snr;
          max_sig_charge = z;
        }
      }
    }

    // determine the final charge ranges based on per charge power.
    // If more than two consecutive charges do not contain any signal peak, the charge range stops at that charge.

    int new_max_abs_charge;
    int new_min_abs_charge;

    new_max_abs_charge = new_min_abs_charge = max_sig_charge;
    double threshold = std::min(max_sig / 10, 1.0);
    // a charge below the threshold is passed over (another species' peaks inside the isotope span at one charge make its SNR collapse);
    // the range ends at the last charge above it before the third weak charge in a row
    int weak = 0;
    for (int z = max_sig_charge; z <= max_abs_charge_; z++)
    {
      if ((bounded && charge_snr_bound(z) < threshold) || get_charge_snr(z) < threshold)
      {
        if (++weak > 2) { break; }
        continue;
      }
      weak = 0;
      new_max_abs_charge = z;
    }

    weak = 0;
    for (int z = max_sig_charge; z >= min_abs_charge_; z--)
    {
      if ((bounded && charge_snr_bound(z) < threshold) || get_charge_snr(z) < threshold)
      {
        if (++weak > 2) { break; }
        continue;
      }
      weak = 0;
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
                                            const std::vector<double>* log_mzs,
                                            const NearestPeakIndex* nearest_peaks)
  {
    const double mul_tol = .8; // not all peaks within tolerance are considered as signal peaks.
    if (log_mzs != nullptr && log_mzs->size() != spec.size()) { log_mzs = nullptr; }
    if (nearest_peaks != nullptr && ! nearest_peaks->usableFor(spec)) { nearest_peaks = nullptr; }
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

    // the signal peaks are noted (peak, isotope index, charge) in scan order and made LogMzPeaks at the end, in a vector grown once
    struct RecruitedPeak
    {
      Size index;
      int isotope_index;
      int abs_charge;
    };
    static thread_local std::vector<RecruitedPeak> signal_recruits; // per-worker scratch
    signal_recruits.clear();
    const Size peak_count = spec.size();

    // scan from the largest to the smallest charges and recruit the raw peaks for this monoisotopic_mass
    for (int c = max_abs_charge_; c >= min_abs_charge_; c--)
    {
      if (c <= 0) { break; }
      double cmz = (mono_mass) / c + FLASHHelperClasses::getChargeMass(is_positive_);
      double left_mz = (mono_mass - (1 - min_negative_isotope_index_) * iso_da_distance_) / c + FLASHHelperClasses::getChargeMass(is_positive_);
      const double first_mz = left_mz * (1 - tol * mul_tol);
      Size index = nearest_peaks != nullptr ? nearest_peaks->findNearest(spec, first_mz) : spec.findNearest(first_mz);
      double iso_delta = iso_da_distance_ / c;

      for (; index < peak_count; index++)
      {
        float pint = spec[index].getIntensity();
        if (pint <= 0) { continue; }
        double pmz = spec[index].getMZ();
        int iso_index = (int)round((pmz - cmz) / iso_delta);
        if (iso_index > max_isotope) { break; }

        if (iso_index < min_isotope) { continue; }

        if (std::abs(pmz - cmz - iso_index * iso_delta) <= pmz * tol * mul_tol)
        {
          if (iso_index < 0)
          {
            auto p = make_peak(index);
            p.isotopeIndex = iso_index;
            p.abs_charge = c;
            negative_iso_peaks_.push_back(p);
          }
          else
          {
            max_signal_isotope = std::max(max_signal_isotope, iso_index);
            if (renew_signal_peaks) signal_recruits.push_back({index, iso_index, c});
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

      if (index >= peak_count) { break; }
    }

    if (renew_signal_peaks)
    {
      logMzpeaks_.reserve(signal_recruits.size());
      for (const auto& r : signal_recruits)
      {
        auto p = make_peak(r.index);
        p.isotopeIndex = r.isotope_index;
        p.abs_charge = r.abs_charge;
        logMzpeaks_.push_back(p);
      }
    }
    noisy_peaks.erase(std::remove_if(noisy_peaks.begin(), noisy_peaks.end(), [max_signal_isotope](const LogMzPeak& p) { return p.isotopeIndex > max_signal_isotope; }),
                      noisy_peaks.end());
  }

  void PeakGroup::NearestPeakIndex::build(const MSSpectrum& spec)
  {
    size_ = spec.size();
    usable_ = size_ > 0;
    for (Size i = 0; usable_ && i < size_; i++)
    {
      const double mz = spec[i].getMZ();
      usable_ = std::isfinite(mz) && (i == 0 || spec[i - 1].getMZ() <= mz);
    }
    if (! usable_)
    {
      first_peaks_.clear();
      return;
    }
    first_mz_ = spec[0].getMZ();
    last_mz_ = spec[size_ - 1].getMZ();
    const Size cell_count = 2 * size_;
    const double range = last_mz_ - first_mz_;
    cells_per_mz_ = range > 0 ? (double)cell_count / range : 0;
    if (! std::isfinite(cells_per_mz_)) { cells_per_mz_ = 0; } // only the start of the search depends on the grid
    first_peaks_.resize(cell_count + 1);
    Size i = 0;
    for (Size k = 0; k <= cell_count; k++)
    {
      const double cell_start = cells_per_mz_ > 0 ? first_mz_ + (double)k / cells_per_mz_ : first_mz_;
      while (i < size_ && spec[i].getMZ() < cell_start) { i++; }
      first_peaks_[k] = i;
    }
  }

  bool PeakGroup::NearestPeakIndex::usableFor(const MSSpectrum& spec) const
  {
    return usable_ && spec.size() == size_;
  }

  Size PeakGroup::NearestPeakIndex::findNearest(const MSSpectrum& spec, const double mz) const
  {
    const Size n = size_;
    // the lower bound: the first peak whose m/z is not below mz (MSSpectrum::MZBegin), unique as the m/z values are sorted
    Size lower;
    if (! (mz > first_mz_)) { lower = 0; } // also for NaN, as MZBegin
    else if (mz > last_mz_) { lower = n; }
    else
    {
      // from the first peak of the cell of mz, step to the lower bound; a binary search on the remaining side after a few steps
      const double cell = (mz - first_mz_) * cells_per_mz_;
      const Size last_cell = first_peaks_.size() - 1;
      lower = first_peaks_[cell < (double)last_cell ? (Size)cell : last_cell];
      const int max_steps = 8;
      if (lower < n && spec[lower].getMZ() < mz)
      {
        for (int step = 0; lower < n && spec[lower].getMZ() < mz; step++, lower++)
        {
          if (step == max_steps)
          {
            lower = spec.MZBegin(spec.begin() + lower, mz, spec.end()) - spec.begin();
            break;
          }
        }
      }
      else
      {
        for (int step = 0; lower > 0 && ! (spec[lower - 1].getMZ() < mz); step++, lower--)
        {
          if (step == max_steps)
          {
            lower = spec.MZBegin(spec.begin(), mz, spec.begin() + lower) - spec.begin();
            break;
          }
        }
      }
    }
    // the nearer of the lower bound and the peak before it, as MSSpectrum::findNearest
    if (lower == 0) { return 0; }
    if (lower == n) { return n - 1; }
    if (std::fabs(spec[lower].getMZ() - mz) < std::fabs(spec[lower - 1].getMZ() - mz)) { return lower; }
    return lower - 1;
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

  void PeakGroup::updateSNR_(float mul_factor, ChargeNoise_* lazy_noise)
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
        float denom = 1e-6f * mul_factor + chargeNoisePower_((int)c, lazy_noise) + (1 - per_charge_cos_squared) * per_charge_sum_signal_squared_[c];

        // A charge with no signal peak has no SNR. Without this the two 1e-6 regularisers cancel and a charge that
        // updateChargeRange_ passed over (it skips up to two weak charges) scores exactly 1.0, above real charges,
        // and becomes the representative charge with no peak to isolate.
        per_charge_snr_[c] = (per_charge_sum_signal_squared_[c] <= 0 || denom <= 0) ? .0f : (nom / denom);

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
