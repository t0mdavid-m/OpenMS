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
// $Maintainer: Tom David Müller$
// $Authors: Kyowon Jeong, Tom David Müller $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/KERNEL/ConsensusMap.h> 
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroupScoring.h>
#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/SpectrumSettings.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <functional>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <unordered_set>
#include <nlohmann/json.hpp>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{

// Anonymous namespace for identification workflow helper structures and functions
namespace
{
  /// Get ion type mass shift for N-terminal (prefix) ions
  /// Uses Residue class methods for consistency with FLASHTnT
  inline double getPrefixIonShift(const String& ion_type)
  {
    if (ion_type == "b") return Residue::getInternalToBIon().getMonoWeight();
    if (ion_type == "a") return Residue::getInternalToAIon().getMonoWeight();
    if (ion_type == "c") return Residue::getInternalToCIon().getMonoWeight();
    return Residue::getInternalToBIon().getMonoWeight(); // default to b
  }

  /// Get ion type mass shift for C-terminal (suffix) ions
  /// Uses Residue class methods for consistency with FLASHTnT
  inline double getSuffixIonShift(const String& ion_type)
  {
    if (ion_type == "y") return Residue::getInternalToYIon().getMonoWeight();
    if (ion_type == "x") return Residue::getInternalToXIon().getMonoWeight();
    if (ion_type == "z") return Residue::getInternalToZIon().getMonoWeight();
    return Residue::getInternalToYIon().getMonoWeight(); // default to y
  }

  /// Map fragmentation method name to ion types (case-insensitive via lowercase)
  inline std::vector<std::string> getIonTypesForFragmentationMethod(const String& method)
  {
    String lower_method = method;
    std::transform(lower_method.begin(), lower_method.end(), lower_method.begin(), ::tolower);

    if (lower_method == "hcd") return {"b", "y"};
    if (lower_method == "cid") return {"b", "y"};
    if (lower_method == "etd") return {"c", "z"};
    if (lower_method == "ethcd") return {"b", "c", "y", "z"};
    if (lower_method == "etcid") return {"b", "c", "y", "z"};
    if (lower_method == "uvpd") return {"a", "b", "c", "x", "y", "z"};
    return {"b", "y"};  // default to HCD
  }

  /// Check if ion type is a prefix (N-terminal) ion
  inline bool isPrefixIon(char ion_type)
  {
    return ion_type == 'a' || ion_type == 'b' || ion_type == 'c';
  }

  // PTMSite is now defined as FLASHIda::PTMSite in the header file
  using PTMSite = FLASHIda::PTMSite;


  /// Calculate theoretical fragment masses with PTM adjustments from FLASHExtender
  /// For each PTM at positions [start, end], fragments containing that region get the PTM mass added
  void calculatePTMAdjustedFragmentMasses(
      const String& sequence,
      const std::vector<FLASHIda::PTMSite>& ptm_sites,
      double b_ion_shift,
      double y_ion_shift,
      std::vector<double>& prefix_masses,
      std::vector<double>& suffix_masses)
  {
    // Clean sequence (replace I with L for mass calculation)
    String clean_seq = sequence;
    std::replace(clean_seq.begin(), clean_seq.end(), 'I', 'L');

    int seq_len = static_cast<int>(clean_seq.size());
    prefix_masses.clear();
    suffix_masses.clear();
    prefix_masses.reserve(seq_len);
    suffix_masses.reserve(seq_len);

    // For each position, calculate cumulative PTM mass contribution
    // A b-ion at position i (b_i) contains residues 1..i
    // For ambiguous PTMs [S,E], use end_position as cutoff to be consistent with y-ions
    // This ensures b_E+ are modified (definitely contain PTM regardless of exact position)
    std::vector<double> ptm_at_or_before(seq_len + 1, 0.0);
    for (const auto& ptm : ptm_sites)
    {
      // PTM affects all b-ions from b(end_position) onwards (definitely contain PTM)
      for (int i = ptm.end_position; i <= seq_len; ++i)
      {
        ptm_at_or_before[i] += ptm.mass_shift;
      }
    }

    // Calculate prefix masses (b-ions) with PTM adjustment
    double cumulative_mass = 0.0;
    for (int i = 0; i < seq_len; ++i)
    {
      char aa = clean_seq[i];
      if (aa != 'X')
      {
        const Residue* res = ResidueDB::getInstance()->getResidue(aa);
        if (res != nullptr)
        {
          cumulative_mass += res->getMonoWeight(Residue::Internal);
        }
      }
      // b_{i+1} covers positions 1..(i+1), add PTM if start <= i+1
      double ptm_contribution = ptm_at_or_before[i + 1];
      prefix_masses.push_back(cumulative_mass + ptm_contribution + b_ion_shift);
    }

    // For y-ions: y_n covers positions (L-n+1)..L
    // y_n definitely contains PTM [S,E] if it covers position S (start_position)
    // That means (L-n+1) <= S, i.e., n >= L - S + 1
    // This is consistent with b-ion logic: only mark as modified if definitely contains PTM
    std::vector<double> ptm_for_suffix(seq_len + 1, 0.0);
    for (const auto& ptm : ptm_sites)
    {
      // y_n includes PTM if n >= L - start_position + 1 (definitely contains PTM)
      int min_n = seq_len - ptm.start_position + 1;
      for (int n = min_n; n <= seq_len; ++n)
      {
        ptm_for_suffix[n] += ptm.mass_shift;
      }
    }

    // Calculate suffix masses (y-ions) with PTM adjustment
    cumulative_mass = 0.0;
    for (int i = seq_len - 1; i >= 0; --i)
    {
      char aa = clean_seq[i];
      if (aa != 'X')
      {
        const Residue* res = ResidueDB::getInstance()->getResidue(aa);
        if (res != nullptr)
        {
          cumulative_mass += res->getMonoWeight(Residue::Internal);
        }
      }
      // This is y_{seq_len - i} (e.g., i=seq_len-1 gives y_1)
      int y_number = seq_len - i;
      double ptm_contribution = ptm_for_suffix[y_number];
      suffix_masses.push_back(cumulative_mass + ptm_contribution + y_ion_shift);
    }
  }


  /// Calculate theoretical fragment masses for multiple ion types with PTM adjustments
  /// @param sequence the protein sequence
  /// @param ptm_sites PTM sites from FLASHExtender
  /// @param ion_types vector of ion type strings (e.g., {"a", "b", "c", "x", "y", "z"})
  /// @param fragment_masses_map output: map from ion type char to vector of masses
  void calculatePTMAdjustedFragmentMassesMulti(
      const String& sequence,
      const std::vector<FLASHIda::PTMSite>& ptm_sites,
      const std::vector<std::string>& ion_types,
      std::map<char, std::vector<double>>& fragment_masses_map)
  {
    fragment_masses_map.clear();

    // Clean sequence (replace I with L for mass calculation)
    String clean_seq = sequence;
    std::replace(clean_seq.begin(), clean_seq.end(), 'I', 'L');
    int seq_len = static_cast<int>(clean_seq.size());

    // Separate prefix and suffix ion types
    std::vector<std::string> prefix_ions, suffix_ions;
    for (const auto& ion : ion_types)
    {
      if (ion == "a" || ion == "b" || ion == "c") prefix_ions.push_back(ion);
      else if (ion == "x" || ion == "y" || ion == "z") suffix_ions.push_back(ion);
    }

    // Pre-calculate PTM contributions for prefix ions (b-ion style)
    std::vector<double> ptm_at_or_before(seq_len + 1, 0.0);
    for (const auto& ptm : ptm_sites)
    {
      for (int i = ptm.end_position; i <= seq_len; ++i)
      {
        ptm_at_or_before[i] += ptm.mass_shift;
      }
    }

    // Pre-calculate PTM contributions for suffix ions (y-ion style)
    std::vector<double> ptm_for_suffix(seq_len + 1, 0.0);
    for (const auto& ptm : ptm_sites)
    {
      int min_n = seq_len - ptm.start_position + 1;
      for (int n = min_n; n <= seq_len; ++n)
      {
        ptm_for_suffix[n] += ptm.mass_shift;
      }
    }

    // Calculate base cumulative mass for prefix direction (N→C)
    std::vector<double> prefix_base_masses;
    prefix_base_masses.reserve(seq_len);
    double cumulative_mass = 0.0;
    for (int i = 0; i < seq_len; ++i)
    {
      char aa = clean_seq[i];
      if (aa != 'X')
      {
        const Residue* res = ResidueDB::getInstance()->getResidue(aa);
        if (res != nullptr)
        {
          cumulative_mass += res->getMonoWeight(Residue::Internal);
        }
      }
      double ptm_contribution = ptm_at_or_before[i + 1];
      prefix_base_masses.push_back(cumulative_mass + ptm_contribution);
    }

    // Calculate base cumulative mass for suffix direction (C→N)
    std::vector<double> suffix_base_masses;
    suffix_base_masses.reserve(seq_len);
    cumulative_mass = 0.0;
    for (int i = seq_len - 1; i >= 0; --i)
    {
      char aa = clean_seq[i];
      if (aa != 'X')
      {
        const Residue* res = ResidueDB::getInstance()->getResidue(aa);
        if (res != nullptr)
        {
          cumulative_mass += res->getMonoWeight(Residue::Internal);
        }
      }
      int y_number = seq_len - i;
      double ptm_contribution = ptm_for_suffix[y_number];
      suffix_base_masses.push_back(cumulative_mass + ptm_contribution);
    }

    // Apply ion-specific shifts for each prefix ion type
    for (const auto& ion : prefix_ions)
    {
      double shift = getPrefixIonShift(ion);
      char ion_char = ion[0];
      std::vector<double> masses;
      masses.reserve(seq_len);
      for (double base_mass : prefix_base_masses)
      {
        masses.push_back(base_mass + shift);
      }
      fragment_masses_map[ion_char] = std::move(masses);
    }

    // Apply ion-specific shifts for each suffix ion type
    for (const auto& ion : suffix_ions)
    {
      double shift = getSuffixIonShift(ion);
      char ion_char = ion[0];
      std::vector<double> masses;
      masses.reserve(seq_len);
      for (double base_mass : suffix_base_masses)
      {
        masses.push_back(base_mass + shift);
      }
      fragment_masses_map[ion_char] = std::move(masses);
    }
  }


} // anonymous namespace

/// optimal window margin
inline const double optimal_window_margin_ = .4;

/// constructor
FLASHIda::FLASHIda(char* arg)
{
  #ifdef _OPENMP
    omp_set_num_threads(4);
  #endif

    // Config must be JSON (legacy space-delimited format removed in Phase 8)
    std::string arg_str(arg);
    if (arg_str.empty() || arg_str[0] != '{')
    {
      throw std::invalid_argument("FLASHIda: config must be JSON (starts with '{')");
    }
    parseJSONConfig_(arg_str);
  }

  bool FLASHIda::isDifferentiallyAbundant(const double* mzs,
                              const double* ints,
                              const int length,
                              const double rt,
                              const int ms_level,
                              const char* name,
                              double reporter_mz_tol,
                              double fold_change_threshold,
                              bool only_one_condition
                              )
  {
    // Create spectrum
    auto spec = makeMSSpectrum_(mzs, ints, length, rt, ms_level, name);

    // Set precursor with HCD activation - neccessary for channel extractor
    OpenMS::Precursor precursor;
    precursor.setActivationMethods({OpenMS::Precursor::ActivationMethod::HCD});
    spec.setPrecursors({precursor});

    // Create experiment
    MSExperiment exp;
    exp.addSpectrum(spec);

    // TODO: Variable channel extractors
    TMTSixPlexQuantitationMethod quant_method;
    IsobaricChannelExtractor channel_extractor(&quant_method);

    // Set parameters
    Param p = channel_extractor.getParameters();
    p.setValue("reporter_mass_shift", reporter_mz_tol);
    channel_extractor.setParameters(p);

    // Extract reporter-ion channels into a ConsensusMap
    ConsensusMap consensus_map_raw;
    channel_extractor.extractChannels(exp, consensus_map_raw);

    // Collect m/z-intensity pairs from the ConsensusFeatures
    std::vector<std::pair<double, double>> mz_int;
    mz_int.reserve(quant_method.getNumberOfChannels());
    for (const auto& cf : consensus_map_raw)
    {
        for (auto& i : cf)
        {
            mz_int.emplace_back(i.getMZ(), i.getIntensity());
        }
    }

    // Something went wrong – bail out early
    if (mz_int.size() != quant_method.getNumberOfChannels())
    {
        std::cout << "FLASHIda - channel extraction failed..." << std::endl;
        return false;
    }

    // Sort by m/z to ensure channel order
    std::sort(mz_int.begin(), mz_int.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    // Extract intensities
    std::vector<double> intensities;
    intensities.reserve(mz_int.size());
    for (const auto& p2 : mz_int) intensities.push_back(p2.second);

    for (auto intensity : intensities) {
      std::cout << intensity << std::endl; 
    }

    // TODO: Make configurable
    if (only_one_condition) {
      bool first_sample_present = std::none_of(
        intensities.begin(), intensities.begin()+3, [](double x){ return x < 1e-3; }
      );
      std::cout << first_sample_present << std::endl;
      bool second_sample_present = std::none_of(
        intensities.begin()+3, intensities.end(), [](double x){ return x < 1e-3; }
      );
      std::cout << second_sample_present << std::endl;
      bool first_sample_missing = std::all_of(
        intensities.begin(), intensities.begin()+3, [](double x){ return x < 1e-3; }
      );
      std::cout << first_sample_missing << std::endl;
      bool second_sample_missing = std::all_of(
        intensities.begin()+3, intensities.end(), [](double x){ return x < 1e-3; }
      );
      std::cout << second_sample_missing << std::endl;
      // No signal
      if (first_sample_missing && second_sample_missing) {
        return false;
      }
      // Both are incomplete
      else if (!first_sample_present && !second_sample_present) {
        return false;
      }
      // One signal
      else if (
        (first_sample_missing || second_sample_missing)
        && (first_sample_present || second_sample_present)
      ) {
        return true;
      }
    }

    // Reject spectra with missing / too-low reporter peaks
    if (!only_one_condition && std::any_of(intensities.begin(), intensities.end(),
                    [](double x){ return x < 1e-3; }))
    {
        return false;
    }

    const double sample1_mean = std::accumulate(intensities.begin(),
                                                intensities.begin() + 3, 0.0) / 3.0;

    const double sample2_mean = std::accumulate(intensities.begin() + 3,
                                                intensities.end(), 0.0) / 3.0;

    const double fold_change = sample1_mean / sample2_mean;

    return (fold_change > fold_change_threshold) || ((1/fold_change) > fold_change_threshold);

  }

  int FLASHIda::getPeakGroups(const double* mzs,
                              const double* ints,
                              const int length,
                              const double rt,
                              const int ms_level,
                              const char* name,
                              const char* cv)
  {
    // int ret[2] = {0,0};
    auto spec = makeMSSpectrum_(mzs, ints, length, rt, ms_level, name);
    if (cv != nullptr) { spec.setMetaValue("filter string", DataValue("cv=" + std::string(cv))); }
    // selected_peak_groups_ = DeconvolvedSpectrum(spec, 1);
    if (ms_level == 1)
    {
      // current_max_mass_ = max_mass;
      // currentChargeRange = chargeRange;
    }
    else
    {
      return 0;
      // TODO precursor infor here
    }

    std::vector<DeconvolvedSpectrum> tmp;
    PeakGroup empty;

    target_masses_.clear();
    excluded_masses_.clear();
    if (targeting_mode_ == 1)  // Unified inclusion mode (merged mode 4 functionality)
    {
      active_targets_.clear();
      target_priority_map_.clear();

      if (!inclusion_targets_.empty())  // TSV-based targets
      {
        for (const auto& target : inclusion_targets_)
        {
          if (target.isActiveAt(rt))
          {
            active_targets_.push_back(&target);
            target_masses_.push_back(target.mass);

            // Build priority map for tie-breaking (use highest priority for each nominal mass)
            int nominal = SpectralDeconvolution::getNominalMass(target.mass);
            if (target_priority_map_.find(nominal) == target_priority_map_.end()
                || target.priority > target_priority_map_[nominal])
            {
              target_priority_map_[nominal] = target.priority;
            }
          }
        }
      }
      else  // Legacy .log/.out targets (existing behavior)
      {
        for (const auto& [mass, rts] : target_mass_rt_map_)
        {
          target_masses_.push_back(mass);
        }
      }

      std::sort(target_masses_.begin(), target_masses_.end());
      fd_.setTargetMasses(target_masses_, false);
    }
    else if (targeting_mode_ == 3)
    {
      for (const auto& [prt, masses] : exclusion_rt_masses_map_)
      {
        if (std::abs(rt - prt) >= rt_window_ && prt != 0) continue;
        for (double mass : masses)
        {
          excluded_masses_.push_back(mass);
        }
      }
      std::sort(excluded_masses_.begin(), excluded_masses_.end());
    }

    selected_peak_groups_.clear();
    deconvolved_spectrum_.clear();

    fd_.performSpectrumDeconvolution(spec, 0, empty);
    deconvolved_spectrum_ = fd_.getDeconvolvedSpectrum();
    // per spec deconvolution
    FLASHIda::filterPeakGroupsUsingMassExclusion_(ms_level, rt);
    // spec.clear(true);
    return (int)selected_peak_groups_.size();
  }

  void FLASHIda::filterPeakGroupsUsingMassExclusion_(const int ms_level, const double rt)
  {
    if (use_idscore_ && consider_all_Charge_states_ && hcd_energy_ < 0) {
      deconvolved_spectrum_.sortByIDScoreAllCharges();
    }
    else if (use_idscore_ && consider_all_Charge_states_) {
      deconvolved_spectrum_.sortByIDScoreAllCharges(hcd_energy_);
    }
    else if (use_idscore_ && !consider_all_Charge_states_ && hcd_energy_ < 0) {
      deconvolved_spectrum_.sortByIDScoreRepresentative();
    }
    else if (use_idscore_ && !consider_all_Charge_states_) {
      deconvolved_spectrum_.sortByIDScoreRepresentative(hcd_energy_);
    }
    else if (!use_idscore_ && consider_all_Charge_states_) {
      deconvolved_spectrum_.sortByQScoreAllCharges();
    }
    else {
      deconvolved_spectrum_.sortByQscore();
    }

    // Apply priority tie-breaking when TSV targets are loaded
    if (targeting_mode_ == 1 && !inclusion_targets_.empty())
    {
      std::stable_sort(deconvolved_spectrum_.begin(), deconvolved_spectrum_.end(),
        [this](const PeakGroup& a, const PeakGroup& b) {
          if (std::abs(a.getQscore() - b.getQscore()) < tie_threshold_)
          {
            int nom_a = SpectralDeconvolution::getNominalMass(a.getMonoMass());
            int nom_b = SpectralDeconvolution::getNominalMass(b.getMonoMass());
            int pri_a = target_priority_map_.count(nom_a) ? target_priority_map_.at(nom_a) : 0;
            int pri_b = target_priority_map_.count(nom_b) ? target_priority_map_.at(nom_b) : 0;
            return pri_a > pri_b;  // Higher priority first
          }
          return a.getQscore() > b.getQscore();
        });
    }

    Size mass_count = (Size)mass_count_[ms_level - 1];
    trigger_charges.clear();
    trigger_hcds.clear();
    trigger_scores.clear();
    trigger_charges.reserve(mass_count);
    trigger_hcds.reserve(mass_count);
    trigger_scores.reserve(mass_count);
    trigger_left_isolation_mzs_.clear();
    trigger_left_isolation_mzs_.reserve(mass_count);
    trigger_right_isolation_mzs_.clear();
    trigger_right_isolation_mzs_.reserve(mass_count);
    trigger_ids_.clear();
    trigger_ids_.reserve(mass_count);
    std::vector<int>* charges = nullptr;

    selected_peak_groups_.reserve(mass_count_.size());
    std::set<double> current_selected_mzs;    // current selected mzs
    std::set<double> current_selected_masses; // current selected mzs

    std::unordered_map<int, double> new_mz_rt_map_;
    std::unordered_map<int, double> new_mass_rt_map_;
    std::unordered_map<int, double> new_all_mass_rt_map_;
    std::unordered_map<int, double> new_mass_score_map_;
    std::unordered_map<int, double> t_mass_score_map_;

    // exclusion mode
    // TODO: Update IDScore bla, currently only qscore
    if (targeting_mode_ == 2)
    {
      for (const auto& [mass, rts] : target_mass_rt_map_)
      {
        int nominal_mass = SpectralDeconvolution::getNominalMass(mass);
        auto qscores = target_mass_qscore_map_[mass];
        for (uint i = 0; i < rts.size(); i++)
        {
          double prt = rts[i];
          double qscore = qscores[i];
          if (std::abs(rt - prt) < rt_window_)
          {
            auto inter = t_mass_score_map_.find(nominal_mass);
            if (inter == t_mass_score_map_.end()) { t_mass_score_map_[nominal_mass] = 1 - qscore; }
            else { t_mass_score_map_[nominal_mass] *= 1 - qscore; }
          }
        }
      }
    }

    // remove expired entries for tqscore_exceeding_mz_rt_map_
    for (const auto& [m, r] : tqscore_exceeding_mz_rt_map_)
    {
      if (rt - r > rt_window_) { continue; }
      new_mz_rt_map_[m] = r;
    }
    new_mz_rt_map_.swap(tqscore_exceeding_mz_rt_map_);
    std::unordered_map<int, double>().swap(new_mz_rt_map_);

    // remove expired entries for tqscore_exceeding_mass_rt_map_
    for (const auto& [m, r] : tqscore_exceeding_mass_rt_map_)
    {
      if (rt - r > rt_window_) { continue; }
      new_mass_rt_map_[m] = r;
    }
    new_mass_rt_map_.swap(tqscore_exceeding_mass_rt_map_);
    std::unordered_map<int, double>().swap(new_mass_rt_map_);

    // remove expired entries for all_mass_rt_map_, mass_qscore_map_
    for (const auto& item : all_mass_rt_map_)
    {
      if (rt - item.second > rt_window_) { continue; }
      new_all_mass_rt_map_[item.first] = item.second;
      new_mass_score_map_[item.first] = mass_qscore_map_[item.first];
    }
    new_all_mass_rt_map_.swap(all_mass_rt_map_);
    std::unordered_map<int, double>().swap(new_all_mass_rt_map_);
    new_mass_score_map_.swap(mass_qscore_map_);
    std::unordered_map<int, double>().swap(new_mass_score_map_);

    const int selection_phase_start = 0;
    const int selection_phase_end = 2; // inclusive
    // When selection_phase == 0, consider only the masses whose tqscore did not exceed total qscore threshold. 
    // when selection_phase == 1, consider all other masses for selection but the same m/z is avoided
    // when selection_phase == 2, consider all.
    // for target inclusive masses, qscore precursor snr threshold is not applied.
    // In all phase, for target exclusive mode, all the exclusive masses are excluded. For target inclusive mode, only the target masses are considered.

    for (int iteration = targeting_mode_ == 2 ? 0 : 1; iteration < 2; iteration++) 
    // for mass exclusion, first collect masses with exclusion list. Then collect without exclusion. This works the best
    {
      for (int selection_phase = selection_phase_start; selection_phase <= selection_phase_end; selection_phase++)
      {
        // Phase 0: targets (for inclusion mode) or tqscore-filtered masses
        // Phase 1: non-targets (only if non-strict inclusion mode and targets exist)
        if (selection_phase > 0)
        {
          // Allow phase 1 only for non-strict inclusion mode with active targets
          if (!(targeting_mode_ == 1 && !strict_inclusion_ && target_masses_.size() > 0 && selection_phase == 1))
          {
            break;
          }
        }

        // Iterate over candidates (sorted by qscore)
        for (const auto& pg : deconvolved_spectrum_)
        {
          // dont acquire the same mass multiple times
          if (selected_peak_groups_.size() >= mass_count) { break; }

          int charge;
          double score;
          int hcd = hcd_energy_;
          
          if (use_idscore_ && consider_all_Charge_states_ && hcd_energy_ < 0) {
            charge = pg.getBestIDScoreCharge();
            score = pg.getBestIDScore();
            hcd = pg.getBestIDScoreHCD();
          }
          else if (use_idscore_ && consider_all_Charge_states_) {
            charge = pg.getBestIDScoreChargeForHCD(hcd_energy_);
            score = pg.getBestIDScoreForHCD(hcd_energy_);
          }
          else if (use_idscore_ && !consider_all_Charge_states_ && hcd_energy_ < 0) {
            charge = pg.getRepAbsCharge();
            score = pg.getBestIDScoreForCharge(charge);
            hcd = pg.getBestHCDForCharge(charge);
          }
          else if (use_idscore_ && !consider_all_Charge_states_) {
            charge = pg.getRepAbsCharge();
            score = pg.getIDScoreForChargeAndHCD(charge, hcd_energy_);
          }
          else if (!use_idscore_ && consider_all_Charge_states_) {
            charge = pg.getBestQScoreCharge();
            score = pg.getBestQScore();
          }
          else {
            charge = pg.getRepAbsCharge();
            score = pg.getQscore();
          }
          
          double mass = pg.getMonoMass();

          auto [mz1, mz2] = pg.getMzRange(charge);

          double center_mz = (mz1 + mz2) / 2.0;

          mz1 -= optimal_window_margin_;
          mz2 += optimal_window_margin_;
          int integer_mz = (int)round(center_mz);

          int nominal_mass = SpectralDeconvolution::getNominalMass(mass);
          bool target_matched = false;
          double snr_threshold = snr_threshold_;
          double qscore_threshold = qscore_threshold_;
          double tqscore_factor_for_exclusion = 1.0;
          
          

          // Only triggered in exclusion mode
          if (iteration == 0)
          {
            auto inter = t_mass_score_map_.find(nominal_mass);
            if (inter != t_mass_score_map_.end()) { tqscore_factor_for_exclusion = t_mass_score_map_[nominal_mass]; }
            if (1 - tqscore_factor_for_exclusion > tqscore_threshold) { continue; }
          }

          // Inclusion mode (unified: supports TSV-based targets and legacy .log/.out targets)
          if (targeting_mode_ == 1 && target_masses_.size() > 0)
          {
            double delta = 2 * tol_[0] * mass * 1e-6;
            auto ub = std::upper_bound(target_masses_.begin(), target_masses_.end(), mass + delta);

            while (!target_matched)
            {
              if (ub != target_masses_.end())
              {
                if (std::abs(*ub - mass) < delta) // target is detected.
                {
                  // Check charge matching for both TSV and legacy modes
                  if (!inclusion_targets_.empty())  // TSV mode: check active_targets_ with charge
                  {
                    for (const auto* t : active_targets_)
                    {
                      if (t->charge < 0) {
                        target_matched = true;
                        break;
                      }
                      auto [min_charge, max_charge] = pg.getAbsChargeRange();
                      if ((t->charge >= min_charge) && (t->charge <= max_charge))
                      {
                        // Update with matched charge
                        charge = t->charge;
                        std::tie(mz1, mz2) = pg.getMzRange(charge);
                        center_mz = (mz1 + mz2) / 2.0;
                        mz1 -= optimal_window_margin_;
                        mz2 += optimal_window_margin_;
                        integer_mz = (int)round(center_mz);
                        
                        target_matched = true;
                        break;
                      }
                    }
                  }
                  else if (!target_mass_charge_map_.empty())  // Legacy .log mode with charge
                  {
                    auto it = target_mass_charge_map_.find(*ub);
                    if (it != target_mass_charge_map_.end())
                    {
                      charges = &it->second;
                      if (std::find(charges->begin(), charges->end(), charge) != charges->end())
                      {
                        target_matched = true;
                      }
                    }
                  }
                  else  // Legacy mode without charge (mass-only matching)
                  {
                    target_matched = true;
                  }

                  if (!target_matched)
                  {
                    if (ub == target_masses_.begin()) { break; }
                    ub--;
                    continue;
                  }
                }
                if (mass - *ub > delta) { break; }
              }
              if (ub == target_masses_.begin()) { break; }
              ub--;
            }

            if (target_matched)
            {
              snr_threshold = 0.0;
              qscore_threshold = 0.0; // stop exclusion for targets. todo tqscore lowest first? charge change.
            }
            else if (selection_phase == 0)
            {
              // Phase 0: only targets (strict behavior)
              continue;
            }
            // Phase 1: non-targets proceed with default thresholds
          }
          else if (targeting_mode_ == 1 && strict_inclusion_)
          {
            // Strict inclusion mode with no active targets - skip all candidates
            continue;
          }
          // deep mode
          else if (targeting_mode_ == 3 && excluded_masses_.size() > 0)
          {
            bool to_exclude = false;
            double delta = 2 * tol_[0] * mass * 1e-6;
            auto ub = std::upper_bound(excluded_masses_.begin(), excluded_masses_.end(), mass + delta);

            while (! to_exclude)
            {
              if (ub != excluded_masses_.end())
              {
                if (std::abs(*ub - mass) < delta) // target is detected.
                {
                  to_exclude = true;
                }
                if (mass - *ub > delta) { break; }
              }
              if (ub == excluded_masses_.begin()) { break; }
              ub--;
            }

            if (to_exclude) { continue; }
          }

          if (score < qscore_threshold) { break; }

          // TODO: Check
          if (pg.getChargeSNR(charge) < snr_threshold) { continue; }

          if (current_selected_mzs.find(center_mz) != current_selected_mzs.end()) // mz has been triggered
          {
            if (selection_phase < selection_phase_end) { continue; }
            if (! target_matched && current_selected_masses.find(pg.getMonoMass()) == current_selected_masses.end()) // but mass is different
            {
              continue;
            }
          }
          
          // selection phase 0, skip masses over tqscore threshold
          if (selection_phase < selection_phase_end - 1)
          { 
            if (tqscore_exceeding_mass_rt_map_.find(nominal_mass) != tqscore_exceeding_mass_rt_map_.end()
                || tqscore_exceeding_mz_rt_map_.find(integer_mz) != tqscore_exceeding_mz_rt_map_.end())
            {
              continue;
            }
          }
          
          // save mass acquisition
          all_mass_rt_map_[nominal_mass] = rt;

          if (!use_idscore_) {
            // Compute total qscore
            auto inter = mass_qscore_map_.find(nominal_mass);
            if (inter == mass_qscore_map_.end()) 
            { 
              mass_qscore_map_[nominal_mass] = score; 
            }
            else {
              // If mass has previously been acquired with higher qscore, skip
              if (score < mass_qscore_map_[nominal_mass]) {
                continue;
              } 
              mass_qscore_map_[nominal_mass] = score; 
            }

            // Add to exclusion list if neccessary
            if (mass_qscore_map_[nominal_mass] > tqscore_threshold)
            {
              tqscore_exceeding_mass_rt_map_[nominal_mass] = rt;
              tqscore_exceeding_mz_rt_map_[integer_mz] = rt;
            }
          }
          else {
            // Compute total qscore
            auto inter = mass_qscore_map_.find(nominal_mass);
            if (inter == mass_qscore_map_.end()) { mass_qscore_map_[nominal_mass] = 1 - score; }
            else { mass_qscore_map_[nominal_mass] *= 1 - score; }

            // Add to exclusion list if neccessary
            if (1 - mass_qscore_map_[nominal_mass] * tqscore_factor_for_exclusion > tqscore_threshold)
            {
              tqscore_exceeding_mass_rt_map_[nominal_mass] = rt;
              tqscore_exceeding_mz_rt_map_[integer_mz] = rt;
            }
          }

          // For legacy .log mode with charge targeting, remove this charge from list
          if (targeting_mode_ == 1 && charges != nullptr) {
            auto it = std::find(charges->begin(), charges->end(), charge);
            if (it != charges->end()) {
              charges->erase(it);
            }
          }

          // Store acquisition
          id_mass_map_[window_id_] = nominal_mass;
          id_mz_map_[window_id_] = integer_mz;
          id_qscore_map_[window_id_] = score;
          trigger_ids_.push_back(window_id_);
          window_id_++;

          selected_peak_groups_.push_back(pg);
          trigger_charges.push_back(charge);
          trigger_hcds.push_back(hcd);
          trigger_scores.push_back(score);

          trigger_left_isolation_mzs_.push_back(mz1);
          trigger_right_isolation_mzs_.push_back(mz2);
          current_selected_masses.insert(pg.getMonoMass());
          current_selected_mzs.insert(center_mz);
        }
      }
    }
  }

  void FLASHIda::removeFromExlusionList(int id)
  {
    // Check if id is valid
    if (id >= window_id_) { return; }

    // Obtain information needed for removal
    int nominal_mass = id_mass_map_[id];
    int integer_mz = id_mz_map_[id];
    double qscore = id_qscore_map_[id];

    // Remove from mass exclusion
    if (tqscore_exceeding_mass_rt_map_.find(nominal_mass) != tqscore_exceeding_mass_rt_map_.end())
    {
      tqscore_exceeding_mass_rt_map_.erase(nominal_mass);
    }

    // Remove from mz exclusion
    if (tqscore_exceeding_mz_rt_map_.find(integer_mz) != tqscore_exceeding_mz_rt_map_.end()) { tqscore_exceeding_mz_rt_map_.erase(integer_mz); }

    // Remove qscore from further calculations
    if (mass_qscore_map_.find(nominal_mass) != mass_qscore_map_.end()) { mass_qscore_map_[nominal_mass] /= 1 - qscore; }
  }

  void FLASHIda::getAllMonoisotopicMasses(double* masses, int length)
  {
    int len = std::min(length, (int)deconvolved_spectrum_.size());
    for (int i = 0; i < len; i++)
    {
      masses[i] = deconvolved_spectrum_[i].getMonoMass();
    }
  }

  int FLASHIda::GetAllPeakGroupSize()
  {
    return deconvolved_spectrum_.size();
  }

  int FLASHIda::deconvolveMS2(const double* mzs,
                              const double* ints,
                              int length,
                              double rt,
                              double precursor_mass,
                              int precursor_charge)
  {
    // Clear previous state
    ms2_deconvolved_spectrum_.clear();
    ms2_deconv_valid_ = false;
    ms2_deconv_rt_ = rt;

    if (length == 0)
    {
      return 0;
    }

    // Create MSSpectrum from input
    auto spec = makeMSSpectrum_(mzs, ints, length, rt, 2, "ms2_spectrum");

    // Create precursor PeakGroup - only set if precursor_mass > 0 AND precursor_charge != 0
    PeakGroup precursor_pg;
    if (precursor_mass > 0 && precursor_charge != 0)
    {
      int abs_charge = std::abs(precursor_charge);
      bool is_positive = precursor_charge > 0;

      // Calculate precursor m/z from mass and charge
      double charge_mass = FLASHHelperClasses::getChargeMass(is_positive);
      double precursor_mz = (precursor_mass + abs_charge * charge_mass) / abs_charge;

      // Set precursor on the MSSpectrum (required for deconvolution mass range calculation)
      Precursor precursor;
      precursor.setMZ(precursor_mz);
      precursor.setCharge(precursor_charge);
      spec.getPrecursors().push_back(precursor);

      // Construct PeakGroup with proper charge range and polarity
      precursor_pg = PeakGroup(abs_charge, abs_charge, is_positive);
      precursor_pg.push_back(FLASHHelperClasses::LogMzPeak());
      precursor_pg.setMonoisotopicMass(precursor_mass);
      precursor_pg.setRepAbsCharge(abs_charge);
      precursor_pg.setQscore(1.0);  // Known precursor from MS1, high confidence
      precursor_pg.setSNR(1.0);
    }

    // Perform deconvolution (empty precursor_pg if mass <= 0 or charge == 0)
    fd_.performSpectrumDeconvolution(spec, 0, precursor_pg);
    ms2_deconvolved_spectrum_ = fd_.getDeconvolvedSpectrum();

    if (ms2_deconvolved_spectrum_.empty())
    {
      return 0;
    }

    // Sort by qscore (highest first) for getBestMS2Masses
    ms2_deconvolved_spectrum_.sortByQscore();
    ms2_deconv_valid_ = true;

    return static_cast<int>(ms2_deconvolved_spectrum_.size());
  }

  int FLASHIda::deconvolveMS2Py(const std::vector<double>& mzs,
                                const std::vector<double>& ints,
                                double rt,
                                double precursor_mass,
                                int precursor_charge)
  {
    if (mzs.empty() || mzs.size() != ints.size())
    {
      return 0;
    }
    return deconvolveMS2(mzs.data(), ints.data(), static_cast<int>(mzs.size()), rt, precursor_mass, precursor_charge);
  }

  int FLASHIda::getBestMS2Masses(int n,
                                 double* masses,
                                 double* qscores,
                                 int* charges,
                                 double* window_starts,
                                 double* window_ends)
  {
    if (!ms2_deconv_valid_ || ms2_deconvolved_spectrum_.empty())
    {
      return 0;
    }

    std::sort(ms2_deconvolved_spectrum_.begin(), ms2_deconvolved_spectrum_.end(),
    [](const PeakGroup& a, const PeakGroup& b) {
      return a.getChargeIntensity(a.getMaxIntensityAbsCharge()) > b.getChargeIntensity(b.getMaxIntensityAbsCharge());
    });

    int output_idx = 0;
    for (Size pg_idx = 0; pg_idx < ms2_deconvolved_spectrum_.size() && output_idx < n; ++pg_idx)
    {
      const auto& pg = ms2_deconvolved_spectrum_[pg_idx];
      double mono_mass = pg.getMonoMass();

      if (ms3_all_charges_)
      {
        // Collect all charges with non-zero intensity and sort by intensity descending
        auto [min_c, max_c] = pg.getAbsChargeRange();
        std::vector<std::pair<float, int>> charge_intensities;  // (intensity, charge)
        for (int c = min_c; c <= max_c; ++c)
        {
          float intensity = pg.getChargeIntensity(c);
          if (intensity > 0)
          {
            charge_intensities.emplace_back(intensity, c);
          }
        }
        std::sort(charge_intensities.begin(), charge_intensities.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        for (const auto& [intensity, c] : charge_intensities)
        {
          if (output_idx >= n) break;
          masses[output_idx] = mono_mass;
          qscores[output_idx] = intensity;
          charges[output_idx] = c;
          auto [mz1, mz2] = pg.getMzRange(c);
          window_starts[output_idx] = mz1 - optimal_window_margin_;
          window_ends[output_idx] = mz2 + optimal_window_margin_;
          ++output_idx;
        }
        if (output_idx >= n) break;
      }
      else
      {
        // Single charge behavior (most intense charge)
        int charge = pg.getMaxIntensityAbsCharge();
        masses[output_idx] = mono_mass;
        qscores[output_idx] = pg.getQscore();
        charges[output_idx] = charge;
        auto [mz1, mz2] = pg.getMzRange(charge);
        window_starts[output_idx] = mz1 - optimal_window_margin_;
        window_ends[output_idx] = mz2 + optimal_window_margin_;
        ++output_idx;
      }
    }

    return output_idx;
  }

  int FLASHIda::getBestMS2MassesPy(int n,
                                   std::vector<double>& masses,
                                   std::vector<double>& qscores,
                                   std::vector<int>& charges,
                                   std::vector<double>& window_starts,
                                   std::vector<double>& window_ends)
  {
    masses.resize(n);
    qscores.resize(n);
    charges.resize(n);
    window_starts.resize(n);
    window_ends.resize(n);

    int count = getBestMS2Masses(n, masses.data(), qscores.data(), charges.data(),
                                 window_starts.data(), window_ends.data());

    masses.resize(count);
    qscores.resize(count);
    charges.resize(count);
    window_starts.resize(count);
    window_ends.resize(count);

    return count;
  }

  bool FLASHIda::hasMS2Deconvolution() const
  {
    return ms2_deconv_valid_;
  }

  int FLASHIda::getMS2PeakGroupCount() const
  {
    return ms2_deconv_valid_ ? static_cast<int>(ms2_deconvolved_spectrum_.size()) : 0;
  }

  void FLASHIda::clearMS2Deconvolution()
  {
    ms2_deconvolved_spectrum_.clear();
    ms2_deconv_valid_ = false;
    ms2_deconv_rt_ = -1.0;
  }

  int FLASHIda::runTagBasedFragmentMatching_(const String& protein_sequence,
                                              std::vector<TagBasedFragmentMatch>& matches,
                                              std::vector<PTMSite>* ptm_sites,
                                              const String& fragmentation_method)
  {
    matches.clear();
    if (ptm_sites != nullptr)
    {
      ptm_sites->clear();
    }

    if (!ms2_deconv_valid_ || ms2_deconvolved_spectrum_.empty() || protein_sequence.empty())
    {
      return 0;
    }

    // 1. Copy and sort deconvolved spectrum
    DeconvolvedSpectrum dspec = ms2_deconvolved_spectrum_;
    dspec.sort();
    double precursor_mass = dspec.getPrecursorPeakGroup().getMonoMass();
    std::cout << "PM=" << precursor_mass << std::endl;
    for (const auto& pg : dspec)
    {
      double mono_mass = pg.getMonoMass();
      std::cout << mono_mass << std::endl;
    }

    // 2. Run FLASHTagger to generate sequence tags
    double ppm_tolerance = tol_[1];
    std::vector<std::string> ion_types_str = getIonTypesForFragmentationMethod(fragmentation_method);
    FLASHTaggerAlgorithm tagger;
    Param tagger_param = tagger.getDefaults();
    tagger_param.setValue("ion_type", ion_types_str);
    tagger.setParameters(tagger_param);
    tagger.run(dspec, ppm_tolerance);

    std::vector<FLASHHelperClasses::Tag> tags;
    tagger.fillTags(tags);

    if (tags.empty())
    {
      // No tags found - return empty
      std::cout << "no-tags" << std::endl;
      return 0;
    }
    std::cout << "tags" << std::endl;

    // 3. Create ProteinHit with proper metadata for matching
    String clean_seq = protein_sequence;
    std::replace(clean_seq.begin(), clean_seq.end(), 'I', 'L');

    ProteinHit hit(0.0, 0, "input_protein", protein_sequence);
    hit.setMetaValue("Scan", 1);
    hit.setMetaValue("FastaIndex", 0);

    // Find tag positions in the protein
    std::vector<int> tag_positions;
    std::set<int> tag_indices_set;
    for (Size i = 0; i < tags.size(); ++i)
    {
      String tag_seq = tags[i].getUppercaseSequence();
      Size pos = clean_seq.find(tag_seq);
      if (pos != String::npos)
      {
        tag_positions.push_back(static_cast<int>(pos));
        tag_indices_set.insert(static_cast<int>(i));
      }
    }

    if (tag_positions.empty())
    {
      // No tags match the protein - return empty
      std::cout << "no-matching-tags" << std::endl;
      return 0;
    }
    std::cout << "matching-tags" << std::endl;

    hit.setMetaValue("TagIndices", std::vector<int>(tag_indices_set.begin(), tag_indices_set.end()));
    hit.setMetaValue("TagPositions", tag_positions);
    hit.setMetaValue("MatchedAA", static_cast<int>(tag_positions.size()));
    hit.setCoverage(static_cast<double>(tag_positions.size()) / protein_sequence.size());

    // Calculate flanking masses for tag matching
    std::set<int> n_flanking_masses_set, c_flanking_masses_set;
    for (Size i = 0; i < tags.size(); ++i)
    {
      String tag_seq = tags[i].getUppercaseSequence();
      Size pos = clean_seq.find(tag_seq);
      if (pos != String::npos)
      {
        double n_mass = 0;
        for (Size j = 0; j < pos; ++j)
        {
          n_mass += ResidueDB::getInstance()->getResidue(clean_seq[j])->getMonoWeight(Residue::Internal);
        }
        double n_diff = static_cast<int>(std::round(n_mass - tags[i].getNtermMass()));
        n_flanking_masses_set.insert(n_diff);

        double c_mass = 0;
        Size tag_end = pos + tag_seq.size();
        for (Size j = tag_end; j < clean_seq.size(); ++j)
        {
          c_mass += ResidueDB::getInstance()->getResidue(clean_seq[j])->getMonoWeight(Residue::Internal);
        }
        double c_diff = static_cast<int>(std::round(c_mass - tags[i].getCtermMass()));
        c_flanking_masses_set.insert(c_diff);
      }
    }
    hit.setMetaValue("NtermFlankingMasses", std::vector<int>(n_flanking_masses_set.begin(), n_flanking_masses_set.end()));
    hit.setMetaValue("CtermFlankingMasses", std::vector<int>(c_flanking_masses_set.begin(), c_flanking_masses_set.end()));

    std::vector<ProteinHit> hits;
    hits.push_back(hit);

    // 4. Create spec_vec (integer masses from spectrum)
    std::vector<int> spec_vec;
    spec_vec.reserve(dspec.size() + 1);
    spec_vec.push_back(0);
    for (const auto& pg : dspec)
    {
      spec_vec.push_back(static_cast<int>(std::round(pg.getMonoMass())));
    }

    // 5. Create vec_pro and rev_vec_pro for the protein
    std::unordered_set<int> vec, rev_vec;
    double mass = 0;
    vec.insert(0);
    rev_vec.insert(0);

    for (Size j = 0; j < clean_seq.size(); ++j)
    {
      mass += ResidueDB::getInstance()->getResidue(clean_seq[j])->getMonoWeight(Residue::Internal);
      vec.insert(static_cast<int>(std::round(mass)));
    }

    mass = 0;
    for (Size j = clean_seq.size(); j > 0; --j)
    {
      mass += ResidueDB::getInstance()->getResidue(clean_seq[j - 1])->getMonoWeight(Residue::Internal);
      rev_vec.insert(static_cast<int>(std::round(mass)));
    }

    std::vector<std::unordered_set<int>> vec_pro = {vec};
    std::vector<std::unordered_set<int>> rev_vec_pro = {rev_vec};

    // 6. Run FLASHTagger matching
    double max_mod_mass = 700.0;
    std::cout << "max_mod_mass=" << max_mod_mass << std::endl;
    std::cout << "ion_types_str=";
    for (const auto& ion : ion_types_str) {
      std::cout << ion << " ";
    }
    std::cout << std::endl;

    FLASHTaggerAlgorithm::runMatching(hits, dspec, spec_vec, vec_pro, rev_vec_pro, max_mod_mass);

    if (hits.empty())
    {
      std::cout << "no-hits" << std::endl;
      return 0;
    }
    std::cout << "matching-hits" << std::endl;

    // 7. Run FLASHExtender for path-based validation
    FLASHExtenderAlgorithm extender;
    Param extender_param = extender.getDefaults();
    extender_param.setValue("ion_type", ion_types_str);
    extender_param.setValue("max_mod_mass", max_mod_mass);
    extender.setParameters(extender_param);
    extender.run(hits, dspec, spec_vec, vec_pro, rev_vec_pro, tags, ppm_tolerance, false);

    std::vector<ProteinHit> proteoform_hits;
    extender.fillProteoforms(proteoform_hits);

    if (proteoform_hits.empty())
    {
      std::cout << "no-proteoform-hits" << std::endl;
      return 0;
    }
    std::cout << "matching-proteoform-hits" << std::endl;

    // === Diagnostic Output: Proteoform and PTM sites ===
    const auto& best_hit = proteoform_hits[0];

    // Print matched proteoform sequence (using FLASHExtender truncation info)
    // Note: FLASHExtender positions are 1-based, convert to 0-based for substr
    int start_pos = -1, end_pos = -1;
    if (best_hit.metaValueExists("StartPosition"))
      start_pos = static_cast<int>(best_hit.getMetaValue("StartPosition")) - 1;
    if (best_hit.metaValueExists("EndPosition"))
      // TODO: Why not end??
      end_pos = static_cast<int>(best_hit.getMetaValue("EndPosition"));

    if (start_pos >= 0 && end_pos >= 0 && end_pos > start_pos)
    {
      // Truncated proteoform detected by FLASHExtender
      String truncated_seq = protein_sequence.substr(start_pos, end_pos - start_pos);
      String display_seq = truncated_seq.size() > 30
          ? truncated_seq.substr(0, 27) + "..."
          : truncated_seq;
      std::cout << "[runTagBasedFragmentMatching_] Matched proteoform: " << display_seq
                << " (residues " << (start_pos + 1) << "-" << end_pos
                << " of " << protein_sequence.size() << ")" << std::endl;
    }
    else
    {
      // Full sequence (no truncation detected)
      String display_seq = protein_sequence.size() > 30
          ? protein_sequence.substr(0, 15) + "..." + protein_sequence.substr(protein_sequence.size() - 12)
          : protein_sequence;
      std::cout << "[runTagBasedFragmentMatching_] Matched proteoform: " << display_seq
                << " (full, " << protein_sequence.size() << " aa)" << std::endl;
    }

    // Extract PTM sites from FLASHExtender metadata
    std::vector<double> mod_masses;
    std::vector<int> mod_starts;
    std::vector<int> mod_ends;

    if (best_hit.metaValueExists("Modifications") &&
        best_hit.metaValueExists("ModificationStarts") &&
        best_hit.metaValueExists("ModificationEnds"))
    {
      mod_masses = best_hit.getMetaValue("Modifications");
      mod_starts = best_hit.getMetaValue("ModificationStarts");
      mod_ends = best_hit.getMetaValue("ModificationEnds");
      // Convert positions - TODO: Why??
      for (auto& pos : mod_starts) pos += 1;
      for (auto& pos : mod_ends) pos += 1;
    }

    // Print diagnostic output for PTM sites
    std::cout << "[runTagBasedFragmentMatching_] PTM sites detected: " << mod_masses.size() << std::endl;
    for (Size i = 0; i < mod_masses.size(); ++i)
    {
      int idx_start = std::max(0, mod_starts[i] - 1);
      int idx_end = std::min(static_cast<int>(protein_sequence.size()), mod_ends[i]);
      String subseq = protein_sequence.substr(idx_start, idx_end - idx_start);
      std::cout << "  Residue " << mod_starts[i] << "-" << mod_ends[i] << ": "
                << std::showpos << mod_masses[i] << std::noshowpos << " Da "
                << "\"" << subseq << "\"" << std::endl;
    }

    // Populate ptm_sites output if requested
    if (ptm_sites != nullptr)
    {
      for (Size i = 0; i < mod_masses.size(); ++i)
      {
        PTMSite site;
        site.start_position = mod_starts[i];
        site.end_position = mod_ends[i];
        site.position = (site.start_position + site.end_position) / 2;
        site.mass_shift = mod_masses[i];
        ptm_sites->push_back(site);
      }
    }

    // 8. Build local PTM sites for fragment matching
    std::vector<PTMSite> local_ptm_sites;
    for (Size i = 0; i < mod_masses.size(); ++i)
    {
      PTMSite site;
      site.start_position = mod_starts[i];
      site.end_position = mod_ends[i];
      site.position = (site.start_position + site.end_position) / 2;
      site.mass_shift = mod_masses[i];
      local_ptm_sites.push_back(site);
    }

    // 9. Handle truncation: use truncated sequence if FLASHExtender detected it
    String matching_sequence = protein_sequence;

    if (start_pos >= 0 && end_pos > start_pos)
    {
      // Use truncated sequence for fragment matching
      matching_sequence = protein_sequence.substr(start_pos, end_pos - start_pos);

      // Adjust PTM positions to be relative to truncated sequence
      for (auto& site : local_ptm_sites)
      {
        site.start_position -= start_pos;
        site.end_position -= start_pos;
        site.position -= start_pos;

        // Clamp to valid range within truncated sequence
        site.start_position = std::max(1, site.start_position);
        site.end_position = std::min(static_cast<int>(matching_sequence.size()), site.end_position);
      }

      std::cout << "[runTagBasedFragmentMatching_] Using truncated sequence for matching ("
                << matching_sequence.size() << " aa, offset " << start_pos << ")" << std::endl;
    }

    // 10. Calculate PTM-adjusted theoretical fragment masses for all configured ion types
    std::map<char, std::vector<double>> fragment_masses_map;
    calculatePTMAdjustedFragmentMassesMulti(matching_sequence, local_ptm_sites, ion_types_str, fragment_masses_map);

    // 11. Match observed masses against PTM-adjusted theoretical masses (all ion types)
    for (Size peak_idx = 0; peak_idx < ms2_deconvolved_spectrum_.size(); ++peak_idx)
    {
      const auto& pg = ms2_deconvolved_spectrum_[peak_idx];
      double observed_mass = pg.getMonoMass();

      // Find best match across all configured ion types
      char best_ion_type = '\0';
      int best_frag_idx = -1;
      double best_ppm = ppm_tolerance;
      double best_theo = 0;

      for (const auto& [ion_char, masses] : fragment_masses_map)
      {
        for (Size i = 0; i < masses.size(); ++i)
        {
          double theo_mass = masses[i];
          if (theo_mass <= 0) continue;
          double ppm_error = std::abs(observed_mass - theo_mass) / theo_mass * 1e6;
          if (ppm_error < best_ppm)
          {
            best_ppm = ppm_error;
            best_ion_type = ion_char;
            best_frag_idx = static_cast<int>(i) + 1;  // 1-based index
            best_theo = theo_mass;
          }
        }
      }

      // Store best match if found
      if (best_ion_type != '\0')
      {
        TagBasedFragmentMatch match;
        match.peak_index = static_cast<int>(peak_idx);
        match.observed_mass = observed_mass;
        match.qscore = pg.getChargeIntensity(pg.getMaxIntensityAbsCharge());
        match.charge = pg.getMaxIntensityAbsCharge();
        match.fragment_index = best_frag_idx;
        match.ion_type = best_ion_type;
        match.theoretical_mass = best_theo;
        match.ppm_error = best_ppm;
        matches.push_back(match);
      }
    }

    // 11. Sort by qscore descending
    std::sort(matches.begin(), matches.end(),
              [](const TagBasedFragmentMatch& a, const TagBasedFragmentMatch& b) {
                return a.qscore > b.qscore;
              });

    // === Diagnostic Output: Matched fragment ions ===
    std::cout << "[runTagBasedFragmentMatching_] Matched " << matches.size() << " fragment ions:" << std::endl;
    Size display_count = std::min(Size(10), matches.size());
    for (Size i = 0; i < display_count; ++i)
    {
      const auto& m = matches[i];
      std::cout << "  " << m.ion_type << m.fragment_index
                << " - " << std::fixed << std::setprecision(2) << m.observed_mass << " Da"
                << " (qscore: " << std::setprecision(2) << m.qscore << ")" << std::endl;
    }
    if (matches.size() > 10)
    {
      std::cout << "  ... (" << (matches.size() - 10) << " more)" << std::endl;
    }

    return static_cast<int>(matches.size());
  }

  int FLASHIda::getTopFragmentMatches(const String& protein_sequence,
                                      int n,
                                      double* masses,
                                      double* qscores,
                                      int* charges,
                                      double* window_starts,
                                      double* window_ends,
                                      char* ion_types,
                                      int* fragment_indices,
                                      const String& fragmentation_method)
  {
    std::cout << "Matching fragments!" << std::endl;
    // Use tag-based matching workflow (FLASHTagger + FLASHExtender)
    std::vector<TagBasedFragmentMatch> matches;
    runTagBasedFragmentMatching_(protein_sequence, matches, nullptr, fragmentation_method);

    int output_idx = 0;
    for (size_t i = 0; i < matches.size() && output_idx < n; ++i)
    {
      const auto& m = matches[i];
      const auto& pg = ms2_deconvolved_spectrum_[m.peak_index];

      if (ms3_all_charges_)
      {
        // Collect all charges with non-zero intensity and sort by intensity descending
        auto [min_c, max_c] = pg.getAbsChargeRange();
        std::vector<std::pair<float, int>> charge_intensities;  // (intensity, charge)
        for (int c = min_c; c <= max_c; ++c)
        {
          float intensity = pg.getChargeIntensity(c);
          if (intensity > 0)
          {
            charge_intensities.emplace_back(intensity, c);
          }
        }
        std::sort(charge_intensities.begin(), charge_intensities.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        for (const auto& [intensity, c] : charge_intensities)
        {
          if (output_idx >= n) break;
          masses[output_idx] = m.observed_mass;
          qscores[output_idx] = intensity;
          charges[output_idx] = c;
          ion_types[output_idx] = m.ion_type;
          fragment_indices[output_idx] = m.fragment_index;
          auto [mz1, mz2] = pg.getMzRange(c);
          window_starts[output_idx] = mz1 - optimal_window_margin_;
          window_ends[output_idx] = mz2 + optimal_window_margin_;
          ++output_idx;
        }
        if (output_idx >= n) break;
      }
      else
      {
        // Single charge behavior
        masses[output_idx] = m.observed_mass;
        qscores[output_idx] = m.qscore;
        charges[output_idx] = m.charge;
        ion_types[output_idx] = m.ion_type;
        fragment_indices[output_idx] = m.fragment_index;
        auto [mz1, mz2] = pg.getMzRange(m.charge);
        window_starts[output_idx] = mz1 - optimal_window_margin_;
        window_ends[output_idx] = mz2 + optimal_window_margin_;
        ++output_idx;
      }
    }

    return output_idx;
  }

  int FLASHIda::getTopFragmentMatchesPy(const String& protein_sequence,
                                        int n,
                                        std::vector<double>& masses,
                                        std::vector<double>& qscores,
                                        std::vector<int>& charges,
                                        std::vector<double>& window_starts,
                                        std::vector<double>& window_ends,
                                        std::vector<int>& is_b_ions,
                                        std::vector<int>& fragment_indices)
  {
    masses.resize(n);
    qscores.resize(n);
    charges.resize(n);
    window_starts.resize(n);
    window_ends.resize(n);

    // Use raw char array for the C-style function
    std::unique_ptr<char[]> ion_types_temp(new char[n]);
    fragment_indices.resize(n);

    int count = getTopFragmentMatches(protein_sequence, n, masses.data(), qscores.data(),
                                      charges.data(), window_starts.data(), window_ends.data(),
                                      ion_types_temp.get(), fragment_indices.data(), "HCD");

    masses.resize(count);
    qscores.resize(count);
    charges.resize(count);
    window_starts.resize(count);
    window_ends.resize(count);
    fragment_indices.resize(count);

    // Convert char array to int vector (1 for prefix ion, 0 for suffix ion)
    is_b_ions.resize(count);
    for (int i = 0; i < count; ++i)
    {
      char c = ion_types_temp[i];
      is_b_ions[i] = (c == 'a' || c == 'b' || c == 'c') ? 1 : 0;
    }

    return count;
  }

  int FLASHIda::getAmbiguityEnclosingIons(const String& protein_sequence,
                                          int n,
                                          double* masses,
                                          double* qscores,
                                          int* charges,
                                          double* window_starts,
                                          double* window_ends,
                                          char* ion_types,
                                          int* fragment_indices,
                                          const String& fragmentation_method)
  {
    // Get fragment matches AND PTM sites from FLASHExtender
    std::vector<TagBasedFragmentMatch> fragment_ion_match;
    std::vector<PTMSite> ptm_sites;
    int match_count = runTagBasedFragmentMatching_(protein_sequence, fragment_ion_match, &ptm_sites, fragmentation_method);

    if (match_count == 0)
    {
      std::cout << "[getAmbiguityEnclosingIons] No fragment matches found" << std::endl;
      return 0;
    }
    std::cout << "[getAmbiguityEnclosingIons] Found " << match_count << " fragment matches" << std::endl;

    if (ptm_sites.empty())
    {
      std::cout << "[getAmbiguityEnclosingIons] No PTM sites detected by FLASHExtender" << std::endl;
      return 0;
    }

    // Debug output for PTM sites
    std::cout << "[getAmbiguityEnclosingIons] FLASHExtender detected " << ptm_sites.size() << " PTM sites:" << std::endl;
    for (const auto& site : ptm_sites)
    {
      int start_idx = std::max(0, site.start_position - 1);
      int end_idx = std::min(static_cast<int>(protein_sequence.size()), site.end_position);
      String subsequence = protein_sequence.substr(start_idx, end_idx - start_idx);
      bool is_ambiguous = site.start_position != site.end_position;
      std::cout << "  PTM [" << site.start_position << "-" << site.end_position << "] "
                << "\"" << subsequence << "\" "
                << std::showpos << site.mass_shift << std::noshowpos << " Da"
                << (is_ambiguous ? " (AMBIGUOUS)" : " (localized)") << std::endl;
    }

    // New data structures for interleaved PTM site output
    struct EnclosingIon {
      float qscore;
      int peak_index;
      char ion_type;
      int fragment_index;
      int n_terminal_pos;  // For boundary checking: (L-fragment_index+1) for suffix, 1 for prefix
      int c_terminal_pos;  // For boundary checking: L for suffix, fragment_index for prefix
    };

    struct PTMBrackets {
      std::vector<EnclosingIon> left_ions;   // All valid left brackets, sorted by distance (closest first)
      std::vector<EnclosingIon> right_ions;  // All valid right brackets, sorted by distance (closest first)
      int ptm_start;
      int ptm_end;
      float priority;
    };

    std::set<int> used_peaks;
    std::vector<PTMBrackets> ptm_brackets;
    int seq_len = static_cast<int>(protein_sequence.size());

    // Collect ALL valid brackets per PTM site
    for (const auto& site : ptm_sites)
    {
      PTMBrackets brackets;
      brackets.ptm_start = site.start_position;
      brackets.ptm_end = site.end_position;

      // Collect ALL valid left brackets (suffix ions), sorted by distance to PTM
      std::vector<std::pair<int, const TagBasedFragmentMatch*>> left_candidates;  // (distance, match)
      for (const auto& m : fragment_ion_match)
      {
        if (!isPrefixIon(m.ion_type))  // suffix ion (y, x, z)
        {
          // suffix ion covers positions (L-n+1)..L where n=fragment_index
          int min_y = seq_len - site.end_position;
          if (m.fragment_index >= min_y)
          {
            int distance = m.fragment_index + 1 - min_y;
            left_candidates.emplace_back(distance, &m);
          }
        }
      }
      std::sort(left_candidates.begin(), left_candidates.end());  // Sort by distance (closest first)

      for (const auto& [dist, m] : left_candidates)
      {
        int n_term = seq_len - m->fragment_index + 1;  // N-terminal position
        brackets.left_ions.push_back({static_cast<float>(m->qscore), m->peak_index, m->ion_type,
                                      m->fragment_index, n_term, seq_len});
      }

      // Collect ALL valid right brackets (prefix ions), sorted by distance to PTM
      std::vector<std::pair<int, const TagBasedFragmentMatch*>> right_candidates;
      for (const auto& m : fragment_ion_match)
      {
        if (isPrefixIon(m.ion_type))  // prefix ion (b, a, c)
        {
          // prefix ion covers positions 1..fragment_index
          if (m.fragment_index >= site.end_position - 1)
          {
            int distance = m.fragment_index - site.end_position;
            right_candidates.emplace_back(distance, &m);
          }
        }
      }
      std::sort(right_candidates.begin(), right_candidates.end());

      for (const auto& [dist, m] : right_candidates)
      {
        brackets.right_ions.push_back({static_cast<float>(m->qscore), m->peak_index, m->ion_type,
                                       m->fragment_index, 1, m->fragment_index});
      }

      // Priority = max qscore of primary brackets
      brackets.priority = 0.0f;
      if (!brackets.left_ions.empty()) brackets.priority = std::max(brackets.priority, brackets.left_ions[0].qscore);
      if (!brackets.right_ions.empty()) brackets.priority = std::max(brackets.priority, brackets.right_ions[0].qscore);

      // Debug output for this PTM site
      std::cout << "[getAmbiguityEnclosingIons] PTM [" << site.start_position << "-" << site.end_position << "]:" << std::endl;
      if (!brackets.left_ions.empty())
      {
        std::cout << "  Left bracket (primary): " << brackets.left_ions[0].ion_type << brackets.left_ions[0].fragment_index
                  << " (qscore: " << brackets.left_ions[0].qscore << ", total candidates: " << brackets.left_ions.size() << ")" << std::endl;
      }
      else
      {
        std::cout << "  Left bracket: none found" << std::endl;
      }
      if (!brackets.right_ions.empty())
      {
        std::cout << "  Right bracket (primary): " << brackets.right_ions[0].ion_type << brackets.right_ions[0].fragment_index
                  << " (qscore: " << brackets.right_ions[0].qscore << ", total candidates: " << brackets.right_ions.size() << ")" << std::endl;
      }
      else
      {
        std::cout << "  Right bracket: none found" << std::endl;
      }

      if (!brackets.left_ions.empty() || !brackets.right_ions.empty())
        ptm_brackets.push_back(brackets);
    }

    // Sort PTM sites by priority (highest first)
    std::sort(ptm_brackets.begin(), ptm_brackets.end(),
              [](const auto& a, const auto& b) { return a.priority > b.priority; });

    // Helper to check if ion covers a PTM
    auto covers_ptm = [&](const EnclosingIon& ion, const PTMBrackets& ptm) -> bool {
      if (isPrefixIon(ion.ion_type)) {
        return ion.c_terminal_pos >= ptm.ptm_start;  // Prefix covers if C-term >= PTM start
      } else {
        return ion.n_terminal_pos <= ptm.ptm_end;    // Suffix covers if N-term <= PTM end
      }
    };

    // Helper to check if secondary covers an ADDITIONAL PTM that primary doesn't
    auto crosses_additional_ptm = [&](const EnclosingIon& secondary,
                                       const EnclosingIon& primary,
                                       const PTMBrackets& current) -> bool {
      for (const auto& other : ptm_brackets) {
        if (&other == &current) continue;
        if (covers_ptm(secondary, other) && !covers_ptm(primary, other)) {
          return true;  // Secondary covers a PTM that primary doesn't
        }
      }
      return false;
    };

    // Helper to output one ion (handles MS3AllCharges)
    int output_idx = 0;
    auto output_ion = [&](const EnclosingIon& ion) {
      if (used_peaks.find(ion.peak_index) != used_peaks.end()) return;  // Already output
      used_peaks.insert(ion.peak_index);

      const auto& pg = ms2_deconvolved_spectrum_[ion.peak_index];
      double mono_mass = pg.getMonoMass();

      if (ms3_all_charges_)
      {
        // Collect all charges with non-zero intensity and sort by intensity descending
        auto [min_c, max_c] = pg.getAbsChargeRange();
        std::vector<std::pair<float, int>> charge_intensities;
        for (int c = min_c; c <= max_c; ++c)
        {
          float intensity = pg.getChargeIntensity(c);
          if (intensity > 0)
          {
            charge_intensities.emplace_back(intensity, c);
          }
        }
        std::sort(charge_intensities.begin(), charge_intensities.end(),
                  [](const auto& a, const auto& b) { return a.first > b.first; });

        for (const auto& [intensity, c] : charge_intensities)
        {
          if (output_idx >= n) break;
          masses[output_idx] = mono_mass;
          qscores[output_idx] = intensity;
          charges[output_idx] = c;
          ion_types[output_idx] = ion.ion_type;
          fragment_indices[output_idx] = ion.fragment_index;
          auto [mz1, mz2] = pg.getMzRange(c);
          window_starts[output_idx] = mz1 - optimal_window_margin_;
          window_ends[output_idx] = mz2 + optimal_window_margin_;
          ++output_idx;
        }
      }
      else
      {
        // Single charge behavior
        int charge = pg.getMaxIntensityAbsCharge();
        masses[output_idx] = mono_mass;
        qscores[output_idx] = ion.qscore;
        charges[output_idx] = charge;
        ion_types[output_idx] = ion.ion_type;
        fragment_indices[output_idx] = ion.fragment_index;
        auto [mz1, mz2] = pg.getMzRange(charge);
        window_starts[output_idx] = mz1 - optimal_window_margin_;
        window_ends[output_idx] = mz2 + optimal_window_margin_;
        ++output_idx;
      }
    };

    std::cout << "[getAmbiguityEnclosingIons] Output Phase 1: Primary brackets for " << ptm_brackets.size() << " PTM sites" << std::endl;

    // Phase 1: Primary brackets (first ion per PTM site)
    for (auto& brackets : ptm_brackets)
    {
      if (output_idx >= n) break;

      // Output primary left bracket
      if (!brackets.left_ions.empty())
      {
        output_ion(brackets.left_ions[0]);
      }
      if (output_idx >= n) break;

      // Output primary right bracket
      if (!brackets.right_ions.empty())
      {
        output_ion(brackets.right_ions[0]);
      }
    }

    std::cout << "[getAmbiguityEnclosingIons] After Phase 1: " << output_idx << " ions output" << std::endl;

    // Phase 2: Secondary brackets (additional ions that don't cover ADDITIONAL PTMs)
    std::cout << "[getAmbiguityEnclosingIons] Output Phase 2: Secondary brackets" << std::endl;
    for (size_t ion_idx = 1; output_idx < n; ++ion_idx)
    {
      bool any_output = false;
      for (auto& brackets : ptm_brackets)
      {
        if (output_idx >= n) break;

        // Output secondary left bracket if available and doesn't cover additional PTM
        if (ion_idx < brackets.left_ions.size() && !brackets.left_ions.empty())
        {
          const auto& secondary = brackets.left_ions[ion_idx];
          const auto& primary = brackets.left_ions[0];
          if (used_peaks.find(secondary.peak_index) == used_peaks.end() &&
              !crosses_additional_ptm(secondary, primary, brackets))
          {
            output_ion(secondary);
            any_output = true;
          }
        }
        if (output_idx >= n) break;

        // Output secondary right bracket if available and doesn't cover additional PTM
        if (ion_idx < brackets.right_ions.size() && !brackets.right_ions.empty())
        {
          const auto& secondary = brackets.right_ions[ion_idx];
          const auto& primary = brackets.right_ions[0];
          if (used_peaks.find(secondary.peak_index) == used_peaks.end() &&
              !crosses_additional_ptm(secondary, primary, brackets))
          {
            output_ion(secondary);
            any_output = true;
          }
        }
      }
      if (!any_output) break;  // No more secondary ions available
    }

    std::cout << "[getAmbiguityEnclosingIons] Total output: " << output_idx << " ions" << std::endl;

    return output_idx;
  }

  int FLASHIda::getAmbiguityEnclosingIonsPy(const String& protein_sequence,
                                            int n,
                                            std::vector<double>& masses,
                                            std::vector<double>& qscores,
                                            std::vector<int>& charges,
                                            std::vector<double>& window_starts,
                                            std::vector<double>& window_ends,
                                            std::vector<int>& is_b_ions,
                                            std::vector<int>& fragment_indices)
  {
    masses.resize(n);
    qscores.resize(n);
    charges.resize(n);
    window_starts.resize(n);
    window_ends.resize(n);

    // Use raw char array for the C-style function
    std::unique_ptr<char[]> ion_types_temp(new char[n]);
    fragment_indices.resize(n);

    int count = getAmbiguityEnclosingIons(protein_sequence, n, masses.data(), qscores.data(),
                                          charges.data(), window_starts.data(), window_ends.data(),
                                          ion_types_temp.get(), fragment_indices.data(), "HCD");

    masses.resize(count);
    qscores.resize(count);
    charges.resize(count);
    window_starts.resize(count);
    window_ends.resize(count);
    fragment_indices.resize(count);

    // Convert char array to int vector (1 for prefix ion, 0 for suffix ion)
    is_b_ions.resize(count);
    for (int i = 0; i < count; ++i)
    {
      char c = ion_types_temp[i];
      is_b_ions[i] = (c == 'a' || c == 'b' || c == 'c') ? 1 : 0;
    }

    return count;
  }

  int FLASHIda::getTerminalFragmentIons(
      const String& protein_sequence,
      int n,
      double* masses,
      double* qscores,
      int* charges,
      double* window_starts,
      double* window_ends,
      char* ion_types,
      int* fragment_indices,
      const String& fragmentation_method)
  {
    // Run fragment matching to get all matches
    std::vector<TagBasedFragmentMatch> matches;
    runTagBasedFragmentMatching_(protein_sequence, matches, nullptr, fragmentation_method);

    if (matches.empty()) return 0;

    // Separate into prefix ions (a, b, c) and suffix ions (x, y, z)
    std::vector<TagBasedFragmentMatch> prefix_ions, suffix_ions;
    for (const auto& m : matches)
    {
      if (isPrefixIon(m.ion_type)) prefix_ions.push_back(m);
      else suffix_ions.push_back(m);
    }

    std::cout << "[getTerminalFragmentIons] Found " << prefix_ions.size() << " prefix ions and "
              << suffix_ions.size() << " suffix ions" << std::endl;

    // Sort prefix ions by fragment_index descending (rightmost first), then qscore descending
    std::sort(prefix_ions.begin(), prefix_ions.end(), [](const auto& a, const auto& b) {
      if (a.fragment_index != b.fragment_index) return a.fragment_index > b.fragment_index;
      return a.qscore > b.qscore;
    });

    // Sort suffix ions by fragment_index descending (leftmost first), then qscore descending
    std::sort(suffix_ions.begin(), suffix_ions.end(), [](const auto& a, const auto& b) {
      if (a.fragment_index != b.fragment_index) return a.fragment_index > b.fragment_index;
      return a.qscore > b.qscore;
    });

    if (!prefix_ions.empty())
    {
      std::cout << "[getTerminalFragmentIons] Top prefix ions (rightmost): ";
      for (size_t i = 0; i < std::min(prefix_ions.size(), size_t(3)); ++i)
      {
        std::cout << prefix_ions[i].ion_type << prefix_ions[i].fragment_index << " ";
      }
      std::cout << std::endl;
    }
    if (!suffix_ions.empty())
    {
      std::cout << "[getTerminalFragmentIons] Top suffix ions (leftmost): ";
      for (size_t i = 0; i < std::min(suffix_ions.size(), size_t(3)); ++i)
      {
        std::cout << suffix_ions[i].ion_type << suffix_ions[i].fragment_index << " ";
      }
      std::cout << std::endl;
    }

    // Interleave output: prefix, suffix, prefix, suffix, ...
    int output_idx = 0;
    size_t idx_prefix = 0, idx_suffix = 0;
    bool select_prefix = true;  // Start with prefix

    while (output_idx < n && (idx_prefix < prefix_ions.size() || idx_suffix < suffix_ions.size()))
    {
      const TagBasedFragmentMatch* selected = nullptr;

      if (select_prefix)  // Prefer prefix ion
      {
        if (idx_prefix < prefix_ions.size()) {
          selected = &prefix_ions[idx_prefix++];
        } else if (idx_suffix < suffix_ions.size()) {
          selected = &suffix_ions[idx_suffix++];
        }
      }
      else  // Prefer suffix ion
      {
        if (idx_suffix < suffix_ions.size()) {
          selected = &suffix_ions[idx_suffix++];
        } else if (idx_prefix < prefix_ions.size()) {
          selected = &prefix_ions[idx_prefix++];
        }
      }

      if (selected)
      {
        const auto& pg = ms2_deconvolved_spectrum_[selected->peak_index];

        if (ms3_all_charges_)
        {
          // Collect all charges with non-zero intensity and sort by intensity descending
          auto [min_c, max_c] = pg.getAbsChargeRange();
          std::vector<std::pair<float, int>> charge_intensities;  // (intensity, charge)
          for (int c = min_c; c <= max_c; ++c)
          {
            float intensity = pg.getChargeIntensity(c);
            if (intensity > 0)
            {
              charge_intensities.emplace_back(intensity, c);
            }
          }
          std::sort(charge_intensities.begin(), charge_intensities.end(),
                    [](const auto& a, const auto& b) { return a.first > b.first; });

          for (const auto& [intensity, c] : charge_intensities)
          {
            if (output_idx >= n) break;
            masses[output_idx] = selected->observed_mass;
            qscores[output_idx] = intensity;
            charges[output_idx] = c;
            ion_types[output_idx] = selected->ion_type;
            fragment_indices[output_idx] = selected->fragment_index;
            auto [mz1, mz2] = pg.getMzRange(c);
            window_starts[output_idx] = mz1 - optimal_window_margin_;
            window_ends[output_idx] = mz2 + optimal_window_margin_;
            ++output_idx;
          }
          if (output_idx >= n) break;
        }
        else
        {
          // Single charge behavior
          masses[output_idx] = selected->observed_mass;
          qscores[output_idx] = selected->qscore;
          charges[output_idx] = selected->charge;
          ion_types[output_idx] = selected->ion_type;
          fragment_indices[output_idx] = selected->fragment_index;
          auto [mz_start, mz_end] = pg.getMzRange(selected->charge);
          window_starts[output_idx] = mz_start - optimal_window_margin_;
          window_ends[output_idx] = mz_end + optimal_window_margin_;
          ++output_idx;
        }
      }
      select_prefix = !select_prefix;  // Alternate between prefix and suffix
    }

    std::cout << "[getTerminalFragmentIons] Returning " << output_idx << " interleaved ions:" << std::endl;
    for (int i = 0; i < output_idx; ++i)
    {
      std::cout << "  " << ion_types[i] << fragment_indices[i] << " mass=" << masses[i]
                << " qscore=" << qscores[i] << " charge=" << charges[i] << std::endl;
    }

    return output_idx;
  }

  int FLASHIda::getTerminalFragmentIonsPy(
      const String& protein_sequence,
      int n,
      std::vector<double>& masses,
      std::vector<double>& qscores,
      std::vector<int>& charges,
      std::vector<double>& window_starts,
      std::vector<double>& window_ends,
      std::vector<int>& is_b_ions,
      std::vector<int>& fragment_indices)
  {
    masses.resize(n);
    qscores.resize(n);
    charges.resize(n);
    window_starts.resize(n);
    window_ends.resize(n);

    // Use raw char array for the C-style function
    std::unique_ptr<char[]> ion_types_temp(new char[n]);
    fragment_indices.resize(n);

    int count = getTerminalFragmentIons(
        protein_sequence, n,
        masses.data(), qscores.data(), charges.data(),
        window_starts.data(), window_ends.data(),
        ion_types_temp.get(), fragment_indices.data(), "HCD");

    // Resize to actual count
    masses.resize(count);
    qscores.resize(count);
    charges.resize(count);
    window_starts.resize(count);
    window_ends.resize(count);
    fragment_indices.resize(count);

    // Convert char array to int vector (1 for prefix ion, 0 for suffix ion)
    is_b_ions.resize(count);
    for (int i = 0; i < count; ++i)
    {
      char c = ion_types_temp[i];
      is_b_ions[i] = (c == 'a' || c == 'b' || c == 'c') ? 1 : 0;
    }

    return count;
  }

  double FLASHIda::getRepresentativeMass()
  {/*
    const int max_count = 10;
    double threshold = 0;
    double mass = 0;
    double intensity_sum = 0;

    if (deconvolved_spectrum_.size() > max_count)
    {
      std::vector<float> intensites;
      intensites.reserve(deconvolved_spectrum_.size());
      for (const auto& pg : deconvolved_spectrum_)
      {
        intensites.push_back(pg.getIntensity());
      }
      std::sort(intensites.rbegin(), intensites.rend());
      threshold = intensites[max_count];
    }

    for (const auto& pg : deconvolved_spectrum_)
    {
      if (pg.getIntensity() < threshold) continue;
      mass += pg.getMonoMass() * pg.getIntensity();
      intensity_sum += pg.getIntensity();
    }
    if (intensity_sum <= 0) return 0;
    return mass / intensity_sum;
    */
    auto filter_str = deconvolved_spectrum_.getOriginalSpectrum().getMetaValue("filter string").toString();
    Size pos = filter_str.find("cv=");
    double cv;

    if (pos != String::npos) // get the preferred mass ranges accding to CV values.
    {
      Size end = filter_str.find(" ", pos);
      if (end == String::npos) end = filter_str.length() - 1;
      cv = std::stod(filter_str.substr(pos + 3, end - pos));
      return cv;
    }
    return -100;
  }

  void FLASHIda::getIsolationWindows(double* wstart,
                                     double* wend,
                                     double* qscores,
                                     int* charges,
                                     int* min_charges,
                                     int* max_charges,
                                     double* mono_masses,
                                     double* chare_cos,
                                     double* charge_snrs,
                                     double* iso_cos,
                                     double* snrs,
                                     double* charge_scores,
                                     double* ppm_errors,
                                     double* precursor_intensities,
                                     double* peakgroup_intensities,
                                     int* hcds,
                                     int* ids)
  {
    // std::sort(selected_peak_groups_.begin(), selected_peak_groups_.end(), QscoreComparator_);

    for (Size i = 0; i < selected_peak_groups_.size(); i++)
    {
      if (trigger_charges[i] == 0) { continue; }
      auto peakgroup = selected_peak_groups_[i];
      charges[i] = trigger_charges[i];
      auto cr = peakgroup.getAbsChargeRange();
      min_charges[i] = std::get<0>(cr);
      max_charges[i] = std::get<1>(cr);

      wstart[i] = trigger_left_isolation_mzs_[i]; // std::get<0>(mz_range) - optimal_window_margin_;
      wend[i] = trigger_right_isolation_mzs_[i];  // std::get<1>(mz_range) + optimal_window_margin_;

      qscores[i] = trigger_scores[i];
      mono_masses[i] = peakgroup.getMonoMass();
      chare_cos[i] = peakgroup.getChargeIsotopeCosine(charges[i]);
      charge_snrs[i] = peakgroup.getChargeSNR(charges[i]);
      iso_cos[i] = peakgroup.getIsotopeCosine();
      snrs[i] = peakgroup.getSNR();
      charge_scores[i] = peakgroup.getChargeScore();
      ppm_errors[i] = peakgroup.getAvgPPMError();
      peakgroup_intensities[i] = peakgroup.getIntensity();
      precursor_intensities[i] = peakgroup.getChargeIntensity(charges[i]);
      hcds[i] = trigger_hcds[i];
      ids[i] = trigger_ids_[i];
    }
  }

  MSSpectrum FLASHIda::makeMSSpectrum_(const double* mzs, const double* ints, const int length, const double rt, const int ms_level, const char* name)
  {
    auto spec = MSSpectrum();
    for (int i = 0; i < length; i++)
    {
      if (ints[i] <= 0) { continue; }
      spec.emplace_back(mzs[i], ints[i]);
    }
    spec.setMSLevel(ms_level);
    spec.setName(name);
    spec.setRT(rt);
    return spec;
  }

  int FLASHIda::getConfigInt(const std::string& key) const
  {
    if (key == "targeting_mode") return targeting_mode_;
    if (key == "hcd_energy") return hcd_energy_;
    return 0;
  }

  double FLASHIda::getConfigDouble(const std::string& key) const
  {
    if (key == "rt_window") return rt_window_;
    return 0.0;
  }

  std::map<int, std::vector<std::vector<float>>> FLASHIda::parseFLASHIdaLog(const String& in_log_file)
  {
    std::map<int, std::vector<std::vector<float>>>
      precursor_map_for_real_time_acquisition; // ms1 scan -> mass, charge ,score, mz range, precursor int, mass int, color

    if (in_log_file.empty()) { return precursor_map_for_real_time_acquisition; }

    std::ifstream f(in_log_file.c_str());
    if (! f.good())
    {
      std::cout << "FLASHIda log file " << in_log_file << " is NOT found. FLASHIda support is not active." << std::endl;
      return precursor_map_for_real_time_acquisition;
    }


    std::cout << "FLASHIda log file used: " << in_log_file << std::endl;
    std::ifstream instream(in_log_file);
    if (instream.good())
    {
      String line;
      int scan;
      float mass, charge, w1, w2, qscore, pint, mint, z1, z2;
      float features[6];
      while (std::getline(instream, line))
      {
        if (line.find("0 targets") != line.npos) { continue; }
        if (line.hasPrefix("MS1"))
        {
          Size st = line.find("MS1 Scan# ") + 10;
          Size ed = line.find(' ', st);
          String n = line.substr(st, ed);
          scan = atoi(n.c_str());
          precursor_map_for_real_time_acquisition[scan]
            = std::vector<std::vector<float>>(); //// ms1 scan -> mass, charge ,score, mz range, precursor int, mass int, color
        }
        if (line.hasPrefix("Mass"))
        {
          Size st = 5;
          Size ed = line.find('\t');
          String n = line.substr(st, ed);
          mass = (float)atof(n.c_str());

          st = line.find("Z=") + 2;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          charge = (float)atof(n.c_str());

          st = line.find("Score=") + 6;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          qscore = (float)atof(n.c_str());

          st = line.find("[") + 1;
          ed = line.find('-', st);
          n = line.substr(st, ed);
          w1 = (float)atof(n.c_str());

          st = line.find('-', ed) + 1;
          ed = line.find(']', st);
          n = line.substr(st, ed);
          w2 = (float)atof(n.c_str());

          st = line.find("PrecursorIntensity=", ed) + 19;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          pint = (float)atof(n.c_str());

          st = line.find("PrecursorMassIntensity=", ed) + 23;
          ed = line.find('\t', st);
          n = line.substr(st, ed);
          mint = (float)atof(n.c_str());

          st = line.find("Features=", ed) + 9;
          // ed = line.find(' ', st);

          st = line.find('[', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[0] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[1] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[2] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[3] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(',', st);
          n = line.substr(st, ed);
          features[4] = (float)atof(n.c_str());

          st = line.find(',', st) + 1;
          ed = line.find(']', st);
          n = line.substr(st, ed);
          features[5] = (float)atof(n.c_str());

          st = line.find("ChargeRange=[", ed) + 13;
          ed = line.find('-', st);
          n = line.substr(st, ed);
          z1 = (float)atof(n.c_str());

          st = line.find("-", ed) + 1;
          ed = line.find(']', st);
          n = line.substr(st, ed);
          z2 = (float)atof(n.c_str());
          std::vector<float> e(15);
          e[0] = mass;
          e[1] = charge;
          e[2] = qscore;
          e[3] = w1;
          e[4] = w2;
          e[5] = pint;
          e[6] = mint;
          e[7] = z1;
          e[8] = z2;
          for (int i = 9; i < 15; i++)
          {
            e[i] = features[i - 9];
          }
          precursor_map_for_real_time_acquisition[scan].push_back(e);
        }
      }
      instream.close();
    }
    else { std::cout << in_log_file << " not found\n"; }
    int mass_cntr = 0;
    for (auto& v : precursor_map_for_real_time_acquisition)
    {
      std::sort(v.second.begin(), v.second.end(), [](const std::vector<float>& left, const std::vector<float>& right) { return left[0] < right[0]; });
      mass_cntr += v.second.size();
    }

    std::cout << "Used precursor size : " << precursor_map_for_real_time_acquisition.size() << " precursor masses : " << mass_cntr << std::endl;

    return precursor_map_for_real_time_acquisition;
  }

  void FLASHIda::parseInclusionListTSV_(const String& filename)
  {
    std::ifstream instream(filename);
    if (!instream.good())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    String line;
    bool header_skipped = false;

    while (std::getline(instream, line))
    {
      line = line.trim();
      if (line.empty() || line.hasPrefix("#")) continue;

      StringList cells;
      line.split('\t', cells);

      // Skip header row (check if first cell looks like "mass" header)
      if (!header_skipped && cells.size() > 0 && cells[0].toLower() == "mass")
      {
        header_skipped = true;
        continue;
      }

      // Expect 5 columns: mass, charge, rt_start, rt_end, priority
      if (cells.size() < 5)
      {
        OPENMS_LOG_WARN << "Skipping malformed line (expected 5 columns): " << line << std::endl;
        continue;
      }

      InclusionTarget target;
      target.mass = cells[0].toDouble();
      // Charge: -1 or empty means "any charge"
      String charge_str = cells[1].trim();
      target.charge = (charge_str.empty() || charge_str == "-1") ? -1 : cells[1].toInt();
      // RT values in minutes, convert to seconds
      target.rt_start = cells[2].toDouble() * 60.0;
      target.rt_end = cells[3].toDouble() * 60.0;
      target.priority = cells[4].toInt();

      inclusion_targets_.push_back(target);
    }

    // Sort by mass for efficient binary search later
    std::sort(inclusion_targets_.begin(), inclusion_targets_.end(),
      [](const InclusionTarget& a, const InclusionTarget& b) { return a.mass < b.mass; });
  }

  void FLASHIda::parseTargetPTMsTSV_(const String& filename)
  {
    std::ifstream instream(filename);
    if (!instream.good())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    String line;
    bool header_skipped = false;

    while (std::getline(instream, line))
    {
      line = line.trim();
      if (line.empty() || line.hasPrefix("#")) continue;

      StringList cells;
      line.split('\t', cells);

      // Skip header row (check if first cell looks like "name" header)
      if (!header_skipped && cells.size() > 0 && cells[0].toLower() == "name")
      {
        header_skipped = true;
        continue;
      }

      // Expect 3 columns: name, mass, max_count
      if (cells.size() < 3)
      {
        OPENMS_LOG_WARN << "Skipping malformed PTM line (expected 3 columns): " << line << std::endl;
        continue;
      }

      TargetPTM ptm;
      ptm.name = cells[0].trim();
      ptm.mass = cells[1].toDouble();
      ptm.max_count = cells[2].toInt();

      target_ptms_.push_back(ptm);
    }
  }

  std::vector<double> FLASHIda::generatePTMCombinations_(double base_mass,
                                                         const std::vector<TargetPTM>& ptms) const
  {
    std::vector<double> result;
    // Unmodified base mass intentionally excluded - only PTM-modified forms

    if (ptms.empty())
    {
      return result;
    }

    // Use iterative approach to generate all combinations
    // For each PTM, generate 0 to max_count occurrences, respecting global max
    std::function<void(Size, double, int)> generate;
    generate = [&](Size ptm_idx, double current_mass, int total_count) {
      if (ptm_idx >= ptms.size())
      {
        if (total_count > 0)  // Only include modified forms
        {
          result.push_back(current_mass);
        }
        return;
      }

      const TargetPTM& ptm = ptms[ptm_idx];
      // Generate both negative and positive mass shifts: -max_count to +max_count
      for (int count = -ptm.max_count; count <= ptm.max_count; ++count)
      {
        int new_total = total_count + std::abs(count);  // Track absolute count for pruning
        if (new_total > max_total_ptm_count_) continue;  // Skip if exceeds global max
        double new_mass = current_mass + count * ptm.mass;
        generate(ptm_idx + 1, new_mass, new_total);
      }
    };

    // Generate all combinations starting from base mass
    generate(0, base_mass, 0);

    // Remove duplicates (masses within 0.01 Da are considered equal)
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end(),
      [](double a, double b) { return std::abs(a - b) < 0.01; }), result.end());

    return result;
  }

  void FLASHIda::addDynamicTargets_(const std::vector<double>& masses,
                                    double rt,
                                    int priority)
  {
    // Add new targets to inclusion list
    // Note: target_masses_ and fd_.setTargetMasses() are managed by getPeakGroups()
    // which clears and rebuilds them from inclusion_targets_ on each MS1 scan.
    // We only need to add to inclusion_targets_ here.
    for (double mass : masses)
    {
      InclusionTarget target;
      target.mass = mass;
      target.charge = -1;  // Any charge
      target.rt_start = rt;
      target.rt_end = rt + rt_window_;
      target.priority = priority;

      inclusion_targets_.push_back(target);
    }

    // Re-sort inclusion_targets_ by mass for efficient lookup
    std::sort(inclusion_targets_.begin(), inclusion_targets_.end(),
      [](const InclusionTarget& a, const InclusionTarget& b) { return a.mass < b.mass; });

    std::cout << "Added " << masses.size() << " dynamic target masses (RT window: "
              << rt << "-" << (rt + rt_window_) << "s)\n";
    std::cout << "  Masses: ";
    for (Size i = 0; i < masses.size(); ++i)
    {
      if (i > 0) std::cout << ", ";
      std::cout << std::fixed << std::setprecision(4) << masses[i];
    }
    std::cout << std::endl;
  }

  bool FLASHIda::processMS2ForTagBasedTargeting(double precursor_mass)
  {
    // Early exit if tag-based targeting not enabled
    if (!tag_based_targeting_enabled_ || target_protein_database_.empty())
    {
      return false;
    }

    // Require deconvolveMS2() to be called first
    if (!ms2_deconv_valid_)
    {
      return false;
    }

    // Use stored MS2 deconvolution
    DeconvolvedSpectrum dspec = ms2_deconvolved_spectrum_;
    if (dspec.empty())
    {
      return false;
    }

    // Sort deconvolved spectrum by mass
    dspec.sort();

    // Create and configure tagger with tag length parameters
    FLASHTaggerAlgorithm tagger;
    Param tagger_param = tagger.getDefaults();
    tagger_param.setValue("min_length", min_tag_length_for_targeting_);
    tagger_param.setValue("max_length", max_tag_length_for_targeting_);
    tagger.setParameters(tagger_param);

    // Run tag generation
    tagger.run(dspec, tag_matching_tolerance_ppm_);

    // Get the generated tags
    std::vector<FLASHHelperClasses::Tag> tags;
    tagger.fillTags(tags);

    if (tags.empty())
    {
      return false;
    }

    // Prepare protein sequences for matching (replace I with L for matching)
    std::vector<String> protein_seqs;
    protein_seqs.reserve(target_protein_database_.size());
    for (const auto& fe : target_protein_database_)
    {
      String seq = fe.sequence;
      std::replace(seq.begin(), seq.end(), 'I', 'L');
      protein_seqs.push_back(seq);
    }

    // Match tags against target protein database
    std::set<Size> matched_protein_indices;
    for (const auto& tag : tags)
    {
      std::cout << tag.getSequence() << std::endl;
      for (Size protein_idx = 0; protein_idx < protein_seqs.size(); ++protein_idx)
      {
        const String& pseq = protein_seqs[protein_idx];

        std::vector<int> positions;
        std::vector<double> flanking_mass_diffs;

        FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(
          positions,
          flanking_mass_diffs,
          max_flanking_mass_diff_,
          pseq,
          tag);

        if (!positions.empty())
        {
          matched_protein_indices.insert(protein_idx);
        }
      }
    }

    if (matched_protein_indices.empty())
    {
      return false;
    }

    // Target protein detected! Expand target masses with PTM combinations
    // Use precursor mass from iAPI as base (not theoretical protein mass)
    std::vector<double> new_targets;
    if (precursor_mass > 0)
    {
      std::vector<double> ptm_masses = generatePTMCombinations_(precursor_mass, target_ptms_);
      new_targets.insert(new_targets.end(), ptm_masses.begin(), ptm_masses.end());
    }

    // Remove duplicates and already-tracked masses
    std::sort(new_targets.begin(), new_targets.end());
    new_targets.erase(std::unique(new_targets.begin(), new_targets.end(),
      [](double a, double b) { return std::abs(a - b) < 0.01; }), new_targets.end());

    // Filter out masses already in expanded set
    std::vector<double> truly_new_targets;
    for (double mass : new_targets)
    {
      int nominal = SpectralDeconvolution::getNominalMass(mass);
      bool already_tracked = false;
      for (double tracked : expanded_target_masses_)
      {
        if (std::abs(SpectralDeconvolution::getNominalMass(tracked) - nominal) <= 1)
        {
          already_tracked = true;
          break;
        }
      }
      if (!already_tracked)
      {
        truly_new_targets.push_back(mass);
        expanded_target_masses_.insert(mass);
      }
    }

    if (truly_new_targets.empty())
    {
      return true; // Target protein was detected, but all masses already tracked
    }

    // Add to dynamic inclusion list
    addDynamicTargets_(truly_new_targets, ms2_deconv_rt_, 100); // High priority

    return true;
  }

  void FLASHIda::parseJSONConfig_(const std::string& json_str)
  {
    try
    {
      using json = nlohmann::json;
      json config = json::parse(json_str);

      // --- deconvolution section ---
      auto deconv = config.value("deconvolution", json::object());
      qscore_threshold_ = deconv.value("score_threshold", 0.0);
      tqscore_threshold = deconv.value("tqscore_threshold", 0.9);
      int min_charge = deconv.value("min_charge", 4);
      int max_charge = deconv.value("max_charge", 50);
      double min_mass = deconv.value("min_mass", 500.0);
      double max_mass = deconv.value("max_mass", 50000.0);

      DoubleList tol_values;
      if (deconv.contains("tol") && deconv["tol"].is_array())
      {
        for (const auto& v : deconv["tol"])
          tol_values.push_back(v.get<double>());
      }
      if (tol_values.empty())
        tol_values = {10.0, 10.0};
      if (tol_values.size() == 1)
        tol_values.push_back(tol_values[0]);
      tol_ = std::vector<double>(tol_values);

      // Use MS2 tolerance for tag matching
      if (tol_.size() >= 2)
        tag_matching_tolerance_ppm_ = tol_[1];
      else if (tol_.size() == 1)
        tag_matching_tolerance_ppm_ = tol_[0];

      // --- precursor_selection section ---
      auto ps = config.value("precursor_selection", json::object());
      auto mass_count_arr = ps.value("max_mass_count", std::vector<int>{1});
      for (int j : mass_count_arr)
        mass_count_.push_back(j);

      rt_window_ = ps.value("RT_window", 180.0);
      targeting_mode_ = ps.value("target_mode", 0);
      use_idscore_ = ps.value("IDScore", false);
      consider_all_Charge_states_ = ps.value("AllCharges", false);
      ms3_all_charges_ = ps.value("MS3AllCharges", false);
      hcd_energy_ = ps.value("HCDEnergy", -1);
      strict_inclusion_ = ps.value("strict_inclusion", false);
      tie_threshold_ = ps.value("tie_threshold", 0.1);

      if (targeting_mode_ == 1)
        std::cout << "Inclusion mode: " << (strict_inclusion_ ? "strict" : "non-strict") << "\n";

      // --- tagging section ---
      auto tagging = config.value("tagging", json::object());
      min_tag_length_for_targeting_ = tagging.value("min_tag_length", 3);
      max_tag_length_for_targeting_ = tagging.value("max_tag_length", 8);
      max_total_ptm_count_ = tagging.value("max_ptm_count", 3);
      max_flanking_mass_diff_ = tagging.value("max_flanking_mass_diff", 50000.0);

      // --- files section ---
      auto files = config.value("files", json::object());

      // Target log files
      std::vector<String> log_files;
      if (files.contains("target_logs") && files["target_logs"].is_array())
      {
        for (const auto& f : files["target_logs"])
          log_files.push_back(f.get<std::string>());
      }

      // FASTA files
      std::vector<String> fasta_files;
      if (files.contains("fasta") && !files["fasta"].get<std::string>().empty())
        fasta_files.push_back(files["fasta"].get<std::string>());

      // TSV inclusion list
      std::vector<String> tsv_files;
      if (files.contains("inclusion_list") && !files["inclusion_list"].get<std::string>().empty())
        tsv_files.push_back(files["inclusion_list"].get<std::string>());

      // PTM list
      std::vector<String> ptm_files;
      if (files.contains("ptm_list") && !files["ptm_list"].get<std::string>().empty())
        ptm_files.push_back(files["ptm_list"].get<std::string>());

      // --- Build SpectralDeconvolution Param (must match legacy path exactly) ---
      Param sd_defaults = SpectralDeconvolution().getDefaults();
      sd_defaults.setValue("min_charge", min_charge);
      sd_defaults.setValue("max_charge", max_charge);
      sd_defaults.setValue("min_mass", min_mass);
      sd_defaults.setValue("max_mass", max_mass);
      sd_defaults.setValue("tol", tol_values);

      // --- Load target log files (same as legacy path) ---
      std::stringstream ss{};
      for (const auto& log_file : log_files)
      {
        ss << log_file << " ";
        std::ifstream instream(log_file);
        if (instream.good())
        {
          String line;
          double rt = .0, mass, qscore;
          int charge;
          while (std::getline(instream, line))
          {
            if (line.find("0 targets") != line.npos)
              continue;
            if (line.hasPrefix("MS1"))
            {
              Size st = line.find("RT ") + 3;
              Size ed = line.find('(') - 2;
              String n = line.substr(st, ed - st + 1);
              rt = atof(n.c_str());
            }
            if (line.hasPrefix("Mass"))
            {
              Size st = 5;
              Size ed = line.find('\t');
              String n = line.substr(st, ed - st + 1);
              mass = atof(n.c_str());

              st = line.find("Score=") + 6;
              ed = line.find('\t', st);
              n = line.substr(st, ed - st + 1);
              qscore = atof(n.c_str());

              st = line.find("Z=") + 2;
              ed = line.find('\t', st);
              n = line.substr(st, ed - st + 1);
              charge = n.toInt();

              if (targeting_mode_ == 1 || targeting_mode_ == 2)
              {
                if (target_mass_rt_map_.find(mass) == target_mass_rt_map_.end())
                  target_mass_rt_map_[mass] = std::vector<double>();
                target_mass_rt_map_[mass].push_back(rt * 60.0);
                if (target_mass_qscore_map_.find(mass) == target_mass_qscore_map_.end())
                  target_mass_qscore_map_[mass] = std::vector<double>();
                target_mass_qscore_map_[mass].push_back(qscore);
                if (target_mass_charge_map_.find(mass) == target_mass_charge_map_.end())
                  target_mass_charge_map_[mass] = std::vector<int>();
                target_mass_charge_map_[mass].push_back(charge);
              }
            }
            else if (line.hasPrefix("AllMass"))
            {
              if (targeting_mode_ == 3)
              {
                Size st = 8;
                Size ed = line.size();
                String n = line.substr(st, ed - st + 1);
                std::stringstream tmp_stream(n);
                String str;
                std::vector<double> results;
                while (getline(tmp_stream, str, ' '))
                  results.push_back(atof(str.c_str()));
                if (exclusion_rt_masses_map_.find(rt * 60.0) == exclusion_rt_masses_map_.end())
                  exclusion_rt_masses_map_[rt * 60.0] = std::vector<double>();
                for (double m : results)
                  exclusion_rt_masses_map_[rt * 60.0].push_back(m);
              }
            }
          }
          instream.close();
        }
      }

      if (targeting_mode_ == 1)
        std::cout << ss.str() << "file(s) is(are) used for inclusion mode\n";
      else if (targeting_mode_ == 2)
        std::cout << ss.str() << "file(s) is(are) used for in-depth mode\n";
      else if (targeting_mode_ == 3)
        std::cout << ss.str() << "file(s) is(are) used for exclusion mode\n";

      // Parse TSV inclusion list files
      if (targeting_mode_ == 1 && !tsv_files.empty())
      {
        for (const auto& tsv_file : tsv_files)
          parseInclusionListTSV_(tsv_file);
        std::cout << inclusion_targets_.size() << " targets loaded from TSV inclusion list\n";
      }

      // Load FASTA files for tag-based targeting
      if (!fasta_files.empty())
      {
        for (const auto& fasta_file : fasta_files)
        {
          std::vector<FASTAFile::FASTAEntry> entries;
          FASTAFile().load(fasta_file, entries);
          target_protein_database_.insert(target_protein_database_.end(), entries.begin(), entries.end());
        }
        std::cout << target_protein_database_.size() << " protein entries loaded for tag-based targeting\n";
        tag_based_targeting_enabled_ = true;
      }

      if (tag_based_targeting_enabled_)
      {
        std::cout << "Tag-based targeting: min_tag_length=" << min_tag_length_for_targeting_
                  << ", max_tag_length=" << max_tag_length_for_targeting_
                  << ", tolerance=" << tag_matching_tolerance_ppm_ << " ppm"
                  << ", max_flanking_mass_diff=" << max_flanking_mass_diff_ << " Da\n";
      }

      // Load PTM TSV files
      for (const auto& ptm_file : ptm_files)
        parseTargetPTMsTSV_(ptm_file);
      if (!target_ptms_.empty())
      {
        std::cout << target_ptms_.size() << " PTM modifications loaded for target expansion (max "
                  << max_total_ptm_count_ << " total per proteoform)\n";
      }

      // --- Store new JSON-only sections for future phases ---

      // ms3
      auto ms3 = config.value("ms3", json::object());
      ms3_enabled_ = ms3.value("enabled", false);
      ms3_mode_ = ms3.value("mode", 0);
      max_ms3_per_ms2_ = ms3.value("max_per_ms2", 4);
      ms3_protein_sequence_ = ms3.value("protein_sequence", "");

      // conditional ms2
      conditional_ms2_enabled_ = config.value("conditional_ms2", false);

      // quantification
      auto quant = config.value("quantification", json::object());
      quant_enabled_ = quant.value("enabled", false);
      reporter_mz_tol_ = quant.value("reporter_mz_tol", 0.002);
      fold_change_threshold_ = quant.value("fold_change_threshold", 1.4);

      // faims
      auto faims = config.value("faims", json::object());
      if (faims.contains("cv_values") && faims["cv_values"].is_array())
      {
        for (const auto& v : faims["cv_values"])
          faims_cv_values_.push_back(v.get<double>());
      }
      max_cv_skip_ = faims.value("max_cv_skip", 0);
      cv_precursor_threshold_ = faims.value("cv_precursor_threshold", 15);

      // Initialize per-CV adaptive skip state
      faims_enabled_ = (faims_cv_values_.size() > 1);
      if (faims_enabled_)
      {
        int n = static_cast<int>(faims_cv_values_.size());
        cv_skip_amount_.resize(n, 0);
        cv_skip_count_.resize(n, 0);
        current_cv_index_ = 0;
      }

      // ms_settings
      auto ms_settings = config.value("ms_settings", json::object());
      auto ms1 = ms_settings.value("ms1", json::object());
      ms1_analyzer_ = ms1.value("analyzer", "");
      ms1_first_mass_ = ms1.value("first_mass", 0.0);
      ms1_last_mass_ = ms1.value("last_mass", 0.0);
      ms1_resolution_ = ms1.value("resolution", 0);
      ms1_agc_target_ = ms1.value("agc_target", 0);
      ms1_max_it_ = ms1.value("max_it", 0.0);

      ms2_configs_.clear();
      if (ms_settings.contains("ms2") && ms_settings["ms2"].is_array())
      {
        for (const auto& m : ms_settings["ms2"])
        {
          MS2ConfigJson cfg;
          cfg.analyzer = m.value("analyzer", "");
          cfg.activation = m.value("activation", "");
          cfg.collision_energy = m.value("collision_energy", 0);
          cfg.resolution = m.value("resolution", 0);
          ms2_configs_.push_back(cfg);
        }
      }

      // scheduling
      auto sched = config.value("scheduling", json::object());
      auto ct = sched.value("cycle_time", json::object());
      cycle_time_enabled_ = ct.value("enabled", false);
      cycle_time_ms_ = ct.value("value_ms", 60000.0);
      auto to = sched.value("scan_timeout", json::object());
      timeout_enabled_ = to.value("enabled", false);
      timeout_ms_ = to.value("value_ms", 30000.0);
      double agc_interval_sec = sched.value("agc_interval_seconds", 30.0);
      agc_interval_ms_ = static_cast<uint64_t>(agc_interval_sec * 1000.0);
      last_agc_time_ = std::chrono::steady_clock::now();

      // --- Phase 7: selection_strategy (required) ---
      if (!config.contains("selection_strategy"))
      {
        throw std::runtime_error("FLASHIda: missing required 'selection_strategy' in JSON config");
      }
      const auto& sel_strategy = config["selection_strategy"];
      for (auto it = sel_strategy.begin(); it != sel_strategy.end(); ++it)
      {
        std::string ms_key = it.key();
        if (ms_key.substr(0, 2) == "ms" && ms_key.size() > 2)
        {
          int level = std::stoi(ms_key.substr(2));
          const auto& level_obj = it.value();
          MSLevelConfig cfg;
          // Selection metric
          std::string sel_str = level_obj.value("selection", level == 1 ? std::string("qscore") : std::string("intensity"));
          if (sel_str == "intensity") cfg.selection = SelectionMetric::Intensity;
          else if (sel_str == "qscore") cfg.selection = SelectionMetric::QScore;
          else if (sel_str == "none") cfg.selection = SelectionMetric::None;
          else cfg.selection = SelectionMetric::Intensity;
          // Max targets (with aliases)
          cfg.max_targets = level_obj.value("max_targets",
              level_obj.value("max_precursors",
              level_obj.value("max_fragments", 10)));
          // Exploration (optional, MS2+ only; guard against JSON null)
          if (level_obj.contains("exploration") && !level_obj["exploration"].is_null() && level > 1)
          {
            const auto& expl_obj = level_obj["exploration"];
            std::string met_str = expl_obj.value("metric", std::string("none"));
            if (met_str == "mass_count") cfg.exploration = ExplorationMetric::MassCount;
            else if (met_str == "remaining_precursor") cfg.exploration = ExplorationMetric::RemainingPrecursor;
            else if (met_str == "fragment_count") cfg.exploration = ExplorationMetric::FragmentCount;
            else cfg.exploration = ExplorationMetric::None;
            cfg.ce_min = expl_obj.value("ce_min", 20.0);
            cfg.ce_max = expl_obj.value("ce_max", 40.0);
            cfg.ce_step = expl_obj.value("ce_step", 5.0);
            cfg.activation = expl_obj.value("activation", std::string("HCD"));
            if (expl_obj.contains("overrides"))
            {
              const auto& ov_obj = expl_obj["overrides"];
              for (auto ov_it = ov_obj.begin(); ov_it != ov_obj.end(); ++ov_it)
                cfg.overrides[ov_it.key()] = ov_it.value().get<std::string>();
            }
          }
          level_configs_[level] = cfg;
        }
      }
      // Compute convenience boolean
      exploration_enabled_ = false;
      for (auto lc_it = level_configs_.begin(); lc_it != level_configs_.end(); ++lc_it)
      {
        if (hasExploration(lc_it->second)) { exploration_enabled_ = true; break; }
      }

      // --- Initialize SpectralDeconvolution (must match legacy path) ---
      snr_threshold_ = 1;
      fd_.setParameters(sd_defaults);
      fd_.calculateAveragine(false);
    }
    catch (const nlohmann::json::exception& e)
    {
      std::cerr << "FLASHIda JSON config parse error: " << e.what() << std::endl;
    }
  }

  // --- Phase 3: Scan command queue and tracking ---

  std::string FLASHIda::encodeBase36_(int value)
  {
    static const char digits[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    char buf[5] = {'0', '0', '0', '0', '\0'};
    for (int i = 3; i >= 0; --i)
    {
      buf[i] = digits[value % 36];
      value /= 36;
    }
    return std::string(buf);
  }

  int FLASHIda::nextTrackingIdInt_()
  {
    int id = tracking_id_counter_++;
    // Wrap at 36^4 - 1 = 1679615 to stay within 4-char base-36 range
    if (tracking_id_counter_ > 1679615)
      tracking_id_counter_ = 0;
    return id;
  }

  bool FLASHIda::needsAGCScan_() const
  {
    auto now = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - last_agc_time_).count();
    return static_cast<uint64_t>(elapsed) > agc_interval_ms_;
  }

  uint64_t FLASHIda::msSinceLastMS1_() const
  {
    auto now = std::chrono::steady_clock::now();
    return static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(now - last_ms1_time_).count());
  }

  ScanCommand FLASHIda::makeMS1Command_() const
  {
    ScanCommand cmd{};
    cmd.msn_level = 1;
    cmd.priority = 3; // lowest priority — MS1 is the fallback
    cmd.is_agc = 0;
    cmd.num_stages = 0;
    cmd.orbitrap_resolution = ms1_resolution_;
    cmd.agc_target = ms1_agc_target_;
    cmd.first_mass = ms1_first_mass_;
    cmd.last_mass = ms1_last_mass_;
    cmd.max_it = ms1_max_it_;

    // Copy analyzer string safely
    std::strncpy(cmd.analyzer, ms1_analyzer_.c_str(), sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';

    std::strncpy(cmd.scan_description, "MS1 survey scan", sizeof(cmd.scan_description) - 1);
    cmd.scan_description[sizeof(cmd.scan_description) - 1] = '\0';

    return cmd;
  }

  ScanCommand FLASHIda::makeAGCCommand_() const
  {
    ScanCommand cmd{};
    cmd.msn_level = 1;
    cmd.priority = 0; // highest priority — AGC is time-critical
    cmd.is_agc = 1;
    cmd.num_stages = 0;
    cmd.orbitrap_resolution = 0;
    cmd.agc_target = ms1_agc_target_;
    cmd.first_mass = ms1_first_mass_;
    cmd.last_mass = ms1_last_mass_;
    cmd.max_it = ms1_max_it_;
    std::strncpy(cmd.analyzer, "IonTrap", sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    std::strncpy(cmd.scan_description, "AGC calibration", sizeof(cmd.scan_description) - 1);
    cmd.scan_description[sizeof(cmd.scan_description) - 1] = '\0';

    return cmd;
  }

  void FLASHIda::cleanupExpiredCommands_()
  {
    if (!timeout_enabled_)
      return;

    auto now_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    auto it = pending_scan_map_.begin();
    while (it != pending_scan_map_.end())
    {
      if (it->second.enqueue_timestamp_ms > 0 &&
          (now_ms - it->second.enqueue_timestamp_ms) > static_cast<uint64_t>(timeout_ms_))
      {
        std::string id_str = encodeBase36_(it->first);
        std::cout << "[TRACK-EXPIRE] id=" << id_str
                  << " age_ms=" << (now_ms - it->second.enqueue_timestamp_ms)
                  << std::endl;
        it = pending_scan_map_.erase(it);
      }
      else
      {
        ++it;
      }
    }
  }

  // ---- Phase 6: FAIMS CV adaptive skip state machine ----
  //
  // Behavioral audit of C# ScanScheduler (reference implementation):
  //
  // Q1. CV advance trigger: updateCV() called after every MS1 deconvolution with
  //     the CV of that scan and precursor count. getFAIMSMS1Scan() does actual cycling.
  // Q2. Threshold comparison: precursors < MassThreshold (strictly less than). Default=15.
  // Q3. Skip counter: CVSkipAmount[pos] doubles (min=1, cap=MaxCVSkip).
  //     CVSkipCount[pos] resets to 0 in BOTH branches of updateCV().
  //     In getFAIMSMS1Scan(), CVSkipCount[i]++ when skipping.
  // Q4. Forced cycle: when CVSkipAmount >= MaxCVSkip, skip amount stops growing.
  //     CV is visited when CVSkipCount >= CVSkipAmount (counter exhausted).
  // Q5. CV transition: getFAIMSMS1Scan() enqueues AGC+MS1 for next CV.
  // Q6. Call order: deconvolve → create MS2s → updateCV(cv, precursors) → getFAIMSMS1Scan(true).
  // Q7. Cycling: forward only (currentCV++), wraps at end. Increment-first (index 0 skipped on first call).
  // Q8. Single CV: Flash.cs uses UnifiedScanProcessor, ScanScheduler with UseFAIMS=false.
  // Q9. Non-FAIMS: faims_cv = 0.0 on all commands.

  void FLASHIda::updateCVSkip_(double cv, int precursor_count)
  {
    if (!faims_enabled_) return;

    // Find position for this CV value
    int pos = -1;
    for (int i = 0; i < static_cast<int>(faims_cv_values_.size()); ++i)
    {
      if (faims_cv_values_[i] == cv) { pos = i; break; }
    }
    if (pos < 0) return;  // unknown CV

    if (precursor_count < cv_precursor_threshold_)  // strictly < (audit Q2)
    {
      if (cv_skip_amount_[pos] < max_cv_skip_)
      {
        cv_skip_amount_[pos] *= 2;                  // double spacing (audit Q3)
        if (cv_skip_amount_[pos] <= 0)
          cv_skip_amount_[pos] = 1;                 // min = 1
        if (cv_skip_amount_[pos] > max_cv_skip_)
          cv_skip_amount_[pos] = max_cv_skip_;      // cap at max
      }
      cv_skip_count_[pos] = 0;                      // reset in BOTH branches (audit Q3)
    }
    else
    {
      cv_skip_amount_[pos] = 0;                     // high precursor count: reset
      cv_skip_count_[pos] = 0;
    }
  }

  double FLASHIda::advanceToNextCV_()
  {
    int n = static_cast<int>(faims_cv_values_.size());
    // Safety bound (C# uses while(true); we bound to n iterations)
    for (int attempts = 0; attempts < n; ++attempts)
    {
      current_cv_index_++;                           // increment-first (audit Q7)
      if (current_cv_index_ >= n)
        current_cv_index_ = 0;                       // wrap (audit Q7)

      if (cv_skip_count_[current_cv_index_] < cv_skip_amount_[current_cv_index_])
      {
        cv_skip_count_[current_cv_index_]++;         // skip this CV (audit Q3)
        OPENMS_LOG_DEBUG << "[FAIMS] Skipping CV=" << faims_cv_values_[current_cv_index_]
                         << " (" << cv_skip_count_[current_cv_index_]
                         << "/" << cv_skip_amount_[current_cv_index_] << ")" << std::endl;
      }
      else
      {
        OPENMS_LOG_DEBUG << "[FAIMS] Changed to CV=" << faims_cv_values_[current_cv_index_] << std::endl;
        return faims_cv_values_[current_cv_index_];  // use this CV
      }
    }
    // Fallback: all CVs being skipped — use current anyway
    return faims_cv_values_[current_cv_index_];
  }

  int FLASHIda::processScan(const double* mzs, const double* ints, int length,
                             double rt_min, int ms_level, const char* scan_description,
                             double faims_cv)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);

    if (ms_level == 1)
    {
      last_ms1_time_ = std::chrono::steady_clock::now();

      // MS1 path: deconvolve, score, filter, select top-N, push MS2 commands
      double parent_cv = faims_enabled_ ? faims_cv : 0.0;

      int n = getPeakGroups(mzs, ints, length, rt_min, 1, "ms1", nullptr);
      int commands_pushed = 0;
      for (int i = 0; i < n; i++)
      {
        ScanCommand cmd = buildMS2Command_(selected_peak_groups_[i],
                                           trigger_charges[i], trigger_hcds[i]);
        cmd.faims_cv = parent_cv;  // MS2 carries parent MS1's CV
        pushCommand_(cmd);
        commands_pushed++;
      }

      // Phase 7: Initiate exploration for selected precursors if MS2 exploration is enabled
      if (hasExploration(getLevelConfig_(2)))
      {
        for (int i = 0; i < n; i++)
        {
          auto [mz1, mz2] = selected_peak_groups_[i].getMzRange(trigger_charges[i]);
          double center_mz = (mz1 + mz2) / 2.0;
          initiateExploration_(2, center_mz,
              selected_peak_groups_[i].getMonoMass(),
              trigger_charges[i], parent_cv);
        }
      }

      // FAIMS CV cycling: update skip policy, advance to next CV, push MS1
      if (faims_enabled_)
      {
        double current_cv = faims_cv_values_[current_cv_index_];
        updateCVSkip_(current_cv, commands_pushed);

        double next_cv = advanceToNextCV_();
        ScanCommand ms1 = makeMS1Command_();
        ms1.faims_cv = next_cv;
        ms1.scan_id = nextTrackingIdInt_();
        ms1.priority = 0;  // priority 0 to send before pending MS2s

        std::string id_str = encodeBase36_(ms1.scan_id);
        std::snprintf(ms1.scan_description, sizeof(ms1.scan_description),
                      "%s|CV transition MS1 CV=%.1f", id_str.c_str(), next_cv);

        auto now = std::chrono::steady_clock::now();
        ms1.enqueue_timestamp_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()).count();

        pushCommand_(ms1);
        std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=cv_transition cv=" << next_cv << std::endl;
      }

      return commands_pushed;
    }
    else if (ms_level == 2)
    {
      return processMS2Path_(mzs, ints, length, rt_min, scan_description);
    }
    return 0;
  }

  int FLASHIda::decodeBase36_(const std::string& s) const
  {
    int result = 0;
    for (char c : s)
    {
      result *= 36;
      if (c >= '0' && c <= '9')
        result += c - '0';
      else if (c >= 'a' && c <= 'z')
        result += 10 + (c - 'a');
    }
    return result;
  }

  ScanCommand FLASHIda::buildMS2Command_(const PeakGroup& pg, int charge, int hcd)
  {
    ScanCommand cmd{};
    int id = nextTrackingIdInt_();
    cmd.scan_id = id;
    cmd.msn_level = 2;
    cmd.priority = 1;
    cmd.is_agc = 0;
    cmd.num_stages = 1;

    // Use first MS2 config for resolution and analyzer
    if (!ms2_configs_.empty())
    {
      cmd.orbitrap_resolution = ms2_configs_[0].resolution;
      std::strncpy(cmd.analyzer, ms2_configs_[0].analyzer.c_str(), sizeof(cmd.analyzer) - 1);
      cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    }
    else
    {
      cmd.orbitrap_resolution = ms1_resolution_;
      std::strncpy(cmd.analyzer, ms1_analyzer_.c_str(), sizeof(cmd.analyzer) - 1);
      cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    }

    cmd.agc_target = ms1_agc_target_;
    cmd.first_mass = ms1_first_mass_;
    cmd.last_mass = ms1_last_mass_;
    cmd.max_it = ms1_max_it_;

    // Isolation window from peak group m/z range
    auto [mz1, mz2] = pg.getMzRange(charge);
    double center_mz = (mz1 + mz2) / 2.0;
    mz1 -= optimal_window_margin_;
    mz2 += optimal_window_margin_;
    double iso_width = mz2 - mz1;

    cmd.stages[0].precursor_mz = center_mz;
    cmd.stages[0].isolation_width = iso_width;
    cmd.stages[0].charge_state = charge;

    // Activation type and collision energy
    std::string activation = "HCD";
    int ce = 0;
    if (!ms2_configs_.empty())
    {
      activation = ms2_configs_[0].activation;
      ce = ms2_configs_[0].collision_energy;
    }
    // For HCD activation, use hcd as fallback only when config CE is unset
    if (activation == "HCD" && ce <= 0)
    {
      ce = (hcd > 0) ? hcd : 29;
    }
    cmd.stages[0].collision_energy = static_cast<double>(ce);
    std::strncpy(cmd.stages[0].activation_type, activation.c_str(), sizeof(cmd.stages[0].activation_type) - 1);
    cmd.stages[0].activation_type[sizeof(cmd.stages[0].activation_type) - 1] = '\0';

    // Scan description with tracking ID — base-36 format: {base36}|{mass:.2f}@{charge}
    std::string id_str = encodeBase36_(id);
    char desc_buf[256];
    std::snprintf(desc_buf, sizeof(desc_buf), "%s|%.2f@%d", id_str.c_str(), pg.getMonoMass(), charge);
    std::strncpy(cmd.scan_description, desc_buf, sizeof(cmd.scan_description) - 1);
    cmd.scan_description[sizeof(cmd.scan_description) - 1] = '\0';

    // Timestamp
    cmd.enqueue_timestamp_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    // Precursor scoring data for diagnostic TSV output
    cmd.qscore = pg.getQscore();
    cmd.mono_mass = pg.getMonoMass();
    cmd.charge_cos = pg.getChargeIsotopeCosine(std::abs(charge));
    cmd.charge_snr = pg.getChargeSNR(std::abs(charge));
    cmd.iso_cos = pg.getIsotopeCosine();
    cmd.snr = pg.getSNR();
    cmd.charge_score = pg.getChargeScore();
    cmd.ppm_error = pg.getAvgPPMError();
    cmd.precursor_intensity = pg.getChargeIntensity(std::abs(charge));
    cmd.peakgroup_intensity = pg.getIntensity();
    cmd.hcd_energy = (hcd > 0) ? hcd : ce;
    cmd.pad2 = 0;

    // Store in pending map for MS2 tracking resolution
    pending_scan_map_[id] = cmd;

    std::cout << "[TRACK-CREATE] id=" << id_str
              << " ms_level=2"
              << " mz=" << center_mz
              << " z=" << charge
              << " mass=" << pg.getMonoMass()
              << std::endl;

    return cmd;
  }

  ScanCommand FLASHIda::buildMS3Command_(const ScanCommand& ms2_ctx, double frag_mz, int frag_charge, double iso_width,
                                          char ion_type, int frag_index)
  {
    ScanCommand cmd{};
    int id = nextTrackingIdInt_();
    cmd.scan_id = id;
    cmd.msn_level = 3;
    cmd.priority = 3;  // lowest priority for MS3
    cmd.is_agc = 0;
    cmd.num_stages = 2;

    // Copy analyzer/resolution from MS2 context
    cmd.orbitrap_resolution = ms2_ctx.orbitrap_resolution;
    std::strncpy(cmd.analyzer, ms2_ctx.analyzer, sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    cmd.agc_target = ms2_ctx.agc_target;
    cmd.first_mass = ms2_ctx.first_mass;
    cmd.last_mass = ms2_ctx.last_mass;
    cmd.max_it = ms2_ctx.max_it;

    // Stage 0: MS2 precursor (from MS2 context)
    cmd.stages[0] = ms2_ctx.stages[0];

    // Stage 1: Fragment target
    cmd.stages[1].precursor_mz = frag_mz;
    cmd.stages[1].isolation_width = iso_width;
    cmd.stages[1].charge_state = frag_charge;
    cmd.stages[1].collision_energy = cmd.stages[0].collision_energy;
    std::strncpy(cmd.stages[1].activation_type, "HCD", sizeof(cmd.stages[1].activation_type) - 1);
    cmd.stages[1].activation_type[sizeof(cmd.stages[1].activation_type) - 1] = '\0';

    // Description — base-36 format with optional ion annotation (modes 3/4)
    std::string id_str = encodeBase36_(id);
    char desc_buf[256];
    if (ion_type != '\0' && frag_index > 0)
      std::snprintf(desc_buf, sizeof(desc_buf), "%s|MS3 %c%d mz=%.2f z=%d",
                    id_str.c_str(), ion_type, frag_index, frag_mz, frag_charge);
    else
      std::snprintf(desc_buf, sizeof(desc_buf), "%s|MS3 mz=%.2f z=%d",
                    id_str.c_str(), frag_mz, frag_charge);
    std::strncpy(cmd.scan_description, desc_buf, sizeof(cmd.scan_description) - 1);
    cmd.scan_description[sizeof(cmd.scan_description) - 1] = '\0';

    cmd.enqueue_timestamp_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    // D6: Store in pending map for future MS3 resolution
    pending_scan_map_[id] = cmd;

    std::cout << "[TRACK-CREATE] id=" << id_str
              << " ms_level=3"
              << " frag_mz=" << frag_mz
              << std::endl;

    return cmd;
  }

  void FLASHIda::pushCommand_(ScanCommand cmd)
  {
    int p = std::clamp(cmd.priority, 0, 3);
    queues_[p].push_back(cmd);
  }

  std::vector<FLASHIda::MS3Target> FLASHIda::selectMS3Targets_()
  {
    std::vector<MS3Target> targets;
    if (!ms3_enabled_ || ms3_mode_ == 0 || !ms2_deconv_valid_)
      return targets;

    const int n = max_ms3_per_ms2_;
    std::vector<double> masses(n), qscores(n), wstarts(n), wends(n);
    std::vector<int> charges(n);
    std::vector<char> ion_types(n, '\0');
    std::vector<int> frag_indices(n, 0);

    int count = 0;

    if (ms3_mode_ == 1 || ms3_mode_ == 2)
    {
      // Mode 1 (SourceCID) and Mode 2 (SPS): Use getBestMS2Masses
      count = getBestMS2Masses(n, masses.data(), qscores.data(), charges.data(),
                               wstarts.data(), wends.data());
    }
    else if (ms3_mode_ == 3 && !ms3_protein_sequence_.empty())
    {
      // Mode 3 (HCD-triggered): Use getTopFragmentMatches
      count = getTopFragmentMatches(ms3_protein_sequence_, n, masses.data(), qscores.data(),
                                    charges.data(), wstarts.data(), wends.data(),
                                    ion_types.data(), frag_indices.data(), "HCD");
    }
    else if (ms3_mode_ == 4 && !ms3_protein_sequence_.empty())
    {
      // Mode 4 (EThcD-triggered): Use getTerminalFragmentIons
      count = getTerminalFragmentIons(ms3_protein_sequence_, n, masses.data(), qscores.data(),
                                      charges.data(), wstarts.data(), wends.data(),
                                      ion_types.data(), frag_indices.data(), "EThcD");
    }

    for (int i = 0; i < count; i++)
    {
      double center_mz = (wstarts[i] + wends[i]) / 2.0;
      double iso_width = wends[i] - wstarts[i];
      targets.push_back({center_mz, charges[i], iso_width, ion_types[i], frag_indices[i]});
    }
    return targets;
  }

  void FLASHIda::pushFollowUpMS2_(const ScanCommand& ctx)
  {
    if (ms2_configs_.size() < 2)
      return;

    ScanCommand cmd = ctx;
    cmd.scan_id = nextTrackingIdInt_();
    cmd.priority = 2;

    // Use second MS2 config for follow-up
    std::strncpy(cmd.analyzer, ms2_configs_[1].analyzer.c_str(), sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    cmd.orbitrap_resolution = ms2_configs_[1].resolution;
    cmd.stages[0].collision_energy = static_cast<double>(ms2_configs_[1].collision_energy);
    std::strncpy(cmd.stages[0].activation_type, ms2_configs_[1].activation.c_str(),
                 sizeof(cmd.stages[0].activation_type) - 1);
    cmd.stages[0].activation_type[sizeof(cmd.stages[0].activation_type) - 1] = '\0';

    std::string id_str = encodeBase36_(cmd.scan_id);
    char desc_buf[256];
    std::snprintf(desc_buf, sizeof(desc_buf), "%s|followup mz=%.2f", id_str.c_str(), cmd.stages[0].precursor_mz);
    std::strncpy(cmd.scan_description, desc_buf, sizeof(cmd.scan_description) - 1);
    cmd.scan_description[sizeof(cmd.scan_description) - 1] = '\0';

    cmd.enqueue_timestamp_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    pending_scan_map_[cmd.scan_id] = cmd;
    pushCommand_(cmd);

    std::cout << "[TRACK-CREATE] id=" << id_str
              << " ms_level=2 type=followup"
              << std::endl;
  }

  void FLASHIda::pushConditionalFollowUp_(const ScanCommand& ctx)
  {
    if (ms2_configs_.size() < 2)
      return;

    ScanCommand cmd = ctx;
    cmd.scan_id = nextTrackingIdInt_();
    cmd.priority = 2;

    // Use second MS2 config for conditional follow-up
    std::strncpy(cmd.analyzer, ms2_configs_[1].analyzer.c_str(), sizeof(cmd.analyzer) - 1);
    cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
    cmd.orbitrap_resolution = ms2_configs_[1].resolution;
    cmd.stages[0].collision_energy = static_cast<double>(ms2_configs_[1].collision_energy);
    std::strncpy(cmd.stages[0].activation_type, ms2_configs_[1].activation.c_str(),
                 sizeof(cmd.stages[0].activation_type) - 1);
    cmd.stages[0].activation_type[sizeof(cmd.stages[0].activation_type) - 1] = '\0';

    std::string id_str = encodeBase36_(cmd.scan_id);
    char desc_buf[256];
    std::snprintf(desc_buf, sizeof(desc_buf), "%s|conditional mz=%.2f", id_str.c_str(), cmd.stages[0].precursor_mz);
    std::strncpy(cmd.scan_description, desc_buf, sizeof(cmd.scan_description) - 1);
    cmd.scan_description[sizeof(cmd.scan_description) - 1] = '\0';

    cmd.enqueue_timestamp_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    pending_scan_map_[cmd.scan_id] = cmd;
    pushCommand_(cmd);

    std::cout << "[TRACK-CREATE] id=" << id_str
              << " ms_level=2 type=conditional"
              << std::endl;
  }

  // --- Phase 7: Exploration engine implementation ---

  std::vector<double> FLASHIda::buildCEVariants_(double ce_min, double ce_max, double ce_step) const
  {
    std::vector<double> ces;
    for (double ce = ce_min; ce <= ce_max + 1e-9; ce += ce_step)
      ces.push_back(ce);
    return ces;
  }

  ScanCommand FLASHIda::buildMS2Command_(double precursor_mz, int charge, double ce,
                                          const std::string& activation)
  {
    ScanCommand cmd{};
    cmd.msn_level = 2;
    cmd.priority = 1;
    cmd.is_agc = 0;
    cmd.num_stages = 1;
    cmd.scan_id = nextTrackingIdInt_();

    if (!ms2_configs_.empty())
    {
      std::strncpy(cmd.analyzer, ms2_configs_[0].analyzer.c_str(), sizeof(cmd.analyzer) - 1);
      cmd.analyzer[sizeof(cmd.analyzer) - 1] = '\0';
      cmd.orbitrap_resolution = ms2_configs_[0].resolution;
    }
    cmd.first_mass = ms1_first_mass_;
    cmd.last_mass = ms1_last_mass_;
    cmd.max_it = ms1_max_it_;
    cmd.agc_target = ms1_agc_target_;

    cmd.stages[0].precursor_mz = precursor_mz;
    cmd.stages[0].isolation_width = 2.0;
    cmd.stages[0].collision_energy = ce;
    cmd.stages[0].charge_state = charge;
    std::strncpy(cmd.stages[0].activation_type, activation.c_str(),
                 sizeof(cmd.stages[0].activation_type) - 1);
    cmd.stages[0].activation_type[sizeof(cmd.stages[0].activation_type) - 1] = '\0';

    cmd.enqueue_timestamp_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    return cmd;
  }

  void FLASHIda::applyOverrides_(ScanCommand& cmd,
      const std::unordered_map<std::string, std::string>& overrides) const
  {
    for (auto ov_it = overrides.begin(); ov_it != overrides.end(); ++ov_it)
    {
      const auto& okey = ov_it->first;
      const auto& oval = ov_it->second;
      if (okey == "agc_target") cmd.agc_target = static_cast<int32_t>(std::stod(oval));
      else if (okey == "max_injection_time_ms") cmd.max_it = std::stod(oval);
      else if (okey == "isolation_width") cmd.stages[0].isolation_width = std::stod(oval);
    }
  }

  void FLASHIda::initiateExploration_(int msn_level, double precursor_mz,
      double precursor_mass, int precursor_charge, double faims_cv)
  {
    const auto& cfg = getLevelConfig_(msn_level);
    if (!hasExploration(cfg)) return;

    std::vector<double> ces = buildCEVariants_(cfg.ce_min, cfg.ce_max, cfg.ce_step);
    if (ces.empty()) return;

    ExplorationGroup group;
    group.group_id = next_exploration_group_id_++;
    group.msn_level = msn_level;
    group.exploration_metric = cfg.exploration;
    group.precursor_mz = precursor_mz;
    group.precursor_mass = precursor_mass;
    group.precursor_charge = precursor_charge;
    group.isolation_width = 2.0;
    group.faims_cv = faims_cv;
    group.start_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    for (int i = 0; i < static_cast<int>(ces.size()); ++i)
    {
      ExplorationVariant v;
      v.variant_index = i;
      v.collision_energy = ces[i];
      v.activation_type = cfg.activation;

      ScanCommand cmd = buildMS2Command_(precursor_mz, precursor_charge, ces[i], cfg.activation);
      cmd.priority = 0;
      cmd.faims_cv = faims_cv;
      applyOverrides_(cmd, cfg.overrides);

      int id_int = cmd.scan_id;
      std::string id_str = encodeBase36_(id_int);
      v.tracking_id = id_str;
      std::snprintf(cmd.scan_description, sizeof(cmd.scan_description),
                   "%s|EXPL CE=%.1f %.2f@%d", id_str.c_str(),
                   ces[i], precursor_mass, precursor_charge);

      group.variants.push_back(v);
      variant_tracking_to_group_[id_int] = {group.group_id, i};
      queues_[0].push_back(cmd);

      std::cout << "[TRACK-CREATE] id=" << id_str
                << " ms_level=" << msn_level << " type=exploration"
                << " CE=" << ces[i] << std::endl;
    }

    active_exploration_groups_[group.group_id] = std::move(group);
  }

  double FLASHIda::computeExplorationScore_(ExplorationMetric metric,
      const DeconvolvedSpectrum& spec) const
  {
    switch (metric)
    {
      case ExplorationMetric::MassCount:
        return computeMassCount_(spec);
      case ExplorationMetric::RemainingPrecursor:
        return computeRemainingPrecursorScore_(spec);
      case ExplorationMetric::FragmentCount:
        return computeFragmentCount_(spec);
      default:
        return computeMassCount_(spec);
    }
  }

  double FLASHIda::computeMassCount_(const DeconvolvedSpectrum& spec) const
  {
    return static_cast<double>(spec.size());
  }

  double FLASHIda::computeRemainingPrecursorScore_(const DeconvolvedSpectrum& spec) const
  {
    if (spec.empty()) return 0.0;
    double tic = 0.0;
    for (Size i = 0; i < spec.size(); ++i)
      tic += spec[i].getIntensity();
    return tic;
  }

  double FLASHIda::computeFragmentCount_(const DeconvolvedSpectrum& spec) const
  {
    return static_cast<double>(spec.size());
  }

  float FLASHIda::computeTICCoverage_(const DeconvolvedSpectrum& spec) const
  {
    if (spec.empty()) return 0.0f;
    float total = 0.0f;
    for (Size i = 0; i < spec.size(); ++i)
      total += static_cast<float>(spec[i].getIntensity());
    return total > 0.0f ? 1.0f : 0.0f;
  }

  void FLASHIda::feedExplorationResult_(int group_id, int variant_index,
      const DeconvolvedSpectrum& ms2_deconv, double rt)
  {
    (void)rt;
    auto git = active_exploration_groups_.find(group_id);
    if (git == active_exploration_groups_.end()) return;
    ExplorationGroup& group = git->second;

    if (variant_index < 0 || variant_index >= static_cast<int>(group.variants.size())) return;
    ExplorationVariant& v = group.variants[variant_index];
    if (v.received) return;

    v.result = ms2_deconv;
    v.score = computeExplorationScore_(group.exploration_metric, ms2_deconv);
    v.tic_coverage = computeTICCoverage_(ms2_deconv);
    v.fragment_count = static_cast<int>(ms2_deconv.size());
    v.received = true;

    auto& meta = v.result.getOrCreateOptimizationMetadata();
    meta.group_id = group.group_id;
    meta.variant_index = variant_index;
    meta.total_variants = static_cast<int>(group.variants.size());
    meta.is_best_variant = false;
    meta.msn_level_optimized = group.msn_level;
    meta.exploration_metric = static_cast<int>(group.exploration_metric);
    meta.collision_energy = v.collision_energy;
    meta.activation_type = v.activation_type;
    meta.precursor_mass = group.precursor_mass;
    meta.precursor_charge = group.precursor_charge;
    meta.fragmentation_quality_score = v.score;
    meta.tic_coverage = v.tic_coverage;
    meta.fragment_count = v.fragment_count;
    meta.start_ms = group.start_ms;
    meta.complete_ms = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());
    meta.exploration_scans = static_cast<int>(group.variants.size());

    bool all_received = std::all_of(group.variants.begin(), group.variants.end(),
                                    [](const ExplorationVariant& x){ return x.received; });
    if (!all_received) return;

    int best_idx = 0;
    double best_score = group.variants[0].score;
    for (int i = 1; i < static_cast<int>(group.variants.size()); ++i)
    {
      if (group.variants[i].score > best_score)
      {
        best_score = group.variants[i].score;
        best_idx = i;
      }
    }
    group.winner_index = best_idx;
    group.complete = true;
    group.variants[best_idx].result.getOrCreateOptimizationMetadata().is_best_variant = true;

    std::cout << "[EXPL-WINNER] group=" << group.group_id
              << " winner_idx=" << best_idx
              << " CE=" << group.variants[best_idx].collision_energy
              << " score=" << best_score << std::endl;

    const auto& level_config = getLevelConfig_(group.msn_level);
    if (!level_config.overrides.empty())
    {
      ScanCommand prod_cmd = buildMS2Command_(
          group.precursor_mz, group.precursor_charge,
          group.variants[best_idx].collision_energy,
          group.variants[best_idx].activation_type);
      prod_cmd.faims_cv = group.faims_cv;
      prod_cmd.priority = 1;
      queues_[1].push_back(prod_cmd);

      std::string prod_id = encodeBase36_(prod_cmd.scan_id);
      std::cout << "[TRACK-CREATE] id=" << prod_id
                << " ms_level=" << group.msn_level << " type=production"
                << std::endl;
    }
    else
    {
      initiateNextLevel_(group.msn_level,
          group.variants[best_idx].result, group.faims_cv);
    }

    for (const auto& vr : group.variants)
      variant_tracking_to_group_.erase(decodeBase36_(vr.tracking_id));

    active_exploration_groups_.erase(git);
  }

  void FLASHIda::initiateNextLevel_(int msn_level, const DeconvolvedSpectrum& result,
      double faims_cv)
  {
    int next_level = msn_level + 1;
    const auto& next_cfg = getLevelConfig_(next_level);
    if (next_cfg.selection == SelectionMetric::None) return;

    std::vector<std::pair<double, double>> targets;  // (mass, intensity)
    for (Size i = 0; i < result.size(); ++i)
      targets.push_back({result[i].getMonoMass(), static_cast<double>(result[i].getIntensity())});

    std::sort(targets.begin(), targets.end(),
              [](const auto& a, const auto& b){ return a.second > b.second; });
    int num_targets = std::min(static_cast<int>(targets.size()), next_cfg.max_targets);

    if (hasExploration(next_cfg))
    {
      for (int ti = 0; ti < num_targets; ++ti)
        initiateExploration_(next_level, targets[ti].first, 0.0, 0, faims_cv);
    }
    else
    {
      for (int ti = 0; ti < num_targets; ++ti)
      {
        ScanCommand cmd = buildMS2Command_(targets[ti].first, 0,
            next_cfg.ce_min, next_cfg.activation);
        cmd.msn_level = next_level;
        cmd.faims_cv = faims_cv;
        cmd.priority = 1;
        pushCommand_(cmd);

        std::string id_str = encodeBase36_(cmd.scan_id);
        std::cout << "[TRACK-CREATE] id=" << id_str
                  << " ms_level=" << next_level << " type=next_level"
                  << std::endl;
      }
    }
  }

  int FLASHIda::processMS2Path_(const double* mzs, const double* ints, int length,
                                 double rt_min, const char* scan_desc)
  {
    int commands_pushed = 0;

    // Step 1: Decode tracking ID from scan_description format: {base36_id}|{payload}
    std::string desc_str = scan_desc ? scan_desc : "";
    if (desc_str.empty())
      return 0;

    Size pipe_pos = desc_str.find('|');
    if (pipe_pos == std::string::npos)
      return 0;

    std::string id_str = desc_str.substr(0, pipe_pos);
    int tracking_id = decodeBase36_(id_str);

    // Phase 7: Check if this is an exploration variant (before pending_scan_map_)
    auto vit = variant_tracking_to_group_.find(tracking_id);
    if (vit != variant_tracking_to_group_.end())
    {
      int gid = vit->second.group_id;
      int vidx = vit->second.variant_index;
      variant_tracking_to_group_.erase(vit);

      // Deconvolve the MS2 result for scoring
      DeconvolvedSpectrum ms2_deconv(tracking_id);
      if (mzs != nullptr && ints != nullptr && length > 0)
      {
        deconvolveMS2(mzs, ints, length, rt_min, 0.0, 0);
        ms2_deconv = ms2_deconvolved_spectrum_;
      }

      feedExplorationResult_(gid, vidx, ms2_deconv, rt_min);
      return commands_pushed;
    }

    // Step 2: Look up pending scan context
    auto it = pending_scan_map_.find(tracking_id);
    if (it == pending_scan_map_.end())
    {
      std::cout << "[TRACK-RESOLVE] id=" << id_str << " status=not_found" << std::endl;
      return 0;
    }
    ScanCommand ctx = it->second;
    pending_scan_map_.erase(it);

    // Step 3: Deconvolve MS2 with precursor context
    double precursor_mass = 0;
    int precursor_charge = 0;
    if (ctx.num_stages > 0)
    {
      precursor_charge = ctx.stages[0].charge_state;
      precursor_mass = ctx.stages[0].precursor_mz * precursor_charge
                       - precursor_charge * FLASHHelperClasses::getChargeMass(true);
    }

    deconvolveMS2(mzs, ints, length, rt_min, precursor_mass, precursor_charge);

    // Step 4: Route by mode
    // Tag-based targeting
    bool tags_found = false;
    if (tag_based_targeting_enabled_ && precursor_mass > 0)
    {
      tags_found = processMS2ForTagBasedTargeting(precursor_mass);
    }

    // Quantification follow-up (independent of tags)
    if (quant_enabled_ && ms2_configs_.size() >= 2)
    {
      if (isDifferentiallyAbundant(mzs, ints, length, rt_min, 2, "ms2_quant",
                                    reporter_mz_tol_, fold_change_threshold_, false))
      {
        pushFollowUpMS2_(ctx);
        commands_pushed++;
      }
    }

    // Conditional MS2 follow-up — only when tags detected
    if (conditional_ms2_enabled_ && ms2_configs_.size() >= 2 && tags_found)
    {
      pushConditionalFollowUp_(ctx);
      commands_pushed++;
    }

    // Step 5: MS3 targeting — uses level_configs_ when exploration is configured,
    // falls back to legacy ms3_enabled_ path for standard MS3 targeting
    if (hasExploration(getLevelConfig_(3)))
    {
      // MS3 exploration: create exploration groups for top fragments
      initiateNextLevel_(2, ms2_deconvolved_spectrum_, ctx.faims_cv);
    }
    else if (ms3_enabled_ && ms3_mode_ > 0)
    {
      // Legacy MS3 targeting (non-exploration)
      auto ms3_targets = selectMS3Targets_();
      for (const auto& t : ms3_targets)
      {
        ScanCommand ms3_cmd = buildMS3Command_(ctx, t.center_mz, t.charge, t.iso_width,
                                                t.ion_type, t.frag_index);
        pushCommand_(ms3_cmd);
        commands_pushed++;
      }
    }
    else if (getLevelConfig_(3).selection != SelectionMetric::None && !ms3_enabled_)
    {
      // New selection_strategy MS3 targeting (no exploration, not legacy)
      initiateNextLevel_(2, ms2_deconvolved_spectrum_, ctx.faims_cv);
    }

    std::cout << "[TRACK-RESOLVE] id=" << id_str
              << " rt=" << rt_min
              << " commands=" << commands_pushed
              << std::endl;

    return commands_pushed;
  }

  int FLASHIda::getNextScanCommand(ScanCommand& out)
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);

    // Step 1: AGC scan if needed
    if (needsAGCScan_())
    {
      out = makeAGCCommand_();
      out.faims_cv = faims_enabled_ ? faims_cv_values_[current_cv_index_] : 0.0;
      out.scan_id = nextTrackingIdInt_();
      last_agc_time_ = std::chrono::steady_clock::now();

      // Scan description with base-36 tracking ID
      std::string id_str = encodeBase36_(out.scan_id);
      char desc_buf[128];
      std::snprintf(desc_buf, sizeof(desc_buf), "%s|AGC calibration", id_str.c_str());
      std::strncpy(out.scan_description, desc_buf, sizeof(out.scan_description) - 1);
      out.scan_description[sizeof(out.scan_description) - 1] = '\0';

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=agc" << std::endl;
      return 1;
    }

    // Step 2: Cycle time — force MS1 if too long since last survey scan
    // Suppressed while any exploration group is active (Phase 7)
    bool exploration_active = !active_exploration_groups_.empty();
    if (cycle_time_enabled_ && !exploration_active
        && msSinceLastMS1_() > static_cast<uint64_t>(cycle_time_ms_))
    {
      out = makeMS1Command_();
      out.faims_cv = faims_enabled_ ? faims_cv_values_[current_cv_index_] : 0.0;
      out.scan_id = nextTrackingIdInt_();
      last_ms1_time_ = std::chrono::steady_clock::now();

      std::string id_str = encodeBase36_(out.scan_id);
      char desc_buf[128];
      std::snprintf(desc_buf, sizeof(desc_buf), "%s|MS1 survey", id_str.c_str());
      std::strncpy(out.scan_description, desc_buf, sizeof(out.scan_description) - 1);
      out.scan_description[sizeof(out.scan_description) - 1] = '\0';

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=cycle_time" << std::endl;
      return 1;
    }

    // Step 3: Cleanup expired commands
    cleanupExpiredCommands_();

    // Step 4: Dequeue by priority (0 = highest → 3 = lowest)
    for (int p = 0; p < 4; ++p)
    {
      if (!queues_[p].empty())
      {
        out = queues_[p].front();
        queues_[p].pop_front();
        // faims_cv already set at creation time (MS2 → parent CV, CV-transition MS1 → next CV)
        return 1;
      }
    }

    // Step 5: Queue empty — no commands available
    // The FAIMS CV-transition MS1 is pushed into the queue by processScan(),
    // so it will be dequeued in Step 4 on the next call. The C# while-loop
    // in ProcessSpectrum depends on returning 0 to stop draining.
    return 0;
  }

  int FLASHIda::getNextTrackingId()
  {
    std::lock_guard<std::mutex> lock(queue_mutex_);
    return nextTrackingIdInt_();
  }

} // namespace OpenMS
