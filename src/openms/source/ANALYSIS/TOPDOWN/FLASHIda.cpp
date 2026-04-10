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

// tracking_alphabet_ moved to ScanCommandQueue.cpp

/// optimal window margin
inline const double optimal_window_margin_ = .4;

/// constructor
FLASHIda::FLASHIda(char* arg) :
    config_(std::string(arg)),
    queue_(config_),
    deconv_(config_),
    selection_(config_, deconv_),
    faims_(config_),
    exploration_(config_)
{
  #ifdef _OPENMP
    omp_set_num_threads(4);
  #endif

    // AGC timing initialized in ScanCommandQueue constructor
    // Targeting/file-loading (log files, FASTA, TSV) moved to PrecursorSelection constructor
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
    MSSpectrum spec;
    for (int i = 0; i < length; i++)
    {
      if (ints[i] > 0) { spec.emplace_back(mzs[i], ints[i]); }
    }
    spec.setMSLevel(ms_level);
    spec.setName(name);
    spec.setRT(rt);

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

  // getPeakGroups() now inline in header, delegates to selection_.filterAndRank()

  // filterPeakGroupsUsingMassExclusion_() moved to PrecursorSelection

  // removeFromExlusionList() now inline in header, delegates to selection_.removeFromExclusionList()

  void FLASHIda::getAllMonoisotopicMasses(double* masses, int length)
  {
    int len = std::min(length, (int)deconv_.deconvolvedMS1().size());
    for (int i = 0; i < len; i++)
    {
      masses[i] = deconv_.deconvolvedMS1()[i].getMonoMass();
    }
  }

  int FLASHIda::GetAllPeakGroupSize()
  {
    return deconv_.deconvolvedMS1().size();
  }

  int FLASHIda::getBestMS2Masses(int n,
                                 double* masses,
                                 double* qscores,
                                 int* charges,
                                 double* window_starts,
                                 double* window_ends)
  {
    if (!deconv_.hasStoredMS2() || deconv_.storedMS2().empty())
    {
      return 0;
    }

    std::sort(deconv_.storedMS2().begin(), deconv_.storedMS2().end(),
    [](const PeakGroup& a, const PeakGroup& b) {
      return a.getChargeIntensity(a.getMaxIntensityAbsCharge()) > b.getChargeIntensity(b.getMaxIntensityAbsCharge());
    });

    int output_idx = 0;
    for (Size pg_idx = 0; pg_idx < deconv_.storedMS2().size() && output_idx < n; ++pg_idx)
    {
      const auto& pg = deconv_.storedMS2()[pg_idx];
      double mono_mass = pg.getMonoMass();

      if (config_.targeting().ms3_all_charges)
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

    if (!deconv_.hasStoredMS2() || deconv_.storedMS2().empty() || protein_sequence.empty())
    {
      return 0;
    }

    // 1. Copy and sort deconvolved spectrum
    DeconvolvedSpectrum dspec = deconv_.storedMS2();
    dspec.sort();
    double precursor_mass = dspec.getPrecursorPeakGroup().getMonoMass();
    std::cout << "PM=" << precursor_mass << std::endl;
    for (const auto& pg : dspec)
    {
      double mono_mass = pg.getMonoMass();
      std::cout << mono_mass << std::endl;
    }

    // 2. Run FLASHTagger to generate sequence tags
    double ppm_tolerance = config_.level(2).tolerance_ppm;
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
    for (Size peak_idx = 0; peak_idx < deconv_.storedMS2().size(); ++peak_idx)
    {
      const auto& pg = deconv_.storedMS2()[peak_idx];
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
      const auto& pg = deconv_.storedMS2()[m.peak_index];

      if (config_.targeting().ms3_all_charges)
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

      const auto& pg = deconv_.storedMS2()[ion.peak_index];
      double mono_mass = pg.getMonoMass();

      if (config_.targeting().ms3_all_charges)
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
        const auto& pg = deconv_.storedMS2()[selected->peak_index];

        if (config_.targeting().ms3_all_charges)
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

    if (deconv_.deconvolvedMS1().size() > max_count)
    {
      std::vector<float> intensites;
      intensites.reserve(deconv_.deconvolvedMS1().size());
      for (const auto& pg : deconv_.deconvolvedMS1())
      {
        intensites.push_back(pg.getIntensity());
      }
      std::sort(intensites.rbegin(), intensites.rend());
      threshold = intensites[max_count];
    }

    for (const auto& pg : deconv_.deconvolvedMS1())
    {
      if (pg.getIntensity() < threshold) continue;
      mass += pg.getMonoMass() * pg.getIntensity();
      intensity_sum += pg.getIntensity();
    }
    if (intensity_sum <= 0) return 0;
    return mass / intensity_sum;
    */
    auto filter_str = deconv_.deconvolvedMS1().getOriginalSpectrum().getMetaValue("filter string").toString();
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

  // getIsolationWindows() now inline in header, delegates to selection_.getIsolationWindows()

  int FLASHIda::getConfigInt(const std::string& key) const
  {
    return config_.getInt(key);
  }

  double FLASHIda::getConfigDouble(const std::string& key) const
  {
    return config_.getDouble(key);
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

  // parseInclusionListTSV_(), parseTargetPTMsTSV_(), generatePTMCombinations_(),
  // addDynamicTargets_(), processMS2ForTagBasedTargeting() moved to PrecursorSelection

  // parseJSONConfig_ removed: parsing logic moved to Config class

  // Phase 3 scan command queue and tracking methods moved to ScanCommandQueue.cpp

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

  int FLASHIda::processScan(const double* mzs, const double* ints, int length,
                             double rt_min, int ms_level, const char* scan_description,
                             double faims_cv)
  {
    std::lock_guard<std::mutex> lock(mutex_);

    if (ms_level == 1)
    {
      queue_.recordMS1Time();

      // Selection=none: skip MS1 precursor selection entirely
      if (config_.level(1).selection == SelectionMetric::None)
        return 0;

      // MS1 path: deconvolve, score, filter, select top-N, push MS2 commands
      double parent_cv = config_.faims().enabled ? faims_cv : 0.0;

      int n = selection_.filterAndRank(mzs, ints, length, rt_min, 1, nullptr);
      const auto& selected = selection_.selectedPeakGroups();
      const auto& sel_charges = selection_.triggerCharges();
      const auto& sel_hcds = selection_.triggerHcds();
      int commands_pushed = 0;
      for (int i = 0; i < n; i++)
      {
        ScanCommand cmd = queue_.buildMS2(selected[i], sel_charges[i], sel_hcds[i]);
        cmd.faims_cv = parent_cv;  // MS2 carries parent MS1's CV
        queue_.push(cmd);
        commands_pushed++;
      }

      // Phase 7: Initiate exploration for selected precursors if MS2 exploration is enabled
      if (config_.hasExploration(2))
      {
        for (int i = 0; i < n; i++)
        {
          auto [mz1, mz2] = selected[i].getMzRange(sel_charges[i]);
          double center_mz = (mz1 + mz2) / 2.0;
          auto cmds = exploration_.initiate(2, center_mz,
              selected[i].getMonoMass(),
              sel_charges[i], parent_cv, queue_);
          for (auto& c : cmds) queue_.push(c);
        }
      }

      // FAIMS CV cycling: update skip policy, advance to next CV, push MS1
      if (faims_.isEnabled())
      {
        double current_cv = faims_.currentCV();
        faims_.updateSkip(current_cv, commands_pushed);

        double next_cv = faims_.advanceToNextCV();
        ScanCommand ms1 = queue_.makeMS1();
        ms1.faims_cv = next_cv;
        ms1.scan_id = queue_.nextTrackingId();
        ms1.priority = 0;  // priority 0 to send before pending MS2s

        std::string id_str = ScanCommandQueue::encode(ms1.scan_id);
        std::snprintf(ms1.scan_description, 16, "%sS", id_str.c_str());

        auto now = std::chrono::steady_clock::now();
        ms1.enqueue_timestamp_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
            now.time_since_epoch()).count();

        queue_.push(ms1);
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

  // decodeTracking_ moved to ScanCommandQueue::decode()

  // buildMS2Command_(PeakGroup, charge, hcd) moved to ScanCommandQueue::buildMS2()

  // buildMS3Command_ moved to ScanCommandQueue::buildMS3()

  // pushCommand_ moved to ScanCommandQueue::push()

  std::vector<FLASHIda::MS3Target> FLASHIda::selectMS3Targets_()
  {
    std::vector<MS3Target> targets;
    if (config_.targeting().ms3_mode == 0 || !deconv_.hasStoredMS2())
      return targets;

    const int n = config_.level(3).max_targets;
    std::vector<double> masses(n), qscores(n), wstarts(n), wends(n);
    std::vector<int> charges(n);
    std::vector<char> ion_types(n, '\0');
    std::vector<int> frag_indices(n, 0);

    int count = 0;

    if (config_.targeting().ms3_mode == 1 || config_.targeting().ms3_mode == 2)
    {
      // Mode 1 (SourceCID) and Mode 2 (SPS): Use getBestMS2Masses
      count = getBestMS2Masses(n, masses.data(), qscores.data(), charges.data(),
                               wstarts.data(), wends.data());
    }
    else if (config_.targeting().ms3_mode == 3 && !config_.targeting().protein_sequence.empty())
    {
      // Mode 3 (HCD-triggered): Use getTopFragmentMatches
      count = getTopFragmentMatches(config_.targeting().protein_sequence, n, masses.data(), qscores.data(),
                                    charges.data(), wstarts.data(), wends.data(),
                                    ion_types.data(), frag_indices.data(), "HCD");
    }
    else if (config_.targeting().ms3_mode == 4 && !config_.targeting().protein_sequence.empty())
    {
      // Mode 4 (EThcD-triggered): Use getTerminalFragmentIons
      count = getTerminalFragmentIons(config_.targeting().protein_sequence, n, masses.data(), qscores.data(),
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

  // pushFollowUpMS2_ and pushConditionalFollowUp_ moved to
  // ScanCommandQueue::buildFollowUpMS2() and ScanCommandQueue::buildConditionalFollowUp()

  // --- Phase 7: Exploration engine moved to FLASHIda/Exploration.h/.cpp ---

  int FLASHIda::processMS2Path_(const double* mzs, const double* ints, int length,
                                 double rt_min, const char* scan_desc)
  {
    int commands_pushed = 0;

    // Step 1: Decode tracking ID from scan_description -- fixed position chars 0-2
    std::string desc_str = scan_desc ? scan_desc : "";
    if (desc_str.size() < 3)
      return 0;

    std::string id_str = desc_str.substr(0, 3);
    int tracking_id = queue_.decode(id_str);

    // Phase 7: Check if this is an exploration variant (before pending scan lookup)
    if (exploration_.isExplorationVariant(tracking_id))
    {
      // Deconvolve the MS2 result for scoring
      DeconvolvedSpectrum ms2_deconv(tracking_id);
      if (mzs != nullptr && ints != nullptr && length > 0)
      {
        deconv_.deconvolveMS2(mzs, ints, length, rt_min, 0.0, 0);
        ms2_deconv = deconv_.storedMS2();
      }

      auto cmds = exploration_.feedResult(tracking_id, ms2_deconv, rt_min, queue_);
      for (auto& c : cmds) queue_.push(c);
      return commands_pushed;
    }

    // Step 2: Look up pending scan context via queue
    auto resolved = queue_.resolvePending(tracking_id);
    if (!resolved.has_value())
    {
      std::cout << "[TRACK-RESOLVE] id=" << id_str << " status=not_found" << std::endl;
      return 0;
    }
    ScanCommand ctx = resolved.value();

    // Step 3: Deconvolve MS2 with precursor context
    double precursor_mass = 0;
    int precursor_charge = 0;
    if (ctx.num_stages > 0)
    {
      precursor_charge = ctx.stages[0].charge_state;
      precursor_mass = ctx.stages[0].precursor_mz * precursor_charge
                       - precursor_charge * FLASHHelperClasses::getChargeMass(true);
    }

    deconv_.deconvolveMS2(mzs, ints, length, rt_min, precursor_mass, precursor_charge);

    // Step 4: Route by mode
    // Tag-based targeting
    bool tags_found = false;
    if (selection_.hasTargetProteinDatabase() && precursor_mass > 0)
    {
      tags_found = selection_.processMS2ForTagBasedTargeting(precursor_mass);
    }

    // Quantification follow-up (independent of tags)
    if (config_.quantification().enabled && config_.level(2).scans.size() >= 2)
    {
      if (isDifferentiallyAbundant(mzs, ints, length, rt_min, 2, "ms2_quant",
                                    config_.quantification().reporter_mz_tol, config_.quantification().fold_change_threshold, false))
      {
        queue_.push(queue_.buildFollowUpMS2(ctx));
        commands_pushed++;
      }
    }

    // Conditional MS2 follow-up -- only when tags detected
    if (config_.level(2).scans.size() >= 2 && tags_found)
    {
      queue_.push(queue_.buildConditionalFollowUp(ctx));
      commands_pushed++;
    }

    // Step 5: MS3 targeting -- uses config levels when exploration is configured,
    // falls back to legacy MS3 targeting path for standard MS3 targeting
    if (config_.hasExploration(3))
    {
      // MS3 exploration: create exploration groups for top fragments
      auto cmds = exploration_.initiateNextLevel(2, deconv_.storedMS2(), ctx.faims_cv, queue_);
      for (auto& c : cmds) queue_.push(c);
    }
    else if (config_.targeting().ms3_mode > 0)
    {
      // Legacy MS3 targeting (non-exploration)
      auto ms3_targets = selectMS3Targets_();
      for (const auto& t : ms3_targets)
      {
        ScanCommand ms3_cmd = queue_.buildMS3(ctx, t.center_mz, t.charge, t.iso_width,
                                               t.ion_type, t.frag_index);
        queue_.push(ms3_cmd);
        commands_pushed++;
      }
    }
    else if (config_.level(3).selection != SelectionMetric::None && !(config_.targeting().ms3_mode > 0))
    {
      // New selection_strategy MS3 targeting (no exploration, not legacy)
      auto cmds = exploration_.initiateNextLevel(2, deconv_.storedMS2(), ctx.faims_cv, queue_);
      for (auto& c : cmds) queue_.push(c);
    }

    std::cout << "[TRACK-RESOLVE] id=" << id_str
              << " rt=" << rt_min
              << " commands=" << commands_pushed
              << std::endl;

    return commands_pushed;
  }

  int FLASHIda::getNextScanCommand(ScanCommand& out)
  {
    std::lock_guard<std::mutex> lock(mutex_);

    // Step 1: AGC scan if needed
    if (queue_.needsAGC())
    {
      out = queue_.makeAGC();
      out.faims_cv = faims_.isEnabled() ? faims_.currentCV() : 0.0;
      out.scan_id = queue_.nextTrackingId();
      queue_.recordAGCTime();

      // Scan description: {3-char ID}A
      std::string id_str = ScanCommandQueue::encode(out.scan_id);
      std::snprintf(out.scan_description, 16, "%sA", id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=agc" << std::endl;
      return 1;
    }

    // Step 2: Cycle time -- force MS1 if too long since last survey scan
    // Suppressed while any exploration group is active (Phase 7)
    bool exploration_active = exploration_.activeGroupCount() > 0;
    if (config_.scheduling().cycle_time_enabled && !exploration_active
        && queue_.msSinceLastMS1() > static_cast<uint64_t>(config_.scheduling().cycle_time_ms))
    {
      out = queue_.makeMS1();
      out.faims_cv = faims_.isEnabled() ? faims_.currentCV() : 0.0;
      out.scan_id = queue_.nextTrackingId();
      queue_.recordMS1Time();

      std::string id_str = ScanCommandQueue::encode(out.scan_id);
      std::snprintf(out.scan_description, 16, "%sS", id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=cycle_time" << std::endl;
      return 1;
    }

    // Step 3: Cleanup expired commands
    queue_.cleanupExpired();

    // Step 4: Dequeue by priority (0 = highest -> 3 = lowest)
    auto dequeued = queue_.dequeue();
    if (dequeued.has_value())
    {
      out = dequeued.value();
      // faims_cv already set at creation time (MS2 -> parent CV, CV-transition MS1 -> next CV)
      return 1;
    }

    // Step 5: Idle cycle -- queue empty, keep the instrument busy with AGC + MS1
    // Create an AGC command (returned immediately) and push an MS1 at priority 0
    // into the queue so the next dequeue returns it before any MS2s (priority 1+).
    {
      // 5a: AGC
      ScanCommand agc_cmd = queue_.makeAGC();
      agc_cmd.faims_cv = faims_.isEnabled() ? faims_.currentCV() : 0.0;
      agc_cmd.scan_id = queue_.nextTrackingId();
      queue_.recordAGCTime();

      std::string agc_id_str = ScanCommandQueue::encode(agc_cmd.scan_id);
      std::snprintf(agc_cmd.scan_description, 16, "%sA", agc_id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << agc_id_str << " ms_level=1 type=idle_agc" << std::endl;

      // 5b: MS1 -- override priority to 0 (makeMS1 defaults to 3)
      ScanCommand ms1_cmd = queue_.makeMS1();
      ms1_cmd.faims_cv = faims_.isEnabled() ? faims_.currentCV() : 0.0;
      ms1_cmd.scan_id = queue_.nextTrackingId();
      ms1_cmd.priority = 0;
      queue_.recordMS1Time();

      std::string ms1_id_str = ScanCommandQueue::encode(ms1_cmd.scan_id);
      std::snprintf(ms1_cmd.scan_description, 16, "%sS", ms1_id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << ms1_id_str << " ms_level=1 type=idle_ms1" << std::endl;

      // Push MS1 into priority-0 queue for next dequeue call
      queue_.push(ms1_cmd);

      out = agc_cmd;
      return 1;
    }
  }

  int FLASHIda::getNextTrackingId()
  {
    return queue_.nextTrackingId();
  }

} // namespace OpenMS
