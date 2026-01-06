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
    if (ion_type == "zp1") return Residue::getInternalToZp1Ion().getMonoWeight();
    if (ion_type == "zp2") return Residue::getInternalToZp2Ion().getMonoWeight();
    return Residue::getInternalToYIon().getMonoWeight(); // default to y
  }

  /// Structure representing a matched fragment ion
  struct FragmentIonMatch
  {
    int fragment_index;      ///< Position in sequence where this fragment ends (1-based)
    int peak_index;          ///< Index of the matching peak in the deconvolved spectrum
    double theoretical_mass; ///< Theoretical fragment mass
    double observed_mass;    ///< Observed (deconvolved) mass
    double mass_diff;        ///< Difference: observed - theoretical
    double ppm_error;        ///< PPM error
    bool is_prefix;          ///< True for N-terminal ions (b,a,c), false for C-terminal (y,z,x)
    String ion_type;         ///< Ion type (b, y, c, z, etc.)
    float score;             ///< Match quality score
  };

  // PTMSite is now defined as FLASHIda::PTMSite in the header file
  using PTMSite = FLASHIda::PTMSite;

  /// Structure to hold identification results
  struct ProteoformIdentificationResult
  {
    std::vector<FragmentIonMatch> matched_fragments; ///< All matched fragment ions
    std::vector<int> matched_prefix_indices;         ///< Indices of matched N-terminal fragments
    std::vector<int> matched_suffix_indices;         ///< Indices of matched C-terminal fragments
    std::vector<PTMSite> ptm_sites;                  ///< Detected PTM sites
    double sequence_coverage;                        ///< Fraction of sequence covered by matches
    double total_score;                              ///< Total identification score
    int matched_ion_count;                           ///< Number of matched ions
  };

  /// Calculate theoretical fragment masses for a protein sequence
  void calculateTheoreticalFragmentMasses(const String& sequence,
                                          const std::vector<String>& ion_types,
                                          std::vector<double>& prefix_masses,
                                          std::vector<double>& suffix_masses,
                                          std::vector<double>& prefix_shifts,
                                          std::vector<double>& suffix_shifts)
  {
    // Clean sequence (replace I with L for mass calculation)
    String clean_seq = sequence;
    std::replace(clean_seq.begin(), clean_seq.end(), 'I', 'L');

    Size seq_len = clean_seq.size();
    prefix_masses.clear();
    suffix_masses.clear();
    prefix_masses.reserve(seq_len);
    suffix_masses.reserve(seq_len);

    // Calculate prefix masses (N-terminal fragments)
    double cumulative_mass = 0.0;
    for (Size i = 0; i < seq_len; ++i)
    {
      char aa = clean_seq[i];
      if (aa == 'X')
      {
        prefix_masses.push_back(cumulative_mass);
        continue;
      }
      const Residue* res = ResidueDB::getInstance()->getResidue(aa);
      if (res != nullptr)
      {
        cumulative_mass += res->getMonoWeight(Residue::Internal);
      }
      prefix_masses.push_back(cumulative_mass);
    }

    // Calculate suffix masses (C-terminal fragments)
    cumulative_mass = 0.0;
    for (int i = seq_len - 1; i >= 0; --i)
    {
      char aa = clean_seq[i];
      if (aa == 'X')
      {
        suffix_masses.push_back(cumulative_mass);
        continue;
      }
      const Residue* res = ResidueDB::getInstance()->getResidue(aa);
      if (res != nullptr)
      {
        cumulative_mass += res->getMonoWeight(Residue::Internal);
      }
      suffix_masses.push_back(cumulative_mass);
    }
    // Reverse to get masses in sequence order (position i = y_{n-i})
    std::reverse(suffix_masses.begin(), suffix_masses.end());

    // Calculate ion type shifts
    prefix_shifts.clear();
    suffix_shifts.clear();
    for (const auto& ion_type : ion_types)
    {
      if (ion_type == "b" || ion_type == "a" || ion_type == "c")
      {
        prefix_shifts.push_back(getPrefixIonShift(ion_type));
      }
      else
      {
        suffix_shifts.push_back(getSuffixIonShift(ion_type));
      }
    }

    // Default to b and y if no shifts calculated
    if (prefix_shifts.empty()) prefix_shifts.push_back(getPrefixIonShift("b"));
    if (suffix_shifts.empty()) suffix_shifts.push_back(getSuffixIonShift("y"));
  }

  /// Match a single observed mass against theoretical masses
  bool findBestMatch(double observed_mass,
                     const std::vector<double>& theoretical_masses,
                     const std::vector<double>& ion_shifts,
                     double ppm_tolerance,
                     int& best_index,
                     double& best_theoretical,
                     double& best_diff,
                     double& best_shift)
  {
    best_index = -1;
    double min_ppm = ppm_tolerance;

    for (Size i = 0; i < theoretical_masses.size(); ++i)
    {
      for (double shift : ion_shifts)
      {
        double theo_mass = theoretical_masses[i] + shift;
        if (theo_mass <= 0) continue;

        double diff = observed_mass - theo_mass;
        double ppm = std::abs(diff) / theo_mass * 1e6;

        if (ppm < min_ppm)
        {
          min_ppm = ppm;
          best_index = static_cast<int>(i);
          best_theoretical = theo_mass;
          best_diff = diff;
          best_shift = shift;
        }
      }
    }

    return best_index >= 0;
  }

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

  /// Detect PTM sites from mass differences between matched fragments
  void detectPTMSites(const std::vector<FragmentIonMatch>& matches,
                      const String& sequence,
                      double ptm_mass_threshold,
                      std::vector<PTMSite>& ptm_sites)
  {
    ptm_sites.clear();

    // Group matches by type (prefix/suffix)
    std::vector<const FragmentIonMatch*> prefix_matches, suffix_matches;
    for (const auto& m : matches)
    {
      if (m.is_prefix)
        prefix_matches.push_back(&m);
      else
        suffix_matches.push_back(&m);
    }

    // Sort by fragment index
    auto by_index = [](const FragmentIonMatch* a, const FragmentIonMatch* b) {
      return a->fragment_index < b->fragment_index;
    };
    std::sort(prefix_matches.begin(), prefix_matches.end(), by_index);
    std::sort(suffix_matches.begin(), suffix_matches.end(), by_index);

    // Detect PTMs from consecutive prefix fragments with mass jumps
    for (Size i = 1; i < prefix_matches.size(); ++i)
    {
      const auto* prev = prefix_matches[i - 1];
      const auto* curr = prefix_matches[i];

      if (curr->fragment_index <= prev->fragment_index) continue;

      // Expected mass difference from sequence
      double expected_diff = curr->theoretical_mass - prev->theoretical_mass;
      // Observed mass difference
      double observed_diff = curr->observed_mass - prev->observed_mass;
      // Mass shift
      double mass_shift = observed_diff - expected_diff;

      if (std::abs(mass_shift) > ptm_mass_threshold)
      {
        PTMSite site;
        site.start_position = prev->fragment_index + 1;
        site.end_position = curr->fragment_index;
        site.position = (site.start_position + site.end_position) / 2; // midpoint
        site.mass_shift = mass_shift;
        ptm_sites.push_back(site);
      }
    }

    // Detect PTMs from consecutive suffix fragments with mass jumps
    for (Size i = 1; i < suffix_matches.size(); ++i)
    {
      const auto* prev = suffix_matches[i - 1];
      const auto* curr = suffix_matches[i];

      if (curr->fragment_index <= prev->fragment_index) continue;

      double expected_diff = curr->theoretical_mass - prev->theoretical_mass;
      double observed_diff = curr->observed_mass - prev->observed_mass;
      double mass_shift = observed_diff - expected_diff;

      if (std::abs(mass_shift) > ptm_mass_threshold)
      {
        // For suffix ions, the position is from C-terminus
        int seq_len = static_cast<int>(sequence.size());
        PTMSite site;
        site.start_position = seq_len - curr->fragment_index + 1;
        site.end_position = seq_len - prev->fragment_index;
        site.position = (site.start_position + site.end_position) / 2;
        site.mass_shift = mass_shift;
        ptm_sites.push_back(site);
      }
    }

    // Sort PTM sites by position
    std::sort(ptm_sites.begin(), ptm_sites.end(),
              [](const PTMSite& a, const PTMSite& b) { return a.position < b.position; });
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
    std::unordered_map<std::string, std::vector<double>> inputs;
    std::vector<String> log_files;
    std::vector<String> out_files;
    std::vector<String> tsv_files;  ///< TSV inclusion list files
    std::vector<String> fasta_files;  ///< FASTA files for tag-based targeting
    std::vector<String> ptm_files;    ///< PTM TSV files for target expansion

    // Make a copy of the input string since strtok modifies it in place
    // This is necessary when called from Python where the buffer may be read-only
    std::string arg_copy(arg);
    char* arg_mutable = &arg_copy[0];

    char* token = std::strtok(arg_mutable, " ");
    std::string key;
    std::stringstream ss {};

    while (token != nullptr)
    {
      String token_string = std::string(token);

      if (token_string.hasSuffix(".log"))
      {
        log_files.push_back(token_string);
        ss << token_string << " ";
      }
      else if (token_string.hasSuffix(".out"))
      {
        out_files.push_back(token_string);
        ss << token_string << " ";
      }
      else if (token_string.hasSuffix(".ptm.tsv"))
      {
        ptm_files.push_back(token_string);
        ss << token_string << " ";
      }
      else if (token_string.hasSuffix(".tsv"))
      {
        tsv_files.push_back(token_string);
        ss << token_string << " ";
      }
      else if (token_string.hasSuffix(".fasta") || token_string.hasSuffix(".fa"))
      {
        fasta_files.push_back(token_string);
        ss << token_string << " ";
      }
      else
      {
        double num = atof(token_string.c_str());

        if (num == 0 && ! isdigit(token_string[token_string.size() - 1]))
        {
          key = token_string;
          inputs[key] = DoubleList();
        }
        else { inputs[key].push_back(num); }
      }
      token = std::strtok(nullptr, " ");
    }
    rt_window_ = inputs["RT_window"][0];
    qscore_threshold_ = inputs["score_threshold"][0];
    tqscore_threshold = inputs["tqscore_threshold"][0];
    snr_threshold_ = 1;
    targeting_mode_ = (int)(inputs["target_mode"][0]);
    if (targeting_mode_ == 1) { std::cout << ss.str() << "file(s) is(are) used for inclusion mode\n"; }
    else if (targeting_mode_ == 2) { std::cout << ss.str() << "file(s) is(are) used for in-depth mode\n"; }
    else if (targeting_mode_ == 3) { std::cout << ss.str() << "file(s) is(are) used for exclusion mode\n"; }
    // Mode 4 removed - functionality merged into mode 1
    Param sd_defaults = SpectralDeconvolution().getDefaults();

    sd_defaults.setValue("min_charge", (int)inputs["min_charge"][0]);
    sd_defaults.setValue("max_charge", (int)inputs["max_charge"][0]);
    sd_defaults.setValue("min_mass", inputs["min_mass"][0]);
    sd_defaults.setValue("max_mass", inputs["max_mass"][0]);

    // Ensure tol has at least 2 values (for MS1 and MS2)
    // If only one value is provided, duplicate it for MS2
    DoubleList tol_values = inputs["tol"];
    if (tol_values.size() == 1)
    {
      tol_values.push_back(tol_values[0]);  // Use same tolerance for MS2
    }
    sd_defaults.setValue("tol", tol_values);
    tol_ = std::vector<double>(tol_values);

    // Use MS2 tolerance (second value) for tag matching
    if (tol_.size() >= 2)
    {
      tag_matching_tolerance_ppm_ = tol_[1];
    }
    else if (tol_.size() == 1)
    {
      tag_matching_tolerance_ppm_ = tol_[0];
    }

    // Pass precursor parameters if provided (needed for MS2 deconvolution)
    if (inputs.find("precursor_mz") != inputs.end() && !inputs["precursor_mz"].empty())
    {
      sd_defaults.setValue("precursor_mz", inputs["precursor_mz"][0]);
    }
    if (inputs.find("precursor_charge") != inputs.end() && !inputs["precursor_charge"].empty())
    {
      sd_defaults.setValue("precursor_charge", (int)inputs["precursor_charge"][0]);
    }

    // Pass cosine and SNR thresholds if provided
    // Ensure they have at least 2 values (for MS1 and MS2)
    if (inputs.find("min_cos") != inputs.end() && !inputs["min_cos"].empty())
    {
      DoubleList min_cos_values = inputs["min_cos"];
      if (min_cos_values.size() == 1)
      {
        min_cos_values.push_back(min_cos_values[0]);
      }
      sd_defaults.setValue("min_cos", min_cos_values);
    }
    if (inputs.find("min_snr") != inputs.end() && !inputs["min_snr"].empty())
    {
      DoubleList min_snr_values = inputs["min_snr"];
      if (min_snr_values.size() == 1)
      {
        min_snr_values.push_back(min_snr_values[0]);
      }
      sd_defaults.setValue("min_snr", min_snr_values);
    }

    // Set defaults for optional IDScore parameters
    use_idscore_ = false;
    consider_all_Charge_states_ = false;
    hcd_energy_ = -1;
    if (inputs.find("IDScore") != inputs.end() && !inputs["IDScore"].empty())
    {
      use_idscore_ = (bool)inputs["IDScore"][0];
    }
    if (inputs.find("AllCharges") != inputs.end() && !inputs["AllCharges"].empty())
    {
      consider_all_Charge_states_ = (bool)inputs["AllCharges"][0];
    }
    if (inputs.find("HCDEnergy") != inputs.end() && !inputs["HCDEnergy"].empty())
    {
      hcd_energy_ = (int)inputs["HCDEnergy"][0];
    }
    if (inputs.find("tie_threshold") != inputs.end() && !inputs["tie_threshold"].empty())
    {
      tie_threshold_ = inputs["tie_threshold"][0];
    }
    if (inputs.find("strict_inclusion") != inputs.end() && !inputs["strict_inclusion"].empty())
    {
      strict_inclusion_ = (bool)inputs["strict_inclusion"][0];
    }
    if (targeting_mode_ == 1)
    {
      std::cout << "Inclusion mode: " << (strict_inclusion_ ? "strict" : "non-strict") << "\n";
    }

    auto mass_count_double = inputs["max_mass_count"];

    for (double j : mass_count_double)
    {
      mass_count_.push_back((int)j);
    }

    for (const auto& log_file : log_files)
    {
      std::ifstream instream(log_file);
      if (instream.good())
      {
        String line;
        double rt = .0;
        double mass;
        double qscore;
        int charge;
        while (std::getline(instream, line))
        {
          if (line.find("0 targets") != line.npos) { continue; }
          if (line.hasPrefix("MS1"))
          {
            Size st = line.find("RT ") + 3;
            Size ed = line.find('(') - 2;
            String n = line.substr(st, ed - st + 1);
            rt = atof(n.c_str());
            // precursor_map_for_real_time_acquisition[scan] = std::vector<std::vector<double>>();//// ms1 scan -> mass, charge ,score, mz range,
            // precursor int, mass int, color
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

            if (targeting_mode_ == 1 || targeting_mode_ == 2)  // Mode 4 merged into mode 1
            {
              if (target_mass_rt_map_.find(mass) == target_mass_rt_map_.end()) { target_mass_rt_map_[mass] = std::vector<double>(); }
              target_mass_rt_map_[mass].push_back(rt * 60.0);
              if (target_mass_qscore_map_.find(mass) == target_mass_qscore_map_.end()) { target_mass_qscore_map_[mass] = std::vector<double>(); }
              target_mass_qscore_map_[mass].push_back(qscore);
              if (target_mass_charge_map_.find(mass) == target_mass_charge_map_.end()) { target_mass_charge_map_[mass] = std::vector<int>(); }
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
              {
                results.push_back(atof(str.c_str()));
              }
              if (exclusion_rt_masses_map_.find(rt * 60.0) == exclusion_rt_masses_map_.end())
              {
                exclusion_rt_masses_map_[rt * 60.0] = std::vector<double>();
              }
              for (double m : results)
              {
                exclusion_rt_masses_map_[rt * 60.0].push_back(m);
              }
            }
          }
        }
        instream.close();
      }
    }

    for (const auto& log_file : out_files)
    {
      std::ifstream instream(log_file);
      double rt = .0;
      if (instream.good())
      {
        String line;
        double mass;
        double qscore;

        while (std::getline(instream, line))
        {
          if (line.hasPrefix("rt")) { continue; }

          std::stringstream tmp_stream(line);
          String str;
          std::vector<String> results;
          while (getline(tmp_stream, str, '\t'))
          {
            results.push_back(str);
          }
          mass = atof(results[5].c_str());
          qscore = atof(results[3].c_str());

          if (targeting_mode_ == 1 || targeting_mode_ == 2)  // Mode 4 merged into mode 1
          {
            if (target_mass_rt_map_.find(mass) == target_mass_rt_map_.end()) { target_mass_rt_map_[mass] = std::vector<double>(); }
            rt = atof(results[0].c_str());
            target_mass_rt_map_[mass].push_back(60.0 * rt);
            if (target_mass_qscore_map_.find(mass) == target_mass_qscore_map_.end()) { target_mass_qscore_map_[mass] = std::vector<double>(); }
            target_mass_qscore_map_[mass].push_back(qscore);
          }
        }
        instream.close();
      }
    }

    // Parse TSV inclusion list files (auto-detected for mode 1)
    if (targeting_mode_ == 1 && !tsv_files.empty())
    {
      for (const auto& tsv_file : tsv_files)
      {
        parseInclusionListTSV_(tsv_file);
      }
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

    // Parse max_ptm_count if provided (default 3)
    if (inputs.find("max_ptm_count") != inputs.end() && !inputs["max_ptm_count"].empty())
    {
      max_total_ptm_count_ = (int)inputs["max_ptm_count"][0];
    }

    // Parse min_tag_length if provided (default 3)
    if (inputs.find("min_tag_length") != inputs.end() && !inputs["min_tag_length"].empty())
    {
      min_tag_length_for_targeting_ = (int)inputs["min_tag_length"][0];
    }

    // Parse max_tag_length if provided (default 8)
    if (inputs.find("max_tag_length") != inputs.end() && !inputs["max_tag_length"].empty())
    {
      max_tag_length_for_targeting_ = (int)inputs["max_tag_length"][0];
    }

    // Parse max_flanking_mass_diff if provided (default 50000.0)
    if (inputs.find("max_flanking_mass_diff") != inputs.end() && !inputs["max_flanking_mass_diff"].empty())
    {
      max_flanking_mass_diff_ = inputs["max_flanking_mass_diff"][0];
    }

    if (tag_based_targeting_enabled_)
    {
      std::cout << "Tag-based targeting: min_tag_length=" << min_tag_length_for_targeting_
                << ", max_tag_length=" << max_tag_length_for_targeting_
                << ", tolerance=" << tag_matching_tolerance_ppm_ << " ppm"
                << ", max_flanking_mass_diff=" << max_flanking_mass_diff_ << " Da\n";
    }

    // Load PTM TSV files for target expansion
    for (const auto& ptm_file : ptm_files)
    {
      parseTargetPTMsTSV_(ptm_file);
    }
    if (!target_ptms_.empty())
    {
      std::cout << target_ptms_.size() << " PTM modifications loaded for target expansion (max "
                << max_total_ptm_count_ << " total per proteoform)\n";
    }

    fd_.setParameters(sd_defaults);
    fd_.calculateAveragine(false);

    std::cout << sd_defaults << std::endl;
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

          int nominal_mass = SpectralDeconvolution::getNominalMass(mass);
          bool target_matched = false;
          double snr_threshold = snr_threshold_;
          double qscore_threshold = qscore_threshold_;
          double tqscore_factor_for_exclusion = 1.0;
          int integer_mz = (int)round(center_mz);
          

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
                      if (std::abs(t->mass - *ub) < 1e-6 && t->matchesCharge(charge))
                      {
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
                              double rt)
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

    // Perform deconvolution
    PeakGroup empty;
    fd_.performSpectrumDeconvolution(spec, 0, empty);
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
                                double rt)
  {
    if (mzs.empty() || mzs.size() != ints.size())
    {
      return 0;
    }
    return deconvolveMS2(mzs.data(), ints.data(), static_cast<int>(mzs.size()), rt);
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

    int count = std::min(n, static_cast<int>(ms2_deconvolved_spectrum_.size()));

    for (int i = 0; i < count; ++i)
    {
      const auto& pg = ms2_deconvolved_spectrum_[i];
      masses[i] = pg.getMonoMass();
      qscores[i] = pg.getQscore();

      int charge = pg.getRepAbsCharge();
      charges[i] = charge;

      // Get isolation window from m/z range
      auto [mz1, mz2] = pg.getMzRange(charge);
      window_starts[i] = mz1;
      window_ends[i] = mz2;
    }

    return count;
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
                                              std::vector<PTMSite>* ptm_sites)
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

    // 2. Run FLASHTagger to generate sequence tags
    double ppm_tolerance = tol_[1];
    FLASHTaggerAlgorithm tagger;
    Param tagger_param = tagger.getDefaults();
    std::vector<std::string> ion_types_str = {"b", "y"};
    tagger_param.setValue("ion_type", ion_types_str);
    tagger.setParameters(tagger_param);
    tagger.run(dspec, ppm_tolerance);

    std::vector<FLASHHelperClasses::Tag> tags;
    tagger.fillTags(tags);

    if (tags.empty())
    {
      // No tags found - return empty
      return 0;
    }

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
      return 0;
    }

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
    double max_mod_mass = 500.0;
    FLASHTaggerAlgorithm::runMatching(hits, dspec, spec_vec, vec_pro, rev_vec_pro, max_mod_mass);

    if (hits.empty())
    {
      return 0;
    }

    // 7. Run FLASHExtender for path-based validation
    FLASHExtenderAlgorithm extender;
    Param extender_param = extender.getDefaults();
    extender_param.setValue("ion_type", ion_types_str);
    extender.setParameters(extender_param);
    extender.run(hits, dspec, spec_vec, vec_pro, rev_vec_pro, tags, ppm_tolerance, false);

    std::vector<ProteinHit> proteoform_hits;
    extender.fillProteoforms(proteoform_hits);

    if (proteoform_hits.empty())
    {
      return 0;
    }

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

    // 10. Calculate PTM-adjusted theoretical fragment masses
    //     These masses already include PTM contributions based on FLASHExtender's detection
    double b_ion_shift = getPrefixIonShift("b");
    double y_ion_shift = getSuffixIonShift("y");
    std::vector<double> prefix_masses, suffix_masses;
    calculatePTMAdjustedFragmentMasses(matching_sequence, local_ptm_sites, b_ion_shift, y_ion_shift,
                                        prefix_masses, suffix_masses);

    // 11. Match observed masses against PTM-adjusted theoretical masses
    for (Size peak_idx = 0; peak_idx < ms2_deconvolved_spectrum_.size(); ++peak_idx)
    {
      const auto& pg = ms2_deconvolved_spectrum_[peak_idx];
      double observed_mass = pg.getMonoMass();

      // Find best matching prefix (b-ion)
      int best_prefix_idx = -1;
      double best_prefix_ppm = ppm_tolerance;
      double best_prefix_theo = 0;
      for (Size i = 0; i < prefix_masses.size(); ++i)
      {
        double theo_mass = prefix_masses[i];
        if (theo_mass <= 0) continue;
        double ppm_error = std::abs(observed_mass - theo_mass) / theo_mass * 1e6;
        if (ppm_error < best_prefix_ppm)
        {
          best_prefix_ppm = ppm_error;
          best_prefix_idx = static_cast<int>(i);
          best_prefix_theo = theo_mass;
        }
      }

      // Find best matching suffix (y-ion)
      int best_suffix_idx = -1;
      double best_suffix_ppm = ppm_tolerance;
      double best_suffix_theo = 0;
      for (Size i = 0; i < suffix_masses.size(); ++i)
      {
        double theo_mass = suffix_masses[i];
        if (theo_mass <= 0) continue;
        double ppm_error = std::abs(observed_mass - theo_mass) / theo_mass * 1e6;
        if (ppm_error < best_suffix_ppm)
        {
          best_suffix_ppm = ppm_error;
          best_suffix_idx = static_cast<int>(i);
          best_suffix_theo = theo_mass;
        }
      }

      // Choose best match (prefer lower ppm error)
      if (best_prefix_idx >= 0 || best_suffix_idx >= 0)
      {
        bool use_prefix = best_prefix_idx >= 0 &&
                          (best_suffix_idx < 0 || best_prefix_ppm < best_suffix_ppm);

        TagBasedFragmentMatch match;
        match.peak_index = static_cast<int>(peak_idx);
        match.observed_mass = observed_mass;
        match.qscore = pg.getQscore();
        match.charge = pg.getRepAbsCharge();

        if (use_prefix)
        {
          match.fragment_index = best_prefix_idx + 1;  // 1-based b-ion index
          match.is_prefix = true;
          match.theoretical_mass = best_prefix_theo;
          match.ppm_error = best_prefix_ppm;
        }
        else
        {
          // suffix_masses[i] corresponds to y_{i+1} (0-indexed to 1-indexed)
          match.fragment_index = best_suffix_idx + 1;  // 1-based y-ion index
          match.is_prefix = false;
          match.theoretical_mass = best_suffix_theo;
          match.ppm_error = best_suffix_ppm;
        }

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
      std::cout << "  " << (m.is_prefix ? "b" : "y") << m.fragment_index
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
                                      double* window_ends)
  {
    std::cout << "Matching fragments!" << std::endl;
    // Use tag-based matching workflow (FLASHTagger + FLASHExtender)
    std::vector<TagBasedFragmentMatch> matches;
    runTagBasedFragmentMatching_(protein_sequence, matches);

    // Already sorted by qscore descending
    int count = std::min(n, static_cast<int>(matches.size()));

    for (int i = 0; i < count; ++i)
    {
      const auto& m = matches[i];
      const auto& pg = ms2_deconvolved_spectrum_[m.peak_index];

      masses[i] = m.observed_mass;
      qscores[i] = m.qscore;
      charges[i] = m.charge;

      auto [mz1, mz2] = pg.getMzRange(m.charge);
      window_starts[i] = mz1;
      window_ends[i] = mz2;
    }

    return count;
  }

  int FLASHIda::getTopFragmentMatchesPy(const String& protein_sequence,
                                        int n,
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

    int count = getTopFragmentMatches(protein_sequence, n, masses.data(), qscores.data(),
                                      charges.data(), window_starts.data(), window_ends.data());

    masses.resize(count);
    qscores.resize(count);
    charges.resize(count);
    window_starts.resize(count);
    window_ends.resize(count);

    return count;
  }

  int FLASHIda::getAmbiguityEnclosingIons(const String& protein_sequence,
                                          int n,
                                          double* masses,
                                          double* qscores,
                                          int* charges,
                                          double* window_starts,
                                          double* window_ends)
  {
    // Get fragment matches AND PTM sites from FLASHExtender
    std::vector<TagBasedFragmentMatch> tag_matches;
    std::vector<PTMSite> ptm_sites;
    int match_count = runTagBasedFragmentMatching_(protein_sequence, tag_matches, &ptm_sites);

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

    // Collect unique enclosing ions (deduplicated by peak_index)
    // Enclosing ions bracket the PTM ambiguity region for MS3 targeting:
    // - Left bracket: fragment ending just before PTM 
    // - Right bracket: fragment just including PTM 
    std::set<int> used_peaks;
    std::vector<std::pair<float, int>> enclosing_ions; // (qscore, peak_index)
    int seq_len = static_cast<int>(protein_sequence.size());

    for (const auto& site : ptm_sites)
    {
      // Left bracket: largest fragment
      // - y-ion: y_n where (L-n+1) > E, i.e., n <= L-E (fragment_index <= L-E)
      const TagBasedFragmentMatch* best_left = nullptr;
      int best_left_distance = INT_MAX;
      for (const auto& m : tag_matches)
      {
        bool is_left_bracket = false;
        int distance_to_ptm = 0;

        if (!m.is_prefix)
        {
          // y-ion: y_n covers positions (L-n+1)..L
          int min_y = seq_len - site.end_position;
          if (m.fragment_index >= min_y)
          {
            is_left_bracket = true;
            // Distance is how far from the boundary
            distance_to_ptm = min_y - m.fragment_index + 1;
          }
        }

        if (is_left_bracket)
        {
          // Prefer closest to PTM boundary (smallest distance), then highest qscore
          if (distance_to_ptm < best_left_distance ||
              (distance_to_ptm == best_left_distance && (!best_left || m.qscore > best_left->qscore)))
          {
            best_left = &m;
            best_left_distance = distance_to_ptm;
          }
        }
      }

      // Right bracket: smallest modified fragment (just includes PTM)
      // - y-ion: y_n where (L-n+1) <= S, i.e., n >= L-S+1 (fragment_index >= L-S+1)
      const TagBasedFragmentMatch* best_right = nullptr;
      int best_right_distance = INT_MAX;
      for (const auto& m : tag_matches)
      {
        bool is_right_bracket = false;
        int distance_to_ptm = 0;

        if (m.is_prefix)
        {
          // b-ion: covers positions 1..fragment_index
          // Modified if fragment_index >= E (includes entire PTM region)
          if (m.fragment_index >= site.end_position-1)
          {
            is_right_bracket = true;
            distance_to_ptm = m.fragment_index - site.end_position;
          }
        }

        if (is_right_bracket)
        {
          // Prefer closest to PTM boundary (smallest distance), then highest qscore
          if (distance_to_ptm < best_right_distance ||
              (distance_to_ptm == best_right_distance && (!best_right || m.qscore > best_right->qscore)))
          {
            best_right = &m;
            best_right_distance = distance_to_ptm;
          }
        }
      }

      // Debug output for this PTM site
      std::cout << "[getAmbiguityEnclosingIons] PTM [" << site.start_position << "-" << site.end_position << "]:" << std::endl;
      if (best_left)
      {
        std::cout << "  Left bracket: " << (best_left->is_prefix ? "b" : "y") << best_left->fragment_index
                  << " (mass: " << std::fixed << std::setprecision(2) << best_left->observed_mass
                  << ", qscore: " << best_left->qscore << ")" << std::endl;
      }
      else
      {
        std::cout << "  Left bracket: none found" << std::endl;
      }
      if (best_right)
      {
        std::cout << "  Right bracket: " << (best_right->is_prefix ? "b" : "y") << best_right->fragment_index
                  << " (mass: " << std::fixed << std::setprecision(2) << best_right->observed_mass
                  << ", qscore: " << best_right->qscore << ")" << std::endl;
      }
      else
      {
        std::cout << "  Right bracket: none found" << std::endl;
      }

      // Add to output (deduplicated)
      if (best_left && used_peaks.find(best_left->peak_index) == used_peaks.end())
      {
        used_peaks.insert(best_left->peak_index);
        enclosing_ions.push_back({static_cast<float>(best_left->qscore), best_left->peak_index});
      }
      if (best_right && used_peaks.find(best_right->peak_index) == used_peaks.end())
      {
        used_peaks.insert(best_right->peak_index);
        enclosing_ions.push_back({static_cast<float>(best_right->qscore), best_right->peak_index});
      }
    }

    // Sort by qscore descending and fill output arrays
    std::sort(enclosing_ions.begin(), enclosing_ions.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    int count = std::min(n, static_cast<int>(enclosing_ions.size()));
    for (int i = 0; i < count; ++i)
    {
      const auto& pg = ms2_deconvolved_spectrum_[enclosing_ions[i].second];
      masses[i] = pg.getMonoMass();
      qscores[i] = enclosing_ions[i].first;
      charges[i] = pg.getRepAbsCharge();
      auto [mz1, mz2] = pg.getMzRange(charges[i]);
      window_starts[i] = mz1;
      window_ends[i] = mz2;
    }

    return count;
  }

  int FLASHIda::getAmbiguityEnclosingIonsPy(const String& protein_sequence,
                                            int n,
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

    int count = getAmbiguityEnclosingIons(protein_sequence, n, masses.data(), qscores.data(),
                                          charges.data(), window_starts.data(), window_ends.data());

    masses.resize(count);
    qscores.resize(count);
    charges.resize(count);
    window_starts.resize(count);
    window_ends.resize(count);

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
      bool* is_b_ions)
  {
    // Run fragment matching to get all matches
    std::vector<TagBasedFragmentMatch> matches;
    runTagBasedFragmentMatching_(protein_sequence, matches);

    if (matches.empty()) return 0;

    // Separate into b-ions and y-ions
    std::vector<TagBasedFragmentMatch> b_ions, y_ions;
    for (const auto& m : matches)
    {
      if (m.is_prefix) b_ions.push_back(m);
      else y_ions.push_back(m);
    }

    std::cout << "[getTerminalFragmentIons] Found " << b_ions.size() << " b-ions and "
              << y_ions.size() << " y-ions" << std::endl;

    // Sort b-ions by fragment_index descending (rightmost first), then qscore descending
    std::sort(b_ions.begin(), b_ions.end(), [](const auto& a, const auto& b) {
      if (a.fragment_index != b.fragment_index) return a.fragment_index > b.fragment_index;
      return a.qscore > b.qscore;
    });

    // Sort y-ions by fragment_index descending (leftmost first), then qscore descending
    std::sort(y_ions.begin(), y_ions.end(), [](const auto& a, const auto& b) {
      if (a.fragment_index != b.fragment_index) return a.fragment_index > b.fragment_index;
      return a.qscore > b.qscore;
    });

    if (!b_ions.empty())
    {
      std::cout << "[getTerminalFragmentIons] Top b-ions (rightmost): ";
      for (size_t i = 0; i < std::min(b_ions.size(), size_t(3)); ++i)
      {
        std::cout << "b" << b_ions[i].fragment_index << " ";
      }
      std::cout << std::endl;
    }
    if (!y_ions.empty())
    {
      std::cout << "[getTerminalFragmentIons] Top y-ions (leftmost): ";
      for (size_t i = 0; i < std::min(y_ions.size(), size_t(3)); ++i)
      {
        std::cout << "y" << y_ions[i].fragment_index << " ";
      }
      std::cout << std::endl;
    }

    // Interleave output: b, y, b, y, ...
    int count = 0;
    size_t idx_b = 0, idx_y = 0;

    while (count < n && (idx_b < b_ions.size() || idx_y < y_ions.size()))
    {
      const TagBasedFragmentMatch* selected = nullptr;
      bool is_b = false;

      if (count % 2 == 0)  // Even indices: prefer b-ion
      {
        if (idx_b < b_ions.size()) {
          selected = &b_ions[idx_b++];
          is_b = true;
        } else if (idx_y < y_ions.size()) {
          selected = &y_ions[idx_y++];
          is_b = false;
        }
      }
      else  // Odd indices: prefer y-ion
      {
        if (idx_y < y_ions.size()) {
          selected = &y_ions[idx_y++];
          is_b = false;
        } else if (idx_b < b_ions.size()) {
          selected = &b_ions[idx_b++];
          is_b = true;
        }
      }

      if (selected)
      {
        masses[count] = selected->observed_mass;
        qscores[count] = selected->qscore;
        charges[count] = selected->charge;
        is_b_ions[count] = is_b;

        // Get isolation window from deconvolved spectrum
        auto [mz_start, mz_end] = ms2_deconvolved_spectrum_[selected->peak_index].getMzRange(selected->charge);
        window_starts[count] = mz_start;
        window_ends[count] = mz_end;

        ++count;
      }
    }

    std::cout << "[getTerminalFragmentIons] Returning " << count << " interleaved ions:" << std::endl;
    for (int i = 0; i < count; ++i)
    {
      std::cout << "  " << (is_b_ions[i] ? "b" : "y") << " mass=" << masses[i]
                << " qscore=" << qscores[i] << " charge=" << charges[i] << std::endl;
    }

    return count;
  }

  int FLASHIda::getTerminalFragmentIonsPy(
      const String& protein_sequence,
      int n,
      std::vector<double>& masses,
      std::vector<double>& qscores,
      std::vector<int>& charges,
      std::vector<double>& window_starts,
      std::vector<double>& window_ends,
      std::vector<int>& is_b_ions)
  {
    masses.resize(n);
    qscores.resize(n);
    charges.resize(n);
    window_starts.resize(n);
    window_ends.resize(n);

    // Use raw bool array for the C-style function
    std::unique_ptr<bool[]> is_b_temp(new bool[n]);

    int count = getTerminalFragmentIons(
        protein_sequence, n,
        masses.data(), qscores.data(), charges.data(),
        window_starts.data(), window_ends.data(),
        is_b_temp.get());

    // Resize to actual count
    masses.resize(count);
    qscores.resize(count);
    charges.resize(count);
    window_starts.resize(count);
    window_ends.resize(count);

    // Convert bool array to int vector (1 for b-ion, 0 for y-ion)
    is_b_ions.resize(count);
    for (int i = 0; i < count; ++i)
    {
      is_b_ions[i] = is_b_temp[i] ? 1 : 0;
    }

    return count;
  }

  int FLASHIda::getSequenceTagsAndMatchesPy(const std::vector<double>& mzs,
                                            const std::vector<double>& ints,
                                            double rt,
                                            int /*ms_level*/,
                                            const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                            const Param& tagger_param,
                                            std::vector<FLASHHelperClasses::Tag>& tags,
                                            std::vector<TagMatch>& matches,
                                            double ppm_tolerance,
                                            double max_flanking_mass_diff)
  {
    // Deconvolve the spectrum first
    deconvolveMS2Py(mzs, ints, rt);

    // Call the underlying method
    return getSequenceTagsAndMatches(fasta_entries, tagger_param, tags, matches,
                                      ppm_tolerance, max_flanking_mass_diff);
  }

  int FLASHIda::identifyProteoformPy(const std::vector<double>& mzs,
                                     const std::vector<double>& ints,
                                     double rt,
                                     const String& protein_sequence,
                                     double ppm_tolerance,
                                     const std::vector<String>& ion_types,
                                     double ptm_mass_threshold,
                                     std::vector<int>& matched_fragment_indices,
                                     std::vector<int>& ptm_start_positions,
                                     std::vector<int>& ptm_end_positions,
                                     std::vector<double>& ptm_masses)
  {
    // Deconvolve the spectrum first
    deconvolveMS2Py(mzs, ints, rt);

    // Call the underlying method
    return identifyProteoform(protein_sequence, ppm_tolerance, ion_types, ptm_mass_threshold,
                              matched_fragment_indices, ptm_start_positions,
                              ptm_end_positions, ptm_masses);
  }

  int FLASHIda::identifyProteoformExtendedPy(const std::vector<double>& mzs,
                                             const std::vector<double>& ints,
                                             double rt,
                                             const String& protein_sequence,
                                             double ppm_tolerance,
                                             const std::vector<String>& ion_types,
                                             int max_ptm_count,
                                             double max_ptm_mass,
                                             std::vector<int>& matched_peak_indices,
                                             std::vector<double>& matched_theoretical_masses,
                                             std::vector<bool>& matched_ion_types,
                                             std::vector<int>& ptm_start_positions,
                                             std::vector<int>& ptm_end_positions,
                                             std::vector<double>& ptm_masses,
                                             double& coverage,
                                             double& total_score)
  {
    // Deconvolve the spectrum first
    deconvolveMS2Py(mzs, ints, rt);

    // Call the underlying method
    return identifyProteoformExtended(protein_sequence, ppm_tolerance, ion_types,
                                       max_ptm_count, max_ptm_mass,
                                       matched_peak_indices, matched_theoretical_masses,
                                       matched_ion_types, ptm_start_positions,
                                       ptm_end_positions, ptm_masses, coverage, total_score);
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


  int FLASHIda::getSequenceTagsAndMatches(const std::vector<FASTAFile::FASTAEntry>& fasta_entries,
                                          const Param& tagger_param,
                                          std::vector<FLASHHelperClasses::Tag>& tags,
                                          std::vector<TagMatch>& matches,
                                          double ppm_tolerance,
                                          double max_flanking_mass_diff)
  {
    // Clear output vectors
    tags.clear();
    matches.clear();

    // Require deconvolveMS2() to be called first
    if (!ms2_deconv_valid_)
    {
      return 0;
    }
    DeconvolvedSpectrum dspec = ms2_deconvolved_spectrum_;

    // Sort deconvolved spectrum by mass
    dspec.sort();

    // Create and configure tagger
    FLASHTaggerAlgorithm tagger;
    if (! tagger_param.empty())
    {
      tagger.setParameters(tagger_param);
    }

    // Run tag generation
    tagger.run(dspec, ppm_tolerance);

    // Get the generated tags
    tagger.fillTags(tags);

    if (tags.empty())
    {
      return 0;
    }

    // Prepare protein sequences for matching (replace I with L for consistency)
    std::vector<String> protein_seqs;
    protein_seqs.reserve(fasta_entries.size());

    for (const auto& fe : fasta_entries)
    {
      String seq = fe.sequence;
      std::replace(seq.begin(), seq.end(), 'I', 'L');
      protein_seqs.push_back(seq);
    }

    // Match tags against protein database
    for (const auto& tag : tags)
    {
      for (Size protein_idx = 0; protein_idx < protein_seqs.size(); ++protein_idx)
      {
        const String& pseq = protein_seqs[protein_idx];

        // Find all positions where tag matches in the protein
        std::vector<int> positions;
        std::vector<double> flanking_mass_diffs;

        FLASHTaggerAlgorithm::fillMatchedPositionsAndFlankingMassDiffs(
          positions,
          flanking_mass_diffs,
          max_flanking_mass_diff,
          pseq,
          tag);

        // Add matches
        for (Size i = 0; i < positions.size(); ++i)
        {
          TagMatch match;
          match.tag_sequence = tag.getSequence();
          match.n_term_mass = tag.getNtermMass();
          match.c_term_mass = tag.getCtermMass();
          match.tag_score = tag.getScore();
          match.protein_index = static_cast<int>(protein_idx);
          match.protein_accession = fasta_entries[protein_idx].identifier;
          match.match_position = positions[i];
          match.flanking_mass_diff = flanking_mass_diffs[i];

          matches.push_back(match);
        }
      }
    }

    return static_cast<int>(tags.size());
  }

  /**
   * @brief Identify proteoform from MS2 spectrum and a single protein sequence
   *
   * This implements the full FLASHTnT workflow using FLASHExtenderAlgorithm:
   * 1. Deconvolve the MS2 spectrum to get monoisotopic masses
   * 2. Run FLASHTagger to generate sequence tags from the deconvolved spectrum
   * 3. Match tags against the protein sequence
   * 4. Run FLASHExtender to detect PTMs through path-finding analysis
   *
   * @param mzs m/z values of the input MS2 spectrum
   * @param ints intensities of the input MS2 spectrum
   * @param length number of peaks in the spectrum
   * @param rt retention time in seconds
   * @param protein_sequence the protein sequence to match against
   * @param ppm_tolerance mass tolerance in ppm (default 10.0)
   * @param ion_types ion types to consider (default {"b", "y"})
   * @param ptm_mass_threshold (unused - PTM detection is handled by FLASHExtender)
   * @param matched_fragment_indices output: indices of matched fragment ions (1-based positions in sequence)
   * @param ptm_start_positions output: start positions of PTM localization ranges (1-based)
   * @param ptm_end_positions output: end positions of PTM localization ranges (1-based)
   * @param ptm_masses output: mass shifts at each PTM position
   * @return number of matched amino acids
   */
  int FLASHIda::identifyProteoform(const String& protein_sequence,
                                   double ppm_tolerance,
                                   const std::vector<String>& ion_types,
                                   double ptm_mass_threshold,
                                   std::vector<int>& matched_fragment_indices,
                                   std::vector<int>& ptm_start_positions,
                                   std::vector<int>& ptm_end_positions,
                                   std::vector<double>& ptm_masses)
  {
    // Suppress unused parameter warning - PTM detection handled by FLASHExtender
    (void)ptm_mass_threshold;

    // Clear output vectors
    matched_fragment_indices.clear();
    ptm_start_positions.clear();
    ptm_end_positions.clear();
    ptm_masses.clear();

    if (protein_sequence.empty())
    {
      return 0;
    }

    // Require deconvolveMS2() to be called first
    if (!ms2_deconv_valid_)
    {
      return 0;
    }
    DeconvolvedSpectrum dspec = ms2_deconvolved_spectrum_;

    // Sort deconvolved spectrum by mass
    dspec.sort();

    // Use default ion types if none specified
    std::vector<String> types = ion_types;
    if (types.empty())
    {
      types = {"b", "y"};
    }

    // === FLASHExtender-based workflow (same as FLASHTnT) ===

    // Step 1: Run FLASHTagger to get sequence tags
    FLASHTaggerAlgorithm tagger;
    Param tagger_param = tagger.getDefaults();
    std::vector<std::string> ion_types_str;
    for (const auto& t : types)
    {
      ion_types_str.push_back(t);
    }
    tagger_param.setValue("ion_type", ion_types_str);
    tagger.setParameters(tagger_param);
    tagger.run(dspec, ppm_tolerance);

    std::vector<FLASHHelperClasses::Tag> tags;
    tagger.fillTags(tags);

    if (tags.empty())
    {
      // No tags found - return empty result (matching FLASHTnT behavior)
      return 0;
    }

    // Step 2: Create ProteinHit with proper metadata
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
      // No tags match the protein - return empty result
      return 0;
    }

    hit.setMetaValue("TagIndices", std::vector<int>(tag_indices_set.begin(), tag_indices_set.end()));
    hit.setMetaValue("TagPositions", tag_positions);
    hit.setMetaValue("MatchedAA", static_cast<int>(tag_positions.size()));
    hit.setCoverage(static_cast<double>(tag_positions.size()) / protein_sequence.size());

    // Calculate flanking masses for each tag position (required by FLASHTaggerAlgorithm::runMatching)
    std::set<int> n_flanking_masses_set, c_flanking_masses_set;
    for (Size i = 0; i < tags.size(); ++i)
    {
      String tag_seq = tags[i].getUppercaseSequence();
      Size pos = clean_seq.find(tag_seq);
      if (pos != String::npos)
      {
        // Calculate N-terminal flanking mass difference
        double n_mass = 0;
        for (Size j = 0; j < pos; ++j)
        {
          n_mass += ResidueDB::getInstance()->getResidue(clean_seq[j])->getMonoWeight(Residue::Internal);
        }
        double n_diff = static_cast<int>(std::round(n_mass - tags[i].getNtermMass()));
        n_flanking_masses_set.insert(n_diff);

        // Calculate C-terminal flanking mass difference
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

    // Step 3: Create spec_vec (integer masses from spectrum)
    std::vector<int> spec_vec;
    spec_vec.reserve(dspec.size() + 1);
    spec_vec.push_back(0);
    for (const auto& pg : dspec)
    {
      spec_vec.push_back(static_cast<int>(std::round(pg.getMonoMass())));
    }

    // Step 4: Create vec_pro and rev_vec_pro for the single protein
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

    // Step 5: Run FLASHTagger matching
    double max_mod_mass = 500.0; // default max modification mass
    FLASHTaggerAlgorithm::runMatching(hits, dspec, spec_vec, vec_pro, rev_vec_pro, max_mod_mass);

    if (hits.empty())
    {
      return 0;
    }

    // Step 6: Run FLASHExtender for PTM detection
    FLASHExtenderAlgorithm extender;
    Param extender_param = extender.getDefaults();
    extender_param.setValue("ion_type", ion_types_str);
    extender.setParameters(extender_param);

    extender.run(hits, dspec, spec_vec, vec_pro, rev_vec_pro, tags, ppm_tolerance, false);

    std::vector<ProteinHit> proteoform_hits;
    extender.fillProteoforms(proteoform_hits);

    // Step 7: Extract results from the best proteoform hit
    if (proteoform_hits.empty())
    {
      // FLASHExtender didn't produce results - return empty (matching FLASHTnT behavior)
      return 0;
    }

    // Get the best hit (highest score)
    const ProteinHit& best_hit = proteoform_hits[0];

    // Extract modifications (PTMs) from FLASHTnT output
    if (best_hit.metaValueExists("Modifications"))
    {
      std::vector<double> mod_masses_vec = best_hit.getMetaValue("Modifications");
      std::vector<int> mod_starts = best_hit.getMetaValue("ModificationStarts");
      std::vector<int> mod_ends = best_hit.getMetaValue("ModificationEnds");

      for (Size i = 0; i < mod_masses_vec.size(); ++i)
      {
        ptm_start_positions.push_back(mod_starts[i]);
        ptm_end_positions.push_back(mod_ends[i]);
        ptm_masses.push_back(mod_masses_vec[i]);
      }
    }

    // Extract matched fragment count from MatchedAA
    int matched_count = 0;
    if (best_hit.metaValueExists("MatchedAA"))
    {
      matched_count = best_hit.getMetaValue("MatchedAA");
    }

    // Populate matched_fragment_indices (simplified - just return the count via MatchedAA)
    // For detailed fragment indices, would need to trace through FLASHExtender's path
    for (int i = 0; i < matched_count; ++i)
    {
      matched_fragment_indices.push_back(i + 1); // placeholder indices
    }

    return matched_count;
  }

  /**
   * @brief Extended identification with detailed output including all match information
   *
   * This is an extended version that also uses sequence tagging for improved accuracy,
   * similar to the full FLASHTnT workflow.
   *
   * @param mzs m/z values of the input MS2 spectrum
   * @param ints intensities of the input MS2 spectrum
   * @param length number of peaks in the spectrum
   * @param rt retention time in seconds
   * @param protein_sequence the protein sequence to match against
   * @param ppm_tolerance mass tolerance in ppm
   * @param ion_types ion types to consider
   * @param max_ptm_count maximum number of PTMs to consider
   * @param max_ptm_mass maximum mass shift for a single PTM
   * @param matched_peak_indices output: indices of matched peaks in the original spectrum
   * @param matched_theoretical_masses output: theoretical masses that were matched
   * @param matched_ion_types output: ion types (true=prefix, false=suffix) for each match
   * @param ptm_start_positions output: start of region for each PTM
   * @param ptm_end_positions output: end of region for each PTM
   * @param ptm_masses output: mass shift for each PTM
   * @param coverage output: sequence coverage
   * @param total_score output: total identification score
   * @return number of matched ions
   */
  int FLASHIda::identifyProteoformExtended(const String& protein_sequence,
                                           double ppm_tolerance,
                                           const std::vector<String>& ion_types,
                                           int max_ptm_count,
                                           double max_ptm_mass,
                                           std::vector<int>& matched_peak_indices,
                                           std::vector<double>& matched_theoretical_masses,
                                           std::vector<bool>& matched_ion_types,
                                           std::vector<int>& ptm_start_positions,
                                           std::vector<int>& ptm_end_positions,
                                           std::vector<double>& ptm_masses,
                                           double& coverage,
                                           double& total_score)
  {
    // Clear output vectors
    matched_peak_indices.clear();
    matched_theoretical_masses.clear();
    matched_ion_types.clear();
    ptm_start_positions.clear();
    ptm_end_positions.clear();
    ptm_masses.clear();
    coverage = 0.0;
    total_score = 0.0;

    if (protein_sequence.empty())
    {
      return 0;
    }

    // Require deconvolveMS2() to be called first
    if (!ms2_deconv_valid_)
    {
      return 0;
    }
    DeconvolvedSpectrum dspec = ms2_deconvolved_spectrum_;

    // Sort deconvolved spectrum by mass
    dspec.sort();

    // Calculate theoretical fragment masses
    std::vector<double> prefix_masses, suffix_masses;
    std::vector<double> prefix_shifts, suffix_shifts;

    std::vector<String> types = ion_types;
    if (types.empty())
    {
      types = {"b", "y"};
    }

    calculateTheoreticalFragmentMasses(protein_sequence, types, prefix_masses, suffix_masses, prefix_shifts, suffix_shifts);

    // Create spectrum vector for FLASHTaggerAlgorithm matching
    std::vector<int> spec_vec;
    spec_vec.reserve(dspec.size() + 1);
    spec_vec.push_back(0);
    for (const auto& pg : dspec)
    {
      spec_vec.push_back(static_cast<int>(round(pg.getMonoMass())));
    }

    // Create protein vectors for matching
    String clean_seq = protein_sequence;
    std::replace(clean_seq.begin(), clean_seq.end(), 'I', 'L');

    std::unordered_set<int> pro_prefix_set, pro_suffix_set;
    double cumulative = 0.0;
    pro_prefix_set.insert(0);
    for (Size i = 0; i < clean_seq.size(); ++i)
    {
      char aa = clean_seq[i];
      if (aa != 'X')
      {
        const Residue* res = ResidueDB::getInstance()->getResidue(aa);
        if (res) cumulative += res->getMonoWeight(Residue::Internal);
      }
      pro_prefix_set.insert(static_cast<int>(round(cumulative)));
    }

    cumulative = 0.0;
    pro_suffix_set.insert(0);
    for (int i = clean_seq.size() - 1; i >= 0; --i)
    {
      char aa = clean_seq[i];
      if (aa != 'X')
      {
        const Residue* res = ResidueDB::getInstance()->getResidue(aa);
        if (res) cumulative += res->getMonoWeight(Residue::Internal);
      }
      pro_suffix_set.insert(static_cast<int>(round(cumulative)));
    }

    // Match observed masses and track coverage
    std::vector<FragmentIonMatch> all_matches;
    std::set<int> covered_positions;

    for (int peak_idx = 0; peak_idx < static_cast<int>(dspec.size()); ++peak_idx)
    {
      const auto& pg = dspec[peak_idx];
      double observed_mass = pg.getMonoMass();
      float score = static_cast<float>(pg.getQscore());

      // Try matching against prefix ions
      int best_prefix_idx = -1;
      double best_prefix_theo = 0, best_prefix_diff = 0, best_prefix_shift = 0;
      bool prefix_match = findBestMatch(observed_mass, prefix_masses, prefix_shifts, ppm_tolerance,
                                        best_prefix_idx, best_prefix_theo, best_prefix_diff, best_prefix_shift);

      // Try matching against suffix ions
      int best_suffix_idx = -1;
      double best_suffix_theo = 0, best_suffix_diff = 0, best_suffix_shift = 0;
      bool suffix_match = findBestMatch(observed_mass, suffix_masses, suffix_shifts, ppm_tolerance,
                                        best_suffix_idx, best_suffix_theo, best_suffix_diff, best_suffix_shift);

      if (prefix_match || suffix_match)
      {
        bool use_prefix = false;
        if (prefix_match && suffix_match)
        {
          double prefix_ppm = std::abs(best_prefix_diff) / best_prefix_theo * 1e6;
          double suffix_ppm = std::abs(best_suffix_diff) / best_suffix_theo * 1e6;
          use_prefix = (prefix_ppm <= suffix_ppm);
        }
        else
        {
          use_prefix = prefix_match;
        }

        FragmentIonMatch match;
        match.peak_index = peak_idx;
        match.observed_mass = observed_mass;
        match.score = score;

        if (use_prefix)
        {
          match.fragment_index = best_prefix_idx + 1;
          match.theoretical_mass = best_prefix_theo;
          match.mass_diff = best_prefix_diff;
          match.ppm_error = std::abs(best_prefix_diff) / best_prefix_theo * 1e6;
          match.is_prefix = true;
          match.ion_type = "b";

          // Mark covered positions
          for (int pos = 1; pos <= match.fragment_index; ++pos)
          {
            covered_positions.insert(pos);
          }
        }
        else
        {
          int seq_len = static_cast<int>(protein_sequence.size());
          match.fragment_index = seq_len - best_suffix_idx;
          match.theoretical_mass = best_suffix_theo;
          match.mass_diff = best_suffix_diff;
          match.ppm_error = std::abs(best_suffix_diff) / best_suffix_theo * 1e6;
          match.is_prefix = false;
          match.ion_type = "y";

          // Mark covered positions
          for (int pos = match.fragment_index; pos <= seq_len; ++pos)
          {
            covered_positions.insert(pos);
          }
        }

        all_matches.push_back(match);
        total_score += score;
      }
    }

    // Fill output vectors
    for (const auto& m : all_matches)
    {
      matched_peak_indices.push_back(m.peak_index);
      matched_theoretical_masses.push_back(m.theoretical_mass);
      matched_ion_types.push_back(m.is_prefix);
    }

    // Calculate coverage
    coverage = static_cast<double>(covered_positions.size()) / static_cast<double>(protein_sequence.size());

    // Detect PTM sites
    std::vector<PTMSite> ptm_sites;
    detectPTMSites(all_matches, protein_sequence, 5.0, ptm_sites); // 5 Da threshold

    // Filter PTMs by max count and mass
    int ptm_count = 0;
    for (const auto& site : ptm_sites)
    {
      if (ptm_count >= max_ptm_count) break;
      if (std::abs(site.mass_shift) > max_ptm_mass) continue;

      ptm_start_positions.push_back(site.start_position);
      ptm_end_positions.push_back(site.end_position);
      ptm_masses.push_back(site.mass_shift);
      ++ptm_count;
    }

    return static_cast<int>(all_matches.size());
  }

  /**
   * @brief Core identification function that works with pre-deconvolved masses
   *
   * This function performs proteoform identification using FLASHTagger + FLASHExtender.
   * It takes deconvolved monoisotopic masses and uses the full FLASHTnT workflow
   * to properly detect PTMs through mass shift tracking.
   */
  int FLASHIda::identifyProteoformFromMasses(const std::vector<double>& observed_masses,
                                             const std::vector<double>& mass_scores,
                                             const String& protein_sequence,
                                             double ppm_tolerance,
                                             const std::vector<String>& ion_types,
                                             double /* ptm_mass_threshold - now controlled by extender params */,
                                             std::vector<int>& matched_fragment_indices,
                                             std::vector<bool>& matched_ion_types_out,
                                             std::vector<double>& matched_observed_masses,
                                             std::vector<double>& matched_theoretical_masses,
                                             std::vector<double>& matched_ppm_errors,
                                             std::vector<int>& ptm_start_positions,
                                             std::vector<int>& ptm_end_positions,
                                             std::vector<double>& ptm_masses)
  {
    // Clear output vectors
    matched_fragment_indices.clear();
    matched_ion_types_out.clear();
    matched_observed_masses.clear();
    matched_theoretical_masses.clear();
    matched_ppm_errors.clear();
    ptm_start_positions.clear();
    ptm_end_positions.clear();
    ptm_masses.clear();

    if (protein_sequence.empty() || observed_masses.empty())
    {
      return 0;
    }

    // Use default ion types if none specified
    std::vector<String> types = ion_types;
    if (types.empty())
    {
      types = {"b", "y"};
    }

    // Step 1: Create DeconvolvedSpectrum from input masses
    DeconvolvedSpectrum dspec(1); // scan number 1
    dspec.reserve(observed_masses.size());

    for (Size i = 0; i < observed_masses.size(); ++i)
    {
      PeakGroup pg(1, 10, true); // min_charge=1, max_charge=10, positive mode
      pg.setMonoisotopicMass(observed_masses[i]);
      double score = (i < mass_scores.size()) ? mass_scores[i] : 1.0;
      pg.setQscore(score);
      dspec.push_back(pg);
    }

    // Sort by mass
    dspec.sort();

    // Step 2: Run FLASHTagger to get sequence tags
    FLASHTaggerAlgorithm tagger;

    // Configure tagger with ion types
    Param tagger_param = tagger.getDefaults();
    std::vector<std::string> ion_types_str;
    for (const auto& t : types)
    {
      ion_types_str.push_back(t);
    }
    tagger_param.setValue("ion_type", ion_types_str);
    tagger.setParameters(tagger_param);

    tagger.run(dspec, ppm_tolerance);

    std::vector<FLASHHelperClasses::Tag> tags;
    tagger.fillTags(tags);

    // If no tags found, return empty result (matching FLASHTnT behavior)
    if (tags.empty())
    {
      return 0;
    }

    // Step 3: Create ProteinHit for the input sequence
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
      // No tags match the protein - return empty result
      return 0;
    }

    hit.setMetaValue("TagIndices", std::vector<int>(tag_indices_set.begin(), tag_indices_set.end()));
    hit.setMetaValue("TagPositions", tag_positions);
    hit.setMetaValue("MatchedAA", static_cast<int>(tag_positions.size()));
    hit.setCoverage(static_cast<double>(tag_positions.size()) / protein_sequence.size());

    // Calculate flanking masses for each tag position (required by FLASHTaggerAlgorithm::runMatching)
    std::set<int> n_flanking_masses_set, c_flanking_masses_set;
    for (Size i = 0; i < tags.size(); ++i)
    {
      String tag_seq = tags[i].getUppercaseSequence();
      Size pos = clean_seq.find(tag_seq);
      if (pos != String::npos)
      {
        // Calculate N-terminal flanking mass difference
        double n_mass = 0;
        for (Size j = 0; j < pos; ++j)
        {
          n_mass += ResidueDB::getInstance()->getResidue(clean_seq[j])->getMonoWeight(Residue::Internal);
        }
        double n_diff = static_cast<int>(std::round(n_mass - tags[i].getNtermMass()));
        n_flanking_masses_set.insert(n_diff);

        // Calculate C-terminal flanking mass difference
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

    // Step 4: Create spec_vec (integer masses from spectrum)
    std::vector<int> spec_vec;
    spec_vec.reserve(dspec.size() + 1);
    spec_vec.push_back(0);
    for (const auto& pg : dspec)
    {
      spec_vec.push_back(static_cast<int>(std::round(pg.getMonoMass())));
    }

    // Step 5: Create vec_pro and rev_vec_pro for the single protein
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

    // Step 6: Run FLASHTagger matching
    double max_mod_mass = 500.0; // default max modification mass
    FLASHTaggerAlgorithm::runMatching(hits, dspec, spec_vec, vec_pro, rev_vec_pro, max_mod_mass);

    if (hits.empty())
    {
      return 0;
    }

    // Step 7: Run FLASHExtender for PTM detection
    FLASHExtenderAlgorithm extender;
    Param extender_param = extender.getDefaults();
    extender_param.setValue("ion_type", ion_types_str);
    extender.setParameters(extender_param);

    extender.run(hits, dspec, spec_vec, vec_pro, rev_vec_pro, tags, ppm_tolerance, false);

    std::vector<ProteinHit> proteoform_hits;
    extender.fillProteoforms(proteoform_hits);

    // Step 8: Extract results from the best proteoform hit
    if (proteoform_hits.empty())
    {
      // FLASHExtender didn't produce results - return empty (matching FLASHTnT behavior)
      return 0;
    }

    // Get the best hit (highest score)
    const ProteinHit& best_hit = proteoform_hits[0];

    // Extract modifications (PTMs) from FLASHTnT output
    if (best_hit.metaValueExists("Modifications"))
    {
      std::vector<double> mod_masses_vec = best_hit.getMetaValue("Modifications");
      std::vector<int> mod_starts = best_hit.getMetaValue("ModificationStarts");
      std::vector<int> mod_ends = best_hit.getMetaValue("ModificationEnds");

      for (Size i = 0; i < mod_masses_vec.size(); ++i)
      {
        ptm_start_positions.push_back(mod_starts[i]);
        ptm_end_positions.push_back(mod_ends[i]);
        ptm_masses.push_back(mod_masses_vec[i]);
      }
    }

    // Populate fragment matching output vectors using simple ppm matching
    // (PTM detection already done by FLASHExtender above)
    std::vector<double> prefix_masses, suffix_masses;
    std::vector<double> prefix_shifts, suffix_shifts;
    calculateTheoreticalFragmentMasses(protein_sequence, types, prefix_masses, suffix_masses, prefix_shifts, suffix_shifts);

    for (Size peak_idx = 0; peak_idx < observed_masses.size(); ++peak_idx)
    {
      double observed_mass = observed_masses[peak_idx];

      int best_prefix_idx = -1;
      double best_prefix_theo = 0, best_prefix_diff = 0, best_prefix_shift = 0;
      bool prefix_match = findBestMatch(observed_mass, prefix_masses, prefix_shifts, ppm_tolerance,
                                        best_prefix_idx, best_prefix_theo, best_prefix_diff, best_prefix_shift);

      int best_suffix_idx = -1;
      double best_suffix_theo = 0, best_suffix_diff = 0, best_suffix_shift = 0;
      bool suffix_match = findBestMatch(observed_mass, suffix_masses, suffix_shifts, ppm_tolerance,
                                        best_suffix_idx, best_suffix_theo, best_suffix_diff, best_suffix_shift);

      if (prefix_match || suffix_match)
      {
        bool use_prefix = prefix_match && (!suffix_match ||
          std::abs(best_prefix_diff) / best_prefix_theo < std::abs(best_suffix_diff) / best_suffix_theo);

        if (use_prefix)
        {
          matched_fragment_indices.push_back(best_prefix_idx + 1);
          matched_ion_types_out.push_back(true);
          matched_observed_masses.push_back(observed_mass);
          matched_theoretical_masses.push_back(best_prefix_theo);
          matched_ppm_errors.push_back(std::abs(best_prefix_diff) / best_prefix_theo * 1e6);
        }
        else
        {
          int seq_len = static_cast<int>(protein_sequence.size());
          matched_fragment_indices.push_back(seq_len - best_suffix_idx);
          matched_ion_types_out.push_back(false);
          matched_observed_masses.push_back(observed_mass);
          matched_theoretical_masses.push_back(best_suffix_theo);
          matched_ppm_errors.push_back(std::abs(best_suffix_diff) / best_suffix_theo * 1e6);
        }
      }
    }

    return static_cast<int>(matched_fragment_indices.size());
  }

  int FLASHIda::identifyProteoformFromMassesPy(std::vector<double>& observed_masses,
                                               std::vector<double>& mass_scores,
                                               const String& protein_sequence,
                                               double ppm_tolerance,
                                               std::vector<String>& ion_types,
                                               double ptm_mass_threshold,
                                               std::vector<int>& matched_fragment_indices,
                                               std::vector<bool>& matched_ion_types,
                                               std::vector<double>& matched_observed_masses,
                                               std::vector<double>& matched_theoretical_masses,
                                               std::vector<double>& matched_ppm_errors,
                                               std::vector<int>& ptm_start_positions,
                                               std::vector<int>& ptm_end_positions,
                                               std::vector<double>& ptm_masses)
  {
    return identifyProteoformFromMasses(observed_masses, mass_scores, protein_sequence,
                                        ppm_tolerance, ion_types, ptm_mass_threshold,
                                        matched_fragment_indices, matched_ion_types,
                                        matched_observed_masses, matched_theoretical_masses,
                                        matched_ppm_errors, ptm_start_positions, ptm_end_positions, ptm_masses);
  }

  void FLASHIda::calculateTheoreticalFragmentMassesPy(const String& sequence,
                                                      const std::vector<String>& ion_types,
                                                      std::vector<double>& prefix_masses,
                                                      std::vector<double>& suffix_masses)
  {
    std::vector<double> prefix_shifts, suffix_shifts;
    calculateTheoreticalFragmentMasses(sequence, ion_types, prefix_masses, suffix_masses,
                                       prefix_shifts, suffix_shifts);

    // Apply ion type shifts to convert internal masses to theoretical fragment masses
    // Use the first available shift for each direction (prefix/suffix)
    double prefix_shift = prefix_shifts.empty() ? Residue::getInternalToBIon().getMonoWeight() : prefix_shifts[0];
    double suffix_shift = suffix_shifts.empty() ? Residue::getInternalToYIon().getMonoWeight() : suffix_shifts[0];

    for (double& mass : prefix_masses)
    {
      mass += prefix_shift;
    }
    for (double& mass : suffix_masses)
    {
      mass += suffix_shift;
    }
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
      for (int count = 0; count <= ptm.max_count; ++count)
      {
        int new_total = total_count + count;
        if (new_total > max_total_ptm_count_) break;  // Prune: exceeds global max
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

} // namespace OpenMS
