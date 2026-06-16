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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FragmentAnalysis.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHExtenderAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHTaggerAlgorithm.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_set>

namespace OpenMS
{

// Anonymous namespace for fragment analysis helper functions
namespace
{
  /// optimal window margin
  inline const double optimal_window_margin_ = .4;

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

  /// Check if ion type is a prefix (N-terminal) ion
  inline bool isPrefixIon(char ion_type)
  {
    return ion_type == 'a' || ion_type == 'b' || ion_type == 'c';
  }

  using PTMSite = FragmentAnalysis::PTMSite;


  /// Calculate theoretical fragment masses with PTM adjustments from FLASHExtender
  /// For each PTM at positions [start, end], fragments containing that region get the PTM mass added
  void calculatePTMAdjustedFragmentMasses(
      const String& sequence,
      const std::vector<PTMSite>& ptm_sites,
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
      const std::vector<PTMSite>& ptm_sites,
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

    // Calculate base cumulative mass for prefix direction (N->C)
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

    // Calculate base cumulative mass for suffix direction (C->N)
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


  std::vector<std::string> FragmentAnalysis::getIonTypesForFragmentationMethod(const String& method)
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

  FragmentAnalysis::FragmentAnalysis(const Config& config)
    : config_(config)
  {
  }

  // Re-compute the fragment precursor's deconvolution scores from the live PeakGroup,
  // mirroring buildMS2 (ScanCommandQueue.cpp) getter-for-getter. Used for the stage-1
  // (fragment) scoring columns of an MS3 scan_commands row.
  FragmentAnalysis::FragmentScores
  FragmentAnalysis::FragmentScores::fromPeakGroup(const PeakGroup& pg, int abs_charge)
  {
    FragmentScores s;
    s.mono_mass            = pg.getMonoMass();
    s.qscore               = pg.getQscore();
    s.charge_cos           = pg.getChargeIsotopeCosine(abs_charge);
    s.charge_snr           = pg.getChargeSNR(abs_charge);
    s.iso_cos              = pg.getIsotopeCosine();
    s.snr                  = pg.getSNR();
    s.charge_score         = pg.getChargeScore();
    s.ppm_error            = pg.getAvgPPMError();
    s.precursor_intensity  = pg.getChargeIntensity(abs_charge);
    s.peakgroup_intensity  = pg.getIntensity();
    return s;
  }

  int FragmentAnalysis::getBestMS2Masses(int n,
                                         double* masses,
                                         double* qscores,
                                         int* charges,
                                         double* window_starts,
                                         double* window_ends,
                                         DeconvolvedSpectrum& stored_ms2)
  {
    if (stored_ms2.empty())
    {
      return 0;
    }

    std::sort(stored_ms2.begin(), stored_ms2.end(),
    [](const PeakGroup& a, const PeakGroup& b) {
      return a.getChargeIntensity(a.getMaxIntensityAbsCharge()) > b.getChargeIntensity(b.getMaxIntensityAbsCharge());
    });

    int output_idx = 0;
    for (Size pg_idx = 0; pg_idx < stored_ms2.size() && output_idx < n; ++pg_idx)
    {
      const auto& pg = stored_ms2[pg_idx];
      double mono_mass = pg.getMonoMass();

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

    return output_idx;
  }

  int FragmentAnalysis::getBestMS2MassesPy(int n,
                                           std::vector<double>& masses,
                                           std::vector<double>& qscores,
                                           std::vector<int>& charges,
                                           std::vector<double>& window_starts,
                                           std::vector<double>& window_ends,
                                           DeconvolvedSpectrum& stored_ms2)
  {
    masses.resize(n);
    qscores.resize(n);
    charges.resize(n);
    window_starts.resize(n);
    window_ends.resize(n);

    int count = getBestMS2Masses(n, masses.data(), qscores.data(), charges.data(),
                                 window_starts.data(), window_ends.data(), stored_ms2);

    masses.resize(count);
    qscores.resize(count);
    charges.resize(count);
    window_starts.resize(count);
    window_ends.resize(count);

    return count;
  }

  int FragmentAnalysis::runTagBasedFragmentMatching_(const String& protein_sequence,
                                                     std::vector<TagBasedFragmentMatch>& matches,
                                                     DeconvolvedSpectrum& stored_ms2,
                                                     ProteoformMatch& result,
                                                     const String& fragmentation_method,
                                                     double tolerance_ppm)
  {
    matches.clear();

    if (stored_ms2.empty() || protein_sequence.empty())
    {
      return 0;
    }

    // 1. Copy and sort deconvolved spectrum
    DeconvolvedSpectrum dspec = stored_ms2;
    dspec.sort();
    double precursor_mass = dspec.getPrecursorPeakGroup().getMonoMass();
    std::cout << "PM=" << precursor_mass << std::endl;
    for (const auto& pg : dspec)
    {
      double mono_mass = pg.getMonoMass();
      std::cout << mono_mass << std::endl;
    }

    // 2. Run FLASHTagger to generate sequence tags
    double ppm_tolerance = (tolerance_ppm > 0.0) ? tolerance_ppm : config_.level(2).tolerance_ppm;
    std::vector<std::string> ion_types_str = FragmentAnalysis::getIonTypesForFragmentationMethod(fragmentation_method);
    FLASHTaggerAlgorithm tagger;
    Param tagger_param = tagger.getDefaults();
    tagger_param.setValue("ion_type", ion_types_str);
    tagger.setParameters(tagger_param);
    tagger.run(dspec, ppm_tolerance);

    std::vector<FLASHHelperClasses::Tag> tags;
    tagger.fillTags(tags);

    // E4: record the FLASHTagger tag count generated during identification so the logged tag_count
    // reflects the real value (reuses the already-generated tags — no extra work, no re-run).
    result.tag_count = static_cast<int>(tags.size());

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
    extender_param.setValue("skip_precursor_inference", "true");
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

    // Populate proteoform region in result
    result.region_start = start_pos;  // 0-based, or -1 if not found
    result.region_end = end_pos;       // 0-based exclusive, or -1
    result.proteoform_sequence = protein_sequence;

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

    // Populate PTM sites in result
    result.ptm_sites.clear();
    for (Size i = 0; i < mod_masses.size(); ++i)
    {
      PTMSite site;
      site.start_position = mod_starts[i];
      site.end_position = mod_ends[i];
      site.position = (site.start_position + site.end_position) / 2;
      site.mass_shift = mod_masses[i];
      result.ptm_sites.push_back(site);
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
    for (Size peak_idx = 0; peak_idx < stored_ms2.size(); ++peak_idx)
    {
      const auto& pg = stored_ms2[peak_idx];
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
        match.intensity = pg.getChargeIntensity(pg.getMaxIntensityAbsCharge());
        match.charge = pg.getMaxIntensityAbsCharge();
        match.fragment_index = best_frag_idx;
        match.ion_type = best_ion_type;
        match.theoretical_mass = best_theo;
        match.ppm_error = best_ppm;
        matches.push_back(match);
      }
    }

    // 11. Sort by intensity descending
    std::sort(matches.begin(), matches.end(),
              [](const TagBasedFragmentMatch& a, const TagBasedFragmentMatch& b) {
                return a.intensity > b.intensity;
              });

    // === Diagnostic Output: Matched fragment ions ===
    std::cout << "[runTagBasedFragmentMatching_] Matched " << matches.size() << " fragment ions (" << dspec.size() << " masses):" << std::endl;
    Size display_count = std::min(Size(10), matches.size());
    for (Size i = 0; i < display_count; ++i)
    {
      const auto& m = matches[i];
      std::cout << "  " << m.ion_type << m.fragment_index
                << " - " << std::fixed << std::setprecision(2) << m.observed_mass << " Da"
                << " (intensity:" << std::setprecision(2) << m.intensity << ")" << std::endl;
    }
    if (matches.size() > 10)
    {
      std::cout << "  ... (" << (matches.size() - 10) << " more)" << std::endl;
    }

    result.total_match_count = static_cast<int>(matches.size());
    result.fragments.clear();
    result.fragments.reserve(matches.size());
    for (const auto& m : matches)
    {
      ProteoformMatch::FragmentMatch fm;
      fm.ion_type = std::string(1, m.ion_type);
      fm.ion_index = m.fragment_index;
      fm.observed_mass = m.observed_mass;
      result.fragments.push_back(std::move(fm));
    }
    return static_cast<int>(matches.size());
  }

  int FragmentAnalysis::getTopFragmentMatches(const String& protein_sequence,
                                              int n,
                                              double* masses,
                                              double* qscores,
                                              int* charges,
                                              double* window_starts,
                                              double* window_ends,
                                              char* ion_types,
                                              int* fragment_indices,
                                              DeconvolvedSpectrum& stored_ms2,
                                              ProteoformMatch& result,
                                              const String& fragmentation_method,
                                              double tolerance_ppm,
                                              FragmentScores* frag_scores)
  {
    std::cout << "Matching fragments!" << std::endl;
    // Use tag-based matching workflow (FLASHTagger + FLASHExtender)
    std::vector<TagBasedFragmentMatch> matches;
    runTagBasedFragmentMatching_(protein_sequence, matches, stored_ms2, result, fragmentation_method, tolerance_ppm);

    int output_idx = 0;
    for (size_t i = 0; i < matches.size() && output_idx < n; ++i)
    {
      const auto& m = matches[i];
      const auto& pg = stored_ms2[m.peak_index];

      // Single charge behavior
      masses[output_idx] = m.observed_mass;
      qscores[output_idx] = m.intensity;
      charges[output_idx] = m.charge;
      ion_types[output_idx] = m.ion_type;
      fragment_indices[output_idx] = m.fragment_index;
      auto [mz1, mz2] = pg.getMzRange(m.charge);
      window_starts[output_idx] = mz1 - optimal_window_margin_;
      window_ends[output_idx] = mz2 + optimal_window_margin_;
      if (frag_scores) frag_scores[output_idx] = FragmentScores::fromPeakGroup(pg, std::abs(m.charge));
      ++output_idx;
    }

    return output_idx;
  }

  int FragmentAnalysis::getTopFragmentMatchesPy(const String& protein_sequence,
                                                int n,
                                                std::vector<double>& masses,
                                                std::vector<double>& qscores,
                                                std::vector<int>& charges,
                                                std::vector<double>& window_starts,
                                                std::vector<double>& window_ends,
                                                std::vector<int>& is_b_ions,
                                                std::vector<int>& fragment_indices,
                                                DeconvolvedSpectrum& stored_ms2)
  {
    masses.resize(n);
    qscores.resize(n);
    charges.resize(n);
    window_starts.resize(n);
    window_ends.resize(n);

    // Use raw char array for the C-style function
    std::unique_ptr<char[]> ion_types_temp(new char[n]);
    fragment_indices.resize(n);

    ProteoformMatch discard;
    int count = getTopFragmentMatches(protein_sequence, n, masses.data(), qscores.data(),
                                      charges.data(), window_starts.data(), window_ends.data(),
                                      ion_types_temp.get(), fragment_indices.data(), stored_ms2, discard, "HCD");

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

  int FragmentAnalysis::getAmbiguityEnclosingIons(const String& protein_sequence,
                                                  int n,
                                                  double* masses,
                                                  double* qscores,
                                                  int* charges,
                                                  double* window_starts,
                                                  double* window_ends,
                                                  char* ion_types,
                                                  int* fragment_indices,
                                                  DeconvolvedSpectrum& stored_ms2,
                                                  ProteoformMatch& result,
                                                  const String& fragmentation_method,
                                                  double tolerance_ppm,
                                                  FragmentScores* frag_scores)
  {
    // Get fragment matches AND PTM sites from FLASHExtender
    std::vector<TagBasedFragmentMatch> fragment_ion_match;
    int match_count = runTagBasedFragmentMatching_(protein_sequence, fragment_ion_match, stored_ms2, result, fragmentation_method, tolerance_ppm);

    if (match_count == 0)
    {
      std::cout << "[getAmbiguityEnclosingIons] No fragment matches found" << std::endl;
      return 0;
    }
    std::cout << "[getAmbiguityEnclosingIons] Found " << match_count << " fragment matches" << std::endl;

    if (result.ptm_sites.empty())
    {
      std::cout << "[getAmbiguityEnclosingIons] No PTM sites detected by FLASHExtender" << std::endl;
      return 0;
    }

    // Debug output for PTM sites
    std::cout << "[getAmbiguityEnclosingIons] FLASHExtender detected " << result.ptm_sites.size() << " PTM sites:" << std::endl;
    for (const auto& site : result.ptm_sites)
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
      float intensity;  // max-charge intensity (sort key, copied from TagBasedFragmentMatch.intensity); NOT a qscore
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
    for (const auto& site : result.ptm_sites)
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
        brackets.left_ions.push_back({static_cast<float>(m->intensity), m->peak_index, m->ion_type,
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
        brackets.right_ions.push_back({static_cast<float>(m->intensity), m->peak_index, m->ion_type,
                                       m->fragment_index, 1, m->fragment_index});
      }

      // Priority = max intensity of primary brackets
      brackets.priority = 0.0f;
      if (!brackets.left_ions.empty()) brackets.priority = std::max(brackets.priority, brackets.left_ions[0].intensity);
      if (!brackets.right_ions.empty()) brackets.priority = std::max(brackets.priority, brackets.right_ions[0].intensity);

      // Debug output for this PTM site
      std::cout << "[getAmbiguityEnclosingIons] PTM [" << site.start_position << "-" << site.end_position << "]:" << std::endl;
      if (!brackets.left_ions.empty())
      {
        std::cout << "  Left bracket (primary): " << brackets.left_ions[0].ion_type << brackets.left_ions[0].fragment_index
                  << " (intensity:" << brackets.left_ions[0].intensity << ", total candidates: " << brackets.left_ions.size() << ")" << std::endl;
      }
      else
      {
        std::cout << "  Left bracket: none found" << std::endl;
      }
      if (!brackets.right_ions.empty())
      {
        std::cout << "  Right bracket (primary): " << brackets.right_ions[0].ion_type << brackets.right_ions[0].fragment_index
                  << " (intensity:" << brackets.right_ions[0].intensity << ", total candidates: " << brackets.right_ions.size() << ")" << std::endl;
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

    // Helper to output one ion (single most-intense charge)
    int output_idx = 0;
    auto output_ion = [&](const EnclosingIon& ion) {
      if (used_peaks.find(ion.peak_index) != used_peaks.end()) return;  // Already output
      used_peaks.insert(ion.peak_index);

      const auto& pg = stored_ms2[ion.peak_index];
      double mono_mass = pg.getMonoMass();

      int charge = pg.getMaxIntensityAbsCharge();
      masses[output_idx] = mono_mass;
      qscores[output_idx] = ion.intensity;
      charges[output_idx] = charge;
      ion_types[output_idx] = ion.ion_type;
      fragment_indices[output_idx] = ion.fragment_index;
      auto [mz1, mz2] = pg.getMzRange(charge);
      window_starts[output_idx] = mz1 - optimal_window_margin_;
      window_ends[output_idx] = mz2 + optimal_window_margin_;
      if (frag_scores) frag_scores[output_idx] = FragmentScores::fromPeakGroup(pg, charge);
      ++output_idx;
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

  int FragmentAnalysis::getAmbiguityEnclosingIonsPy(const String& protein_sequence,
                                                    int n,
                                                    std::vector<double>& masses,
                                                    std::vector<double>& qscores,
                                                    std::vector<int>& charges,
                                                    std::vector<double>& window_starts,
                                                    std::vector<double>& window_ends,
                                                    std::vector<int>& is_b_ions,
                                                    std::vector<int>& fragment_indices,
                                                    DeconvolvedSpectrum& stored_ms2)
  {
    masses.resize(n);
    qscores.resize(n);
    charges.resize(n);
    window_starts.resize(n);
    window_ends.resize(n);

    // Use raw char array for the C-style function
    std::unique_ptr<char[]> ion_types_temp(new char[n]);
    fragment_indices.resize(n);

    ProteoformMatch discard;
    int count = getAmbiguityEnclosingIons(protein_sequence, n, masses.data(), qscores.data(),
                                          charges.data(), window_starts.data(), window_ends.data(),
                                          ion_types_temp.get(), fragment_indices.data(), stored_ms2, discard, "HCD");

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

  int FragmentAnalysis::getTerminalFragmentIons(
      const String& protein_sequence,
      int n,
      double* masses,
      double* qscores,
      int* charges,
      double* window_starts,
      double* window_ends,
      char* ion_types,
      int* fragment_indices,
      DeconvolvedSpectrum& stored_ms2,
      ProteoformMatch& result,
      const String& fragmentation_method,
      double tolerance_ppm,
      FragmentScores* frag_scores)
  {
    // Run fragment matching to get all matches
    std::vector<TagBasedFragmentMatch> matches;
    runTagBasedFragmentMatching_(protein_sequence, matches, stored_ms2, result, fragmentation_method, tolerance_ppm);

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

    // Sort prefix ions by fragment_index descending (rightmost first), then intensity descending
    std::sort(prefix_ions.begin(), prefix_ions.end(), [](const auto& a, const auto& b) {
      if (a.fragment_index != b.fragment_index) return a.fragment_index > b.fragment_index;
      return a.intensity > b.intensity;
    });

    // Sort suffix ions by fragment_index descending (leftmost first), then intensity descending
    std::sort(suffix_ions.begin(), suffix_ions.end(), [](const auto& a, const auto& b) {
      if (a.fragment_index != b.fragment_index) return a.fragment_index > b.fragment_index;
      return a.intensity > b.intensity;
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
        const auto& pg = stored_ms2[selected->peak_index];

        // Single charge behavior
        masses[output_idx] = selected->observed_mass;
        qscores[output_idx] = selected->intensity;
        charges[output_idx] = selected->charge;
        ion_types[output_idx] = selected->ion_type;
        fragment_indices[output_idx] = selected->fragment_index;
        auto [mz_start, mz_end] = pg.getMzRange(selected->charge);
        window_starts[output_idx] = mz_start - optimal_window_margin_;
        window_ends[output_idx] = mz_end + optimal_window_margin_;
        if (frag_scores) frag_scores[output_idx] = FragmentScores::fromPeakGroup(pg, std::abs(selected->charge));
        ++output_idx;
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

  int FragmentAnalysis::getTerminalFragmentIonsPy(
      const String& protein_sequence,
      int n,
      std::vector<double>& masses,
      std::vector<double>& qscores,
      std::vector<int>& charges,
      std::vector<double>& window_starts,
      std::vector<double>& window_ends,
      std::vector<int>& is_b_ions,
      std::vector<int>& fragment_indices,
      DeconvolvedSpectrum& stored_ms2)
  {
    masses.resize(n);
    qscores.resize(n);
    charges.resize(n);
    window_starts.resize(n);
    window_ends.resize(n);

    // Use raw char array for the C-style function
    std::unique_ptr<char[]> ion_types_temp(new char[n]);
    fragment_indices.resize(n);

    ProteoformMatch discard;
    int count = getTerminalFragmentIons(
        protein_sequence, n,
        masses.data(), qscores.data(), charges.data(),
        window_starts.data(), window_ends.data(),
        ion_types_temp.get(), fragment_indices.data(), stored_ms2, discard, "HCD");

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

  double FragmentAnalysis::windowSnr(const MSSpectrum& source, double lo, double hi, double signal_intensity)
  {
    // I2: SNR over the ACTUAL commanded isolation window. The "signal" (the selected charge's intensity,
    // = PeakGroup::getChargeIntensity, already computed by the engine) is passed in; here we measure the
    // total raw intensity in [lo, hi] and treat the remainder as co-isolation noise. No deconvolution math.
    if (signal_intensity <= 0.0 || hi <= lo) return 0.0;
    double total = 0.0;
    for (const auto& peak : source)
    {
      const double mz = peak.getMZ();
      if (mz >= lo && mz <= hi) total += peak.getIntensity();
    }
    const double noise = std::max(total - signal_intensity, 0.0);
    // Floor noise at 0.1% of signal so a perfectly pure window yields a bounded high SNR (~1000) rather
    // than dividing by ~0. (D5: the ε in "signal / (noise + ε)", chosen relative so magnitudes stay sane.)
    const double denom = std::max(noise, 1e-3 * signal_intensity);
    return signal_intensity / denom;
  }

  std::string FragmentAnalysis::toProForma(const std::string& sequence,
                                           const std::vector<PTMSite>& ptm_sites)
  {
    if (ptm_sites.empty()) return sequence;

    std::vector<PTMSite> sorted_sites = ptm_sites;
    std::sort(sorted_sites.begin(), sorted_sites.end(),
              [](const PTMSite& a, const PTMSite& b) { return a.start_position > b.start_position; });

    std::string result = sequence;
    for (const auto& site : sorted_sites)
    {
      std::ostringstream mass_ss;
      mass_ss << std::fixed << std::setprecision(4);
      if (site.mass_shift >= 0)
        mass_ss << "[+" << site.mass_shift << "]";
      else
        mass_ss << "[" << site.mass_shift << "]";
      std::string mass_str = mass_ss.str();

      if (site.start_position == site.end_position)
      {
        // Localized: insert after residue at 1-based position
        int insert_pos = site.start_position;  // 1-based position -> insert at this index
        if (insert_pos >= 0 && insert_pos <= static_cast<int>(result.size()))
          result.insert(insert_pos, mass_str);
      }
      else
      {
        // Ambiguous: wrap region in parentheses
        int end_idx = site.end_position;  // insert ')' + mass after this 1-based position
        int start_idx = site.start_position - 1;  // insert '(' before this 1-based position (0-based)
        if (end_idx >= 0 && end_idx <= static_cast<int>(result.size()))
          result.insert(end_idx, ")" + mass_str);
        if (start_idx >= 0 && start_idx <= static_cast<int>(result.size()))
          result.insert(start_idx, "(");
      }
    }
    return result;
  }

} // namespace OpenMS
