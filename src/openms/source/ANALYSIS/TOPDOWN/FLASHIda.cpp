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


#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <functional>
#include <map>
#include <sstream>
#ifdef _OPENMP
  #include <omp.h>
#endif

namespace OpenMS
{

  DoubleList FLASHIda::buildToleranceList_(const Config& config)
  {
    DoubleList tol;
    for (const auto& [lvl, cfg] : config.levels())
      tol.push_back(cfg.tolerance_ppm);
    return tol;
  }

/// constructor
FLASHIda::FLASHIda(char* arg) :
    config_(std::string(arg)),
    queue_(config_),
    deconv_(config_, buildToleranceList_(config_)),
    fragments_(config_),
    selection_(config_, deconv_),
    quant_(config_),
    faims_(config_),
    exploration_(config_, fragments_)
{
  #ifdef _OPENMP
    omp_set_num_threads(4);
  #endif

    engine_start_time_ = std::chrono::steady_clock::now();

    // Initialize FAIMS CV atomic for getNextScanCommand reads.
    // relaxed is safe: no other thread can observe this object before construction completes.
    current_faims_cv_.store(faims_.isEnabled() ? faims_.currentCV() : 0.0, std::memory_order_relaxed);

    const auto& rt_cfg = config_.runtime();
    if (!rt_cfg.ida_log_path.empty())
    {
      ida_log_stream_.open(rt_cfg.ida_log_path, std::ios::app);
    }
    if (!rt_cfg.scan_commands_path.empty())
    {
      commands_tsv_stream_.open(rt_cfg.scan_commands_path, std::ios::app);
      if (commands_tsv_stream_.is_open())
      {
        commands_tsv_stream_ << "tracking_id\tms_level\tscan_type\tenqueue_ts\tpriority\t"
                             << "faims_cv\tmono_mass\tcharge\tprecursor_mz\tisolation_width\t"
                             << "collision_energy\tactivation\tqscore\tcharge_cos\tcharge_snr\t"
                             << "iso_cos\tsnr\tcharge_score\tppm_error\tprecursor_intensity\t"
                             << "peakgroup_intensity\thcd_energy\tparent_tracking_id\t"
                             << "ion_type\tion_index\t"
                             << "reaction_time\treagent_max_it\treagent_agc_target\n";
        commands_tsv_stream_.flush();
      }
    }
    if (!rt_cfg.scan_results_path.empty())
    {
      results_tsv_stream_.open(rt_cfg.scan_results_path, std::ios::app);
      if (results_tsv_stream_.is_open())
      {
        results_tsv_stream_ << "tracking_id\tresolve_ts\tduration_ms\treceived_ts\tduration_received_ms\trt\t"
                            << "mass_count\tcommands_pushed\tchild_ids\t"
                            << "tag_count\tmatched_protein\tproteoform_sequence\t"
                            << "tic_coverage\tfragment_count\t"
                            << "exploration_group_id\texploration_metric\t"
                            << "variant_index\ttotal_variants\t"
                            << "collision_energy\texploration_score\tremaining_ratio\t"
                            << "activation_type\treaction_time\t"
                            << "deconv_masses\tdeconv_intensities\tdeconv_min_charge\tdeconv_max_charge\tparent_tracking_id\t"
                            << "dequeue_ts\tqueue_duration_ms\tinstrument_duration_ms\tprocessing_duration_ms\n";
        results_tsv_stream_.flush();
      }
    }
    if (!rt_cfg.identification_path.empty())
    {
      identification_tsv_stream_.open(rt_cfg.identification_path, std::ios::app);
      if (identification_tsv_stream_.is_open())
      {
        identification_tsv_stream_ << "tracking_id\tproteoform\tstart_pos\tend_pos\t"
                                   << "ppm_offset\tcorrection_factor\t"
                                   << "ms1_precursor_mass\tms1_precursor_mz\tms1_precursor_charge\t"
                                   << "ms2_precursor_ion\tms2_precursor_mass\tms2_precursor_mz\tms2_precursor_charge\t"
                                   << "ms2_fragments\tms2_fragment_masses\t"
                                   << "ms3_fragments\tms3_fragment_masses\n";
        identification_tsv_stream_.flush();
      }
    }
  }

  FLASHIda::~FLASHIda()
  {
    if (ida_log_stream_.is_open()) ida_log_stream_.close();
    if (commands_tsv_stream_.is_open()) commands_tsv_stream_.close();
    if (results_tsv_stream_.is_open()) results_tsv_stream_.close();
    if (identification_tsv_stream_.is_open()) identification_tsv_stream_.close();
  }

  // Fragment analysis: C-pointer overloads delegate to fragments_ + deconv_.storedMS2()

  int FLASHIda::getBestMS2Masses(int n, double* masses, double* qscores, int* charges,
                                 double* window_starts, double* window_ends)
  {
    if (!deconv_.hasStoredMS2()) return 0;
    return fragments_.getBestMS2Masses(n, masses, qscores, charges, window_starts, window_ends,
                                       deconv_.storedMS2());
  }

  int FLASHIda::getTopFragmentMatches(const String& protein_sequence, int n,
                                      double* masses, double* qscores, int* charges,
                                      double* window_starts, double* window_ends,
                                      char* ion_types, int* fragment_indices,
                                      const String& fragmentation_method)
  {
    if (!deconv_.hasStoredMS2()) return 0;
    return fragments_.getTopFragmentMatches(protein_sequence, n, masses, qscores, charges,
                                            window_starts, window_ends, ion_types, fragment_indices,
                                            deconv_.storedMS2(), fragmentation_method);
  }

  int FLASHIda::getTerminalFragmentIons(const String& protein_sequence, int n,
                                        double* masses, double* qscores, int* charges,
                                        double* window_starts, double* window_ends,
                                        char* ion_types, int* fragment_indices,
                                        const String& fragmentation_method)
  {
    if (!deconv_.hasStoredMS2()) return 0;
    return fragments_.getTerminalFragmentIons(protein_sequence, n, masses, qscores, charges,
                                              window_starts, window_ends, ion_types, fragment_indices,
                                              deconv_.storedMS2(), fragmentation_method);
  }

  int FLASHIda::getAmbiguityEnclosingIons(const String& protein_sequence, int n,
                                          double* masses, double* qscores, int* charges,
                                          double* window_starts, double* window_ends,
                                          char* ion_types, int* fragment_indices,
                                          const String& fragmentation_method)
  {
    if (!deconv_.hasStoredMS2()) return 0;
    return fragments_.getAmbiguityEnclosingIons(protein_sequence, n, masses, qscores, charges,
                                                window_starts, window_ends, ion_types, fragment_indices,
                                                deconv_.storedMS2(), fragmentation_method);
  }

  int FLASHIda::getConfigInt(const std::string& key) const
  {
    return config_.getInt(key);
  }

  double FLASHIda::getConfigDouble(const std::string& key) const
  {
    return config_.getDouble(key);
  }

  // ---- Logging writer implementations ----

  void FLASHIda::writeIDALogEntry_(double rt,
                                    const std::string& tracking_id,
                                    const std::vector<ScanCommand>& ms2_commands,
                                    const DeconvolvedSpectrum& all_peak_groups)
  {
    if (!ida_log_stream_.is_open()) return;

    // MS1 header line
    ida_log_stream_ << "MS1 Scan# 0"
                    << " RT " << std::fixed << std::setprecision(4) << rt
                    << " (Access ID " << tracking_id << ") - "
                    << ms2_commands.size() << " targets\n";

    for (const auto& cmd : ms2_commands)
    {
      double w1 = 0, w2 = 0;
      int charge = 0;
      if (cmd.num_stages > 0)
      {
        double center = cmd.stages[0].precursor_mz;
        double half_width = cmd.stages[0].isolation_width / 2.0;
        w1 = center - half_width;
        w2 = center + half_width;
        charge = cmd.stages[0].charge_state;
      }

      ida_log_stream_ << "Mass=" << std::defaultfloat << cmd.mono_mass
                      << "\tZ=" << charge
                      << "\tScore=" << std::fixed << std::setprecision(5) << cmd.qscore
                      << "\tWindow=[" << std::setprecision(4) << w1 << "-" << w2 << "]"
                      << "\tPrecursorIntensity=" << std::setprecision(5) << cmd.precursor_intensity
                      << "\tPrecursorMassIntensity=" << std::setprecision(5) << cmd.peakgroup_intensity
                      << "\tFeatures=["
                        << std::setprecision(6) << cmd.charge_cos << ","
                        << cmd.charge_snr << ","
                        << cmd.iso_cos << ","
                        << cmd.snr << ","
                        << cmd.charge_score << ","
                        << cmd.ppm_error << "]"
                      << "\tChargeRange=[" << charge << "-" << charge << "]"
                      << "\tHCD=" << cmd.hcd_energy << "\n";
    }

    // AllMass line
    ida_log_stream_ << "AllMass=";
    for (size_t i = 0; i < all_peak_groups.size(); i++)
    {
      if (i > 0) ida_log_stream_ << " ";
      ida_log_stream_.precision(4);
      ida_log_stream_ << std::fixed << all_peak_groups[i].getMonoMass();
    }
    ida_log_stream_ << "\n";
    ida_log_stream_.flush();
  }

  std::string FLASHIda::scanTypeFromDescription_(const ScanCommand& cmd)
  {
    if (std::strlen(cmd.scan_description) < 4)
      return "unknown";
    switch (cmd.scan_description[3])
    {
      case 'S': return "survey";
      case 'A': return "agc";
      case 'R': return "recording";
      case 'F': return "followup";
      case 'C': return "conditional";
      case 'E': return "exploration";
      default: return "unknown";
    }
  }

  void FLASHIda::writeScanCommandRow_(const ScanCommand& cmd)
  {
    std::lock_guard<std::mutex> lock(analysis_mutex_);
    if (!commands_tsv_stream_.is_open()) return;

    std::string id_str = ScanCommandQueue::encode(cmd.scan_id);
    std::string scan_type = scanTypeFromDescription_(cmd);

    std::string charges, precursor_mzs, iso_widths, col_energies, activations;
    std::string reaction_times, reagent_max_its, reagent_agc_targets;
    for (int i = 0; i < cmd.num_stages; ++i)
    {
      if (i > 0) { charges += ";"; precursor_mzs += ";"; iso_widths += ";"; col_energies += ";"; activations += ";"; reaction_times += ";"; reagent_max_its += ";"; reagent_agc_targets += ";"; }
      charges += std::to_string(cmd.stages[i].charge_state);
      std::ostringstream mz_os, iw_os, ce_os;
      mz_os << cmd.stages[i].precursor_mz;
      iw_os << cmd.stages[i].isolation_width;
      ce_os << cmd.stages[i].collision_energy;
      precursor_mzs += mz_os.str();
      iso_widths += iw_os.str();
      col_energies += ce_os.str();
      activations += cmd.stages[i].activation_type;
      std::ostringstream rt_os, rmi_os;
      rt_os << cmd.stages[i].reaction_time;
      rmi_os << cmd.stages[i].reagent_max_it;
      reaction_times += rt_os.str();
      reagent_max_its += rmi_os.str();
      reagent_agc_targets += std::to_string(cmd.stages[i].reagent_agc_target);
    }

    std::string parent_id(cmd.parent_scan_id);
    std::string ion_type;
    int ion_index = 0;
    if (cmd.msn_level == 3 && std::strlen(cmd.scan_description) > 4)
    {
      std::string desc(cmd.scan_description);
      auto at_pos = desc.find('@');
      if (at_pos != std::string::npos)
      {
        size_t pos = at_pos + 1;
        while (pos < desc.size() && (std::isdigit(desc[pos]) || desc[pos] == '-'))
          pos++;
        if (pos < desc.size() && std::isalpha(desc[pos]))
        {
          ion_type = std::string(1, desc[pos]);
          pos++;
          if (pos < desc.size())
            ion_index = std::atoi(desc.c_str() + pos);
        }
      }
    }

    commands_tsv_stream_ << id_str << "\t"
                         << cmd.msn_level << "\t"
                         << scan_type << "\t"
                         << cmd.enqueue_timestamp_ms << "\t"
                         << cmd.priority << "\t"
                         << cmd.faims_cv << "\t"
                         << cmd.mono_mass << "\t"
                         << charges << "\t"
                         << precursor_mzs << "\t"
                         << iso_widths << "\t"
                         << col_energies << "\t"
                         << activations << "\t"
                         << cmd.qscore << "\t"
                         << cmd.charge_cos << "\t"
                         << cmd.charge_snr << "\t"
                         << cmd.iso_cos << "\t"
                         << cmd.snr << "\t"
                         << cmd.charge_score << "\t"
                         << cmd.ppm_error << "\t"
                         << cmd.precursor_intensity << "\t"
                         << cmd.peakgroup_intensity << "\t"
                         << cmd.hcd_energy << "\t"
                         << parent_id << "\t"
                         << ion_type << "\t"
                         << ion_index << "\t"
                         << reaction_times << "\t"
                         << reagent_max_its << "\t"
                         << reagent_agc_targets << "\n";
    commands_tsv_stream_.flush();
  }

  void FLASHIda::writeScanResultRow_(const std::string& tracking_id, double rt,
                                      int mass_count, int commands_pushed,
                                      const std::vector<std::string>& child_ids,
                                      int tag_count, const std::string& matched_protein,
                                      const std::string& proteoform_sequence,
                                      uint64_t enqueue_ts, uint64_t dequeue_ts, uint64_t received_ts,
                                      const DeconvolvedSpectrum* deconv_spectrum,
                                      const std::string& parent_tracking_id,
                                      float tic_coverage, int fragment_count,
                                      int exploration_group_id, int exploration_metric,
                                      int variant_index, int total_variants,
                                      double collision_energy, double exploration_score,
                                      double remaining_ratio,
                                      const std::string& activation_type,
                                      double reaction_time)
  {
    if (!results_tsv_stream_.is_open()) return;

    auto now = std::chrono::steady_clock::now();
    uint64_t resolve_ts = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()).count();
    uint64_t duration = (enqueue_ts > 0) ? (resolve_ts - enqueue_ts) : 0;
    uint64_t duration_received = (enqueue_ts > 0 && received_ts > 0) ? (received_ts - enqueue_ts) : 0;
    uint64_t queue_duration = (dequeue_ts > 0 && enqueue_ts > 0) ? (dequeue_ts - enqueue_ts) : 0;
    uint64_t instrument_duration = (received_ts > 0 && dequeue_ts > 0) ? (received_ts - dequeue_ts) : 0;
    uint64_t processing_duration = (received_ts > 0) ? (resolve_ts - received_ts) : 0;

    std::string child_str;
    for (size_t i = 0; i < child_ids.size(); i++)
    {
      if (i > 0) child_str += ";";
      child_str += child_ids[i];
    }

    results_tsv_stream_ << tracking_id << "\t"
                        << resolve_ts << "\t"
                        << duration << "\t"
                        << received_ts << "\t"
                        << duration_received << "\t"
                        << rt << "\t"
                        << mass_count << "\t"
                        << commands_pushed << "\t"
                        << child_str << "\t"
                        << tag_count << "\t"
                        << matched_protein << "\t"
                        << proteoform_sequence << "\t"
                        << tic_coverage << "\t"
                        << fragment_count << "\t"
                        << exploration_group_id << "\t"
                        << exploration_metric << "\t"
                        << variant_index << "\t"
                        << total_variants << "\t"
                        << collision_energy << "\t"
                        << exploration_score << "\t"
                        << remaining_ratio << "\t";
    results_tsv_stream_ << activation_type << "\t"
                        << reaction_time << "\t";

    // Deconvolved masses and intensities (semicolon-delimited)
    if (deconv_spectrum != nullptr && deconv_spectrum->size() > 0)
    {
      for (size_t i = 0; i < deconv_spectrum->size(); i++)
      {
        if (i > 0) results_tsv_stream_ << ";";
        results_tsv_stream_ << (*deconv_spectrum)[i].getMonoMass();
      }
      results_tsv_stream_ << "\t";
      for (size_t i = 0; i < deconv_spectrum->size(); i++)
      {
        if (i > 0) results_tsv_stream_ << ";";
        results_tsv_stream_ << (*deconv_spectrum)[i].getIntensity();
      }
      results_tsv_stream_ << "\t";
      for (size_t i = 0; i < deconv_spectrum->size(); i++)
      {
        if (i > 0) results_tsv_stream_ << ";";
        results_tsv_stream_ << std::get<0>((*deconv_spectrum)[i].getAbsChargeRange());
      }
      results_tsv_stream_ << "\t";
      for (size_t i = 0; i < deconv_spectrum->size(); i++)
      {
        if (i > 0) results_tsv_stream_ << ";";
        results_tsv_stream_ << std::get<1>((*deconv_spectrum)[i].getAbsChargeRange());

      }
    }
    else
    {
      results_tsv_stream_ << "\t";
    }
    results_tsv_stream_ << "\t" << parent_tracking_id
                        << "\t" << dequeue_ts
                        << "\t" << queue_duration
                        << "\t" << instrument_duration
                        << "\t" << processing_duration << "\n";
    results_tsv_stream_.flush();
  }

  void FLASHIda::writeIdentificationRow_(
    const std::string& tracking_id,
    const Exploration::MS2Context& ctx,
    const MS3FragmentMatcher::MatchResult& result)
  {
    // No lock — called from processScan() which already holds analysis_mutex_
    if (!identification_tsv_stream_.is_open()) return;
    if (result.matches.empty()) return;

    std::ostringstream ms3_frags, ms3_masses, ms2_frags, ms2_masses;
    ms3_frags << std::fixed << std::setprecision(4);
    ms3_masses << std::fixed << std::setprecision(4);
    ms2_frags << std::fixed << std::setprecision(4);
    ms2_masses << std::fixed << std::setprecision(4);

    for (size_t i = 0; i < result.matches.size(); ++i)
    {
      const auto& fm = result.matches[i];
      if (i > 0)
      {
        ms3_frags << ";";
        ms3_masses << ";";
        ms2_frags << ";";
        ms2_masses << ";";
      }
      ms3_frags << fm.ms3_ion_type << fm.ms3_ion_index;
      ms3_masses << fm.observed_mass;
      ms2_frags << fm.ms2_equiv_type << fm.ms2_equiv_index;
      ms2_masses << fm.adjusted_mass;
    }

    std::string precursor_ion;
    if (ctx.fragment_ion_type != '\0')
      precursor_ion = std::string(1, ctx.fragment_ion_type) + std::to_string(ctx.fragment_ion_index);

    identification_tsv_stream_
      << tracking_id << "\t"
      << ctx.proteoform_sequence << "\t"
      << ctx.start_pos << "\t"
      << ctx.end_pos << "\t"
      << std::fixed << std::setprecision(2) << result.ppm_offset << "\t"
      << std::setprecision(8) << result.correction_factor << "\t"
      << std::setprecision(4) << ctx.ms1_precursor_mass << "\t"
      << ctx.ms1_precursor_mz << "\t"
      << ctx.ms1_precursor_charge << "\t"
      << precursor_ion << "\t"
      << ctx.fragment_mass << "\t"
      << ctx.fragment_mz << "\t"
      << ctx.fragment_charge << "\t"
      << ms2_frags.str() << "\t"
      << ms2_masses.str() << "\t"
      << ms3_frags.str() << "\t"
      << ms3_masses.str() << "\n";
    identification_tsv_stream_.flush();
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

  // ---- FAIMS CV adaptive skip state machine ----
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
    std::lock_guard<std::mutex> lock(analysis_mutex_);

    // Centralized tracking ID extraction — early return for instrument-triggered scans
    std::string desc_str = scan_description ? std::string(scan_description) : "";
    if (desc_str.size() < 3) return 0;
    std::string id_str = desc_str.substr(0, 3);
    int tracking_id = queue_.decode(id_str);

    // Retrieve timestamps from pending map (enqueue set by push(), dequeue set by dequeue())
    uint64_t enqueue_ts = 0;
    uint64_t dequeue_ts = 0;
    {
      auto peeked = queue_.peekPending(tracking_id);
      if (peeked.has_value())
      {
        enqueue_ts = peeked->enqueue_timestamp_ms;
        dequeue_ts = peeked->dequeue_timestamp_ms;
      }
    }

    // Stamp received timestamp (instrument → C++ handoff)
    uint64_t received_ts = static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now().time_since_epoch()).count());

    if (ms_level == 1)
    {
      // Selection=none: skip MS1 precursor selection entirely
      if (config_.level(1).selection == SelectionMetric::None)
        return 0;

      // MS1 path: deconvolve, score, filter, select top-N, push MS2 commands
      int n = selection_.filterAndRank(mzs, ints, length, rt_min, 1, faims_cv);
      const auto& selected = selection_.selectedPeakGroups();
      const auto& sel_charges = selection_.triggerCharges();
      const auto& sel_hcds = selection_.triggerHcds();
      int commands_pushed = 0;
      std::vector<ScanCommand> ms2_commands;
      
      if (config_.hasExploration(2))
      {
        // Exploration path: initiate CE sweep variants INSTEAD of regular MS2
        for (int i = 0; i < n; i++)
        {
          ScanCommand ms1_ctx{};
          ms1_ctx.scan_id = tracking_id;
          auto cmds = exploration_.initiate(2, selected[i], sel_charges[i], faims_cv, queue_, &ms1_ctx);
          for (auto& c : cmds)
          {
            std::string ms1_enc = ScanCommandQueue::encode(tracking_id);
            std::strncpy(c.parent_scan_id, ms1_enc.c_str(), 3);
            c.parent_scan_id[3] = '\0';
            queue_.push(c);
            ms2_commands.push_back(c);
            commands_pushed++;
          }
        }
      }
      else
      {
        // Normal path: push MS2 for each precursor, for each scan config
        for (int i = 0; i < n; i++)
        {
          for (const auto& sc : config_.level(2).scans)
          {
            ScanConfig ms2_config = sc;
            if (config_.targeting().use_idscore)
              ms2_config.collision_energy = sel_hcds[i];
            ScanCommand cmd = queue_.buildMS2(selected[i], sel_charges[i], ms2_config, 2, tracking_id);
            cmd.faims_cv = faims_cv;  // MS2 carries parent MS1's CV
            queue_.push(cmd);
            ms2_commands.push_back(cmd);
            commands_pushed++;
          }
        }
      }

      // IDA log entry
      writeIDALogEntry_(rt_min, id_str, ms2_commands, selection_.deconvolvedMS1());

      // Results TSV entry for MS1
      std::vector<std::string> child_ids;
      for (const auto& c : ms2_commands)
        child_ids.push_back(ScanCommandQueue::encode(c.scan_id));
      int all_mass_count = static_cast<int>(selection_.deconvolvedMS1().size());
      writeScanResultRow_(id_str, rt_min, all_mass_count, commands_pushed,
                          child_ids, 0, "", "", enqueue_ts, dequeue_ts, received_ts,
                          &selection_.deconvolvedMS1(), "");

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

        std::string next_ms1_id = ScanCommandQueue::encode(ms1.scan_id);
        std::snprintf(ms1.scan_description, 16, "%sS", next_ms1_id.c_str());

        queue_.push(ms1);
        std::cout << "[TRACK-CREATE] id=" << next_ms1_id << " ms_level=1 type=cv_transition cv=" << next_cv << std::endl;
      }

      // Update atomics for lock-free reads by getNextScanCommand
      exploration_active_.store(exploration_.activeGroupCount() > 0, std::memory_order_release);
      current_faims_cv_.store(faims_.isEnabled() ? faims_.currentCV() : 0.0, std::memory_order_release);

      return commands_pushed;
    }
    else if (ms_level == 2)
    {
      int commands_pushed = 0;

      // Check if this is an exploration variant (before pending scan lookup)
      if (exploration_.isExplorationVariant(tracking_id))
      {
        auto info = exploration_.feedResult(tracking_id, mzs, ints, length, rt_min, queue_);
        for (auto& c : info.commands)
        {
          std::strncpy(c.parent_scan_id, info.parent_scan_id, 4);
          queue_.push(c);
        }
        
        const bool has_expl_ms2 = exploration_.exploration_deconv_ != nullptr &&
                exploration_.exploration_deconv_->hasStoredMS2();
        int expl_mass_count = has_expl_ms2 ? static_cast<int>(exploration_.exploration_deconv_->storedMS2().size()) : 0;
        const DeconvolvedSpectrum* ms2_spec = has_expl_ms2 ? &exploration_.exploration_deconv_->storedMS2() : nullptr;
        std::string parent_id(info.parent_scan_id);
        writeScanResultRow_(id_str, rt_min, expl_mass_count, static_cast<int>(info.commands.size()),
                            {}, 0, info.matched_protein, info.proteoform_sequence, enqueue_ts, dequeue_ts, received_ts,
                            ms2_spec, parent_id,
                            info.tic_coverage, info.fragment_count,
                            info.group_id, info.exploration_metric,
                            info.variant_index, info.total_variants,
                            info.collision_energy, info.score, info.remaining_ratio,
                            info.activation_type, info.reaction_time);

        exploration_active_.store(exploration_.activeGroupCount() > 0, std::memory_order_release);
        return static_cast<int>(info.commands.size());
      }

      auto resolved = queue_.resolvePending(tracking_id);
      if (!resolved.has_value())
      {
        std::cout << "[TRACK-RESOLVE] id=" << id_str << " status=not_found" << std::endl;
        return 0;
      }
      ScanCommand ctx = resolved.value();
      std::vector<std::string> child_ids;

      // Step 3: Deconvolve MS2 with precursor context
      double precursor_mass = 0;
      int precursor_charge = 0;
      if (ctx.num_stages > 0)
      {
        precursor_charge = ctx.stages[0].charge_state;
        precursor_mass = ctx.mono_mass;
      }

      deconv_.deconvolveMSn(mzs, ints, length, rt_min, precursor_mass, precursor_charge);

      // Step 4: Route by mode
      // Tag-based targeting
      bool tags_found = false;
      if (selection_.hasTargetProteinDatabase() && precursor_mass > 0)
      {
        tags_found = selection_.processMS2ForTagBasedTargeting(precursor_mass);
      }

      // Quantification follow-up (independent of tags)
      if (config_.quantification().enabled && !config_.quantification().follow_up_scan.activation.empty())
      {
        if (quant_.isDifferentiallyAbundant(mzs, ints, length, rt_min, 2, "ms2_quant",
                                            config_.quantification().reporter_mz_tol, config_.quantification().fold_change_threshold, false))
        {
          auto followup = queue_.buildFollowUp(ctx, config_.quantification().follow_up_scan, 'F');
          queue_.push(followup);
          child_ids.push_back(ScanCommandQueue::encode(followup.scan_id));
          commands_pushed++;
        }
      }

      // Conditional MS2 follow-up -- only when tags detected AND explicitly enabled
      if (config_.targeting().conditional_ms2_enabled && tags_found)
      {
        auto cond = queue_.buildFollowUp(ctx, config_.targeting().tagging_follow_up_scan, 'C');
        queue_.push(cond);
        child_ids.push_back(ScanCommandQueue::encode(cond.scan_id));
        commands_pushed++;
      }

      // Step 5: MS3 targeting via selection_strategy
      Exploration::NextLevelResult nlr;
      if (config_.level(2).selection != SelectionMetric::None)
      {
        std::cout << "flashida call site" << std::endl;
        nlr = exploration_.initiateNextLevel(2, deconv_.storedMS2(), ctx.faims_cv, queue_, &ctx);
        for (auto& c : nlr.commands)
        {
          std::string ms2_enc = ScanCommandQueue::encode(ctx.scan_id);
          std::strncpy(c.parent_scan_id, ms2_enc.c_str(), 3);
          c.parent_scan_id[3] = '\0';
          queue_.push(c);
          child_ids.push_back(ScanCommandQueue::encode(c.scan_id));
          commands_pushed++;
        }

        // Cache MS2 context for each MS3 command (for non-exploration identification.tsv)
        for (size_t ci = 0; ci < nlr.commands.size(); ++ci)
        {
          if (nlr.commands[ci].msn_level >= 3 && ci < nlr.ms3_contexts.size())
          {
            ms2_context_cache_[nlr.commands[ci].scan_id] = nlr.ms3_contexts[ci];
          }
        }
      }

      int ms2_mass_count = deconv_.hasStoredMS2() ? static_cast<int>(deconv_.storedMS2().size()) : 0;
      int tag_count = tags_found ? 1 : 0;
      const DeconvolvedSpectrum* ms2_spec = deconv_.hasStoredMS2() ? &deconv_.storedMS2() : nullptr;
      std::string parent_id(ctx.parent_scan_id);

      writeScanResultRow_(id_str, rt_min, ms2_mass_count, commands_pushed,
                          child_ids, tag_count, nlr.matched_protein, nlr.proteoform_sequence,
                          enqueue_ts, dequeue_ts, received_ts,
                          ms2_spec, parent_id,
                          nlr.tic_coverage, nlr.fragment_count);

      std::cout << "[TRACK-RESOLVE] id=" << id_str
                << " rt=" << rt_min
                << " commands=" << commands_pushed
                << std::endl;

      // Update atomic for lock-free reads by getNextScanCommand
      exploration_active_.store(exploration_.activeGroupCount() > 0, std::memory_order_release);

      return commands_pushed;
    }
    // MS3 (or higher)
    {
      // Route exploration variants to feedResult for scoring/winner selection
      if (exploration_.isExplorationVariant(tracking_id))
      {
        auto info = exploration_.feedResult(tracking_id, mzs, ints, length, rt_min, queue_);
        for (auto& c : info.commands)
        {
          std::strncpy(c.parent_scan_id, info.parent_scan_id, 4);
          queue_.push(c);
        }
        
        const bool has_expl_ms3 = exploration_.exploration_deconv_ != nullptr &&
                exploration_.exploration_deconv_->hasStoredMS2();
        int expl_mass_count = has_expl_ms3 ? static_cast<int>(exploration_.exploration_deconv_->storedMS2().size()) : 0;
        const DeconvolvedSpectrum* ms3_spec = has_expl_ms3 ? &exploration_.exploration_deconv_->storedMS2() : nullptr;
        std::string parent_id(info.parent_scan_id);
        writeScanResultRow_(id_str, rt_min, expl_mass_count,
                            static_cast<int>(info.commands.size()),
                            {}, 0, info.matched_protein, info.proteoform_sequence, enqueue_ts, dequeue_ts, received_ts,
                            ms3_spec, parent_id,
                            info.tic_coverage, info.fragment_count,
                            info.group_id, info.exploration_metric,
                            info.variant_index, info.total_variants,
                            info.collision_energy, info.score, info.remaining_ratio,
                            info.activation_type, info.reaction_time);

        if (!info.identification_result.matches.empty())
        {
          writeIdentificationRow_(id_str, info.ms2_context, info.identification_result);
        }

        exploration_active_.store(exploration_.activeGroupCount() > 0,
                                  std::memory_order_release);
        return static_cast<int>(info.commands.size());
      }

      // Non-exploration MS3: resolve pending, deconvolve with precursor context
      auto resolved = queue_.resolvePending(tracking_id);

      double precursor_mass = 0.0;
      int precursor_charge = 0;
      if (resolved.has_value() && resolved->num_stages >= 2)
      {
        precursor_charge = resolved->stages[1].charge_state;
        precursor_mass = resolved->mono_mass;
      }

      int ms3_mass_count = 0;
      if (mzs != nullptr && ints != nullptr && length > 0)
      {
        deconv_.deconvolveMSn(mzs, ints, length, rt_min, precursor_mass, precursor_charge);
        ms3_mass_count = deconv_.hasStoredMS2()
            ? static_cast<int>(deconv_.storedMS2().size()) : 0;
      }

      const DeconvolvedSpectrum* ms3_spec = deconv_.hasStoredMS2() ? &deconv_.storedMS2() : nullptr;

      // Identification: look up MS2 context from cache and run fragment matching
      {
        auto cache_it = ms2_context_cache_.find(tracking_id);
        if (cache_it != ms2_context_cache_.end() && ms3_spec != nullptr && !ms3_spec->empty())
        {
          const auto& mc = cache_it->second;
          MS3FragmentMatcher::ProteoformContext proto_ctx;
          proto_ctx.region_start = mc.start_pos;
          proto_ctx.region_end = mc.end_pos;
          proto_ctx.ptm_sites = mc.ptm_sites;

          std::vector<const DeconvolvedSpectrum*> spectra = {ms3_spec};
          std::vector<MS3FragmentMatcher::MatchResult> detailed;
          MS3FragmentMatcher::calibrateAndScore(
            spectra,
            config_.targeting().protein_sequence,
            proto_ctx,
            mc.fragment_ion_type,
            mc.fragment_ion_index,
            MS3FragmentMatcher::LOOSE_TOLERANCE_PPM,
            config_.level(3).tolerance_ppm,
            &detailed);

          if (!detailed.empty() && !detailed[0].matches.empty())
          {
            writeIdentificationRow_(id_str, mc, detailed[0]);
          }

          ms2_context_cache_.erase(cache_it);
        }
      }

      std::string parent_id;
      if (resolved.has_value())
        parent_id = std::string(resolved->parent_scan_id);
      writeScanResultRow_(id_str, rt_min, ms3_mass_count, 0,
                          {}, 0, "", "", enqueue_ts, dequeue_ts, received_ts,
                          ms3_spec, parent_id);
      return 0;
    }
  }

  int FLASHIda::getNextScanCommand(ScanCommand& out)
  {
    // No analysis_mutex_ acquired — queue methods lock internally, exploration/FAIMS via atomics
    double faims_cv = faims_.isEnabled() ? current_faims_cv_.load(std::memory_order_acquire) : 0.0;

    // Step 1: AGC scan if needed
    if (queue_.needsAGC())
    {
      out = queue_.makeAGC();
      out.faims_cv = faims_cv;
      out.scan_id = queue_.nextTrackingId();
      queue_.recordAGCTime();

      uint64_t now_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::steady_clock::now().time_since_epoch()).count());
      out.enqueue_timestamp_ms = now_ms;
      out.dequeue_timestamp_ms = now_ms;

      // Scan description: {3-char ID}A
      std::string id_str = ScanCommandQueue::encode(out.scan_id);
      std::snprintf(out.scan_description, 16, "%sA", id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=agc" << std::endl;
      queue_.registerPending(out);
      writeScanCommandRow_(out);
      return 1;
    }

    // Step 2: Cycle time -- force MS1 if too long since last survey scan
    // Suppressed while any exploration group is active.
    // Queued at priority 0 (not returned immediately) so it goes through normal dequeue.
    if (config_.scheduling().cycle_time_enabled 
        && queue_.msSinceLastMS1() > static_cast<uint64_t>(config_.scheduling().cycle_time_ms)
    )
    {
      ScanCommand ms1_cmd = queue_.makeMS1();
      ms1_cmd.faims_cv = faims_cv;
      ms1_cmd.scan_id = queue_.nextTrackingId();
      ms1_cmd.priority = 0;
      ms1_cmd.scan_description[3] = 'C';
      queue_.recordMS1Time();

      std::string id_str = ScanCommandQueue::encode(ms1_cmd.scan_id);
      std::snprintf(ms1_cmd.scan_description, 16, "%sS", id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << id_str << " ms_level=1 type=cycle_time" << std::endl;
      queue_.push(ms1_cmd);
      // Fall through to Step 3 (cleanup) and Step 4 (dequeue)
    }

    // Step 3: Cleanup expired commands
    queue_.cleanupExpired();

    // Step 4: Dequeue by priority (0 = highest -> 3 = lowest)
    auto dequeued = queue_.dequeue();
    if (dequeued.has_value())
    {
      out = dequeued.value();
      if (out.msn_level == 1 && out.is_agc == 0)
        queue_.recordMS1Time();
      // faims_cv already set at creation time (MS2 -> parent CV, CV-transition MS1 -> next CV)
      writeScanCommandRow_(out);
      return 1;
    }

    // Step 5: Idle cycle -- queue empty, keep the instrument busy with AGC + MS1
    // Create an AGC command (returned immediately) and push an MS1 at priority 3
    // into the queue as a fallback scan (lowest priority, behind follow-ups/MS3/MS2).
    {
      // 5a: AGC
      ScanCommand agc_cmd = queue_.makeAGC();
      agc_cmd.faims_cv = faims_cv;
      agc_cmd.scan_id = queue_.nextTrackingId();
      queue_.recordAGCTime();

      uint64_t now_ms = static_cast<uint64_t>(
        std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::steady_clock::now().time_since_epoch()).count());
      agc_cmd.enqueue_timestamp_ms = now_ms;
      agc_cmd.dequeue_timestamp_ms = now_ms;

      std::string agc_id_str = ScanCommandQueue::encode(agc_cmd.scan_id);
      std::snprintf(agc_cmd.scan_description, 16, "%sA", agc_id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << agc_id_str << " ms_level=1 type=idle_agc" << std::endl;

      // 5b: MS1 -- use default priority 3 (lowest, behind follow-ups/MS3/MS2)
      ScanCommand ms1_cmd = queue_.makeMS1();
      ms1_cmd.faims_cv = faims_cv;
      ms1_cmd.scan_id = queue_.nextTrackingId();
      ms1_cmd.priority = 3;

      std::string ms1_id_str = ScanCommandQueue::encode(ms1_cmd.scan_id);
      std::snprintf(ms1_cmd.scan_description, 16, "%sS", ms1_id_str.c_str());

      std::cout << "[TRACK-CREATE] id=" << ms1_id_str << " ms_level=1 type=idle_ms1" << std::endl;

      // Push MS1 into priority-3 queue for next dequeue call
      queue_.push(ms1_cmd);

      out = agc_cmd;
      queue_.registerPending(out);
      writeScanCommandRow_(out);
      return 1;
    }
  }

  int FLASHIda::getNextTrackingId()
  {
    return queue_.nextTrackingId();
  }

} // namespace OpenMS
