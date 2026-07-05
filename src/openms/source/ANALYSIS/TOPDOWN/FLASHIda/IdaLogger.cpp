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


#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/IdaLogger.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace OpenMS
{

  IdaLogger::IdaLogger(const Config& config) :
      config_(config)
  {
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
                             << "reaction_time\treagent_max_it\treagent_agc_target\t"
                             << "scan_description\t"  // E6: raw descriptor, asserted == what's sent to instrument
                             // P5: per-MS1-selection precursor identity (plain decimal); appended LAST so the
                             // existing positional column indices don't shift. 0 for MS1/AGC/untracked.
                             << "precursor_id\n";
        commands_tsv_stream_.flush();
      }
    }
    if (!rt_cfg.scan_results_path.empty())
    {
      results_tsv_stream_.open(rt_cfg.scan_results_path, std::ios::app);
      if (results_tsv_stream_.is_open())
      {
        results_tsv_stream_ << "tracking_id\tms_level\tresolve_ts\tduration_ms\treceived_ts\tduration_received_ms\trt\t"
                            << "mass_count\tcommands_pushed\tchild_ids\t"
                            << "tag_count\tmatched_protein\tproteoform_sequence\t"
                            << "tic_coverage\tfragment_count\t"
                            << "exploration_group_id\texploration_metric\t"
                            << "variant_index\ttotal_variants\t"
                            << "collision_energy\texploration_score\tremaining_ratio\t"
                            << "activation_type\treaction_time\t"
                            << "deconv_masses\tdeconv_intensities\tdeconv_min_charge\tdeconv_max_charge\tparent_tracking_id\t"
                            << "dequeue_ts\tqueue_duration_ms\tinstrument_duration_ms\tprocessing_duration_ms\t"
                            << "winner_tracking_id\n";  // F5: trailing winner pointer (33 -> 34 columns)
        results_tsv_stream_.flush();
      }
    }
    if (!rt_cfg.identification_path.empty())
    {
      identification_tsv_stream_.open(rt_cfg.identification_path, std::ios::app);
      if (identification_tsv_stream_.is_open())
      {
        identification_tsv_stream_ << "ms_level\tscan_mode\t"
                                   << "tracking_id\tproteoform\tstart_pos\tend_pos\t"
                                   << "ppm_offset\tcorrection_factor\t"
                                   << "ms1_precursor_mass\tms1_precursor_mz\tms1_precursor_charge\t"
                                   << "ms2_precursor_ion\tms2_precursor_mass\tms2_precursor_mz\tms2_precursor_charge\t"
                                   << "ms2_fragments\tms2_fragment_masses\t"
                                   << "ms3_fragments\tms3_fragment_masses\t"
                                   // I2: isolation-window width / over-window SNR / selected-charge intensity,
                                   // for the MS2 precursor and (on MS3 rows) the MS3 fragment precursor.
                                   << "ms2_isolation_width\tms2_window_snr\tms2_charge_intensity\t"
                                   << "ms3_isolation_width\tms3_window_snr\tms3_charge_intensity\t"
                                   // P5: per-MS1-selection precursor identity (plain decimal).
                                   << "precursor_id\t"
                                   // Fragment-mass table (MS2 & MS3 rows): per-scan theoretical + residual; appended LAST.
                                   << "theoretical_masses\tdiff_da\tdiff_ppm\t"
                                   // C2: MS3 per-ion fragment coverage (distinct backbone bonds / (L-1)); appended LAST, -1 on MS2.
                                   << "ms3_fragment_coverage\n";
        identification_tsv_stream_.flush();
      }
    }
    if (!rt_cfg.pooled_identification_path.empty())
    {
      pooled_stream_.open(rt_cfg.pooled_identification_path, std::ios::app);
      if (pooled_stream_.is_open())
      {
        pooled_stream_ << "nominal_mass\tmono_mass\tproteoform\tflash_extender_score\t"
                       << "coverage_pct\tn_fragments\tlocalized_mods\tambiguous_mods\t"
                       << "contributing_scan_ids\t"
                       // Fragment-mass table (grouped): masses + ion labels + measured/theoretical/residual,
                       // all aligned index-for-index with combined_ms2_frame_masses.
                       << "combined_ms2_frame_masses\tcombined_ms2_fragment_ions\t"
                       << "combined_measured\tcombined_theoretical\tcombined_diff_da\tcombined_diff_ppm\t"
                       << "update_index\t"
                       // P5: per-MS1-selection precursor identity (the model key; plain decimal).
                       << "precursor_id\t"
                       // P6: trajectory columns — trigger source and the tracking-id of the driving scan.
                       << "trigger\ttrigger_scan_id\n";
        pooled_stream_.flush();
      }
    }
  }

  void IdaLogger::writeIDALogEntry(double rt, int scan_number,
                                    const std::string& tracking_id,
                                    const std::vector<ScanCommand>& ms2_commands,
                                    const DeconvolvedSpectrum& all_peak_groups)
  {
    if (!ida_log_stream_.is_open()) return;

    // MS1 header line. scan_number = the decoded tracking_id (a distinct, increasing per-command
    // base-94 counter) so multi-scan logs no longer collapse to map key 0 (B1). NOTE: it is NOT the
    // instrument scan index (that would need a bridge change); it is unique-per-scan, which is the fix.
    ida_log_stream_ << "MS1 Scan# " << scan_number
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

  std::string IdaLogger::scanTypeFromDescription_(const ScanCommand& cmd)
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

  void IdaLogger::writeScanCommandRow(const ScanCommand& cmd, int precursor_id)
  {
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
    // Stage-less commands (MS1 / AGC) populate no per-stage fields. Emit explicit
    // placeholders so the row still serializes the full column count: a tab-split
    // parser drops a trailing empty field, which would leave such rows one column
    // short. Mirrors the empty-field guard in the ScanResults row writer.
    if (cmd.num_stages == 0)
    {
      charges = precursor_mzs = iso_widths = col_energies = reaction_times =
          reagent_max_its = reagent_agc_targets = "0";
      activations = "none";
    }

    std::string parent_id(cmd.parent_scan_id);
    std::string ion_type;
    int ion_index = 0;
    if (cmd.msn_level == 3 && std::strlen(cmd.scan_description) > 4)
    {
      std::string desc(cmd.scan_description);
      // rfind: the charge delimiter is the LAST '@'. The 3-char base-94 tracking id can itself contain '@'
      // (0x40, e.g. id "!#@"); find() would lock onto the id's '@' and mis-parse the ion (logging the marker
      // 'R' as ion_type + a mass digit as ion_index). The mass token has no '@', so rfind('@') is exactly the
      // charge delimiter. (Matches the test-side decodeTrailingIon, which also uses rfind.)
      auto at_pos = desc.rfind('@');
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

    // MS3 scoring columns are 2-stage 'MS2-precursor;fragment'; sc/sci append ';v1' only for MS3.
    // MS1/MS2/followup rows (msn_level != 3) emit stage-0 only -> byte-identical to before.
    auto sc = [&](double v0, double v1) {
      std::ostringstream os; os << v0; if (cmd.msn_level == 3) { os << ';' << v1; } return os.str();
    };
    auto sci = [&](int v0, int v1) {
      std::ostringstream os; os << v0; if (cmd.msn_level == 3) { os << ';' << v1; } return os.str();
    };
    commands_tsv_stream_ << id_str << "\t"
                         << cmd.msn_level << "\t"
                         << scan_type << "\t"
                         << cmd.enqueue_timestamp_ms << "\t"
                         << cmd.priority << "\t"
                         << cmd.faims_cv << "\t"
                         << sc(cmd.mono_mass, cmd.mono_mass_s1) << "\t"
                         << charges << "\t"
                         << precursor_mzs << "\t"
                         << iso_widths << "\t"
                         << col_energies << "\t"
                         << activations << "\t"
                         << sc(cmd.qscore, cmd.qscore_s1) << "\t"
                         << sc(cmd.charge_cos, cmd.charge_cos_s1) << "\t"
                         << sc(cmd.charge_snr, cmd.charge_snr_s1) << "\t"
                         << sc(cmd.iso_cos, cmd.iso_cos_s1) << "\t"
                         << sc(cmd.snr, cmd.snr_s1) << "\t"
                         << sc(cmd.charge_score, cmd.charge_score_s1) << "\t"
                         << sc(cmd.ppm_error, cmd.ppm_error_s1) << "\t"
                         << sc(cmd.precursor_intensity, cmd.precursor_intensity_s1) << "\t"
                         << sc(cmd.peakgroup_intensity, cmd.peakgroup_intensity_s1) << "\t"
                         << sci(cmd.hcd_energy, cmd.hcd_energy_s1) << "\t"
                         << parent_id << "\t"
                         << ion_type << "\t"
                         << ion_index << "\t"
                         << reaction_times << "\t"
                         << reagent_max_its << "\t"
                         << reagent_agc_targets << "\t"
                         << cmd.scan_description << "\t"  // E6: raw descriptor (tab/newline-free by construction)
                         << precursor_id << "\n";         // P5: per-MS1-selection precursor identity (appended LAST)
    commands_tsv_stream_.flush();
  }

  void IdaLogger::writeScanResultRow(const ScanRowDescriptor& row)
  {
    if (!results_tsv_stream_.is_open()) return;

    // Alias the descriptor fields so the row-building body below is unchanged (byte-identical output).
    const std::string& tracking_id = row.tracking_id;
    const int ms_level = row.ms_level;
    const double rt = row.rt;
    const int mass_count = row.mass_count;
    const int commands_pushed = row.commands_pushed;
    const std::vector<std::string>& child_ids = row.child_ids;
    const int tag_count = row.tag_count;
    const std::string& matched_protein = row.matched_protein;
    const std::string& proteoform_sequence = row.proteoform_sequence;
    const uint64_t enqueue_ts = row.enqueue_ts;
    const uint64_t dequeue_ts = row.dequeue_ts;
    const uint64_t received_ts = row.received_ts;
    const DeconvolvedSpectrum* deconv_spectrum = row.deconv_spectrum;
    const std::string& parent_tracking_id = row.parent_tracking_id;
    const float tic_coverage = row.tic_coverage;
    const int fragment_count = row.fragment_count;
    const int exploration_group_id = row.exploration_group_id;
    const int exploration_metric = row.exploration_metric;
    const int variant_index = row.variant_index;
    const int total_variants = row.total_variants;
    const std::string& collision_energy = row.collision_energy;
    const double exploration_score = row.exploration_score;
    const double remaining_ratio = row.remaining_ratio;
    const std::string& activation_type = row.activation_type;
    const std::string& reaction_time = row.reaction_time;
    const std::string& winner_tracking_id = row.winner_tracking_id;

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
      // The child_ids separator MUST be outside the 0x21-0x7E tracking-id alphabet
      // (encoded ids may contain ';', ',', '"', etc.); space (0x20) is the only printable
      // char the alphabet excludes, so ids can never collide with it. Readers split on ' '.
      if (i > 0) child_str += " ";
      child_str += child_ids[i];
    }

    results_tsv_stream_ << tracking_id << "\t"
                        << ms_level << "\t"
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
      results_tsv_stream_ << "\t\t\t";
    }
    results_tsv_stream_ << "\t" << parent_tracking_id
                        << "\t" << dequeue_ts
                        << "\t" << queue_duration
                        << "\t" << instrument_duration
                        << "\t" << processing_duration
                        // F5: trailing winner pointer — non-empty ONLY on the group-completing exploration
                        // row (the winning variant's encoded id); "" everywhere else. Appended LAST so the
                        // C# comparer's fixed column indices don't shift (scan_results 33 -> 34 columns).
                        << "\t" << winner_tracking_id << "\n";
    results_tsv_stream_.flush();
  }

  void IdaLogger::writeIdentificationRow(const IdRowDescriptor& row)
  {
    if (!identification_tsv_stream_.is_open()) return;

    // Alias the descriptor fields so the row-building body below is unchanged (byte-identical output).
    const std::string& tracking_id = row.tracking_id;
    const int ms_level = row.ms_level;
    const char scan_mode = row.scan_mode;
    const Exploration::MS2Context& ctx = row.ctx;
    const FragmentAnalysis::ProteoformMatch& match = row.match;

    if (match.fragments.empty()) return;

    // CONTRACT: `match` = the identified species (MS2 = the MS2 proteoform; MS3 = the FRAGMENT
    // sub-sequence + its clipped mods, filled by calibrateAndScore), `ctx` = the acquisition context
    // (MS1/MS2/MS3 precursor identity + isolation). Source proteoform/PTMs/region from `match`; fall back
    // to the parent `ctx` only if `match` is unpopulated (defensive — a path that produced no match).
    const bool match_has_proteoform = !match.proteoform_sequence.empty();
    const std::string& id_proteoform = match_has_proteoform ? match.proteoform_sequence : ctx.proteoform_sequence;
    const std::vector<FragmentAnalysis::PTMSite>& id_ptm_sites = match_has_proteoform ? match.ptm_sites : ctx.ptm_sites;
    // region: the fragment sub-range (MS3) / proteoform region (MS2) from match.region_start/end (0-based,
    // exclusive end); fall back to the parent ctx range only if match has no region (region_start < 0).
    const bool use_match_region = (match.region_start >= 0);
    const int id_region_start = use_match_region ? match.region_start : ctx.start_pos;
    const int id_region_end   = use_match_region ? match.region_end   : ctx.end_pos;

    // ISSUE(N): the MS3 identification leaf carries THIS scan's OWN PTM narrowing (per-scan evidence),
    // not the parent-wide clip. Narrow a LOCAL copy of the PTM sites from the match's own matched
    // sub-fragments; this writer emits ONLY identification.tsv, so pooled (ProteoformTracker, seeded from
    // the MS2 winner) and scan_results (rendered from the parent MS2 context) are untouched. This single
    // site covers all three MS3 identification-row sinks (regular 'R' + exploration-'E' primary + winner).
    // MS2 rows keep their proteoform range as-is (narrowing is an MS3-only, fragment-frame notion).
    std::vector<FragmentAnalysis::PTMSite> narrowed_sites;
    const std::vector<FragmentAnalysis::PTMSite>* eff_ptm_sites = &id_ptm_sites;
    if (ms_level == 3 && match_has_proteoform)
    {
      narrowed_sites = FragmentAnalysis::narrowFragmentPTMSites(
          match.ptm_sites, static_cast<int>(match.proteoform_sequence.size()), match.fragments);
      eff_ptm_sites = &narrowed_sites;
    }

    std::string proforma = FragmentAnalysis::toProForma(id_proteoform, *eff_ptm_sites);

    std::ostringstream ms2_frags, ms2_masses;
    ms2_frags << std::fixed << std::setprecision(4);
    ms2_masses << std::fixed << std::setprecision(4);

    if (ms_level == 2)
    {
      for (size_t i = 0; i < match.fragments.size(); ++i)
      {
        if (i > 0) { ms2_frags << ";"; ms2_masses << ";"; }
        ms2_frags << match.fragments[i].ion_type << match.fragments[i].ion_index;
        ms2_masses << match.fragments[i].observed_mass;
      }
    }
    else
    {
      for (size_t i = 0; i < match.fragments.size(); ++i)
      {
        if (i > 0) { ms2_frags << ";"; ms2_masses << ";"; }
        ms2_frags << match.fragments[i].equiv_type << match.fragments[i].equiv_index;
        ms2_masses << match.fragments[i].adjusted_mass;
      }
    }

    std::ostringstream ms3_frags, ms3_masses;
    ms3_frags << std::fixed << std::setprecision(4);
    ms3_masses << std::fixed << std::setprecision(4);

    if (ms_level == 3)
    {
      for (size_t i = 0; i < match.fragments.size(); ++i)
      {
        if (i > 0) { ms3_frags << ";"; ms3_masses << ";"; }
        ms3_frags << match.fragments[i].ion_type << match.fragments[i].ion_index;   // raw MS3 sub-fragment (measured)
        ms3_masses << match.fragments[i].observed_mass;
      }
    }

    // Fragment-mass table theoretical + residual, per fragment, for BOTH levels (match.fragments): MS2 rows
    // carry the MS2 matcher's best_theo + diff, MS3 rows the equiv-frame theoretical + diff (calibrateAndScore).
    std::ostringstream frag_theo, frag_diff_da, frag_diff_ppm;
    frag_theo << std::fixed << std::setprecision(4);
    frag_diff_da << std::fixed << std::setprecision(4);
    frag_diff_ppm << std::fixed << std::setprecision(2);
    for (size_t i = 0; i < match.fragments.size(); ++i)
    {
      if (i > 0) { frag_theo << ";"; frag_diff_da << ";"; frag_diff_ppm << ";"; }
      frag_theo << match.fragments[i].theoretical_mass;
      frag_diff_da << match.fragments[i].diff_da;
      frag_diff_ppm << match.fragments[i].diff_ppm;
    }

    std::string precursor_ion;
    if (ms_level == 3 && ctx.fragment_ion_type != '\0')
      precursor_ion = std::string(1, ctx.fragment_ion_type)
                      + std::to_string(ctx.fragment_ion_index);

    identification_tsv_stream_
      << ms_level << "\t"
      << scan_mode << "\t"
      << tracking_id << "\t"
      << proforma << "\t"
      << id_region_start << "\t"
      << id_region_end << "\t"
      << std::fixed << std::setprecision(2) << match.ppm_offset << "\t"
      << std::setprecision(8) << match.correction_factor << "\t"
      << std::setprecision(4) << ctx.ms1_precursor_mass << "\t"
      << ctx.ms1_precursor_mz << "\t"
      << ctx.ms1_precursor_charge << "\t"
      << precursor_ion << "\t"
      << (ms_level == 3 ? ctx.fragment_mass : 0.0) << "\t"
      << (ms_level == 3 ? ctx.fragment_mz : 0.0) << "\t"
      << (ms_level == 3 ? ctx.fragment_charge : 0) << "\t"
      << ms2_frags.str() << "\t"
      << ms2_masses.str() << "\t"
      << ms3_frags.str() << "\t"
      << ms3_masses.str() << "\t"
      // I2: isolation-window reporting (std::fixed setprecision(4) still in effect). MS2 triplet always
      // written; MS3 triplet only on MS3 rows (0.0 otherwise), mirroring the fragment_* columns above.
      << ctx.ms2_isolation_width << "\t"
      << ctx.ms2_window_snr << "\t"
      << ctx.ms2_charge_intensity << "\t"
      << (ms_level == 3 ? ctx.ms3_isolation_width : 0.0) << "\t"
      << (ms_level == 3 ? ctx.ms3_window_snr : 0.0) << "\t"
      << (ms_level == 3 ? ctx.ms3_charge_intensity : 0.0) << "\t"
      // P5: per-MS1-selection precursor identity (plain decimal).
      << row.precursor_id << "\t"
      // Fragment-mass table: per-scan theoretical + residual (Da, ppm) for MS2 AND MS3 fragments; appended LAST.
      << frag_theo.str() << "\t" << frag_diff_da.str() << "\t" << frag_diff_ppm.str() << "\t"
      // C2: MS3 per-ion fragment coverage (distinct backbone bonds / (L-1)); -1 on MS2. std::fixed setprecision(4) in effect.
      << (ms_level == 3 ? match.ms3_fragment_coverage : -1.0f) << "\n";
    identification_tsv_stream_.flush();
  }

  void IdaLogger::writePooledModelRow(const PooledModelDescriptor& r)
  {
    if (!pooled_stream_.is_open()) return;

    // localized_mods and ambiguous_mods joined with ';' (no tracking-id collision risk: these are
    // human-readable strings, not encoded ids).
    std::string loc_str;
    for (size_t i = 0; i < r.localized_mods.size(); ++i)
    {
      if (i > 0) loc_str += ";";
      loc_str += r.localized_mods[i];
    }
    std::string amb_str;
    for (size_t i = 0; i < r.ambiguous_mods.size(); ++i)
    {
      if (i > 0) amb_str += ";";
      amb_str += r.ambiguous_mods[i];
    }
    // contributing_scan_ids joined with SPACE — the established delimiter that avoids the
    // tracking-id-alphabet ';' collision (same precedent as child_ids in writeScanResultRow).
    // Each id is base-94 encoded (ScanCommandQueue::encode) for consistency with every other
    // tracking-id in every stream — including the pooled trigger_scan_id (col 18) in this same row
    // and encode(cmd.scan_id) in writeScanCommandRow. Encoded ids are exactly 3 chars from the
    // 0x21-0x7E alphabet, so the space (0x20) delimiter never collides.
    std::string scan_ids_str;
    for (size_t i = 0; i < r.contributing_scan_ids.size(); ++i)
    {
      if (i > 0) scan_ids_str += " ";
      scan_ids_str += ScanCommandQueue::encode(r.contributing_scan_ids[i]);
    }
    // Fragment-mass table (grouped) — the aligned lists joined with ';'; ion labels are strings.
    auto joinDoubles = [](const std::vector<double>& v, int precision) -> std::string {
      std::ostringstream ss;
      ss << std::fixed << std::setprecision(precision);
      for (size_t i = 0; i < v.size(); ++i) { if (i > 0) ss << ";"; ss << v[i]; }
      return ss.str();
    };
    std::string ions_str;
    for (size_t i = 0; i < r.combined_ms2_fragment_ions.size(); ++i)
    { if (i > 0) ions_str += ";"; ions_str += r.combined_ms2_fragment_ions[i]; }
    const std::string masses_str   = joinDoubles(r.combined_masses, 4);      // combined_ms2_frame_masses (adjusted)
    const std::string measured_str = joinDoubles(r.combined_measured, 4);
    const std::string theo_str     = joinDoubles(r.combined_theoretical, 4);
    const std::string diff_da_str  = joinDoubles(r.combined_diff_da, 4);
    const std::string diff_ppm_str = joinDoubles(r.combined_diff_ppm, 2);

    pooled_stream_ << r.nominal_mass << "\t"
                   << std::fixed << std::setprecision(4) << r.mono_mass << "\t"
                   << r.proforma << "\t"
                   << std::setprecision(4) << r.score << "\t"
                   << std::setprecision(4) << r.coverage_pct << "\t"
                   << r.n_fragments << "\t"
                   << loc_str << "\t"
                   << amb_str << "\t"
                   << scan_ids_str << "\t"
                   // Fragment-mass table (grouped): masses | ion labels | measured | theoretical | diff_da | diff_ppm.
                   << masses_str << "\t" << ions_str << "\t"
                   << measured_str << "\t" << theo_str << "\t" << diff_da_str << "\t" << diff_ppm_str << "\t"
                   << r.update_index << "\t"
                   // P5: per-MS1-selection precursor identity (the model key).
                   << r.precursor_id << "\t"
                   // P6: trajectory columns — trigger source and driving scan id.
                   << r.trigger << "\t" << r.trigger_scan_id << "\n";
    pooled_stream_.flush();
  }

  std::map<int, std::vector<std::vector<float>>> IdaLogger::parseFLASHIdaLog(const String& in_log_file)
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

}
