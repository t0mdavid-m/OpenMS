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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/MS3FragmentMatcher.h>
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
    // Empty log_dir = open nothing. Otherwise log_dir is the already-created, already-resolved
    // run folder and every stream is <log_dir>/<fixed basename>. The host owns the filesystem:
    // this constructor never creates a directory, so a folder that does not exist leaves every
    // stream closed and every writer a no-op.
    const std::string dir = rt_cfg.log_dir.empty() ? std::string{} : rt_cfg.log_dir + "/";
    if (!dir.empty())
    {
      ida_log_stream_.open(dir + kIdaLogName, std::ios::app);
    }
    if (!dir.empty())
    {
      commands_tsv_stream_.open(dir + kScanCommandsName, std::ios::app);
      if (commands_tsv_stream_.is_open())
      {
        commands_tsv_stream_ << "tracking_id\tscan_type\tms_level\tparent_tracking_id\t"
                             // P5: per-MS1-selection precursor identity (plain decimal); 0 for MS1/AGC/untracked.
                             << "precursor_id\t"
                             << "priority\tmono_mass\tcharge\tprecursor_mz\tisolation_width\t"
                             << "qscore\tcharge_cos\tcharge_snr\tiso_cos\tsnr\tcharge_score\t"
                             << "activation\tcollision_energy\thcd_energy\t"
                             << "reaction_time\treagent_max_it\treagent_agc_target\t"
                             << "ppm_error\tprecursor_intensity\tpeakgroup_intensity\t"
                             << "ion_type\tion_index\t"
                             // The wide clipped b/y fragment ProForma of the MS3 target being fired.
                             // "" on MS1/MS2/AGC (no fragment sub-sequence).
                             << "ms3_proteoform\t"
                             << "scan_description\t"  // E6: raw descriptor, asserted == what's sent to instrument
                             << "faims_cv\t"
                             // faims_cv alone is ambiguous: 0.0 is both "no FAIMS" and a
                             // legitimate compensation voltage, and ADR-0012 made cv_values: [0]
                             // -- FAIMS on at CV 0 -- expressible. Without this column a reader
                             // cannot tell that run from a FAIMS-off one. Sits between faims_cv
                             // and enqueue_ts deliberately: enqueue_ts must stay the trailing
                             // column (FLASHIda_LoggingFields_test asserts headers.back()), and
                             // every scan_commands column index pinned by that test is < 29.
                             << "faims_enabled\t"
                             // ADR-0026 decision 6: the acquired scan range. Without it a scan whose range was BOUND
                             // to its isolation window and one whose range was hand-overridden in ms_settings emit
                             // byte-identical rows, so the escape hatch is undetectable from a run folder. Logged for
                             // EVERY command, not only exploration ones. 0 is a real value here -- "unset, use the
                             // instrument method default" -- and is deliberately not suppressed.
                             // Placed between faims_enabled and enqueue_ts for the same forced reason as faims_enabled
                             // above: FLASHIda_LoggingFields_test pins headers.back() == "enqueue_ts" and pins six
                             // scan_commands column indices, all <= 30, so slots 31/32 are the only pair of positions
                             // that satisfies both.
                             << "first_mass\tlast_mass\t"
                             << "enqueue_ts\n";
        commands_tsv_stream_.flush();
      }
    }
    if (!dir.empty())
    {
      results_tsv_stream_.open(dir + kScanResultsName, std::ios::app);
      if (results_tsv_stream_.is_open())
      {
        // scan_results is a pure acquisition-event log: identification payload
        // (tag_count/matched_protein/proteoform_sequence/fragment_count) moved to
        // identification.tsv (proteoform -> scan_commands.ms3_proteoform); tic_coverage
        // moved to identification.tsv; tag_count dropped. (34 -> 29 columns.)
        results_tsv_stream_ << "tracking_id\tms_level\tparent_tracking_id\t"
                            << "commands_pushed\tchild_ids\trt\t"
                            << "duration_ms\tduration_received_ms\tqueue_duration_ms\tinstrument_duration_ms\tprocessing_duration_ms\t"
                            << "mass_count\tremaining_ratio\t"
                            // Per-scan identification YIELD, beside the other "how much did this scan
                            // return" columns. -1 = no tagger ran (MS1, MS3); 0 = ran, read nothing;
                            // >0 = real. tic_coverage also appears on identification.tsv, per ID row.
                            << "tag_count\tfragment_count\ttic_coverage\t"
                            << "exploration_group_id\texploration_metric\texploration_score\t"
                            << "variant_index\ttotal_variants\t"
                            << "winner_tracking_id\t"  // F5: group-completing exploration winner pointer
                            << "activation_type\tcollision_energy\treaction_time\t"
                            << "deconv_masses\tdeconv_qscores\tdeconv_charges\tdeconv_intensities\t"
                            << "resolve_ts\treceived_ts\tdequeue_ts\t"
                            // ADR-0038. Isobaric quantification, measured on the 'Q' quantification
                            // scan and nowhere else; every other row carries the sentinels. APPENDED,
                            // deliberately: every scan_results colIndex pinned in
                            // FLASHIda_LoggingFields_test's schema_column_counts is <= 21, so appending
                            // is the one placement that adds columns without invalidating a single pin
                            // -- only the count and headers.back() move.
                            << "quant_channels\tquant_condition_means\tquant_fold_change\tquant_verdict\n";
        results_tsv_stream_.flush();
      }
    }
    if (!dir.empty())
    {
      identification_tsv_stream_.open(dir + kIdentificationName, std::ios::app);
      if (identification_tsv_stream_.is_open())
      {
        identification_tsv_stream_ << "tracking_id\tscan_mode\tms_level\t"
                                   // P5: per-MS1-selection precursor identity (plain decimal).
                                   << "precursor_id\t"
                                   << "ms1_precursor_mass\tms2_precursor_ion\tproteoform\t"
                                   // C: FLASHExtender score of this scan's OWN match; -1 = winner-re-matched row (no own ID).
                                   << "flash_extender_score\t"
                                   // Sequence tags FLASHTagger read off THIS spectrum, beside the other
                                   // "how strong is this ID" scalar. -1 on every MS3 row: that is a
                                   // POLICY, not an absence of tagging (MS3 exploration variants are
                                   // tagged via computeFragmentMatch_). An MS3 tag count would measure
                                   // the sub-fragment spectrum, not the precursor's identifiability.
                                   << "tag_count\t"
                                   << "ms2_fragments\tms2_fragment_masses\t"
                                   << "ppm_offset\tcorrection_factor\t"
                                   << "ms1_precursor_mz\tms1_precursor_charge\t"
                                   << "ms2_precursor_mass\tms2_precursor_mz\tms2_precursor_charge\t"
                                   << "start_pos\tend_pos\t"
                                   << "ms3_fragments\tms3_fragment_masses\t"
                                   // I2: isolation-window width / over-window SNR / selected-charge intensity,
                                   // for the MS2 precursor and (on MS3 rows) the MS3 fragment precursor.
                                   << "ms2_isolation_width\tms2_window_snr\tms2_charge_intensity\t"
                                   << "ms3_isolation_width\tms3_window_snr\tms3_charge_intensity\t"
                                   // Fragment-mass table (MS2 & MS3 rows): per-scan theoretical + residual.
                                   << "theoretical_masses\tdiff_da\tdiff_ppm\t"
                                   // One PeakGroup::getQscore() per MATCHED fragment, index-aligned with
                                   // the block above and with ms2_fragments/ms3_fragments -- all five
                                   // lists iterate the same match.fragments vector. Distinct from
                                   // scan_results' deconv_qscores, which reports EVERY deconvolved mass
                                   // of the spectrum; this is the matched subset.
                                   << "fragment_qscores\t"
                                   // C2: MS3 per-ion fragment coverage (distinct backbone bonds / (L-1)); -1 on MS2.
                                   << "ms3_fragment_coverage\t"
                                   // Per-scan TIC / matched-fragment coverage (moved here from scan_results).
                                   << "tic_coverage\n";
        identification_tsv_stream_.flush();
      }
    }
    if (!dir.empty())
    {
      pooled_stream_.open(dir + kPooledIdentificationName, std::ios::app);
      if (pooled_stream_.is_open())
      {
        // P6: trajectory columns lead — trigger source and the tracking-id of the driving scan.
        pooled_stream_ << "trigger_scan_id\ttrigger\t"
                       // P5: per-MS1-selection precursor identity (the model key; plain decimal).
                       << "precursor_id\tupdate_index\t"
                       << "mono_mass\tproteoform\tflash_extender_score\t"
                       << "coverage_pct\tn_fragments\tlocalized_mods\tambiguous_mods\t"
                       // Fragment-mass table (grouped): masses + ion labels + measured/theoretical/residual,
                       // all aligned index-for-index with combined_ms2_frame_masses.
                       << "combined_ms2_frame_masses\tcombined_ms2_fragment_ions\t"
                       << "combined_measured_raw\tcombined_theoretical\tcombined_diff_da\tcombined_diff_ppm\t"
                       << "contributing_scan_ids\tnominal_mass\n";
        pooled_stream_.flush();
      }
    }
  }

  void IdaLogger::writeIDALogEntry(double rt, int tracking_id, int instrument_scan_number,
                                    const std::vector<ScanCommand>& ms2_commands,
                                    const std::vector<const PeakGroup*>& ms2_sources,
                                    const DeconvolvedSpectrum& all_peak_groups)
  {
    if (!ida_log_stream_.is_open()) return;
    // After the is_open() check, so a run with no log_dir pays nothing at all.
    std::lock_guard<std::mutex> lk(ida_log_mutex_);

    // The access-id token stays the TRACKING id: it is the join key to scan_commands.tsv,
    // scan_results.tsv, identification.tsv and pooled_identification.tsv, and nothing else on this
    // line carries it.
    const std::string tracking_id_str = ScanCommandQueue::encode(tracking_id);

    // "MS1 Scan#" is the INSTRUMENT's scan number -- what FLASHDeconv matches against the mzML
    // native id, and what the pre-port C# writer put here (ADR-0035). It used to be the tracking id,
    // which made the join structurally unsatisfiable rather than merely wrong: tracking ids count
    // 1,2,3... while instrument scan numbers run into the thousands, so past ~50 scans
    // findPrecursorPeakGroupsFormIdaLog_'s `iter->first < scan_number - 50` cutoff sits above every
    // key in the map and it returns before attempting a single isolation-window match.
    //
    // Falling back to the tracking id (rather than writing 0, or refusing the entry) keeps an
    // offline or synthetically-driven run producing exactly the log it produced before.
    const int logged_scan_number = (instrument_scan_number > 0) ? instrument_scan_number : tracking_id;
    if (instrument_scan_number <= 0 && !warned_missing_scan_number_)
    {
      warned_missing_scan_number_ = true;
      std::cout << "No instrument scan number supplied to processScan; ida.log 'MS1 Scan#' will "
                   "carry the tracking id instead. FLASHDeconv coupling on this log will not resolve."
                << std::endl;
    }

    ida_log_stream_ << "MS1 Scan# " << logged_scan_number
                    << " RT " << std::fixed << std::setprecision(4) << rt
                    << " (Access ID " << tracking_id_str << ") - "
                    << ms2_commands.size() << " targets\n";

    // Indexed, not range-for: ms2_sources is index-parallel to ms2_commands.
    for (size_t ci = 0; ci < ms2_commands.size(); ++ci)
    {
      const ScanCommand& cmd = ms2_commands[ci];
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

      // ChargeRange is the SPECIES' measured envelope (ADR-0035 decision 6), read off the PeakGroup
      // this command was actually built from -- ms2_sources is index-parallel to ms2_commands, filled
      // at push time in FLASHIda::processScan. Not a mass lookup into all_peak_groups: several
      // PeakGroups routinely share one mass within a survey (ADR-0036 split envelopes -- 48 of the
      // 1324 committed golden target lines sit on a species with 2-4 of them), and each carries a
      // different charge subset, so a lookup would silently pick one of several answers.
      // It reaches FLASHDeconvFeatureFile as columns 10-11 (z, Z) of *_ms2.feature, so a degenerate
      // [z-z] tells TopPIC every ida.log-sourced feature is single-charge.
      // Missing source => fall back to the trigger charge, i.e. the pre-ADR-0035 output rather than
      // a fabricated range. FLASHIda_LoggingFields_test's any_wider assertion is what catches that
      // fallback becoming the norm.
      // The well-formedness test is not defensive noise: PeakGroup's members default to
      // min_abs_charge_ = 0, max_abs_charge_ = -1, so a group whose range was never set would print
      // "ChargeRange=[0--1]" -- an inverted range that both readers would parse without complaint.
      // Every selected PeakGroup comes out of the deconvolution with a real range, so this arm is
      // unreachable today; it exists so the fallback means the same thing on every path.
      int z_lo = charge, z_hi = charge;
      if (ci < ms2_sources.size() && ms2_sources[ci] != nullptr)
      {
        auto [lo, hi] = ms2_sources[ci]->getAbsChargeRange();
        if (lo > 0 && hi >= lo) { z_lo = lo; z_hi = hi; }
      }

      // Mass= sets its OWN format. It used to say std::defaultfloat with no precision, inheriting
      // whatever the previous field left -- and under defaultfloat precision counts SIGNIFICANT
      // DIGITS, not decimals. The first target of an entry inherited the header's setprecision(4)
      // (-> "1.235e+04") and every later one the previous line's Features setprecision(6)
      // (-> "12383.3"), against an AllMass= of 12351.3933 / 12383.3180. std::fixed is redundant on
      // the second and later targets and load-bearing on the first, so it is set unconditionally.
      ida_log_stream_ << "Mass=" << std::fixed << std::setprecision(4) << cmd.mono_mass
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
                      << "\tChargeRange=[" << z_lo << "-" << z_hi << "]"
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
      // ADR-0038: 'F'/"followup" is retired. It named the scan quantification BOUGHT, which was
      // never measured; 'Q' names the scan quantification MEASURES.
      case 'Q': return "quantification";
      case 'C': return "conditional";
      case 'E': return "exploration";
      default: return "unknown";
    }
  }

  void IdaLogger::writeScanCommandRow(const ScanCommand& cmd, int precursor_id, const std::string& ms3_proteoform)
  {
    if (!commands_tsv_stream_.is_open()) return;
    // The one writer called from the instrument event thread. Its own lock, so it never queues
    // behind a processScan-side row and its flush.
    std::lock_guard<std::mutex> lk(commands_tsv_mutex_);

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

      // Co-isolation notches, ','-joined INSIDE this stage's ';'-group -- the same two-axis grammar
      // the wire uses, so scan_commands.tsv records what was actually commanded rather than only the
      // anchor window. Only the charge-dependent geometry gains the axis: collision energy and
      // activation are per fragmentation event and all notches of a stage share one.
      // No notches (every mode but multiplexed) emits nothing here, so these columns stay
      // byte-identical to the pre-notch format.
      auto stage_notches = notchesForStage(cmd, i);
      for (int n = 0; n < stage_notches.second; ++n)
      {
        std::ostringstream nmz_os, niw_os;
        nmz_os << stage_notches.first[n].precursor_mz;
        niw_os << stage_notches.first[n].isolation_width;
        charges += "," + std::to_string(stage_notches.first[n].charge_state);
        precursor_mzs += "," + nmz_os.str();
        iso_widths += "," + niw_os.str();
      }
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
    // Scan-range bounds, formatted through a LOCAL stream like sc()/sci() above and for the same reason:
    // std::fixed/setprecision are STICKY, and commands_tsv_stream_ is deliberately left at the stream default,
    // so setting them on it here would reformat faims_cv on every SUBSEQUENT row (0 -> 0.0000) and move goldens
    // for all 22 modes. The default is not sufficient on its own either: defaultfloat gives 6 SIGNIFICANT
    // figures, so an isolation-bound edge of 1236.625 prints as "1236.62" -- 0.005 Th gone, which is exactly
    // the resolution a reader needs to tell a bound window from a hand-override. fixed+4 holds 4 decimals at
    // any m/z in range.
    auto mz = [](double v) { std::ostringstream os; os << std::fixed << std::setprecision(4) << v; return os.str(); };
    commands_tsv_stream_ << id_str << "\t"
                         << scan_type << "\t"
                         << cmd.msn_level << "\t"
                         << parent_id << "\t"
                         << precursor_id << "\t"          // P5: per-MS1-selection precursor identity
                         << cmd.priority << "\t"
                         << sc(cmd.mono_mass, cmd.mono_mass_s1) << "\t"
                         << charges << "\t"
                         << precursor_mzs << "\t"
                         << iso_widths << "\t"
                         << sc(cmd.qscore, cmd.qscore_s1) << "\t"
                         << sc(cmd.charge_cos, cmd.charge_cos_s1) << "\t"
                         << sc(cmd.charge_snr, cmd.charge_snr_s1) << "\t"
                         << sc(cmd.iso_cos, cmd.iso_cos_s1) << "\t"
                         << sc(cmd.snr, cmd.snr_s1) << "\t"
                         << sc(cmd.charge_score, cmd.charge_score_s1) << "\t"
                         << activations << "\t"
                         << col_energies << "\t"
                         << sci(cmd.hcd_energy, cmd.hcd_energy_s1) << "\t"
                         << reaction_times << "\t"
                         << reagent_max_its << "\t"
                         << reagent_agc_targets << "\t"
                         << sc(cmd.ppm_error, cmd.ppm_error_s1) << "\t"
                         << sc(cmd.precursor_intensity, cmd.precursor_intensity_s1) << "\t"
                         << sc(cmd.peakgroup_intensity, cmd.peakgroup_intensity_s1) << "\t"
                         << ion_type << "\t"
                         << ion_index << "\t"
                         << ms3_proteoform << "\t"        // wide MS3 fragment ProForma ("" for MS1/MS2/AGC)
                         << cmd.scan_description << "\t"  // E6: raw descriptor (tab/newline-free by construction)
                         << cmd.faims_cv << "\t"
                         << cmd.faims_enabled << "\t"     // 0/1; disambiguates CV 0 from no FAIMS
                         << mz(cmd.first_mass) << "\t"    // ADR-0026 d6: 0 = unset (instrument method default), a real value
                         << mz(cmd.last_mass) << "\t"
                         << cmd.enqueue_timestamp_ms << "\n";
    commands_tsv_stream_.flush();
  }

  void IdaLogger::writeScanResultRow(const ScanRowDescriptor& row)
  {
    if (!results_tsv_stream_.is_open()) return;
    std::lock_guard<std::mutex> lk(results_tsv_mutex_);

    // Alias the descriptor fields so the row-building body below is unchanged (byte-identical output).
    const std::string& tracking_id = row.tracking_id;
    const int ms_level = row.ms_level;
    const double rt = row.rt;
    const int mass_count = row.mass_count;
    const int commands_pushed = row.commands_pushed;
    const std::vector<std::string>& child_ids = row.child_ids;
    const uint64_t enqueue_ts = row.enqueue_ts;
    const uint64_t dequeue_ts = row.dequeue_ts;
    const uint64_t received_ts = row.received_ts;
    const DeconvolvedSpectrum* deconv_spectrum = row.deconv_spectrum;
    const std::string& parent_tracking_id = row.parent_tracking_id;
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
                        << parent_tracking_id << "\t"
                        << commands_pushed << "\t"
                        << child_str << "\t"
                        << rt << "\t"
                        << duration << "\t"
                        << duration_received << "\t"
                        << queue_duration << "\t"
                        << instrument_duration << "\t"
                        << processing_duration << "\t"
                        << mass_count << "\t"
                        << remaining_ratio << "\t"
                        // Identification yield. Emitted raw: results_tsv_stream_ is the one IdaLogger
                        // stream deliberately never given a precision, and both fixed/setprecision are
                        // sticky, so setting one here would reformat rt, remaining_ratio,
                        // exploration_score and the next row's deconv block. Two are ints anyway.
                        << row.tag_count << "\t"
                        << row.fragment_count << "\t"
                        << row.tic_coverage << "\t"
                        << exploration_group_id << "\t"
                        << exploration_metric << "\t"
                        << exploration_score << "\t"
                        << variant_index << "\t"
                        << total_variants << "\t"
                        // F5: winner pointer — non-empty ONLY on the group-completing exploration row
                        // (the winning variant's encoded id); "" everywhere else.
                        << winner_tracking_id << "\t";
    results_tsv_stream_ << activation_type << "\t"
                        << collision_energy << "\t"
                        << reaction_time << "\t";

    // Deconvolved masses, then the per-mass score, then the per-charge breakdown: one ';'-group per
    // PeakGroup in each column, with that group's observed charges and their OWN intensities ','-joined
    // inside, index-aligned.
    //
    //   deconv_masses       12358.31;15234.11
    //   deconv_qscores      0.83;0.44
    //   deconv_charges      12,14;9,10
    //   deconv_intensities  3e5,8e5;1e5,2e5
    //
    // deconv_qscores is the PER-MASS value: PeakGroupScoring's logistic evaluated at the representative
    // charge (getRepAbsCharge(), the max-ChargeSNR charge). It is therefore an ELEMENT of the per-charge
    // qscore set, not an aggregate of it -- which is also why it needs no ','-axis while charges and
    // intensities do -- and for a mass that was selected it equals that command's scan_commands.qscore.
    //
    // It is written unconditionally, on every MS level, and deliberately is NOT the config-dependent
    // selection score (consider_all_charges swaps in getBestQScore() at the selection site). One
    // definition everywhere means a reader can compare the column across rows, runs and configs without
    // first knowing which selection knobs were set. The one caveat is provenance, not definition:
    // exploration-variant rows are deconvolved by Exploration's SECOND Deconvolution instance, built
    // from a different tolerance list, so their scores are not directly comparable to production rows.
    //
    // deconv_intensities used to be getIntensity() -- the PeakGroup TOTAL, summed across charge states
    // -- and the charge information was only the [min,max] range, so the log could say how much signal
    // a mass carried but never how it was distributed across the envelope. That distribution is the
    // input to charge-state co-isolation: it is what decides which charges clear the SNR gate.
    // min/max are dropped as derivable from the charge list, so that reshape left three columns where
    // there had been four; deconv_qscores makes it four again for the unrelated reason above.
    //
    // A charge is listed only when it actually carries intensity, so the group is the OBSERVED
    // envelope rather than every integer in the range.
    if (deconv_spectrum != nullptr && deconv_spectrum->size() > 0)
    {
      for (size_t i = 0; i < deconv_spectrum->size(); i++)
      {
        if (i > 0) results_tsv_stream_ << ";";
        results_tsv_stream_ << (*deconv_spectrum)[i].getMonoMass();
      }
      results_tsv_stream_ << "\t";
      // No std::fixed / setprecision here, deliberately: results_tsv_stream_ is the one IdaLogger stream
      // that is never given a precision, and both are STICKY -- setting them for this column would silently
      // reformat rt, remaining_ratio, exploration_score and the NEXT row's deconv_masses.
      for (size_t i = 0; i < deconv_spectrum->size(); i++)
      {
        if (i > 0) results_tsv_stream_ << ";";
        results_tsv_stream_ << (*deconv_spectrum)[i].getQscore();
      }
      results_tsv_stream_ << "\t";
      for (size_t i = 0; i < deconv_spectrum->size(); i++)
      {
        if (i > 0) results_tsv_stream_ << ";";
        auto [lo, hi] = (*deconv_spectrum)[i].getAbsChargeRange();
        bool first = true;
        for (int z = lo; z <= hi; ++z)
        {
          if ((*deconv_spectrum)[i].getChargeIntensity(z) <= 0) continue;
          if (!first) results_tsv_stream_ << ",";
          results_tsv_stream_ << z;
          first = false;
        }
      }
      results_tsv_stream_ << "\t";
      for (size_t i = 0; i < deconv_spectrum->size(); i++)
      {
        if (i > 0) results_tsv_stream_ << ";";
        auto [lo, hi] = (*deconv_spectrum)[i].getAbsChargeRange();
        bool first = true;
        for (int z = lo; z <= hi; ++z)
        {
          const float ci = (*deconv_spectrum)[i].getChargeIntensity(z);
          if (ci <= 0) continue;
          if (!first) results_tsv_stream_ << ",";
          results_tsv_stream_ << ci;
          first = false;
        }
      }
    }
    else
    {
      // One tab per gap BETWEEN the four deconv columns; the trailing tab before resolve_ts is emitted
      // below for both branches alike. Four columns => three tabs (it was two while there were three).
      results_tsv_stream_ << "\t\t\t";
    }
    // ADR-0038 quant block. Joined HERE rather than by the producer, so processScan hands over
    // values and never learns the delimiter. Built in a LOCAL ostringstream: std::fixed /
    // setprecision are sticky, and results_tsv_stream_ is deliberately left at the stream default
    // (see the notes above), which is also what renders the -1.0 sentinels as "-1".
    // Empty vector -> empty cell, which is the "not measured" sentinel for every non-'Q' row.
    auto joinQuant = [](const std::vector<double>& v) {
      std::ostringstream os;
      for (size_t i = 0; i < v.size(); ++i) { if (i) os << ';'; os << v[i]; }
      return os.str();
    };
    const std::string quant_channels_str = joinQuant(row.quant_channels);
    const std::string quant_means_str    = joinQuant(row.quant_condition_means);

    results_tsv_stream_ << "\t" << resolve_ts
                        << "\t" << received_ts
                        << "\t" << dequeue_ts
                        // ADR-0038. No stream manipulator here on purpose: std::fixed/setprecision are
                        // STICKY on this stream (see the notes above), and results_tsv_stream_ is
                        // deliberately left at the default, which is what renders the -1.0 sentinel as
                        // "-1" exactly like remaining_ratio and exploration_score already do.
                        << "\t" << quant_channels_str
                        << "\t" << quant_means_str
                        << "\t" << row.quant_fold_change
                        << "\t" << row.quant_verdict << "\n";
    results_tsv_stream_.flush();
  }

  void IdaLogger::writeIdentificationRow(const IdRowDescriptor& row)
  {
    if (!identification_tsv_stream_.is_open()) return;
    std::lock_guard<std::mutex> lk(identification_tsv_mutex_);

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

    // FRESH PER-SCAN ambiguity (MS3): the identification leaf shows what THIS scan's ion looks like. Its
    // localizer (narrowFragmentPTMSites) brackets each mod over THIS scan's EQUIVALENT (full-protein) ions --
    // seeded WIDE over the fragment region [1,L] and tightened inward with the flip/mod-aware verdict shared
    // with pooled Pass B -- then MERGES any mods this scan cannot separate into one summed shift. So the leaf
    // reflects the per-scan evidence exactly: its ranges can be WIDER than the cumulative pooled log AND wider
    // than the scan_commands a-priori base (which is rendered pre-acquisition and cannot know the MS3 ions --
    // there is no longer a leaf-subset-of-scan_commands guarantee). `wide_sites` (from the shared
    // fragmentProFormaSites, region-frame RENDER context reconstructed from row.ctx) is passed only as the
    // a-priori classification reference, never as an output clamp. This writer emits ONLY identification.tsv;
    // pooled (ProteoformTracker) + scan_commands are untouched. MS2 rows keep their range as-is (narrowing is
    // an MS3-only, equivalent-frame notion).
    std::vector<FragmentAnalysis::PTMSite> narrowed_sites;
    const std::vector<FragmentAnalysis::PTMSite>* eff_ptm_sites = &id_ptm_sites;
    std::string id_render_seq = id_proteoform;   // default: the fragment sub-sequence carried by match
    if (ms_level == 3 && match_has_proteoform)
    {
      MS3FragmentMatcher::ProteoformContext render_ctx;
      render_ctx.region_start = ctx.start_pos;   // render context (region frame) from the cached MS2Context
      render_ctx.region_end = ctx.end_pos;
      render_ctx.ptm_sites = ctx.ptm_sites;
      std::string wide_seq;
      std::vector<FragmentAnalysis::PTMSite> wide_sites;
      if (MS3FragmentMatcher::fragmentProFormaSites(config_.characterization().protein_sequence, render_ctx,
              ctx.fragment_ion_type, ctx.fragment_ion_index, wide_seq, wide_sites))
      {
        narrowed_sites = FragmentAnalysis::narrowFragmentPTMSites(
            wide_sites, static_cast<int>(wide_seq.size()), match.fragments);   // THIS scan's ions only
        eff_ptm_sites = &narrowed_sites;
        id_render_seq = wide_seq;                // render the SAME fragment scan_commands renders
      }
      else
      {
        // Defensive fallback (render ctx not populated, e.g. a degenerate row): previous behaviour.
        narrowed_sites = FragmentAnalysis::narrowFragmentPTMSites(
            match.ptm_sites, static_cast<int>(match.proteoform_sequence.size()), match.fragments);
        eff_ptm_sites = &narrowed_sites;
      }
    }

    std::string proforma = FragmentAnalysis::toProForma(id_render_seq, *eff_ptm_sites);

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
    // One qscore per MATCHED fragment, built in this same loop so it is aligned with the mass table by
    // construction rather than by convention. Its OWN precision is set explicitly: this stream is
    // fixed/setprecision(4) but that state is re-established several times around the output statement,
    // and a qscore is a [0,1] probability, so 4 decimals is both sufficient and what its neighbours use.
    std::ostringstream frag_qscores;
    frag_qscores << std::fixed << std::setprecision(4);
    for (size_t i = 0; i < match.fragments.size(); ++i)
    {
      if (i > 0) { frag_theo << ";"; frag_diff_da << ";"; frag_diff_ppm << ";"; frag_qscores << ";"; }
      frag_theo << match.fragments[i].theoretical_mass;
      frag_diff_da << match.fragments[i].diff_da;
      frag_diff_ppm << match.fragments[i].diff_ppm;
      frag_qscores << match.fragments[i].qscore;
    }

    std::string precursor_ion;
    if (ms_level == 3 && ctx.fragment_ion_type != '\0')
      precursor_ion = std::string(1, ctx.fragment_ion_type)
                      + std::to_string(ctx.fragment_ion_index);

    identification_tsv_stream_
      << tracking_id << "\t"
      << scan_mode << "\t"
      << ms_level << "\t"
      // P5: per-MS1-selection precursor identity (plain decimal).
      << row.precursor_id << "\t"
      // First float on this stream: establish std::fixed + setprecision(4) here (KEEP std::fixed explicit).
      << std::fixed << std::setprecision(4) << ctx.ms1_precursor_mass << "\t"
      << precursor_ion << "\t"
      << proforma << "\t"
      // C: FLASHExtender score of this scan's own match (-1 = winner-re-matched row, no own ID). fixed+4 in effect.
      << row.flash_extender_score << "\t"
      // -1 on every MS3 row by POLICY, not because nothing tagged: MS3 exploration variants ARE tagged
      // (computeFragmentMatch_), but a tag count there measures the sub-fragment spectrum rather than
      // the precursor's identifiability, which is not what this column means. On MS2 rows 0 is a real
      // reading -- the tagger ran and read nothing -- and stays distinct from -1.
      << (ms_level == 3 ? -1 : match.tag_count) << "\t"
      << ms2_frags.str() << "\t"
      << ms2_masses.str() << "\t"
      << std::setprecision(2) << match.ppm_offset << "\t"
      << std::setprecision(8) << match.correction_factor << "\t"
      // Re-establish setprecision(4) (sticky for every float below).
      << std::setprecision(4) << ctx.ms1_precursor_mz << "\t"
      << ctx.ms1_precursor_charge << "\t"
      << (ms_level == 3 ? ctx.fragment_mass : 0.0) << "\t"
      << (ms_level == 3 ? ctx.fragment_mz : 0.0) << "\t"
      << (ms_level == 3 ? ctx.fragment_charge : 0) << "\t"
      << id_region_start << "\t"
      << id_region_end << "\t"
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
      // Fragment-mass table: per-scan theoretical + residual (Da, ppm) for MS2 AND MS3 fragments.
      << frag_theo.str() << "\t" << frag_diff_da.str() << "\t" << frag_diff_ppm.str() << "\t"
      // Per-matched-fragment deconvolution confidence, aligned with the three columns above.
      << frag_qscores.str() << "\t"
      // C2: MS3 per-ion fragment coverage (distinct backbone bonds / (L-1)); -1 on MS2. std::fixed setprecision(4) in effect.
      << (ms_level == 3 ? match.ms3_fragment_coverage : -1.0f) << "\t"
      // Per-scan TIC / matched-fragment coverage (moved from scan_results); the scan's actual value at this row.
      << row.tic_coverage << "\n";
    identification_tsv_stream_.flush();
  }

  void IdaLogger::writePooledModelRow(const PooledModelDescriptor& r)
  {
    if (!pooled_stream_.is_open()) return;
    std::lock_guard<std::mutex> lk(pooled_mutex_);

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

    pooled_stream_ << r.trigger_scan_id << "\t"
                   << r.trigger << "\t"
                   // P5: per-MS1-selection precursor identity (the model key).
                   << r.precursor_id << "\t"
                   << r.update_index << "\t"
                   << std::fixed << std::setprecision(4) << r.mono_mass << "\t"
                   << r.proforma << "\t"
                   << std::setprecision(4) << r.score << "\t"
                   << std::setprecision(4) << r.coverage_pct << "\t"
                   << r.n_fragments << "\t"
                   << loc_str << "\t"
                   << amb_str << "\t"
                   // Fragment-mass table (grouped): masses | ion labels | measured | theoretical | diff_da | diff_ppm.
                   << masses_str << "\t" << ions_str << "\t"
                   << measured_str << "\t" << theo_str << "\t" << diff_da_str << "\t" << diff_ppm_str << "\t"
                   << scan_ids_str << "\t"
                   << r.nominal_mass << "\n";
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
      int scan = -1;   // no MS1 header seen yet; a Mass= line before one is not attributable
      float mass, charge, w1, w2, qscore, pint, mint, z1, z2;
      float features[6];
      while (std::getline(instream, line))
      {
        // " 0 targets" WITH THE LEADING SPACE. The bare substring also matches "10 targets" and
        // "20 targets", and command fan-out (one MS2 config x CE-sweep variant x charge under
        // `separate`) makes counts >= 10 routine -- the committed separate_charges golden alone
        // carries eight "- 10 targets" headers and one "- 20 targets". Each was silently dropped,
        // and its Mass= rows then landed under the PREVIOUS header's scan key.
        if (line.find(" 0 targets") != line.npos) { continue; }
        if (line.hasPrefix("MS1"))
        {
          Size st = line.find("MS1 Scan# ") + 10;
          Size ed = line.find(' ', st);
          String n = line.substr(st, ed);
          scan = atoi(n.c_str());
          precursor_map_for_real_time_acquisition[scan]
            = std::vector<std::vector<float>>(); //// ms1 scan -> mass, charge ,score, mz range, precursor int, mass int, color
        }
        // A Mass= line before ANY header is not attributable to a scan. `scan` used to be
        // uninitialised here, so such a row was filed under whatever the stack held. Must sit AFTER
        // the MS1 block (which assigns `scan`) and BEFORE the Mass block. "AllMass=" does not have
        // the prefix "Mass", so this cannot swallow the AllMass line.
        if (line.hasPrefix("Mass") && scan < 0) { continue; }
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
