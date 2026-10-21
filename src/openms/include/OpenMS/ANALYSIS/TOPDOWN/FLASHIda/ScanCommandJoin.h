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

#pragma once

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace OpenMS
{
  /**
   * @brief Joins a converted data file back to the FLASHIda acquisition that produced it (ADR-0046).
   *
   * A scan FLASHIda commanded is recognised in the mzML by its TRACKING ID -- the first three characters
   * of the scan description, which msconvert preserves as the spectrum meta value "scan description".
   * That id finds the scan's row in the run's scan_commands.tsv, and the row answers the two questions
   * an offline analysis cannot answer from the spectrum alone: WHICH mass was commanded, and -- through
   * parent_tracking_id -- WHICH survey the command was decided from. Commands queue, so that survey is
   * routinely not the one acquired just before the scan.
   *
   * Pure on purpose: parse() reads one file, join() maps data to data. Neither knows FLASHDeconv, which is
   * what lets both be tested from the FLASHIda side. MS2 only; rows of other levels keep id, parent and
   * level, because decision 5 checks the level of every joined scan and the parent walk needs the parents.
   */
  struct ScanCommandJoin
  {
    /// Locator tolerance: 10 ppm, plus 0.06 Da so that a mono_mass logged at six SIGNIFICANT digits --
    /// every log written before ADR-0046, 12351.4 for 12351.3933 -- still locates. For one species in one
    /// survey the engine's and FLASHDeconv's masses differ by a median of 0.003 ppm, so this is generous.
    static constexpr double LOCATOR_ABS_TOL_DA = 0.06;
    static constexpr double LOCATOR_REL_TOL_PPM = 10.0;
    /// How far FLASHDeconv's monoisotopic call may sit from the engine's, in isotopes.
    static constexpr int MAX_ISOTOPE_OFFSET = 2;
    /// A row's anchor precursor_mz against the spectrum's isolation target. Worst measured deviation
    /// 0.00505 over 13,873 MS2; a FOREIGN file agrees by chance for 0.09 % of scans.
    static constexpr double ANCHOR_MZ_TOL = 0.01;
    /// A follow-up MS2's parent is the MS2 that triggered it, not the survey. Depth is one by design;
    /// the bound only stops a malformed file from looping.
    static constexpr int MAX_PARENT_HOPS = 8;

    /// One scan_commands.tsv row, reduced to what a consumer of the acquisition needs.
    struct Row
    {
      std::string tracking_id;
      std::string parent_tracking_id;   ///< "" for a root command
      int ms_level = 0;
      // Everything below is filled for ms_level == 2 rows ONLY (ADR-0046 decision 8).
      double mono_mass = 0;             ///< the commanded precursor's monoisotopic mass
      int charge = 0;                   ///< the ANCHOR charge: first ',' value of a multiplexed cell
      double precursor_mz = 0;          ///< the anchor's isolation centre
      double qscore = 0, charge_cos = 0, charge_snr = 0, iso_cos = 0, snr = 0, charge_score = 0, ppm_error = 0;
      double precursor_intensity = 0, peakgroup_intensity = 0;
    };

    /// What the join needs to know about one spectrum of the data file.
    struct Scan
    {
      int scan_number = 0;
      int ms_level = 0;
      std::string scan_description;           ///< "" when the spectrum carries none
      std::vector<double> isolation_targets;  ///< one per precursor element (several under MSX)
    };

    /// The join's answer for one commanded MS2.
    struct Located
    {
      int parent_scan_number = -1;   ///< the survey the command was decided from; -1 = not in this data file
      Row row;
    };

    /// The tracking id of a scan description, or "" when it is too short to carry one.
    static std::string trackingIdOf(const std::string& scan_description)
    {
      return scan_description.size() >= 3 ? scan_description.substr(0, 3) : std::string();
    }

    /// 0 = @p candidate is the commanded mass; k = it sits k isotopes away; -1 = neither.
    /// @p residual receives |candidate - nearest accepted isotope of commanded|, the tie-break within a rank.
    static int massRank(double candidate, double commanded, double& residual)
    {
      const double tol = LOCATOR_ABS_TOL_DA + commanded * LOCATOR_REL_TOL_PPM * 1e-6;
      for (int k = 0; k <= MAX_ISOTOPE_OFFSET; ++k)
      {
        for (int sign = 1; sign >= -1; sign -= 2)
        {
          const double d = std::abs(candidate - (commanded + sign * k * Constants::ISOTOPE_MASSDIFF_55K_U));
          if (d <= tol) { residual = d; return k; }
          if (k == 0) { break; }   // +0 and -0 are the same shift
        }
      }
      return -1;
    }

    /// Read a scan_commands.tsv. Columns are resolved BY HEADER NAME, so a reorder is free and every file the
    /// engine has ever written parses. Throws FileNotFound, ParseError (a missing required column, a short or
    /// long row, a non-numeric cell) or InvalidValue (a duplicate tracking_id).
    static std::map<std::string, Row> parse(const std::string& path)
    {
      std::ifstream in(path);
      if (!in.good()) { throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path); }

      std::string line;
      if (!std::getline(in, line)) { throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, "the scan_commands file is empty"); }
      const std::vector<std::string> header = splitTabs_(line);
      std::map<std::string, size_t> col;
      for (size_t i = 0; i < header.size(); ++i) { col[header[i]] = i; }
      for (const char* name : {"tracking_id", "ms_level", "parent_tracking_id", "mono_mass", "charge", "precursor_mz", "qscore",
                               "charge_cos", "charge_snr", "iso_cos", "snr", "charge_score", "ppm_error", "precursor_intensity",
                               "peakgroup_intensity"})
      {
        if (col.find(name) == col.end())
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, String("the scan_commands file has no '") + name + "' column");
        }
      }

      std::map<std::string, Row> rows;
      size_t line_no = 1;
      while (std::getline(in, line))
      {
        ++line_no;
        const std::vector<std::string> cell = splitTabs_(line);
        if (cell.size() == 1 && cell[0].empty()) { continue; }   // blank line
        if (cell.size() != header.size())
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path,
                                      String("line ") + line_no + " has " + cell.size() + " cells, the header has " + header.size());
        }
        Row r;
        r.tracking_id = cell[col["tracking_id"]];
        r.parent_tracking_id = cell[col["parent_tracking_id"]];
        try
        {
          r.ms_level = String(cell[col["ms_level"]]).toInt();
          if (r.ms_level == 2)
          {
            r.mono_mass = String(cell[col["mono_mass"]]).toDouble();
            r.charge = String(anchorOf_(cell[col["charge"]])).toInt();
            r.precursor_mz = String(anchorOf_(cell[col["precursor_mz"]])).toDouble();
            r.qscore = String(cell[col["qscore"]]).toDouble();
            r.charge_cos = String(cell[col["charge_cos"]]).toDouble();
            r.charge_snr = String(cell[col["charge_snr"]]).toDouble();
            r.iso_cos = String(cell[col["iso_cos"]]).toDouble();
            r.snr = String(cell[col["snr"]]).toDouble();
            r.charge_score = String(cell[col["charge_score"]]).toDouble();
            r.ppm_error = String(cell[col["ppm_error"]]).toDouble();
            r.precursor_intensity = String(cell[col["precursor_intensity"]]).toDouble();
            r.peakgroup_intensity = String(cell[col["peakgroup_intensity"]]).toDouble();
          }
        }
        catch (const Exception::ConversionError& e)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, path, String("line ") + line_no + ": " + e.what());
        }
        if (!rows.emplace(r.tracking_id, r).second)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "duplicate tracking_id in the scan_commands file -- the id counter wrapped (it does at 830,583 commands) "
                                        "or two runs were concatenated; either way an id no longer names one scan",
                                        r.tracking_id);
        }
      }
      return rows;
    }

    /// Join the spectra of a data file to the rows of its acquisition. Returns the commanded MS2 scans, by
    /// scan number. Throws InvalidValue when the file does not belong to the data (see the class comment).
    static std::map<int, Located> join(const std::vector<Scan>& scans, const std::map<std::string, Row>& rows)
    {
      // Which spectrum carries which COMMANDED id. An id found in no row is an uncommanded scan -- the
      // instrument's own, which has been seen to carry a three-blank description -- and is ignored here,
      // duplicates included: thousands of placeholders may share one such non-id.
      std::map<std::string, const Scan*> scan_of;
      for (const Scan& s : scans)
      {
        const std::string tid = trackingIdOf(s.scan_description);
        if (tid.empty() || rows.find(tid) == rows.end()) { continue; }
        auto ins = scan_of.emplace(tid, &s);
        if (!ins.second)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        String("scans ") + ins.first->second->scan_number + " and " + s.scan_number + " carry the same tracking id", tid);
        }
      }
      if (scan_of.empty())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "no spectrum carries a tracking id found in the scan_commands file: either the converter dropped the "
                                      "'scan description' of each spectrum (msconvert keeps it) or the file belongs to another run",
                                      String(rows.size()) + " rows");
      }

      std::map<int, Located> located;
      String offenders;
      Size n_offenders = 0;
      for (const auto& entry : scan_of)
      {
        const Scan& s = *entry.second;
        const Row& r = rows.at(entry.first);
        bool ok = (r.ms_level == s.ms_level);
        if (ok && s.ms_level == 2)
        {
          ok = false;
          for (double target : s.isolation_targets) { if (std::abs(target - r.precursor_mz) <= ANCHOR_MZ_TOL) { ok = true; } }
        }
        if (!ok)
        {
          if (n_offenders++ < 5)
          {
            offenders += String("\n  scan ") + s.scan_number + " (tracking id " + entry.first + "): the row says MS" + r.ms_level + " at m/z "
                       + r.precursor_mz + ", the spectrum is MS" + s.ms_level + " isolating "
                       + (s.isolation_targets.empty() ? String("nothing") : String(s.isolation_targets.front()));
          }
          continue;
        }
        if (s.ms_level != 2) { continue; }
        Located l;
        l.row = r;
        l.parent_scan_number = surveyScanNumber_(r, rows, scan_of);
        located[s.scan_number] = l;
      }
      if (n_offenders > 0)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      String("this scan_commands file does not belong to this data file: ") + n_offenders + " of " + scan_of.size()
                                        + " joined scans disagree with their rows. Tracking ids restart in every run, so a foreign file joins "
                                          "every scan; first offenders:" + offenders,
                                      String(n_offenders));
      }
      return located;
    }

  private:
    static std::vector<std::string> splitTabs_(std::string line)
    {
      if (!line.empty() && line.back() == '\r') { line.pop_back(); }   // a Windows-written log read under gcc
      std::vector<std::string> out;
      size_t st = 0;
      while (true)
      {
        const size_t ed = line.find('\t', st);
        if (ed == std::string::npos) { out.push_back(line.substr(st)); break; }
        out.push_back(line.substr(st, ed - st));
        st = ed + 1;
      }
      return out;
    }

    /// "17,16,15" -> "17": the anchor comes first, its co-isolated notches after (ADR-0016).
    static std::string anchorOf_(const std::string& cell) { return cell.substr(0, cell.find(',')); }

    /// Walk parent_tracking_id up the ROWS to the MS1 survey, then ask which spectrum carries that id.
    static int surveyScanNumber_(const Row& row, const std::map<std::string, Row>& rows, const std::map<std::string, const Scan*>& scan_of)
    {
      std::string parent = row.parent_tracking_id;
      for (int hop = 0; hop < MAX_PARENT_HOPS && !parent.empty(); ++hop)
      {
        auto r = rows.find(parent);
        if (r == rows.end()) { return -1; }
        if (r->second.ms_level == 1)
        {
          auto s = scan_of.find(parent);
          return s == scan_of.end() ? -1 : s->second->scan_number;
        }
        parent = r->second.parent_tracking_id;
      }
      return -1;
    }
  };
} // namespace OpenMS
