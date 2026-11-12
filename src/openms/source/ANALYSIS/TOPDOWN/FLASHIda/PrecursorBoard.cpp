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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorBoard.h>

#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>

#include <algorithm>
#include <set>

namespace OpenMS
{
  namespace
  {
    bool hasResult(PrecursorState s)
    {
      return s == PrecursorState::NoId || s == PrecursorState::Identified || s == PrecursorState::ClosedOut;
    }

    bool hasProteoform(PrecursorState s)
    {
      return s == PrecursorState::Identified || s == PrecursorState::ClosedOut;
    }
  } // namespace

  std::string PrecursorBoard::proformaOf(const FragmentAnalysis::ProteoformMatch& match)
  {
    return FragmentAnalysis::regionProForma(match);   // the matched region, as identification.tsv's MS2 row
  }

  BoardRow& PrecursorBoard::rowFor_(double mass)
  {
    const int nominal = SpectralDeconvolution::getNominalMass(mass);
    auto it = rows_.find(nominal);
    if (it == rows_.end())
    {
      it = rows_.emplace(nominal, BoardRow{}).first;
      it->second.nominal = nominal;
      it->second.mass = mass;
    }
    it->second.last_rt_s = std::max(it->second.last_rt_s, rt_s_);
    return it->second;
  }

  void PrecursorBoard::setRt(double rt_s)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    rt_s_ = std::max(rt_s_, rt_s);
  }

  void PrecursorBoard::setLogDir(const std::string& log_dir)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    log_dir_ = log_dir;
  }

  void PrecursorBoard::onSent(const ScanCommand& cmd, int precursor_id)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    if (cmd.msn_level <= 1)
    {
      if (cmd.is_agc == 0) ++ms1_;
      return;
    }
    if (cmd.msn_level == 2)
    {
      ++ms2_;
      BoardRow& row = rowFor_(cmd.mono_mass);
      ++row.ms2;
      if (cmd.num_stages > 0) row.charge = cmd.stages[0].charge_state;
      if (!hasResult(row.state)) row.state = PrecursorState::Acquiring;
      if (precursor_id > 0) nominal_by_precursor_[precursor_id] = row.nominal;
      return;
    }
    ++ms3_;
    const auto it = nominal_by_precursor_.find(precursor_id);
    if (it == nominal_by_precursor_.end()) return;
    auto rit = rows_.find(it->second);
    if (rit == rows_.end()) return;
    ++rit->second.ms3;
    rit->second.last_rt_s = std::max(rit->second.last_rt_s, rt_s_);
  }

  void PrecursorBoard::onPlanned(int units)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    planned_ = units;
  }

  void PrecursorBoard::onReturned(double commanded_mass, bool searched)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    BoardRow& row = rowFor_(commanded_mass);
    if (row.state == PrecursorState::ClosedOut) return;
    if (searched) row.state = PrecursorState::Searching;
    else if (!hasResult(row.state)) row.state = PrecursorState::Acquired;
  }

  void PrecursorBoard::onResult(double commanded_mass, const FragmentAnalysis::ProteoformMatch& match, bool closed_out,
                                const std::optional<double>& threshold)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    BoardRow& row = rowFor_(commanded_mass);
    const bool identified = !match.proteoform_sequence.empty() && !match.fragments.empty();
    if (!identified)
    {
      if (!hasProteoform(row.state)) row.state = PrecursorState::NoId;   // an earlier proteoform stands
      return;
    }
    if (!hasProteoform(row.state)) id_rts_s_.push_back(rt_s_);   // the Precursor's first identification
    const double score = match.search_identification_score;
    if (score >= row.score || !hasProteoform(row.state))
    {
      row.score = score;
      row.accession = match.matched_protein;
      row.proforma = proformaOf(match);
    }
    if (closed_out || row.state == PrecursorState::ClosedOut)
    {
      row.state = PrecursorState::ClosedOut;
      row.below_threshold = false;
    }
    else
    {
      row.state = PrecursorState::Identified;
      row.below_threshold = threshold.has_value() && row.score < *threshold;
      row.threshold = threshold.value_or(0.0);
    }
  }

  void PrecursorBoard::onDropped(double commanded_mass)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    BoardRow& row = rowFor_(commanded_mass);
    if (!hasProteoform(row.state)) row.state = PrecursorState::NoId;
  }

  void PrecursorBoard::onSearchWait(int ms)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    search_wait_ms_.push_back(ms);
  }

  void PrecursorBoard::onSearchTime(int ms)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    search_ms_.push_back(ms);
  }

  BoardSnapshot PrecursorBoard::snapshot(double rt_window_s, double duration_min, size_t search_waiting,
                                         size_t search_capacity) const
  {
    std::lock_guard<std::mutex> lock(mutex_);
    BoardSnapshot s;
    s.rt_s = rt_s_;
    s.duration_min = duration_min;
    s.ms1 = ms1_;
    s.ms2 = ms2_;
    s.ms3 = ms3_;
    s.planned = planned_;
    s.search_waiting = search_waiting;
    s.search_capacity = search_capacity;
    s.log_dir = log_dir_;
    s.id_rts_s = id_rts_s_;
    s.search_wait_ms = search_wait_ms_;
    s.search_ms = search_ms_;
    std::set<std::string> accessions;
    for (const auto& kv : rows_)
    {
      const BoardRow& r = kv.second;
      ++s.precursors;
      if (r.state == PrecursorState::ClosedOut) ++s.closed_out;
      if (r.state == PrecursorState::Identified) ++s.identified;   // exclusive: a closed-out row is counted there
      if (r.state == PrecursorState::NoId) ++s.no_id;
      if (hasProteoform(r.state))
      {
        if (!r.accession.empty()) accessions.insert(r.accession);
        if (r.score > s.best_score)
        {
          s.best_score = r.score;
          s.best_accession = r.accession;
        }
      }
      if (r.last_rt_s >= rt_s_ - rt_window_s) s.live.push_back(r);
    }
    s.proteins = static_cast<int>(accessions.size());
    std::stable_sort(s.live.begin(), s.live.end(),
                     [](const BoardRow& a, const BoardRow& b) { return a.last_rt_s > b.last_rt_s; });
    return s;
  }
} // namespace OpenMS
