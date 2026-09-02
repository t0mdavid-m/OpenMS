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
// $Authors: Kyowon Jeong, Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Quantification.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTEighteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTElevenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <sstream>
#include <memory>
#include <numeric>
#include <vector>

namespace OpenMS
{
  namespace
  {
    /// An assigned channel at or below this intensity counts as absent. Same 1e-3 floor the
    /// pre-ADR-0038 completeness gate used; kept so a config that merely gains `conditions`
    /// classifies its spectra the same way.
    constexpr double kEmptyChannel = 1e-3;

    /// Build the OpenMS method named by @p labelling. nullptr for an unknown name -- Config
    /// validates the name at load, so a nullptr here means the two lists drifted.
    std::unique_ptr<IsobaricQuantitationMethod> makeMethod(const std::string& labelling)
    {
      if (labelling == "tmt6plex")   return std::make_unique<TMTSixPlexQuantitationMethod>();
      if (labelling == "tmt10plex")  return std::make_unique<TMTTenPlexQuantitationMethod>();
      if (labelling == "tmt11plex")  return std::make_unique<TMTElevenPlexQuantitationMethod>();
      if (labelling == "tmt16plex")  return std::make_unique<TMTSixteenPlexQuantitationMethod>();
      if (labelling == "tmt18plex")  return std::make_unique<TMTEighteenPlexQuantitationMethod>();
      if (labelling == "itraq4plex") return std::make_unique<ItraqFourPlexQuantitationMethod>();
      if (labelling == "itraq8plex") return std::make_unique<ItraqEightPlexQuantitationMethod>();
      return nullptr;
    }

    /// True when every entry of an authored correction matrix records NO cross-talk, i.e. the
    /// matrix builds to the identity and correction would be a no-op.
    ///
    /// This has to be detected HERE rather than delegated. IsobaricIsotopeCorrector THROWS on an
    /// identity matrix ("...leading to no correction. Please provide a valid isotope_correction
    /// matrix as it was provided with the sample kit!"), so passing one through would abort the
    /// measurement of every scan rather than skip the correction. An all-zero matrix is the
    /// documented way to turn correction off in this schema, so the engine honours that by not
    /// running the corrector at all -- which is also why FLASHDeconv exposes correction as a
    /// boolean instead.
    bool isNoCorrection(const std::vector<std::string>& matrix)
    {
      for (const auto& row : matrix)
      {
        std::string tok;
        std::istringstream rs(row);
        while (std::getline(rs, tok, '/'))
        {
          // The same three spellings IsobaricQuantitationMethod treats as "no contribution".
          std::string t;
          for (char c : tok) if (!std::isspace(static_cast<unsigned char>(c))) t += std::toupper(static_cast<unsigned char>(c));
          if (t.empty() || t == "NA" || t == "-1") continue;
          try { if (std::stod(t) != 0.0) return false; }
          catch (...) { return false; }   // unparseable -> let OpenMS produce the real error
        }
      }
      return true;
    }
  } // namespace

  Quantification::Quantification(const Config& config) :
      config_(config)
  {
  }

  const std::vector<std::string>& Quantification::labellingNames()
  {
    // TopDownIsobaricQuantification.cpp's setValidStrings("type", ...) MINUS "none": that tool has
    // no separate on/off key, this section does (quantification.enabled), so accepting "none" here
    // would make `enabled: true, labelling: "none"` a contradiction someone has to resolve.
    static const std::vector<std::string> names {
      "itraq4plex", "itraq8plex", "tmt6plex", "tmt10plex", "tmt11plex", "tmt16plex", "tmt18plex"};
    return names;
  }

  std::vector<std::string> Quantification::channelNames(const std::string& labelling)
  {
    auto method = makeMethod(labelling);
    if (method == nullptr) return {};
    std::vector<std::string> names;
    for (const auto& ch : method->getChannelInformation())
      names.push_back(ch.name);
    return names;
  }

  const char* Quantification::verdictName(Result::Verdict v)
  {
    switch (v)
    {
      case Result::Verdict::Differential:       return "differential";
      case Result::Verdict::NotDifferential:    return "not_differential";
      case Result::Verdict::IncompleteChannels: return "incomplete_channels";
      case Result::Verdict::ExtractionFailed:   return "extraction_failed";
    }
    return "extraction_failed";
  }

  // ---- the quantification group and its vote (ADR-0040) ----------------------------------------

  void Quantification::addMember(const int group_id, const ScanCommand& cmd, const int precursor_id)
  {
    GroupMember m;
    m.tracking_id  = cmd.scan_id;
    m.precursor_id = precursor_id;
    m.charge       = cmd.num_stages > 0 ? std::abs(cmd.stages[0].charge_state) : 0;
    m.window_snr   = cmd.window_snr;
    m.intensity    = cmd.precursor_intensity;
    m.cmd          = cmd;
    groups_[group_id].push_back(m);
    group_by_tracking_[cmd.scan_id] = group_id;
  }

  int Quantification::deposit(const int tracking_id, const Result& result)
  {
    auto it = group_by_tracking_.find(tracking_id);
    if (it == group_by_tracking_.end()) return 0;
    auto git = groups_.find(it->second);
    if (git == groups_.end()) return 0;
    for (GroupMember& m : git->second)
    {
      if (m.tracking_id != tracking_id) continue;
      m.received = true;
      m.result   = result;
      return it->second;
    }
    return 0;
  }

  bool Quantification::isComplete(const int group_id) const
  {
    auto it = groups_.find(group_id);
    if (it == groups_.end()) return false;
    // Exploration.cpp:643's all_of(variants, received), deliberately -- including its lack of a
    // timeout. A 'Q' the instrument never returns leaves the group open forever, exactly as a lost
    // exploration variant does today (ADR-0040 decision 8).
    return std::all_of(it->second.begin(), it->second.end(),
                       [](const GroupMember& m) { return m.received; });
  }

  void Quantification::closeGroup(const int group_id)
  {
    auto it = groups_.find(group_id);
    if (it == groups_.end()) return;
    for (const GroupMember& m : it->second) group_by_tracking_.erase(m.tracking_id);
    groups_.erase(it);
  }

  Quantification::Consensus Quantification::decide(const int group_id) const
  {
    Consensus out;
    auto git = groups_.find(group_id);
    if (git == groups_.end()) return out;
    const std::vector<GroupMember>& members = git->second;

    // ---- ballots -------------------------------------------------------------------------------
    // A member that measured nothing usable ABSTAINS; it is evidence of nothing, not evidence of
    // disagreement, and counting it as dissent would manufacture chimericity out of a weak charge.
    std::vector<const GroupMember*> ballots;
    bool any_returned = false;
    for (const GroupMember& m : members)
    {
      if (m.received) any_returned = true;
      if (!m.received) continue;
      if (m.result.verdict == Result::Verdict::Differential
          || m.result.verdict == Result::Verdict::NotDifferential)
        ballots.push_back(&m);
    }

    const int min_ballots = std::max(1, config_.targeting().min_charge_states);
    if (static_cast<int>(ballots.size()) < min_ballots)
    {
      // Not a verdict about the SPECIES but about the MEASUREMENT, and it reads as one. Folding it
      // into IncompleteChannels follows ADR-0039's own reasoning that the two failure values "differ
      // only in why" -- and each member's own why is already on its own scan_results row.
      out.result.verdict = any_returned ? Result::Verdict::IncompleteChannels
                                        : Result::Verdict::ExtractionFailed;
      for (const GroupMember* b : ballots) { out.ballot_charges.push_back(b->charge);
                                             out.ballot_agreed.push_back(false); }
      out.total_ballots = static_cast<int>(ballots.size());
      return out;
    }

    // ---- the DIRECTIONAL key -------------------------------------------------------------------
    // Voting on the verdict ALONE is insufficient and fails silently: fold changes 2.50 / 2.50 /
    // 0.42 are all Differential and vote 3-0 unanimous while disagreeing about WHICH condition is
    // enriched -- precisely the interference this feature exists to detect.
    //
    // The direction is read from condition_means, NEVER from fold_change: a wholly-absent condition
    // is Differential with fold_change == -1.0, a SENTINEL, so any fold_change-based test silently
    // drops the strongest results in the experiment (ADR-0039).
    auto ballotKey = [](const GroupMember* m) -> std::pair<int, int> {
      const int verdict = static_cast<int>(m->result.verdict);
      if (m->result.verdict != Result::Verdict::Differential) return {verdict, -1};
      const int enriched = m->result.condition_means[0] > m->result.condition_means[1] ? 0 : 1;
      return {verdict, enriched};
    };

    std::map<std::pair<int, int>, std::vector<const GroupMember*>> blocs;
    for (const GroupMember* b : ballots) blocs[ballotKey(b)].push_back(b);

    // ---- pick the winning bloc -----------------------------------------------------------------
    // Majority; a tie goes to the bloc holding the highest-window_snr ballot. Breaking it on
    // window_snr rather than intensity is what makes the tie-breaking ballot the ID charge BY
    // CONSTRUCTION -- the same member wins both comparisons, so the two can never name different
    // charges (ADR-0040 decision 5).
    auto better = [](const GroupMember* a, const GroupMember* b) {
      if (a->window_snr != b->window_snr) return a->window_snr > b->window_snr;
      if (a->intensity  != b->intensity)  return a->intensity  > b->intensity;
      return a->charge < b->charge;    // total and reproducible: window_snr saturates at 1000.0
    };

    const std::vector<const GroupMember*>* winner_bloc = nullptr;
    const GroupMember* winner = nullptr;
    for (const auto& kv : blocs)
    {
      const std::vector<const GroupMember*>& bloc = kv.second;
      const GroupMember* best = bloc.front();
      for (const GroupMember* m : bloc) if (better(m, best)) best = m;

      if (winner_bloc == nullptr
          || bloc.size() > winner_bloc->size()
          || (bloc.size() == winner_bloc->size() && better(best, winner)))
      {
        winner_bloc = &bloc;
        winner      = best;
      }
    }

    // ---- the aggregate, over the AGREEING members only ------------------------------------------
    // Intensity-WEIGHTED, by summing condition_means and ratioing the sums rather than averaging
    // the per-charge ratios: summing IS the weighting, and ratios of small numbers are noisy. The
    // vote above is deliberately UNWEIGHTED -- a weak charge's opinion is still evidence of
    // interference. SNR admits, intensity ranks; the same split selectNotches already makes.
    std::array<double, 2> sums {{0.0, 0.0}};
    std::vector<double> channels;
    for (const GroupMember* m : *winner_bloc)
    {
      sums[0] += m->result.condition_means[0];
      sums[1] += m->result.condition_means[1];
      if (channels.size() < m->result.channels.size()) channels.resize(m->result.channels.size(), 0.0);
      for (size_t i = 0; i < m->result.channels.size(); ++i) channels[i] += m->result.channels[i];
    }

    out.result.verdict         = winner_bloc->front()->result.verdict;
    out.result.condition_means = sums;
    out.result.channels        = channels;
    // Inherit ADR-0038's sentinel rules: no finite ratio when a condition is wholly absent.
    out.result.fold_change = (sums[0] > 0.0 && sums[1] > 0.0) ? sums[0] / sums[1] : -1.0;

    out.has_winner         = true;
    out.winner_charge      = winner->charge;
    out.winner_precursor_id= winner->precursor_id;
    out.winner_cmd         = winner->cmd;
    out.agreeing           = static_cast<int>(winner_bloc->size());
    out.total_ballots      = static_cast<int>(ballots.size());

    for (const GroupMember* b : ballots)
    {
      out.ballot_charges.push_back(b->charge);
      out.ballot_agreed.push_back(
          std::find(winner_bloc->begin(), winner_bloc->end(), b) != winner_bloc->end());
    }
    return out;
  }

  Quantification::Result Quantification::measure(const double* mzs, const double* ints,
                                                 const int length, const double rt,
                                                 const int ms_level, const char* name) const
  {
    Result out;
    const auto& cfg = config_.quantification();

    auto method = makeMethod(cfg.labelling);
    if (method == nullptr) return out;  // ExtractionFailed; Config rejects this at load

    // Optional authored correction matrix. Empty leaves the method's stock matrix in place, which
    // is the manufacturer's typical lot rather than this kit's -- an override exists for anyone
    // holding their data sheet, and an all-zero matrix means "no correction" (handled below).
    const bool correction_off = !cfg.correction_matrix.empty() && isNoCorrection(cfg.correction_matrix);
    if (!cfg.correction_matrix.empty() && !correction_off)
    {
      Param mp = method->getParameters();
      mp.setValue("correction_matrix", cfg.correction_matrix);
      method->setParameters(mp);  // throws on a wrong entry/column count -- OpenMS validates it for us
    }

    // ---- the spectrum, as the extractor wants it ----------------------------------------------
    MSSpectrum spec;
    for (int i = 0; i < length; i++)
    {
      if (ints[i] > 0) { spec.emplace_back(mzs[i], ints[i]); }
    }
    spec.setMSLevel(ms_level);
    spec.setName(name);
    spec.setRT(rt);

    // The precursor is fabricated with HCD purely to defeat IsobaricChannelExtractor's own
    // `select_activation: auto` filter -- it is not a claim about how this scan was acquired. Note
    // the consequence: the extractor never rejects on activation, so a scan whose activation cannot
    // release the reporter yields empty channels rather than ExtractionFailed. That surfaces as
    // IncompleteChannels, which is the honest report and is why the verdict is a four-way enum.
    Precursor precursor;
    precursor.setActivationMethods({Precursor::ActivationMethod::HCD});
    spec.setPrecursors({precursor});

    MSExperiment exp;
    exp.addSpectrum(spec);

    // ---- extract, then correct ----------------------------------------------------------------
    IsobaricChannelExtractor extractor(method.get());
    Param ep = extractor.getParameters();
    ep.setValue("reporter_mass_shift", cfg.reporter_mz_tol);
    extractor.setParameters(ep);

    ConsensusMap raw;
    extractor.extractChannels(exp, raw);

    // IsobaricQuantifier's stock parameters correct isotope impurity and do NOTHING else:
    // `isotope_correction` defaults true, `normalization` defaults false. Leaving normalization
    // alone is load-bearing -- median-of-ratios against a reference channel would silently rescale
    // every fold change.
    IsobaricQuantifier quantifier(method.get());
    if (correction_off)
    {
      Param qp = quantifier.getParameters();
      qp.setValue("isotope_correction", "false");
      quantifier.setParameters(qp);
    }
    ConsensusMap corrected;
    quantifier.quantify(raw, corrected);

    // ---- read intensities BY CHANNEL ORDINAL --------------------------------------------------
    // IsobaricChannelExtractor inserts each channel with map_index == its ordinal in
    // getChannelInformation(). That identity is stated, so read it; the pre-ADR-0038 code discarded
    // it and re-derived the ordering with an m/z sort, which only agreed because every method class
    // happens to declare its channels in ascending m/z.
    const size_t n = method->getNumberOfChannels();
    out.channels.assign(n, 0.0);
    size_t seen = 0;
    for (const auto& cf : corrected)
    {
      for (auto it = cf.begin(); it != cf.end(); ++it)
      {
        const size_t idx = static_cast<size_t>(it->getMapIndex());
        if (idx < n) { out.channels[idx] = it->getIntensity(); ++seen; }
      }
    }
    if (seen != n) return out;  // ExtractionFailed: no usable channel set for this spectrum

    // ---- condition means, over ASSIGNED channels only -----------------------------------------
    // A channel named in no condition is unassigned: ignored here, still reported in `channels`.
    // That is what makes a bridge channel, or a kit run below capacity, legal -- the old gate read
    // every channel in the scheme, so four samples in six-plex chemistry rejected every spectrum.
    if (cfg.conditions.size() != 2) return out;  // Config rejects this at load

    bool any_assigned_empty = false;
    std::array<bool, 2> wholly_absent {{false, false}};
    for (size_t c = 0; c < 2; ++c)
    {
      const auto& chans = cfg.conditions[c].channels;
      if (chans.empty()) return out;  // Config rejects this at load

      double sum = 0.0;
      size_t empties = 0;
      for (size_t idx : chans)
      {
        const double v = (idx < n) ? out.channels[idx] : 0.0;
        if (v < kEmptyChannel) ++empties;
        sum += v;
      }
      out.condition_means[c] = sum / static_cast<double>(chans.size());
      if (empties == chans.size()) { wholly_absent[c] = true; }
      else if (empties > 0)        { any_assigned_empty = true; }
    }

    // A condition that is WHOLLY absent is not missing data, it is the measurement: the species is
    // present in one condition and not the other, which is the strongest result the experiment can
    // produce. There is no finite ratio, so fold_change stays -1 and condition_means carries the
    // truth. Both absent means no signal at all.
    if (wholly_absent[0] && wholly_absent[1])
    {
      out.verdict = Result::Verdict::IncompleteChannels;
      return out;
    }
    if (wholly_absent[0] || wholly_absent[1])
    {
      out.verdict = Result::Verdict::Differential;
      return out;
    }
    // A PARTIALLY empty condition is untrustworthy rather than informative: the zero is folded into
    // the mean and biases the ratio, and there is no honest number to report.
    if (any_assigned_empty)
    {
      out.verdict = Result::Verdict::IncompleteChannels;
      return out;
    }

    out.fold_change = out.condition_means[0] / out.condition_means[1];
    const double t = cfg.fold_change_threshold;
    out.verdict = (out.fold_change > t || (1.0 / out.fold_change) > t)
                      ? Result::Verdict::Differential
                      : Result::Verdict::NotDifferential;
    return out;
  }

} // namespace OpenMS
