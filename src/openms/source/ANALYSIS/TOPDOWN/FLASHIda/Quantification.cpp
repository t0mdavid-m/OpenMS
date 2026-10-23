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
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>

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

#include <boost/math/distributions/complement.hpp>
#include <boost/math/distributions/students_t.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
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
    /// Two-sided Welch p for two samples of log2 values (ADR-0056): Welch-Satterthwaite degrees of
    /// freedom, the Student t distribution, exactly scipy.stats.ttest_ind(equal_var=False). -1 when
    /// either sample has fewer than two values. With no variance on either side t is undefined, and
    /// the answer is scipy's: p 0 for different means, p 1 for equal ones.
    double welchP(const std::vector<double>& a, const std::vector<double>& b)
    {
      if (a.size() < 2 || b.size() < 2) return -1.0;
      auto meanVar = [](const std::vector<double>& x, double& mean, double& var) {
        mean = std::accumulate(x.begin(), x.end(), 0.0) / static_cast<double>(x.size());
        double ss = 0.0;
        for (double v : x) ss += (v - mean) * (v - mean);
        var = ss / static_cast<double>(x.size() - 1);
      };
      double ma, va, mb, vb;
      meanVar(a, ma, va);
      meanVar(b, mb, vb);
      const double sa = va / static_cast<double>(a.size());
      const double sb = vb / static_cast<double>(b.size());
      const double se2 = sa + sb;
      if (!(se2 > 0.0)) return (ma != mb) ? 0.0 : 1.0;
      const double t  = (ma - mb) / std::sqrt(se2);
      const double nu = se2 * se2 / (sa * sa / static_cast<double>(a.size() - 1)
                                     + sb * sb / static_cast<double>(b.size() - 1));
      const boost::math::students_t dist(nu);
      return std::min(1.0, 2.0 * boost::math::cdf(boost::math::complement(dist, std::fabs(t))));
    }

    /// Median of a non-empty sample; an even count takes the mean of the two middle values.
    double medianOf(std::vector<double> v)
    {
      std::sort(v.begin(), v.end());
      const size_t n = v.size();
      return (n % 2 == 1) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
    }

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

  void Quantification::addMember(const int group_id, const ScanCommand& cmd, const int precursor_id, const int survey_id)
  {
    GroupMember m;
    m.tracking_id  = cmd.scan_id;
    m.precursor_id = precursor_id;
    m.charge       = cmd.num_stages > 0 ? std::abs(cmd.stages[0].charge_state) : 0;
    m.window_snr   = cmd.window_snr;
    m.intensity    = cmd.precursor_intensity;
    m.cmd          = cmd;
    m.isolated_charges = ScanCommandQueue::isolatedCharges(cmd);
    m.survey_id    = survey_id;
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

  Quantification::GroupDecision Quantification::decide(const int group_id) const
  {
    GroupDecision out;
    auto git = groups_.find(group_id);
    if (git == groups_.end()) return out;
    const std::vector<GroupMember>& members = git->second;
    const auto& cfg = config_.quantification();
    out.members = static_cast<int>(members.size());

    // ---- who can speak ------------------------------------------------------------------------
    // USABLE: every assigned channel measured (the pipeline's COMPLETE_ONLY). ABSENT: wholly absent
    // in one condition against a full other, read as Differential by evaluate() only under
    // empty_condition "differential" -- under "incomplete" such a member is IncompleteChannels and
    // so is neither. Everything else that returned is neither, and is never averaged in.
    std::vector<const GroupMember*> received, usable, absent;
    for (const GroupMember& m : members)
    {
      if (!m.received) continue;
      received.push_back(&m);
      if (m.result.complete) usable.push_back(&m);
      else if (m.result.verdict == Result::Verdict::Differential) absent.push_back(&m);
    }
    out.usable = static_cast<int>(usable.size());
    for (const GroupMember* m : usable) out.usable_charges.push_back(m->charge);

    // The charges the group's Q scans covered -- every one that returned, usable or not, since each isolated its
    // charges whatever it measured: the only ones the bought R may move to (survey re-planning).
    for (const GroupMember* m : received)
      out.q_covered_charges.insert(out.q_covered_charges.end(), m->isolated_charges.begin(), m->isolated_charges.end());
    std::sort(out.q_covered_charges.begin(), out.q_covered_charges.end());
    out.q_covered_charges.erase(std::unique(out.q_covered_charges.begin(), out.q_covered_charges.end()), out.q_covered_charges.end());

    // The identification window: the purest (window_snr), then the brightest, then the lower
    // charge -- total and reproducible, since window_snr saturates at 1000.0.
    auto better = [](const GroupMember* a, const GroupMember* b) {
      if (a->window_snr != b->window_snr) return a->window_snr > b->window_snr;
      if (a->intensity  != b->intensity)  return a->intensity  > b->intensity;
      return a->charge < b->charge;
    };
    auto pick = [&](const std::vector<const GroupMember*>& set) -> const GroupMember* {
      const GroupMember* best = nullptr;
      for (const GroupMember* m : set) if (best == nullptr || better(m, best)) best = m;
      return best;
    };

    const size_t floor = static_cast<size_t>(std::max(1, cfg.min_usable_scans));
    const GroupMember* winner = nullptr;

    if (!usable.empty() && usable.size() >= floor)
    {
      // ---- the aggregate (ADR-0056 decision 5) --------------------------------------------------
      // Each member scaled to its share of its own ASSIGNED reporter signal -- the channels that
      // enter the test, as the pipeline's SCAN_NORMALIZE divides by its QUANT_COLS -- so a bright
      // charge state does not dominate; then one value per channel by the authored aggregate.
      std::vector<size_t> assigned;
      for (const auto& cond : cfg.conditions) for (size_t idx : cond.channels) assigned.push_back(idx);
      const size_t n = usable.front()->result.channels.size();
      std::vector<std::vector<double>> per_channel(n);
      for (const GroupMember* m : usable)
      {
        const auto& ch = m->result.channels;
        double scale = 1.0;
        if (cfg.scan_normalize)
        {
          double total = 0.0;
          for (size_t idx : assigned) if (idx < ch.size()) total += ch[idx];
          if (total > 0.0) scale = 1.0 / total;
        }
        for (size_t i = 0; i < n && i < ch.size(); ++i) per_channel[i].push_back(ch[i] * scale);
      }
      std::vector<double> aggregate(n, 0.0);
      for (size_t i = 0; i < n; ++i)
      {
        const auto& v = per_channel[i];
        if (v.empty()) continue;
        const double sum = std::accumulate(v.begin(), v.end(), 0.0);
        switch (cfg.aggregate)
        {
          case QuantAggregate::Median: aggregate[i] = medianOf(v); break;
          case QuantAggregate::Mean:   aggregate[i] = sum / static_cast<double>(v.size()); break;
          case QuantAggregate::Sum:    aggregate[i] = sum; break;
        }
      }
      // Floor 0: a share can legitimately sit below 1e-3, and usability already guaranteed every
      // assigned channel was measured.
      out.result = evaluate(aggregate, cfg, 0.0);
      winner = pick(usable);
    }
    else if (usable.empty() && !absent.empty() && absent.size() >= floor)
    {
      // ---- the empty condition decides (ADR-0056 decisions 4 and 6) ---------------------------
      // Only when no member was complete, and -- owner, 2026-10-15 -- only when at least the floor's
      // worth of members say so. Its condition_means carry the direction enriched_in reads.
      winner = pick(absent);
      out.result = winner->result;
    }
    else
    {
      // ---- undecided -------------------------------------------------------------------------
      // A verdict about the MEASUREMENT, not the species. The window is still named (ADR-0056
      // decision 9) so `identify: all` buys it; quantified/differential refuse it on the verdict.
      out.result.verdict = received.empty() ? Result::Verdict::ExtractionFailed
                                            : Result::Verdict::IncompleteChannels;
      winner = pick(received);
    }

    if (winner != nullptr)
    {
      out.has_winner          = true;
      out.winner_charge       = winner->charge;
      out.winner_precursor_id = winner->precursor_id;
      out.winner_cmd          = winner->cmd;
      out.winner_survey_id    = winner->survey_id;
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

    return evaluate(out.channels, cfg, kEmptyChannel);
  }

  Quantification::Result Quantification::evaluate(const std::vector<double>& channels,
                                                   const QuantConfig& cfg, const double empty_floor)
  {
    Result out;
    out.channels = channels;
    const size_t n = channels.size();

    // ---- condition by condition, over ASSIGNED channels only ----------------------------------
    // A channel named in no condition is unassigned: ignored here, still reported in `channels`.
    // That is what makes a bridge channel, or a kit run below capacity, legal.
    if (cfg.conditions.size() != 2) return out;  // Config rejects this at load

    std::array<std::vector<double>, 2> log2v;
    std::array<bool, 2> wholly_absent {{false, false}};
    bool any_assigned_empty = false;
    for (size_t c = 0; c < 2; ++c)
    {
      const auto& chans = cfg.conditions[c].channels;
      if (chans.empty()) return out;  // Config rejects this at load

      size_t empties = 0;
      for (size_t idx : chans)
      {
        const double v = (idx < n) ? channels[idx] : 0.0;
        if (!(v > 0.0) || v < empty_floor) { ++empties; continue; }
        // ADR-0057: emptiness above is judged on the RAW value; the authored channel
        // normalisation enters only what the test reads. `out.channels` stays raw -- it is what
        // scan_results.tsv logs, and the offline pipeline normalises those columns itself.
        const double f = cfg.channel_normalization.empty() ? 1.0 : cfg.channel_normalization[idx];
        log2v[c].push_back(std::log2(v) - std::log2(f));
      }
      if (empties == chans.size()) { wholly_absent[c] = true; }
      else if (empties > 0)        { any_assigned_empty = true; }
      // GEOMETRIC mean (ADR-0056), defined only for a fully measured condition; 0 otherwise -- a
      // geometric mean of a 0 is undefined, and 0 is what an absent condition has always logged.
      if (empties == 0)
        out.condition_means[c] = std::exp2(std::accumulate(log2v[c].begin(), log2v[c].end(), 0.0)
                                           / static_cast<double>(log2v[c].size()));
    }

    // Both absent means no signal at all.
    if (wholly_absent[0] && wholly_absent[1])
    {
      out.verdict = Result::Verdict::IncompleteChannels;
      return out;
    }
    // A PARTIALLY empty condition is untrustworthy rather than informative, and is checked BEFORE
    // the empty-condition rule (ADR-0038, amended 2026-10-15): absence is read only against a fully
    // measured condition.
    if (any_assigned_empty)
    {
      out.verdict = Result::Verdict::IncompleteChannels;
      return out;
    }
    // A WHOLLY absent condition against one measured in every channel: what it means is authored
    // (ADR-0056 decision 4), and it is read BEFORE any test -- a test has nothing to compute on an
    // empty side. Under "differential" the species is present on one side and not the other, the
    // strongest result the experiment can produce; fold_change and p stay -1, condition_means carry
    // the direction.
    if (wholly_absent[0] || wholly_absent[1])
    {
      out.verdict = (cfg.empty_condition == EmptyCondition::Differential)
                        ? Result::Verdict::Differential
                        : Result::Verdict::IncompleteChannels;
      return out;
    }

    // ---- complete: both numbers, under every test, then the authored test reads one or both ----
    out.complete    = true;
    out.fold_change = out.condition_means[0] / out.condition_means[1];
    out.p_value     = welchP(log2v[0], log2v[1]);

    const double t = cfg.fold_change_threshold;
    const bool ratio_hit = out.fold_change > t || (1.0 / out.fold_change) > t;
    const bool p_hit     = out.p_value >= 0.0 && out.p_value <= cfg.p_value_threshold;
    bool hit = ratio_hit;
    switch (cfg.test)
    {
      case QuantTest::FoldChange:         hit = ratio_hit; break;
      case QuantTest::Welch:              hit = p_hit; break;
      case QuantTest::WelchAndFoldChange: hit = ratio_hit && p_hit; break;
    }
    out.verdict = hit ? Result::Verdict::Differential : Result::Verdict::NotDifferential;
    return out;
  }

} // namespace OpenMS
