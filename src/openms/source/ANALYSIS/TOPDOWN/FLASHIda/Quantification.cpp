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
