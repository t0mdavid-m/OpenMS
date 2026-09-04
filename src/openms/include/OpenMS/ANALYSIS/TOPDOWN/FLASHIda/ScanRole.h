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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommand.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>

namespace OpenMS
{

  /**
   * @brief What a scan is FOR, decoded from the 4th character of its scan_description.
   *
   * The marker has always been on the wire -- `{3-char tracking id}{marker}...` -- but it used to be
   * decoded by four hand-written `char` comparisons in two files plus a six-case switch, so "which
   * roles drive acquisition decisions?" was a question you answered by grepping. This is the single
   * decoder, and the traits table below is what makes a new role declare its own answer (ADR-0042).
   */
  enum class ScanRole : uint8_t
  {
    Survey = 0,          ///< 'S' the MS1 acquisition decisions are made from
    Agc,                 ///< 'A' the flux prescan; discarded on return, never deconvolved
    Monitor,             ///< 'M' observed and never acted upon (ADR-0042)
    Identification,      ///< 'R' the MS2/MS3 acquired to sequence a proteoform
    Quantification,      ///< 'Q' the MS2 acquired to measure reporter ions (ADR-0038)
    FollowUp,            ///< 'C' a tagging follow-up, bought by tag detection; depth one
    ExplorationVariant,  ///< 'E' one point of a CE/reaction-time sweep
    _Count
  };

  /// One row per ScanRole, indexed BY that role. See kScanRoleTraits.
  struct ScanRoleTraits
  {
    ScanRole role;
    char marker;           ///< scan_description[3]
    const char* log_name;  ///< scan_commands.tsv `scan_type` value
    bool observes;         ///< does the engine read this scan's peaks at all?
    bool decides;          ///< does processing it drive an acquisition decision?
  };

  /**
   * @brief The role table. Rows MUST be in enumerator order -- scanRoleTableIsTotal() asserts it.
   *
   * `log_name` IS NOT A NAMING OPPORTUNITY. These seven strings are the exact values
   * IdaLogger::scanTypeFromDescription_ has always returned, and they are emitted into
   * scan_commands.tsv's `scan_type` column, which is a golden column. "recording" is a poor name for
   * an identification MS2 and it MUST NOT be improved here: renaming it revalues 28 golden files for
   * no observable gain. Pinned by FLASHIda_exploration_test's
   * roleName_matches_the_legacy_scan_type_strings.
   *
   * `observes` is false for exactly one row today, and `decides` for exactly two. Neither is a rich
   * policy -- the point is that adding an eighth role cannot compile until someone has answered both
   * questions for it.
   */
  inline constexpr std::array<ScanRoleTraits, static_cast<size_t>(ScanRole::_Count)> kScanRoleTraits {{
    {ScanRole::Survey,             'S', "survey",         true,  true },
    {ScanRole::Agc,                'A', "agc",            false, false},
    {ScanRole::Monitor,            'M', "monitor",        true,  false},  // ADR-0042
    {ScanRole::Identification,     'R', "recording",      true,  true },
    {ScanRole::Quantification,     'Q', "quantification", true,  true },
    {ScanRole::FollowUp,           'C', "conditional",    true,  true },
    {ScanRole::ExplorationVariant, 'E', "exploration",    true,  true },
  }};

  /// Every enumerator has a row, at its own index, with a unique marker and a non-null name.
  constexpr bool scanRoleTableIsTotal()
  {
    for (size_t i = 0; i < kScanRoleTraits.size(); ++i)
    {
      if (static_cast<size_t>(kScanRoleTraits[i].role) != i) return false;
      if (kScanRoleTraits[i].log_name == nullptr) return false;
      if (kScanRoleTraits[i].marker == '\0') return false;
      for (size_t j = i + 1; j < kScanRoleTraits.size(); ++j)
        if (kScanRoleTraits[i].marker == kScanRoleTraits[j].marker) return false;
    }
    return true;
  }

  // Adding an enumerator without adding its row leaves the array's tail value-initialised, so `role`
  // reads Survey(0) at a non-zero index and THIS FAILS TO COMPILE. That is deliberately a hard error
  // rather than a warning: ENABLE_GCC_WERROR is OFF by default and MSVC's C4062 (unhandled enum in a
  // switch) is off-by-default, so a switch-without-default would only produce a diagnostic nobody reads.
  static_assert(scanRoleTableIsTotal(),
                "kScanRoleTraits: one row per ScanRole enumerator, in enumerator order, with unique "
                "non-null markers. Did you add an enumerator without its row?");

  constexpr const ScanRoleTraits& scanRoleTraits(ScanRole r)
  {
    return kScanRoleTraits[static_cast<size_t>(r)];
  }

  constexpr char roleMarker(ScanRole r)      { return scanRoleTraits(r).marker; }
  constexpr const char* roleName(ScanRole r) { return scanRoleTraits(r).log_name; }
  constexpr bool roleObserves(ScanRole r)    { return scanRoleTraits(r).observes; }
  constexpr bool roleDecides(ScanRole r)     { return scanRoleTraits(r).decides; }

  /**
   * @brief Decode a role from a scan description. nullopt when it is too short or the marker is unknown.
   *
   * @param len bounded length of @p d (the caller's strnlen), NOT the buffer size.
   */
  constexpr std::optional<ScanRole> roleFromDescription(const char* d, size_t len)
  {
    if (d == nullptr || len < 4) return std::nullopt;
    for (size_t i = 0; i < kScanRoleTraits.size(); ++i)
      if (kScanRoleTraits[i].marker == d[3]) return kScanRoleTraits[i].role;
    return std::nullopt;
  }

  /// Bounded strlen over scan_description, which is a fixed 256-byte buffer that a truncating
  /// snprintf can in principle leave unterminated.
  inline size_t scanDescriptionLength(const char* d, size_t cap)
  {
    if (d == nullptr) return 0;
    size_t n = 0;
    while (n < cap && d[n] != '\0') ++n;
    return n;
  }

  inline std::optional<ScanRole> roleOf(const ScanCommand& cmd)
  {
    return roleFromDescription(cmd.scan_description,
                               scanDescriptionLength(cmd.scan_description, sizeof(cmd.scan_description)));
  }

  // An UNDECODABLE description -- too short, or a marker no row claims -- resolves as a survey would:
  // observes AND decides both true. That is what keeps the migrated read sites byte-identical with
  // the `char` comparisons they replace:
  //   FLASHIda.cpp:147  `desc.size() >= 4 && desc[3] == 'A'`  ==  !roleObserves(echo_role)
  //   the MS1 fork      an unknown MS1 has always taken the survey arm
  // Do NOT "fail closed" here by returning false: that would silently stop processing every scan
  // whose description the engine did not mint, which is not this table's decision to make.
  inline bool roleObserves(std::optional<ScanRole> r) { return r ? roleObserves(*r) : true; }
  inline bool roleDecides(std::optional<ScanRole> r)  { return r ? roleDecides(*r)  : true; }

  /// scan_commands.tsv `scan_type`. "unknown" for an undecodable description -- the same string, for
  /// the same inputs, that scanTypeFromDescription_'s `default:` has always returned.
  inline const char* roleNameOf(const ScanCommand& cmd)
  {
    const std::optional<ScanRole> r = roleOf(cmd);
    return r ? roleName(*r) : "unknown";
  }

} // namespace OpenMS
