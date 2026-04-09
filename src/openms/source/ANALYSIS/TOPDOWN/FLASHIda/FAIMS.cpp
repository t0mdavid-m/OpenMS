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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/FAIMS.h>

#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  FAIMS::FAIMS(const Config& config)
  {
    const auto& fc = config.faims();
    cv_values_ = fc.cv_values;
    max_cv_skip_ = fc.max_cv_skip;
    precursor_threshold_ = fc.precursor_threshold;
    enabled_ = fc.enabled;
    if (enabled_)
    {
      int n = static_cast<int>(cv_values_.size());
      cv_skip_amount_.resize(n, 0);
      cv_skip_count_.resize(n, 0);
      current_cv_index_ = 0;
    }
  }

  bool FAIMS::isEnabled() const
  {
    return enabled_;
  }

  double FAIMS::currentCV() const
  {
    return cv_values_[current_cv_index_];
  }

  int FAIMS::cvSkipAmount(size_t index) const
  {
    return cv_skip_amount_[index];
  }

  void FAIMS::updateSkip(double cv, int precursor_count)
  {
    if (!enabled_) return;

    // Find position for this CV value
    int pos = -1;
    for (int i = 0; i < static_cast<int>(cv_values_.size()); ++i)
    {
      if (cv_values_[i] == cv) { pos = i; break; }
    }
    if (pos < 0) return;  // unknown CV

    if (precursor_count < precursor_threshold_)  // strictly < (audit Q2)
    {
      if (cv_skip_amount_[pos] < max_cv_skip_)
      {
        cv_skip_amount_[pos] *= 2;                  // double spacing (audit Q3)
        if (cv_skip_amount_[pos] <= 0)
          cv_skip_amount_[pos] = 1;                 // min = 1
        if (cv_skip_amount_[pos] > max_cv_skip_)
          cv_skip_amount_[pos] = max_cv_skip_;      // cap at max
      }
      cv_skip_count_[pos] = 0;                      // reset in BOTH branches (audit Q3)
    }
    else
    {
      cv_skip_amount_[pos] = 0;                     // high precursor count: reset
      cv_skip_count_[pos] = 0;
    }
  }

  double FAIMS::advanceToNextCV()
  {
    int n = static_cast<int>(cv_values_.size());
    // Safety bound (C# uses while(true); we bound to n iterations)
    for (int attempts = 0; attempts < n; ++attempts)
    {
      current_cv_index_++;                           // increment-first (audit Q7)
      if (current_cv_index_ >= n)
        current_cv_index_ = 0;                       // wrap (audit Q7)

      if (cv_skip_count_[current_cv_index_] < cv_skip_amount_[current_cv_index_])
      {
        cv_skip_count_[current_cv_index_]++;         // skip this CV (audit Q3)
        OPENMS_LOG_DEBUG << "[FAIMS] Skipping CV=" << cv_values_[current_cv_index_]
                         << " (" << cv_skip_count_[current_cv_index_]
                         << "/" << cv_skip_amount_[current_cv_index_] << ")" << std::endl;
      }
      else
      {
        OPENMS_LOG_DEBUG << "[FAIMS] Changed to CV=" << cv_values_[current_cv_index_] << std::endl;
        return cv_values_[current_cv_index_];  // use this CV
      }
    }
    // Fallback: all CVs being skipped — use current anyway
    return cv_values_[current_cv_index_];
  }

} // namespace OpenMS
