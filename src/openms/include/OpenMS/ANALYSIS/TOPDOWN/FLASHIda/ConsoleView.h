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

#include <OpenMS/config.h>

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/PrecursorBoard.h>

#include <atomic>
#include <functional>
#include <memory>
#include <string>
#include <thread>
#include <vector>

namespace OpenMS
{
  /**
   * @brief The operator view (ADR-0065 push 2): the Precursor board, the identification timeline and the status line
   * at the bottom of the console, under a VT scroll region the plain lines scroll in; and the run summary.
   *
   * The renderers are pure functions of a BoardSnapshot and a size, and are what the tests check. The driver exists
   * only when stdout is a Windows console that accepts virtual-terminal processing (start() returns nullptr on a pipe,
   * a file, or any other platform): it redraws on its own thread, writing through EventLog::paint so a redraw never
   * interleaves with an engine line.
   */
  class OPENMS_DLLAPI ConsoleView
  {
  public:
    /// The board, the status line and the timeline: exactly @p height lines, none wider than @p width columns.
    static std::vector<std::string> renderPanel(const BoardSnapshot& snap, int width, int height, bool utf8, bool ansi);

    /// The end-of-run block (no leading spaces; each line becomes "[OP-SUMMARY] <line>" in the event log).
    static std::vector<std::string> renderSummary(const BoardSnapshot& snap, bool utf8);

    /// The identification timeline's @p cells cells: shaded by identifications per cell, relative to the busiest one,
    /// blank past the current retention time. Spans global.duration, or the retention time so far when it is 0.
    static std::string timelineCells(const BoardSnapshot& snap, int cells, bool utf8);

    /// Display columns of @p s: UTF-8 code points, ANSI escape sequences excluded.
    static int columns(const std::string& s);

    /// Start the view when stdout is a VT-capable console; nullptr otherwise. @p snapshot is called on the render
    /// thread and must not take analysis_mutex_.
    static std::unique_ptr<ConsoleView> start(std::function<BoardSnapshot()> snapshot);

    /// Print the run summary through EventLog (console, events.log, sinks), in any mode.
    static void printSummary(const BoardSnapshot& snap, bool utf8);

    /// Stop redrawing and release the scroll region; the last frame stays on screen. Idempotent.
    void stop();

    /// Whether the console takes UTF-8 glyphs (else ASCII).
    bool utf8() const { return utf8_; }

    /// Stops, then restores the console's mode and code page.
    ~ConsoleView();

    ConsoleView(const ConsoleView&) = delete;
    ConsoleView& operator=(const ConsoleView&) = delete;

  private:
    ConsoleView() = default;
    void layout_();   // query the window, push the old text up, set the scroll region
    void redraw_();

    std::function<BoardSnapshot()> snapshot_;
    std::thread thread_;
    std::atomic<bool> running_{false};
    bool utf8_ = false;
    bool stopped_ = false;
    void* handle_ = nullptr;        // HANDLE, kept opaque so this header needs no <windows.h>
    unsigned long old_mode_ = 0;
    unsigned int old_cp_ = 0;
    int width_ = 0, height_ = 0, panel_ = 0;
  };
} // namespace OpenMS
