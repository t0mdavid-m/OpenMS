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

#include <functional>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

namespace OpenMS
{
  /// Severity of an engine event. Warn and Error also reach the console (ADR-0065 d.1).
  enum class EventLevel
  {
    Info,
    Warn,
    Error
  };

  /// The operator line's kind: its column on the console, and its `[OP-…]` tag in the event log.
  enum class OperatorKind
  {
    Ms1,   ///< one decision survey
    Id,    ///< one live identification
    Ms3,   ///< one MS3 batch
    Info   ///< startup and status prose
  };

  /**
   * @brief The engine's one writer of text output (ADR-0065): the event log and the operator lines.
   *
   * Process-wide, as std::cout was before it: the emitters live in Config, ScanCommandQueue, Exploration,
   * ProteoformTracker, FragmentAnalysis, MS3FragmentMatcher, PrecursorSelection, SearchPool and FLASHIda, and many
   * of them are built standalone in tests. Every line is composed whole and written in one call under one mutex,
   * so lines from the analysis thread and the search pool never tear.
   *
   * Destinations:
   * - **events.log** (`<runtime.log_dir>/events.log`, opened by IdaLogger): every event and every operator line,
   *   each prefixed `HH:mm:ss.fff  rt <minutes>  LEVEL  `. Not a golden stream. Lines emitted before the file
   *   opens (Config warns while it parses) are kept, at most kPendingCap of them, and written when it opens;
   *   openFile("") (no log_dir) discards them.
   * - **Sinks** (tests): the bare line, without the prefix, so `[TAG] key=value` assertions read as before.
   * - **The console**: operator lines and Warn/Error events only, as `HH:mm:ss.fff  KIND   text`. An Info event
   *   never reaches it.
   *
   * Hazard: two engines alive at once share one file; the second openFile replaces the first.
   */
  class OPENMS_DLLAPI EventLog
  {
  public:
    using Sink = std::function<void(const std::string& line)>;

    /// Lines kept for a file that is not open yet; the oldest are dropped past this.
    static constexpr size_t kPendingCap = 1000;

    /// An event: `line` is "[TAG] key=value …". File and sinks; the console too when level is Warn or Error.
    static void record(EventLevel level, const std::string& line);

    /// An operator line: console, file (as "[OP-<KIND>] text") and sinks.
    static void say(OperatorKind kind, const std::string& text);

    /// The run summary (ADR-0065 push 2): each line to the console as it is, and to the file and the sinks as
    /// "[OP-SUMMARY] <line>".
    static void summary(const std::vector<std::string>& lines);

    /// Raw bytes to the console under the same lock as every line (the operator view's redraw: escape sequences
    /// included). Not recorded in the file or the sinks.
    static void paint(const std::string& raw);

    /// Colour the console's KIND column (the operator view sets it while it runs on a VT console).
    static void setAnsi(bool ansi);

    /// The last retention time seen, in seconds (the engine's unit); the file prefix prints it in minutes.
    static void setRt(double rt_seconds);

    /// Open @p path in append mode (closing any open file) and write the pending lines; resets the RT. An empty
    /// path discards the pending lines when no file is open, and leaves an open file alone.
    static void openFile(const std::string& path);

    /// Close the file (IdaLogger's destructor) and forget the RT. Later lines are pending again.
    static void closeFile();

    /// Register a sink; returns its handle for removeSink.
    static int addSink(Sink sink);
    static void removeSink(int handle);

    /// Where operator lines go: std::cout by default; nullptr silences the console.
    static void setConsole(std::ostream* console);

    /// The tag operator lines carry in the file and the sinks, e.g. "[OP-ID]".
    static const char* operatorTag(OperatorKind kind);
  };

  /**
   * @brief One event line, composed with `<<` and emitted whole when it goes out of scope.
   *
   * `EventLine() << "[REPLAN] cv=" << cv << " kept=" << k;` replaces `std::cout << … << std::endl;`. Each line
   * formats on its own stream, so no manipulator (std::fixed, std::setprecision) leaks into the next line the way
   * it did on the shared std::cout. A named one (`EventLine ev; ev << …; if (x) ev << …;`) emits at its scope's end.
   */
  class EventLine
  {
  public:
    explicit EventLine(EventLevel level = EventLevel::Info) : level_(level) {}
    EventLine(const EventLine&) = delete;
    EventLine& operator=(const EventLine&) = delete;
    ~EventLine()
    {
      try { EventLog::record(level_, os_.str()); }
      catch (...) {}   // a log line is never worth terminating the acquisition for
    }
    template<typename T> EventLine& operator<<(const T& value)
    {
      os_ << value;
      return *this;
    }
    /// Manipulators (std::fixed, std::showpos); std::endl and std::flush are no-ops, the line ends at scope exit.
    EventLine& operator<<(std::ostream& (*manip)(std::ostream&))
    {
      if (manip != static_cast<std::ostream& (*)(std::ostream&)>(std::endl<char, std::char_traits<char>>)
          && manip != static_cast<std::ostream& (*)(std::ostream&)>(std::flush<char, std::char_traits<char>>))
        manip(os_);
      return *this;
    }
    EventLine& operator<<(std::ios_base& (*manip)(std::ios_base&))
    {
      manip(os_);
      return *this;
    }

  private:
    EventLevel level_;
    std::ostringstream os_;
  };
} // namespace OpenMS
