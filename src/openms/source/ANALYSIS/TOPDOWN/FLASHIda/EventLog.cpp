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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/EventLog.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <mutex>

namespace OpenMS
{
  namespace
  {
    struct State
    {
      std::mutex mutex;
      std::ofstream file;
      std::deque<std::string> pending;   // prefixed file lines waiting for openFile
      std::map<int, EventLog::Sink> sinks;
      int next_handle = 1;
      std::ostream* console = &std::cout;
      double rt_seconds = std::nan("");
      bool ansi = false;   // the console takes SGR colour (the operator view is running)
    };

    State& state()
    {
      // Deliberately leaked: an engine destroyed during static teardown still finds its log alive.
      static State* s = new State;
      return *s;
    }

    /// "HH:mm:ss.fff", local wall-clock time.
    std::string clock()
    {
      const auto now = std::chrono::system_clock::now();
      const std::time_t t = std::chrono::system_clock::to_time_t(now);
      const long ms = static_cast<long>(std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count() % 1000);
      std::tm tm{};
#ifdef _WIN32
      localtime_s(&tm, &t);
#else
      localtime_r(&t, &tm);
#endif
      char buf[16];
      std::snprintf(buf, sizeof(buf), "%02d:%02d:%02d.%03ld", tm.tm_hour, tm.tm_min, tm.tm_sec, ms);
      return buf;
    }

    const char* levelName(EventLevel level)
    {
      switch (level)
      {
        case EventLevel::Warn: return "WARN";
        case EventLevel::Error: return "ERROR";
        default: return "INFO";
      }
    }

    const char* kindName(OperatorKind kind)
    {
      switch (kind)
      {
        case OperatorKind::Ms1: return "MS1";
        case OperatorKind::Id: return "ID";
        case OperatorKind::Ms3: return "MS3";
        default: return "INFO";
      }
    }

    /// "<name>" padded to 5, then two spaces: the column every console and file line shares.
    std::string column(const char* name)
    {
      char buf[16];
      std::snprintf(buf, sizeof(buf), "%-5s  ", name);
      return buf;
    }

    /// Under the mutex: one line to the file (or the pending buffer), then every sink.
    void record_(State& s, const std::string& stamp, EventLevel level, const std::string& line)
    {
      std::string rt = "rt -";
      if (!std::isnan(s.rt_seconds))
      {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "rt %.2f", s.rt_seconds / 60.0);
        rt = buf;
      }
      std::string out = stamp + "  " + rt + "  " + column(levelName(level)) + line + "\n";
      if (s.file.is_open())
      {
        s.file.write(out.data(), static_cast<std::streamsize>(out.size()));
        s.file.flush();
      }
      else
      {
        s.pending.push_back(std::move(out));
        if (s.pending.size() > EventLog::kPendingCap) s.pending.pop_front();
      }
      for (const auto& kv : s.sinks) kv.second(line);
    }

    /// The KIND column's colour when the console takes SGR (ADR-0065 push 2).
    const char* colourOf(const char* name)
    {
      const std::string n(name);
      if (n == "ID") return "\x1b[32m";
      if (n == "MS3") return "\x1b[36m";
      if (n == "MS1") return "\x1b[2m";
      if (n == "WARN") return "\x1b[33m";
      if (n == "ERROR") return "\x1b[31m";
      return "";
    }

    void console_(State& s, const std::string& stamp, const char* name, const std::string& text)
    {
      if (s.console == nullptr) return;
      const char* colour = s.ansi ? colourOf(name) : "";
      const std::string out = stamp + "  " + colour + column(name) + (*colour ? "\x1b[0m" : "") + text + "\n";
      s.console->write(out.data(), static_cast<std::streamsize>(out.size()));
      s.console->flush();
    }
  } // namespace

  const char* EventLog::operatorTag(OperatorKind kind)
  {
    switch (kind)
    {
      case OperatorKind::Ms1: return "[OP-MS1]";
      case OperatorKind::Id: return "[OP-ID]";
      case OperatorKind::Ms3: return "[OP-MS3]";
      default: return "[OP-INFO]";
    }
  }

  void EventLog::record(EventLevel level, const std::string& line)
  {
    State& s = state();
    const std::string stamp = clock();
    std::lock_guard<std::mutex> lock(s.mutex);
    record_(s, stamp, level, line);
    if (level != EventLevel::Info) console_(s, stamp, levelName(level), line);
  }

  void EventLog::say(OperatorKind kind, const std::string& text)
  {
    State& s = state();
    const std::string stamp = clock();
    std::lock_guard<std::mutex> lock(s.mutex);
    record_(s, stamp, EventLevel::Info, std::string(operatorTag(kind)) + " " + text);
    console_(s, stamp, kindName(kind), text);
  }

  void EventLog::summary(const std::vector<std::string>& lines)
  {
    State& s = state();
    const std::string stamp = clock();
    std::lock_guard<std::mutex> lock(s.mutex);
    for (const std::string& l : lines)
    {
      record_(s, stamp, EventLevel::Info, "[OP-SUMMARY] " + l);
      if (s.console != nullptr)
      {
        const std::string out = l + "\n";
        s.console->write(out.data(), static_cast<std::streamsize>(out.size()));
      }
    }
    if (s.console != nullptr) s.console->flush();
  }

  void EventLog::paint(const std::string& raw)
  {
    State& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (s.console == nullptr) return;
    s.console->write(raw.data(), static_cast<std::streamsize>(raw.size()));
    s.console->flush();
  }

  void EventLog::setAnsi(bool ansi)
  {
    State& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.ansi = ansi;
  }

  void EventLog::setRt(double rt_seconds)
  {
    State& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.rt_seconds = rt_seconds;
  }

  void EventLog::openFile(const std::string& path)
  {
    State& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (path.empty())
    {
      // No run folder: the pending lines have nowhere to go. A file another engine opened stays open.
      if (!s.file.is_open()) s.pending.clear();
      return;
    }
    if (s.file.is_open()) s.file.close();
    s.rt_seconds = std::nan("");
    s.file.open(path, std::ios::app);
    if (s.file.is_open())
    {
      for (const std::string& l : s.pending) s.file.write(l.data(), static_cast<std::streamsize>(l.size()));
      s.file.flush();
    }
    s.pending.clear();
  }

  void EventLog::closeFile()
  {
    State& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    if (s.file.is_open()) s.file.close();
    s.rt_seconds = std::nan("");   // the run is over: a line before the next one carries no retention time
  }

  int EventLog::addSink(Sink sink)
  {
    State& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    const int h = s.next_handle++;
    s.sinks.emplace(h, std::move(sink));
    return h;
  }

  void EventLog::removeSink(int handle)
  {
    State& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.sinks.erase(handle);
  }

  void EventLog::setConsole(std::ostream* console)
  {
    State& s = state();
    std::lock_guard<std::mutex> lock(s.mutex);
    s.console = console;
  }
} // namespace OpenMS
