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

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ConsoleView.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/EventLog.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>

#ifdef _WIN32
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
  #ifndef WIN32_LEAN_AND_MEAN
    #define WIN32_LEAN_AND_MEAN
  #endif
  #include <windows.h>
#endif

// Every non-ASCII glyph is written as escaped UTF-8 bytes, so the source compiles the same under any code page.
namespace OpenMS
{
  namespace
  {
    const char* const kLevelsUtf8[] = {" ", "\xE2\x96\x81", "\xE2\x96\x82", "\xE2\x96\x83", "\xE2\x96\x84",
                                       "\xE2\x96\x85", "\xE2\x96\x86", "\xE2\x96\x87", "\xE2\x96\x88"};
    const char* const kLevelsAscii[] = {" ", ":", "-", "=", "+", "*", "#", "%", "@"};
    constexpr int kLevels = 9;

    const char* future(bool utf8) { return utf8 ? "\xC2\xB7" : "."; }
    const char* sep(bool utf8) { return utf8 ? " \xC2\xB7 " : " | "; }
    const char* rule(bool utf8) { return utf8 ? "\xE2\x94\x80" : "-"; }

    std::string fmtd(const char* f, double v)
    {
      char buf[64];
      std::snprintf(buf, sizeof(buf), f, v);
      return buf;
    }

    const char* stateName(PrecursorState s)
    {
      switch (s)
      {
        case PrecursorState::Acquiring: return "ACQUIRING";
        case PrecursorState::Searching: return "SEARCHING";
        case PrecursorState::NoId: return "NO-ID";
        case PrecursorState::Identified: return "IDENTIFIED";
        case PrecursorState::ClosedOut: return "CLOSED-OUT";
        default: return "ACQUIRED";
      }
    }

    const char* stateLower(PrecursorState s)
    {
      switch (s)
      {
        case PrecursorState::Acquiring: return "acquiring";
        case PrecursorState::Searching: return "searching";
        case PrecursorState::NoId: return "no-id";
        case PrecursorState::Identified: return "identified";
        case PrecursorState::ClosedOut: return "closed out";
        default: return "acquired";
      }
    }

    const char* stateColour(PrecursorState s)
    {
      switch (s)
      {
        case PrecursorState::Acquiring:
        case PrecursorState::Searching: return "\x1b[33m";
        case PrecursorState::NoId: return "\x1b[2m";
        case PrecursorState::Identified: return "\x1b[36m";
        case PrecursorState::ClosedOut: return "\x1b[32m";
        default: return "\x1b[2m";
      }
    }

    /// The prefix of @p s that fits @p cols display columns (escape sequences kept whole, counted as zero).
    std::string clip(const std::string& s, int cols)
    {
      std::string out;
      int used = 0;
      for (size_t i = 0; i < s.size();)
      {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        if (c == 0x1b)
        {
          size_t j = i + 1;
          if (j < s.size() && s[j] == '[')
            for (++j; j < s.size() && !(s[j] >= '@' && s[j] <= '~'); ++j) {}
          out.append(s, i, j + 1 - i);
          i = j + 1;
          continue;
        }
        size_t len = 1;
        if (c >= 0xF0) len = 4;
        else if (c >= 0xE0) len = 3;
        else if (c >= 0xC0) len = 2;
        if (used + 1 > cols) break;
        out.append(s, i, len);
        ++used;
        i += len;
      }
      return out;
    }

    /// The nearest-rank @p p-th percentile of @p v (sorted copy); -1 for an empty list.
    int percentile(std::vector<int> v, double p)
    {
      if (v.empty()) return -1;
      std::sort(v.begin(), v.end());
      size_t rank = static_cast<size_t>(std::ceil(p / 100.0 * static_cast<double>(v.size())));
      rank = std::max<size_t>(rank, 1);
      return v[std::min(rank, v.size()) - 1];
    }

    std::string statusLine(const BoardSnapshot& s, bool utf8)
    {
      const char* d = sep(utf8);
      std::string out = "rt " + fmtd("%.1f", s.rt_s / 60.0);
      if (s.duration_min > 0) out += " / " + fmtd("%.0f", s.duration_min) + " min";
      else out += " min";
      out += "  MS1 " + std::to_string(s.ms1) + d + "MS2 " + std::to_string(s.ms2) + d + "MS3 " + std::to_string(s.ms3);
      out += "   live " + std::to_string(s.live.size()) + d + "identified " + std::to_string(s.identified) + d + "closed out "
             + std::to_string(s.closed_out) + d + "plan " + std::to_string(s.planned);
      out += "   backlog " + std::to_string(s.backlog);
      return out;
    }
  } // namespace

  int ConsoleView::columns(const std::string& s)
  {
    int n = 0;
    for (size_t i = 0; i < s.size(); ++i)
    {
      const unsigned char c = static_cast<unsigned char>(s[i]);
      if (c == 0x1b)
      {
        size_t j = i + 1;
        if (j < s.size() && s[j] == '[')
          for (++j; j < s.size() && !(s[j] >= '@' && s[j] <= '~'); ++j) {}
        i = j;
        continue;
      }
      if ((c & 0xC0) != 0x80) ++n;   // not a UTF-8 continuation byte
    }
    return n;
  }

  std::string ConsoleView::timelineCells(const BoardSnapshot& snap, int cells, bool utf8)
  {
    if (cells <= 0) return std::string();
    const double now_min = snap.rt_s / 60.0;
    const double span = snap.duration_min > 0 ? snap.duration_min : std::max(now_min, 1.0);
    const double per = span / cells;
    std::vector<int> count(cells, 0);
    for (double rt : snap.id_rts_s)
    {
      const int c = static_cast<int>(std::floor(rt / 60.0 / per));
      if (c >= 0 && c < cells) ++count[c];
    }
    const int busiest = *std::max_element(count.begin(), count.end());
    const int now_cell = static_cast<int>(std::floor(now_min / per));
    std::string out;
    for (int i = 0; i < cells; ++i)
    {
      if (i > now_cell)
      {
        out += future(utf8);
        continue;
      }
      int level = 0;
      if (count[i] > 0 && busiest > 0) level = 1 + (count[i] * (kLevels - 2) + busiest - 1) / busiest;   // 1..8, the busiest 8
      level = std::min(level, kLevels - 1);
      out += utf8 ? kLevelsUtf8[level] : kLevelsAscii[level];
    }
    return out;
  }

  std::vector<std::string> ConsoleView::renderPanel(const BoardSnapshot& snap, int width, int height, bool utf8, bool ansi)
  {
    std::vector<std::string> lines;
    if (width <= 0 || height <= 0) return lines;

    // 1. the status line, under a rule
    std::string status = std::string(rule(utf8)) + rule(utf8) + " " + statusLine(snap, utf8) + " ";
    while (columns(status) < width) status += rule(utf8);
    lines.push_back(clip(status, width));

    // 2. the identification timeline
    const int cells = std::max(10, width - 24);
    lines.push_back(clip(std::string("ids ") + (utf8 ? "\xE2\x96\x95" : "|") + timelineCells(snap, cells, utf8)
                             + (utf8 ? "\xE2\x96\x8F" : "|") + " " + std::to_string(snap.identified) + " identified",
                         width));

    // 3. the header
    lines.push_back(clip(" nominal       mass   z  state        score  MS2  MS3      rt  protein", width));

    // 4. the rows that fit whole, most recent first; the rest go to the footer
    const int wrap = std::max(10, width - 9);
    const int room = height - static_cast<int>(lines.size()) - 1;   // one line for the footer
    int used = 0;
    size_t shown = 0;
    for (; shown < snap.live.size(); ++shown)
    {
      const BoardRow& r = snap.live[shown];
      const int pf_lines = r.proforma.empty() ? 0 : static_cast<int>((r.proforma.size() + wrap - 1) / wrap);
      if (used + 1 + pf_lines > room) break;
      char head[160];
      std::snprintf(head, sizeof(head), "%8d %10.2f %3d  ", r.nominal, r.mass, r.charge);
      std::string state = stateName(r.state);
      state.resize(11, ' ');
      std::string score = r.score >= 0 ? fmtd("%6.0f", r.score) : std::string("     -");
      char tail[96];
      std::snprintf(tail, sizeof(tail), "  %4d %4d  %6.2f  ", r.ms2, r.ms3, r.last_rt_s / 60.0);
      std::string row = std::string(head) + (ansi ? stateColour(r.state) : "") + state + (ansi ? "\x1b[0m" : "") + score + tail
                        + r.accession;
      if (r.below_threshold) row += "  (below " + fmtd("%.0f", r.threshold) + ")";
      lines.push_back(clip(row, width));
      for (int k = 0; k < pf_lines; ++k)
        lines.push_back(clip(std::string(9, ' ') + r.proforma.substr(static_cast<size_t>(k) * wrap, wrap), width));
      used += 1 + pf_lines;
    }

    // 5. pad, then the footer on the last line
    while (static_cast<int>(lines.size()) < height - 1) lines.push_back(std::string());
    std::string footer;
    if (shown < snap.live.size())
    {
      std::map<PrecursorState, int> by_state;
      for (size_t i = shown; i < snap.live.size(); ++i) ++by_state[snap.live[i].state];
      footer = " + " + std::to_string(snap.live.size() - shown) + " more live (";
      bool first = true;
      for (const auto& kv : by_state)
      {
        footer += std::string(first ? "" : sep(utf8)) + std::to_string(kv.second) + " " + stateLower(kv.first);
        first = false;
      }
      footer += ")";
    }
    lines.push_back(clip(footer, width));
    lines.resize(static_cast<size_t>(height));
    return lines;
  }

  std::vector<std::string> ConsoleView::renderSummary(const BoardSnapshot& snap, bool utf8)
  {
    const char* d = sep(utf8);
    std::vector<std::string> out;
    std::string head = std::string(rule(utf8)) + rule(utf8) + " run summary" + d + "rt " + fmtd("%.1f", snap.rt_s / 60.0) + " min ";
    for (int i = 0; i < 30; ++i) head += rule(utf8);
    out.push_back(head);
    out.push_back("scans       MS1 " + std::to_string(snap.ms1) + d + "MS2 " + std::to_string(snap.ms2) + d + "MS3 "
                  + std::to_string(snap.ms3));
    out.push_back("precursors  " + std::to_string(snap.precursors) + d + "identified " + std::to_string(snap.identified) + d
                  + "closed out " + std::to_string(snap.closed_out) + d + "no-id " + std::to_string(snap.no_id));
    std::string proteins = "proteins    " + std::to_string(snap.proteins);
    if (snap.best_score >= 0) proteins += std::string(d) + "best score " + fmtd("%.0f", snap.best_score) + " (" + snap.best_accession + ")";
    out.push_back(proteins);
    if (snap.search_wait_ms.empty())
      out.push_back("search      none");
    else
      out.push_back("search      wait median " + std::to_string(percentile(snap.search_wait_ms, 50)) + " ms" + d + "p95 "
                    + std::to_string(percentile(snap.search_wait_ms, 95)) + " ms" + d + "n "
                    + std::to_string(snap.search_wait_ms.size()));
    out.push_back(std::string("timeline    ") + (utf8 ? "\xE2\x96\x95" : "|") + timelineCells(snap, 60, utf8)
                  + (utf8 ? "\xE2\x96\x8F" : "|") + " " + std::to_string(snap.identified) + " identified");
    if (!snap.log_dir.empty()) out.push_back("logs        " + snap.log_dir);
    return out;
  }

  void ConsoleView::printSummary(const BoardSnapshot& snap, bool utf8)
  {
    EventLog::summary(renderSummary(snap, utf8));
  }

  std::unique_ptr<ConsoleView> ConsoleView::start(std::function<BoardSnapshot()> snapshot)
  {
#ifdef _WIN32
    HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    if (h == INVALID_HANDLE_VALUE || h == nullptr || !GetConsoleMode(h, &mode)) return nullptr;   // a pipe or a file
    if (!SetConsoleMode(h, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING)) return nullptr;            // no VT: plain lines
    std::unique_ptr<ConsoleView> v(new ConsoleView());
    v->handle_ = h;
    v->old_mode_ = mode;
    v->old_cp_ = GetConsoleOutputCP();
    v->utf8_ = SetConsoleOutputCP(CP_UTF8) != 0;
    v->snapshot_ = std::move(snapshot);
    EventLog::setAnsi(true);
    v->layout_();
    v->running_ = true;
    ConsoleView* self = v.get();
    v->thread_ = std::thread([self] {
      while (self->running_.load())
      {
        self->redraw_();
        for (int i = 0; i < 5 && self->running_.load(); ++i) std::this_thread::sleep_for(std::chrono::milliseconds(50));
      }
    });
    return v;
#else
    (void)snapshot;
    return nullptr;
#endif
  }

  void ConsoleView::layout_()
  {
#ifdef _WIN32
    CONSOLE_SCREEN_BUFFER_INFO info;
    if (!GetConsoleScreenBufferInfo(static_cast<HANDLE>(handle_), &info)) return;
    width_ = info.srWindow.Right - info.srWindow.Left + 1;
    height_ = info.srWindow.Bottom - info.srWindow.Top + 1;
    panel_ = height_ >= 16 ? height_ / 2 : 0;   // too small a window: plain lines only
    std::string raw = "\x1b[r";
    if (panel_ > 0)
    {
      raw += std::string(static_cast<size_t>(panel_), '\n');   // push what is on screen up, out of the panel's rows
      raw += "\x1b[1;" + std::to_string(height_ - panel_) + "r";
      raw += "\x1b[" + std::to_string(height_ - panel_) + ";1H";
    }
    EventLog::paint(raw);
#endif
  }

  void ConsoleView::redraw_()
  {
#ifdef _WIN32
    CONSOLE_SCREEN_BUFFER_INFO info;
    if (GetConsoleScreenBufferInfo(static_cast<HANDLE>(handle_), &info))
    {
      const int w = info.srWindow.Right - info.srWindow.Left + 1;
      const int h = info.srWindow.Bottom - info.srWindow.Top + 1;
      if (w != width_ || h != height_) layout_();
    }
    if (panel_ <= 0) return;
    const std::vector<std::string> lines = renderPanel(snapshot_(), width_, panel_, utf8_, true);
    std::string raw = "\x1b" "7";
    for (size_t i = 0; i < lines.size(); ++i)
      raw += "\x1b[" + std::to_string(height_ - panel_ + 1 + static_cast<int>(i)) + ";1H" + lines[i] + "\x1b[0m\x1b[K";
    raw += "\x1b" "8";
    EventLog::paint(raw);
#endif
  }

  void ConsoleView::stop()
  {
    if (stopped_) return;
    stopped_ = true;
    running_ = false;
    if (thread_.joinable()) thread_.join();
#ifdef _WIN32
    redraw_();   // the final frame
    EventLog::paint("\x1b[r\x1b[" + std::to_string(std::max(height_, 1)) + ";1H\n");
#endif
    EventLog::setAnsi(false);
  }

  ConsoleView::~ConsoleView()
  {
    stop();
#ifdef _WIN32
    if (handle_ != nullptr)
    {
      SetConsoleMode(static_cast<HANDLE>(handle_), old_mode_);
      if (old_cp_ != 0) SetConsoleOutputCP(old_cp_);
    }
#endif
  }
} // namespace OpenMS
