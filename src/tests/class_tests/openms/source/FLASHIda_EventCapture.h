// FLASHIda test helper: capture the engine's event lines (ADR-0065).
//
// The engine no longer prints to std::cout; its events go to EventLog's sinks. EventCapture is one sink for the
// lifetime of the object: str() is every line emitted meanwhile, one per '\n', as the line was emitted (no
// events.log prefix), operator lines included as "[OP-<KIND>] text". Thread-safe: the search pool emits from its
// workers. Test scaffolding, so it lives here and not in the engine headers.

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/EventLog.h>

#include <mutex>
#include <string>

class EventCapture
{
public:
  EventCapture() :
      handle_(OpenMS::EventLog::addSink([this](const std::string& line) {
        std::lock_guard<std::mutex> lock(mutex_);
        text_ += line;
        text_ += '\n';
      }))
  {
  }
  EventCapture(const EventCapture&) = delete;
  EventCapture& operator=(const EventCapture&) = delete;
  ~EventCapture() { OpenMS::EventLog::removeSink(handle_); }

  std::string str() const
  {
    std::lock_guard<std::mutex> lock(mutex_);
    return text_;
  }

private:
  mutable std::mutex mutex_;
  std::string text_;
  int handle_;   // last: initialised after the buffer the sink writes into
};
