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


#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/SearchPool.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/EventLog.h>

#include <algorithm>
#include <exception>
#include <iostream>

namespace OpenMS
{
  SearchPool::SearchPool(int threads, size_t backlog, std::function<void()> on_thread_start) :
      backlog_(std::max<size_t>(1, backlog)),
      on_thread_start_(std::move(on_thread_start))
  {
    const int n = std::max(1, threads);
    workers_.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i)
      workers_.emplace_back([this] { run_(); });
  }

  SearchPool::~SearchPool()
  {
    stop();
  }

  std::unique_ptr<SearchPool::Job> SearchPool::submit(std::unique_ptr<Job> job)
  {
    std::unique_ptr<Job> dropped;
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (stopping_) return job;   // after stop() nothing runs: the caller drops it like a backlog drop
      if (waiting_.size() >= backlog_)
      {
        dropped = std::move(waiting_.front());
        waiting_.pop_front();
      }
      waiting_.push_back(std::move(job));
    }
    work_cv_.notify_one();
    return dropped;
  }

  int SearchPool::waitIdle()
  {
    std::unique_lock<std::mutex> lock(mutex_);
    idle_cv_.wait(lock, [this] { return (waiting_.empty() && in_flight_ == 0) || stopping_; });
    return pushed_since_wait_.exchange(0);
  }

  void SearchPool::stop()
  {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      if (stopping_ && workers_.empty()) return;
      stopping_ = true;
      waiting_.clear();
    }
    work_cv_.notify_all();
    idle_cv_.notify_all();
    for (auto& w : workers_)
      if (w.joinable()) w.join();
    workers_.clear();
  }

  size_t SearchPool::backlog() const
  {
    std::lock_guard<std::mutex> lock(mutex_);
    return backlog_;
  }

  size_t SearchPool::waiting() const
  {
    std::lock_guard<std::mutex> lock(mutex_);
    return waiting_.size();
  }

  void SearchPool::setBacklog(size_t backlog)
  {
    std::lock_guard<std::mutex> lock(mutex_);
    backlog_ = std::max<size_t>(1, backlog);
  }

  void SearchPool::run_()
  {
    if (on_thread_start_) on_thread_start_();
    for (;;)
    {
      std::unique_ptr<Job> job;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        work_cv_.wait(lock, [this] { return stopping_ || !waiting_.empty(); });
        if (stopping_) return;
        job = std::move(waiting_.front());
        waiting_.pop_front();
        ++in_flight_;
      }

      // A throw must not escape a worker (std::terminate would take the instrument down with it). The
      // job's scan then identifies nothing, like a dropped one; the message says which step threw.
      int pushed = 0;
      try
      {
        if (job->search) job->search();
        searches_finished_.fetch_add(1);
        if (job->apply) pushed = job->apply();
      }
      catch (const std::exception& e)
      {
        EventLine(EventLevel::Error) << "[SEARCH-ERROR] " << e.what();
      }
      catch (...)
      {
        EventLine(EventLevel::Error) << "[SEARCH-ERROR] unknown exception";
      }
      applies_finished_.fetch_add(1);
      pushed_since_wait_.fetch_add(pushed);

      {
        std::lock_guard<std::mutex> lock(mutex_);
        --in_flight_;
        if (waiting_.empty() && in_flight_ == 0) idle_cv_.notify_all();
      }
    }
  }
} // namespace OpenMS
