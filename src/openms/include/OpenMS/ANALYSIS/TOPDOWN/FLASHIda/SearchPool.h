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

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace OpenMS
{
  /**
   * @brief A fixed pool of search workers with a bounded backlog (ADR-0062).
   *
   * A job is three steps, and the pool runs the first two on a worker:
   *   - search(): the expensive part, run with NO lock held. It may read only what the job copied and
   *     what is immutable after load (the Config).
   *   - apply():  the result's effect on the engine. It takes analysis_mutex_ ITSELF -- the engine's
   *     second acquisition site -- and returns the number of commands it pushed.
   *   - drop():   run instead of both when the backlog is full and this is the oldest waiting job. It is
   *     run by the SUBMITTER, which already holds analysis_mutex_ (submit() hands the job back for that).
   *
   * The pool knows nothing about spectra or the engine; FLASHIda builds the jobs.
   *
   * Lock hierarchy: analysis_mutex_ -> pool mutex (a leaf). A worker never holds the pool mutex while it
   * runs search() or apply(), so apply() taking analysis_mutex_ cannot invert the edge. waitIdle() blocks
   * on the pool and must never be called with analysis_mutex_ held: a worker's apply() would wait for the
   * caller forever.
   */
  class OPENMS_DLLAPI SearchPool
  {
  public:
    struct Job
    {
      std::function<void()> search;
      std::function<int()> apply;
      std::function<void()> drop;
    };

    /// @p threads workers (at least 1) and a backlog of @p backlog WAITING jobs (running ones excluded).
    /// @p on_thread_start runs once on each worker before its first job (FLASHIda sets OpenMP's per-thread
    /// team size there, which a new thread does not inherit).
    SearchPool(int threads, size_t backlog, std::function<void()> on_thread_start = {});

    /// stop()
    ~SearchPool();

    SearchPool(const SearchPool&) = delete;
    SearchPool& operator=(const SearchPool&) = delete;

    /// Queue @p job. When the backlog is already full the OLDEST waiting job is taken out and returned, and
    /// the caller runs its drop() (ADR-0062 d.5: the oldest is the least useful -- its window is nearest
    /// expiry). Returns nullptr when nothing was dropped.
    std::unique_ptr<Job> submit(std::unique_ptr<Job> job);

    /// Block until no job is waiting or running. Returns the commands the applies pushed since the previous
    /// call. Harnesses only (WaitForSearchIdle): production never waits on the pool.
    int waitIdle();

    /// Stop accepting work, discard the waiting jobs, let a running job finish (its apply included) and join
    /// every worker. Idempotent.
    void stop();

    size_t backlog() const;
    void setBacklog(size_t backlog);   ///< test seam (FLASHIdaTestAccess::setSearchBacklog)
    int threadCount() const { return static_cast<int>(workers_.size()); }
    int searchesFinished() const { return searches_finished_.load(); }
    int appliesFinished() const { return applies_finished_.load(); }

  private:
    void run_();

    mutable std::mutex mutex_;
    std::condition_variable work_cv_;
    std::condition_variable idle_cv_;
    std::deque<std::unique_ptr<Job>> waiting_;
    int in_flight_ = 0;
    size_t backlog_;
    bool stopping_ = false;
    std::function<void()> on_thread_start_;
    std::atomic<int> pushed_since_wait_{0};
    std::atomic<int> searches_finished_{0};
    std::atomic<int> applies_finished_{0};
    std::vector<std::thread> workers_;   // last: the workers start once everything above exists
  };
} // namespace OpenMS
