// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Tom David Mueller $
// $Authors: Tom David Mueller $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/ScanCommandQueue.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda/Config.h>

#include <atomic>
#include <set>
#include <thread>
#include <vector>

using namespace OpenMS;

namespace
{
  const char* minimal_config = R"({
    "deconvolution": {
      "min_charge": 4, "max_charge": 50,
      "min_mass": 500, "max_mass": 50000,
      "tol": [10, 10]
    },
    "precursor_selection": { "RT_window": 180, "target_mode": 0 },
    "tagging": {},
    "quantification": { "enabled": false },
    "faims": { "cv_values": [-50] },
    "ms_settings": {
      "ms1": {
        "analyzer": "Orbitrap", "first_mass": 500, "last_mass": 2000,
        "resolution": 120000, "agc_target": 800000, "max_it": 246
      },
      "ms2": [{ "analyzer": "Orbitrap", "activation": "HCD", "collision_energy": 29, "resolution": 120000 }]
    },
    "scheduling": {
      "cycle_time": { "enabled": false },
      "scan_timeout": { "enabled": true, "value_ms": 100 },
      "agc_interval_seconds": 30
    },
    "files": {},
    "selection_strategy": {
      "ms1": { "selection": "qscore", "max_precursors": 10 },
      "ms2": { "selection": "intensity" }
    }
  })";
}

START_TEST(ScanCommandQueue_Concurrent, "$Id$")

/////////////////////////////////////////////////////////////

// T1: Concurrent push/dequeue — no lost commands
START_SECTION(concurrent_push_dequeue)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  const int N_PRODUCERS = 4;
  const int CMDS_PER_PRODUCER = 250;
  const int TOTAL = N_PRODUCERS * CMDS_PER_PRODUCER;
  std::atomic<int> dequeued_count{0};

  // Producers push commands
  auto producer = [&](int thread_id)
  {
    for (int i = 0; i < CMDS_PER_PRODUCER; ++i)
    {
      ScanCommand cmd{};
      cmd.msn_level = 2;
      cmd.priority = 1;
      cmd.scan_id = thread_id * CMDS_PER_PRODUCER + i;
      queue.push(cmd);
    }
  };

  // Consumer drains the queue
  auto consumer = [&]()
  {
    while (dequeued_count.load() < TOTAL)
    {
      auto cmd = queue.dequeue();
      if (cmd.has_value())
        dequeued_count.fetch_add(1);
      else
        std::this_thread::yield();  // avoid starving producers
    }
  };

  std::vector<std::thread> threads;
  // Start consumer first so it's draining while producers push
  threads.emplace_back(consumer);
  for (int t = 0; t < N_PRODUCERS; ++t)
    threads.emplace_back(producer, t);

  for (auto& th : threads)
    th.join();

  // Drain any remaining
  while (auto cmd = queue.dequeue())
    dequeued_count.fetch_add(1);

  TEST_EQUAL(dequeued_count.load(), TOTAL)
}
END_SECTION

// T2: Concurrent tracking ID uniqueness
START_SECTION(concurrent_tracking_id_uniqueness)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  const int N_THREADS = 4;
  const int IDS_PER_THREAD = 250;
  const int TOTAL = N_THREADS * IDS_PER_THREAD;

  std::vector<std::vector<int>> per_thread_ids(N_THREADS);

  auto worker = [&](int thread_id)
  {
    per_thread_ids[thread_id].reserve(IDS_PER_THREAD);
    for (int i = 0; i < IDS_PER_THREAD; ++i)
    {
      per_thread_ids[thread_id].push_back(queue.nextTrackingId());
    }
  };

  std::vector<std::thread> threads;
  for (int t = 0; t < N_THREADS; ++t)
    threads.emplace_back(worker, t);
  for (auto& th : threads)
    th.join();

  std::set<int> all_ids;
  for (const auto& ids : per_thread_ids)
    for (int id : ids)
      all_ids.insert(id);

  TEST_EQUAL(static_cast<int>(all_ids.size()), TOTAL)
}
END_SECTION

// T3: Concurrent build + resolve — every built command is resolvable exactly once
START_SECTION(concurrent_build_resolve)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  const int N = 100;
  std::vector<std::atomic<int>> built_ids(N);
  for (int i = 0; i < N; ++i) built_ids[i].store(-1, std::memory_order_relaxed);
  std::atomic<int> resolved_count{0};
  std::atomic<int> double_resolve_count{0};

  // Producer: build MS2 commands concurrently (each writes to pending_scan_map_)
  auto builder = [&]()
  {
    for (int i = 0; i < N; ++i)
    {
      ScanCommand cmd = queue.buildMS2(500.0 + i, 10, 29.0, "HCD");
      built_ids[i].store(cmd.scan_id, std::memory_order_release);
    }
  };

  // Consumer: resolve pending commands by ID concurrently with builder
  auto resolver = [&]()
  {
    int local_resolved = 0;
    while (local_resolved < N)
    {
      bool made_progress = false;
      for (int i = 0; i < N; ++i)
      {
        int id = built_ids[i].load(std::memory_order_acquire);
        if (id < 0) continue;  // not yet built or already resolved
        auto result = queue.resolvePending(id);
        if (result.has_value())
        {
          local_resolved++;
          made_progress = true;
          // Try resolving again — must return nullopt (exactly-once guarantee)
          auto duplicate = queue.resolvePending(id);
          if (duplicate.has_value())
            double_resolve_count.fetch_add(1);
          built_ids[i].store(-2, std::memory_order_relaxed);  // mark as resolved
        }
      }
      if (!made_progress)
        std::this_thread::yield();  // avoid starving builder
    }
    resolved_count.store(local_resolved);
  };

  // Run concurrently: builder writes to pending_scan_map_ via buildMS2,
  // resolver reads+erases via resolvePending. Both acquire queue_mutex_ internally.
  std::thread build_thread(builder);
  std::thread resolve_thread(resolver);
  build_thread.join();
  resolve_thread.join();

  TEST_EQUAL(resolved_count.load(), N)
  TEST_EQUAL(double_resolve_count.load(), 0)
}
END_SECTION

// T4: Concurrent push + cleanupExpired — no crashes
START_SECTION(concurrent_push_cleanup)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  const int N = 200;
  std::atomic<bool> done{false};

  // Producer pushes commands with timestamp 0 (will be expired immediately)
  auto pusher = [&]()
  {
    for (int i = 0; i < N; ++i)
    {
      ScanCommand cmd{};
      cmd.msn_level = 2;
      cmd.priority = 1;
      cmd.scan_id = i;
      cmd.enqueue_timestamp_ms = 1;  // old timestamp -> will expire
      queue.push(cmd);
      // Also register in pending map so cleanupExpired has something to clean
      queue.registerPending(i, cmd);
    }
    done.store(true);
  };

  // Cleaner runs cleanupExpired concurrently
  auto cleaner = [&]()
  {
    while (!done.load())
    {
      queue.cleanupExpired();
      std::this_thread::yield();  // avoid starving pusher
    }
    // Final cleanup
    queue.cleanupExpired();
  };

  std::thread push_thread(pusher);
  std::thread clean_thread(cleaner);
  push_thread.join();
  clean_thread.join();

  // If we get here without crash/hang, the test passes
  TEST_EQUAL(true, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
