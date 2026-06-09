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
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>

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
      "ms1": { "selection": "qscore", "max_targets": 10 },
      "ms2": { "selection": "none" }
    }
  })";
}

START_TEST(ScanCommandQueue_Concurrent, "$Id$")

/////////////////////////////////////////////////////////////

// T1: Concurrent push from multiple threads, then concurrent dequeue.
// Phase-based: no spin-waits, no live consumer during push.
START_SECTION(concurrent_push_dequeue)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  const int N_THREADS = 4;
  const int CMDS_PER_THREAD = 250;
  const int TOTAL = N_THREADS * CMDS_PER_THREAD;

  // Phase 1: 4 threads push concurrently
  {
    std::vector<std::thread> threads;
    for (int t = 0; t < N_THREADS; ++t)
    {
      threads.emplace_back([&, t]()
      {
        for (int i = 0; i < CMDS_PER_THREAD; ++i)
        {
          ScanCommand cmd{};
          cmd.msn_level = 2;
          cmd.priority = 1;
          cmd.scan_id = t * CMDS_PER_THREAD + i;
          queue.push(cmd);
        }
      });
    }
    for (auto& th : threads) th.join();
  }

  // Phase 2: 4 threads dequeue concurrently
  std::atomic<int> dequeued_count{0};
  {
    std::vector<std::thread> threads;
    for (int t = 0; t < N_THREADS; ++t)
    {
      threads.emplace_back([&]()
      {
        while (true)
        {
          auto cmd = queue.dequeue();
          if (cmd.has_value())
            dequeued_count.fetch_add(1);
          else
            break;  // queue empty, done
        }
      });
    }
    for (auto& th : threads) th.join();
  }

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

// T3: Concurrent resolve — build first, then multiple resolvers race.
// Tests that each pending command is resolved exactly once.
START_SECTION(concurrent_build_resolve)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  const int N = 100;

  // Phase 1: build and push all commands (single-threaded)
  ScanConfig sc;
  sc.collision_energy = 29;
  sc.activation = "HCD";
  std::vector<int> built_ids(N);
  for (int i = 0; i < N; ++i)
  {
    PeakGroup pg(10, 10, true);
    pg.setMonoisotopicMass(500.0 + i);
    FLASHHelperClasses::LogMzPeak lp;
    lp.mz = 500.0 + i;
    lp.abs_charge = 10;
    pg.push_back(lp);
    ScanCommand cmd = queue.buildMS2(pg, 10, sc);
    built_ids[i] = cmd.scan_id;
    queue.push(cmd);
  }

  // Phase 1b: dequeue all commands so they enter pending_scan_map_
  // (push() only enqueues; dequeue() stamps dequeue_timestamp_ms and inserts into map)
  for (int i = 0; i < N; ++i)
  {
    auto dequeued = queue.dequeue();
    TEST_EQUAL(dequeued.has_value(), true)
  }

  // Phase 2: 4 resolver threads race to resolve
  const int N_RESOLVERS = 4;
  std::vector<std::atomic<int>> per_thread_resolved(N_RESOLVERS);
  for (int t = 0; t < N_RESOLVERS; ++t) per_thread_resolved[t].store(0);

  {
    std::vector<std::thread> threads;
    for (int t = 0; t < N_RESOLVERS; ++t)
    {
      threads.emplace_back([&, t]()
      {
        int local = 0;
        // Each thread tries to resolve ALL ids — only one thread will succeed per id
        for (int i = 0; i < N; ++i)
        {
          auto result = queue.resolvePending(built_ids[i]);
          if (result.has_value())
            local++;
        }
        per_thread_resolved[t].store(local);
      });
    }
    for (auto& th : threads) th.join();
  }

  // Total resolved across all threads must equal N (each resolved exactly once)
  int total_resolved = 0;
  for (int t = 0; t < N_RESOLVERS; ++t)
    total_resolved += per_thread_resolved[t].load();

  TEST_EQUAL(total_resolved, N)
}
END_SECTION

// T4: Concurrent push + cleanupExpired — no crashes, no data corruption.
// Phase-based: pusher and cleaner run concurrently but pusher is bounded.
START_SECTION(concurrent_push_cleanup)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  const int N = 200;
  std::atomic<bool> done{false};

  // Pusher pushes commands with current timestamps (exercises concurrent push + cleanup)
  auto pusher = [&]()
  {
    for (int i = 0; i < N; ++i)
    {
      ScanCommand cmd{};
      cmd.msn_level = 2;
      cmd.priority = 1;
      cmd.scan_id = i;
      queue.push(cmd);
    }
    done.store(true, std::memory_order_release);
  };

  // Cleaner runs cleanupExpired concurrently
  auto cleaner = [&]()
  {
    while (!done.load(std::memory_order_acquire))
    {
      queue.cleanupExpired();
      std::this_thread::yield();
    }
    queue.cleanupExpired();
  };

  std::thread push_thread(pusher);
  std::thread clean_thread(cleaner);
  push_thread.join();
  clean_thread.join();

  // No command reaches the 100 ms timeout during this sub-millisecond test, so the
  // concurrent cleaner must not have dropped any of the 200 pushes — verifies no loss /
  // no corruption under the push+cleanup race (not just "didn't crash").
  TEST_EQUAL(queue.queueSize(1), N)
}
END_SECTION

// T5: cancelByScanIds removes matching queued commands and returns their ids.
START_SECTION(cancelByScanIds_removes_queued_and_returns_ids)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  for (int id : {10, 20, 30, 40})
  {
    ScanCommand cmd{};
    cmd.msn_level = 2;
    cmd.priority = 1;
    cmd.scan_id = id;
    queue.push(cmd);
  }
  TEST_EQUAL(static_cast<int>(queue.queueSize(1)), 4)

  auto removed = queue.cancelByScanIds({20, 40});
  std::set<int> removed_set(removed.begin(), removed.end());
  TEST_EQUAL(removed_set.size() == 2 && removed_set.count(20) == 1 && removed_set.count(40) == 1, true)
  TEST_EQUAL(static_cast<int>(queue.queueSize(1)), 2)   // 10, 30 remain

  auto removed2 = queue.cancelByScanIds({10, 30});
  TEST_EQUAL(static_cast<int>(removed2.size()), 2)
  TEST_EQUAL(static_cast<int>(queue.queueSize(1)), 0)
}
END_SECTION

// T6: cancelByScanIds also removes in-flight (pending_scan_map_) commands.
START_SECTION(cancelByScanIds_removes_inflight_pending)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  ScanCommand cmd{};
  cmd.msn_level = 2;
  cmd.priority = 1;
  cmd.scan_id = 77;
  queue.push(cmd);
  auto dq = queue.dequeue();                            // 77: queues_ -> pending_scan_map_ (in-flight)
  TEST_EQUAL(dq.has_value(), true)
  TEST_EQUAL(queue.peekPending(77).has_value(), true)

  auto removed = queue.cancelByScanIds({77});
  TEST_EQUAL(static_cast<int>(removed.size()), 1)
  TEST_EQUAL(removed.empty() ? -1 : removed[0], 77)
  TEST_EQUAL(queue.peekPending(77).has_value(), false)  // also cleared from pending_scan_map_
}
END_SECTION

// T7: cancelByScanIds with no matching id is a harmless no-op.
START_SECTION(cancelByScanIds_no_match_is_noop)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  ScanCommand cmd{};
  cmd.msn_level = 2;
  cmd.priority = 1;
  cmd.scan_id = 5;
  queue.push(cmd);

  auto removed = queue.cancelByScanIds({9999});
  TEST_EQUAL(removed.empty(), true)
  TEST_EQUAL(static_cast<int>(queue.queueSize(1)), 1)   // untouched
}
END_SECTION

// T8: cancelByScanIds spans all priority queues; only targeted ids are removed.
START_SECTION(cancelByScanIds_spans_priorities)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  for (int p = 0; p < 4; ++p)
  {
    ScanCommand cmd{};
    cmd.msn_level = 2;
    cmd.priority = p;
    cmd.scan_id = 100 + p;                              // ids 100,101,102,103
    queue.push(cmd);
  }

  auto removed = queue.cancelByScanIds({100, 102});     // priorities 0 and 2
  std::set<int> removed_set(removed.begin(), removed.end());
  TEST_EQUAL(removed_set.size() == 2 && removed_set.count(100) == 1 && removed_set.count(102) == 1, true)
  TEST_EQUAL(static_cast<int>(queue.queueSize(0)), 0)
  TEST_EQUAL(static_cast<int>(queue.queueSize(1)), 1)   // 101 intact
  TEST_EQUAL(static_cast<int>(queue.queueSize(2)), 0)
  TEST_EQUAL(static_cast<int>(queue.queueSize(3)), 1)   // 103 intact
}
END_SECTION

// T9: Concurrent cancelByScanIds over disjoint id ranges while pushing — no races, each id removed once.
START_SECTION(cancelByScanIds_concurrent)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  const int N = 1000;
  for (int i = 0; i < N; ++i)
  {
    ScanCommand cmd{};
    cmd.msn_level = 2;
    cmd.priority = 1;
    cmd.scan_id = i;
    queue.push(cmd);
  }

  const int N_CANCELERS = 4;
  const int PER = N / N_CANCELERS;                      // 250
  std::vector<std::atomic<int>> per_thread_removed(N_CANCELERS);
  for (int t = 0; t < N_CANCELERS; ++t) per_thread_removed[t].store(0);

  // One extra pusher adds ids 1000..1249 concurrently (never targeted by any canceler).
  auto pusher = [&]()
  {
    for (int i = N; i < N + 250; ++i)
    {
      ScanCommand cmd{};
      cmd.msn_level = 2;
      cmd.priority = 1;
      cmd.scan_id = i;
      queue.push(cmd);
    }
  };

  auto canceler = [&](int t)
  {
    std::vector<int> ids;
    ids.reserve(PER);
    for (int i = t * PER; i < (t + 1) * PER; ++i) ids.push_back(i);
    auto removed = queue.cancelByScanIds(ids);
    per_thread_removed[t].store(static_cast<int>(removed.size()));
  };

  std::vector<std::thread> threads;
  threads.emplace_back(pusher);
  for (int t = 0; t < N_CANCELERS; ++t) threads.emplace_back(canceler, t);
  for (auto& th : threads) th.join();

  int total_removed = 0;
  for (int t = 0; t < N_CANCELERS; ++t) total_removed += per_thread_removed[t].load();

  TEST_EQUAL(total_removed, N)                          // each of the 1000 original ids removed exactly once
  TEST_EQUAL(static_cast<int>(queue.queueSize(1)), 250) // the 250 newly pushed (untargeted) ids remain
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
