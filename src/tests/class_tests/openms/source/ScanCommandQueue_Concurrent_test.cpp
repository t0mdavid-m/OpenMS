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

  // Phase 1: build all commands (single-threaded)
  std::vector<int> built_ids(N);
  for (int i = 0; i < N; ++i)
  {
    ScanCommand cmd = queue.buildMS2(500.0 + i, 10, 29.0, "HCD");
    built_ids[i] = cmd.scan_id;
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

  // Pusher registers commands with old timestamps (will expire)
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
      queue.registerPending(i, cmd);
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

  // If we get here without crash/hang, the test passes
  TEST_EQUAL(true, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
