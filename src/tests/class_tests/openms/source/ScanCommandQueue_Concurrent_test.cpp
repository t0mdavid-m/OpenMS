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

#include "FLASHIda_TestAccess.h"   // the queue's plan registry, read without a public getter

#include <atomic>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace OpenMS;

namespace
{
  const char* minimal_config = R"({
  "deconvolution": {
    "min_charge": 4,
    "max_charge": 50,
    "min_mass": 500,
    "max_mass": 50000,
    "tol": [
      10,
      10,
      10
    ]
  },
  "faims": {
    "cv_values": []
  },
  "scheduling": {
    "cycle_time": {
      "enabled": false
    },
    "scan_timeout": {
      "enabled": true,
      "value_ms": 100
    },
    "agc_interval_seconds": 30
  },
  "files": {},
  "precursor_selection": {
    "rt_window": 180,
    "targeting": "none",
    "rank_by": "qscore",
    "max_precursors": 10
  },
  "characterization": {
    "mode": "off",
    "max_targets": 10
  },
  "ms_settings": {
    "ms1": {
      "analyzer": "Orbitrap",
      "first_mass": 500,
      "last_mass": 2000,
      "resolution": 120000,
      "agc_target": 800000,
      "max_it": 246
    },
    "ms2": {
      "analyzer": "Orbitrap",
      "activation": "HCD",
      "collision_energy": 29,
      "resolution": 120000
    }
  },
  "tagging": {},
  "quantification": {
    "enabled": false
  }
}
)";

  // ---- Survey re-planning: plan units built by hand -------------------------------------------------

  /// The minimal config with scan_timeout at @p ms, for the expiry section.
  std::string configWithTimeout(int ms)
  {
    std::string s(minimal_config);
    const std::string from = "\"value_ms\": 100";
    s.replace(s.find(from), from.size(), "\"value_ms\": " + std::to_string(ms));
    return s;
  }

  /// An MS2 at @p charge isolating @p mz, parented to nothing (a survey's MS2s carry the survey's id;
  /// the queue reads the parent only for the causal invariant, which none of these sections exercises).
  ScanCommand ms2(int scan_id, int charge, double mz = 951.0, int priority = 2)
  {
    ScanCommand cmd{};
    cmd.scan_id = scan_id;
    cmd.msn_level = 2;
    cmd.priority = priority;
    cmd.num_stages = 1;
    cmd.stages[0].charge_state = charge;
    cmd.stages[0].precursor_mz = mz;
    cmd.qscore = 0.9;
    return cmd;
  }

  MemoryStamp stamp(MemoryStamp::Kind kind, int key, int charge = 0, double value = 60.0)
  {
    MemoryStamp s;
    s.kind = kind;
    s.key = key;
    s.charge = charge;
    s.value = value;
    return s;
  }

  /// One species of one survey: a command per charge, a species-level TqMassRt stamp and a per-charge
  /// AllMassRt stamp keyed on the charge, so which stamps were committed is readable off their charges.
  PlanUnit unit(int nominal, const std::vector<ScanCommand>& cmds, double cv, int decision_survey_id = 1)
  {
    PlanUnit u;
    u.faims_cv = cv;
    u.nominal_mass = nominal;
    u.decision_survey_id = decision_survey_id;
    u.commands = cmds;
    for (const auto& c : cmds)
    {
      u.purity.push_back(0.5);
      u.stamps.per_charge[c.stages[0].charge_state].push_back(stamp(MemoryStamp::Kind::AllMassRt, nominal, c.stages[0].charge_state));
    }
    u.stamps.nominal_mass = nominal;
    u.stamps.species_level.push_back(stamp(MemoryStamp::Kind::TqMassRt, nominal));
    return u;
  }

  std::vector<int> idsOf(const PlanUnit& u)
  {
    std::vector<int> ids;
    for (const auto& c : u.commands) ids.push_back(c.scan_id);
    return ids;
  }

  /// Dequeue until empty; the scan ids in dequeue order.
  std::vector<int> drainIds(ScanCommandQueue& q)
  {
    std::vector<int> ids;
    while (auto c = q.dequeue()) ids.push_back(c->scan_id);
    return ids;
  }

  /// (kind, key, charge) of each stamp, so a committed sequence compares as a value.
  std::vector<std::string> signature(const std::vector<MemoryStamp>& stamps)
  {
    std::vector<std::string> out;
    for (const auto& s : stamps)
    {
      out.push_back(std::to_string(static_cast<int>(s.kind)) + ":" + std::to_string(s.key) + ":" + std::to_string(s.charge));
    }
    return out;
  }

  uint64_t steadyNowMs()
  {
    return static_cast<uint64_t>(
      std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now().time_since_epoch()).count());
  }
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
    ScanCommand cmd = queue.buildMS2(pg, 10, sc, 2, 0);
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

// ---- Survey re-planning: the plan registry -----------------------------------------------------------
//
// A survey hands the queue its whole plan for one CV (replacePlan). A plan unit stays UNCOMMITTED, and so
// droppable by the next survey at that CV, until the drain dequeues one of its commands; the dequeue commits
// it and moves its stamps into the inbox the analysis thread folds (takeCommitted).

START_SECTION(replacePlan_drops_uncommitted_units_at_its_cv_only)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  const PlanUnit u1 = unit(12351, {ms2(1, 13), ms2(2, 14)}, -40.0);
  const PlanUnit u2 = unit(13000, {ms2(3, 13, 1001.0)}, -60.0);
  const PlanUnit u3 = unit(14000, {ms2(4, 13, 1078.0)}, -40.0);

  q.replacePlan(-40.0, {u1}, {}, {});
  q.replacePlan(-60.0, {u2}, {}, {});
  const ReplanResult r = q.replacePlan(-40.0, {u3}, {}, {});

  std::set<int> dropped(r.dropped_scan_ids.begin(), r.dropped_scan_ids.end());
  TEST_EQUAL(dropped == std::set<int>({1, 2}), true)
  TEST_EQUAL(r.added, 1)
  const std::vector<int> drained = drainIds(q);
  TEST_EQUAL(std::set<int>(drained.begin(), drained.end()) == std::set<int>({3, 4}), true)
  TEST_EQUAL(static_cast<int>(drained.size()), 2)
}
END_SECTION

START_SECTION(dequeue_commits_the_unit_and_keeps_it)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  const PlanUnit u1 = unit(12351, {ms2(11, 13), ms2(12, 14)}, 0.0);
  q.replacePlan(0.0, {u1}, {}, {});

  const auto first = q.dequeue();
  ABORT_IF(!first.has_value())
  TEST_EQUAL(first->scan_id, 11)

  // The species is committed: a newer survey at its CV that no longer plans it does not drop its remainder.
  const ReplanResult r = q.replacePlan(0.0, {}, {}, {});
  TEST_EQUAL(r.dropped_scan_ids.empty(), true)

  const auto second = q.dequeue();
  ABORT_IF(!second.has_value())
  TEST_EQUAL(second->scan_id, 12)

  std::vector<MemoryStamp> expected = u1.stamps.species_level;
  expected.insert(expected.end(), u1.stamps.per_charge.at(13).begin(), u1.stamps.per_charge.at(13).end());
  expected.insert(expected.end(), u1.stamps.per_charge.at(14).begin(), u1.stamps.per_charge.at(14).end());
  TEST_EQUAL(signature(q.takeCommitted()) == signature(expected), true)
}
END_SECTION

START_SECTION(per_charge_stamps_commit_with_their_own_command)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  PlanUnit u = unit(12351, {ms2(21, 13), ms2(22, 14), ms2(23, 16)}, 0.0);
  u.separate = true;
  q.replacePlan(0.0, {u}, {}, {});

  const auto first = q.dequeue();
  ABORT_IF(!first.has_value())
  TEST_EQUAL(first->scan_id, 21)

  // The species-level stamps and charge 13's own; 14 and 16 commit only when their commands go out.
  std::vector<MemoryStamp> expected = u.stamps.species_level;
  expected.insert(expected.end(), u.stamps.per_charge.at(13).begin(), u.stamps.per_charge.at(13).end());
  TEST_EQUAL(signature(q.takeCommitted()) == signature(expected), true)
}
END_SECTION

START_SECTION(reconcile_refuses_a_unit_committed_since_the_fold)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  // U1 commits between the newer survey's plan and its swap; the analysis thread has not folded it yet,
  // so that survey selected A again. The reconcile step reads the commit and refuses the new unit.
  q.replacePlan(0.0, {unit(12351, {ms2(31, 13)}, 0.0)}, {}, {});
  ABORT_IF(!q.dequeue().has_value())

  PlanUnit u2 = unit(12351, {ms2(32, 13)}, 0.0, 2);
  u2.unit_id = 77;
  const ReplanResult r = q.replacePlan(0.0, {u2}, {}, {});

  TEST_EQUAL(r.reconciled_unit_ids == std::vector<int>({77}), true)
  TEST_EQUAL(r.added, 0)
  TEST_EQUAL(drainIds(q).empty(), true)   // U2's command never reached the queue
}
END_SECTION

START_SECTION(unregistered_commands_are_never_dropped)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  // A follow-up and a monitor MS1 belong to no survey's plan, so no replacePlan may touch them.
  ScanCommand follow_up = ms2(41, 13, 951.0, 0);
  ScanCommand monitor{};
  monitor.scan_id = 42;
  monitor.msn_level = 1;
  monitor.priority = 0;
  q.push(follow_up);
  q.push(monitor);

  const ReplanResult r = q.replacePlan(0.0, {}, {}, {});
  TEST_EQUAL(r.dropped_scan_ids.empty(), true)
  TEST_EQUAL(drainIds(q) == std::vector<int>({41, 42}), true)
}
END_SECTION

START_SECTION(adjustable_is_rewritten_never_dropped)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  q.push(ms2(51, 13, 951.0, 0));   // the bought R
  q.registerAdjustable(51, {13, 14}, 5);
  std::this_thread::sleep_for(std::chrono::milliseconds(5));
  const uint64_t before_adjust = steadyNowMs();

  const ReplanResult r = q.replacePlan(0.0, {}, {{51, ms2(51, 14, 882.0, 0), 9, 0.7}}, {});
  TEST_EQUAL(r.adjusted, 1)
  TEST_EQUAL(r.dropped_scan_ids.empty(), true)

  const ReplanResult r2 = q.replacePlan(0.0, {}, {}, {});
  TEST_EQUAL(r2.dropped_scan_ids.empty(), true)

  const auto cmd = q.dequeue();
  ABORT_IF(!cmd.has_value())
  TEST_EQUAL(cmd->scan_id, 51)
  TEST_EQUAL(cmd->stages[0].charge_state, 14)
  TEST_EQUAL(cmd->enqueue_timestamp_ms >= before_adjust, true)   // reset by the adjust, not left at the push

  const auto info = q.takeDecisionInfo(51);
  ABORT_IF(!info.has_value())
  TEST_EQUAL(info->decision_survey_id, 9)
  TEST_REAL_SIMILAR(info->isolation_purity, 0.7)
  TEST_EQUAL(q.takeDecisionInfo(51).has_value(), false)   // a take, not a peek
}
END_SECTION

START_SECTION(expired_member_leaves_registry)
{
  Config cfg{configWithTimeout(1)};
  ScanCommandQueue q(cfg);

  q.replacePlan(0.0, {unit(12351, {ms2(61, 13), ms2(62, 14)}, 0.0)}, {}, {});
  std::this_thread::sleep_for(std::chrono::milliseconds(5));
  q.cleanupExpired();

  // Never acquired, so never excluded: an uncommitted unit that empties is discarded with its stamps.
  TEST_EQUAL(q.takeCommitted().empty(), true)
  TEST_EQUAL(FLASHIdaTestAccess::planUnitCount(q), 0)
  TEST_EQUAL(drainIds(q).empty(), true)
}
END_SECTION

START_SECTION(enqueue_timestamp_reset_on_insert_and_adjust)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  // A stale timestamp on the way in is overwritten on insert ...
  ScanCommand planned = ms2(71, 13);
  planned.enqueue_timestamp_ms = 1;
  const uint64_t before_insert = steadyNowMs();
  q.replacePlan(0.0, {unit(12351, {planned}, 0.0)}, {}, {});
  const auto inserted = q.dequeue();
  ABORT_IF(!inserted.has_value())
  TEST_EQUAL(inserted->enqueue_timestamp_ms >= before_insert, true)

  // ... and on adjust, where the rewritten command carries one too.
  q.push(ms2(72, 13, 951.0, 0));
  q.registerAdjustable(72, {13, 14}, 5);
  std::this_thread::sleep_for(std::chrono::milliseconds(5));
  ScanCommand rewritten = ms2(72, 14, 882.0, 0);
  rewritten.enqueue_timestamp_ms = 1;
  const uint64_t before_adjust = steadyNowMs();
  q.replacePlan(0.0, {}, {{72, rewritten, 9, 0.7}}, {});
  const auto adjusted = q.dequeue();
  ABORT_IF(!adjusted.has_value())
  TEST_EQUAL(adjusted->enqueue_timestamp_ms >= before_adjust, true)
}
END_SECTION

// A multiplexed command isolates its anchor AND its stage-0 notches, so its dequeue commits the per-charge
// stamps of every charge it isolates -- under an authored charge set those are what keep the notch charges
// from being selected again (ADR-0028).
START_SECTION(multiplexed_command_commits_every_notch_charge)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  ScanCommand cmd = ms2(81, 13);
  cmd.stage0_notch_count = 2;
  cmd.notches[0].charge_state = 14;
  cmd.notches[1].charge_state = 16;
  PlanUnit u = unit(12351, {cmd}, 0.0);
  u.stamps.per_charge[14].push_back(stamp(MemoryStamp::Kind::AuthoredAcquired, 12351, 14));
  u.stamps.per_charge[16].push_back(stamp(MemoryStamp::Kind::AuthoredAcquired, 12351, 16));
  q.replacePlan(0.0, {u}, {}, {});
  ABORT_IF(!q.dequeue().has_value())

  std::vector<MemoryStamp> expected = u.stamps.species_level;
  for (int z : {13, 14, 16}) expected.insert(expected.end(), u.stamps.per_charge.at(z).begin(), u.stamps.per_charge.at(z).end());
  TEST_EQUAL(signature(q.takeCommitted()) == signature(expected), true)
}
END_SECTION

// A cancelled member leaves the registry exactly as an expired one does: its unit empties when the rest is sent,
// and nothing derived from it is mistaken for the child of a command that never went out.
START_SECTION(cancelled_member_leaves_registry)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  q.replacePlan(0.0, {unit(12351, {ms2(91, 13), ms2(92, 14)}, 0.0)}, {}, {});
  TEST_EQUAL(q.cancelByScanIds({91}) == std::vector<int>({91}), true)
  const auto second = q.dequeue();
  ABORT_IF(!second.has_value())
  TEST_EQUAL(second->scan_id, 92)
  TEST_EQUAL(FLASHIdaTestAccess::planUnitCount(q), 0)

  ScanCommand child = ms2(93, 13, 951.0, 0);
  const std::string parent = ScanCommandQueue::encode(91);
  std::copy(parent.begin(), parent.end(), child.parent_scan_id);
  std::ostringstream captured;
  std::streambuf* const real_cout = std::cout.rdbuf(captured.rdbuf());
  q.push(child);
  std::cout.rdbuf(real_cout);
  TEST_EQUAL(captured.str().find("[REPLAN-INVARIANT]") == std::string::npos, true)
}
END_SECTION

// A reconciled unit is reported with everything its caller must undo: its commands (never queued, so their
// tracking is forgotten) and the quant and sweep groups planned for it.
START_SECTION(reconciled_unit_is_reported_in_full)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  q.replacePlan(0.0, {unit(12351, {ms2(101, 13)}, 0.0)}, {}, {});
  ABORT_IF(!q.dequeue().has_value())

  PlanUnit u2 = unit(12351, {ms2(102, 13), ms2(103, 14)}, 0.0, 2);
  u2.unit_id = 78;
  u2.quant_group_ids = {5};
  u2.sweep_group_ids = {7};
  const ReplanResult r = q.replacePlan(0.0, {u2}, {}, {});

  TEST_EQUAL(r.reconciled_unit_ids == std::vector<int>({78}), true)
  TEST_EQUAL(r.reconciled_scan_ids == idsOf(u2), true)
  TEST_EQUAL(r.reconciled_quant_groups == std::vector<int>({5}), true)
  TEST_EQUAL(r.reconciled_sweep_groups == std::vector<int>({7}), true)
}
END_SECTION

// Review focus 1: a `separate` remainder is re-planned from a waitingAt snapshot, and the drain may send one of its
// charges before the swap. The replacement inherits what was sent AT THE SWAP, so that charge is not sent twice.
START_SECTION(remainder_replacement_skips_a_charge_sent_after_waitingAt)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  PlanUnit u1 = unit(12351, {ms2(111, 13), ms2(112, 14), ms2(113, 16)}, 0.0);
  u1.separate = true;
  q.replacePlan(0.0, {u1}, {}, {});
  std::vector<int> sent_charges;
  const auto first = q.dequeue();   // z13: the species commits
  ABORT_IF(!first.has_value())
  sent_charges.push_back(first->stages[0].charge_state);

  const ScanCommandQueue::Waiting w = q.waitingAt(0.0);
  ABORT_IF(w.committed_units.size() != 1)
  TEST_EQUAL(idsOf(w.committed_units[0]) == std::vector<int>({112, 113}), true)

  const auto raced = q.dequeue();   // z14 goes out between the snapshot and the swap
  ABORT_IF(!raced.has_value())
  sent_charges.push_back(raced->stages[0].charge_state);

  PlanUnit remainder = unit(12351, {ms2(114, 14), ms2(115, 17)}, 0.0, 2);
  remainder.separate = true;
  remainder.committed = true;
  const ReplanResult r = q.replacePlan(0.0, {remainder}, {}, {w.committed_units[0].unit_id});

  TEST_EQUAL(r.dropped_scan_ids == std::vector<int>({113}), true)
  TEST_EQUAL(r.reconciled_scan_ids == std::vector<int>({114}), true)
  std::vector<int> rest;
  while (auto c = q.dequeue())
  {
    rest.push_back(c->scan_id);
    sent_charges.push_back(c->stages[0].charge_state);
  }
  TEST_EQUAL(rest == std::vector<int>({115}), true)
  TEST_EQUAL(std::set<int>(sent_charges.begin(), sent_charges.end()).size(), sent_charges.size())
}
END_SECTION

// ...and when the drain sent the whole remainder before the swap, the named unit is gone: there is nothing left to
// replace, so the replacement is refused rather than acquiring the species past its plan.
START_SECTION(remainder_replacement_refused_when_old_unit_emptied)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  PlanUnit u1 = unit(12351, {ms2(121, 13), ms2(122, 14), ms2(123, 16)}, 0.0);
  u1.separate = true;
  q.replacePlan(0.0, {u1}, {}, {});
  ABORT_IF(!q.dequeue().has_value())

  const ScanCommandQueue::Waiting w = q.waitingAt(0.0);
  ABORT_IF(w.committed_units.size() != 1)
  ABORT_IF(!q.dequeue().has_value())
  ABORT_IF(!q.dequeue().has_value())

  PlanUnit remainder = unit(12351, {ms2(124, 14), ms2(125, 17)}, 0.0, 2);
  remainder.unit_id = 88;
  remainder.separate = true;
  remainder.committed = true;
  const ReplanResult r = q.replacePlan(0.0, {remainder}, {}, {w.committed_units[0].unit_id});

  TEST_EQUAL(r.reconciled_unit_ids == std::vector<int>({88}), true)
  TEST_EQUAL(r.reconciled_scan_ids == idsOf(remainder), true)
  TEST_EQUAL(r.added, 0)
  TEST_EQUAL(drainIds(q).empty(), true)
  TEST_EQUAL(FLASHIdaTestAccess::planUnitCount(q), 0)
}
END_SECTION

// The race at queue level: four drains dequeue while one analysis thread re-plans and pushes unregistered
// work. Nothing may be dequeued twice, nothing unregistered lost, and a unit is either sent (from its first
// dequeue on, never dropped) or dropped (none of it ever sent).
START_SECTION(replacePlan_vs_concurrent_dequeue_conservation)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  const int N_DRAINS = 4;
  std::atomic<bool> stop{false};
  std::vector<std::vector<int>> dequeued(N_DRAINS + 1);
  auto drain = [&](int t)
  {
    while (!stop.load(std::memory_order_acquire))
    {
      if (auto c = q.dequeue()) dequeued[t].push_back(c->scan_id);
      else std::this_thread::yield();
    }
  };

  std::set<int> unregistered, dropped, reconciled_ids, registered;
  std::map<int, int> unit_of;   // scan id -> unit id
  std::set<int> dropped_units;
  auto replanner = [&]()
  {
    std::mt19937 rng(20261022);
    int next_id = 1000, next_unit = 1;
    for (int i = 0; i < 500; ++i)
    {
      const double cv = (rng() % 2 == 0) ? -40.0 : -60.0;
      std::vector<PlanUnit> units;
      const int n_units = static_cast<int>(rng() % 4);
      for (int k = 0; k < n_units; ++k)
      {
        const int nominal = 10000 + static_cast<int>(rng() % 20);
        std::vector<ScanCommand> cmds;
        const int n_cmds = 1 + static_cast<int>(rng() % 3);
        for (int c = 0; c < n_cmds; ++c) cmds.push_back(ms2(next_id++, 10 + c, 900.0 + (nominal - 10000)));
        PlanUnit u = unit(nominal, cmds, cv);
        u.unit_id = next_unit++;
        for (int id : idsOf(u))
        {
          unit_of[id] = u.unit_id;
          registered.insert(id);
        }
        units.push_back(u);
      }
      std::map<int, std::vector<int>> ids_of_unit;
      for (const auto& u : units) ids_of_unit[u.unit_id] = idsOf(u);

      const ReplanResult r = q.replacePlan(cv, std::move(units), {}, {});
      for (int id : r.dropped_scan_ids)
      {
        dropped.insert(id);
        dropped_units.insert(unit_of[id]);
      }
      for (int uid : r.reconciled_unit_ids)
      {
        for (int id : ids_of_unit[uid]) reconciled_ids.insert(id);
      }

      if (rng() % 2 == 0)
      {
        unregistered.insert(next_id);
        q.push(ms2(next_id++, 13, 951.0, 0));
      }
      if (i % 50 == 49) q.takeCommitted();   // the analysis thread's fold
    }
  };

  std::vector<std::thread> threads;
  for (int t = 0; t < N_DRAINS; ++t) threads.emplace_back(drain, t);
  std::thread planner(replanner);
  planner.join();
  stop.store(true, std::memory_order_release);
  for (auto& th : threads) th.join();
  dequeued[N_DRAINS] = drainIds(q);

  std::map<int, int> times;
  for (const auto& v : dequeued)
    for (int id : v) ++times[id];

  bool none_twice = true, unregistered_once = true, dropped_never_sent = true, dropped_unit_never_sent = true,
       reconciled_never_sent = true, conserved = true;
  for (const auto& [id, n] : times)
  {
    if (n != 1) none_twice = false;
    if (dropped.count(id)) dropped_never_sent = false;
    if (reconciled_ids.count(id)) reconciled_never_sent = false;
    if (registered.count(id) && dropped_units.count(unit_of[id])) dropped_unit_never_sent = false;
  }
  for (int id : unregistered)
    if (times[id] != 1) unregistered_once = false;
  for (int id : registered)   // every planned command is sent, dropped or was never inserted -- exactly one
  {
    const int fates = (times.count(id) ? 1 : 0) + (dropped.count(id) ? 1 : 0) + (reconciled_ids.count(id) ? 1 : 0);
    if (fates != 1) conserved = false;
  }

  TEST_EQUAL(none_twice, true)
  TEST_EQUAL(unregistered_once, true)
  TEST_EQUAL(dropped_never_sent, true)
  TEST_EQUAL(dropped_unit_never_sent, true)
  TEST_EQUAL(reconciled_never_sent, true)
  TEST_EQUAL(conserved, true)
  TEST_EQUAL(FLASHIdaTestAccess::planUnitCount(q), 0)   // every unit sent in full or dropped: none left behind
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
