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
#include "FLASHIda_EventCapture.h"

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

// ADR-0061 d.9 / ADR-0062 d.3: a close-out withdraws the species' wholly unsent plan units at EVERY CV in one
// queue_mutex_ step and puts its stamps in the inbox. A committed unit of the species (a scan already sent) keeps
// its remaining commands, and another species is untouched.
START_SECTION(close_out_withdraws_unsent_units_at_every_cv)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  q.replacePlan(-40.0, {unit(12351, {ms2(51, 13), ms2(52, 14)}, -40.0)}, {}, {});   // committed below
  ABORT_IF(!q.dequeue().has_value())                                                  // sends 51: the unit commits
  (void)q.takeCommitted();                                                            // the analysis thread folded it
  q.replacePlan(-60.0, {unit(12351, {ms2(53, 13)}, -60.0)}, {}, {});                   // unsent, other CV
  q.replacePlan(-80.0, {unit(12351, {ms2(54, 13)}, -80.0), unit(13000, {ms2(55, 13, 1001.0)}, -80.0)}, {}, {});

  const std::vector<MemoryStamp> stamps = {stamp(MemoryStamp::Kind::CloseOutMassRt, 12351, 0, 70.0),
                                           stamp(MemoryStamp::Kind::CloseOutMzRt, 951, 0, 70.0)};
  const ReplanResult r = q.closeOut(12351, stamps);

  std::set<int> dropped(r.dropped_scan_ids.begin(), r.dropped_scan_ids.end());
  TEST_EQUAL(dropped == std::set<int>({53, 54}), true)   // both CVs' unsent units
  TEST_EQUAL(r.kept, 1)                                   // the committed one
  TEST_EQUAL(signature(q.takeCommitted()) == signature(stamps), true)
  const std::vector<int> drained = drainIds(q);
  TEST_EQUAL(std::set<int>(drained.begin(), drained.end()) == std::set<int>({52, 55}), true)   // its remainder; another species
  TEST_EQUAL(static_cast<int>(drained.size()), 2)
}
END_SECTION

// The close-out's stamps also reach committed_since_fold_: a survey that planned the species before the analysis
// thread folded the close-out has its fresh unit reconciled away, the same way a commit in that gap is.
START_SECTION(close_out_stamps_reconcile_a_racing_plan)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  q.closeOut(12351, {stamp(MemoryStamp::Kind::CloseOutMassRt, 12351, 0, 70.0)});
  PlanUnit u = unit(12351, {ms2(61, 13, 1100.0)}, 0.0, 2);   // another isolation m/z: the mass bar alone refuses it
  u.unit_id = 88;
  const ReplanResult r = q.replacePlan(0.0, {u, unit(13000, {ms2(62, 13, 1001.0)}, 0.0, 2)}, {}, {});

  TEST_EQUAL(r.reconciled_unit_ids == std::vector<int>({88}), true)
  TEST_EQUAL(r.added, 1)
  TEST_EQUAL(drainIds(q) == std::vector<int>({62}), true)
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
  EventCapture cap;
  q.push(child);
  TEST_EQUAL(cap.str().find("[REPLAN-INVARIANT]") == std::string::npos, true)
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

// reisolate is buildMS2's geometry applied to an existing command: whatever the command isolated before, the
// re-isolated one carries the same stage-0 window, notches and scan_description suffix as buildMS2 builds for
// the same PeakGroup and charge set -- and keeps its own scan id. Multiplexed, so the notch write is exercised.
START_SECTION(reisolate_matches_buildMS2_geometry)
{
  std::string multiplexed(minimal_config);
  const std::string from = "\"rank_by\": \"qscore\",";
  multiplexed.replace(multiplexed.find(from), from.size(), from + " \"precursor_charges\": \"multiplexed\",");
  Config cfg{multiplexed};
  ScanCommandQueue q(cfg);

  // One species at charges 12..16, three isotopes each, so every charge has a measured m/z range.
  const double mass = 12351.3933;
  PeakGroup pg(12, 16, true);
  pg.setMonoisotopicMass(mass);
  for (int z = 12; z <= 16; ++z)
  {
    for (int iso = 0; iso < 3; ++iso)
    {
      FLASHHelperClasses::LogMzPeak lp;
      lp.mz = (mass + iso * 1.003355) / z + 1.007276;
      lp.abs_charge = z;
      lp.isotopeIndex = iso;
      lp.intensity = 1000.0f * z;
      pg.push_back(lp);
    }
  }
  ScanConfig sc;
  sc.collision_energy = 29;
  sc.activation = "HCD";

  // The command to re-isolate: another charge with another notch set, as a bought R a survey ago would be.
  const std::vector<int> old_set = {13, 15};
  const ScanCommand r = q.buildMS2(pg, 13, sc, 0, 0, &old_set, 'R');
  ABORT_IF(r.stage0_notch_count != 1)

  for (const std::vector<int>& set : std::vector<std::vector<int>>{{14, 12, 16}, {16}})
  {
    const ScanCommand expected = q.buildMS2(pg, set.front(), sc, 0, 0, &set, 'R');
    std::vector<std::pair<int, const PeakGroup*>> picks;
    for (int z : set) picks.emplace_back(z, &pg);
    const ScanCommand got = ScanCommandQueue::reisolate(r, picks, optimal_window_margin_, 'R');

    TEST_EQUAL(got.scan_id, r.scan_id)
    TEST_EQUAL(got.stages[0].charge_state, expected.stages[0].charge_state)
    TEST_REAL_SIMILAR(got.stages[0].precursor_mz, expected.stages[0].precursor_mz)
    TEST_REAL_SIMILAR(got.stages[0].isolation_width, expected.stages[0].isolation_width)
    TEST_EQUAL(got.stage0_notch_count, expected.stage0_notch_count)
    TEST_EQUAL(got.stage0_notch_count, static_cast<int>(set.size()) - 1)
    for (int j = 0; j < MAX_NOTCHES_PER_STAGE; ++j) // the used slots match, and no stale slot survives
    {
      TEST_EQUAL(got.notches[j].charge_state, expected.notches[j].charge_state)
      TEST_EQUAL(got.notches[j].precursor_mz == expected.notches[j].precursor_mz, true)
      TEST_EQUAL(got.notches[j].isolation_width == expected.notches[j].isolation_width, true)
    }
    TEST_EQUAL(std::string(got.scan_description).substr(3), std::string(expected.scan_description).substr(3))
    TEST_EQUAL(std::string(got.scan_description).substr(0, 3), ScanCommandQueue::encode(r.scan_id))
    TEST_REAL_SIMILAR(got.mono_mass, expected.mono_mass)
    TEST_REAL_SIMILAR(got.precursor_intensity, expected.precursor_intensity)
  }
}
END_SECTION

// A bought R may move only among the charges its group's Q scans covered (owner 2026-10-23). With none, there
// is nothing it may move to: the registration is refused, so no survey can re-choose it -- an empty set must
// never read as "no restriction" -- and the R goes out exactly as bought.
START_SECTION(bought_r_empty_covered_set_is_not_adjustable)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  q.push(ms2(91, 13, 951.0, 0));   // the bought R
  q.registerAdjustable(91, {}, 5, 12351);
  TEST_EQUAL(FLASHIdaTestAccess::adjustableResolvedCharges(q, 91).has_value(), false)
  TEST_EQUAL(q.waitingAt(0.0).adjustable.empty(), true)

  const ReplanResult r = q.replacePlan(0.0, {}, {{91, ms2(91, 14, 882.0, 0), 9, 0.7}}, {});
  TEST_EQUAL(r.adjusted, 0)
  const auto cmd = q.dequeue();
  ABORT_IF(!cmd.has_value())
  TEST_EQUAL(cmd->scan_id, 91)
  TEST_EQUAL(cmd->stages[0].charge_state, 13)
  TEST_REAL_SIMILAR(cmd->stages[0].precursor_mz, 951.0)
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

// T10: an in-flight command whose scan never comes back is dropped once MAX_IN_FLIGHT later commands
// have been sent -- and not before. Until this, pending_scan_map_ was emptied only by an answer or a
// cancel, so every unanswered command (~2 KB) stayed for the life of the engine.
//
// The drive mimics a host: each later command is dequeued (which registers it) and answered at once.
// The sweep runs every IN_FLIGHT_SWEEP_EVERY registrations, so the lost command (registration 1) must
// survive every sweep up to the first one at which it is MORE than MAX_IN_FLIGHT registrations old.
START_SECTION(an_unanswered_command_is_dropped_after_max_in_flight_later_commands)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  const uint64_t max_in_flight = ScanCommandQueue::MAX_IN_FLIGHT;
  const uint64_t every = ScanCommandQueue::IN_FLIGHT_SWEEP_EVERY;
  // The first sweep at which registration 1 is more than max_in_flight registrations old.
  const uint64_t deciding = ((max_in_flight + 1) / every + 1) * every;
  TEST_EQUAL(deciding > max_in_flight + 1 && deciding % every == 0 && deciding - every <= max_in_flight + 1, true)

  ScanCommand lost{};
  lost.msn_level = 2;
  lost.scan_id = queue.nextTrackingId();
  queue.registerPending(lost);                                   // registration 1, never answered

  int recent_orphan = -1;
  bool lost_survived_until_deciding = true;
  for (uint64_t r = 2; r <= deciding; ++r)
  {
    ScanCommand cmd{};
    cmd.msn_level = 2;
    cmd.priority = 2;
    cmd.scan_id = queue.nextTrackingId();
    queue.push(cmd);
    auto sent = queue.dequeue();                                 // registration r
    TEST_EQUAL(sent.has_value() && sent->scan_id == cmd.scan_id, true)
    // One more command that never comes back, but recent: it must outlive the deciding sweep.
    if (r == deciding - 100)
    {
      recent_orphan = cmd.scan_id;
      continue;
    }
    TEST_EQUAL(queue.resolvePending(cmd.scan_id).has_value(), true)
    if (r < deciding && !queue.peekPending(lost.scan_id).has_value()) lost_survived_until_deciding = false;
  }

  TEST_EQUAL(lost_survived_until_deciding, true)                 // not dropped early
  TEST_EQUAL(queue.peekPending(lost.scan_id).has_value(), false) // dropped at the deciding sweep
  TEST_EQUAL(queue.resolvePending(lost.scan_id).has_value(), false)
  TEST_EQUAL(queue.peekPending(recent_orphan).has_value(), true) // a recent one is still awaited
  TEST_EQUAL(static_cast<int>(queue.pendingScanMapSize()), 1)    // and nothing else is left
}
END_SECTION

// T11: the age is counted from an entry's LATEST registration. Registering an id again (it is
// re-sent) restarts its count, so the sweep must not drop it on the strength of the first one.
START_SECTION(a_re_registered_command_is_aged_from_its_latest_registration)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue queue(cfg);

  const uint64_t max_in_flight = ScanCommandQueue::MAX_IN_FLIGHT;
  const uint64_t every = ScanCommandQueue::IN_FLIGHT_SWEEP_EVERY;
  const uint64_t deciding = ((max_in_flight + 1) / every + 1) * every;

  ScanCommand resent{};
  resent.msn_level = 2;
  resent.scan_id = queue.nextTrackingId();
  queue.registerPending(resent);                                 // registration 1

  for (uint64_t r = 2; r <= deciding; ++r)
  {
    if (r == deciding - 10)
    {
      queue.registerPending(resent);                             // registration r: sent again
      continue;
    }
    ScanCommand cmd{};
    cmd.msn_level = 2;
    cmd.scan_id = queue.nextTrackingId();
    queue.registerPending(cmd);
    queue.resolvePending(cmd.scan_id);
  }

  TEST_EQUAL(queue.peekPending(resent.scan_id).has_value(), true)
  TEST_EQUAL(static_cast<int>(queue.pendingScanMapSize()), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
// ===== Per-proteoform MS3 plans (ADR-0064) =========================================================================
namespace
{
  Ms3ProteoformKey pkey(double mass, const std::string& acc = "P1", int start = 0, int end = 104)
  {
    Ms3ProteoformKey k;
    k.accession = acc;
    k.region_start = start;
    k.region_end = end;
    k.commanded_mass = mass;
    return k;
  }

  /// One single-command item: an MS3 (priority 1) on cleavage (@p type, @p index) at @p charge, chosen at @p intensity.
  Ms3PlanItem ms3Item(int scan_id, const std::string& type, int index, int charge, double intensity, int owner = 1)
  {
    Ms3PlanItem it;
    it.cleavage = {type, index};
    it.charges = {{charge, intensity}};
    ScanCommand c = ms2(scan_id, charge, 900.0 + index, 1);
    c.msn_level = 3;
    it.commands = {c};
    it.owner_precursor_id = owner;
    return it;
  }

  std::set<int> asSet(const std::vector<int>& v) { return std::set<int>(v.begin(), v.end()); }
}

// The pure rules: one proteoform across selections is the same protein and region at the same commanded mass, within
// the MS1 tolerance or 1-3 isotope steps / a charge harmonic (sameSpeciesMass); a sent fragment needs 2x.
START_SECTION(ms3_registry_rules_same_proteoform_and_admits)
{
  const double tol = 10.0;
  TEST_EQUAL(sameProteoform(pkey(12351.39), pkey(12351.40), tol), true)              // within 10 ppm
  TEST_EQUAL(sameProteoform(pkey(12351.39), pkey(12351.39 + 1.00235), tol), true)    // one isotope step
  TEST_EQUAL(sameProteoform(pkey(12351.39), pkey(12351.39 * 2.0), tol), true)        // a harmonic
  TEST_EQUAL(sameProteoform(pkey(12351.39), pkey(12351.39 + 57.02146), tol), false)  // another species
  TEST_EQUAL(sameProteoform(pkey(12351.39), pkey(12351.39, "P2"), tol), false)       // another protein
  TEST_EQUAL(sameProteoform(pkey(12351.39), pkey(12351.39, "P1", 1, 104), tol), false) // another region

  Ms3RegistrySnapshot snap;
  snap.sent[{"b", 50}][8] = 1000.0;
  TEST_EQUAL(ms3Admits(snap, {"y", 29}, 4, 1.0), true)        // never sent
  TEST_EQUAL(ms3Admits(snap, {"b", 50}, 8, 1999.0), false)    // same charge, 1.999x
  TEST_EQUAL(ms3Admits(snap, {"b", 50}, 8, 2000.0), true)     // same charge, 2x
  TEST_EQUAL(ms3Admits(snap, {"b", 50}, 6, 1500.0), false)    // another charge of a sent cleavage: vs its max
  TEST_EQUAL(ms3Admits(snap, {"b", 50}, 6, 2000.0), true)
  snap.sent[{"b", 50}][6] = 3000.0;
  TEST_EQUAL(ms3Admits(snap, {"b", 50}, 7, 4000.0), false)    // a third charge vs the cleavage's highest (3000)
  TEST_EQUAL(ms3NewCleavageRoom(snap, 3), 2)
}
END_SECTION

// A newer plan replaces the still-unsent items of the same proteoform (matched across a one-isotope mass step) and
// keeps a committed one; a dequeue commits the item's fragment per charge at its chosen intensity; another
// proteoform's entry is untouched; unregistered MS3s are never touched.
START_SECTION(ms3_plan_replaces_unsent_and_commits_on_dequeue)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  ScanCommand loose = ms2(70, 4, 700.0, 1);   // an unregistered MS3-priority command
  loose.msn_level = 3;
  q.push(loose);

  Ms3ReplanResult r1 = q.replaceMs3Plan(pkey(12351.39), {ms3Item(71, "b", 50, 8, 1000.0), ms3Item(72, "y", 29, 4, 800.0)}, 0, 3);
  TEST_EQUAL(asSet(r1.inserted_scan_ids) == std::set<int>({71, 72}), true)
  Ms3ReplanResult other = q.replaceMs3Plan(pkey(8000.0, "P2"), {ms3Item(80, "b", 10, 2, 100.0)}, 0, 3);
  TEST_EQUAL(other.entry_id != r1.entry_id, true)

  // The drain sends 70 (unregistered), then 71: b50 z8 commits at 1000.
  auto d1 = q.dequeue();
  auto d2 = q.dequeue();
  ABORT_IF(!d1.has_value() || !d2.has_value())
  TEST_EQUAL(d1->scan_id, 70)
  TEST_EQUAL(d2->scan_id, 71)
  const Ms3RegistrySnapshot snap = q.ms3Snapshot(pkey(12351.39 + 1.00235));   // one isotope up: the same entry
  TEST_EQUAL(snap.found, true)
  TEST_EQUAL(snap.sent.size(), 1)
  TEST_REAL_SIMILAR(snap.sent.at({"b", 50}).at(8), 1000.0)

  // A newer identification: 72 (unsent) is withdrawn, 71 stays committed; its plan holds y29 again and a new b60.
  Ms3ReplanResult r2 = q.replaceMs3Plan(pkey(12351.40), {ms3Item(73, "y", 29, 4, 900.0), ms3Item(74, "b", 60, 9, 700.0)},
                                        snap.generation, 3);
  TEST_EQUAL(r2.entry_id, r1.entry_id)
  TEST_EQUAL(asSet(r2.withdrawn_scan_ids) == std::set<int>({72}), true)
  TEST_EQUAL(asSet(r2.inserted_scan_ids) == std::set<int>({73, 74}), true)
  TEST_EQUAL(r2.committed_in_gap.empty(), true)
  const std::vector<int> drained = drainIds(q);
  TEST_EQUAL(std::set<int>(drained.begin(), drained.end()) == std::set<int>({73, 74, 80}), true)   // 72 never goes out
}
END_SECTION

// The cleavage cap counts DISTINCT cleavages per proteoform per run (charges of one cleavage count once); a sent
// fragment below 2x is refused at the swap; a commit made after the planner's snapshot reconciles a stale item and
// is reported; a refused item's proteoform stash is dropped.
START_SECTION(ms3_plan_cap_two_x_and_commit_in_gap)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);

  // b50 at two charges + y29: two cleavages; then b60 makes three; b70 would be a fourth.
  Ms3ReplanResult r1 = q.replaceMs3Plan(pkey(12351.39), {ms3Item(1, "b", 50, 8, 1000.0), ms3Item(2, "b", 50, 6, 500.0),
                                                         ms3Item(3, "y", 29, 4, 800.0), ms3Item(4, "b", 60, 9, 700.0),
                                                         ms3Item(5, "b", 70, 9, 600.0)}, 0, 3);
  TEST_EQUAL(asSet(r1.inserted_scan_ids) == std::set<int>({1, 2, 3, 4}), true)
  TEST_EQUAL(asSet(r1.reconciled_scan_ids) == std::set<int>({5}), true)   // the fourth distinct cleavage

  const Ms3RegistrySnapshot before = q.ms3Snapshot(pkey(12351.39));
  ABORT_IF(!q.dequeue().has_value())   // 1 commits b50 z8 at 1000 AFTER the snapshot
  Ms3ReplanResult r2 = q.replaceMs3Plan(pkey(12351.39), {ms3Item(6, "b", 50, 8, 1500.0), ms3Item(7, "b", 50, 8, 2500.0)},
                                        before.generation, 3);
  TEST_EQUAL(r2.committed_in_gap.size(), 1)
  TEST_EQUAL(asSet(r2.reconciled_scan_ids) == std::set<int>({6}), true)   // 1.5x of the sent 1000
  TEST_EQUAL(asSet(r2.inserted_scan_ids) == std::set<int>({7}), true)     // 2.5x: a re-acquisition, not a new cleavage
  TEST_EQUAL(asSet(r2.withdrawn_scan_ids) == std::set<int>({2, 3, 4}), true)
  TEST_EQUAL(q.takeMS3Proteoform(6).empty(), true)
}
END_SECTION

// An expired command of an item that never went out withdraws the WHOLE item into the inbox the analysis thread
// undoes; an expired command of a committed item only leaves the registry.
START_SECTION(ms3_plan_expired_unsent_item_reaches_the_inbox)
{
  Config cfg{std::string(minimal_config)};   // scan_timeout 100 ms
  ScanCommandQueue q(cfg);
  Ms3PlanItem group = ms3Item(11, "b", 40, 5, 900.0);   // a two-variant exploration group
  ScanCommand v2 = group.commands[0];
  v2.scan_id = 12;
  group.commands.push_back(v2);
  group.sweep_group_id = 7;
  q.replaceMs3Plan(pkey(12351.39), {group}, 0, 3);
  std::this_thread::sleep_for(std::chrono::milliseconds(150));
  q.cleanupExpired();
  const Ms3Expired expired = q.takeMs3Expired();
  TEST_EQUAL(expired.items.size(), 1)
  ABORT_IF(expired.items.empty())
  TEST_EQUAL(expired.items[0].sweep_group_id, 7)
  TEST_EQUAL(expired.committed_owner_nominals.empty(), true)       // nothing of the entry was committed
  TEST_EQUAL(drainIds(q).empty(), true)
  TEST_EQUAL(q.ms3Snapshot(pkey(12351.39)).sent.empty(), true)   // nothing was sent
  TEST_EQUAL(q.takeMs3Expired().items.empty(), true)               // taken once
}
END_SECTION

// An expired unsent item arrives with its entry's committed (owner, nominal) pairs, so its undo keeps the stamp of a
// mass whose sibling (another charge, separate mode) went out -- including a sibling sent after the expiry.
START_SECTION(ms3_plan_expired_item_carries_its_entry_committed_nominals)
{
  Config cfg{std::string(minimal_config)};   // scan_timeout 100 ms
  ScanCommandQueue q(cfg);
  Ms3PlanItem z5 = ms3Item(21, "u", 3456, 5, 900.0, /*owner=*/4);
  z5.nominal = 3456;
  Ms3PlanItem z6 = ms3Item(22, "u", 3456, 6, 800.0, /*owner=*/4);
  z6.nominal = 3456;
  q.replaceMs3Plan(pkey(12351.39), {z5, z6}, 0, 3);
  std::this_thread::sleep_for(std::chrono::milliseconds(60));
  auto first = q.dequeue();                   // z5 (priority 1) goes out first and commits (4, 3456)
  ABORT_IF(!first.has_value())
  TEST_EQUAL(first->scan_id, 21)
  std::this_thread::sleep_for(std::chrono::milliseconds(100));
  q.cleanupExpired();                         // z6 never went out: its item reaches the inbox
  const Ms3Expired expired = q.takeMs3Expired();
  TEST_EQUAL(expired.items.size(), 1)
  ABORT_IF(expired.items.empty())
  TEST_EQUAL(expired.items[0].commands.front().scan_id, 22)
  TEST_EQUAL(expired.committed_owner_nominals.count({4, 3456}), 1)
}
END_SECTION

// Conservation under a racing drain: every MS3 id a swap was handed ends exactly once as dequeued, withdrawn,
// reconciled or still queued, and none is dequeued twice.
START_SECTION(ms3_plan_conservation_under_a_racing_drain)
{
  Config cfg{std::string(minimal_config)};
  ScanCommandQueue q(cfg);
  std::atomic<bool> stop{false};
  std::vector<int> dequeued;
  std::thread drain([&] {
    while (!stop.load())
    {
      if (auto c = q.dequeue()) dequeued.push_back(c->scan_id);
      else std::this_thread::yield();
    }
  });
  std::vector<int> built, withdrawn, reconciled;
  int id = 1000;
  for (int round = 0; round < 200; ++round)
  {
    const Ms3RegistrySnapshot snap = q.ms3Snapshot(pkey(12351.39));
    std::vector<Ms3PlanItem> items;
    for (int k = 0; k < 3; ++k)
    {
      items.push_back(ms3Item(++id, "b", 10 + (round + k) % 6, 4 + k, 100.0 * (1 + round % 7)));
      built.push_back(id);
    }
    const Ms3ReplanResult r = q.replaceMs3Plan(pkey(12351.39), std::move(items), snap.generation, 3);
    withdrawn.insert(withdrawn.end(), r.withdrawn_scan_ids.begin(), r.withdrawn_scan_ids.end());
    reconciled.insert(reconciled.end(), r.reconciled_scan_ids.begin(), r.reconciled_scan_ids.end());
  }
  stop = true;
  drain.join();
  const std::vector<int> queued = drainIds(q);
  std::map<int, int> seen;
  for (const auto& v : {dequeued, withdrawn, reconciled, queued})
    for (int x : v) ++seen[x];
  int missing = 0, twice = 0;
  for (int x : built)
  {
    if (seen[x] == 0) ++missing;
    if (seen[x] > 1) ++twice;
  }
  TEST_EQUAL(missing, 0)
  TEST_EQUAL(twice, 0)
}
END_SECTION

END_TEST
